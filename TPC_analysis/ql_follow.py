#!/usr/bin/env python3
#
# QL 用 PC に同期した run のファイルを見張り、間引いたイベントを ql_tpc.py / ql_psa.py へ流す。
#
# ql_sync.sh が Storage Node から run を同期する。書き込み中のファイル (*.tmp) も、書き足された分
# ずつ届く。この送信器はファイルごとにどこまで読んだかを覚え、増えた分のうち書き終わった枠
# だけを読み、ql_replay.py と同じ形式で、つないできた全員へ送る。ql_tpc.py と ql_psa.py は
# --host 127.0.0.1 で同時につないでよい。
#
# 書き込み中のファイルは、末尾の書きかけの枠の手前で止め、続きは次に届いたときに読む。
# X.raw.tmp が閉じて X.raw になっても同じファイルとして続きから読み、同じ枠を二度送らない。
# 途中に壊れた箇所があったときに読み飛ばすのは、閉じたファイルだけにする。書き込み中の
# ファイルでは、書きかけの枠と見分けがつかないため。
#
# いちばん新しい run (YYYYMMDD_runNN) を追い、より新しい run が届いたらそちらへ移る。
# 起動したときにすでに届いている分は読み飛ばし、そこから先を読む (長い run の途中で起動しても
# すぐに表示が始まる)。run の最初から読むときは --from-start を付ける。起動したあとに始まった
# run は最初から読む。
# 途中からつないだ受け手にも、その run でそれまでに送った枠を先に送る。
#
# ql_tpc.py は、全ボードを通した N_trg の最大値から大きく戻ると新しい run とみなす。ボードごとに
# ファイルの閉じる時刻がずれ、calltt と killua の分も別々に届くので、TPC の枠は次のように送る。
#   1. ql_sync.sh が全ノードを 1 周同期するたびに更新する .ql_sync_done を見て、更新されたとき
#      だけ新しいファイルを読む。起動する前の .ql_sync_done は、同期の途中かもしれないので使わず、
#      起動したあとの 1 周目が終わるまで待つ。ql_sync.sh を使わずにファイルを置いたときは、
#      置き終えてから touch .ql_sync_done する
#   2. 全ボードの分がそろった N_trg までを、N_trg の順に送る。前の run にあったボードは、
#      この run でまだ届いていなくても待つ
#   3. 止まったボードを待ち続けないよう、最新のファイルが他より --stale 秒以上古いボードと、
#      run が始まって --stale 秒たっても届かないボードは待たない。新しいファイルが --stale 秒
#      届かなければ、残りをすべて送る
# HA は ql_psa.py が HA ごとに判定するので、届いたファイルから順に送る。
#

import argparse
import asyncio
import glob
import os
import re
import time
from collections import deque

import numpy as np
from natsort import natsorted

from check_tpc import (
    FADC_MAX,
    FADC_SIZE,
    FILE_PATTERN,
    FOOTER,
    HEAD_SIZE,
    MIN_FRAME_SIZE,
    N_CLOCKS,
    N_FADC_CH,
    N_FADC_SAMPLES,
    N_STRIPS,
    RECORD_SIZE,
    decode_file,
)
from ql_replay import HA_FOOTER, HA_HEADER, HA_PATTERN
from ql_tpc import encode_packet


RUN_PATTERN = re.compile(r"\d{8}_run\d+$")

# ql_sync.sh が全ノードを 1 周同期するたびに更新するファイル
SYNC_MARKER = ".ql_sync_done"

# Storage Node が書き込み中のファイルに付ける拡張子
WRITING_SUFFIX = ".tmp"


def parse_args():

    parser = argparse.ArgumentParser(description="Send thinned events from synced run files")
    parser.add_argument("-d", "--data-dir", default=os.path.expanduser("~/smile_ics_ql"),
                        help="directory synced by ql_sync.sh (storage_g0/, storage_g1/, ...)")
    parser.add_argument("--run", default=None,
                        help="run ID to send (default: follow the newest run)")
    parser.add_argument("--port", type=int, default=10012, help="TCP port to listen on")
    parser.add_argument("--bind", default="127.0.0.1", help="address to listen on")
    parser.add_argument("--from-start", action="store_true",
                        help="read the run running at startup from its beginning "
                             "(default: skip what has already arrived)")
    parser.add_argument("--thin", type=int, default=450,
                        help="send triggers whose N_trg is a multiple of this [triggers]")
    parser.add_argument("--poll", type=float, default=5.0, help="interval to look for new files [s]")
    parser.add_argument("--stale", type=float, default=120.0,
                        help="do not wait for boards whose files are this much older [s]")
    parser.add_argument("--history-mb", type=float, default=256.0,
                        help="packets kept for viewers that connect later [MB]")

    return parser.parse_args()


def newest_run(data_dir):
    """data_dir/storage_*/ の下でいちばん新しい run ID を返す。無ければ None。"""

    names = [os.path.basename(path) for path in glob.glob(os.path.join(data_dir, "*", "*"))
             if os.path.isdir(path) and RUN_PATTERN.match(os.path.basename(path))]

    return max(names) if names else None


def tpc_chunk(side, data, closed, thin):
    """TPC のファイルの続き data から、書き終わった枠を読む。

    間引いた [(N_trg, 枠のバイト列)]、読めた最後の N_trg、読み終えたバイト数を返す。
    """

    buf = np.frombuffer(data, dtype=np.uint8)
    capacity = len(buf) // MIN_FRAME_SIZE + 1

    trig = np.zeros(capacity, dtype=np.uint32)
    moment = np.zeros(capacity, dtype=np.uint32)
    charge = np.zeros((capacity, N_FADC_CH), dtype=np.float32)
    n_records = np.zeros(capacity, dtype=np.int32)
    offsets = np.zeros(capacity, dtype=np.uint32)
    # decode_file が積算する図の分。ここでは使わない
    image = np.zeros((N_CLOCKS, N_STRIPS), dtype=np.uint32)
    baseline = np.zeros((N_FADC_CH, FADC_MAX + 1), dtype=np.uint32)
    wave_map = np.zeros((FADC_MAX + 1, N_FADC_SAMPLES), dtype=np.uint32)

    n_events, _ = decode_file(buf, image, baseline, wave_map, trig, moment, charge, n_records, offsets,
                              30, 20.0, side == "c")

    frames = []
    last = None
    position = 0

    for event in range(n_events):

        start = int(offsets[event])

        # 枠はすき間なく並ぶ。すき間は、閉じたファイルなら壊れた箇所なので読み飛ばし、
        # 書き込み中なら書きかけの枠かもしれないので、そこで止める
        if start != position and not closed:
            break

        body = start + HEAD_SIZE

        if int.from_bytes(data[body:body + 4], "big") == FOOTER:
            stop = body + 4
        else:
            stop = body + FADC_SIZE + 4 + RECORD_SIZE * int(n_records[event]) + 4

        if stop > len(buf):
            break

        position = stop
        number = int(trig[event])
        last = number if last is None else max(last, number)

        if number % thin == 0:
            frames.append((number, data[start:stop]))

    return frames, last, position


def ha_chunk(data, closed, thin):
    """HA のファイルの続き data から、書き終わった枠を読む。

    間引いた [(trig, 枠のバイト列)] と、読み終えたバイト数を返す。
    """

    frames = []
    position = 0

    while position < len(data):

        if data[position:position + 4] != HA_HEADER:
            # 閉じたファイルなら壊れた箇所なので次のヘッダまで飛ばし、書き込み中なら止める
            if not closed:
                break
            following = data.find(HA_HEADER, position + 1)
            if following < 0:
                break
            position = following
            continue

        stop = data.find(HA_FOOTER, position + 4)

        if stop < 0:
            break

        frame = data[position:stop + 4]
        position = stop + 4

        if len(frame) >= 12 and int.from_bytes(frame[4:8], "big") % thin == 0:
            frames.append((int.from_bytes(frame[4:8], "big"), frame))

    return frames, position


def read_more(path, offset):
    """path の offset から後ろを読む。"""

    with open(path, "rb") as handle:
        handle.seek(offset)
        return handle.read()


class Follower:
    """run のファイルを読み、送る順に並べた枠を返す。"""

    def __init__(self, args):

        self.args = args
        self.run_id = None
        self.expected = set()       # 前の run にあったボード。この run でも待つ
        self.marker_mtime = None
        self.started = time.time()  # これより前の .ql_sync_done は使わない
        self.waiting = False        # 1 周目の同期を待っていることを表示したか
        self.skip_existing = not args.from_start   # 起動したときの run は、届いている分を読み飛ばす
        self.reset()

    def reset(self):

        self.offsets = {}           # {X.raw のパス: 読み終えたバイト数}。X.raw.tmp もこのパスで覚える
        self.resync = set()         # 読み飛ばした位置が枠の途中かもしれないので、次の枠の頭を探すファイル
        self.pending = []           # まだ送らない TPC の枠 [(N_trg, 側, ボード番号, 枠)]
        self.boards = {}            # {(側, ボード番号): (届いた最後の N_trg, 最新ファイルの mtime)}
        self.last_new_file = time.monotonic()
        self.run_started = time.monotonic()

    def scan(self):
        """新しいファイルを読み、送ってよい枠を [(側, 番号, 枠)] で返す。run が変われば True も返す。"""

        # 起動したときに届いている分は、すぐに位置だけ覚えて読み飛ばす。1 周目の同期が終われば、
        # 起動したあとに増えた分から表示できる。同期先が空なら、1 周目の後で読み飛ばす
        # ql_sync.sh は転送の前に run のディレクトリを作るので、中身が無ければ空とみなす
        if self.skip_existing:
            run_id = self.args.run or newest_run(self.args.data_dir)
            directories = sorted(glob.glob(os.path.join(self.args.data_dir, "*", run_id or "")))
            if run_id is not None and self.received_bytes(directories) > 0:
                self.skip_existing = False
                self.run_id = run_id
                print(f"[RUN] None -> {run_id}", flush=True)
                self.skip_to_end(directories, run_id)

        # 同期の途中で読むと、送られてきた順 (ファイルの名前の順) に 1 つずつ表示されたり、片方の
        # ノードの分だけを先に送ったりする。1 周終わるまで読まない
        marker = os.path.join(self.args.data_dir, SYNC_MARKER)
        mtime = os.path.getmtime(marker) if os.path.exists(marker) else None

        if mtime is None or mtime < self.started:
            if not self.waiting:
                self.waiting = True
                print(f"[WAIT] waiting for ql_sync.sh to finish a round ({marker})", flush=True)
            return [], False

        if mtime == self.marker_mtime:
            return self.release_tpc(), False

        self.marker_mtime = mtime

        run_id = self.args.run or newest_run(self.args.data_dir)
        switched = run_id != self.run_id

        if switched:
            print(f"[RUN] {self.run_id} -> {run_id}", flush=True)
            self.run_id = run_id
            self.expected = set(self.boards) or self.expected
            self.reset()

        if run_id is None:
            return [], switched

        directories = sorted(glob.glob(os.path.join(self.args.data_dir, "*", run_id)))
        ready = []
        n_files = 0

        if self.skip_existing:
            self.skip_existing = False
            self.skip_to_end(directories, run_id)
            return [], switched

        for directory in directories:
            for name in natsorted(self.data_files(directory)):

                path = os.path.join(directory, name)
                # 閉じたあとの名前があればそちらを読む。無ければ書き込み中の .tmp を読む
                closed = os.path.exists(path)
                actual = path if closed else path + WRITING_SUFFIX
                offset = self.offsets.get(path, 0)

                try:
                    if os.path.getsize(actual) <= offset:
                        continue
                    data = read_more(actual, offset)
                    mtime = os.path.getmtime(actual)
                except FileNotFoundError:
                    # 読もうとした間に .tmp が閉じて名前が変わった。次に読む
                    continue

                tpc = FILE_PATTERN.search(name)

                if tpc:
                    key = (tpc.group(1), int(tpc.group(2)))
                    frames, last, used = tpc_chunk(key[0], data, closed or path in self.resync,
                                                   self.args.thin)
                    self.pending.extend((trig, key[0], key[1], frame) for trig, frame in frames)
                    if last is not None:
                        previous = self.boards.get(key, (last, 0.0))[0]
                        self.boards[key] = (max(previous, last), mtime)
                else:
                    index = int(HA_PATTERN.search(name).group(1))
                    frames, used = ha_chunk(data, closed or path in self.resync, self.args.thin)
                    ready.extend(("h", index, frame) for _, frame in frames)

                if used > 0:
                    self.offsets[path] = offset + used
                    self.resync.discard(path)
                    n_files += 1

        if n_files > 0:
            self.last_new_file = time.monotonic()
            print(f"[FILES] {run_id}: {n_files} files grew", flush=True)

        ready.extend(self.release_tpc())

        return ready, switched

    def skip_to_end(self, directories, run_id):
        """届いている分を読み飛ばし、次からはその続きを読む。"""

        n_files = 0
        n_bytes = 0

        for directory in directories:
            for name in self.data_files(directory):
                path = os.path.join(directory, name)
                actual = path if os.path.exists(path) else path + WRITING_SUFFIX
                try:
                    self.offsets[path] = os.path.getsize(actual)
                except FileNotFoundError:
                    continue
                # 書き込み中のファイルの末尾は枠の途中かもしれない。次は枠の頭を探して読む
                self.resync.add(path)
                n_files += 1
                n_bytes += self.offsets[path]

        print(f"[SKIP] {run_id}: skip {n_files} files ({n_bytes / 1e6:.1f} MB) already received. "
              "use --from-start to read them", flush=True)

    def received_bytes(self, directories):
        """directories に届いている TPC と HA のファイルの合計の大きさ。"""

        total = 0

        for directory in directories:
            for name in os.listdir(directory):
                base = name[:-len(WRITING_SUFFIX)] if name.endswith(WRITING_SUFFIX) else name
                if FILE_PATTERN.search(base) or HA_PATTERN.search(base):
                    total += os.path.getsize(os.path.join(directory, name))

        return total

    @staticmethod
    def data_files(directory):
        """directory の TPC と HA のファイルの名前を、書き込み中の .tmp を外した形で返す。"""

        names = set()

        for name in os.listdir(directory):
            base = name[:-len(WRITING_SUFFIX)] if name.endswith(WRITING_SUFFIX) else name
            if FILE_PATTERN.search(base) or HA_PATTERN.search(base):
                names.add(base)

        return names

    def release_tpc(self):
        """全ボードの分がそろった N_trg までの TPC の枠を、N_trg の順に取り出す。"""

        if len(self.pending) == 0:
            return []

        now = time.monotonic()
        quiet = now - self.last_new_file

        if quiet >= self.args.stale:
            limit = float("inf")
        elif self.expected - self.boards.keys() and now - self.run_started < self.args.stale:
            return []
        else:
            newest = max(mtime for _, mtime in self.boards.values())
            limit = min(last for last, mtime in self.boards.values()
                        if newest - mtime < self.args.stale)

        self.pending.sort(key=lambda item: item[0])
        count = 0

        while count < len(self.pending) and self.pending[count][0] <= limit:
            count += 1

        released, self.pending = self.pending[:count], self.pending[count:]

        return [(side, index, frame) for _, side, index, frame in released]


class Broadcaster:
    """つないできた受け手全員へ送り、途中からの受け手のために送った分を持っておく。"""

    def __init__(self, history_bytes):

        self.history_bytes = history_bytes
        self.clients = set()
        self.clear()

    def clear(self):

        self.history = deque()
        self.n_bytes = 0

    async def send(self, packets):

        for packet in packets:
            self.history.append(packet)
            self.n_bytes += len(packet)

        while self.n_bytes > self.history_bytes and len(self.history) > 0:
            self.n_bytes -= len(self.history.popleft())

        for writer in list(self.clients):
            try:
                for packet in packets:
                    writer.write(packet)
                await writer.drain()
            except ConnectionError:
                self.clients.discard(writer)

    async def handle(self, reader, writer):

        peer = writer.get_extra_info("peername")
        print(f"[CONNECT] {peer}, sending {len(self.history)} packets so far", flush=True)

        try:
            for packet in list(self.history):
                writer.write(packet)
            await writer.drain()
            self.clients.add(writer)
            # 受け手は何も送ってこない。切れるまで待つ
            await reader.read()
        except ConnectionError:
            pass
        finally:
            self.clients.discard(writer)
            writer.close()
            print(f"[DISCONNECT] {peer}", flush=True)


async def run(args):

    follower = Follower(args)
    broadcaster = Broadcaster(int(args.history_mb * 1024 * 1024))
    server = await asyncio.start_server(broadcaster.handle, args.bind, args.port)

    print(f"[FOLLOW] {args.data_dir}, listening on {args.bind}:{args.port}", flush=True)

    async with server:
        while True:

            # 読むのに時間がかかるので、受け手への送信を止めないよう別スレッドで読む
            frames, switched = await asyncio.to_thread(follower.scan)

            if switched:
                broadcaster.clear()

            if len(frames) > 0:
                await broadcaster.send([encode_packet(side, index, frame)
                                        for side, index, frame in frames])

            await asyncio.sleep(args.poll)


def main():

    try:
        asyncio.run(run(parse_args()))
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
