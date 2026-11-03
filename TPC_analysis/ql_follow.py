#!/usr/bin/env python3
#
# QL 用 PC に同期した run のファイルを見張り、間引いたイベントを ql_tpc.py / ql_psa.py へ流す。
#
# ql_sync.sh が Storage Node から run を同期し、閉じたファイル (tpc_enc_*.raw, haNN_*.dat) から
# 順に届く。この送信器は新しく届いたファイルを読み、ql_replay.py と同じ形式で、つないできた
# 全員へ送る。ql_tpc.py と ql_psa.py は --host 127.0.0.1 で同時につないでよい。
#
# いちばん新しい run (YYYYMMDD_runNN) を追い、より新しい run が届いたらそちらへ移る。
# 途中からつないだ受け手にも、その run でそれまでに送った枠を先に送る。
#
# ql_tpc.py は、全ボードを通した N_trg の最大値から大きく戻ると新しい run とみなす。ボードごとに
# ファイルの閉じる時刻がずれ、calltt と killua の分も別々に届くので、TPC の枠は次のように送る。
#   1. ql_sync.sh が全ノードを 1 周同期するたびに更新する .ql_sync_done を見て、更新されたとき
#      だけ新しいファイルを読む。.ql_sync_done が無ければ、新しいファイルが --settle 秒届かなく
#      なってから送る
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

from check_tpc import FADC_SIZE, FILE_PATTERN, FOOTER, HEAD_SIZE, RECORD_SIZE, load_board
from ql_replay import HA_PATTERN, split_ha_frames
from ql_tpc import encode_packet


RUN_PATTERN = re.compile(r"\d{8}_run\d+$")

# ql_sync.sh が全ノードを 1 周同期するたびに更新するファイル
SYNC_MARKER = ".ql_sync_done"


def parse_args():

    parser = argparse.ArgumentParser(description="Send thinned events from synced run files")
    parser.add_argument("-d", "--data-dir", default=os.path.expanduser("~/smile_ics_ql"),
                        help="directory synced by ql_sync.sh (storage_g0/, storage_g1/, ...)")
    parser.add_argument("--run", default=None,
                        help="run ID to send (default: follow the newest run)")
    parser.add_argument("--port", type=int, default=10012, help="TCP port to listen on")
    parser.add_argument("--bind", default="127.0.0.1", help="address to listen on")
    parser.add_argument("--thin", type=int, default=450,
                        help="send triggers whose N_trg is a multiple of this [triggers]")
    parser.add_argument("--poll", type=float, default=5.0, help="interval to look for new files [s]")
    parser.add_argument("--settle", type=float, default=10.0,
                        help="wait this long after the last new file before sending TPC frames [s]")
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


def tpc_frames(key, path, thin):
    """1 ファイルから、間引いた [(N_trg, 枠のバイト列)] と最後の N_trg を返す。"""

    analysis = argparse.Namespace(base_window=30, threshold=20.0)
    board = load_board(key, [path], analysis)

    if len(board["trig"]) == 0:
        return [], None

    buf = np.fromfile(path, dtype=np.uint8)
    frames = []

    for event in np.flatnonzero(board["trig"] % thin == 0):

        start = int(board["offsets"][event])
        body = start + HEAD_SIZE

        if int.from_bytes(buf[body:body + 4].tobytes(), "big") == FOOTER:
            stop = body + 4
        else:
            stop = body + FADC_SIZE + 4 + RECORD_SIZE * int(board["n_records"][event]) + 4

        frames.append((int(board["trig"][event]), buf[start:stop].tobytes()))

    return frames, int(board["trig"].max())


def ha_frames(path, thin):
    """HA の 1 ファイルから、間引いた [(trig, 枠のバイト列)] を返す。"""

    with open(path, "rb") as handle:
        frames = split_ha_frames(handle.read())

    selected = []

    for frame in frames:

        if len(frame) < 12:
            continue

        trig = int.from_bytes(frame[4:8], "big")

        if trig % thin == 0:
            selected.append((trig, frame))

    return selected


class Follower:
    """run のファイルを読み、送る順に並べた枠を返す。"""

    def __init__(self, args):

        self.args = args
        self.run_id = None
        self.expected = set()       # 前の run にあったボード。この run でも待つ
        self.marker_mtime = None
        self.synced = False         # .ql_sync_done を見て読むか
        self.reset()

    def reset(self):

        self.seen = set()
        self.pending = []           # まだ送らない TPC の枠 [(N_trg, 側, ボード番号, 枠)]
        self.boards = {}            # {(側, ボード番号): (届いた最後の N_trg, 最新ファイルの mtime)}
        self.last_new_file = time.monotonic()
        self.run_started = time.monotonic()

    def scan(self):
        """新しいファイルを読み、送ってよい枠を [(側, 番号, 枠)] で返す。run が変われば True も返す。"""

        # 同期の途中で読むと、片方のノードの分だけを先に送ってしまう。1 周終わるまで読まない
        marker = os.path.join(self.args.data_dir, SYNC_MARKER)
        self.synced = os.path.exists(marker)

        if self.synced:
            mtime = os.path.getmtime(marker)
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

        for directory in directories:
            for path in natsorted(glob.glob(os.path.join(directory, "*"))):

                if path in self.seen:
                    continue

                name = os.path.basename(path)
                tpc = FILE_PATTERN.search(name)
                ha = HA_PATTERN.search(name)

                if tpc is None and ha is None:
                    continue

                self.seen.add(path)
                n_files += 1

                if tpc:
                    key = (tpc.group(1), int(tpc.group(2)))
                    frames, last = tpc_frames(key, path, self.args.thin)
                    self.pending.extend((trig, key[0], key[1], frame) for trig, frame in frames)
                    if last is not None:
                        previous = self.boards.get(key, (last, 0.0))[0]
                        self.boards[key] = (max(previous, last), os.path.getmtime(path))
                else:
                    index = int(ha.group(1))
                    ready.extend(("h", index, frame) for _, frame in ha_frames(path, self.args.thin))

        if n_files > 0:
            self.last_new_file = time.monotonic()
            print(f"[FILES] {run_id}: {n_files} new files", flush=True)

        ready.extend(self.release_tpc())

        return ready, switched

    def release_tpc(self):
        """全ボードの分がそろった N_trg までの TPC の枠を、N_trg の順に取り出す。"""

        if len(self.pending) == 0:
            return []

        now = time.monotonic()
        quiet = now - self.last_new_file

        if not self.synced and quiet < self.args.settle:
            return []

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
