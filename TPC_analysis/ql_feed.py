#!/usr/bin/env python3
#
# Storage Node で動かし、いちばん新しい run のファイルを追いかけて、間引いた枠を標準出力に流す。
# QL 用 PC の ql_stream.py が ssh で起動する。スクリプトは ssh の標準入力で渡すので、ノードに
# 置かなくてよい。ノードの python3 の標準ライブラリだけで動く。
#
# 全データを QL 用 PC へ送ると、トリガレートが高いときに転送が追いつかない。ここで間引けば、
# 送るのは --thin 分の 1 で済む。
#
# 流す形式は ql_replay.py と同じパケット (マジック, 側, 番号, 予約, 長さ, 枠, 枠の CRC-32)。
# ほかに、ql_stream.py が 2 台のノードの TPC の枠を N_trg の順にそろえるための知らせを流す。
#   側 'p' : ボードの進み。枠の代わりに、側 ('a' か 'c') 1 byte と、読み終えた最後の N_trg 4 byte
#   側 'r' : 追いかける run が変わった。枠の代わりに run ID
#
# 起動したときの run は、それまでに書かれた分を読み飛ばし、そこから先を流す。そのあとに始まった
# run は最初から流す。書き込み中のファイル (*.tmp) は、書き終わった枠まで流し、続きは次に読む。
# X.raw.tmp が閉じて X.raw になっても、同じファイルとして続きから読む。
#

import argparse
import os
import re
import struct
import sys
import time
import zlib


MAGIC = 0xeb905154
PACKET_HEAD = struct.Struct(">IBBHI")
PACKET_TAIL = struct.Struct(">I")

# TPC エンコーダと HA の枠の頭と尻 (どちらも同じ)
HEADER = bytes.fromhex("eb901964")
FOOTER = bytes.fromhex("75504943")

# TPC エンコーダの枠: ヘッダ + N_trg + N_time + N_gso, FADC 512 サンプル x 8 byte, 1 語,
# 飛跡の記録 20 byte x 最大 1024, フッタ。ヒットの無い枠は、3 語の直後にフッタが来る
TPC_HEAD_SIZE = 16
TPC_FADC_SIZE = 8 * 512
TPC_RECORD_SIZE = 20
TPC_MAX_RECORDS = 1024

# HA の枠の大きさ (ヒット無し, ヒットあり)
HA_FRAME_SIZES = (40, 840)

TPC_PATTERN = re.compile(r"tpc_enc_([ac])(\d+)_\d+\.raw$")
HA_PATTERN = re.compile(r"ha(\d+)_\d+\.dat$")
RUN_PATTERN = re.compile(r"\d{8}_run\d+$")
WRITING_SUFFIX = ".tmp"


def parse_args():

    parser = argparse.ArgumentParser(description="Stream thinned frames of the newest run")
    parser.add_argument("--root", required=True, help="storage directory (/data/smile_ics/storage_gX)")
    parser.add_argument("--thin", type=int, default=450,
                        help="send triggers whose N_trg is a multiple of this [triggers]")
    parser.add_argument("--poll", type=float, default=1.0, help="interval to read the files [s]")
    parser.add_argument("--from-start", action="store_true",
                        help="send the run running at startup from its beginning")

    return parser.parse_args()


def packet(side, index, payload):

    return PACKET_HEAD.pack(MAGIC, ord(side), index, 0, len(payload)) + payload + \
        PACKET_TAIL.pack(zlib.crc32(payload))


def natural_key(name):

    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name)]


def tpc_frame_end(data, start):
    """start から始まる TPC の枠の終わりを返す。書きかけなら None、形が合わなければ -1。"""

    body = start + TPC_HEAD_SIZE

    if body + 4 > len(data):
        return None

    if data[body:body + 4] == FOOTER:
        return body + 4

    first_record = body + TPC_FADC_SIZE + 4
    search = first_record

    # 飛跡の記録は 20 byte ごと。記録の中にフッタと同じ並びが出ることがあるので、位置もそろえる
    while True:
        found = data.find(FOOTER, search)
        if found < 0:
            if len(data) - first_record > TPC_RECORD_SIZE * TPC_MAX_RECORDS + 4:
                return -1
            return None
        if (found - first_record) % TPC_RECORD_SIZE == 0:
            if (found - first_record) // TPC_RECORD_SIZE > TPC_MAX_RECORDS:
                return -1
            return found + 4
        search = found + 1


def ha_frame_end(data, start):
    """start から始まる HA の枠の終わりを返す。書きかけなら None、形が合わなければ -1。"""

    for size in HA_FRAME_SIZES:
        if start + size > len(data):
            return None
        if data[start + size - 4:start + size] == FOOTER:
            return start + size

    return -1


def walk(data, frame_end, resync):
    """data の枠を先頭から順にたどり、[(始まり, 終わり)] と読み終えた位置を返す。

    書きかけの枠の手前で止める。形の合わない箇所は次のヘッダまで飛ばす。resync のとき (読み飛ばした
    位置が枠の途中かもしれない) は、見つけた枠の次がヘッダで始まることも確かめてから使う。
    """

    frames = []
    position = 0

    while position < len(data):

        if data[position:position + 4] != HEADER:
            if len(data) - position < 4:
                break
            following = data.find(HEADER, position + 1)
            if following < 0:
                position = max(position, len(data) - 3)
                break
            position = following
            continue

        stop = frame_end(data, position)

        if stop is None:
            break

        if stop < 0:
            position += 1
            continue

        if resync:
            if stop + 4 > len(data):
                break
            if data[stop:stop + 4] != HEADER:
                position += 1
                continue
            resync = False

        frames.append((position, stop))
        position = stop

    return frames, position, resync


class Feed:

    def __init__(self, args):

        self.args = args
        self.out = sys.stdout.buffer
        self.run_id = None
        self.skip_existing = not args.from_start
        self.reset()

    def reset(self):

        self.offsets = {}   # {X.raw のパス: 読み終えたバイト数}。X.raw.tmp もこのパスで覚える
        self.resync = set()
        self.last = {}      # {(側, ボード番号): 読み終えた最後の N_trg}
        # 閉じて読み終えたファイル (.tmp を外した名前)。もう増えないので、次からは開かない。
        # 1 つの run に何万ファイルもあるので、毎回すべてを開き直すとノードの負担になる
        self.done = set()
        self.kinds = {}     # {ディレクトリの中の名前: .tmp を外した名前。TPC と HA 以外は None}

    def newest_run(self):

        names = [name for name in os.listdir(self.args.root)
                 if RUN_PATTERN.match(name) and os.path.isdir(os.path.join(self.args.root, name))]

        return max(names) if names else None

    def data_files(self, directory):
        """まだ読み終えていない TPC と HA のファイルの名前を、.tmp を外した形で名前の順に返す。"""

        names = os.listdir(directory)

        for name in names:
            if name not in self.kinds:
                base = name[:-len(WRITING_SUFFIX)] if name.endswith(WRITING_SUFFIX) else name
                self.kinds[name] = base if TPC_PATTERN.search(base) or HA_PATTERN.search(base) else None

        bases = {self.kinds[name] for name in names} - {None} - self.done

        return sorted(bases, key=natural_key)

    def poll(self):

        run_id = self.newest_run()

        if run_id is None:
            return

        directory = os.path.join(self.args.root, run_id)

        if run_id != self.run_id:
            self.run_id = run_id
            self.reset()
            self.out.write(packet("r", 0, run_id.encode()))
            if self.skip_existing:
                for name in self.data_files(directory):
                    path = os.path.join(directory, name)
                    if os.path.exists(path):
                        self.done.add(name)
                        continue
                    try:
                        self.offsets[path] = os.path.getsize(path + WRITING_SUFFIX)
                    except FileNotFoundError:
                        continue
                    self.resync.add(path)
        self.skip_existing = False

        tpc = []
        progress = set()

        for name in self.data_files(directory):

            path = os.path.join(directory, name)
            closed = os.path.exists(path)
            actual = path if closed else path + WRITING_SUFFIX
            offset = self.offsets.get(path, 0)

            try:
                size = os.path.getsize(actual)
                if size <= offset:
                    if closed:
                        self.done.add(name)
                    continue
                with open(actual, "rb") as handle:
                    handle.seek(offset)
                    data = handle.read()
            except FileNotFoundError:
                # 読もうとした間に .tmp が閉じて名前が変わった。次に読む
                continue

            match = TPC_PATTERN.search(name)
            side, index = (match.group(1), int(match.group(2))) if match else \
                ("h", int(HA_PATTERN.search(name).group(1)))
            frames, used, resync = walk(data, tpc_frame_end if match else ha_frame_end,
                                        path in self.resync)

            if not resync:
                self.resync.discard(path)

            self.offsets[path] = offset + used

            for start, stop in frames:

                trig = int.from_bytes(data[start + 4:start + 8], "big")

                if match:
                    self.last[(side, index)] = max(self.last.get((side, index), trig), trig)
                    progress.add((side, index))
                    if trig % self.args.thin == 0:
                        tpc.append((trig, side, index, data[start:stop]))
                elif trig % self.args.thin == 0:
                    self.out.write(packet("h", index, data[start:stop]))

        # TPC は N_trg の順に流し、そのあとで各ボードの進みを知らせる
        for trig, side, index, frame in sorted(tpc, key=lambda item: item[0]):
            self.out.write(packet(side, index, frame))

        for side, index in sorted(progress):
            self.out.write(packet("p", index, side.encode() + self.last[(side, index)].to_bytes(4, "big")))

        self.out.flush()


def main():

    args = parse_args()
    feed = Feed(args)
    parent = os.getppid()

    try:
        # ssh が切れても、書き出すものが無いと気づけない。親 (ssh のセッション) がいなくなったら終わる
        while os.getppid() == parent:
            feed.poll()
            time.sleep(args.poll)
    except (BrokenPipeError, KeyboardInterrupt):
        # QL 用 PC 側が終わった
        pass


if __name__ == "__main__":
    main()
