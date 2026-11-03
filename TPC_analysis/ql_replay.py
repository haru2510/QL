#!/usr/bin/env python3
#
# 保存済みの raw を間引いて TCP で流す、ql_tpc.py / ql_psa.py の試験用の送信器。
# zeno の中継が無くても QL を確認できる。接続を待ち、run を最後まで流したら終わる。
#
# トリガ番号が thin の倍数のイベントを、全ボードの分だけ記録時の間隔で送る。
# run ディレクトリに HA の haNN_*.dat があれば、それも側 'h' で一緒に送る。
# TPC は N_time、HA は 1 MHz の clk を、それぞれ最初のイベントからの秒に直して並べる。
#

import argparse
import asyncio
import glob
import os
import re
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from natsort import natsorted

from check_tpc import FADC_SIZE, FOOTER, HEAD_SIZE, RECORD_SIZE, collect_boards, load_board
from ql_tpc import encode_packet


HA_PATTERN = re.compile(r"ha(\d+)_\d+\.dat$")
HA_HEADER = bytes.fromhex("eb901964")
HA_FOOTER = bytes.fromhex("75504943")

# HA の clk は 1 MHz の 32 bit カウンタ
HA_CLOCK = 1e6
HA_CLOCK_WRAP = 2 ** 32


def parse_args():

    parser = argparse.ArgumentParser(description="Replay a saved run as thinned events")
    parser.add_argument("-i", "--input", required=True,
                        help="run directory (tpc_enc_*.raw and/or haNN_*.dat)")
    parser.add_argument("--port", type=int, default=10012, help="TCP port to listen on")
    parser.add_argument("--thin", type=int, default=450,
                        help="send triggers whose N_trg is a multiple of this [triggers]")
    parser.add_argument("--speed", type=float, default=1.0, help="replay speed [times real time]")
    parser.add_argument("--time-tick", type=float, default=100e-6, help="N_time tick [s]")
    parser.add_argument("--workers", type=int, default=4, help="boards decoded at the same time")

    return parser.parse_args()


def load_tpc_events(directory, thin, args):
    """TPC の [(秒, N_trg, 側, ボード番号, 枠のバイト列)] を返す。"""

    files = collect_boards(directory, None)

    if len(files) == 0:
        return []

    analysis = argparse.Namespace(base_window=30, threshold=20.0)

    with ThreadPoolExecutor(args.workers) as pool:
        boards = list(pool.map(lambda item: load_board(item[0], item[1], analysis), files.items()))

    origin = min(int(board["time"].min()) for board in boards if len(board["time"]) > 0)
    events = []

    for (side, index), board in zip(files.keys(), boards):

        buffers = [np.fromfile(path, dtype=np.uint8) for path in board["paths"]]

        for event in np.flatnonzero(board["trig"] % thin == 0):

            buf = buffers[board["file_index"][event]]
            start = int(board["offsets"][event])
            body = start + HEAD_SIZE

            if int.from_bytes(buf[body:body + 4].tobytes(), "big") == FOOTER:
                stop = body + 4
            else:
                stop = body + FADC_SIZE + 4 + RECORD_SIZE * int(board["n_records"][event]) + 4

            seconds = (int(board["time"][event]) - origin) * args.time_tick
            events.append((seconds, int(board["trig"][event]), side, index, buf[start:stop].tobytes()))

    return events


def split_ha_frames(data):
    """HA の .dat をヘッダからフッタまでの枠に分ける。"""

    frames = []
    start = data.find(HA_HEADER)

    while start >= 0:

        stop = data.find(HA_FOOTER, start + 4)

        if stop < 0:
            break

        frames.append(data[start:stop + 4])
        start = data.find(HA_HEADER, stop + 4)

    return frames


def load_ha_events(directory, thin):
    """HA の [(秒, trig, 'h', HA 番号, 枠のバイト列)] を返す。"""

    paths = {}

    for path in glob.glob(os.path.join(directory, "ha*_*.dat")):
        match = HA_PATTERN.search(os.path.basename(path))
        if match:
            paths.setdefault(int(match.group(1)), []).append(path)

    events = []

    for index, files in sorted(paths.items()):

        origin = None
        previous = None
        elapsed = 0

        for path in natsorted(files):

            with open(path, "rb") as handle:
                frames = split_ha_frames(handle.read())

            for frame in frames:

                if len(frame) < 12:
                    continue

                trig = int.from_bytes(frame[4:8], "big")
                clock = int.from_bytes(frame[8:12], "big")

                # clk の 1 周 (約 72 分) をまたいでも時刻が戻らないようにする
                if previous is not None:
                    elapsed += (clock - previous) % HA_CLOCK_WRAP
                previous = clock

                if origin is None:
                    origin = elapsed

                if trig % thin == 0:
                    events.append(((elapsed - origin) / HA_CLOCK, trig, "h", index, frame))

    return events


async def run(args):

    tpc = load_tpc_events(args.input, args.thin, args)
    ha = load_ha_events(args.input, args.thin)
    events = sorted(tpc + ha, key=lambda item: item[:2])

    print(f"[REPLAY] TPC {len(tpc)} frames, HA {len(ha)} frames")

    if len(events) == 0:
        raise FileNotFoundError(f"no tpc_enc_*.raw or haNN_*.dat in {args.input}")

    finished = asyncio.Event()

    async def send(_reader, writer):

        previous = None

        try:
            for seconds, _, side, index, frame in events:

                if previous is not None and seconds > previous:
                    await asyncio.sleep((seconds - previous) / args.speed)

                previous = seconds
                writer.write(encode_packet(side, index, frame))
                await writer.drain()
        except ConnectionError as exc:
            print(f"[REPLAY] connection lost: {exc}")
        finally:
            writer.close()
            finished.set()

    server = await asyncio.start_server(send, "127.0.0.1", args.port)

    async with server:
        await finished.wait()

    print("[REPLAY] finished")


def main():

    asyncio.run(run(parse_args()))


if __name__ == "__main__":

    main()
