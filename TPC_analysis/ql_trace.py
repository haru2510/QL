#!/usr/bin/env python3
#
# QL に流れる枠を記録し、trig (TPC は N_trg) が --threshold 以上戻ったところを、直前の様子と一緒に出す。
# QL で「trig went back」が出るときに、どこで順番が崩れているかを見分ける。
#
#   python3 TPC_analysis/ql_trace.py --port 10012 --seconds 120
#       動いている ql_stream.py につなぎ、ビューアに流れているものを記録する
#   python3 TPC_analysis/ql_trace.py --node calltt:storage_g0 --seconds 120
#       ql_stream.py を通さず、ノードで ql_feed.py を直接動かして記録する。直前に読んだファイルも出す
#
# --port で戻りが出て --node で出なければ、ql_stream.py (2 台の合わせ方) の問題。--node でも出れば
# ql_feed.py (ノードでの読み方) の問題。
#

import argparse
import asyncio
import struct
import time
from collections import defaultdict, deque
from pathlib import Path

from ql_stream import die_with_parent


PACKET_HEAD = struct.Struct(">IBBHI")
PACKET_TAIL = struct.Struct(">I")
FEED_SCRIPT = Path(__file__).resolve().parent / "ql_feed.py"


def parse_args():

    parser = argparse.ArgumentParser(description="Record QL frames and report trig going back")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--port", type=int, help="TCP port of ql_stream.py on localhost")
    source.add_argument("--node", help="ssh host and storage process id (calltt:storage_g0)")
    parser.add_argument("--root", default="/data/smile_ics", help="storage root on the node")
    parser.add_argument("--seconds", type=float, default=120.0, help="how long to record [s]")
    parser.add_argument("--threshold", type=int, default=1000, help="report going back more than this")

    return parser.parse_args()


async def open_source(args):
    """(reader, 後始末) を返す。"""

    if args.port:
        reader, writer = await asyncio.open_connection("127.0.0.1", args.port)
        return reader, writer.close

    node, process_id = args.node.rsplit(":", 1)
    process = await asyncio.create_subprocess_exec(
        "ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", node,
        f"nice -n 19 ionice -c3 python3 - --root {args.root}/{process_id} --debug",
        stdin=asyncio.subprocess.PIPE, stdout=asyncio.subprocess.PIPE, preexec_fn=die_with_parent)
    process.stdin.write(FEED_SCRIPT.read_bytes())
    await process.stdin.drain()
    process.stdin.close()

    return process.stdout, process.kill


async def run(args):

    reader, close = await open_source(args)
    started = time.monotonic()
    latest = {}
    counts = defaultdict(int)
    recent = defaultdict(lambda: deque(maxlen=4))   # {枠の種類: 直前の (時刻, trig)}
    reads = defaultdict(lambda: deque(maxlen=8))     # {枠の種類: そのファイルを直前に読んだ様子} (--node のとき)
    drops = 0

    print(f"[TRACE] {args.port and f'port {args.port}' or args.node}, {args.seconds:g} s", flush=True)

    try:
        while time.monotonic() - started < args.seconds:

            head = await asyncio.wait_for(reader.readexactly(PACKET_HEAD.size),
                                          args.seconds - (time.monotonic() - started))
            _, side, index, _, size = PACKET_HEAD.unpack(head)
            payload = await reader.readexactly(size)
            await reader.readexactly(PACKET_TAIL.size)
            side = chr(side)
            now = time.monotonic() - started

            if side == "d":
                name = payload.decode().split()[0]
                tpc = name.startswith("tpc_enc_")
                kind = f"{name[8]}{name[9:11]}" if tpc else f"h{name[2:4]}"
                reads[kind].append(f"{now:6.1f} s  {payload.decode()}")
                continue
            if side == "r":
                print(f"{now:6.1f} s  run {payload.decode()}", flush=True)
                latest.clear()
                continue
            if side not in "ach":
                continue

            key = f"{side}{index:02d}"
            trig = int.from_bytes(payload[4:8], "big")
            counts[key] += 1

            if key in latest and trig + args.threshold < latest[key]:
                drops += 1
                if drops <= 10:
                    history = ", ".join(f"{t:.1f} s {v}" for t, v in recent[key])
                    print(f"{now:6.1f} s  {key}: {latest[key]} -> {trig} (差 {latest[key] - trig})"
                          f"  直前: {history}", flush=True)
                    for line in reads[key]:
                        print(f"           {key} のファイルを読んだ: {line}  (名前 読み始め 読んだ量 閉じていたか)",
                              flush=True)

            latest[key] = max(latest.get(key, trig), trig)
            recent[key].append((now, trig))

    except (asyncio.TimeoutError, asyncio.IncompleteReadError):
        pass
    finally:
        close()

    print(f"[TRACE] 枠 {sum(counts.values())} ({len(counts)} 種類), {args.threshold} 以上の戻り {drops} 回",
          flush=True)


def main():

    asyncio.run(run(parse_args()))


if __name__ == "__main__":
    main()
