#!/usr/bin/env python3
#
# 各 Storage Node で ql_feed.py を ssh で動かし、間引いた枠を受けて ql_tpc.py / ql_psa.py へ流す。
#
# ql_feed.py がノード側で間引くので、QL 用 PC へ届くのは --thin 分の 1 だけで済む。全データの
# 同期 (ql_sync.sh) を待たないので、トリガレートが高くても数秒の遅れで表示できる。
# 流す形式は ql_replay.py / ql_follow.py と同じで、ql_tpc.py と ql_psa.py は --host 127.0.0.1 で
# 同時につないでよい。途中からつないだ受け手にも、その run でそれまでに送った枠を先に送る。
#
# ql_tpc.py は、全ボードを通した N_trg の最大値から大きく戻ると新しい run とみなす。2 台のノードの
# 枠は別々に届くので、TPC の枠は全ボードの進みがそろった N_trg まで待ち、N_trg の順に送る。
# 進みが --stale 秒止まったボードは待たない。run が始まって全ノードの進みがそろうまでも待つが、
# --stale 秒たてば待たない。HA は ql_psa.py が HA ごとに判定するので、届いた順に送る。
#
# ノードの ql_feed.py は、CPU とディスクの読み出しをどちらも最低の優先度 (nice 19, ionice idle) で
# 動かし、DAQ を妨げない。ノードが忙しいときは QL の方が遅れる。
# ノードには ~/.ssh/config の Host 名でつなぐ (ql_sync.sh と同じ。magallan を ProxyJump で経由する)。
# 接続が切れたら --reconnect 秒ごとにつなぎ直す。
#

import argparse
import asyncio
import ctypes
import signal
import struct
import time
from pathlib import Path

from ql_follow import Broadcaster


MAGIC = 0xeb905154
PACKET_HEAD = struct.Struct(">IBBHI")
PACKET_TAIL = struct.Struct(">I")

FEED_SCRIPT = Path(__file__).resolve().parent / "ql_feed.py"


def die_with_parent():
    """起動する ssh に、親 (このプロセス) が終わったら自分も終わるよう Linux に頼む。

    ssh が残ると、ノードの ql_feed.py も動き続ける。Ctrl-C だけでなく、kill されたり落ちたりした
    ときにも ssh が残らないようにする。
    """

    PR_SET_PDEATHSIG = 1
    ctypes.CDLL(None, use_errno=True).prctl(PR_SET_PDEATHSIG, signal.SIGTERM)


def parse_args():

    parser = argparse.ArgumentParser(description="Merge thinned frames streamed from storage nodes")
    parser.add_argument("--nodes", default="calltt:storage_g0 killua:storage_g1",
                        help="ssh host and storage process id pairs (same as QL_NODES of ql_sync.sh)")
    parser.add_argument("--root", default="/data/smile_ics", help="storage root on the nodes")
    parser.add_argument("--thin", type=int, default=450,
                        help="send triggers whose N_trg is a multiple of this [triggers]")
    parser.add_argument("--poll", type=float, default=1.0, help="interval the nodes read the files [s]")
    parser.add_argument("--from-start", action="store_true",
                        help="send the run running at startup from its beginning "
                             "(default: from the moment of startup)")
    parser.add_argument("--stale", type=float, default=30.0,
                        help="do not wait for boards or nodes silent this long [s]")
    parser.add_argument("--reconnect", type=float, default=5.0, help="reconnect interval [s]")
    parser.add_argument("--port", type=int, default=10012, help="TCP port to listen on")
    parser.add_argument("--bind", default="127.0.0.1", help="address to listen on")
    parser.add_argument("--history-mb", type=float, default=256.0,
                        help="packets kept for viewers that connect later [MB]")

    return parser.parse_args()


class Merger:
    """ノードごとの流れを合わせ、TPC の枠を N_trg の順にそろえて送る。"""

    def __init__(self, args, broadcaster, nodes):

        self.args = args
        self.broadcaster = broadcaster
        self.nodes = nodes
        self.node_run = {}          # {ノード: そのノードが今追いかけている run}
        self.run_id = None
        self.reset()

    def reset(self):

        self.pending = []           # まだ送らない TPC の枠 [(N_trg, パケット)]
        self.boards = {}            # {(側, ボード番号): (進んだ N_trg, 知らせを受けた時刻)}
        self.ready_nodes = set()    # この run で進みを知らせてきたノード
        self.run_started = time.monotonic()

    async def handle(self, node, side, index, payload, raw):

        if side == "r":
            run_id = payload.decode()
            self.node_run[node] = run_id
            if self.run_id is None or run_id > self.run_id:
                print(f"[RUN] {self.run_id} -> {run_id} ({node})", flush=True)
                self.run_id = run_id
                self.reset()
                self.broadcaster.clear()
            return

        # 古い run を追いかけているノードの枠は捨てる (すぐに新しい run へ移ってくる)
        if self.node_run.get(node) != self.run_id:
            return

        if side == "h":
            await self.broadcaster.send([raw])
        elif side in "ac":
            self.pending.append((int.from_bytes(payload[4:8], "big"), raw))
        elif side == "p":
            board = (chr(payload[0]), index)
            self.boards[board] = (int.from_bytes(payload[1:5], "big"), time.monotonic())
            self.ready_nodes.add(node)

    def release(self):
        """全ボードの進みがそろった N_trg までの TPC の枠を、N_trg の順に取り出す。"""

        if len(self.pending) == 0:
            return []

        now = time.monotonic()
        waiting_nodes = now - self.run_started < self.args.stale and self.ready_nodes != set(self.nodes)
        active = [trig for trig, seen in self.boards.values() if now - seen < self.args.stale]

        if waiting_nodes:
            return []

        limit = min(active) if active else float("inf")

        self.pending.sort(key=lambda item: item[0])
        count = 0

        while count < len(self.pending) and self.pending[count][0] <= limit:
            count += 1

        released, self.pending = self.pending[:count], self.pending[count:]

        return [raw for _, raw in released]

    async def release_loop(self):

        while True:
            await asyncio.sleep(0.5)
            packets = self.release()
            if packets:
                await self.broadcaster.send(packets)


async def follow_node(node, process_id, merger, args, script):
    """ノードで ql_feed.py を動かし、届いたパケットを merger に渡す。切れたらつなぎ直す。"""

    command = (f"nice -n 19 ionice -c3 python3 - --root {args.root}/{process_id} --thin {args.thin} "
               f"--poll {args.poll}" + (" --from-start" if args.from_start else ""))
    # つなぎ直したときに最初から送り直さないよう、--from-start は最初の 1 回だけにする
    first = True

    while True:

        process = await asyncio.create_subprocess_exec(
            "ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "-o", "ServerAliveInterval=10",
            node, command if first else command.replace(" --from-start", ""),
            stdin=asyncio.subprocess.PIPE, stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE, preexec_fn=die_with_parent)
        first = False

        process.stdin.write(script)
        await process.stdin.drain()
        process.stdin.close()
        print(f"[FEED] {node}:{process_id} started", flush=True)

        try:
            while True:
                head = await process.stdout.readexactly(PACKET_HEAD.size)
                magic, side, index, _, size = PACKET_HEAD.unpack(head)
                if magic != MAGIC:
                    raise ConnectionError(f"wrong packet header: {head.hex()}")
                payload = await process.stdout.readexactly(size)
                tail = await process.stdout.readexactly(PACKET_TAIL.size)
                await merger.handle(node, chr(side), index, payload, head + payload + tail)
        except (asyncio.IncompleteReadError, ConnectionError) as exc:
            if process.returncode is None:
                process.kill()
            error = (await process.stderr.read()).decode(errors="replace").strip().splitlines()
            print(f"[FEED] {node}:{process_id} stopped ({exc.__class__.__name__}). "
                  + (error[-1] if error else ""), flush=True)

        await process.wait()
        await asyncio.sleep(args.reconnect)


async def run(args):

    nodes = [entry.rsplit(":", 1) for entry in args.nodes.split()]
    broadcaster = Broadcaster(int(args.history_mb * 1024 * 1024))
    merger = Merger(args, broadcaster, [node for node, _ in nodes])
    script = FEED_SCRIPT.read_bytes()
    server = await asyncio.start_server(broadcaster.handle, args.bind, args.port)

    print(f"[STREAM] {args.nodes}, listening on {args.bind}:{args.port}", flush=True)

    async with server:
        await asyncio.gather(merger.release_loop(),
                             *(follow_node(node, process_id, merger, args, script)
                               for node, process_id in nodes))


def main():

    try:
        asyncio.run(run(parse_args()))
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
