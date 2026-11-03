#!/usr/bin/env python3
#
# zeno から届く間引きイベントを受け、check_tpc.py と同じ確認図を積算しながら更新する。
#
# zeno が無いときは ql_replay.py が同じ形式で送る。起動するとブラウザが開き、
# 図は自動で更新される。別の PC から見るときは --http-bind 0.0.0.0 を付ける。
#
# ブラウザには次の 3 枚を、モニターに収まるように並べる。
#   1 イベントの図 : 飛跡と FADC 波形。--event-update ごとに最新のイベントへ更新
#   レートの図     : カウントレート、FADC ベースライン、dE/dx (全 uPIC の飛跡長とアノード電荷の 2 次元図)
#   積算の図       : 飛跡の積算画像。上の段が uPIC 0-3 の anode、下の段が cathode。画面の大半を使う
#
# 1 イベントの図の色は、x が N 側 (uPIC 0, 3) 赤・S 側 (1, 2) 青、
# y が E 側 (uPIC 0, 1) 赤・W 側 (2, 3) 青。FADC は 1 ch あたり 32 ストリップなので、
# その中央に波形を縦向きに重ねる。縦軸は 1 サンプル = 2 クロックとしている。
#
# 受ける形式 (ビッグエンディアン。TPC エンコーダの 1 枠を包んだもの):
#   マジック 4 byte, 側 1 byte ('a' か 'c'), ボード番号 1 byte, 予約 2 byte,
#   枠の長さ 4 byte, 枠, 枠の CRC-32 4 byte
# HA の枠 (側 'h') も同じ流れで届くが、ここでは読み飛ばす (ql_psa.py が使う)。
#
# N_trg が大きく戻ったら新しい run とみなし、積算をやり直す。
#

import argparse
import asyncio
import http.server
import os
import shutil
import struct
import subprocess
import threading
import time
import traceback
import webbrowser
import zlib
from collections import OrderedDict

import numpy as np

from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

from check_tpc import (
    FADC_MAX,
    FADC_SIZE,
    FOOTER,
    HEAD_SIZE,
    MIN_FRAME_SIZE,
    N_BOARDS_PER_UPIC,
    N_CLOCKS,
    N_FADC_CH,
    N_FADC_SAMPLES,
    N_STRIPS,
    RECORD_SIZE,
    STRIP_AXES,
    add_analysis_args,
    decode_file,
    N_UPIC,
    draw_rate,
    measure_tracks,
    save_figure,
    track_lengths,
)


MAGIC = 0xeb905154
PACKET_HEAD = struct.Struct(">IBBHI")
PACKET_TAIL = struct.Struct(">I")
MAX_FRAME_SIZE = 64 * 1024 * 1024

# 枠を持っておくイベント数。1 イベント表示に使う
KEEP_FRAMES = 32

# N_trg がこれ以上戻ったら新しい run とする
RUN_RESET_TRIGGERS = 1000

RANGE_NAMES = ("strip_min", "strip_max", "clock_min", "clock_max")

# FADC は 20 ns、飛跡は 10 ns ごとなので、1 サンプルが 2 クロックにあたる
CLOCKS_PER_SAMPLE = 2
STRIPS_PER_FADC_CH = 32

# 1 イベントの図の色。軸ごとに {uPIC: 色}
EVENT_COLORS = {
    "x": {0: "red", 3: "red", 1: "blue", 2: "blue"},
    "y": {0: "red", 1: "red", 2: "blue", 3: "blue"},
}
EVENT_TITLES = {"x": "x (red: N, blue: S)", "y": "y (red: E, blue: W)"}

VIEW_PAGE = """<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<title>TPC QL</title>
<style>
html, body { margin: 0; height: 100%; overflow: hidden; background: #fff; }
body { display: grid; height: 100vh; grid-template-columns: 42fr 58fr; grid-template-rows: 26vh 74vh; }
img { width: 100%; height: 100%; min-height: 0; object-fit: contain; }
#summary { grid-column: 1 / 3; }
#status { position: fixed; right: 8px; bottom: 4px; font: 12px sans-serif; color: #555; }
</style>
</head>
<body>
<img id="event" alt="">
<img id="stats" alt="">
<img id="summary" alt="">
<div id="status">図を待っています</div>
<script>
const statusLine = document.getElementById("status");
function watch(name, interval) {
  const figure = document.getElementById(name);
  function refresh() {
    const next = new Image();
    next.onload = () => {
      figure.src = next.src;
      statusLine.textContent = "読み込み " + new Date().toLocaleTimeString();
    };
    next.src = "/" + name + ".png?t=" + Date.now();
  }
  refresh();
  setInterval(refresh, interval);
}
watch("event", EVENT_MS);
watch("stats", SUMMARY_MS);
watch("summary", SUMMARY_MS);
</script>
</body>
</html>
"""


def parse_args():

    parser = argparse.ArgumentParser(description="TPC quick look for thinned events from zeno")
    parser.add_argument("--host", default="zeno", help="host sending the thinned events")
    parser.add_argument("--port", type=int, default=10012, help="TCP port of the thinned events")
    parser.add_argument("--update", type=float, default=10.0,
                        help="update interval of the accumulated figures [s]")
    parser.add_argument("--event-update", type=float, default=1.0,
                        help="update interval of the one-event figure [s]")
    parser.add_argument("--fadc-scale", type=float, default=0.2,
                        help="width of the FADC waveform in the one-event figure [mm/ADU]")
    parser.add_argument("--reconnect", type=float, default=3.0, help="reconnect interval [s]")
    parser.add_argument("--name", default="ql_tpc", help="name of the output PNG")
    parser.add_argument("--http-bind", default="127.0.0.1", help="address of the viewer")
    parser.add_argument("--http-port", type=int, default=8090, help="TCP port of the viewer")
    parser.add_argument("--no-open", action="store_true", help="do not open a browser")
    add_analysis_args(parser)
    parser.set_defaults(rate_bin=60.0)

    return parser.parse_args()


def encode_packet(side, index, frame):

    head = PACKET_HEAD.pack(MAGIC, ord(side), index, 0, len(frame))

    return head + frame + PACKET_TAIL.pack(zlib.crc32(frame))


def thinned_rate(all_boards, time_origin, args):
    """間引いたデータから、飛跡のあるトリガのレートを推定する。

    区間ごとに トリガレート (N_trg の増え方) x 届いたイベントのうち飛跡のある割合 とする。
    """

    trig = np.concatenate([board["trig"] for board in all_boards]).astype(np.int64)
    times = np.concatenate([board["time"] for board in all_boards]).astype(np.float64)
    tracked = np.concatenate([board["n_records"] > 0 for board in all_boards])

    unique, first = np.unique(trig, return_index=True)
    moment = times[first] * args.time_tick - time_origin
    has_track = np.isin(unique, trig[tracked])

    if len(unique) < 2:
        return None

    n_bins = int(moment.max() // args.rate_bin) + 1
    edges = np.arange(n_bins + 1) * args.rate_bin
    rate = np.full(n_bins, np.nan)

    for index in range(n_bins):

        selected = (moment >= edges[index]) & (moment < edges[index + 1])

        if selected.sum() < 2:
            continue

        span = moment[selected].max() - moment[selected].min()

        if span > 0:
            trigger_rate = (unique[selected].max() - unique[selected].min()) / span
            rate[index] = trigger_rate * has_track[selected].mean()

    return edges, rate


class ViewHandler(http.server.BaseHTTPRequestHandler):
    """ブラウザに表示用のページと最新の PNG を返す。"""

    def do_GET(self):

        path = self.path.split("?")[0]

        if path == "/":
            self.reply("text/html; charset=utf-8", self.server.page)
        elif path in self.server.images:
            try:
                with open(self.server.images[path], "rb") as handle:
                    self.reply("image/png", handle.read())
            except FileNotFoundError:
                self.send_error(404)
        else:
            self.send_error(404)

    def reply(self, content_type, body):

        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *_):

        pass


def start_viewer(args, paths, page=VIEW_PAGE):
    """表示用の HTTP サーバーを別スレッドで動かし、URL を返す。"""

    server = http.server.ThreadingHTTPServer((args.http_bind, args.http_port), ViewHandler)
    server.daemon_threads = True
    server.images = {f"/{name}.png": path for name, path in paths.items()}
    server.page = (page
                   .replace("EVENT_MS", str(int(args.event_update * 1000)))
                   .replace("SUMMARY_MS", str(int(args.update * 1000)))
                   .encode("utf-8"))

    threading.Thread(target=server.serve_forever, daemon=True).start()

    host = "127.0.0.1" if args.http_bind == "0.0.0.0" else args.http_bind

    return f"http://{host}:{args.http_port}/"


def parse_event(frame):
    """1 枠から (N_time, 波形 [ch, sample], 飛跡の (strip, clock) の配列) を返す。"""

    moment = int.from_bytes(frame[8:12], "big")

    if int.from_bytes(frame[HEAD_SIZE:HEAD_SIZE + 4], "big") == FOOTER:
        return moment, None, np.zeros((0, 2), dtype=np.int64)

    fadc = np.frombuffer(frame, dtype=">u2", count=N_FADC_CH * N_FADC_SAMPLES, offset=HEAD_SIZE)
    wave = (fadc.reshape(N_FADC_SAMPLES, N_FADC_CH).T & 0x3ff).astype(np.float64)

    first_record = HEAD_SIZE + FADC_SIZE + 4
    count = (len(frame) - first_record - 4) // RECORD_SIZE
    records = np.frombuffer(frame, dtype=">u4", count=5 * count,
                            offset=first_record).reshape(count, 5)

    hits = []

    for record in records:

        clock = int(record[0]) & 0x7ff
        if clock >= N_CLOCKS:
            continue

        for ch in range(N_FADC_CH):
            strips = int(record[4 - ch])
            for bit in range(32):
                if (strips >> bit) & 1:
                    hits.append((32 * ch + bit, clock))

    return moment, wave, np.array(hits, dtype=np.int64).reshape(-1, 2)


def draw_event_figure(frames, trig, out_path, args):
    """1 イベントの飛跡を描き、飛跡のあるボードの FADC 波形を対応するストリップに重ねる。"""

    fig = Figure(figsize=(12, 3.6), constrained_layout=True)
    ax_x, ax_y = fig.subplots(1, 2)
    axes = {"x": ax_x, "y": ax_y}

    extent = N_BOARDS_PER_UPIC * N_STRIPS * args.pitch_strip * 1.1
    moment = None
    n_hit_boards = 0

    for (side, index), frame in frames.items():

        moment, wave, hits = parse_event(frame)

        if wave is None or len(hits) == 0:
            continue

        n_hit_boards += 1

        upic = index // N_BOARDS_PER_UPIC
        offset = index % N_BOARDS_PER_UPIC
        axis, sign, constant = STRIP_AXES[upic][0 if side == "a" else 1]
        color = EVENT_COLORS[axis][upic]
        ax = axes[axis]

        position = args.pitch_strip * (sign * (hits[:, 0] + offset * N_STRIPS) + constant)
        ax.scatter(position, hits[:, 1], s=4, color=color, marker="s", linewidths=0)

        clocks = np.arange(N_FADC_SAMPLES) * CLOCKS_PER_SAMPLE

        for ch in range(N_FADC_CH):

            height = wave[ch] - wave[ch, :args.base_window].mean()

            # カソードは反転して、アノードと同じ向きに出す
            if side == "c":
                height = -height

            strip = offset * N_STRIPS + STRIPS_PER_FADC_CH * ch + (STRIPS_PER_FADC_CH - 1) / 2
            center = args.pitch_strip * (sign * strip + constant)

            ax.plot(center + args.fadc_scale * height, clocks, color=color, linewidth=0.6)

    for name, ax in axes.items():
        ax.set_xlim(-extent, extent)
        ax.set_ylim(0, N_CLOCKS)
        ax.set_title(EVENT_TITLES[name])
        ax.set_xlabel(f"{name} [mm]")
        ax.set_ylabel("clock")

    fig.suptitle(f"Event ID (N_trg): {trig}   Time cnt: {moment}   Hit boards: {n_hit_boards}",
                 fontsize=10)

    save_figure(fig, out_path)


def draw_dedx(ax, boards, args):
    """全 uPIC の dE/dx (飛跡長とアノード電荷の 2 次元図) を描く。check_tpc.py の 5 段目と同じ計算。"""

    lengths = []
    charges = []

    for upic in range(N_UPIC):
        selected = {
            side: [(offset, boards[(side, upic * N_BOARDS_PER_UPIC + offset)])
                   for offset in range(N_BOARDS_PER_UPIC)
                   if (side, upic * N_BOARDS_PER_UPIC + offset) in boards]
            for side in "ac"
        }
        _, length, charge = track_lengths(selected, args)
        lengths.append(length)
        charges.append(charge)

    length = np.concatenate(lengths)
    charge = np.concatenate(charges)
    counts, x_edges, y_edges = np.histogram2d(charge, length, bins=100,
                                              range=[[0, args.dedx_q_max], [0, args.length_max]])

    if counts.max() > 0:
        ax.pcolormesh(x_edges, y_edges, np.ma.masked_equal(counts.T, 0), norm=LogNorm(vmin=1),
                      cmap="viridis")

    ax.set_xlim(0, args.dedx_q_max)
    ax.set_ylim(0, args.length_max)
    ax.set_title(f"dE/dx (all uPIC, {len(length)} tracks)")
    ax.set_xlabel("anode charge [ADU x sample]")
    ax.set_ylabel("track length [mm]")


def draw_baselines(ax_cathode, ax_anode, boards):
    """FADC の ch ごとのベースラインの中央値と標準偏差を、cathode (上) と anode (下) に分けて描く。

    ベースラインは 1-2 ADC の幅に集まるので、ADC ごとのヒストグラムにすると線が 1 画素より細くなり、
    多くの ch が描かれずに抜けて見える。そのため ch ごとの代表値で描く。値の無い ch は灰色の x にする。
    """

    n_ch = N_UPIC * N_BOARDS_PER_UPIC * N_FADC_CH
    adc = np.arange(FADC_MAX + 1)

    for ax, side, name in ((ax_cathode, "c", "cathode"), (ax_anode, "a", "anode")):

        medians = np.full(n_ch, np.nan)
        spreads = np.full(n_ch, np.nan)

        for (board_side, index), board in boards.items():
            if board_side != side or index >= N_UPIC * N_BOARDS_PER_UPIC:
                continue
            for ch in range(N_FADC_CH):
                counts = board["baseline"][ch].astype(np.float64)
                if counts.sum() == 0:
                    continue
                cumulative = np.cumsum(counts)
                medians[index * N_FADC_CH + ch] = adc[np.searchsorted(cumulative, cumulative[-1] / 2)]
                mean = (adc * counts).sum() / counts.sum()
                spreads[index * N_FADC_CH + ch] = np.sqrt((counts * (adc - mean) ** 2).sum() / counts.sum())

        positions = np.arange(n_ch)
        found = ~np.isnan(medians)
        ax.errorbar(positions[found], medians[found], yerr=spreads[found], fmt="o", markersize=2.5,
                    capsize=1.5, linewidth=0.8, color="tab:orange" if side == "c" else "tab:blue")

        if found.any():
            low, high = np.nanmin(medians - spreads), np.nanmax(medians + spreads)
            margin = max(5.0, 0.15 * (high - low))
            ax.set_ylim(low - margin, high + margin)
        ax.plot(positions[~found], np.full((~found).sum(), ax.get_ylim()[0]), "x", color="0.6",
                markersize=4, clip_on=False)

        # uPIC の境目に線を引き、目盛りは uPIC ごとに 1 つにする
        per_upic = N_BOARDS_PER_UPIC * N_FADC_CH
        for upic in range(1, N_UPIC):
            ax.axvline(upic * per_upic - 0.5, color="0.8", linewidth=0.8)
        ax.set_xlim(-0.5, n_ch - 0.5)
        ax.set_xticks([upic * per_upic + (per_upic - 1) / 2 for upic in range(N_UPIC)],
                      [f"uPIC{upic}" for upic in range(N_UPIC)])
        ax.set_ylabel(f"{name} [ADC]")

    ax_cathode.set_title(f"FADC baseline (median, std; x: no data, {N_FADC_CH} ch x 3 boards per uPIC)",
                         fontsize=10)
    ax_cathode.tick_params(labelbottom=False)


def draw_stats_figure(boards, title, out_path, args):
    """カウントレート、FADC ベースライン、dE/dx を描く。"""

    fig = Figure(figsize=(16, 3.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 3)
    ax_rate = fig.add_subplot(grid[:, 0])
    ax_cathode = fig.add_subplot(grid[0, 1])
    ax_anode = fig.add_subplot(grid[1, 1], sharex=ax_cathode)
    ax_dedx = fig.add_subplot(grid[:, 2])

    draw_rate(ax_rate, boards, args, rate_curve=thinned_rate)
    draw_baselines(ax_cathode, ax_anode, boards)
    draw_dedx(ax_dedx, boards, args)
    fig.suptitle(title, fontsize=10)

    save_figure(fig, out_path)


def draw_summary_figure(boards, out_path, args):
    """飛跡の積算画像を描く。上の段が uPIC 0-3 の anode、下の段が cathode。"""

    # 画面の下段いっぱいに大きく出すので、軸の名前は外側の図だけに付けて余白を減らす
    fig = Figure(figsize=(20, 8.3), constrained_layout=True)
    axes = fig.subplots(2, N_UPIC, sharex=True, sharey=True)

    for row, side in enumerate("ac"):
        for upic in range(N_UPIC):

            ax = axes[row, upic]
            image = np.zeros((N_CLOCKS, N_BOARDS_PER_UPIC * N_STRIPS), dtype=np.float64)

            # 1 つの uPIC の 1 つの側は、ボード 3 枚 (128 ストリップずつ) を横に並べる
            for offset in range(N_BOARDS_PER_UPIC):
                key = (side, upic * N_BOARDS_PER_UPIC + offset)
                if key in boards:
                    image[:, offset * N_STRIPS:(offset + 1) * N_STRIPS] = boards[key]["image"]

            if image.max() > 0:
                ax.imshow(np.ma.masked_equal(image, 0), origin="lower", aspect="auto",
                          norm=LogNorm(vmin=1), cmap="viridis", interpolation="nearest")

            # データの無い uPIC も、ほかと同じ目盛りにする
            ax.set_xlim(-0.5, N_BOARDS_PER_UPIC * N_STRIPS - 0.5)
            ax.set_ylim(-0.5, N_CLOCKS - 0.5)
            ax.set_title(f"uPIC{upic} {'anode' if side == 'a' else 'cathode'}")
            if row == 1:
                ax.set_xlabel("strip")
            if upic == 0:
                ax.set_ylabel("clock")

    save_figure(fig, out_path)


def open_browser(url):

    # WSL では Windows 側のブラウザで開く
    if shutil.which("explorer.exe"):
        subprocess.Popen(["explorer.exe", url])
    else:
        webbrowser.open(url)


def prepare_view(args, paths, page=VIEW_PAGE):
    """前回の図を消し、表示用のサーバーを動かしてブラウザを開く。"""

    os.makedirs(args.out_dir, exist_ok=True)

    # 前回の図を今の図と見間違えないよう、最初の更新までは表示しない
    for path in paths.values():
        if os.path.exists(path):
            os.remove(path)

    url = start_viewer(args, paths, page)
    print(f"[VIEW] {url}", flush=True)

    if not args.no_open:
        open_browser(url)


async def receive_packets(reader, sides, on_frame, on_bad):
    """パケットを読み、sides の側の枠を on_frame(側, 番号, 枠) に渡す。CRC が合わなければ on_bad()。"""

    while True:

        head = await reader.readexactly(PACKET_HEAD.size)
        magic, side, index, _, size = PACKET_HEAD.unpack(head)

        if magic != MAGIC or size > MAX_FRAME_SIZE:
            raise ConnectionError(f"wrong packet header: {head.hex()}")

        frame = await reader.readexactly(size)
        crc, = PACKET_TAIL.unpack(await reader.readexactly(PACKET_TAIL.size))

        # 別の検出器の枠は読み飛ばす
        if chr(side) not in sides:
            continue

        if zlib.crc32(frame) != crc:
            on_bad()
            continue

        on_frame(chr(side), index, frame)


async def receive_loop(args, sides, on_frame, on_bad):
    """送信元へつなぎ、切れたら args.reconnect 秒ごとにつなぎ直す。"""

    while True:

        writer = None

        try:
            reader, writer = await asyncio.open_connection(args.host, args.port)
            print(f"[CONNECT] {args.host}:{args.port}", flush=True)
            await receive_packets(reader, sides, on_frame, on_bad)
        except (OSError, asyncio.IncompleteReadError, ConnectionError) as exc:
            print(f"[RETRY] {args.host}:{args.port}: {exc}", flush=True)
        finally:
            if writer is not None:
                writer.close()

        await asyncio.sleep(args.reconnect)


class BoardState:
    """1 ボードの積算。check_tpc の描画関数に渡す形へ変換できる。"""

    def __init__(self, side):

        self.side = side
        self.image = np.zeros((N_CLOCKS, N_STRIPS), dtype=np.uint32)
        self.baseline = np.zeros((N_FADC_CH, FADC_MAX + 1), dtype=np.uint32)
        self.wave_map = np.zeros((FADC_MAX + 1, N_FADC_SAMPLES), dtype=np.uint32)
        self.rows = {name: [] for name in ("trig", "time", "n_records") + RANGE_NAMES}
        self.charge = []
        self.frames = OrderedDict()

    def add(self, frame, args):
        """1 枠を積算する。読めた枠の N_trg を返す。読めなければ None。"""

        buf = np.frombuffer(frame, dtype=np.uint8)

        trig = np.zeros(1, dtype=np.uint32)
        moment = np.zeros(1, dtype=np.uint32)
        charge = np.zeros((1, N_FADC_CH), dtype=np.float32)
        n_records = np.zeros(1, dtype=np.int32)
        offsets = np.zeros(1, dtype=np.uint32)

        if len(buf) < MIN_FRAME_SIZE:
            return None

        n_events, _ = decode_file(buf, self.image, self.baseline, self.wave_map, trig, moment,
                                  charge, n_records, offsets, args.base_window, args.threshold,
                                  self.side == "c")

        if n_events != 1:
            return None

        ranges = {name: np.full(1, -1, dtype=np.int32) for name in RANGE_NAMES}

        measure_tracks(buf, offsets, n_records, ranges["strip_min"], ranges["strip_max"],
                       ranges["clock_min"], ranges["clock_max"])

        self.rows["trig"].append(int(trig[0]))
        self.rows["time"].append(int(moment[0]))
        self.rows["n_records"].append(int(n_records[0]))
        self.charge.append(charge[0].copy())

        for name in RANGE_NAMES:
            self.rows[name].append(int(ranges[name][0]))

        self.frames[int(trig[0])] = bytes(frame[int(offsets[0]):])

        while len(self.frames) > KEEP_FRAMES:
            self.frames.popitem(last=False)

        return int(trig[0])

    def snapshot(self):
        """描画中も受信を続けられるよう、今の積算の写しを返す。"""

        board = {
            "trig": np.array(self.rows["trig"], dtype=np.uint32),
            "time": np.array(self.rows["time"], dtype=np.uint32),
            "n_records": np.array(self.rows["n_records"], dtype=np.int32),
            "charge": np.array(self.charge, dtype=np.float32).reshape(-1, N_FADC_CH),
            "offsets": np.zeros(len(self.rows["trig"]), dtype=np.uint32),
            "image": self.image.copy(),
            "baseline": self.baseline.copy(),
            "wave_map": self.wave_map.copy(),
            "frames": dict(self.frames),
        }

        for name in RANGE_NAMES:
            board[name] = np.array(self.rows[name], dtype=np.int32)

        return board


class QuickLook:

    def __init__(self, args):

        self.args = args
        self.paths = {name: os.path.join(args.out_dir, f"{args.name}_{name}.png")
                      for name in ("event", "stats", "summary")}
        self.bad_packets = 0
        self.reset()

    def reset(self):

        self.boards = {}
        self.latest = -1
        self.latest_tracked = None
        self.n_packets = 0

    def add(self, side, index, frame):

        board = self.boards.get((side, index))

        if board is None:
            board = BoardState(side)
            self.boards[(side, index)] = board

        trig = board.add(frame, self.args)

        if trig is None:
            self.bad_packets += 1
            return

        if trig + RUN_RESET_TRIGGERS < self.latest:
            print(f"[RUN] N_trg went back ({self.latest} -> {trig}). restart the accumulation", flush=True)
            self.reset()
            self.add(side, index, frame)
            return

        self.latest = max(self.latest, trig)
        self.n_packets += 1

        if board.rows["n_records"][-1] > 0 and trig in board.frames:
            self.latest_tracked = trig

    def snapshot(self):
        """描画に渡す積算の写しと、図の題を返す。"""

        boards = {key: board.snapshot() for key, board in self.boards.items()}
        n_events = len(set(trig for board in boards.values() for trig in board["trig"].tolist()))

        now = time.strftime("%H:%M:%S")
        title = (f"QL {now}: {len(boards)} boards, {n_events} thinned events, "
                 f"latest N_trg {self.latest}, broken packets {self.bad_packets}")

        return boards, title

    def event_frames(self):
        """最新の飛跡イベントの枠を、ボードごとに返す。"""

        trig = self.latest_tracked

        return trig, {key: board.frames[trig] for key, board in self.boards.items()
                      if trig in board.frames}

    async def draw_event(self):

        trig, frames = self.event_frames()
        await asyncio.to_thread(draw_event_figure, frames, trig, self.paths["event"], self.args)

    async def draw_accumulated(self):

        boards, title = self.snapshot()

        await asyncio.to_thread(draw_stats_figure, boards, title, self.paths["stats"], self.args)
        await asyncio.to_thread(draw_summary_figure, boards, self.paths["summary"], self.args)

        print(f"[DRAW] {title}", flush=True)

    async def update_loop(self):
        """図を 1 枚ずつ順に描く。matplotlib を同時に使わないようにするため。"""

        drawn_event = None
        next_accumulated = time.monotonic() + self.args.update

        while True:

            await asyncio.sleep(self.args.event_update)

            # 描画に失敗しても受信は止めない
            try:
                if self.latest_tracked is not None and self.latest_tracked != drawn_event:
                    drawn_event = self.latest_tracked
                    await self.draw_event()

                if len(self.boards) > 0 and time.monotonic() >= next_accumulated:
                    next_accumulated = time.monotonic() + self.args.update
                    await self.draw_accumulated()
            except Exception:
                traceback.print_exc()

    def count_bad(self):

        self.bad_packets += 1


async def run(args):

    quick_look = QuickLook(args)

    prepare_view(args, quick_look.paths)

    await asyncio.gather(receive_loop(args, "ac", quick_look.add, quick_look.count_bad),
                         quick_look.update_loop())


def main():

    args = parse_args()

    try:
        asyncio.run(run(args))
    except KeyboardInterrupt:
        print("[STOP]")


if __name__ == "__main__":

    main()
