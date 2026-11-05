#!/usr/bin/env python3
#
# SMILE-ICS が保存した TPC エンコーダの raw (tpc_enc_a00_00000.raw など) を読み、
# データが正しそうかを見るための確認図を 1 枚の PNG に描く。
#
#   1 段目 : uPIC ごとのカウントレート、FADC ベースライン (全 ch)
#   2 段目 : uPIC ごとの飛跡の積算画像 (anode, cathode)
#   3 段目 : uPIC ごとの FADC 波形の積算図 (anode, cathode)。カソードも生の値
#   4 段目 : uPIC ごとのスペクトル (anode, cathode)。エネルギーには換算していない
#   5 段目 : 1 イベントの飛跡、uPIC ごとの dE/dx (飛跡長と電荷の 2 次元図)
#
# 電荷の計算は EricaNext の CalcFADC、座標の向きは TOT_Skewness と同じ。
# イベントの選別はしない。飛跡長は TOT 補正と 3 次元コインシデンスをせず、
# ストリップとクロックの広がりから出す。
# 画面には出さず PNG に保存するだけなので、XLaunch は不要。
#

import argparse
import glob
import os
import re
from concurrent.futures import ThreadPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from natsort import natsorted
from numba import njit

# このリポジトリ (analysis_sm3_FM) の場所。既定のパスはここからの相対で決める
REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


DEFAULT_OUTPUT_DIR = os.path.join(REPO_DIR, "check_output")

HEADER = 0xeb901964
FOOTER = 0x75504943

N_UPIC = 4
N_BOARDS_PER_UPIC = 3
N_FADC_CH = 4          # 1 ボードの FADC チャンネル数
N_FADC_SAMPLES = 512   # 1 チャンネルのサンプル数 (1 サンプル 20 ns)
N_STRIPS = 128         # 1 ボードのストリップ数
N_CLOCKS = 1024        # 飛跡の時間方向のクロック数
FADC_MAX = 1023        # 10 bit

HEAD_SIZE = 16         # ヘッダ + N_trg + N_time + N_gso [byte]
FADC_SIZE = 8 * N_FADC_SAMPLES
RECORD_SIZE = 20       # clk + ストリップ 32 bit x 4 [byte]
MIN_FRAME_SIZE = 20    # ヒットのない枠 (ヘッダ + 3 語 + フッタ) [byte]

UPIC_COLORS = ("tab:blue", "tab:orange", "tab:green", "tab:red")

# uPIC ごとの座標 = pitch x (符号 x strip + 定数)。TOT_Skewness::Get_TPC_XYZ と同じ
# {uPIC: (anode の軸, anode の符号, 定数), (cathode の軸, cathode の符号, 定数)}
STRIP_AXES = {
    0: (("x", 1, 0.5), ("y", 1, -383.5)),
    1: (("y", 1, 0.5), ("x", -1, 383.5)),
    2: (("x", -1, -0.5), ("y", -1, 383.5)),
    3: (("y", -1, -0.5), ("x", 1, -383.5)),
}

FILE_PATTERN = re.compile(r"tpc_enc_([ac])(\d+)_\d+\.raw$")


def parse_args():

    parser = argparse.ArgumentParser(description="Draw sanity-check plots for TPC encoder raw files")
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="run directories (tpc_enc_*.raw). give all nodes of a 4-node run")
    parser.add_argument("--tag", default=None, help="name of the output PNG (default: input directory name)")
    parser.add_argument("--max-files", type=int, default=None, help="files to read per board (default: all)")
    parser.add_argument("--workers", type=int, default=4, help="boards decoded at the same time")
    add_analysis_args(parser)

    return parser.parse_args()


def add_analysis_args(parser):
    """check_tpc.py と ql_tpc.py で共通の引数を足す。"""

    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-window", type=int, default=30,
                        help="samples for the baseline mean [samples, 1 sample = 20 ns]")
    parser.add_argument("--threshold", type=float, default=20.0,
                        help="height above the baseline to integrate [ADU]")
    parser.add_argument("--rate-bin", type=float, default=1.0, help="count rate bin [s]")
    parser.add_argument("--time-tick", type=float, default=100e-6,
                        help="N_time tick [s]. 100 us was measured from run_20260623T072649Z")
    parser.add_argument("--q-max", type=float, default=20000.0,
                        help="upper edge of the spectrum [ADU x sample]")
    parser.add_argument("--q-bins", type=int, default=400, help="bins of the spectrum")
    parser.add_argument("--event", type=int, default=None,
                        help="N_trg of the event to draw (default: the first event with a track)")
    parser.add_argument("--pitch-strip", type=float, default=0.4, help="strip pitch [mm/strip]")
    parser.add_argument("--pitch-clock", type=float, default=0.5, help="drift length [mm/clock]")
    parser.add_argument("--length-max", type=float, default=300.0, help="dE/dx y range [mm]")
    parser.add_argument("--dedx-q-max", type=float, default=5000.0,
                        help="dE/dx x range [ADU x sample]")


@njit(cache=True, nogil=True)
def read_word(buf, pos):

    return ((int(buf[pos]) << 24) | (int(buf[pos + 1]) << 16)
            | (int(buf[pos + 2]) << 8) | int(buf[pos + 3]))


@njit(cache=True, nogil=True)
def decode_file(buf, image, baseline, wave_map, trig, time, charge, n_records, offsets,
                base_window, threshold, invert):
    """1 ファイルを読む。image, baseline, wave_map に積算し、読んだ枠の数と壊れた枠の数を返す。"""

    n_bytes = len(buf)
    wave = np.zeros((N_FADC_CH, N_FADC_SAMPLES), dtype=np.int32)

    pos = 0
    n_events = 0
    n_bad = 0

    while pos + MIN_FRAME_SIZE <= n_bytes:

        if read_word(buf, pos) != HEADER:
            pos += 1
            continue

        event_trig = read_word(buf, pos + 4)
        event_time = read_word(buf, pos + 8)
        offsets[n_events] = pos
        body = pos + HEAD_SIZE

        # ヒットのない枠は、3 語の直後にフッタが来る
        if read_word(buf, body) == FOOTER:
            trig[n_events] = event_trig
            time[n_events] = event_time
            n_events += 1
            pos = body + 4
            continue

        # FADC の後に 1 語あり、その後に飛跡の記録がフッタまで続く
        first_record = body + FADC_SIZE + 4
        record = first_record
        count = 0
        found = False
        while count < N_CLOCKS and record + 4 <= n_bytes:
            if read_word(buf, record) == FOOTER:
                found = True
                break
            if record + RECORD_SIZE > n_bytes:
                break
            record += RECORD_SIZE
            count += 1

        if not found and count < N_CLOCKS:
            # ファイルの途中で切れている、または形式が合わない
            n_bad += 1
            pos += 4
            continue

        trig[n_events] = event_trig
        time[n_events] = event_time
        n_records[n_events] = count

        if count > 0:

            for sample in range(N_FADC_SAMPLES):
                word0 = read_word(buf, body + 8 * sample)
                word1 = read_word(buf, body + 8 * sample + 4)
                wave[0, sample] = (word0 >> 16) & 0x3ff
                wave[1, sample] = word0 & 0x3ff
                wave[2, sample] = (word1 >> 16) & 0x3ff
                wave[3, sample] = word1 & 0x3ff
                for ch in range(N_FADC_CH):
                    wave_map[wave[ch, sample], sample] += 1

            for ch in range(N_FADC_CH):

                mean = 0.0
                for sample in range(base_window):
                    mean += wave[ch, sample]
                mean /= base_window

                # ベースラインの図はカソードも生の値のまま描く
                baseline[ch, min(int(mean), FADC_MAX)] += 1

                total = 0.0
                for sample in range(base_window, N_FADC_SAMPLES):
                    if invert:
                        height = mean - wave[ch, sample]
                    else:
                        height = wave[ch, sample] - mean
                    if height > threshold:
                        total += height
                charge[n_events, ch] = total

            for index in range(count):
                start = first_record + RECORD_SIZE * index
                clock = read_word(buf, start) & 0x7ff
                if clock >= N_CLOCKS:
                    continue
                # 語の並びは ch3, ch2, ch1, ch0 の順、語の中は LSB が小さいストリップ（Erica の reverse_TPC: false と同じ）
                for ch in range(N_FADC_CH):
                    strips = read_word(buf, start + 4 * (4 - ch))
                    for bit in range(32):
                        if (strips >> bit) & 1:
                            image[clock, 32 * ch + bit] += 1

        n_events += 1
        pos = record + 4 if found else record

    return n_events, n_bad


@njit(cache=True, nogil=True)
def measure_tracks(buf, offsets, n_records, strip_min, strip_max, clock_min, clock_max):
    """イベントごとのストリップとクロックの範囲を求める。"""

    for event in range(len(offsets)):

        first_record = int(offsets[event]) + HEAD_SIZE + FADC_SIZE + 4

        for index in range(n_records[event]):

            start = first_record + RECORD_SIZE * index
            clock = read_word(buf, start) & 0x7ff
            if clock >= N_CLOCKS:
                continue

            for ch in range(N_FADC_CH):
                strips = read_word(buf, start + 4 * (4 - ch))
                for bit in range(32):
                    strip = 32 * ch + bit
                    if (strips >> bit) & 1:
                        if strip_min[event] < 0 or strip < strip_min[event]:
                            strip_min[event] = strip
                        strip_max[event] = max(strip_max[event], strip)
                        if clock_min[event] < 0 or clock < clock_min[event]:
                            clock_min[event] = clock
                        clock_max[event] = max(clock_max[event], clock)


def collect_boards(directories, max_files):
    """1 つまたは複数の run ディレクトリから、ボードごとの raw のリストを返す。

    4 ノード構成では 1 つの run がノードごとのディレクトリに分かれるので、まとめて渡す。
    """

    if isinstance(directories, str):
        directories = [directories]

    boards = {}

    for path in (path for directory in directories
                 for path in glob.glob(os.path.join(directory, "tpc_enc_*.raw"))):

        match = FILE_PATTERN.search(os.path.basename(path))

        if match:
            boards.setdefault((match.group(1), int(match.group(2))), []).append(path)

    return {key: natsorted(paths)[:max_files] for key, paths in sorted(boards.items())}


def load_board(key, paths, args):

    invert = key[0] == "c"

    image = np.zeros((N_CLOCKS, N_STRIPS), dtype=np.uint32)
    baseline = np.zeros((N_FADC_CH, FADC_MAX + 1), dtype=np.uint32)
    wave_map = np.zeros((FADC_MAX + 1, N_FADC_SAMPLES), dtype=np.uint32)

    parts = {"trig": [], "time": [], "charge": [], "n_records": [], "offsets": [], "file_index": []}
    n_bad = 0

    for file_index, path in enumerate(paths):

        buf = np.fromfile(path, dtype=np.uint8)
        capacity = len(buf) // MIN_FRAME_SIZE + 1

        trig = np.zeros(capacity, dtype=np.uint32)
        time = np.zeros(capacity, dtype=np.uint32)
        charge = np.zeros((capacity, N_FADC_CH), dtype=np.float32)
        n_records = np.zeros(capacity, dtype=np.int32)
        offsets = np.zeros(capacity, dtype=np.uint32)

        n_events, bad = decode_file(buf, image, baseline, wave_map, trig, time, charge, n_records, offsets,
                                    args.base_window, args.threshold, invert)

        n_bad += bad
        parts["trig"].append(trig[:n_events].copy())
        parts["time"].append(time[:n_events].copy())
        parts["charge"].append(charge[:n_events].copy())
        parts["n_records"].append(n_records[:n_events].copy())
        parts["offsets"].append(offsets[:n_events].copy())
        parts["file_index"].append(np.full(n_events, file_index, dtype=np.int32))

    data = {name: np.concatenate(value) for name, value in parts.items()}
    data.update(image=image, baseline=baseline, wave_map=wave_map, n_bad=n_bad, paths=paths)

    return data


def measure_board(board):
    """イベントごとの飛跡の範囲を board に追加する。"""

    n_events = len(board["trig"])
    ranges = {name: np.full(n_events, -1, dtype=np.int32)
              for name in ("strip_min", "strip_max", "clock_min", "clock_max")}

    for file_index, path in enumerate(board["paths"]):

        selected = np.flatnonzero(board["file_index"] == file_index)
        part = {name: value[selected].copy() for name, value in ranges.items()}

        measure_tracks(np.fromfile(path, dtype=np.uint8), board["offsets"][selected],
                       board["n_records"][selected],
                       part["strip_min"], part["strip_max"],
                       part["clock_min"], part["clock_max"])

        for name in ranges:
            ranges[name][selected] = part[name]

    board.update(ranges)


def side_ranges(boards):
    """uPIC の片側のボードをまとめ、トリガごとのストリップとクロックの範囲を返す。"""

    parts = []

    for offset, board in boards:

        valid = board["strip_min"] >= 0
        parts.append(np.column_stack([
            board["trig"][valid],
            board["strip_min"][valid] + offset * N_STRIPS,
            board["strip_max"][valid] + offset * N_STRIPS,
            board["clock_min"][valid],
            board["clock_max"][valid],
        ]))

    table = np.concatenate(parts).astype(np.int64)
    trig, inverse = np.unique(table[:, 0], return_inverse=True)

    low = np.full((len(trig), 2), np.iinfo(np.int64).max, dtype=np.int64)
    high = np.full((len(trig), 2), -1, dtype=np.int64)

    np.minimum.at(low, inverse, table[:, [1, 3]])
    np.maximum.at(high, inverse, table[:, [2, 4]])

    return trig, low, high


def track_lengths(selected, args):
    """飛跡のあるトリガの (trig, 長さ [mm], アノード電荷) を返す。片側だけの飛跡も含める。"""

    tables = {side: side_ranges(selected[side]) for side in "ac" if len(selected[side]) > 0}

    if len(tables) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0), np.zeros(0)

    trig = np.unique(np.concatenate([table[0] for table in tables.values()]))

    widths = np.zeros((len(trig), 2))
    clock_low = np.full(len(trig), np.iinfo(np.int64).max, dtype=np.int64)
    clock_high = np.full(len(trig), -1, dtype=np.int64)

    for column, (side_trig, low, high) in enumerate(tables.values()):

        index = np.searchsorted(trig, side_trig)
        widths[index, column] = (high[:, 0] - low[:, 0]) * args.pitch_strip
        clock_low[index] = np.minimum(clock_low[index], low[:, 1])
        clock_high[index] = np.maximum(clock_high[index], high[:, 1])

    depth = (clock_high - clock_low) * args.pitch_clock
    length = np.sqrt(widths[:, 0] ** 2 + widths[:, 1] ** 2 + depth ** 2)

    charge = np.zeros(len(trig))

    if len(selected["a"]) > 0:
        trig_q, total = trig_charge([board for _, board in selected["a"]])
        index = np.searchsorted(trig_q, trig)
        found = index < len(trig_q)
        found[found] = trig_q[index[found]] == trig[found]
        charge[found] = total[index[found]]

    return trig, length, charge


def read_event(board, trig):
    """1 枠を読み直し、飛跡の (strip, clock) の配列を返す。"""

    where = np.flatnonzero(board["trig"] == trig)

    if len(where) == 0 or board["n_records"][where[0]] == 0:
        return None

    event = where[0]
    offset = int(board["offsets"][event])
    count = int(board["n_records"][event])
    size = HEAD_SIZE + FADC_SIZE + 4 + RECORD_SIZE * count

    if "frames" in board:
        # ql_tpc.py は最近の枠をメモリに持っている
        words = board["frames"].get(trig)
        if words is None:
            return None
    else:
        words = np.fromfile(board["paths"][board["file_index"][event]], dtype=np.uint8,
                            count=size, offset=offset).tobytes()

    records = np.frombuffer(words, dtype=">u4", count=5 * count,
                            offset=HEAD_SIZE + FADC_SIZE + 4).reshape(count, 5)

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

    return np.array(hits, dtype=np.int64).reshape(-1, 2)


def draw_event(fig, grid, boards, trig, args):

    axes = {"x": fig.add_subplot(grid[4, :2]), "y": fig.add_subplot(grid[4, 2:4])}
    extent = N_BOARDS_PER_UPIC * N_STRIPS * args.pitch_strip * 1.1

    for (side, index), board in boards.items():

        hits = read_event(board, trig)

        if hits is None or len(hits) == 0:
            continue

        upic = index // N_BOARDS_PER_UPIC
        offset = index % N_BOARDS_PER_UPIC
        axis, sign, constant = STRIP_AXES[upic][0 if side == "a" else 1]

        position = args.pitch_strip * (sign * (hits[:, 0] + offset * N_STRIPS) + constant)
        axes[axis].scatter(position, hits[:, 1], s=4, color=UPIC_COLORS[upic],
                           marker="s" if side == "a" else "o", linewidths=0)

    for name, ax in axes.items():
        for upic in range(N_UPIC):
            ax.scatter([], [], s=12, color=UPIC_COLORS[upic], label=f"uPIC{upic}")
        ax.legend(fontsize=8)
        ax.set_xlim(-extent, extent)
        ax.set_ylim(0, N_CLOCKS)
        ax.set_title(f"N_trg {trig}: {name}")
        ax.set_xlabel(f"{name} [mm]")
        ax.set_ylabel("clock")


def trig_charge(boards):
    """同じトリガの電荷をボードをまたいで足し、(trig, 電荷) を返す。"""

    trig = np.concatenate([board["trig"] for board in boards])
    charge = np.concatenate([board["charge"].sum(axis=1) for board in boards])

    unique, inverse = np.unique(trig, return_inverse=True)

    return unique, np.bincount(inverse, weights=charge)


def tracked_times(boards, time_tick):
    """飛跡のあるトリガの時刻 [s] を、トリガごとに 1 つ返す。"""

    trig = np.concatenate([board["trig"][board["n_records"] > 0] for board in boards])
    time = np.concatenate([board["time"][board["n_records"] > 0] for board in boards])

    _, first = np.unique(trig, return_index=True)

    return time[first].astype(np.float64) * time_tick


def saved_rate(all_boards, time_origin, args):
    """全イベントがあるとき: 飛跡のあるトリガを時間で数える。(edges, rate) か None を返す。"""

    times = tracked_times(all_boards, args.time_tick) - time_origin

    if len(times) == 0:
        return None

    # 最後の区間は途中で切れているので描かない
    n_bins = int(times.max() // args.rate_bin)

    if n_bins == 0:
        return None

    edges = np.arange(n_bins + 1) * args.rate_bin
    counts, _ = np.histogram(times, bins=edges)

    return edges, counts / args.rate_bin


def main():

    args = parse_args()

    tag = args.tag or os.path.basename(os.path.normpath(args.input[0]))

    files = collect_boards(args.input, args.max_files)

    if len(files) == 0:
        raise FileNotFoundError(f"No tpc_enc_*.raw files found: {args.input}")

    print(f"[READ] {len(files)} boards, {sum(len(paths) for paths in files.values())} files")

    with ThreadPoolExecutor(args.workers) as pool:
        results = list(pool.map(lambda item: load_board(item[0], item[1], args), files.items()))

    boards = dict(zip(files.keys(), results))

    for (side, index), board in boards.items():

        trig = board["trig"]
        n_events = len(trig)
        n_tracked = int((board["n_records"] > 0).sum())
        n_missing = int((np.diff(trig.astype(np.int64)) - 1).clip(0).sum())

        print(f"[BOARD {side}{index:02d}] events {n_events}, tracked {n_tracked}, "
              f"missing trig {n_missing}, broken frames {board['n_bad']}")

    with ThreadPoolExecutor(args.workers) as pool:
        list(pool.map(measure_board, boards.values()))

    n_events = max(len(board["trig"]) for board in boards.values())

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f"{tag}_tpc_check.png")

    draw_figure(boards, f"{tag}: {len(files)} boards, {n_events} events", out_path, args)

    print(f"[SAVE] {out_path}")


def select_upics(boards):
    """uPIC ごとに {'a': [ボード], 'c': [ボード]} をまとめる。"""

    return [
        {
            side: [boards[(side, upic * N_BOARDS_PER_UPIC + offset)]
                   for offset in range(N_BOARDS_PER_UPIC)
                   if (side, upic * N_BOARDS_PER_UPIC + offset) in boards]
            for side in "ac"
        }
        for upic in range(N_UPIC)
    ]


def draw_rate(ax, boards, args, rate_curve=saved_rate, log=lambda *_: None):

    time_origin = min(
        board["time"][0] for board in boards.values() if len(board["time"]) > 0
    ) * args.time_tick

    for upic, selected in enumerate(select_upics(boards)):

        all_boards = selected["a"] + selected["c"]

        if len(all_boards) == 0:
            continue

        curve = rate_curve(all_boards, time_origin, args)

        if curve is None:
            continue

        edges, rate = curve

        ax.step(edges[:-1], rate, where="post", label=f"uPIC{upic}")
        log(f"[UPIC {upic}] mean rate {np.nanmean(rate):.1f} /s")

    ax.set_title("count rate (triggers with a track)")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("count / s")
    ax.legend()


def draw_baseline(fig, ax, boards):

    n_boards = N_UPIC * N_BOARDS_PER_UPIC
    base_map = np.zeros((FADC_MAX + 1, 2 * n_boards * N_FADC_CH), dtype=np.float64)

    for (side, index), board in boards.items():

        if index >= n_boards:
            continue

        start = (index + (n_boards if side == "c" else 0)) * N_FADC_CH
        base_map[:, start:start + N_FADC_CH] = board["baseline"].T

    mesh = ax.pcolormesh(np.arange(base_map.shape[1] + 1) - 0.5,
                         np.arange(FADC_MAX + 2) - 0.5,
                         np.ma.masked_equal(base_map, 0),
                         norm=LogNorm(vmin=1), cmap="viridis")
    fig.colorbar(mesh, ax=ax)
    ax.axvline(n_boards * N_FADC_CH - 0.5, color="gray", linewidth=0.8)
    ax.set_title("FADC baseline (left: anode, right: cathode)")
    ax.set_xlabel("FADC ch")
    ax.set_ylabel("baseline [ADC]")


def draw_images(fig, grid, row, boards):
    """飛跡の積算画像を row の段に描く。"""

    for upic in range(N_UPIC):

        for column, side in enumerate("ac"):

            ax = fig.add_subplot(grid[row, 2 * upic + column])

            image = np.zeros((N_CLOCKS, N_BOARDS_PER_UPIC * N_STRIPS), dtype=np.float64)

            for offset in range(N_BOARDS_PER_UPIC):
                key = (side, upic * N_BOARDS_PER_UPIC + offset)
                if key in boards:
                    image[:, offset * N_STRIPS:(offset + 1) * N_STRIPS] = boards[key]["image"]

            if image.max() > 0:
                ax.imshow(np.ma.masked_equal(image, 0), origin="lower", aspect="auto",
                          norm=LogNorm(vmin=1), cmap="viridis", interpolation="nearest")

            name = "anode" if side == "a" else "cathode"
            ax.set_title(f"uPIC{upic} {name}")
            ax.set_xlabel("strip")
            ax.set_ylabel("clock")


def draw_wave_maps(fig, grid, row, boards):
    """FADC 波形の積算図を row の段に描く。"""

    for upic in range(N_UPIC):

        for column, side in enumerate("ac"):

            ax = fig.add_subplot(grid[row, 2 * upic + column])

            wave_map = np.zeros((FADC_MAX + 1, N_FADC_SAMPLES), dtype=np.float64)

            for offset in range(N_BOARDS_PER_UPIC):
                key = (side, upic * N_BOARDS_PER_UPIC + offset)
                if key in boards:
                    wave_map += boards[key]["wave_map"]

            if wave_map.max() > 0:
                ax.imshow(np.ma.masked_equal(wave_map, 0), origin="lower", aspect="auto",
                          norm=LogNorm(vmin=1), cmap="viridis", interpolation="nearest")

            name = "anode" if side == "a" else "cathode"
            ax.set_title(f"uPIC{upic} {name} FADC")
            ax.set_xlabel("sample (20 ns)")
            ax.set_ylabel("FADC [ADC]")


def draw_spectra(fig, grid, row, boards, args):
    """uPIC ごとのスペクトルを row の段に描く。"""

    for upic, selected in enumerate(select_upics(boards)):

        ax = fig.add_subplot(grid[row, 2 * upic:2 * upic + 2])

        for side, name in (("a", "anode"), ("c", "cathode")):

            if len(selected[side]) == 0:
                continue

            _, total = trig_charge(selected[side])
            total = total[total > 0]
            n_over = int((total > args.q_max).sum())

            ax.hist(total, bins=args.q_bins, range=(0, args.q_max), histtype="step",
                    label=f"{name} ({len(total)} events, over {n_over})")

        ax.set_yscale("log")
        ax.set_title(f"uPIC{upic} spectrum")
        ax.set_xlabel("charge [ADU x sample]")
        ax.set_ylabel("events")
        ax.legend()


def save_figure(fig, out_path):
    """保存中の不完全な PNG が見えないよう、書き終えてから置き換える。"""

    temporary = f"{out_path}.tmp"
    fig.savefig(temporary, dpi=80, format="png")
    os.replace(temporary, out_path)


def draw_figure(boards, title, out_path, args):
    """確認図を描いて PNG に保存する。"""

    fig = plt.figure(figsize=(26, 25), constrained_layout=True)
    grid = fig.add_gridspec(5, 2 * N_UPIC)

    draw_rate(fig.add_subplot(grid[0, :N_UPIC]), boards, args, log=print)
    draw_baseline(fig, fig.add_subplot(grid[0, N_UPIC:]), boards)
    draw_images(fig, grid, 1, boards)
    draw_wave_maps(fig, grid, 2, boards)
    draw_spectra(fig, grid, 3, boards, args)

    # 5 段目: dE/dx
    tracked = []

    for upic in range(N_UPIC):

        selected = {
            side: [(offset, boards[(side, upic * N_BOARDS_PER_UPIC + offset)])
                   for offset in range(N_BOARDS_PER_UPIC)
                   if (side, upic * N_BOARDS_PER_UPIC + offset) in boards]
            for side in "ac"
        }

        trig, length, charge = track_lengths(selected, args)
        tracked.append(trig)
        print(f"[UPIC {upic}] {len(trig)} events with a track")

        ax = fig.add_subplot(grid[4, N_UPIC + upic])

        counts, x_edges, y_edges = np.histogram2d(
            charge, length, bins=100,
            range=[[0, args.dedx_q_max], [0, args.length_max]])

        if counts.max() > 0:
            ax.pcolormesh(x_edges, y_edges, np.ma.masked_equal(counts.T, 0),
                          norm=LogNorm(vmin=1), cmap="viridis")

        ax.set_title(f"uPIC{upic} dE/dx")
        ax.set_xlabel("anode charge [ADU x sample]")
        ax.set_ylabel("track length [mm]")

    # 1 イベント表示: 指定がなければ、飛跡のある最初のトリガを描く
    event = args.event
    tracked = np.concatenate(tracked)

    if event is None and len(tracked) > 0:
        event = int(tracked.min())

    if event is not None:
        print(f"[EVENT] N_trg {event}")
        draw_event(fig, grid, boards, event, args)

    fig.suptitle(title)

    save_figure(fig, out_path)
    plt.close(fig)


if __name__ == "__main__":

    main()
