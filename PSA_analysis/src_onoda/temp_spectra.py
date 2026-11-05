#!/usr/bin/env python3
#
# psa_decode の .h5 を読み、イベントを温度で区切って、ch ごとに温度別の hg スペクトルを描く。
#
#   図   : <out-dir>/<tag>_temp_spectra.png   ch ごとに 1 枚、温度の区間ごとのスペクトルを重ねる (面積で規格化)
#   CSV  : <out-dir>/<tag>_temp_peaks.csv     温度の区間、ch ごとのヒット数と光電ピークの位置
#
# 温度はイベントごとの 6 か所 (TL, TR, ML, MR, BL, BR) の平均。ゲインは温度で大きく変わるので
# (2026-11-03 の ha02 で 1 degC あたり約 -6 %)、温度の違う run を合わせるときはこれで区切って見る。
# トリガを全ボードで共有するモード (HA_cal など) では、他のボードのトリガで記録した枠に、時刻のずれた
# 信号が入る。その信号は積分の窓が波形の端で切れて hg が小さく出るので、既定では peak_id がトリガの
# 時刻に当たるヒットだけを使う (範囲は on_time_range で HA ごとにデータから決める)。
#
# 光電ピークの位置は、平滑化したスペクトルの --peak-min より上で、十分な高さがあって最もエネルギーの
# 高い山とする (目安)。ch のスペクトルは 64 ピクセルの和なので、ピクセルごとのゲインの違いで山が
# 肩になることがある。値は図と見比べて使う。
#
# -i に HA のディレクトリ (decode_runs.py の h5_runs/ha06 など) を渡す。
# h5_runs のように HA のディレクトリを並べた親を渡すと、HA ごとに 1 組ずつ作る。
#

import argparse
import glob
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted

# このリポジトリ (analysis_sm3_FM) の場所。既定のパスはここからの相対で決める
REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


DEFAULT_OUTPUT_DIR = os.path.join(REPO_DIR, "check_output", "temp_spectra")

N_CHANNELS = 6

# --all-timing のときは grid_spect.py と同じく、peak_id がこれより大きいヒットを使う
PEAK_ID_MIN = 25


def parse_args():

    parser = argparse.ArgumentParser(description="Draw hg spectra per temperature range")
    parser.add_argument("-i", "--input", required=True,
                        help="directory of .h5 files of one HA, or a parent of HA directories (h5_runs)")
    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--temp-width", type=float, default=2.0, help="width of a temperature range [degC]")
    parser.add_argument("--min-hits", type=int, default=5000,
                        help="ranges with fewer hits in a channel are not drawn [hits]")
    parser.add_argument("--hg-max", type=float, default=100e3, help="upper edge of the spectrum [ADC x sample]")
    parser.add_argument("--hg-bins", type=int, default=250, help="bins of the spectrum")
    parser.add_argument("--peak-min", type=float, default=15e3,
                        help="search the photopeak above this [ADC x sample]")
    parser.add_argument("--all", action="store_true", help="include saturated hits (default: remove them)")
    parser.add_argument("--peak-id-window", type=int, default=None,
                        help="keep hits with |peak_id - most frequent peak_id| <= this [samples] "
                             "(default: extend from the most frequent peak_id while the count is above "
                             "twice the off-time background)")
    parser.add_argument("--all-timing", action="store_true",
                        help=f"keep all hits with peak_id > {PEAK_ID_MIN} (includes off-time pulses)")

    return parser.parse_args()


def find_ha_directories(path):
    """.h5 が直接あればそのディレクトリだけ、無ければ .h5 を持つ子ディレクトリを返す。"""

    if glob.glob(os.path.join(path, "*.h5")):
        return [path]

    return natsorted(
        child for child in glob.glob(os.path.join(path, "*"))
        if os.path.isdir(child) and glob.glob(os.path.join(child, "*.h5"))
    )


def load(directory, args):

    parts = {"hg": [], "ch": [], "temp": [], "peak_id": []}

    for path in natsorted(glob.glob(os.path.join(directory, "*.h5"))):

        if path.endswith(("_adc.h5", "_pixel.h5", "gridcorr_events.h5")):
            continue

        with h5py.File(path, "r") as h5:

            selected = h5["peak_id"][:] > PEAK_ID_MIN
            if not args.all:
                selected &= (h5["satu_h"][:] == 0) & (h5["satu_l"][:] == 0)

            parts["hg"].append(h5["hg"][:][selected])
            parts["ch"].append(h5["ch"][:][selected])
            parts["temp"].append(h5["temp"][:][selected].mean(axis=1))
            parts["peak_id"].append(h5["peak_id"][:][selected])

    data = {name: np.concatenate(value) for name, value in parts.items()}

    if args.all_timing or len(data["peak_id"]) == 0:
        return data, None

    low, high = on_time_range(data["peak_id"], args.peak_id_window)
    on_time = (data["peak_id"] >= low) & (data["peak_id"] <= high)

    return {name: value[on_time] for name, value in data.items()}, (low, high)


def on_time_range(peak_id, window):
    """トリガの時刻に当たる peak_id の範囲 (両端を含む) を返す。

    トリガの時刻は DAQ やボードで違い、温度や時間とともに 1 サンプルほどずれることもあるので、
    HA ごとにデータから求める。window を指定しないときは、最頻値から外側へ、
    時刻のずれた信号の背景の 2 倍を超えるビンまでを含める。背景はトリガの前後で違う
    (後ろの方が多い) ので、前側は最頻値の 20-8 サンプル前、後ろ側は 5-17 サンプル後の中央値を使う。
    """

    counts = np.bincount(peak_id, minlength=50)
    mode = int(counts.argmax())

    if window is not None:
        return mode - window, mode + window

    background_before = np.median(counts[max(0, mode - 20):max(1, mode - 7)])
    background_after = np.median(counts[min(mode + 5, 48):min(mode + 18, 49)])
    low, high = mode, mode
    while low - 1 > PEAK_ID_MIN and counts[low - 1] > 2 * background_before:
        low -= 1
    while high + 1 < len(counts) - 1 and counts[high + 1] > 2 * background_after:
        high += 1

    return low, high


def photopeak(counts, centers, peak_min, min_height=0.3):
    """--peak-min より上の山のうち、最も高い山の min_height 倍以上あって、最もエネルギーが高い山。

    光電ピークはコンプトン部分のこぶより低いことがあるので、単に最も高いビンは使わない。
    """

    # 統計のゆらぎで 1 ビンだけ高い所を拾わないよう、5 ビンで平滑化する
    smoothed = np.convolve(counts, np.ones(5) / 5, mode="same")
    above = centers >= peak_min

    if not np.any(smoothed[above] > 0):
        return np.nan

    local_max = np.zeros(len(smoothed), dtype=bool)
    local_max[1:-1] = (smoothed[1:-1] >= smoothed[:-2]) & (smoothed[1:-1] > smoothed[2:])
    candidates = np.flatnonzero(local_max & above & (smoothed >= min_height * smoothed[above].max()))

    if len(candidates) == 0:
        return centers[above][np.argmax(smoothed[above])]

    return centers[candidates[-1]]


def process(directory, args):

    tag = os.path.basename(os.path.normpath(directory))
    data, on_time = load(directory, args)

    if len(data["hg"]) == 0:
        print(f"[WARN] {tag}: no hits")
        return

    timing = f"peak_id {on_time[0]}-{on_time[1]}" if on_time is not None else f"peak_id > {PEAK_ID_MIN}"

    width = args.temp_width
    lows = np.arange(np.floor(data["temp"].min() / width) * width, data["temp"].max(), width)
    edges = np.linspace(0, args.hg_max, args.hg_bins + 1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    # 低温を青、高温を黄にする (中間の色も白に近くならない色の並び)
    colors = plt.cm.plasma(np.linspace(0, 0.85, max(len(lows), 2)))

    print(f"[READ] {tag}: {len(data['hg'])} hits ({timing}), "
          f"temp {data['temp'].min():.1f} - {data['temp'].max():.1f} degC")

    fig, axes = plt.subplots(2, 3, figsize=(16, 8), constrained_layout=True)
    rows = []

    for ich, ax in enumerate(axes.ravel()):

        in_channel = data["ch"] == ich

        for low, color in zip(lows, colors):

            selected = in_channel & (data["temp"] >= low) & (data["temp"] < low + width)
            n_hits = int(selected.sum())

            if n_hits < args.min_hits:
                continue

            counts, _ = np.histogram(data["hg"][selected], bins=edges)
            peak = photopeak(counts, centers, args.peak_min)
            rows.append((ich, low, low + width, n_hits, peak))

            ax.step(edges[:-1], counts / n_hits, where="post", color=color,
                    label=f"{low:.0f}-{low + width:.0f} degC: peak {peak / 1e3:.1f}k ({n_hits})")

        ax.set_yscale("log")
        ax.set_xlim(0, args.hg_max)
        ax.set_title(f"{tag} ch{ich}")
        ax.set_xlabel("hg [ADC x sample]")
        ax.set_ylabel("fraction / bin")
        if ax.get_lines():
            ax.legend(fontsize=7)
        else:
            ax.set_title(f"{tag} ch{ich} (no range with {args.min_hits} hits)")

    fig.suptitle(f"{tag}: hg spectrum per temperature range (mean of 6 sensors, {width:g} degC steps, {timing})")

    os.makedirs(args.out_dir, exist_ok=True)
    fig_path = os.path.join(args.out_dir, f"{tag}_temp_spectra.png")
    fig.savefig(fig_path, dpi=80)
    plt.close(fig)

    csv_path = os.path.join(args.out_dir, f"{tag}_temp_peaks.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("ch,temp_low,temp_high,hits,photopeak\n")
        for ich, low, high, n_hits, peak in rows:
            f.write(f"{ich},{low:g},{high:g},{n_hits},{peak:.6g}\n")

    for ich in range(N_CHANNELS):
        peaks = ", ".join(f"{low:.0f}degC {peak / 1e3:.1f}k" for c, low, _, _, peak in rows if c == ich)
        print(f"  ch{ich}: {peaks}")

    print(f"[SAVE] {fig_path}, {csv_path}")


def main():

    args = parse_args()

    directories = find_ha_directories(args.input)

    if len(directories) == 0:
        raise FileNotFoundError(f"No .h5 files found: {args.input}")

    for directory in directories:
        process(directory, args)


if __name__ == "__main__":

    main()
