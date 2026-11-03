#!/usr/bin/env python3
#
# psa_decode.py の .h5 をまとめて読み、データが正しそうかを見るための確認図を 1 枚の PNG に描く。
#
#   1-3 段目（ch ごと）: バタフライマップ (x_hg, y_hg)、hg の分布、peak_id の分布
#   4 段目            : トリガの時間間隔、4 信号ごとのペデスタル（HG0-3, LG0-3）
#
# 画面には出さず PNG に保存するだけなので、XLaunch は不要。
#

import argparse
import glob
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from natsort import natsorted


DEFAULT_OUTPUT_DIR = "/home/onoda/PSA_analysis/check_output"

N_CHANNELS = 6
N_SAMPLES = 50

# peak_search.py と同じバタフライマップの範囲 [x_hg, y_hg は無次元の比]
X_MIN, X_MAX, N_BINS_X = -0.32, 0.32, 250
Y_MIN, Y_MAX, N_BINS_Y = -0.16, 0.16, 300

# hg の範囲 [ADC カウント x サンプル]
HG_MIN, HG_MAX, HG_BINS = 0.0, 400e3, 400

KEYS = ["trig", "event_time", "ch", "hg", "x_hg", "y_hg", "peak_id", "ped", "satu_h"]


def parse_args():

    parser = argparse.ArgumentParser(description="Draw sanity-check plots for psa_decode .h5 files")
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help=".h5 files, directories or glob patterns (psa_decode output)")
    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--tag", default=None, help="name of the output PNG (default: input directory name)")

    return parser.parse_args()


def collect_inputs(inputs):

    files = []

    for item in inputs:

        if os.path.isdir(item):
            files += [os.path.join(item, name) for name in os.listdir(item)]
        elif os.path.isfile(item):
            files.append(item)
        else:
            files += glob.glob(item)

    return natsorted(set(
        path for path in files
        if path.endswith(".h5") and not path.endswith(("_adc.h5", "_pixel.h5", "gridcorr_events.h5"))
    ))


def load(files):

    data = {key: [] for key in KEYS}
    intervals = []
    trig_first = []

    for path in files:

        with h5py.File(path, "r") as h5:

            part = {key: h5[key][:] for key in KEYS}

        for key in KEYS:
            data[key].append(part[key])

        # 時間間隔はファイルの中だけで計算する（event_time はファイルごとに 0 始まり）
        trig, first = np.unique(part["trig"], return_index=True)
        intervals.append(np.diff(part["event_time"][first]) * 1e6)
        trig_first.append(trig)

    data = {key: np.concatenate(value) for key, value in data.items()}

    return data, np.concatenate(intervals), trig_first


def main():

    args = parse_args()

    files = collect_inputs(args.input)

    if len(files) == 0:
        raise FileNotFoundError(f"No psa_decode .h5 files found: {args.input}")

    tag = args.tag or os.path.basename(os.path.normpath(args.input[0]))

    data, intervals, trigs = load(files)

    n_events = sum(len(t) for t in trigs)
    n_missing = sum(int((np.diff(t) - 1).clip(0).sum()) for t in trigs)

    print(f"[READ] {len(files)} files, {n_events} events, {len(data['ch'])} hits")
    print(f"[TRIG] missing inside files: {n_missing} ({n_missing / (n_events + n_missing) * 100:.3f} %)")
    print(f"[TIME] interval median {np.median(intervals):.1f} us")

    fig, axes = plt.subplots(4, N_CHANNELS, figsize=(24, 15), constrained_layout=True)

    for ich in range(N_CHANNELS):

        selected = data["ch"] == ich
        n_hits = int(selected.sum())

        print(f"[CH {ich}] hits {n_hits}")

        ax = axes[0, ich]
        if n_hits > 0:
            ax.hist2d(data["x_hg"][selected], data["y_hg"][selected],
                      bins=[N_BINS_X, N_BINS_Y], range=[[X_MIN, X_MAX], [Y_MIN, Y_MAX]],
                      norm=LogNorm(vmin=1), cmap="viridis")
        ax.set_title(f"ch{ich} butterfly ({n_hits} hits)")
        ax.set_xlabel("x_hg")
        ax.set_ylabel("y_hg")

        ax = axes[1, ich]
        ax.hist(data["hg"][selected], bins=HG_BINS, range=(HG_MIN, HG_MAX), histtype="step")
        ax.set_yscale("log")
        ax.set_title(f"ch{ich} hg (saturated {data['satu_h'][selected].mean() * 100 if n_hits else 0:.1f} %)")
        ax.set_xlabel("hg [ADC x sample]")

        ax = axes[2, ich]
        ax.hist(data["peak_id"][selected], bins=N_SAMPLES, range=(-0.5, N_SAMPLES - 0.5), histtype="step")
        ax.set_title(f"ch{ich} peak_id")
        ax.set_xlabel("peak_id [sample]")

    for ich in range(N_CHANNELS):
        axes[3, ich].remove()

    grid = axes[3, 0].get_gridspec()

    ax = fig.add_subplot(grid[3, 0:3])
    ax.hist(intervals, bins=200, range=(0, max(1.0, np.percentile(intervals, 99.5))), histtype="step")
    ax.set_yscale("log")
    ax.set_title(f"trigger interval (median {np.median(intervals):.1f} us, missing trig {n_missing})")
    ax.set_xlabel("interval [us]")

    ax = fig.add_subplot(grid[3, 3:6])
    labels = ["HG0", "HG1", "HG2", "HG3", "LG0", "LG1", "LG2", "LG3"]
    for ich in range(N_CHANNELS):
        selected = data["ch"] == ich
        if selected.sum() < 10:
            continue
        ax.plot(labels, np.median(data["ped"][selected], axis=0), marker="o", label=f"ch{ich}")
    ax.set_title("pedestal median per signal")
    ax.set_ylabel("pedestal [ADC]")
    ax.legend()

    fig.suptitle(f"{tag}: {len(files)} files, {n_events} events")

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f"{tag}_check.png")
    fig.savefig(out_path, dpi=80)

    print(f"[SAVE] {out_path}")


if __name__ == "__main__":

    main()
