#!/usr/bin/env python3
#
# .dat の波形を ch ごとに重ねて描く（積算波形図）。
#
#   上段: HG0-3 の平均、下段: LG0-3 の平均（どちらもペデスタルを引く前の生の ADC 値）
#   色: 全イベントを重ねた密度（2 次元ヒストグラム、対数）、白線: 中央値の波形
#   青の帯: ペデスタルの区間（先頭 23 サンプル）
#   橙の帯: 典型的な積分窓（最頻の peak_id - 4 〜 + 9）
#
# 波形は .h5 に入っていないので、.dat を sub_scr/base_decode で直接デコードする。
# メモリを抑えるため、先頭から --max-events イベントだけを使う。
#

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from natsort import natsorted

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sub_scr import base_decode

# このリポジトリ (analysis_sm3_FM) の場所。既定のパスはここからの相対で決める
REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


DEFAULT_OUTPUT_DIR = os.path.join(REPO_DIR, "PSA_analysis", "check_output")

N_CHANNELS = 6
N_SAMPLES = 50

# 温度は描かないので既定値でよい [LSB/degC], [LSB]
TEMP_A = np.ascontiguousarray(np.full(6, 3.925820834, dtype=np.float32))
TEMP_B = np.ascontiguousarray(np.full(6, 1230.29635, dtype=np.float32))


def parse_args():

    parser = argparse.ArgumentParser(description="Overlay raw waveforms of .dat files per channel")
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help=".dat files, directories or glob patterns")
    parser.add_argument("-n", "--max-events", type=int, default=5000,
                        help="number of events to use, from the first file (default 5000)")
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

    return natsorted(set(path for path in files if path.lower().endswith(".dat")))


def load_waveforms(files, max_events):

    waves, hits, peaks = [], [], []
    n_total = 0

    for path in files:

        events_buf, event_sizes, n_events = base_decode.read_file(path)

        if n_events == 0:
            continue

        n_use = min(n_events, max_events - n_total)

        out = base_decode.decode_events(
            events_buf[:n_use], event_sizes[:n_use], False, TEMP_A, TEMP_B, True
        )

        hits.append(out[2])
        peaks.append(out[13])
        waves.append(out[17])

        n_total += n_use

        if n_total >= max_events:
            break

    return np.concatenate(waves), np.concatenate(hits), np.concatenate(peaks), n_total


def draw_panel(ax, wave, peak_id, title):

    # wave: (イベント数, 50) の ADC 値
    samples = np.tile(np.arange(N_SAMPLES), len(wave))
    values = wave.ravel()

    lo, hi = np.percentile(values, [0.1, 99.9])
    margin = 0.05 * (hi - lo + 1)

    ax.hist2d(samples, values, bins=[N_SAMPLES, 200],
              range=[[-0.5, N_SAMPLES - 0.5], [lo - margin, hi + margin]],
              norm=LogNorm(vmin=1), cmap="viridis")
    ax.plot(np.arange(N_SAMPLES), np.median(wave, axis=0), color="white", linewidth=1.5, label="median")

    ax.axvspan(-0.5, base_decode.PEDESTAL_RANGE - 0.5, color="tab:blue", alpha=0.15, label="pedestal")

    if len(peak_id) > 0:
        mode = int(np.bincount(peak_id, minlength=N_SAMPLES).argmax())
        ws = max(0, mode - base_decode.WAVE_FRONT)
        we = min(N_SAMPLES - 1, mode + base_decode.WAVE_BACK)
        ax.axvspan(ws - 0.5, we + 0.5, color="tab:orange", alpha=0.2, label=f"signal ({ws}-{we})")

    ax.set_title(title)
    ax.set_xlabel("sample")
    ax.set_ylabel("ADC")
    ax.legend(fontsize=7, loc="upper right")


def main():

    args = parse_args()

    files = collect_inputs(args.input)

    if len(files) == 0:
        raise FileNotFoundError(f"No .dat files found: {args.input}")

    tag = args.tag or os.path.basename(os.path.normpath(args.input[0]))

    waves, hits, peaks, n_total = load_waveforms(files, args.max_events)

    print(f"[READ] {n_total} events from {files[0]} ...")

    fig, axes = plt.subplots(2, N_CHANNELS, figsize=(24, 8), constrained_layout=True)

    for out_ch in range(N_CHANNELS):

        ich = 5 - out_ch   # .h5 の出力 ch = 5 - 内部 ch

        selected = hits[:, ich] == 1
        n_hit = int(selected.sum())

        print(f"[CH {out_ch}] {n_hit} events")

        if n_hit == 0:
            axes[0, out_ch].set_title(f"ch{out_ch} HG (no events)")
            axes[1, out_ch].set_title(f"ch{out_ch} LG (no events)")
            continue

        wave = waves[selected, ich]          # (イベント数, 50, 10)
        hg = wave[:, :, 2:6].mean(axis=2)    # HG0-3 の平均
        lg = wave[:, :, 6:10].mean(axis=2)   # LG0-3 の平均
        peak_id = peaks[selected, ich]

        draw_panel(axes[0, out_ch], hg, peak_id, f"ch{out_ch} HG mean of 4 ({n_hit} events)")
        draw_panel(axes[1, out_ch], lg, peak_id, f"ch{out_ch} LG mean of 4")

    fig.suptitle(f"{tag}: raw waveforms of the first {n_total} events")

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f"{tag}_waveform.png")
    fig.savefig(out_path, dpi=80)

    print(f"[SAVE] {out_path}")


if __name__ == "__main__":

    main()
