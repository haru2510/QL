#!/usr/bin/env python3
#
# grid_spect.py が保存した位置補正（grid_correction.npz）を新しいデータの .h5 に当てはめ、
# 各ヒットがどのピクセル (n, m) に入ったかを求める。
#
# 補正の関数は、grid_correction.npz の対応点（observed_chN → ideal_chN）から
# grid_spect.py と同じ RBFInterpolator（薄板スプライン）で作り直すので、
# grid_spect.py を実行したときと同じ変換になる。
#
# ピクセルの判定は grid_spect.py と同じ半開区間:
#   n - 0.5 <= x_corrected < n + 0.5,  m - 0.5 <= y_corrected < m + 0.5
# 8x8 の外に出たヒットは pix_n = pix_m = -1。
#

import argparse
import glob
import os

import h5py
import numpy as np
from natsort import natsorted
from scipy.interpolate import RBFInterpolator

# このリポジトリ (analysis_sm3_FM) の場所。既定のパスはここからの相対で決める
REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


DEFAULT_CALIBRATION = os.path.join(REPO_DIR, "PSA_analysis", "grid_output", "grid_correction.npz")
DEFAULT_OUTPUT_DIR = os.path.join(REPO_DIR, "PSA_analysis", "pixel_store")

N_GRID = 8
N_CHANNELS = 6

# psa_decode.py の .h5 からそのまま写すデータセット
COPY_KEYS = [
    "trig", "event_time", "ch", "hitnum",
    "hg", "lg", "x_hg", "y_hg",
    "peak_id", "satu_h", "satu_l",
]

# 新しく作るデータセットの説明と単位
NEW_DATASETS = {
    "x_corrected": ("補正後の x。ピクセル n の中心が整数 n になる", "pixel"),
    "y_corrected": ("補正後の y。ピクセル m の中心が整数 m になる", "pixel"),
    "pix_n": ("ピクセルの列番号 1-8（x 方向）。8x8 の外は -1", "index"),
    "pix_m": ("ピクセルの行番号 1-8（y 方向）。8x8 の外は -1", "index"),
}


def parse_args():

    parser = argparse.ArgumentParser(
        description="Assign each hit in psa_decode .h5 files to a pixel using grid_correction.npz"
    )
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help=".h5 files, directories or glob patterns (psa_decode output)")
    parser.add_argument("-c", "--calibration", default=DEFAULT_CALIBRATION,
                        help="grid_correction.npz made by grid_spect.py")
    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUTPUT_DIR)

    return parser.parse_args()


def is_event_h5(path):

    name = os.path.basename(path)

    return (
        name.endswith(".h5")
        and not name.endswith(("_adc.h5", "_pixel.h5", "gridcorr_events.h5"))
    )


def collect_inputs(inputs):

    files = []

    for item in inputs:

        if os.path.isdir(item):
            files += [os.path.join(item, name) for name in os.listdir(item)]
        elif os.path.isfile(item):
            files.append(item)
        else:
            files += glob.glob(item)

    return natsorted(set(path for path in files if is_event_h5(path)))


def load_warps(calibration_path):

    calibration = np.load(calibration_path)

    kernel = str(calibration["kernel"])
    smoothing = float(calibration["smoothing"])

    warps = {}

    for ich in range(N_CHANNELS):

        observed = calibration[f"observed_ch{ich}"].astype(float)
        ideal = calibration[f"ideal_ch{ich}"].astype(float)

        warps[ich] = RBFInterpolator(observed, ideal, kernel=kernel, smoothing=smoothing)

    return warps


def to_pixel(corrected):

    # n - 0.5 <= x < n + 0.5  <=>  n = floor(x + 0.5)
    index = np.floor(corrected + 0.5).astype(np.int16)
    index[(index < 1) | (index > N_GRID)] = -1

    return index


def process_file(path, warps, out_dir, calibration_path):

    with h5py.File(path, "r") as h5:
        data = {key: h5[key][:] for key in COPY_KEYS}

    n_rows = len(data["trig"])

    x_corrected = np.full(n_rows, np.nan)
    y_corrected = np.full(n_rows, np.nan)

    for ich in range(N_CHANNELS):

        selected = data["ch"] == ich

        if not np.any(selected):
            continue

        corrected = warps[ich](
            np.column_stack([data["x_hg"][selected], data["y_hg"][selected]]).astype(float)
        )
        x_corrected[selected] = corrected[:, 0]
        y_corrected[selected] = corrected[:, 1]

    pix_n = to_pixel(x_corrected)
    pix_m = to_pixel(y_corrected)

    # x か y のどちらかが外なら両方 -1
    outside = (pix_n < 0) | (pix_m < 0)
    pix_n[outside] = -1
    pix_m[outside] = -1

    new = {
        "x_corrected": x_corrected.astype(np.float32),
        "y_corrected": y_corrected.astype(np.float32),
        "pix_n": pix_n,
        "pix_m": pix_m,
    }

    stem = os.path.splitext(os.path.basename(path))[0]
    out_path = os.path.join(out_dir, f"{stem}_pixel.h5")

    with h5py.File(out_path, "w") as out:

        for key, value in data.items():
            out.create_dataset(key, data=value, compression="gzip")

        for key, value in new.items():
            dataset = out.create_dataset(key, data=value, compression="gzip")
            dataset.attrs["description"] = NEW_DATASETS[key][0]
            dataset.attrs["unit"] = NEW_DATASETS[key][1]

        out.attrs["source_h5"] = os.path.abspath(path)
        out.attrs["calibration"] = os.path.abspath(calibration_path)
        out.attrs["pixel_rule"] = "n - 0.5 <= x_corrected < n + 0.5 (same for m); outside 8x8 -> -1"

    return out_path, n_rows, data["ch"], outside


def main():

    args = parse_args()

    files = collect_inputs(args.input)

    if len(files) == 0:
        raise FileNotFoundError(f"No psa_decode .h5 files found: {args.input}")

    os.makedirs(args.out_dir, exist_ok=True)

    print(f"[CAL] {args.calibration}")
    warps = load_warps(args.calibration)

    for idx, path in enumerate(files, start=1):

        out_path, n_rows, ch, outside = process_file(path, warps, args.out_dir, args.calibration)

        inside_ratio = ", ".join(
            f"ch{ich} {100 * np.mean(~outside[ch == ich]):.1f}%"
            for ich in range(N_CHANNELS) if np.any(ch == ich)
        )

        print(f"[{idx}/{len(files)}] {out_path} : {n_rows} hits, in 8x8: {inside_ratio}")


if __name__ == "__main__":

    main()
