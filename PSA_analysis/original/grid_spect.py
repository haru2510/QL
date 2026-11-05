#!/usr/bin/env python3

import argparse
import glob

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.interpolate import RBFInterpolator


parser = argparse.ArgumentParser()
parser.add_argument(
    "--all",
    action="store_true",
    help="include saturated data (default: remove saturated data)",
)
args = parser.parse_args()


INPUT_FILES = "./*.h5"
POSITION_CSV = "peak_output/position_.csv"

N_GRID = 8
N_CHANNELS = 6
PEAK_ID_MIN = 25

HG_MIN, HG_MAX = 0.0, 400e3
HG_BINS = 400

N_BINS_X = 300
X_MIN, X_MAX = -0.32, 0.32
N_BINS_Y = 300
Y_MIN, Y_MAX = -0.16, 0.16

N_BINS_CORRECTED = 160
CORRECTED_X_MIN, CORRECTED_X_MAX = -2.0, 11.0
CORRECTED_Y_MIN, CORRECTED_Y_MAX = -2.0, 11.0


position_data = np.genfromtxt(
    POSITION_CSV,
    delimiter=",",
    names=True,
    dtype=None,
    encoding="utf-8",
)

files = sorted(path for path in glob.glob(INPUT_FILES) if not path.endswith("gridcorr_events.h5"))
if not files:
    raise FileNotFoundError(f"No HDF5 files found: {INPUT_FILES}")

x_all = []
y_all = []
ch_all = []
hg_all = []
saturated_all = []

for path in files:
    print(f"[READ] {path}")
    with h5py.File(path, "r") as h5:
        peak_id = h5["peak_id"][:]
        satu_h = h5["satu_h"][:]
        satu_l = h5["satu_l"][:]

        saturated = (satu_h != 0) | (satu_l != 0)
        selected = peak_id > PEAK_ID_MIN
        if not args.all:
            selected &= ~saturated

        x_all.append(h5["x_hg"][:][selected])
        y_all.append(h5["y_hg"][:][selected])
        ch_all.append(h5["ch"][:][selected])
        hg_all.append(h5["hg"][:][selected])
        saturated_all.append(saturated[selected])

x = np.concatenate(x_all).astype(float)
y = np.concatenate(y_all).astype(float)
ch = np.concatenate(ch_all)
hg = np.concatenate(hg_all).astype(float)
saturated = np.concatenate(saturated_all)

corrected_x = np.empty(len(hg), dtype=float)
corrected_y = np.empty(len(hg), dtype=float)
before_maps = {}
after_maps = {}
pixel_spectra = {}
correction_data = {}

fig, axes = plt.subplots(
    2,
    N_CHANNELS,
    figsize=(18, 6),
    constrained_layout=True,
)

for ich in range(N_CHANNELS):
    rows = position_data[position_data["ch"] == ich]

    if len(rows) != N_GRID * N_GRID:
        raise ValueError(f"ch{ich} has {len(rows)} points, expected 64")

    peak_x = rows["x"].astype(float)
    peak_y = rows["y"].astype(float)

    n = np.empty(len(rows), dtype=int)
    m = np.empty(len(rows), dtype=int)

    x_order = np.argsort(peak_x)
    x_groups = np.array_split(x_order, N_GRID)

    # 各グループに(n, m)を割り当てる
    for group, group_indices in enumerate(x_groups, start=1):
        n[group_indices] = group

        group_y = peak_y[group_indices]
        y_order = np.argsort(group_y)
        sorted_indices = group_indices[y_order]
        m[sorted_indices] = np.arange(1, N_GRID + 1)

    warp = RBFInterpolator(
        np.column_stack([peak_x, peak_y]),
        np.column_stack([n, m]),
        kernel="thin_plate_spline",
        smoothing=0,
    )

    selected = ch == ich
    corrected = warp(np.column_stack([x[selected], y[selected]]))
    hg_ch = hg[selected]
    corrected_x[selected] = corrected[:, 0]
    corrected_y[selected] = corrected[:, 1]

    before_map, x_edges, y_edges = np.histogram2d(
        x[selected],
        y[selected],
        bins=[N_BINS_X, N_BINS_Y],
        range=[[X_MIN, X_MAX], [Y_MIN, Y_MAX]],
    )
    after_map, corrected_x_edges, corrected_y_edges = np.histogram2d(
        corrected[:, 0],
        corrected[:, 1],
        bins=N_BINS_CORRECTED,
        range=[
            [CORRECTED_X_MIN, CORRECTED_X_MAX],
            [CORRECTED_Y_MIN, CORRECTED_Y_MAX],
        ],
    )
    before_maps[f"before_ch{ich}"] = before_map
    after_maps[f"after_ch{ich}"] = after_map

    spectra = np.zeros((N_GRID, N_GRID, HG_BINS), dtype=np.int64)

    axes[0, ich].pcolormesh(
        x_edges,
        y_edges,
        before_map.T,
        cmap="viridis",
        norm=LogNorm(vmin=1),
        shading="auto",
    )
    axes[0, ich].set_title(f"ch{ich} before")
    axes[0, ich].set_xlabel("x_hg")
    axes[0, ich].set_ylabel("y_hg")

    axes[1, ich].pcolormesh(
        corrected_x_edges,
        corrected_y_edges,
        after_map.T,
        cmap="viridis",
        norm=LogNorm(vmin=1),
        shading="auto",
    )
    axes[1, ich].set_title(f"ch{ich} after")
    axes[1, ich].set_xlabel("corrected x")
    axes[1, ich].set_ylabel("corrected y")
    axes[1, ich].set_xticks(range(1, N_GRID + 1))
    axes[1, ich].set_yticks(range(1, N_GRID + 1))
    axes[1, ich].set_xticks(np.arange(0.5, N_GRID + 1.5, 1), minor=True)
    axes[1, ich].set_yticks(np.arange(0.5, N_GRID + 1.5, 1), minor=True)
    axes[1, ich].grid(True, which="minor")

    correction_data[f"observed_ch{ich}"] = np.column_stack([peak_x, peak_y])
    correction_data[f"ideal_ch{ich}"] = np.column_stack([n, m])

    # 各ピクセルのhgスペクトル
    for n_value in range(1, N_GRID + 1):
        for m_value in range(1, N_GRID + 1):
            in_pixel = (
                (corrected[:, 0] >= n_value - 0.5)
                & (corrected[:, 0] < n_value + 0.5)
                & (corrected[:, 1] >= m_value - 0.5)
                & (corrected[:, 1] < m_value + 0.5)
            )
            spectra[n_value - 1, m_value - 1] = np.histogram(
                hg_ch[in_pixel],
                bins=HG_BINS,
                range=(HG_MIN, HG_MAX),
            )[0]

    pixel_spectra[f"spectra_ch{ich}"] = spectra

    axes[0, ich].hist2d(
        x[selected],
        y[selected],
        bins=[N_BINS_X, N_BINS_Y],
        range=[[X_MIN, X_MAX], [Y_MIN, Y_MAX]],
        cmap="viridis",
        norm=LogNorm(vmin=1),
    )
    axes[0, ich].set_title(f"ch{ich} before")
    axes[0, ich].set_xlabel("x_hg")
    axes[0, ich].set_ylabel("y_hg")

    axes[1, ich].hist2d(
        corrected[:, 0],
        corrected[:, 1],
        bins=N_BINS_CORRECTED,
        range=[
            [CORRECTED_X_MIN, CORRECTED_X_MAX],
            [CORRECTED_Y_MIN, CORRECTED_Y_MAX],
        ],
        cmap="viridis",
        norm=LogNorm(vmin=1),
    )
    axes[1, ich].set_title(f"ch{ich} after")
    axes[1, ich].set_xlabel("corrected x")
    axes[1, ich].set_ylabel("corrected y")
    axes[1, ich].set_xticks(range(1, N_GRID + 1))
    axes[1, ich].set_yticks(range(1, N_GRID + 1))
    axes[1, ich].set_xticks(np.arange(0.5, N_GRID + 1.5, 1), minor=True)
    axes[1, ich].set_yticks(np.arange(0.5, N_GRID + 1.5, 1), minor=True)
    axes[1, ich].grid(True, which="minor")

    spectrum_fig, spectrum_axes = plt.subplots(
        N_GRID,
        N_GRID,
        figsize=(12, 12),
        constrained_layout=True,
    )

    for n_value in range(1, N_GRID + 1):
        for m_value in range(1, N_GRID + 1):
            spectrum_ax = spectrum_axes[N_GRID - m_value, n_value - 1]

            in_pixel = (
                (corrected[:, 0] >= n_value - 0.5)
                & (corrected[:, 0] < n_value + 0.5)
                & (corrected[:, 1] >= m_value - 0.5)
                & (corrected[:, 1] < m_value + 0.5)
            )

            spectrum_ax.hist(
                hg_ch[in_pixel],
                bins=HG_BINS,
                range=(HG_MIN, HG_MAX),
                histtype="step",
                linewidth=0.8,
            )
            spectrum_ax.set_title(f"({n_value},{m_value})", fontsize=8)
            spectrum_ax.set_xlim(HG_MIN, HG_MAX)
            spectrum_ax.tick_params(labelsize=6)

    spectrum_fig.suptitle(f"ch{ich} pixel spectra")

np.savez_compressed(
    "gridcorr_histograms.npz",
    x_edges=x_edges,
    y_edges=y_edges,
    corrected_x_edges=corrected_x_edges,
    corrected_y_edges=corrected_y_edges,
    hg_edges=np.linspace(HG_MIN, HG_MAX, HG_BINS + 1),
    **before_maps,
    **after_maps,
    **pixel_spectra,
)

with h5py.File("gridcorr_events.h5", "w") as output:
    output.create_dataset("x_hg", data=x.astype(np.float32), compression="gzip")
    output.create_dataset("y_hg", data=y.astype(np.float32), compression="gzip")
    output.create_dataset(
        "x_corrected", data=corrected_x.astype(np.float32), compression="gzip"
    )
    output.create_dataset(
        "y_corrected", data=corrected_y.astype(np.float32), compression="gzip"
    )
    output.create_dataset("hg", data=hg.astype(np.float32), compression="gzip")
    output.create_dataset("ch", data=ch, compression="gzip")
    output.create_dataset("saturated", data=saturated, compression="gzip")
    output.attrs["peak_id_min"] = PEAK_ID_MIN
    output.attrs["satu_h_condition"] = 0 if not args.all else "not applied"
    output.attrs["satu_l_condition"] = 0 if not args.all else "not applied"
    output.attrs["correction"] = "thin_plate_spline"

np.savez_compressed(
    "grid_correction.npz",
    **correction_data,
    kernel="thin_plate_spline",
    smoothing=0.0,
)

print("[SAVE] gridcorr_histograms.npz")
print("[SAVE] gridcorr_events.h5")
print("[SAVE] grid_correction.npz")

fig.suptitle("Butterfly map before and after correction")
plt.show()

