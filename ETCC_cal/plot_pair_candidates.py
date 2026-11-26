#!/usr/bin/env python3
"""tpc_pair_modulation.py の分類（SS・DS・SD・DD）ごとに、散乱点まわりのXZ・YZのrawストリップと、fitした直線を描く。

保存済みの <result_dir>/with_z0_cut/tpc_pair_branches_r<radius>/<run>.npz から分類し、各分類の先頭 --count 件の
raw を読み直す。灰色は --context mm 以内の全ストリップ、色付きは fit に使う --radius mm 以内（散乱点より3 mm以上下）。
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from settings import load_config, repo_path
from stage_coincidence_xy import scan_run
from tpc_pair_modulation import classify
from tpc_projected_pair import projected_hits
from vertex_sphere import take_local_cells

NAMES = {1: "SS", 2: "DS", 3: "SD", 4: "DD"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).parent / "config/20261124_run00_07_horizontal.json")
    parser.add_argument("--run", default="20261124_run06")
    parser.add_argument("--classes", nargs="+", default=["DD", "DS", "SD", "SS"])
    parser.add_argument("--count", type=int, default=6)
    parser.add_argument("--skip", type=int, default=0, help="各分類の先頭からこの件数を飛ばす")
    parser.add_argument("--radius", type=float, default=40.0)
    parser.add_argument("--context", type=float, default=100.0, help="灰色で描く範囲 [mm]")
    parser.add_argument("--min-span", type=float, default=10.0)
    parser.add_argument("--min-branch-span", type=float, default=25.0)
    parser.add_argument("--min-separation", type=float, default=6.0)
    parser.add_argument("--min-inlier-fraction", type=float, default=0.7)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    folder = repo_path(config["result_dir"]) / f"with_z0_cut/tpc_pair_branches_r{args.radius:g}"
    saved = np.load(folder / f"{args.run}.npz")
    ids, vertex, single, double = saved["ids"], saved["vertex"], saved["single"], saved["double"]
    cls, psi = classify(single, double, args)
    wanted = {k for k, v in NAMES.items() if v in args.classes}
    picks = {k: np.flatnonzero(cls == k)[args.skip:args.skip + args.count] for k in sorted(wanted)}
    index = np.concatenate(list(picks.values()))
    lookup = {int(ids[i]): i for i in index}
    calibration = tuple(config[k] for k in ("velocity", "clock_offset", "z_gem"))
    hits = {}

    def process(event, upic, occupancy, occupied, numbers):
        i = lookup.get(event)
        if i is None:
            take_local_cells(occupancy, occupied, numbers, np.zeros(3), 0.0, *calibration)
            return 0
        cells = take_local_cells(occupancy, occupied, numbers, vertex[i], args.context, *calibration)
        h = projected_hits(cells, vertex[i], args.context, *calibration)
        if h:
            hits.setdefault(i, []).append(np.asarray(h))
        return len(h)

    scan_run(args.run, config, [np.sort(ids[index])], event_processor=process)
    out = args.out or folder / "event_displays"
    out.mkdir(parents=True, exist_ok=True)
    colors = ["#2a78d6", "#eb6834"]
    for k, rows in picks.items():
        fig, axes = plt.subplots(len(rows), 2, figsize=(11, 3.3 * len(rows)), squeeze=False, constrained_layout=True)
        for r, i in enumerate(rows):
            h = np.unique(np.concatenate(hits[i]), axis=0) if i in hits else np.zeros((0, 3))
            for axis, name in enumerate("XY"):
                ax = axes[r, axis]
                p = h[h[:, 0] == axis]
                near = np.hypot(p[:, 1], p[:, 2]) <= args.radius
                ax.scatter(p[~near, 1], p[~near, 2], s=3, color="#b8b7b2", lw=0)
                ax.scatter(p[near, 1], p[near, 2], s=4, color="black", lw=0)
                zz = np.linspace(0, args.radius, 2)
                two = double[i, axis]
                is_two = (k == 4) or (k == 2 and axis == 0) or (k == 3 and axis == 1)
                if is_two:
                    for b in range(2):
                        ax.plot(zz, two[b] * zz, color=colors[b], lw=1.4, alpha=0.8)
                    label = (f"2 branches: {np.degrees(np.arctan(two[0])):.1f}°, {np.degrees(np.arctan(two[1])):.1f}° "
                             f"(sep {np.degrees(np.arctan(two[1]) - np.arctan(two[0])):.1f}°), "
                             f"in {int(two[2])}+{int(two[3])}/{int(single[i, axis, 2])}")
                else:
                    ax.plot(zz, single[i, axis, 0] * zz, color="#1baf7a", lw=1.4, alpha=0.8)
                    label = f"1 line: {np.degrees(np.arctan(single[i, axis, 0])):.1f}°, in {int(single[i, axis, 1])}/{int(single[i, axis, 2])}"
                ax.add_patch(plt.Circle((0, 0), args.radius, fill=False, ls=":", color="0.5"))
                ax.plot(0, 0, "r*", ms=9)
                lim = args.context
                ax.set(xlim=(-5, lim), ylim=(-lim / 2, lim / 2), xlabel="Z − Z_vertex [mm] (drift, beam direction)",
                       ylabel=f"{name} − {name}_vertex [mm]")
                ax.set_aspect("equal")
                ax.set_title(f"{name}Z: {label}", fontsize=8)
            extra = f", ψ = {psi[i]:.0f}°" if k == 4 else ""
            axes[r, 0].text(0.01, 0.97, f"{args.run} event {ids[i] >> 32}:{ids[i] & 0xffffffff}\n"
                            f"vertex ({vertex[i, 0]:.1f}, {vertex[i, 1]:.1f}, {vertex[i, 2]:.1f}) mm{extra}",
                            transform=axes[r, 0].transAxes, va="top", fontsize=7)
        fig.suptitle(f"{NAMES[k]} events ({args.run}); black: strips used (r ≤ {args.radius:g} mm), grey: up to {args.context:g} mm")
        target = out / f"{args.run}_{NAMES[k]}_{args.skip}.png"
        fig.savefig(target, dpi=110)
        plt.close(fig)
        print("Saved", target)


if __name__ == "__main__":
    main()
