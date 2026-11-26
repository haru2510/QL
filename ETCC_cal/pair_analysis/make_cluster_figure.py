#!/usr/bin/env python3
"""README 用: 枝を追うときの「塊」の分け方と割り当てを、実際の DD イベント（1つの投影）で描く。raw を読む（SSD が要る）。

  python3 ETCC_cal/pair_analysis/make_cluster_figure.py [--run 20261124_run06] [--event 1181] [--axis 0] [--zoom 0 70]
  → ETCC_cal/pair_analysis/figures/pairing_clusters.png

左: 投影の全体（灰色の点 = ストリップ、色の線 = 枝に割り当てた塊の重心、四角 = 拡大する範囲）。
中: 拡大。各クロックの塊を縦の線で描き、青 = 枝0、橙 = 枝1 に割り当てた塊、黒 = どちらにも割り当てなかった塊。
    細い点線 = 各枝の予測位置、薄い帯 = 予測位置 ± window_mm。灰色の網かけ = 2本の予測位置が shared_mm より近いクロック（使わない）。
右: 1クロックの切り出し。ストリップ（縦の目盛り）を cluster_gap_mm より離れたところで塊に分け、予測位置に近い塊を割り当てる。
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pipeline.branch_pairing import follow_branches
from pipeline.config import Pipeline
from pipeline.stage5_event_display import read_hits

COLORS = ("#2a78d6", "#eb6834")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).parent / "config/20261124_25_polarization.json")
    parser.add_argument("--dataset", default="horizontal")
    parser.add_argument("--run", default="20261124_run06")
    parser.add_argument("--event", type=int, default=1181, help="trig（epoch 0）")
    parser.add_argument("--axis", type=int, default=0, help="0 = XZ, 1 = YZ")
    parser.add_argument("--zoom", type=float, nargs=2, default=[0.0, 70.0], help="拡大する Z の範囲 [mm]")
    args = parser.parse_args()

    pipeline = Pipeline(args.config)
    settings = pipeline.pairing
    saved = np.load(pipeline.stage_dir(2, args.dataset) / f"{args.run}.npz")
    row = int(np.flatnonzero(saved["event"] == args.event)[0])
    hits = read_hits(pipeline, args.dataset, args.run, saved["event"][[row]], saved["vertex"][[row]], settings["context_mm"])[0]
    points = hits[hits[:, 0] == args.axis][:, 1:]
    trace = []
    mm_per_clock = pipeline.etcc_config(args.dataset)["velocity"]
    _, ends, _, paths = follow_branches(points, saved["double"][row, args.axis, :2], settings, mm_per_clock,
                                        return_paths=True, trace=trace)
    name = "XY"[args.axis]
    window = settings["window_mm"]

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.6), constrained_layout=True, gridspec_kw={"width_ratios": [1.1, 1.6, 0.7]})
    ax = axes[0]
    ax.scatter(points[:, 0], points[:, 1], s=2, color="#b8b7b2", lw=0)
    for k in range(2):
        ax.plot(paths[k][:, 0], paths[k][:, 1], color=COLORS[k], lw=1.5, label=f"branch {k} (end at Z = {ends[k]:.0f} mm)")
    ax.plot(0, 0, "r*", ms=10)
    lo, hi = args.zoom
    visible = (points[:, 0] >= lo) & (points[:, 0] <= hi)
    ylo, yhi = points[visible, 1].min() - 2, points[visible, 1].max() + 2
    ax.add_patch(plt.Rectangle((lo, ylo), hi - lo, yhi - ylo, fill=False, ls="--", color="0.3"))
    ax.set(xlabel="Z − Z_vertex [mm]", ylabel=f"{name} − {name}_vertex [mm]", title=f"{args.run} event {args.event}, {name}Z: followed branches")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    shown = [r for r in trace if lo <= r["dz"] <= hi]
    example = None
    for r in shown:
        z = r["dz"]
        if r["shared"]:
            ax.axvspan(z - mm_per_clock / 2, z + mm_per_clock / 2, color="0.88", lw=0, zorder=0)
        owner = {g: k for k, g in r["assigned"].items()}
        for g, group in enumerate(r["groups"]):
            color = COLORS[owner[g]] if g in owner else "black"
            ax.plot([z, z], [group.min() - 0.2, group.max() + 0.2], color=color, lw=1.6, solid_capstyle="butt")
        for k in range(2):
            if r["alive"][k]:
                ax.plot(z, r["predicted"][k], ".", color=COLORS[k], ms=2)
                ax.fill_between([z - mm_per_clock / 2, z + mm_per_clock / 2], r["predicted"][k] - window, r["predicted"][k] + window,
                                color=COLORS[k], alpha=0.08, lw=0)
        if example is None and not r["shared"] and len(r["groups"]) >= 3 and z > lo + 0.4 * (hi - lo):
            example = r
    if example is None:
        example = next(r for r in shown if not r["shared"])
    ax.axvline(example["dz"], color="#7b4fd0", lw=1, ls=":")
    ax.plot(0, 0, "r*", ms=10)
    ax.set(xlim=(lo, hi), ylim=(ylo, yhi), xlabel="Z − Z_vertex [mm] (one column = one clock)", ylabel=f"{name} − {name}_vertex [mm]",
           title="zoom: clusters per clock (blue / orange = assigned to branch 0 / 1, black = not assigned)\n"
                 f"dots = predicted position, band = ± {window:g} mm, grey = branches closer than {settings['shared_mm']:g} mm (not used)")
    ax.title.set_fontsize(9)

    ax = axes[2]  # 1クロックの切り出し
    r = example
    owner = {g: k for k, g in r["assigned"].items()}
    ax.scatter(np.zeros(len(r["strips"])), r["strips"], marker="_", s=300, color="0.4")
    for g, group in enumerate(r["groups"]):
        color = COLORS[owner[g]] if g in owner else "black"
        ax.plot([0.25, 0.25], [group.min() - 0.15, group.max() + 0.15], color=color, lw=5, solid_capstyle="butt")
        ax.text(0.33, group.mean(), f"cluster: {len(group)} strips" + (f" → branch {owner[g]}" if g in owner else " (not assigned)"),
                va="center", fontsize=8, color=color)
    for k in range(2):
        if r["alive"][k]:
            p = r["predicted"][k]
            ax.plot([-0.35, -0.05], [p, p], color=COLORS[k], lw=1.5)
            ax.fill_between([-0.35, -0.05], p - window, p + window, color=COLORS[k], alpha=0.12, lw=0)
            ax.text(-0.37, p, f"branch {k}\npredicted", ha="right", va="center", fontsize=8, color=COLORS[k])
    ax.set(xlim=(-0.9, 1.6), xticks=[], ylabel=f"{name} − {name}_vertex [mm]",
           title=f"one clock (Z = {r['dz']:.1f} mm, purple line):\nstrips (grey ticks, 0.4 mm pitch) split at gaps > {settings['cluster_gap_mm']:g} mm")
    ax.title.set_fontsize(9)
    fig.suptitle("Following two branches clock by clock: split the strips of each clock into clusters and give each branch the nearest cluster")
    target = Path(__file__).parent / "figures/pairing_clusters.png"
    fig.savefig(target, dpi=120)
    print("saved", target)


if __name__ == "__main__":
    main()
