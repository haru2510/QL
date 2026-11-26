#!/usr/bin/env python3
"""README 用: 枝を追うときの「塊」の分け方と割り当てを、実際の DD イベント（1つの投影）で描く。raw を読む（SSD が要る）。

  python3 ETCC_cal/pair_analysis/make_cluster_figure.py [--run 20261124_run06] [--event 1181] [--axis 0] [--zoom 0 70]
  → ETCC_cal/pair_analysis/figures/pairing_clusters.png（塊の分け方と割り当て）
    ETCC_cal/pair_analysis/figures/pairing_ranges.png（円の範囲と、予測に使う Z 方向の帯）

左: 投影の全体（灰色の点 = ストリップ、色の線 = 枝に割り当てた塊の中心、四角 = 拡大する範囲）。
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
from pipeline.plot_style import use_japanese_font
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
    parser.add_argument("--predict-at", type=float, default=100.0, help="予測の例を描く Z [mm]（pairing_ranges.png）")
    parser.add_argument("--predict-branch", type=int, default=0, help="予測の例を描く枝")
    args = parser.parse_args()
    use_japanese_font()

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
        ax.plot(paths[k][:, 0], paths[k][:, 1], color=COLORS[k], lw=1.5, label=f"枝{k}（端の Z = {ends[k]:.0f} mm）")
    ax.plot(0, 0, "r*", ms=10)
    lo, hi = args.zoom
    visible = (points[:, 0] >= lo) & (points[:, 0] <= hi)
    ylo, yhi = points[visible, 1].min() - 2, points[visible, 1].max() + 2
    ax.add_patch(plt.Rectangle((lo, ylo), hi - lo, yhi - ylo, fill=False, ls="--", color="0.3"))
    ax.set(xlabel="Z - Z_散乱点 [mm]", ylabel=f"{name} - {name}_散乱点 [mm]", title=f"{args.run} event {args.event}、{name}Z: 追った枝")
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
    ax.set(xlim=(lo, hi), ylim=(ylo, yhi), xlabel="Z - Z_散乱点 [mm]（縦の線1本 = 1クロックの1つの塊）", ylabel=f"{name} - {name}_散乱点 [mm]",
           title="拡大: クロックごとの塊（青・橙 = 枝0・枝1 に割り当て、黒 = 割り当てなし）\n"
                 f"点線 = 予測位置、薄い帯 = ± {window:g} mm、灰色 = 2本の予測位置が {settings['shared_mm']:g} mm より近い（使わない）")
    ax.title.set_fontsize(9)

    ax = axes[2]  # 1クロックの切り出し
    r = example
    owner = {g: k for k, g in r["assigned"].items()}
    ax.scatter(np.zeros(len(r["strips"])), r["strips"], marker="_", s=300, color="0.4")
    for g, group in enumerate(r["groups"]):
        color = COLORS[owner[g]] if g in owner else "black"
        ax.plot([0.25, 0.25], [group.min() - 0.15, group.max() + 0.15], color=color, lw=5, solid_capstyle="butt")
        ax.text(0.33, group.mean(), f"塊: {len(group)} 本" + (f" → 枝{owner[g]}" if g in owner else "（割り当てなし）"),
                va="center", fontsize=8, color=color)
    for k in range(2):
        if r["alive"][k]:
            p = r["predicted"][k]
            ax.plot([-0.35, -0.05], [p, p], color=COLORS[k], lw=1.5)
            ax.fill_between([-0.35, -0.05], p - window, p + window, color=COLORS[k], alpha=0.12, lw=0)
            ax.text(-0.37, p, f"枝{k}の\n予測位置", ha="right", va="center", fontsize=8, color=COLORS[k])
    ax.set(xlim=(-0.9, 1.6), xticks=[], ylabel=f"{name} - {name}_散乱点 [mm]",
           title=f"1クロックの切り出し（Z = {r['dz']:.1f} mm、紫の線）:\n灰色 = 鳴ったストリップ（0.4 mm 間隔）、\n{settings['cluster_gap_mm']:g} mm より空いたところで塊に分ける")
    ax.title.set_fontsize(9)
    fig.suptitle("2本の枝をクロックごとに追う: 各クロックのストリップを塊に分け、予測位置に近い塊を各枝に割り当てる")
    target = Path(__file__).parent / "figures/pairing_clusters.png"
    fig.savefig(target, dpi=120)
    print("saved", target)
    plot_ranges(args, pipeline, saved, row, points, paths, name)


def plot_ranges(args, pipeline, saved, row, points, paths, name):
    """円の範囲（scan.radius_mm、pairing.context_mm）と、予測に使う Z 方向の帯（refit_length_mm）を描く。"""
    settings = pipeline.pairing
    radius, context, refit = pipeline.scan["radius_mm"], settings["context_mm"], settings["refit_length_mm"]
    window = settings["window_mm"]
    fig, axes = plt.subplots(1, 2, figsize=(17, 5.8), constrained_layout=True, gridspec_kw={"width_ratios": [1, 1.25]})

    ax = axes[0]  # 円の範囲
    ax.scatter(points[:, 0], points[:, 1], s=2, color="#b8b7b2", lw=0)
    theta = np.linspace(-np.pi / 2, np.pi / 2, 200)
    ax.plot(radius * np.cos(theta), radius * np.sin(theta), color="#1baf7a", lw=1.5,
            label=f"scan.radius_mm = {radius:g} mm（円）: 出発点の2本 fit に使うストリップ")
    ax.plot(context * np.cos(theta), context * np.sin(theta), color="#7b4fd0", lw=1.5, ls="--",
            label=f"pairing.context_mm = {context:g} mm（円）: 枝を追うのに使うストリップ")
    z = np.array([0, context])
    for k in range(2):
        slope = saved["double"][row, args.axis, k]
        ax.plot(z, slope * z, color=COLORS[k], lw=1, ls=":", label=f"枝{k}の出発点の直線（2本 fit）")
    ax.plot(0, 0, "r*", ms=10)
    ax.set(xlim=(-5, context + 10), ylim=(-context / 2, context / 2), xlabel="Z - Z_散乱点 [mm]", ylabel=f"{name} - {name}_散乱点 [mm]",
           title="円の範囲: 投影の平面での散乱点からの距離 √(dz² + t²)")
    ax.set_aspect("equal")
    ax.legend(fontsize=7, frameon=False, loc="lower left")

    ax = axes[1]  # Z 方向の帯で予測
    k = args.predict_branch
    path = paths[k]
    z0 = args.predict_at
    before = path[path[:, 0] < z0]
    last = before[-1, 0]
    used = before[before[:, 0] >= last - refit]
    slope, intercept = np.polyfit(used[:, 0], used[:, 1], 1)
    predicted = slope * z0 + intercept
    lo, hi = last - refit - 25, z0 + 15
    near = (points[:, 0] >= lo) & (points[:, 0] <= hi)
    ax.scatter(points[near, 0], points[near, 1], s=3, color="#b8b7b2", lw=0, label="ストリップ")
    ax.axvspan(last - refit, last, color=COLORS[k], alpha=0.10, lw=0,
               label=f"refit_length_mm = {refit:g} mm: Z 方向の帯（dz >= 最後の中心 - {refit:g} mm）、円ではない")
    other = before[before[:, 0] < last - refit]
    ax.plot(other[:, 0], other[:, 1], "o", ms=3, color=COLORS[k], alpha=0.35, label=f"枝{k}の前の方の中心（使わない）")
    ax.plot(used[:, 0], used[:, 1], "o", ms=3.5, color=COLORS[k], label=f"fit に使う中心（{len(used)} 個）")
    zz = np.linspace(last - refit, z0, 50)
    ax.plot(zz, slope * zz + intercept, color="black", lw=1.2, label="最小二乗の直線 t = 傾き·dz + 切片")
    ax.errorbar([z0], [predicted], yerr=[window], fmt="s", color="black", ms=6, capsize=5,
                label=f"Z = {z0:g} mm での予測位置 ± window_mm（{window:g} mm）")
    ax.axvline(last, color="0.4", lw=0.8, ls=":")
    ax.text(last, used[:, 1].max() + 3, " 最後の中心", fontsize=8, color="0.3")
    ax.set(xlim=(lo, hi), xlabel="Z - Z_散乱点 [mm]", ylabel=f"{name} - {name}_散乱点 [mm]",
           title=f"枝{k}の Z = {z0:g} mm での予測: その枝の塊の中心のうち、Z で直近 {refit:g} mm のものに直線を当てる")
    ax.legend(fontsize=7, frameon=False, loc="best")
    fig.suptitle(f"{args.run} event {args.event}、{name}Z: 枝を追うときに使う範囲")
    target = Path(__file__).parent / "figures/pairing_ranges.png"
    fig.savefig(target, dpi=120)
    print("saved", target)


if __name__ == "__main__":
    main()
