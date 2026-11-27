#!/usr/bin/env python3
"""README 用: 組み方を決めた DD の φ_pair の結果（水平・45°・垂直、分母の日による違い、同じ偏光同士）を1枚の日本語の図にまとめる。

保存済みの stage4（cut set default）の pair_phi.json だけを読む。設定の comparisons に、次の名前の比較が要る:
  <horizontal|diagonal45|vertical>_over_circular、..._over_circular_1124、..._over_circular_1126、
  null_circular_1126_over_1124、null_diagonal45_run04_over_run06

  python3 ETCC_cal/pair_analysis/make_summary_figure.py   → ETCC_cal/pair_analysis/figures/summary_pair_phi.png
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pipeline.config import Pipeline
from pipeline.plot_style import use_japanese_font

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--config", type=Path, default=Path(__file__).parent / "config/20261124_25_polarization.json")
parser.add_argument("--cut-set", default="default")
parser.add_argument("--out", type=Path, default=Path(__file__).parent / "figures/summary_pair_phi.png")
args = parser.parse_args()
BASE = Pipeline(args.config).output_dir / "stage4_modulation" / args.cut_set
use_japanese_font()


def load(name):
    d = json.loads((BASE / name / "pair_phi.json").read_text())
    a, b = (np.array(c, float) for c in d["counts"])
    edges = np.array(d["phi_edges_deg"])
    r = (a / a.sum()) / (b / b.sum())
    e = r * np.sqrt(1 / a + 1 / b)
    return 0.5 * (edges[1:] + edges[:-1]), r, e, d["fit"], d["events"], 0.5 * np.diff(edges)


def model(x, f):
    return f["A"] * (1 + f["mu"] * np.cos(2 * np.radians(x - f["phi0_deg"])))


def wrap(x):
    return (x + 90) % 180 - 90


POL = [("horizontal", "水平", 0.0), ("diagonal45", "45° のデータ（検出器から −45°）", 135.0), ("vertical", "垂直", 90.0)]
DEN = [("circular_1124", "分母 = 円偏光 1124 だけ", "#2a78d6", "{}_over_circular_1124"),
       ("circular_1126", "分母 = 円偏光 1126 run03 だけ", "#eb6834", "{}_over_circular_1126"),
       ("circular", "分母 = 合わせた円偏光（本命）", "#222222", "{}_over_circular")]
xx = np.linspace(-15, 165, 361)

fig, axes = plt.subplots(2, 3, figsize=(18, 10.5), constrained_layout=True)
for ax, (key, label, expect) in zip(axes[0], POL):
    x, r, e, f, ev, hw = load(f"{key}_over_circular")
    ax.errorbar(x, r, xerr=hw, yerr=e, fmt="o", color="#222222", ms=5, capsize=0, label="比（それぞれ合計で規格化）")
    ax.plot(xx, model(xx, f), color="#e34948", lw=2,
            label=f"fit: μ = {f['mu']:.3f} ± {f['mu_err']:.3f}\nφ0 = {f['phi0_deg'] % 180:.0f} ± {f['phi0_err_deg']:.0f}°（{f['significance_sigma']:.1f}σ）")
    ax.axvline(expect, color="#1baf7a", ls="--", lw=2, label=f"予想の φ0 = {expect:.0f}°")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set(xlim=(-15, 165), ylim=(0.6, 1.4), xticks=range(0, 166, 15), xlabel="φ_pair [deg]（対の面の向き。0 = X、90 = Y）",
           ylabel="比", title=f"{label} ÷ 合わせた円偏光（組み方が決まった DD {ev[0]['decided']:,} / {ev[1]['decided']:,}）")
    ax.title.set_fontsize(10)
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    if key == "horizontal":
        ax.text(0.98, 0.03, "0°・90° のまわりは受け入れが小さい（片方の投影で枝が重なる）ので 30° のビン。\n"
                "90° のビンは fit より 2σ 上。外しても φ0 は 1° も変わらない",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="0.3")

# 下段左: φ0 − 予想（分母ごと）
ax = axes[1, 0]
for j, (_, dlabel, color, pattern) in enumerate(DEN):
    diffs = []
    for i, (key, label, expect) in enumerate(POL):
        f = load(pattern.format(key))[3]
        d, s = wrap(f["phi0_deg"] - expect), f["phi0_err_deg"]
        diffs.append((d, s))
        ax.errorbar(i + (j - 1) * 0.2, d, yerr=s, fmt="o", color=color, ms=7, capsize=3)
    w = np.array([1 / s ** 2 for _, s in diffs]); v = np.array([d for d, _ in diffs])
    m, me = (w * v).sum() / w.sum(), w.sum() ** -0.5
    chi2 = float(((v - m) ** 2 * w).sum())
    ax.errorbar(3 + (j - 1) * 0.2, m, yerr=me, fmt="D", color=color, ms=7, capsize=3,
                label=f"{dlabel}: δ = {m:+.1f} ± {me:.1f}°（χ² = {chi2:.1f}/2）")
ax.axhline(0, color="0.5", lw=0.8)
ax.set(xticks=range(4), xticklabels=["水平", "45°", "垂直", "共通のずれ δ\n（重み付き平均）"], ylabel="測った φ0 − 予想の φ0 [deg]",
       title="偏光の向きのずれ: 分母の日によって水平と垂直が逆に動くが、δ はほぼ同じ", ylim=(-15, 25))
ax.title.set_fontsize(10)
ax.legend(fontsize=8, frameon=False, loc="upper left")

# 下段中: μ（分母ごと）
ax = axes[1, 1]
for j, (_, dlabel, color, pattern) in enumerate(DEN):
    for i, (key, label, expect) in enumerate(POL):
        f = load(pattern.format(key))[3]
        ax.errorbar(i + (j - 1) * 0.2, f["mu"], yerr=f["mu_err"], fmt="o", color=color, ms=7, capsize=3, label=dlabel if i == 0 else None)
ax.set(xticks=range(3), xticklabels=["水平", "45°", "垂直"], ylabel="μ（組み方が決まった DD）", ylim=(0, 0.32),
       title="モジュレーションの大きさ μ: 45° は分母で 0.20 と 0.12 に分かれる")
ax.title.set_fontsize(10)
ax.legend(fontsize=8, frameon=False, loc="upper right")

# 下段右: 同じ偏光同士の確認
ax = axes[1, 2]
for name, label, color, dx in [("null_circular_1126_over_1124", "円偏光 1126 run03 ÷ 1124", "#7b4fd0", -1.5),
                               ("null_diagonal45_run04_over_run06", "45° 1126 run04 ÷ 1125 run06", "#1baf7a", 1.5)]:
    x, r, e, f, ev, hw = load(name)
    ax.errorbar(x + dx, r, yerr=e, fmt="o", color=color, ms=5, capsize=0)
    ax.plot(xx, model(xx, f), color=color, lw=1.8,
            label=f"{label}: μ = {f['mu']:.3f} ± {f['mu_err']:.3f}、φ0 = {f['phi0_deg'] % 180:.0f}°（{f['significance_sigma']:.1f}σ）")
ax.axhline(1, color="0.6", lw=0.8)
ax.set(xlim=(-15, 165), ylim=(0.6, 1.4), xticks=range(0, 166, 15), xlabel="φ_pair [deg]", ylabel="比",
       title="同じ偏光同士（平らのはず）: 45° は平ら、円偏光は 2.6σ ずれる")
ax.title.set_fontsize(10)
ax.legend(fontsize=8, frameon=False, loc="upper left")

fig.suptitle("6.6 MeV 電子対の偏光（統計を増やした再解析）: 円偏光 = 1124 run08-12 + 1126 run03、45° = 1125 run06 + 1126 run04。"
             "誤差は統計だけ。0°・90° のまわり（受け入れが小さい）は 30° のビン", fontsize=13)
fig.savefig(args.out, dpi=100)
print("saved", args.out)
