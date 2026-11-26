#!/usr/bin/env python3
"""README 用: DD の枝の組み方に使う相関係数の分布を、保存済みの stage2・stage3 の結果から描く。

  python3 ETCC_cal/pair_analysis/make_pairing_stats.py [--config ...] [--cut-set default]
  → ETCC_cal/pair_analysis/figures/pairing_correlation.png と、標準出力に数値

左: 選んだ組み方の相関 max(r_A, r_B) と、選ばなかった組み方の相関 min(r_A, r_B) の分布（全データセットの DD）。
右: r_A − r_B の分布。|r_A − r_B| >= min_r_difference なら相関で決める（灰色の帯の外）。
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pipeline.config import Pipeline
from pipeline.plot_style import use_japanese_font
from pipeline.stage3_classify import DD, classify, load_sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).parent / "config/20261124_25_polarization.json")
    parser.add_argument("--cut-set", default="default")
    args = parser.parse_args()
    use_japanese_font()
    pipeline = Pipeline(args.config)
    cuts = pipeline.cut_sets[args.cut_set]
    threshold = pipeline.pairing["min_r_difference"]
    r_a, r_b = [], []
    for name in pipeline.datasets:
        data = load_sample(pipeline, pipeline.samples[name])
        base = data["beam_distance"] <= cuts["beam_radius_mm"]
        if cuts["z_veto"]:
            base &= ~data["z_vetoed"]
        cls, _ = classify(data["single"], data["double"], cuts)
        r = data["pair_r"][base & (cls == DD)]
        with np.errstate(invalid="ignore"):
            r_a.append(np.nanmean(np.stack([r[:, 0, 0], r[:, 1, 1]]), axis=0))
            r_b.append(np.nanmean(np.stack([r[:, 0, 1], r[:, 1, 0]]), axis=0))
    r_a, r_b = np.concatenate(r_a), np.concatenate(r_b)
    ok = np.isfinite(r_a) & np.isfinite(r_b)
    better, worse, diff = np.maximum(r_a, r_b)[ok], np.minimum(r_a, r_b)[ok], (r_a - r_b)[ok]
    decided = np.abs(diff) >= threshold
    print(f"DD (all datasets): {len(r_a):,}, both r_A and r_B available: {ok.sum():,}")
    print(f"  chosen pairing r: mean {better.mean():.3f}, median {np.median(better):.3f}")
    print(f"  other pairing r : mean {worse.mean():.3f}, median {np.median(worse):.3f}")
    print(f"  |r_A - r_B|: median {np.median(np.abs(diff)):.3f};  >= {threshold:g}: {decided.sum():,} ({100 * decided.mean():.0f}%)")
    print(f"  A chosen: {np.sum(diff >= threshold):,}, B chosen: {np.sum(diff <= -threshold):,}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    bins = np.linspace(-1, 1, 41)
    ax = axes[0]
    ax.hist(better, bins=bins, histtype="step", lw=2, color="#2a78d6", label=f"選んだ組み方: 平均 {better.mean():.2f}")
    ax.hist(worse, bins=bins, histtype="step", lw=2, color="#7a7a7a", label=f"選ばなかった組み方: 平均 {worse.mean():.2f}")
    ax.set(xlabel="クロックごとのストリップ数の相関（2組の平均: r_A または r_B）", ylabel="DD のイベント数",
           title=f"クロックごとのストリップ数の相関（DD {ok.sum():,} 件）")
    ax.legend(frameon=False, fontsize=9)
    ax = axes[1]
    ax.hist(diff, bins=np.linspace(-2, 2, 41), color="#2a78d6", alpha=0.85)
    ax.axvspan(-threshold, threshold, color="0.85", zorder=0)
    ax.set(xlabel="r_A - r_B（正: 組み方 A = a0-b0, a1-b1、負: 組み方 B = a0-b1, a1-b0）", ylabel="DD のイベント数",
           title=f"相関で決まったもの（|r_A - r_B| >= {threshold:g}、灰色の帯の外）: "
                 f"{decided.sum():,} / {ok.sum():,} 件（{100 * decided.mean():.0f}%）")
    ax.title.set_fontsize(10)
    target = Path(__file__).parent / "figures/pairing_correlation.png"
    fig.savefig(target, dpi=120)
    print("saved", target)


if __name__ == "__main__":
    main()
