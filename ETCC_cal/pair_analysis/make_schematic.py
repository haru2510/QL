#!/usr/bin/env python3
"""README 用の模式図: XZ・YZ の2本ずつの枝から、組み方 A と B で e+e- 面の向きが φ と −φ に分かれることを描く。

  python3 ETCC_cal/pair_analysis/make_schematic.py   → ETCC_cal/pair_analysis/figures/pairing_schematic.png
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pipeline.plot_style import use_japanese_font

# 例: 粒子1 と 粒子2 の傾き (dx/dz, dy/dz)。本当の組み方は A とする。
P1 = (-0.10, -0.06)
P2 = (0.15, 0.12)
BLUE, ORANGE, GREY, RED = "#2a78d6", "#eb6834", "#7a7a7a", "#e34948"


def main():
    use_japanese_font()
    a = sorted([P1[0], P2[0]])          # XZ の2本（傾きの小さい順）: a0, a1
    b = sorted([P1[1], P2[1]])          # YZ の2本: b0, b1
    z = np.array([0, 100])
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    ax = axes[0]
    for k, (slope, color) in enumerate(zip(a, (BLUE, ORANGE))):
        ax.plot(z, slope * z, color=color, lw=2.5, label=f"a{k}  (dx/dz = {slope:+.2f})")
    ax.plot(0, 0, "r*", ms=12)
    ax.set(xlabel="Z（ドリフト方向、ビーム方向）[mm]", ylabel="X [mm]", title="XZ（アノードのストリップ）: 2本の枝 a0, a1", ylim=(-20, 20))
    ax.legend(frameon=False)

    ax = axes[1]
    for k, slope in enumerate(b):
        ax.plot(z, slope * z, color=GREY, lw=2.5, ls="-" if k == 0 else "--", label=f"b{k}  (dy/dz = {slope:+.2f})")
    ax.plot(0, 0, "r*", ms=12)
    ax.set(xlabel="Z（ドリフト方向、ビーム方向）[mm]", ylabel="Y [mm]", title="YZ（カソードのストリップ）: 2本の枝 b0, b1\na0 と同じ粒子はどちら?", ylim=(-20, 20))
    ax.legend(frameon=False)

    ax = axes[2]  # ビームに垂直な面（傾き = 横方向の向き）
    dx, dy = a[1] - a[0], b[1] - b[0]
    for label, sign, color, style in (("組み方 A: a0-b0, a1-b1", +1, RED, "-"), ("組み方 B: a0-b1, a1-b0", -1, GREY, "--")):
        phi = np.degrees(np.arctan2(sign * dy, dx))
        ends = np.array([[a[0], b[0] if sign > 0 else b[1]], [a[1], b[1] if sign > 0 else b[0]]])
        ax.plot(ends[:, 0], ends[:, 1], style, color=color, lw=2)
        ax.scatter(ends[:, 0], ends[:, 1], color=color, s=60, zorder=3)
        ax.annotate(f"{label}\nφ = atan2({'+' if sign > 0 else '-'}Δy', Δx') = {phi:+.0f}°", ends[1], textcoords="offset points",
                    xytext=(8, -4 if sign > 0 else -30), color=color, fontsize=9)
    ax.axhline(0, color="0.85", lw=0.8)
    ax.axvline(0, color="0.85", lw=0.8)
    ax.set(xlabel="dx/dz（ビームに垂直な向き）", ylabel="dy/dz", title="ビームに垂直な面: e+e- 面は φ（A）か -φ（B）",
           xlim=(-0.2, 0.5), ylim=(-0.2, 0.2))
    ax.set_aspect("equal")
    ax.text(0.02, 0.02, "組み方が分からないと |φ|（折り返した ψ、0〜90°）しか分からない:\n"
            "1 + μ cos 2(φ - φ0)  →  1 + μ cos 2φ0 · cos 2ψ（φ0 = 45° なら平ら）",
            transform=ax.transAxes, fontsize=8, color="0.3")
    fig.suptitle("XZ と YZ の枝の組み方が曖昧な理由（2本の飛跡を2つの2次元の投影で読み出す）")
    target = Path(__file__).parent / "figures/pairing_schematic.png"
    fig.savefig(target, dpi=120)
    print("saved", target)


if __name__ == "__main__":
    main()
