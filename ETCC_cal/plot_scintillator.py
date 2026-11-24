#!/usr/bin/env python3
"""保存済みの全rawシンチヒット件数からブロック別XY分布を描く。"""
import argparse  # 設定JSONの指定を受け取る。
import json  # ブロック中心とヒット件数を読む。
from pathlib import Path  # 設定ファイルのパスを扱う。
import numpy as np  # 位置と件数を配列にする。
import matplotlib  # 保存用の描画環境を用意する。
matplotlib.use("Agg")  # GUIを使わず画像を保存する。
import matplotlib.pyplot as plt  # ブロック位置の散布図を描く。
from settings import DEFAULT_CONFIG, load_config, input_file, output_file  # 入出力は他のrun06・07コードと同じJSONを使う。


def main():  # JSONを読んで全rawヒット分布を再描画する。
    parser = argparse.ArgumentParser(description=__doc__)  # 引数の説明を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 入出力の設定を受け取る。
    args = parser.parse_args()  # 外部設定の指定を読む。
    config = load_config(args.config)  # 共通JSONを読む。
    data = json.loads(input_file(config, "scintillator_counts").read_text())  # 保存済み全ヒット件数を取り出す。
    blocks = [block for block in data["blocks"] if block["count"] > 0]  # 一度もヒットしていないブロックを色の集計から除く。
    positions = np.array([block["position_mm"] for block in blocks])  # ブロック中心のX/Y/Zを取り出す。
    counts = np.array([block["count"] for block in blocks])  # rawヒット数を取り出す。
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)  # 位置分布の図を用意する。
    image = ax.scatter(positions[:, 0], positions[:, 1], c=counts, s=40 + 160 * (positions[:, 2] > 330), cmap="viridis", marker="s", edgecolor="black", lw=0.3)  # 底面を大きい四角、側面を小さい四角で示す。
    fig.colorbar(image, ax=ax, label="all raw HA hits")  # 色が全rawヒット件数であることを明示する。
    ax.set(xlabel="X [mm] (E)", ylabel="Y [mm] (S)", title=", ".join(config["runs"]) + ": all raw HA hits\nNo multiplicity / HG / peak_id cuts; large squares: bottom")  # 2ヒット限定の図と区別する。
    ax.invert_yaxis()  # Yが増える方向を下向きにする。
    ax.set_aspect("equal", adjustable="datalim")  # XとYの距離の比を保つ。
    destination = output_file(config, "scintillator_figure")  # 出力先を外部JSONから読む。
    fig.savefig(destination, dpi=150)  # ブロック別分布を保存する。
    plt.close(fig)  # 描画メモリを解放する。
    print(f"Saved {destination}: {counts.sum():,} raw HA hits")  # 保存先と件数を報告する。


if __name__ == "__main__":  # importだけでは再描画しない。
    main()  # コマンド実行時にだけ処理する。
