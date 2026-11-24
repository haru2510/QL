#!/usr/bin/env python3
"""保存した頂点をXYで選別せず投影する。座標の単位はmm。"""
import argparse  # 入出力パスをコマンドラインで受け取る。
from pathlib import Path  # OSに依存しないパス操作を使う。
import numpy as np  # CSVの座標とヒストグラムを扱う。
import matplotlib  # ファイル保存用の描画環境を設定する。
matplotlib.use("Agg")  # 画面を開かず画像ファイルに描画する。
import matplotlib.pyplot as plt  # XY・XZ・YZの3枚を描く。
from matplotlib.colors import LogNorm  # 少数と多数のビンを対数の色で表示する。
from settings import DEFAULT_CONFIG, load_config, output_file  # 入出力を外部JSONから読む。


def main():  # コマンド実行時の処理をまとめる。
    parser = argparse.ArgumentParser(description=__doc__)  # 入力と出力の説明を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 共通の外部JSONを指定する。
    parser.add_argument("--input", type=Path)  # 頂点CSVのパスを受け取る。
    parser.add_argument("--out", type=Path)  # 保存する画像のパスを受け取る。
    parser.add_argument("--selection-label")  # 別の選別条件のCSVを読む場合は図の条件説明も明示する。
    args = parser.parse_args()  # 実際に指定されたパスを読み取る。
    config = load_config(args.config)  # JSONから入力・出力の対応を読む。
    args.input = args.input or output_file(config, "vertices_selected")  # 選別済みCSVを既定の入力にする。
    args.out = args.out or output_file(config, "vertex_figure")  # JSONで指定した図の保存先を使う。
    vertices = np.loadtxt(args.input, delimiter=",", skiprows=2, usecols=(5, 6, 7))  # コメントと見出しを除きX・Y・Z列を読む。
    vertices = np.atleast_2d(vertices)  # 1イベントだけでも行と列を持つ配列にする。
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)  # 3方向の投影を横に並べる。
    projections = [(0, 1, "XY"), (0, 2, "XZ"), (1, 2, "YZ")]  # 各図で使う座標列を指定する。
    labels = ["X [mm] (E)", "Y [mm] (S)", "Z [mm] (down)"]  # 座標軸の向きと単位を明示する。
    for ax, (i, j, name) in zip(axes, projections):  # 各投影に同じ全頂点を渡す。
        x_edges = np.linspace(-164, 164, 83)  # 横方向を検出器全体より広い表示範囲にする。
        y_edges = np.linspace(-164, 164, 83) if j == 1 else np.linspace(-50, 365, 84)  # Zは負座標から高Z端まで表示する。
        counts, _, _ = np.histogram2d(vertices[:, i], vertices[:, j], bins=(x_edges, y_edges))  # 追加の座標選別をせず全頂点を数える。
        masked = np.ma.masked_less(counts.T, 1)  # 空ビンは対数色に割り当てず白くする。
        image = ax.pcolormesh(x_edges, y_edges, masked, cmap="turbo", norm=LogNorm(vmin=1))  # 1点以上のビンを対数色で描く。
        fig.colorbar(image, ax=ax, label="vertices / bin")  # 色が頂点数であることを示す。
        ax.set(xlabel=labels[i], ylabel=labels[j], title=f"{name}: all input vertices")  # XYの追加条件がないことをタイトルに示す。
        ax.invert_yaxis()  # 従来図と同様にYおよびZが増える方向を下向きにする。
        ax.set_aspect("equal", adjustable="box")  # 縦横の1 mmを同じ長さにする。
    selection = args.selection_label or f"Event veto: any raw TPC hit at Z <= {config['z_veto']:g} mm; no additional XY or high-Z selection"  # 入力CSVに対応する条件説明を使う。
    fig.suptitle(", ".join(config["runs"]) + f": N = {len(vertices):,}\n{selection}")  # 事象除外と表示条件を区別する。
    args.out.parent.mkdir(parents=True, exist_ok=True)  # 保存先がなければ作る。
    fig.savefig(args.out, dpi=150)  # 投影図をPNGで保存する。
    plt.close(fig)  # 描画用メモリを解放する。
    print(f"Saved {args.out}: {len(vertices):,} vertices")  # 保存先と使用した頂点数を報告する。


if __name__ == "__main__":  # importだけでは処理を実行しない。
    main()  # コマンドとして呼ばれたときだけ描画する。
