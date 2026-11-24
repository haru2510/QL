#!/usr/bin/env python3
"""保存した候補点群JSONから、TPC投影とシンチ位置のPNG・回転できるHTMLを描く。"""
import argparse  # 共通JSONまたは個別の入力・出力を指定する。
import json  # 保存した点群とヒット情報を読む。
from pathlib import Path  # 入出力ファイル名を扱う。
import numpy as np  # 点群と座標を配列として扱う。
import matplotlib  # 保存用の描画環境を指定する。
matplotlib.use("Agg")  # GUIを使わずPNGへ描画する。
import matplotlib.pyplot as plt  # 3Dと投影を同じ図に並べる。
import plotly.graph_objects as go  # ブラウザで回転できる3D図も保存する。
from settings import DEFAULT_CONFIG, load_config, input_file, output_file  # 入出力を共通JSONから読む。


def main():  # 保存済み候補の再描画を行う。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行時の指定を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 他のrun06・07解析と同じJSONを受け取る。
    parser.add_argument("--input", type=Path)  # 必要なら別候補のJSONで上書きする。
    parser.add_argument("--out", type=Path)  # 必要なら別の保存名で上書きする。
    args = parser.parse_args()  # 実際の指定を読み取る。
    config = load_config(args.config)  # 共通設定を読む。
    source = args.input or input_file(config, "event_points")  # 入力する保存済み点群を決める。
    prefix = args.out or output_file(config, "event_figure_prefix")  # PNGとHTMLの共通保存名を決める。
    prefix.parent.mkdir(parents=True, exist_ok=True)  # 保存先を作る。
    event = json.loads(source.read_text())  # トリガーとTPC点群・シンチヒットを読む。
    points = np.array(event["tpc_combinations_mm"])  # 同一クロックのXY組み合わせ点を取り出す。
    vertex = np.array(event["vertex_mm"])  # 最初の同時計数による推定頂点を取り出す。
    hits = np.array([hit["position_mm"] for hit in event["ha_hits"]])  # シンチブロック中心の座標を取り出す。
    labels = [hit["label"] for hit in event["ha_hits"]]  # 現行chと旧出力chの表示名を取り出す。
    colors = ["red", "darkorange"]  # 2個のシンチヒットを同じ色で対応づける。
    fig = plt.figure(figsize=(16, 6))  # 3枚のパネルを描く領域を用意する。
    ax = fig.add_subplot(131, projection="3d")  # 左側に3D表示を置く。
    ax.scatter(*points.T, s=3, c=points[:, 2], cmap="viridis", alpha=0.4)  # 曖昧さを残したTPC点群をZの色で描く。
    for hit, color, label in zip(hits, colors, labels):  # 同じトリガーのシンチを重ねる。
        ax.scatter(*hit, s=100, c=color, marker="s", label=label)  # シンチ位置はブロック中心として描く。
    ax.scatter(*vertex, s=70, c="magenta", marker="*", label="Estimated vertex")  # 推定頂点を星印で描く。
    ax.set(xlabel="X [mm]", ylabel="Y [mm]", zlabel="Z [mm]", xlim=(-210, 210), ylim=(-210, 210), zlim=(365, 0))  # 側面シンチまで含む範囲を使う。
    ax.view_init(elev=20, azim=-65)  # 固定PNGで見やすい視点を選ぶ。
    ax.legend(fontsize=7)  # シンチと頂点の意味を示す。
    for index, axis in enumerate(("x", "y")):  # 直接観測されたXZとYZも並べる。
        ax = fig.add_subplot(1, 3, index + 2)  # 中央または右パネルを選ぶ。
        projection = np.array(event["projections_mm"][axis])  # rawストリップの投影点を取り出す。
        ax.scatter(projection[:, 0], projection[:, 1], s=2, color="steelblue", alpha=0.4)  # XYの組み合わせを使わず投影を描く。
        ax.scatter(vertex[index], vertex[2], c="magenta", s=70, marker="*")  # 同じ頂点を投影する。
        for hit, color, label in zip(hits, colors, labels):  # シンチヒットも同じ平面へ投影する。
            ax.scatter(hit[index], hit[2], c=color, s=100, marker="s")  # シンチの位置を四角で示す。
            ax.annotate(label, (hit[index], hit[2]), xytext=(-100, -10), textcoords="offset points", fontsize=7, color=color)  # ブロック番号を付ける。
        ax.set(xlabel=axis.upper() + " [mm]", ylabel="Z [mm] (down)", xlim=(-210, 210), ylim=(365, 0), title=axis.upper() + "Z raw projection")  # 投影の方向と単位を示す。
        ax.grid(alpha=0.2)  # 位置関係を見やすくする。
    title = f"{event['run']} / trigger {event['trig']}: two-track topology candidate"  # 対生成の確定ではなく形状候補と明示する。
    fig.suptitle(title + "\nSame-trigger TPC + HA; 3D includes X/Y ambiguities; HA positions are block centres")  # 3D対応の限界を図に残す。
    fig.tight_layout()  # 軸ラベルの重なりを抑える。
    fig.savefig(prefix.with_suffix(".png"), dpi=140)  # 固定視点のPNGを保存する。
    plt.close(fig)  # 描画メモリを解放する。
    interactive = go.Figure()  # 回転可能な3D図を用意する。
    interactive.add_trace(go.Scatter3d(x=points[:, 0], y=points[:, 1], z=points[:, 2], mode="markers", name="TPC combinations (ambiguous)", marker=dict(size=2, color=points[:, 2], colorscale="Viridis", opacity=0.4)))  # 曖昧な点を確定した2本の線にはしない。
    interactive.add_trace(go.Scatter3d(x=hits[:, 0], y=hits[:, 1], z=hits[:, 2], mode="markers+text", text=labels, name="HA block centres", marker=dict(size=8, color=colors, symbol="square")))  # シンチの番号をホバーと文字で確認できるようにする。
    interactive.update_layout(title=title, scene=dict(xaxis_title="X [mm]", yaxis_title="Y [mm]", zaxis_title="Z [mm]", zaxis=dict(autorange="reversed"), aspectmode="data"))  # 3軸の単位と正負の向きをそろえる。
    interactive.write_html(prefix.with_suffix(".html"), include_plotlyjs=True)  # ネット接続なしで開けるHTMLとして保存する。
    print(f"Saved {prefix}.png and .html")  # 保存先を報告する。


if __name__ == "__main__":  # importだけで図を作り直さない。
    main()  # コマンドで実行したときだけ描画する。
