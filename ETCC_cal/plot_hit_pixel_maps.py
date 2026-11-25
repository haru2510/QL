#!/usr/bin/env python3
"""保存済みの1/2ヒットCSVからXY・XZ・YZのピクセル分布を描き、ビーム軸を重ねる。"""
import argparse  # 外部JSONと描画するヒット種別を指定する。
import csv  # 保存済みイベント・ペア一覧を読む。
import json  # 集計とビーム軸の係数を読む。
from pathlib import Path  # 入出力パスを扱う。
import numpy as np  # 投影位置が重なったピクセルの件数を合算する。
import matplotlib  # ファイル出力用の描画環境を設定する。
matplotlib.use("Agg")  # GUIを開かず画像を保存する。
import matplotlib.pyplot as plt  # 3方向の位置分布を描く。
from matplotlib.colors import LogNorm  # 件数の色を対数で表示する。
from settings import DEFAULT_CONFIG, load_config, repo_path  # 全コードと同じ外部JSONを使う。


def beam_reference(config):  # 保存済みfitのビーム直線と底面での通過点を求める。
    options = config["beam_axis_overlay"]  # 表示だけの条件をイベント選別と分けて読む。
    source = repo_path(config["inputs"]["polarization_beam_axis"])  # 他コードで求めた局所ピーク軸を選ぶ。
    saved = json.loads(source.read_text())  # 軸位置と傾きを読む。
    reference_z = float(saved.get("line_reference_z_mm", 155))  # 保存された基準Z [mm]を使う。
    point = np.array([saved["x_line"]["at_z155_mm"], saved["y_line"]["at_z155_mm"], reference_z])  # 軸上の基準点を作る。
    slope = np.array([saved["x_line"]["slope_mm_per_mm"], saved["y_line"]["slope_mm_per_mm"], 1.])  # Zあたりの座標変化を表す。
    marker = point + (options["marker_z_mm"] - reference_z) * slope  # 指定したZ平面を通るビーム位置を求める。
    low, high = options["z_range_mm"]  # 軸を表示するZ範囲を得る。
    z = np.linspace(low, high, 500)  # 投影線をなめらかに描く。
    line = point + (z[:, None] - reference_z) * slope  # 直線上の座標 [mm]を作る。
    slices = saved.get("slices", [])  # 軸fitに実際に使ったZ範囲を確認する。
    fit_low = min(row["z_low_mm"] for row in slices) if slices else low  # fit範囲の下端を得る。
    fit_high = max(row["z_high_mm"] for row in slices) if slices else high  # fit範囲の上端を得る。
    return line, marker, (fit_low, fit_high), str(source)  # 表示と記録に使う情報を返す。


def draw_pixel_map(positions, config, title, count_label, destination, color_max=None, marker_size=12):  # 1/2ヒットで共通の表示方法を使う。
    positions = np.asarray(positions, dtype=float).reshape(-1, 3)  # 空の入力も3列にそろえる。
    projections = []  # 各平面の位置と件数を集める。
    for i, j in ((0, 1), (0, 2), (1, 2)):  # XY・XZ・YZを作る。
        locations, counts = np.unique(positions[:, [i, j]], axis=0, return_counts=True)  # 投影で重なる別深さのピクセルも合算する。
        projections.append((locations, counts))  # 描画と件数確認に使う。
    maximum = max([int(counts.max()) for _, counts in projections if len(counts)] + [2])  # 3方向で同じ色範囲を使う。
    maximum = max(float(color_max), 2) if color_max is not None else maximum  # 段階比較ではカット前の色上限に揃えられるようにする。
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)  # 1ヒット版と同じ3枚構成にする。
    options = config["beam_axis_overlay"]  # ビーム軸の表示設定を読む。
    overlay = beam_reference(config) if options["enabled"] else None  # 表示OFFなら軸ファイルを要求しない。
    metadata = {"input_points": len(positions), "color_max": maximum, "projected_totals": [], "beam_axis_overlay": options, "colormap": config.get("pixel_map", {}).get("colormap", "turbo"), "beam_marker_alpha": options.get("marker_alpha", 0.45)}  # 描画で件数が変わっていないことを残す。
    for ax, (i, j), (locations, counts) in zip(axes, ((0, 1), (0, 2), (1, 2)), projections):  # 各投影を描く。
        image = ax.scatter(locations[:, 0], locations[:, 1], c=counts, marker="s", s=marker_size, cmap=config.get("pixel_map", {}).get("colormap", "turbo"), norm=LogNorm(vmin=1, vmax=maximum))  # ピクセル件数を対数色で表す。
        fig.colorbar(image, ax=ax, label=count_label)  # 1ヒットのイベント数と2ヒットの使用回数を区別する。
        if overlay is not None:  # 保存ビーム軸の印を追加する。
            line, marker, fit_range, source = overlay  # 軸の投影に必要な情報を取り出す。
            if j == 2:  # XZ/YZには直線を重ねる。
                ax.plot(line[:, i], line[:, j], color="magenta", ls="--", lw=1.3, label="Beam axis (extrapolated)")  # fit範囲外は外挿として破線で表示する。
                inside = (line[:, 2] >= fit_range[0]) & (line[:, 2] <= fit_range[1])  # 実際にfitしたZ範囲を選ぶ。
                ax.plot(line[inside, i], line[inside, j], color="magenta", lw=1.8, label="Beam axis (fit Z range)")  # fit範囲内を実線で表示する。
            label = f"Beam at Z={marker[2]:.2f} mm"  # XY上の印がどのZの位置か明示する。
            ax.scatter(marker[i], marker[j], marker="*", s=160, c="magenta", edgecolors="black", linewidths=0.7, label=label, zorder=5, alpha=options.get("marker_alpha", 0.45))  # シンチ位置との比較点を星印で示す。
            ax.legend(fontsize=7, loc="best")  # 星と線の意味を図に書く。
            metadata.update(beam_marker_mm=marker.tolist(), beam_fit_z_range_mm=list(fit_range), beam_axis_file=source)  # 外挿点の数値も保存する。
        ax.set(xlabel=("X", "Y", "Z")[i]+" [mm]", ylabel=("X", "Y", "Z")[j]+" [mm]", title=("X", "Y", "Z")[i]+("X", "Y", "Z")[j])  # 平面名と単位を示す。
        ax.set_xlim(-215, 215)  # 1/2ヒットの比較で同じ横範囲を使う。
        ax.set_ylim((-215, 215) if j == 1 else (-5, 375))  # シンチとビーム投影を両方含む範囲にする。
        ax.invert_yaxis()  # Y/Zが増える方向を下にする。
        ax.set_aspect("equal", adjustable="box")  # 四角の隙間が縦横の縮尺差で縦筋に見えるのを抑える。
        metadata["projected_totals"].append(int(counts.sum()))  # 投影後も入力と同じ件数か記録する。
    fig.suptitle((config["data_label"] + "\n" if config.get("data_label") else "") + title)  # イベント条件・gain・件数を表示する。
    destination = Path(destination)  # 出力先をPathへ変換する。
    destination.parent.mkdir(parents=True, exist_ok=True)  # 保存先を作る。
    np.savez_compressed(destination.with_suffix(".points.npz"), positions=positions)  # 色や印だけ変える際に位置較正を繰り返さず再描画できるようにする。
    fig.savefig(destination, dpi=150)  # 3方向のピクセル図を保存する。
    plt.close(fig)  # 描画メモリを解放する。
    destination.with_suffix(".plot.json").write_text(json.dumps(metadata, indent=2)+"\n")  # 表示条件とビーム位置を保存する。
    return metadata  # 呼び出し側で件数を検証できるようにする。


def main():  # 保存済みCSVからSSDなしで1/2ヒット図を作り直す。
    parser = argparse.ArgumentParser(description=__doc__)  # 再描画の指定を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 全解析と同じJSONを受け取る。
    parser.add_argument("--kind", choices=("single", "pair", "both"), default="both")  # どの分布を再描画するか選ぶ。
    args = parser.parse_args()  # 指定を読み取る。
    config = load_config(args.config)  # 入出力と軸表示を読む。
    kinds = ("single", "pair") if args.kind == "both" else (args.kind,)  # 要求された図だけを処理する。
    for kind in kinds:  # 1ヒットと採用ペアを別々に扱う。
        key = "single_hit_pixel_prefix" if kind == "single" else "pixel_polarization_prefix"  # JSONから対応する保存名を選ぶ。
        prefix = repo_path(config["outputs"][key])  # CSVと図の共通名を得る。
        for gain in ("hg", "lg"):  # 両gainを同じ表示方法にする。
            points = []  # CSVから位置だけを集める。
            suffix = "events.csv" if kind == "single" else "pairs.csv"  # 入力する一覧形式を選ぶ。
            source = Path(str(prefix)+f"_{gain}_{suffix}")  # 入力CSVのパスを決める。
            with source.open() as stream:  # 位置一覧を読む。
                for row in csv.DictReader(stream):  # 全採用イベントまたはペアを使う。
                    if kind == "single":  # 1ヒットは1行1点を使う。
                        points.append([float(row[k]) for k in ("x_mm", "y_mm", "z_mm")])  # 単独のピクセル中心を得る。
                    else:  # ペアは両ヒットを積算する。
                        for side in ("1", "2"):  # 片側だけを落とさず2点とも使う。
                            points.append([float(row[k+side+"_mm"]) for k in ("x", "y", "z")])  # ペアCSVのピクセル中心を得る。
            summary = json.loads(Path(str(prefix)+f"_{gain}.json").read_text())  # 実際に適用したイベント条件を読む。
            if kind == "single":  # 1ヒットの件数をタイトルにする。
                title = f"{gain.upper()}: exactly one raw HA hit; N={len(points):,}\nHG position calibration (LG provisional); beam cut={summary['options']['beam_vertex_cut_enabled']}; z0 veto={summary['options'].get('z0_track_veto_enabled', False)}"  # 保存時の条件を表示する。
                label = "single-hit events / projected pixel"  # 1点が1イベントであることを示す。
                destination = Path(str(prefix)+f"_{gain}.png")  # 1ヒット図をビーム印付きに更新する。
            else:  # 2ヒット条件や多重ペアの違いをタイトルに明記する。
                pairs = len(points)//2  # CSVのペア数を得る。
                title = f"{gain.upper()}: raw2={summary['options']['require_exactly_two_raw_hits']}; {pairs:,} pairs / {len(points):,} hit uses\nHG position calibration (LG provisional); beam cut={summary['options']['beam_vertex_cut_enabled']}; z0 veto={summary['options'].get('z0_track_veto_enabled', False)}"  # 2点ずつ積算したことを明示する。
                label = "hit uses in selected pairs / projected pixel"  # イベント数と2倍のヒット使用数を混同させない。
                destination = Path(str(prefix)+f"_{gain}_pixels.png")  # 方位角の図を残して専用の3投影図を追加する。
            result = draw_pixel_map(points, config, title, label, destination)  # 共通の描画とビーム印を使う。
            assert result["projected_totals"] == [len(points)]*3  # 投影でヒット数が増減していないか確認する。
            print(f"Saved {destination}", flush=True)  # 図の保存先を知らせる。


if __name__ == "__main__":  # importだけではファイルを読まない。
    main()  # コマンドとして呼ばれたときだけ再描画する。
