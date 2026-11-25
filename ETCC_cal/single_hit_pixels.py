#!/usr/bin/env python3
"""rawフラグがちょうど1ヒットのイベントについて、HG/LG別の推定ピクセル分布を作る。"""
import argparse  # 共通の外部JSONを指定する。
import csv  # ピクセルごとのイベントと件数を保存する。
import json  # 設定と採用件数を保存する。
from pathlib import Path  # 入出力ファイルを扱う。
import numpy as np  # rawヒット数とピクセル座標を計算する。
from settings import DEFAULT_CONFIG, load_config, repo_path  # 他のrun06・07解析と同じ外部JSONを使う。
from pixel_polarization import rows, pixel_positions, fitted_axis, vertex_table, select_events  # 偏光解析と同じ位置較正・距離判定を使う。
from plot_hit_pixel_maps import draw_pixel_map  # 1/2ヒットで共通の投影とビーム印を使う。


def main():  # 1ヒットイベントのピクセル分布を作る。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行時の指定を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 入力キャッシュ・較正・保存先をJSONから指定する。
    args = parser.parse_args()  # 指定された設定を読む。
    config = load_config(args.config)  # 較正値と対象ランを読む。
    options = config["single_hit_pixels"]  # この解析専用の任意の距離条件を取り出す。
    prefix = repo_path(config["outputs"]["single_hit_pixel_prefix"])  # 2gainのCSV・JSON・図の共通保存名を読む。
    prefix.parent.mkdir(parents=True, exist_ok=True)  # 保存先を作る。
    geometry = {(int(row["board"]), int(row["ch"])): row for row in rows(repo_path(config["inputs"]["scintillator_geometry"]))}  # ピクセルを検出器座標へ変換する配置を読む。
    calibration_dir = repo_path(config["inputs"]["ha_position_calibration_dir"])  # HG位置較正の保存先を読む。
    point, direction = np.zeros(3), np.array([0., 0., 1.])  # 距離選別を使わない場合の初期値を置く。
    vertices = {}  # 既定ではTPCの散乱点を要求しない。
    if options["beam_vertex_cut_enabled"]:  # 1ヒット分布でも任意でビーム直線条件を使えるようにする。
        point, direction = fitted_axis(repo_path(config["inputs"]["polarization_beam_axis"]))  # 保存されたビーム直線を読む。
        vertices = vertex_table(repo_path(config["inputs"]["polarization_vertices"]))  # 同じトリガーの散乱点を読む。
    results = {gain: {"runs": [], "positions": [], "pixels": []} for gain in ("hg", "lg")}  # 同じrawイベントを2gainで比較する。
    for run in config["runs"]:  # ランごとにトリガーを照合する。
        cache = repo_path(config["inputs"]["single_hit_cache_dir"]) / f"{run}.npz"  # 保存済みの全rawヒットを入力にする。
        with np.load(cache) as stored:  # キャッシュを読み終えたら閉じる。
            metadata = json.loads(str(stored["metadata"]))  # キャッシュの対象範囲を読む。
            if metadata["max_files_per_board"] != 0:  # 試験用の一部ファイルを全ランの結果と混同しない。
                raise ValueError(f"Full-run cache required: {cache}")  # 全rawを読んだキャッシュを要求する。
            events, first, count = np.unique(stored["event"], return_index=True, return_counts=True)  # ピクセルを落とす前のrawヒット数を数える。
            selected = count == 1  # ちょうど1フラグのイベントだけを選ぶ。
            cut_options = dict(options, require_exactly_two_raw_hits=False, required_raw_hits=1)  # 共通距離判定では2ヒット条件を使わない。
            veto = np.load(repo_path(config["outputs"]["veto_dir"]) / f"z0_track_veto_{run}.npz")["event"] if options.get("z0_track_veto_enabled", False) else None  # Z側除外を使うときだけraw-TPCの判定を読む。
            allowed, distances, cutflow = select_events(events, count, vertices.get(run, np.zeros((0, 4))), point, direction, cut_options, veto, return_cutflow=True)  # 任意のビーム条件を適用する。
            indices = first[selected & allowed]  # rawで1ヒットのイベントの唯一のヒットを取り出す。
            hits = {key: stored[key][indices] for key in stored.files if key != "metadata"}  # 必要な1ヒット分だけをメモリへ移す。
        distance = distances[selected & allowed]  # 各1ヒットイベントの距離を出力へ対応づける。
        mapped = {gain: pixel_positions(hits, gain, geometry, calibration_dir) for gain in results}  # HGとLGで独立にピクセル位置を求める。
        common = mapped["hg"]["valid"] & mapped["lg"]["valid"]  # 公平な比較用に両gainで位置が有効なイベントを求める。
        for gain, product in mapped.items():  # 両gainを別々に保存する。
            valid = product["valid"].copy()  # このgainで非飽和・8×8内のヒットを選ぶ。
            if options["common_events_only"]:  # 同じイベントだけで比較したい場合の切替を設ける。
                valid &= common  # HG/LG両方で有効なイベントに揃える。
            summary = {"run": run, "raw_single_hit_events": int(selected.sum()), "events_after_optional_beam_cut": len(indices), "valid_pixel_events": int(valid.sum()), "cutflow": cutflow}  # raw選別とピクセル選別の件数を区別する。
            results[gain]["runs"].append(summary)  # ラン別の集計を残す。
            positions = product["position"][valid]  # 有効なピクセル中心を取り出す。
            local = np.column_stack((hits["board"][valid], hits["ch"][valid], product["pixel"][valid]))  # ボード・ch・n・mをまとめる。
            results[gain]["positions"].append(positions)  # 3方向の位置分布用に保存する。
            results[gain]["pixels"].append(local)  # ブロック内のピクセル件数を後で集計する。
            path = Path(str(prefix) + f"_{gain}_events.csv")  # gainごとのイベント一覧を決める。
            mode = "w" if run == config["runs"][0] else "a"  # 最初のランで新規作成し、後続を追記する。
            with path.open(mode) as stream:  # 出力CSVを開く。
                writer = csv.writer(stream)  # CSVの列を書き出す。
                if mode == "w":  # 最初だけ見出しを付ける。
                    writer.writerow(["run", "epoch", "trig", "raw_hits", "board", "ch", "output_ch", "n", "m", "x_mm", "y_mm", "z_mm", "hg", "lg", "saturated_hg", "saturated_lg", "beam_distance_mm"])  # ピクセルと元ヒットを追跡できる列を保存する。
                for index in np.flatnonzero(valid):  # 採用した1ヒットイベントを書き出す。
                    writer.writerow([run, int(hits["epoch"][index]), int(hits["trig"][index]), 1, int(hits["board"][index]), int(hits["ch"][index]), 5-int(hits["ch"][index]), *product["pixel"][index].tolist(), *product["position"][index].tolist(), float(hits["hg"][index]), float(hits["lg"][index]), int(hits["satu_h"][index]), int(hits["satu_l"][index]), float(distance[index]) if np.isfinite(distance[index]) else ""])  # raw数1を明示し、位置と飽和状態を保存する。
            print(f"{run} {gain.upper()}: {valid.sum():,} pixel events / {selected.sum():,} raw single-hit events", flush=True)  # 進捗と採用率を確認する。
    combined = {}  # 両gainの比較表を作る。
    for gain, result in results.items():  # 各gainの位置分布を描く。
        positions = np.vstack(result["positions"])  # 全ランのピクセル位置を合算する。
        pixel_rows = np.vstack(result["pixels"])  # 全ランのボード・ch・n・mを合算する。
        unique_pixels, counts = np.unique(pixel_rows, axis=0, return_counts=True)  # 同じブロック内ピクセルの件数を求める。
        with Path(str(prefix) + f"_{gain}_pixel_counts.csv").open("w") as stream:  # ピクセル別件数を保存する。
            writer = csv.writer(stream)  # CSV出力を用意する。
            writer.writerow(["board", "ch", "output_ch", "n", "m", "events"])  # 旧出力chも併記する。
            for (board, channel, n, m), count in zip(unique_pixels, counts):  # ピクセルごとのヒット数を並べる。
                writer.writerow([int(board), int(channel), 5-int(channel), int(n), int(m), int(count)])  # 1ヒットなので件数はイベント数と同じになる。
        title = f"{gain.upper()}: exactly one raw HA hit; N={len(positions):,}\nHG position calibration (LG provisional); beam cut={options['beam_vertex_cut_enabled']}; common events={options['common_events_only']}; z0 veto={options.get('z0_track_veto_enabled', False)}"  # 位置選別と表示条件を示す。
        draw_pixel_map(positions, config, title, "single-hit events / projected pixel", Path(str(prefix)+f"_{gain}.png"))  # ビーム軸付きの3投影図を保存する。
        record = {"gain": gain, "selection": "Exactly one raw DataCount hit before pixel cuts", "options": options, "runs": result["runs"], "valid_pixel_events": len(positions), "calibration": "HG position calibration; LG is provisional", "caveat": "A single raw channel hit does not establish single-pixel energy deposition."}  # 単一chと単一ピクセルの仮定を区別する。
        Path(str(prefix)+f"_{gain}.json").write_text(json.dumps(record, indent=2)+"\n")  # 条件と採用件数を保存する。
        combined[gain] = record  # HG/LG比較用にもまとめる。
    Path(str(prefix)+"_comparison.json").write_text(json.dumps(combined, indent=2)+"\n")  # 両版を一覧で確認できるようにする。


if __name__ == "__main__":  # importだけでは解析しない。
    main()  # コマンドとして呼ばれたときだけ実行する。
