#!/usr/bin/env python3
"""同じrawイベントから、カット前→ビーム軸→Z側→raw1/2の累積条件を段階別に描く。"""
from cutflow_report import write_cutflow  # 全run合計と日本語の説明を先頭に保存する。
import argparse  # 共通JSONを指定する。
import json  # 段階ごとの分母と採用数を保存する。
from pathlib import Path  # 出力フォルダを分ける。
import numpy as np  # イベントの対応とピクセル位置を計算する。
from settings import DEFAULT_CONFIG, load_config, repo_path  # 全解析と同じ設定を使う。
from pixel_polarization import rows, vertex_table, fitted_axis, select_events, pixel_positions, pair_indices, project_pairs, legacy  # 既存の定義を再使用する。
from plot_hit_pixel_maps import draw_pixel_map  # HG/LG共通の3投影とビーム印を使う。


def main():  # 累積カットの各段階を独立したフォルダへ保存する。
    parser = argparse.ArgumentParser(description=__doc__)  # 引数を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 入力・軸・距離・保存先を外部JSONから読む。
    args = parser.parse_args()  # 指定を読む。
    config = load_config(args.config)  # 共通設定を読む。
    root = repo_path(config["outputs"]["stage_pixel_map_root"])  # 条件別結果の親を指定する。
    point, direction = fitted_axis(repo_path(config["inputs"]["polarization_beam_axis"]))  # 第1段に使うビーム軸を読む。
    vertices = vertex_table(repo_path(config["inputs"]["vertices_all"]))  # 全散乱点を読み、Z側条件を先にかけない。
    geometry = {(int(row["board"]), int(row["ch"])): row for row in rows(repo_path(config["inputs"]["scintillator_geometry"]))}  # ピクセルの配置を読む。
    calibration = repo_path(config["inputs"]["ha_position_calibration_dir"])  # HG位置較正を読む。
    paths = {  # 累積条件の保存先を明示する。
        "00_before_cuts": root / "without_z0_cut/stages/00_before_cuts",  # 全rawヒットイベント。
        "01_beam_axis": root / "without_z0_cut/stages/01_beam_axis",  # ビーム軸だけで選別。
        "02_after_z0": root / "with_z0_cut/stages/02_beam_axis_then_z0",  # その後にイベント全体のZ側除外を適用。
        "03_raw1": root / "with_z0_cut/stages/04_single_hit_quality",  # さらにraw1ヒットを要求。
        "03_raw2": root / "with_z0_cut/stages/04_two_hits_quality",  # さらにraw2ヒットかつ有効ペアを要求。
    }  # 2種類のraw数は同じ第3段から分岐する。
    collected = {stage: {gain: [] for gain in ("hg", "lg")} for stage in paths}  # 描画点群を段階ごとに集める。
    summaries = []  # 各ランのカットフローを残す。
    options = dict(config["polarization"], require_exactly_two_raw_hits=False, beam_vertex_cut_enabled=True, z0_track_veto_enabled=False)  # まずビーム条件だけを作る。
    e1, e2 = legacy.axes(direction, np.asarray(options["reference_axis"], dtype=float))  # 最終2ヒットの方位角が定義可能かも既存と同じ条件にする。
    for run in config["runs"]:  # ランを混ぜずに判定する。
        cache = repo_path(config["outputs"]["polarization_hit_cache_dir"]) / f"{run}.npz"  # 全raw重心キャッシュを指定する。
        with np.load(cache) as saved:  # SSDのrawを読み直さない。
            meta = json.loads(str(saved["metadata"]))  # キャッシュの対象範囲を確認する。
            if meta["max_files_per_board"] != 0:  # 一部ファイルを全データとして扱わない。
                raise ValueError(f"Full-run cache required: {cache}")  # 入力範囲不一致を明示する。
            hits = {key: saved[key] for key in saved.files if key != "metadata"}  # rawヒットを読む。
        order = np.argsort(hits["event"], kind="stable")  # 同じトリガーを連続させる。
        hits = {key: value[order] for key, value in hits.items()}  # 全列を同じ順番にする。
        events, starts, counts = np.unique(hits["event"], return_index=True, return_counts=True)  # rawヒット数を先に確定する。
        beam, distances = select_events(events, counts, vertices.get(run, np.zeros((0, 4))), point, direction, options)  # 第1段の距離条件だけを適用する。
        veto = np.load(repo_path(config["outputs"]["veto_dir"]) / f"z0_track_veto_{run}.npz")["event"]  # 全TPCヒットから決めたZ側除外集合を読む。
        after_z0 = beam & ~np.isin(events, veto)  # ビーム条件を通ったイベントに第2段を適用する。
        masks = {"00_before_cuts": np.ones(len(events), bool), "01_beam_axis": beam, "02_after_z0": after_z0, "03_raw1": after_z0 & (counts == 1), "03_raw2": after_z0 & (counts == 2)}  # 第3段でraw1/2を選ぶ。
        assert np.all(~after_z0 | beam)  # Z側除外後が第1段の部分集合か確認する。
        run_summary = {"run": run, "events": {stage: int(mask.sum()) for stage, mask in masks.items()}, "displayed_hit_uses": {}}  # ピクセル有効性より前の分母を記録する。
        maps = {gain: pixel_positions(hits, gain, geometry, calibration, reject_saturation=False) for gain in ("hg", "lg")}  # 全段階で同じピクセル推定を使う。
        for gain, product in maps.items():  # gainによる有効性は独立に扱う。
            raw_positions = product["position"].copy()  # 品質選別前の図では全rawフラグを保持する。
            for index in np.flatnonzero(~product["valid"]):  # ピクセル位置を求められないヒットを調べる。
                row = geometry[(int(hits["board"][index]), int(hits["ch"][index]))]  # 既知のch配置を読む。
                raw_positions[index] = [float(row[key]) for key in ("x0", "y0", "z0")]  # 未決定のヒットだけch中心へ置く。
            quality_valid = product["valid"] & (hits["satu_h" if gain == "hg" else "satu_l"] == 0)  # stage04で初めてgain別の飽和も除く。
            run_summary["displayed_hit_uses"][gain] = {}  # gain別の採用数を用意する。
            for stage, event_mask in masks.items():  # 各累積条件の位置を取り出す。
                if stage == "03_raw2":  # 最終2ヒット図は両ピクセル有効なペアとする。
                    pairs = pair_indices(starts, counts, event_mask, quality_valid, True)  # raw数と有効ピクセル数を混同しない。
                    pairs, _ = project_pairs(pairs, product["position"], e1, e2, options["min_pair_distance_mm"])  # 偏光解析と同じ投影長条件を使う。
                    indices = pairs.ravel()  # 採用ペアの2点を積算する。
                else:  # 第0〜2段は全raw多重度、第3段raw1は1点を使う。
                    use = np.repeat(event_mask, counts) & (quality_valid if stage == "03_raw1" else np.ones(len(hits["event"]), bool))  # イベント条件の後で有効ピクセルだけを描く。
                    indices = np.flatnonzero(use)  # 描画するrawヒット番号を得る。
                collected[stage][gain].append((product["position"] if stage.startswith("03_") else raw_positions)[indices])  # 同じ段階の全ランを合算する。
                run_summary["displayed_hit_uses"][gain][stage] = len(indices)  # 表示された点数を記録する。
        summaries.append(run_summary)  # ラン別のカットフローを残す。
        print(run,run_summary["events"],flush=True)  # 順番どおり残ったイベント数を表示する。
    labels = {"00_before_cuts": "Before event cuts; all raw multiplicities", "01_beam_axis": "After beam-axis distance cut; all raw multiplicities", "02_after_z0": "After beam-axis cut THEN Z-side veto; all raw multiplicities", "03_raw1": "Beam-axis THEN Z-side veto THEN exactly one raw HA hit", "03_raw2": "Beam-axis THEN Z-side veto THEN exactly two raw HA hits; valid pairs"}  # 図から段階を読み取れるようにする。
    color_max = {gain: None for gain in ("hg", "lg")}  # 同じgainの段階間で色範囲を固定する。
    for stage, directory in paths.items():  # 各段階を別フォルダに保存する。
        directory.mkdir(parents=True, exist_ok=True)  # フォルダを作る。
        metadata = {"stage": stage, "selection": labels[stage], "cut_order": ["beam_axis", "z0_track_veto", "raw_multiplicity"], "beam_max_distance_mm": options["beam_max_distance_mm"], "z_veto_mm": config["z_veto"], "runs": summaries}  # 順序と分母を保存する。
        (directory / "cutflow.json").write_text(json.dumps(metadata, indent=2)+"\n")  # 各フォルダだけでも条件を確認できるようにする。
        for gain in ("hg", "lg"):  # HG/LGを同じ条件で描く。
            points = np.vstack(collected[stage][gain])  # 全ランの位置をまとめる。
            title = f"{gain.upper()}: {labels[stage]}; N={len(points):,} hit uses\nHG calibration (LG provisional); unresolved raw hits at channel centres before quality cuts; beam radius={options['beam_max_distance_mm']} mm"  # 実際の表示点数を示す。
            result = draw_pixel_map(points, config, title, ("accepted hit uses / projected pixel" if stage.startswith("03_") else "raw hit flags / projected position"), directory/f"{gain}.png", color_max=color_max[gain])  # 各段階で同じビーム印と投影を使う。
            if stage == "00_before_cuts":  # カット前の色上限を基準にする。
                color_max[gain] = result["color_max"]  # 後続段階に同じ上限を渡す。
            assert result["projected_totals"] == [len(points)]*3  # 投影による件数変化がないことを確認する。
            print('Saved',directory/f"{gain}.png",flush=True)  # 保存先を知らせる。
    raw_paths = dict(paths, **{"03_raw1": root / "with_z0_cut/stages/03_raw_single_hit", "03_raw2": root / "with_z0_cut/stages/03_raw_two_hits"})  # stage03は全rawフラグ、stage04は品質選別後へ分離する。
    write_cutflow(root / "stage_cutflow.json", {"cut_order": ["beam_axis", "z0_track_veto", "raw_multiplicity", "gain_specific_quality"], "paths": {k:str(v) for k,v in raw_paths.items()}, "quality_paths": {k:str(v) for k,v in paths.items() if k.startswith("03_")}, "runs": summaries}, "scint")  # 全run合計とラン別件数を一緒に保存する。


if __name__ == "__main__":  # importだけでは処理しない。
    main()  # コマンドとして呼ばれたときだけ実行する。
