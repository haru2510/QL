#!/usr/bin/env python3
"""シンチのHG/LGピクセル位置による方位角解析。入力・出力・選別条件は共通JSONで指定する。"""
import argparse  # 外部JSONとキャッシュからの実行を指定する。
import contextlib  # 詳細なraw読み込みログを別ファイルへ保存する。
import csv  # ピクセル較正・散乱点・ペア一覧を読み書きする。
import itertools  # 2ヒット制限なしの場合に全組み合わせを作る。
import json  # 条件と集計を保存する。
import sys  # 既存のHAデコーダを読み込む。
from pathlib import Path  # 入出力パスを扱う。
import numpy as np  # 重心補正・距離・方位角を計算する。
from scipy.interpolate import RBFInterpolator  # HG用位置較正をthin-plate splineで適用する。
import matplotlib  # ファイル出力用の描画環境を設定する。
matplotlib.use("Agg")  # GUIを開かない。
import matplotlib.pyplot as plt  # HG/LG別の分布を描く。
from settings import DEFAULT_CONFIG, REPO, load_config, repo_path  # 共通外部JSONを使う。
from plot_hit_pixel_maps import draw_pixel_map, beam_reference  # 採用ペアのピクセル位置を3投影でも描く。
sys.path.insert(0, str(REPO / "ETCC_analysis"))  # 既存のHA処理を検索する。
import pair_polarization as legacy  # rawヒット読み取りと変調fitを再使用する。


def rows(path):  # コメント付きCSVを読む。
    with Path(path).open() as source:  # 読み終えたファイルを閉じる。
        return list(csv.DictReader(line for line in source if line.strip() and not line.startswith("#")))  # コメントを除いて列名で読む。


def fitted_axis(path):  # 他コードの局所ピークfitから3次元直線を作る。
    saved = json.loads(Path(path).read_text())  # 保存済みの係数を読む。
    z = float(saved.get("line_reference_z_mm", 155))  # 軸位置の基準Z [mm]を読む。
    point = np.array([saved["x_line"]["at_z155_mm"], saved["y_line"]["at_z155_mm"], z])  # 直線上の1点を得る。
    direction = np.array([saved["x_line"]["slope_mm_per_mm"], saved["y_line"]["slope_mm_per_mm"], 1.0])  # Zが増える向きの方向ベクトルを作る。
    direction /= np.linalg.norm(direction)  # 距離計算用の単位ベクトルにする。
    return point, direction  # ビーム直線を返す。


def vertex_table(path):  # 散乱点をラン・イベント番号で照合できるように読む。
    table = {}  # ランごとに散乱点をまとめる。
    with Path(path).open() as source:  # 大きいCSVを逐次読む。
        next(source)  # 座標の説明行を読み飛ばす。
        for row in csv.DictReader(source):  # 保存された全頂点を読む。
            event = (int(row["epoch"]) << 32) + int(row["trig"])  # HAと同じイベント識別子にする。
            table.setdefault(row["run"], []).append((event, float(row["x"]), float(row["y"]), float(row["z"])))  # 同じトリガーの散乱点を保存する。
    return {run: np.array(sorted(values), dtype=float) for run, values in table.items()}  # イベント順に並べ検索に使う。


def select_events(events, counts, vertices, point, direction, options, veto_events=None, return_cutflow=False):  # ビーム軸→Z側→rawヒット数の順に選別する。
    selected = np.ones(len(events), dtype=bool)  # rawヒットがある全イベントから始める。
    flow = {"initial_events": len(events)}  # カット前の分母を保存する。
    distances = np.full(len(events), np.nan)  # 頂点なしの距離は未定にする。
    if options["beam_vertex_cut_enabled"]:  # 最初に散乱点とビーム直線の距離を判定する。
        valid = np.zeros(len(events), dtype=bool)  # 頂点がないイベントは通さない。
        if len(vertices):  # 同じランの頂点があるか確認する。
            keys = vertices[:, 0].astype(np.int64)  # 頂点イベント番号を整数にする。
            indices = np.searchsorted(keys, events)  # 同じトリガーを検索する。
            matched = np.flatnonzero(indices < len(keys))  # 配列範囲内だけを調べる。
            matched = matched[keys[indices[matched]] == events[matched]]  # 番号が一致する頂点だけを使う。
            delta = vertices[indices[matched], 1:4] - point  # ビーム直線上の基準点からの変位を求める。
            distances[matched] = np.linalg.norm(np.cross(delta, direction), axis=1)  # 直線との最短3D距離を求める。
            valid[matched] = distances[matched] <= options["beam_max_distance_mm"]  # 許容距離以内なら採用する。
        selected &= valid  # 第1段の選別を適用する。
    flow["after_beam_axis"] = int(selected.sum())  # ビーム軸条件だけの残存数を記録する。
    if options.get("z0_track_veto_enabled", False):  # 第2段としてイベント全体のZ側除外を切り替える。
        if veto_events is None:  # 除外リストなしでカットしたことにしない。
            raise ValueError("Z-side veto enabled but no raw-TPC veto list supplied")  # 入力不足を明示する。
        selected &= ~np.isin(events, veto_events)  # どこかにZ<=15 mmのTPCヒットがあるイベントを除外する。
    flow["after_z0_track_veto"] = int(selected.sum())  # 第2段後の数を記録する。
    required = options.get("required_raw_hits", 2 if options["require_exactly_two_raw_hits"] else None)  # 最後にrawヒット数の条件を選ぶ。
    if required is not None:  # 1/2ヒット条件を必要な場合だけ適用する。
        selected &= counts == required  # ピクセル無効ヒットを落とす前の数を使う。
    flow["after_raw_multiplicity"] = int(selected.sum())  # rawヒット数条件後の数を記録する。
    return (selected, distances, flow) if return_cutflow else (selected, distances)  # 従来の2戻り値にも対応する。


def pixel_positions(hits, gain, geometry, calibration_dir, reject_saturation=True):  # 各rawヒットをHGまたはLGのピクセルへ割り当てる。
    size = len(hits["event"])  # rawフラグで数えたヒット数を取得する。
    positions = np.full((size, 3), np.nan)  # 無効な位置をNaNで表す。
    pixels = np.full((size, 2), -1, dtype=np.int16)  # 無効ピクセルを-1で表す。
    corrected = np.full((size, 2), np.nan)  # 補正後の連続ピクセル座標も保存する。
    valid = np.zeros(size, dtype=bool)  # ピクセルを特定できたヒットだけを有効にする。
    channel_keys = np.column_stack((hits["board"], hits["ch"]))  # ボード・chの組を作る。
    for board, channel in np.unique(channel_keys, axis=0):  # 同じ較正表を使うヒットをまとめて計算する。
        indices = np.flatnonzero((hits["board"] == board) & (hits["ch"] == channel))  # このブロックの全rawヒットを取り出す。
        key = (int(board), int(channel))  # 配置表の検索キーを作る。
        path = calibration_dir / f"position_ha{board:02d}.csv"  # HG用位置較正表を選ぶ。
        if key not in geometry or not path.exists():  # 配置か較正がなければピクセルを捏造しない。
            continue  # rawヒット数には残し、位置だけ無効にする。
        calibration = [row for row in rows(path) if int(row["ch"]) == channel]  # 現行ADC番号順のchだけを選ぶ。
        if len(calibration) != 64:  # 64ピクセルの較正が揃っているか確認する。
            continue  # 不完全な較正を位置推定に使わない。
        observed = np.array([[float(row["x"]), float(row["y"])] for row in calibration])  # 重心平面のHGピーク位置を読む。
        ideal = np.array([[int(row["n"]), int(row["m"])] for row in calibration])  # 1〜8のピクセル番号を読む。
        warp = RBFInterpolator(observed, ideal, kernel="thin_plate_spline", smoothing=0)  # EricaNextと同じ位置補間を使う。
        candidates = indices[hits[gain][indices] > 0]  # 重心の計算に使える正積分を要求し、求められない位置は呼出側で処理する。
        if reject_saturation:  # 品質選別後の既定動作ではこのgainの飽和を除く。
            candidates = candidates[hits["satu_h" if gain == "hg" else "satu_l"][candidates] == 0]  # フラグのみの図では飽和を理由に位置推定を止めない。
        centroid = np.column_stack((hits["x_" + gain][candidates], hits["y_" + gain][candidates]))  # 指定gainの重心を得る。
        finite = np.all(np.isfinite(centroid), axis=1)  # NaN重心を較正へ入れない。
        candidates, centroid = candidates[finite], centroid[finite]  # 有限な重心だけに絞る。
        for start in range(0, len(candidates), 20000):  # 補間の一時メモリを小さく抑える。
            subset = candidates[start:start + 20000]  # この塊のヒット番号を得る。
            corrected[subset] = warp(centroid[start:start + 20000])  # HG較正に指定gainの重心を入力する。
        finite = np.all(np.isfinite(corrected[candidates]), axis=1)  # 補正後も有限であることを確認する。
        candidates = candidates[finite]  # 非有限な補正値を除外する。
        inside = np.all((corrected[candidates] >= 0.5) & (corrected[candidates] < 8.5), axis=1)  # 整数化前に8×8領域内か確認する。
        candidates = candidates[inside]  # 範囲外を端のピクセルへ丸め込まない。
        pixels[candidates] = np.floor(corrected[candidates] + 0.5).astype(np.int16)  # n,mをEricaNextと同じ丸め規則で決める。
        row = geometry[key]  # ブロックの向きとピッチを取り出す。
        center = np.array([float(row[name]) for name in ("x0", "y0", "z0")])  # ブロック中心 [mm]を得る。
        u = np.array([float(row[name]) for name in ("ux", "uy", "uz")])  # nが増える方向を得る。
        v = np.array([float(row[name]) for name in ("vx", "vy", "vz")])  # mが増える方向を得る。
        offset = (pixels[candidates, 0, None] - 4.5) * u + (pixels[candidates, 1, None] - 4.5) * v  # ピクセル中心からの変位をピッチ単位で求める。
        positions[candidates] = center + float(row["pitch"]) * offset  # 検出器座標 [mm]へ変換する。
        valid[candidates] = True  # ピクセルが確定できたヒットを採用する。
    return {"position": positions, "pixel": pixels, "corrected": corrected, "valid": valid}  # 方位角と出力に必要な情報を返す。


def pair_indices(starts, counts, selected, valid, exactly_two):  # rawイベント内で有効ピクセルのペアを作る。
    if exactly_two:  # 既定の2ヒット条件では配列でまとめて処理する。
        first = starts[selected & (counts == 2)]  # rawでちょうど2ヒットだったイベントの先頭を選ぶ。
        pairs = np.column_stack((first, first + 1))  # 同じイベントの2ヒットを組にする。
        return pairs[np.all(valid[pairs], axis=1)]  # ピクセルが両方有効な組だけを返す。
    groups = []  # 多重ヒットではイベント内の全ペアを保存する。
    for start, count in zip(starts[selected], counts[selected]):  # イベントをまたがないよう処理する。
        indices = np.flatnonzero(valid[start:start + count]) + start  # このイベントの有効ピクセルを取り出す。
        groups.extend(itertools.combinations(indices, 2))  # 全ての無順序ペアを作る。
    return np.array(groups, dtype=np.int64).reshape(-1, 2)  # 0件でも2列の配列を返す。


def project_pairs(pairs, positions, e1, e2, min_distance):  # ペアの方位角を求め、投影長0を除外する。
    delta = positions[pairs[:, 1]] - positions[pairs[:, 0]]  # 2ピクセルを結ぶベクトルを作る。
    u, v = delta @ e1, delta @ e2  # ビーム軸に垂直な面へ射影する。
    keep = (np.linalg.norm(delta, axis=1) >= min_distance) & (u * u + v * v > 1e-12)  # 距離条件と方位角の定義可能性を確認する。
    return pairs[keep], np.degrees(np.arctan2(v[keep], u[keep])) % 180  # 線の方向を0〜180度に畳む。


def pair_weights(events, weighting):  # 多重ヒットイベントの重みを定義する。
    if weighting == "per_pair":  # 各ペアを1票にする指定を受け付ける。
        return np.ones(len(events))  # 多重ヒットほど多く数える方式。
    if weighting != "per_event":  # 未知の重み付けを黙って適用しない。
        raise ValueError("pair_weighting must be per_event or per_pair")  # 設定ミスを明示する。
    _, inverse, counts = np.unique(events, return_inverse=True, return_counts=True)  # 同じイベントの採用ペア数を数える。
    return 1 / counts[inverse]  # 1イベントの重みの合計が1になるようにする。


def main():  # 同じrawからHG版とLG版を作る。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行方法の説明を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 他の解析と同じJSONを受け取る。
    parser.add_argument("--from-cache", action="store_true")  # 散乱点条件をSSDなしで切り替えられるようにする。
    parser.add_argument("--max-files", type=int)  # 動作確認用のボードごとのファイル制限を受け取る。
    parser.add_argument("--beam-cut", action=argparse.BooleanOptionalAction, default=None)  # ビーム直線条件を上書きできるようにする。
    parser.add_argument("--require-two-hits", action=argparse.BooleanOptionalAction, default=None)  # rawの2ヒット条件を上書きできるようにする。
    args = parser.parse_args()  # コマンド引数を読む。
    config = load_config(args.config)  # raw入力・較正・結果パスを外部JSONから読む。
    options = dict(config["polarization"])  # この解析固有の切替条件を取り出す。
    if args.max_files is not None:  # 試験用の制限が明示されたか確認する。
        options["max_files_per_board"] = args.max_files  # JSONよりコマンド指定を優先する。
    if args.beam_cut is not None:  # ビーム条件が上書きされたか確認する。
        options["beam_vertex_cut_enabled"] = args.beam_cut  # 個別の有効・無効を反映する。
    if args.require_two_hits is not None:  # raw2ヒット条件が上書きされたか確認する。
        options["require_exactly_two_raw_hits"] = args.require_two_hits  # 独立した切替を反映する。
    if options["beam_max_distance_mm"] < 0:  # 負の距離を許さない。
        raise ValueError("beam_max_distance_mm must be nonnegative")  # 設定ミスを明示する。
    prefix = repo_path(config["outputs"]["pixel_polarization_prefix"])  # HG/LG図・CSV・JSONの共通保存名を読む。
    prefix.parent.mkdir(parents=True, exist_ok=True)  # 保存先を用意する。
    cache_dir = repo_path(config["outputs"]["polarization_hit_cache_dir"])  # rawヒットのキャッシュ保存先を読む。
    cache_dir.mkdir(parents=True, exist_ok=True)  # キャッシュの親を作る。
    source = options["beam_direction_source"]  # 方位角の軸を保存fitか任意指定から選ぶ。
    if source not in ("fitted_line", "configured"):  # 誤った文字列を黙って解釈しない。
        raise ValueError("beam_direction_source must be fitted_line or configured")  # 選べる値を示す。
    cut_point, cut_direction = np.zeros(3), np.array([0., 0., 1.])  # 距離条件も保存軸も使わない場合の初期値を置く。
    if options["beam_vertex_cut_enabled"] or source == "fitted_line":  # 必要な場合だけ保存済みビーム軸を読む。
        cut_point, cut_direction = fitted_axis(repo_path(config["inputs"]["polarization_beam_axis"]))  # 他コードの局所ピークfitを使用する。
    direction = cut_direction.copy() if source == "fitted_line" else np.asarray(options["beam_direction"], dtype=float)  # 方位角を定義する向きを選ぶ。
    if not np.all(np.isfinite(direction)) or np.linalg.norm(direction) == 0:  # 不正な方向からNaN角度を作らない。
        raise ValueError("beam_direction must be finite and nonzero")  # 入力ミスを明示する。
    direction /= np.linalg.norm(direction)  # 単位方向に直す。
    e1, e2 = legacy.axes(direction, np.asarray(options["reference_axis"], dtype=float))  # 方位角の基準面を作る。
    vertices = vertex_table(repo_path(config["inputs"]["polarization_vertices"])) if options["beam_vertex_cut_enabled"] else {}  # ビーム条件を切った場合は頂点CSVを読まない。
    geometry = {(int(row["board"]), int(row["ch"])): row for row in rows(repo_path(config["inputs"]["scintillator_geometry"]))}  # ピクセルの検出器座標を定義する。
    calibration_dir = repo_path(config["inputs"]["ha_position_calibration_dir"])  # HG位置較正の保存先を読む。
    edges = np.linspace(0, 180, options["phi_bins"] + 1)  # 方位角のビン境界を作る。
    accum = {gain: {"counts": np.zeros(len(edges) - 1), "sumw2": np.zeros(len(edges) - 1), "runs": [], "xy": []} for gain in ("hg", "lg")}  # HG/LGを別々に集計する。
    header = ["run", "epoch", "trig", "raw_hits", "beam_distance_mm", "board1", "ch1", "n1", "m1", "x1_mm", "y1_mm", "z1_mm", "board2", "ch2", "n2", "m2", "x2_mm", "y2_mm", "z2_mm", "phi_deg", "weight"]  # どのピクセルを使ったか追跡できる列を用意する。
    with contextlib.ExitStack() as stack:  # CSVとログを正常終了・例外の両方で閉じる。
        writers = {}  # gainごとのCSV出力先を作る。
        for gain in accum:  # HGとLGを別ファイルで保存する。
            handle = stack.enter_context(Path(str(prefix) + f"_{gain}_pairs.csv").open("w"))  # ピクセル付きペア一覧を開く。
            writers[gain] = csv.writer(handle)  # CSVへの書き出しを用意する。
            writers[gain].writerow(header)  # 列名を保存する。
        for run in config["runs"]:  # ランをまたいでトリガーを混ぜない。
            cache = cache_dir / f"{run}.npz"  # 再選別用の全rawヒットキャッシュを指定する。
            if args.from_cache:  # SSDなしで条件だけ変える場合を扱う。
                stored = np.load(cache)  # HG/LG重心付きのrawヒットを読む。
                saved = json.loads(str(stored["metadata"]))  # 読み込み範囲の出典を得る。
                if saved["max_files_per_board"] != options["max_files_per_board"]:  # 試験用キャッシュを全データと取り違えない。
                    raise ValueError("Cache file limit differs from configured max_files_per_board")  # 範囲不一致なら止める。
                hits = {key: stored[key] for key in stored.files if key != "metadata"}  # 解析列を取り出す。
            else:  # rawから両gainを同時に読む。
                directories = [str(repo_path(config["data_root"]) / storage / run) for storage in config["storages"]]  # 入力storageを指定する。
                boards = legacy.collect_files(directories, options["include_tmp"], options["max_files_per_board"])  # フラグを数える前にHG/LGの選別をしない。
                if not boards:  # 空の入力を解析成功にしない。
                    raise FileNotFoundError(f"No HA raw files for {run}")  # 入力先を確認できるエラーを返す。
                print(f"Reading {run}: {sum(map(len, boards.values()))} files", flush=True)  # 進捗を知らせる。
                with Path(str(prefix) + f"_{run}_raw.log").open("w") as log, contextlib.redirect_stdout(log):  # ファイル単位の詳細をログへ分ける。
                    hits, frames, resets = legacy.read_hits(boards, include_centroids=True)  # rawフラグとHG/LG重心をまとめて読む。
                saved = {"run": run, "max_files_per_board": options["max_files_per_board"], "frames": frames, "resets": resets, "files": boards}  # キャッシュの対象範囲を記録する。
                np.savez_compressed(cache, **hits, metadata=json.dumps(saved))  # 条件変更時のraw再読込を省けるようにする。
            order = np.argsort(hits["event"], kind="stable")  # 同じイベントのヒットを並べる。
            hits = {key: value[order] for key, value in hits.items()}  # 全列に同じ並べ替えを適用する。
            events, starts, counts = np.unique(hits["event"], return_index=True, return_counts=True)  # rawヒット数を先に確定する。
            veto = np.load(repo_path(config["outputs"]["veto_dir"]) / f"z0_track_veto_{run}.npz")["event"] if options.get("z0_track_veto_enabled", False) else None  # Z側条件を使う場合だけ全TPCヒットの除外集合を読む。
            selected, distances, cutflow = select_events(events, counts, vertices.get(run, np.zeros((0, 4))), cut_point, cut_direction, options, veto, return_cutflow=True)  # 任意のイベント条件を適用する。
            maps = {gain: pixel_positions(hits, gain, geometry, calibration_dir) for gain in accum}  # 同じヒットをHGとLGで別々にピクセルへ割り当てる。
            products = {}  # 各gainの採用ペアと方位角を用意する。
            for gain, mapped in maps.items():  # ピクセルの有効性はgainごとに判定する。
                pairs = pair_indices(starts, counts, selected, mapped["valid"], options["require_exactly_two_raw_hits"])  # rawの2ヒット判定を位置有効性で置き換えない。
                products[gain] = project_pairs(pairs, mapped["position"], e1, e2, options["min_pair_distance_mm"])  # ピクセル中心の線から方位角を求める。
            if options["common_pairs_only"]:  # HG/LGを全く同じペアで比較する場合を扱う。
                common = np.intersect1d(products["hg"][0][:, 0] * len(order) + products["hg"][0][:, 1], products["lg"][0][:, 0] * len(order) + products["lg"][0][:, 1])  # 両gainで有効なヒット番号の組だけを残す。
                for gain, (pairs, phi) in products.items():  # gainごとの図を共通ペアへ揃える。
                    keep = np.isin(pairs[:, 0] * len(order) + pairs[:, 1], common)  # 両gainで同じrawヒットの組か確認する。
                    products[gain] = pairs[keep], phi[keep]  # 方位角の値は各gainで計算したものを使う。
            for gain, (pairs, phi) in products.items():  # HG/LGを独立に保存する。
                mapped = maps[gain]  # 対応するピクセルと座標を使う。
                pair_events = hits["event"][pairs[:, 0]]  # ペアが属するイベント番号を得る。
                weights = pair_weights(pair_events, options["pair_weighting"])  # 多重ヒットの重みを明示する。
                accum[gain]["counts"] += np.histogram(phi, edges, weights=weights)[0]  # 重み付き方位角分布を足す。
                accum[gain]["sumw2"] += np.histogram(phi, edges, weights=weights**2)[0]  # 表示用の重み二乗和を保存する。
                indices = np.searchsorted(events, pair_events)  # ペアからrawヒット数とビーム距離へ戻る。
                for k, ((i, j), angle, weight) in enumerate(zip(pairs, phi, weights)):  # ピクセルを明記したCSVを作る。
                    left = [int(hits["board"][i]), int(hits["ch"][i]), *mapped["pixel"][i].tolist(), *mapped["position"][i].tolist()]  # 1個目のピクセルと位置を出力する。
                    right = [int(hits["board"][j]), int(hits["ch"][j]), *mapped["pixel"][j].tolist(), *mapped["position"][j].tolist()]  # 2個目も同じ形式で出力する。
                    distance = distances[indices[k]]  # オプションの距離を監査できるようにする。
                    writers[gain].writerow([run, int(hits["epoch"][i]), int(hits["trig"][i]), int(counts[indices[k]]), float(distance) if np.isfinite(distance) else "", *left, *right, float(angle), float(weight)])  # イベントと両ピクセルを1行にまとめる。
                used = pairs.ravel()  # 採用ペアに使ったシンチ位置を集める。
                accum[gain]["xy"].append(mapped["position"][used])  # 位置分布を後で描くため保存する。
                information = {"run": run, "raw_hits": len(order), "raw_hit_events": len(events), "events_after_event_cuts": int(selected.sum()), "cutflow": cutflow, "valid_pixel_hits": int(mapped["valid"].sum()), "pairs": len(pairs), "pair_events": len(np.unique(pair_events)), "weighted_sum": float(weights.sum()), "raw_multiplicity": np.bincount(counts).tolist()}  # 条件ごとの減少を追跡する。
                accum[gain]["runs"].append(information)  # gain別の集計へ残す。
                print(f"{run} {gain.upper()}: {len(pairs):,} pairs / {information['pair_events']:,} events", flush=True)  # 違いを確認できる進捗を出す。
    centers = (edges[1:] + edges[:-1]) / 2  # fitに使うビン中心を求める。
    summary = {}  # 両gainの比較結果を保存する。
    for gain, result in accum.items():  # 各gainの図と数値を作る。
        errors = np.sqrt(result["sumw2"])  # 重み付きの表示用誤差を計算する。
        fit = legacy.fit_modulation(centers, result["counts"], errors)  # 既存のcos2phiモデルをfitする。
        record = {"gain": gain, "options": options, "config": str(args.config.resolve()), "beam_point_mm": cut_point.tolist(), "beam_direction": direction.tolist(), "phi_edges_deg": edges.tolist(), "weighted_counts": result["counts"].tolist(), "sumw2": result["sumw2"].tolist(), "runs": result["runs"], "fit": fit, "calibration": "HG position calibration; LG application is provisional", "caveat": "Raw azimuth modulation is not a calibrated polarization fraction. Pair correlations are not included in sumw2 errors."}  # 暫定較正と統計上の制約も残す。
        Path(str(prefix) + f"_{gain}.json").write_text(json.dumps(record, indent=2) + "\n")  # HG/LG別の集計を保存する。
        summary[gain] = record  # 両gainの比較用にもまとめる。
        fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)  # 方位角・ピクセル位置・rawヒット数を並べる。
        axes[0].errorbar(centers, result["counts"], yerr=errors, fmt="o", ms=4, color="black")  # 重み付き方位角分布を描く。
        if fit:  # 十分な方位角ビンがあるときだけfitを描く。
            x = np.linspace(0, 180, 361)  # なめらかなモデル線の角度を作る。
            axes[0].plot(x, fit["A"] * (1 + fit["mu"] * np.cos(np.radians(2 * (x - fit["phi0_deg"])))), "r-", label=f"mu={fit['mu']:.3f}; phi0={fit['phi0_deg']:.1f} deg")  # 生の変調fitを重ねる。
            axes[0].legend(fontsize=8)  # 係数を図に残す。
        axes[0].set(xlabel="azimuth [deg]", ylabel="weighted pairs / bin", title="Raw azimuth (no acceptance correction)", xlim=(0, 180))  # 偏光度と断定しない。
        points = np.vstack(result["xy"]) if result["xy"] else np.zeros((0, 3))  # 採用ピクセル位置を合算する。
        if len(points):  # ペアがある場合に位置図を描く。
            positions, count = np.unique(points, axis=0, return_counts=True)  # 同じピクセルの出現数を数える。
            image = axes[1].scatter(positions[:, 0], positions[:, 1], c=count, s=12, marker="s", cmap=config.get("pixel_map", {}).get("colormap", "turbo"))  # ブロック中心でなくピクセル中心を描く。
            fig.colorbar(image, ax=axes[1], label="hit uses in selected pairs")  # ペアへの使用回数であることを示す。
        axes[1].set(xlabel="X [mm]", ylabel="Y [mm]", title=f"{gain.upper()} assigned pixel centres")  # 使用したgainを示す。
        if config.get("beam_axis_overlay", {}).get("enabled", False):  # 方位角図の中央にも保存ビーム軸の通過点を示す。
            _, marker, _, _ = beam_reference(config)  # 底面シンチ平面でのビーム位置を得る。
            axes[1].scatter(marker[0], marker[1], marker="*", s=160, c="magenta", edgecolors="black", linewidths=0.7, alpha=config["beam_axis_overlay"].get("marker_alpha", 0.45), label=f"Beam at Z={marker[2]:.2f} mm", zorder=5)  # 半透明で下の分布を見えるようにする。
            axes[1].legend(fontsize=7)  # 星印の意味を表示する。
        axes[1].invert_yaxis()  # 既存図と同じ向きにする。
        multiplicity = np.zeros(max([len(row["raw_multiplicity"]) for row in result["runs"]] + [1]))  # 全rawヒット数の分布を用意する。
        for row in result["runs"]:  # ランごとのraw分布を足す。
            multiplicity[:len(row["raw_multiplicity"])] += row["raw_multiplicity"]  # ピクセルで落とす前の数を描く。
        axes[2].bar(np.arange(len(multiplicity)), multiplicity, color="0.4")  # 2ヒット以外も含む全分布を表示する。
        axes[2].set(xlabel="raw HA hit flags / event", ylabel="events", yscale="log", title="Raw multiplicity before pixel cuts")  # 有効ピクセル数と混同させない。
        fig.suptitle((config["data_label"] + "\n" if config.get("data_label") else "") + f"{gain.upper()} pixel polarization; raw2={options['require_exactly_two_raw_hits']}; vertex-line cut={options['beam_vertex_cut_enabled']}; z0 veto={options.get('z0_track_veto_enabled', False)}; R={options['beam_max_distance_mm']} mm\nHG position calibration; LG version is provisional; weighting={options['pair_weighting']}")  # 条件を図から確認できるようにする。
        fig.savefig(str(prefix) + f"_{gain}.png", dpi=140)  # gain別に画像を保存する。
        plt.close(fig)  # 描画メモリを解放する。
        pair_count = sum(row["pairs"] for row in result["runs"])  # 採用ペアの総数を示す。
        title = f"{gain.upper()}: raw2={options['require_exactly_two_raw_hits']}; {pair_count:,} pairs / {len(points):,} hit uses\nHG position calibration (LG provisional); beam cut={options['beam_vertex_cut_enabled']}; z0 veto={options.get('z0_track_veto_enabled', False)}"  # 2点ずつ積算することを明記する。
        draw_pixel_map(points, config, title, "hit uses in selected pairs / projected pixel", Path(str(prefix)+f"_{gain}_pixels.png"))  # 方位角図とは別の3投影を保存する。
    Path(str(prefix) + "_comparison.json").write_text(json.dumps(summary, indent=2) + "\n")  # 両gainの条件・件数・fitを1か所にも保存する。


if __name__ == "__main__":  # import時はrawを読まない。
    main()  # コマンドとして呼ばれたときだけ解析する。
