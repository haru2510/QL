#!/usr/bin/env python3
"""候補イベントのHA rawからLG重心を求め、HG用位置較正を暫定流用してピクセルへ割り当てる。"""
import argparse  # 共通JSONと個別候補を指定する。
import csv  # ピクセル較正とブロック配置を読む。
import hashlib  # 使用した較正表をハッシュで識別する。
import json  # 更新したイベント点群を保存する。
import re  # 保存されたHAラベルからボードとchを読む。
import sys  # 既存のHAデコーダを読み込む。
from pathlib import Path  # 入力rawと出力先を扱う。
import numpy as np  # 波形と重心・座標を計算する。
from scipy.interpolate import RBFInterpolator  # EricaNextと同じthin-plate splineを使う。
from settings import DEFAULT_CONFIG, REPO, load_config, repo_path  # 入出力を外部JSONから読む。
sys.path.insert(0, str(REPO / "ETCC_analysis"))  # 既存のファイル収集処理を探せるようにする。
import pair_polarization as psa  # ADC番号順のraw読み取りを再利用する。


def read_rows(path):  # コメント付きCSVを辞書の行に直す。
    with path.open() as source:  # 読み終えたファイルを閉じる。
        return list(csv.DictReader(line for line in source if line.strip() and not line.startswith("#")))  # 説明コメントを除いて読む。


def decode_targets(run, targets, config):  # 必要なボードとトリガーだけをrawから取り出す。
    directories = [str(repo_path(config["data_root"]) / storage / run) for storage in config["storages"]]  # 両storageを検索する。
    boards = psa.collect_files(directories)  # 完成HA rawをボードごとにまとめる。
    found = {}  # (board,trig,ch)ごとの波形解析結果を保存する。
    for board in sorted({key[0] for key in targets}):  # 候補に使うボードだけを読む。
        pending = {key for key in targets if key[0] == board}  # このボードで探すヒットを選ぶ。
        for filename in boards.get(board, []):  # rawファイルを連番順に読む。
            buf, sizes, count = psa.read_file(filename)  # 生データのフレームを取り出す。
            if buf is None or count == 0:  # 空のファイルは解析しない。
                continue  # 次の入力に進む。
            decoded = psa.decode_events(buf, sizes, False, psa.DEFAULT_TEMP_A, psa.DEFAULT_TEMP_B, False)  # 既存のピーク周辺積分でHG/LGを計算する。
            trig, _, hit = decoded[:3]  # トリガー番号とrawヒットフラグを取り出す。
            for key in list(pending):  # 未取得の候補を探す。
                _, target, channel = key  # この候補のトリガーとchを読む。
                indices = np.flatnonzero((trig == target) & (hit[:, channel] == 1))  # 同じトリガーのrawヒットだけを選ぶ。
                if not len(indices):  # このファイルに対象がなければ続ける。
                    continue  # 次の候補を確認する。
                if len(indices) != 1:  # 同じトリガーが複数あれば曖昧なまま選ばない。
                    raise RuntimeError(f"Ambiguous trigger: {key} in {filename}")  # 一意に同定できない入力を明示する。
                index = int(indices[0])  # 対象フレームの配列位置を得る。
                found[key] = {  # 位置推定と飽和確認に使う値を保存する。
                    "raw_file": str(filename),  # どのrawから読んだか記録する。
                    "hg": float(decoded[7][index, channel]),  # HGの波形積分値 [ADC sample]。
                    "lg": float(decoded[8][index, channel]),  # LGの波形積分値 [ADC sample]。
                    "x_hg": float(decoded[9][index, channel]),  # 比較用のHG重心 [無次元]。
                    "y_hg": float(decoded[10][index, channel]),  # 比較用のHG重心 [無次元]。
                    "x_lg": float(decoded[11][index, channel]),  # 使用するLGのX重心 [無次元]。
                    "y_lg": float(decoded[12][index, channel]),  # 使用するLGのY重心 [無次元]。
                    "peak_id": int(decoded[13][index, channel]),  # LG和で決まるピークサンプル番号。
                    "saturated_hg": bool(decoded[15][index, channel]),  # HG ADCの飽和フラグ。
                    "saturated_lg": bool(decoded[16][index, channel]),  # LG ADCの飽和フラグ。
                }  # 対象ヒットの読み取り結果を完成させる。
                pending.remove(key)  # 見つかった候補を未取得リストから外す。
            if not pending:  # このボードの必要なヒットが全部得られたら終える。
                break  # 不要な後続ファイルを読まない。
        if pending:  # 対象が見つからない場合をブロック中心で黙って代用しない。
            raise RuntimeError(f"Raw hits not found: {sorted(pending)}")  # どの入力が不足したか伝える。
    return found  # 候補ヒットのLG/HG解析値を返す。


def main():  # HG較正を流用したLG位置を求める。
    parser = argparse.ArgumentParser(description=__doc__)  # 入力の説明を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 他のrun06・07解析と同じJSONを使う。
    args = parser.parse_args()  # 指定された設定を読む。
    config = load_config(args.config)  # raw入力と較正・保存先を取得する。
    sources = [repo_path(path) for path in config["inputs"]["event_candidates"]]  # 更新する保存済み候補JSONを選ぶ。
    events = [(path, json.loads(path.read_text())) for path in sources]  # 元のTPC点群を保持したまま読む。
    geometry_path = repo_path(config["inputs"]["scintillator_geometry"])  # ピクセル位置を検出器座標へ変換する配置を読む。
    geometry = {(int(row["board"]), int(row["ch"])): row for row in read_rows(geometry_path)}  # ボードとADC番号順chで配置を検索する。
    targets = {}  # ランごとに必要なrawヒットをまとめる。
    for _, event in events:  # 保存された各候補を確認する。
        if event.get("epoch", 0) != 0:  # この抽出では巻き戻り後のトリガーを未対応とする。
            raise ValueError("Nonzero epoch requires explicit raw event matching")  # 異なるepochを取り違えない。
        for hit in event["ha_hits"]:  # 同じトリガーのシンチヒットを読む。
            match = re.search(r"ha(\d+) ch(\d+)", hit["label"])  # ラベルから現行ボード・chを取得する。
            board, channel = map(int, match.groups())  # ADC番号順の整数に変換する。
            hit["board"], hit["ch"] = board, channel  # 今後はラベル解析なしでも参照できるよう残す。
            targets.setdefault(event["run"], set()).add((board, int(event["trig"]), channel))  # 同じrawを何度も読まないようまとめる。
    raw = {run: decode_targets(run, keys, config) for run, keys in targets.items()}  # 各ランの必要なLG重心を取得する。
    destination = repo_path(config["outputs"]["lg_event_dir"])  # 元のブロック中心図とは別の保存先を使う。
    destination.mkdir(parents=True, exist_ok=True)  # LG版のフォルダを作る。
    results = []  # 位置と飽和状態を一覧にする。
    for source, event in events:  # 各イベントのHA位置を更新する。
        for hit in event["ha_hits"]:  # 保存済みのTPC点群は変更せずHAだけ処理する。
            board, channel = hit["board"], hit["ch"]  # 現行のボード・chを取り出す。
            values = raw[event["run"]][(board, int(event["trig"]), channel)]  # このシンチのraw解析値を取得する。
            calibration = repo_path(config["inputs"]["ha_position_calibration_dir"]) / f"position_ha{board:02d}.csv"  # HG用の位置較正表を選ぶ。
            rows = [row for row in read_rows(calibration) if int(row["ch"]) == channel]  # 現行chが一致する64ピークを使う。
            if len(rows) != 64:  # 8×8の対応表が揃っているか確認する。
                raise ValueError(f"Expected 64 peaks: {calibration}, ch{channel}")  # 不完全な較正で位置を求めない。
            observed = np.array([[float(row["x"]), float(row["y"])] for row in rows])  # HG重心平面のピーク座標を得る。
            ideal = np.array([[int(row["n"]), int(row["m"])] for row in rows])  # 各ピークのピクセル番号を得る。
            warp = RBFInterpolator(observed, ideal, kernel="thin_plate_spline", smoothing=0)  # EricaNextと同じ補間でHG位置較正を再現する。
            assert np.max(np.abs(warp(observed) - ideal)) < 1e-6  # 較正の64点を正しいピクセルへ戻せることを確認する。
            centroid = np.array([[values["x_lg"], values["y_lg"]]])  # 飽和しにくいLG重心をHG用較正へ入力する。
            corrected = warp(centroid)[0]  # 連続値のピクセル座標を求める。
            pixel = np.floor(corrected + 0.5).astype(int)  # EricaNextと同じ最近傍の整数ピクセルに割り当てる。
            valid = bool(values["lg"] > 0 and np.all(np.isfinite(corrected)) and np.all((pixel >= 1) & (pixel <= 8)) and not values["saturated_lg"])  # LG飽和・無信号・範囲外を位置確定に使わない。
            row = geometry[(board, channel)]  # ブロック中心とピクセルの向きを取得する。
            center = np.array([float(row[key]) for key in ("x0", "y0", "z0")])  # 元のブロック中心 [mm]を得る。
            u = np.array([float(row[key]) for key in ("ux", "uy", "uz")])  # nが増える方向の単位ベクトルを得る。
            v = np.array([float(row[key]) for key in ("vx", "vy", "vz")])  # mが増える方向の単位ベクトルを得る。
            estimate = center + float(row["pitch"]) * ((pixel[0] - 4.5) * u + (pixel[1] - 4.5) * v)  # ピクセル中心を検出器座標 [mm]へ変換する。
            hit.update(values)  # HG/LGの積分・重心・飽和フラグを保存する。
            hit["block_center_mm"] = center.tolist()  # 元の位置を比較できるように残す。
            hit["position_mm"] = (estimate if valid else center).tolist()  # 位置不確定ならブロック中心へ戻す。
            hit["pixel_nm"] = pixel.tolist() if valid else None  # 有効なときだけピクセル番号を採用する。
            hit["corrected_lg_nm"] = corrected.tolist()  # 範囲外の結果も診断用に残す。
            hit["position_valid"] = valid  # ブロック中心への代用を明示する。
            hit["position_method"] = "LG centroid + HG position calibration (provisional)" if valid else "block centre fallback: LG saturated / nonpositive / outside 8x8"  # 暫定較正と代用理由を記録する。
            hit["calibration_file"] = str(calibration)  # 使用したHG位置較正の場所を残す。
            hit["calibration_sha256"] = hashlib.sha256(calibration.read_bytes()).hexdigest()  # 表の中身を後から照合できるようにする。
            base_label = f"ha{board:02d} ch{channel} (output ch{5-channel})"  # 現行chと旧出力番号を併記する。
            hit["label"] = base_label + (f" LG pixel {tuple(pixel)}" if valid else " LG invalid; block centre")  # 図のヒット位置の意味を表示する。
            results.append({"run": event["run"], "trig": event["trig"], "board": board, "ch": channel, "pixel_nm": hit["pixel_nm"], "saturated_hg": values["saturated_hg"], "saturated_lg": values["saturated_lg"], "position_mm": hit["position_mm"]})  # 主要な結果を一覧化する。
        event["ha_position_method"] = "LG centroid with HG position calibration (provisional); invalid LG uses block centre"  # 図全体に較正流用の説明を付ける。
        event["source_event_json"] = str(source)  # 元の候補とTPC点群の出典を残す。
        target = destination / (source.stem + "_lg.json")  # LG版を別名で保存する。
        target.write_text(json.dumps(event, indent=2) + "\n")  # SSDなしで再描画できる形で保存する。
        print(f"Saved {target}", flush=True)  # 更新したイベントの場所を報告する。
    (destination / "lg_position_summary.json").write_text(json.dumps(results, indent=2) + "\n")  # 両候補のピクセルと飽和状態をまとめる。
    print(json.dumps(results, indent=2), flush=True)  # 実行結果を確認できるよう表示する。


if __name__ == "__main__":  # importではrawを読み込まない。
    main()  # コマンドとして呼ばれたときだけ処理する。
