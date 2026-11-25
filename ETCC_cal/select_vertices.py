#!/usr/bin/env python3
"""保存済み頂点へイベント単位のZ側除外を適用し、全トリガーに対する残存率を計算する。"""
import argparse  # 設定ファイルのパスを受け取る。
import csv  # 頂点を1行ずつ読み、大きなCSVを丸ごとメモリに持たない。
import json  # 選別の条件と件数を保存する。
from pathlib import Path  # 入出力パスを扱う。
import sys  # 既存のHA rawファイル収集・フレーム切り出しを読み込む検索先を追加する。
import numpy as np  # 除外イベント集合をNPZから読む。
from settings import DEFAULT_CONFIG, REPO, load_config, repo_path, input_file, output_file  # 較正と保存先の解釈を積算と共通にする。
sys.path.insert(0, str(REPO / "ETCC_analysis"))  # pair_polarization.pyの場所を検索先に加える。
import pair_polarization as psa  # HA rawのボード別ファイル収集と、ヘッダ・フッタでのフレーム切り出しを再使用する。


def scan_ha_hits(name, config):  # 1ランのHA rawからイベントごとのヒットフラグ数を数える（波形は読まない）。
    directories = [str(repo_path(config["data_root"]) / storage / name) for storage in config["storages"]]  # 両storageのランを探す。
    boards = psa.collect_files(directories)  # 完成した.datだけをボードごとにまとめる（.dat.tmpは使わない）。
    if not boards:  # HAの入力が無いのにヒット0として保存しない。
        raise FileNotFoundError(f"HA raw files not found: {directories}")  # 見つからない入力先を明示する。
    keys, counts, frames, resets = [], [], {}, {}  # ボードごとのイベント番号・ヒット数・件数を集める。
    for board, paths in sorted(boards.items()):  # ボードを順に処理する。
        seen = set()  # 同じボードの重複フレームを1回だけ数える。
        epoch, previous = 0, None  # trigの巻き戻りの状態を初期化する。
        frames[f"ha{board:02d}"], resets[f"ha{board:02d}"] = 0, 0  # ボードの読み取り件数を初期化する。
        for path in paths:  # rawを1ファイルずつ読む。
            buf, sizes, n = psa.read_file(path)  # フレームを切り出す（decode_eventsと同じ切り出し）。
            if buf is None or n == 0:  # 空ファイルはヒットなしイベントと混同しない。
                print(f"[WARN] no frames: {path}", flush=True)  # 入力の問題をログに残す。
                continue  # 次のファイルへ進む。
            valid = sizes >= 0x28  # decode_eventsが読まない短いフレームを除く。
            b = buf[valid].astype(np.int64)  # ヘッダのバイトを整数で扱う。
            trig = (b[:, 4] << 24) | (b[:, 5] << 16) | (b[:, 6] << 8) | b[:, 7]  # 0x04のTCUトリガー番号（big endian）。
            datacnt = (((b[:, 13] << 8) | b[:, 14]) >> 6) & 0x3F  # 0x0Dの6ビット: データのあるADC（= 単純なヒットフラグ）。
            n_hit = np.zeros(len(trig), np.int64)  # フレームごとのヒットch数を数える。
            for bit in range(6):  # 6 chのフラグを足す。
                n_hit += (datacnt >> bit) & 1  # decode_eventsのhitと同じビットを数える。
            back = np.r_[previous is not None and trig[0] < previous, trig[1:] < trig[:-1]]  # trigが減った場所を調べる。
            epochs = epoch + np.cumsum(back)  # 巻き戻り後を別のイベント番号空間に分ける（TPC・pair_polarizationと同じ）。
            resets[f"ha{board:02d}"] += int(back.sum())  # 巻き戻りの回数を残す。
            epoch, previous = int(epochs[-1]), int(trig[-1])  # 次のファイルへ状態を引き継ぐ。
            events = epochs * (1 << 32) + trig  # (epoch,trig)を1個の整数にまとめる。
            keep = np.array([e not in seen for e in events.tolist()], dtype=bool)  # このボードで初めてのフレームだけを使う。
            _, first = np.unique(events, return_index=True)  # 同一ファイル内の重複も除く。
            keep &= np.isin(np.arange(len(events)), first)  # 最初の出現だけを残す。
            seen.update(events[keep].tolist())  # 読んだイベントを記録する。
            frames[f"ha{board:02d}"] += int(keep.sum())  # 有効フレーム数を加える。
            keys.append(events[keep])  # イベント番号を保存する。
            counts.append(n_hit[keep])  # そのボードのヒットch数を保存する。
        print(f"{name} ha{board:02d}: {frames[f'ha{board:02d}']:,} frames", flush=True)  # 読み込みの進捗を知らせる。
    event, inverse = np.unique(np.concatenate(keys), return_inverse=True)  # ボードをまたいだイベントの和集合を作る。
    n_hits = np.zeros(len(event), np.int64)  # イベントごとの全ボードのヒットch数。
    np.add.at(n_hits, inverse, np.concatenate(counts))  # 各ボードのヒット数を足す。
    return event, n_hits, frames, resets  # HAが記録した全イベントとヒット数を返す。


def load_ha_hits(name, config, rebuild):  # HAのヒット数をキャッシュから読み、無ければrawから作る。
    path = repo_path(config["outputs"]["veto_dir"]) / f"ha_hit_flags_{name}.npz"  # 除外リストの隣に保存する。
    if rebuild or not path.exists():  # 初回か指定時だけrawを読む。
        event, n_hits, frames, resets = scan_ha_hits(name, config)  # HA rawを走査する。
        np.savez_compressed(path, event=event, n_hits=n_hits, info=json.dumps({"frames_per_board": frames, "trig_resets_per_board": resets}))  # 走査結果を保存する。
    data = np.load(path)  # 保存済みのヒット数を読む。
    return data["event"], data["n_hits"], json.loads(str(data["info"]))  # イベント番号・ヒット数・ボード情報を返す。


def main():  # 選別と件数集計をコマンドから実行する。
    parser = argparse.ArgumentParser(description=__doc__)  # 入力条件の説明を作る。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 共通の解析条件を読む。
    parser.add_argument("--rebuild-ha", action="store_true")  # HA rawを読み直してヒットフラグのキャッシュを作り直す。
    args = parser.parse_args()  # 指定された設定パスを取り出す。
    config = load_config(args.config)  # 較正値とラン・保存先を読み取る。
    output = repo_path(config["result_dir"])  # 整理後の結果ディレクトリを求める。
    scan = json.loads(output_file(config, "track_metadata").read_text())  # raw走査から得た正確な全イベント数を使う。
    rows = {run["run"]: dict(run) for run in scan["runs"]}  # ラン別の記録トリガー数を引き継ぐ。
    veto = {}  # ランごとの除外イベント集合を用意する。
    for name in config["runs"]:  # ランをまたいでトリガー番号を混ぜない。
        data = np.load(repo_path(config["outputs"]["veto_dir"]) / f"z0_track_veto_{name}.npz")  # rawで判定したイベントを読む。
        veto[name] = set(map(int, data["event"]))  # 高速に所属判定できる集合にする。
        rows[name]["vertices_before"] = 0  # 入力頂点の件数を初期化する。
        rows[name]["vertices_after"] = 0  # カット後の頂点件数を初期化する。
        rows[name]["selected_events"] = []  # 頂点が残ったイベント番号を集める。
    source_path = input_file(config, "vertices_all")  # 全頂点を保存した入力を使う。
    target_path = output_file(config, "vertices_selected")  # 選別後を別名で保存する。
    with source_path.open() as source, target_path.open("w") as target:  # 入力を消さず出力だけを作り直す。
        next(source)  # 先頭の説明コメントを読み飛ばす。
        reader = csv.DictReader(source)  # 列名で頂点情報を読み取る。
        fields = ["z0_veto_pass" if field == "fiducial" else field for field in reader.fieldnames]  # 古い全方向fiducialフラグと新条件を区別する。
        target.write("# Selection: no raw TPC hit at Z <= 15 mm; no XY / high-Z veto. Coordinates [mm].\n")  # 選別の意味と座標単位を先頭に書く。
        writer = csv.DictWriter(target, fieldnames=fields)  # 入力と同じ座標列の順番で出力する。
        writer.writeheader()  # 列名を保存する。
        for vertex in reader:  # 全イベントの頂点を順に確認する。
            name = vertex["run"]  # ラン名を取り出す。
            event = (int(vertex["epoch"]) << 32) + int(vertex["trig"])  # raw走査と同じ一意イベント番号を作る。
            rows[name]["vertices_before"] += 1  # 選別前の頂点を数える。
            if event in veto[name]:  # 全ボードのどこかに低Zヒットがあったイベントか調べる。
                continue  # 頂点自体のZに関係なくイベント全体を除外する。
            if float(vertex["z"]) <= config["z_veto"]:  # 生存頂点が低Z側なら選別データの不整合になる。
                raise RuntimeError("Veto cache and vertex CSV are inconsistent")  # 不整合のある結果を正常として保存しない。
            vertex.pop("fiducial")  # 旧来の全方向fiducial情報を出力条件と混同させない。
            vertex["z0_veto_pass"] = "1"  # このファイルの全イベントがZ側除外を通過したことを示す。
            vertex["dist_to_line"] = ""  # 旧PCA軸からの距離を新解析へ引き継がない。
            writer.writerow(vertex)  # 生存イベントの頂点を保存する。
            rows[name]["vertices_after"] += 1  # 残った頂点を数える。
            rows[name]["selected_events"].append(event)  # HAのヒットと照合するためにイベント番号を残す。
    summary = []  # レビューしやすい小さな選別表を作る。
    for name, row in rows.items():  # 各ランの全イベント数を分母にする。
        total = row["all_events"]  # ヒットなしも含むraw記録トリガー数を使う。
        if max(row["frames_per_board"].values()) != total:  # TPCのイベントがtrig 1..Nに抜けなく並ぶ前提を確かめる。
            raise RuntimeError(f"{name}: TPC events are not trig 1..{total}; cannot match HA by range")  # 前提が崩れたら数えない。
        veto_sorted = np.array(sorted(veto[name]), dtype=np.int64)  # 除外集合を検索用の配列にする。
        if len(veto_sorted) and (veto_sorted[0] < 1 or veto_sorted[-1] > total):  # 除外イベントもtrig 1..Nの中にあるはず。
            raise RuntimeError(f"{name}: veto events outside trig 1..{total}")  # 前提が崩れたら数えない。
        ha_event, ha_hits, ha_info = load_ha_hits(name, config, args.rebuild_ha)  # HAのイベントごとのヒットフラグ数を得る。
        in_tpc = (ha_event >> 32 == 0) & ((ha_event & 0xFFFFFFFF) >= 1) & ((ha_event & 0xFFFFFFFF) <= total)  # TPCの全イベント（epoch 0、trig 1..N）に入るか。
        passing = in_tpc & ~np.isin(ha_event, veto_sorted)  # Z側除外を通ったイベントか。
        vertex_pass = np.isin(ha_event, np.array(row["selected_events"], dtype=np.int64))  # 除外を通り頂点もあるイベントか。
        hit = ha_hits > 0  # シンチのヒットフラグが1個以上あるイベント。
        def multiplicity(mask, denominator):  # ヒットフラグ数ごとのイベント数（0ヒットは分母から引く）。
            n = [int(np.sum(mask & (ha_hits == k))) for k in (1, 2)] + [int(np.sum(mask & (ha_hits >= 3)))]  # 1、2、3以上を数える。
            return {"0": denominator - sum(n), "1": n[0], "2": n[1], "3+": n[2]}  # HAに記録の無いイベントは0ヒットに入る。
        scint = {  # シンチの単純なヒットフラグ（ADCのデータあり）での集計。
            "definition": "HA frame header DataCount bit (ch with data) of any board; matched to TPC by (epoch, trig); no HG / peak_id / timing cut",  # 何を1ヒットと数えるか。
            "events_with_hit_all": int(np.sum(in_tpc & hit)),  # 全イベントのうちHAヒットが1個以上あるイベント数。
            "events_with_hit_all_percent": 100 * np.sum(in_tpc & hit) / total,  # 全イベントに対する割合 [%]。
            "events_with_hit_passing": int(np.sum(passing & hit)),  # 除外を通ったイベントのうちHAヒットがあるもの。
            "events_with_hit_passing_percent_of_passing": 100 * np.sum(passing & hit) / row["passing_events"],  # 通過イベントに対する割合 [%]。
            "events_with_hit_passing_percent_of_all": 100 * np.sum(passing & hit) / total,  # 全イベントに対する割合 [%]。
            "events_with_hit_vertex": int(np.sum(vertex_pass & hit)),  # 除外を通り頂点もあるイベントのうちHAヒットがあるもの。
            "events_with_hit_vertex_percent_of_vertex": 100 * np.sum(vertex_pass & hit) / row["vertices_after"],  # 頂点が残ったイベントに対する割合 [%]。
            "hits_all": int(ha_hits[in_tpc].sum()),  # 全イベントのヒットフラグの総数。
            "hits_passing": int(ha_hits[passing].sum()),  # 除外を通ったイベントのヒットフラグの総数。
            "hits_passing_percent_of_all_hits": 100 * ha_hits[passing].sum() / max(int(ha_hits[in_tpc].sum()), 1),  # 除外後に残ったヒットの割合 [%]。
            "hits_vertex": int(ha_hits[vertex_pass].sum()),  # 頂点が残ったイベントのヒットフラグの総数。
            "multiplicity_all": multiplicity(in_tpc, total),  # 全イベントのヒット数の分布。
            "multiplicity_passing": multiplicity(passing, row["passing_events"]),  # 除外後のヒット数の分布。
            "ha_events_outside_tpc": int(np.sum(~in_tpc)),  # TPCに対応するイベントの無いHAイベント（epoch>0など）。
            "ha_frames_per_board": ha_info["frames_per_board"],  # ボードごとの読み取り件数。
            "ha_trig_resets_per_board": ha_info["trig_resets_per_board"],  # ボードごとのtrigの巻き戻り回数。
        }  # シンチの集計を完成させる。
        summary.append({  # トリガー通過と頂点成功を別々に保存する。
            "run": name,  # 対象のラン。
            "all_events": total,  # 全イベントの分母。
            "events_with_hits": row["events_with_hits"],  # TPCが1ヒット以上あるイベント数。
            "veto_events": row["veto_events"],  # Z側条件で除外したイベント数。
            "passing_events": row["passing_events"],  # ヒットなしも含むカット通過イベント数。
            "passing_percent": 100 * row["passing_events"] / total,  # 全イベントに対する通過率 [%]。
            "vertices_before": row["vertices_before"],  # 選別前に頂点を求められたイベント数。
            "vertices_after": row["vertices_after"],  # 選別を通過して頂点もあるイベント数。
            "vertex_percent": 100 * row["vertices_after"] / total,  # 全イベントに対する頂点残存率 [%]。
            "scintillator": scint,  # シンチのヒットフラグでの集計。
        })  # ラン別の集計を完成させる。
    result = {"selection": "Any raw anode OR cathode hit Z <= 15 mm vetoes the whole event; no XY/high-Z veto", "runs": summary}  # 選別の定義を添える。
    output_file(config, "selection_summary").write_text(json.dumps(result, indent=2) + "\n")  # 小さな件数表を保存する。
    print(json.dumps(result, indent=2))  # 実行直後にも件数を確認できるよう表示する。


if __name__ == "__main__":  # importでファイルを書き換えない。
    main()  # コマンドとして呼ばれた場合だけ実行する。
