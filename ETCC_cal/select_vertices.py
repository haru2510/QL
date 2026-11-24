#!/usr/bin/env python3
"""保存済み頂点へイベント単位のZ側除外を適用し、全トリガーに対する残存率を計算する。"""
import argparse  # 設定ファイルのパスを受け取る。
import csv  # 頂点を1行ずつ読み、大きなCSVを丸ごとメモリに持たない。
import json  # 選別の条件と件数を保存する。
from pathlib import Path  # 入出力パスを扱う。
import numpy as np  # 除外イベント集合をNPZから読む。
from settings import DEFAULT_CONFIG, load_config, repo_path, input_file, output_file  # 較正と保存先の解釈を積算と共通にする。


def main():  # 選別と件数集計をコマンドから実行する。
    parser = argparse.ArgumentParser(description=__doc__)  # 入力条件の説明を作る。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 共通の解析条件を読む。
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
    summary = []  # レビューしやすい小さな選別表を作る。
    for name, row in rows.items():  # 各ランの全イベント数を分母にする。
        total = row["all_events"]  # ヒットなしも含むraw記録トリガー数を使う。
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
        })  # ラン別の集計を完成させる。
    result = {"selection": "Any raw anode OR cathode hit Z <= 15 mm vetoes the whole event; no XY/high-Z veto", "runs": summary}  # 選別の定義を添える。
    output_file(config, "selection_summary").write_text(json.dumps(result, indent=2) + "\n")  # 小さな件数表を保存する。
    print(json.dumps(result, indent=2))  # 実行直後にも件数を確認できるよう表示する。


if __name__ == "__main__":  # importでファイルを書き換えない。
    main()  # コマンドとして呼ばれた場合だけ実行する。
