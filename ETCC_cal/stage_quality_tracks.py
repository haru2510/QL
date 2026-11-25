#!/usr/bin/env python3
"""stage04のHG/LG別・raw1/2別の採用イベントでTPCのXZ/YZと3D同時ヒットXYを作る。"""
from cutflow_report import write_cutflow  # 全run合計と日本語の説明を先頭に保存する。
import argparse  # 解析設定を指定する。
import csv  # シンチ品質条件を通過したイベント一覧を読む。
import hashlib  # 選別イベント集合をキャッシュへ記録する。
import json  # 条件と件数を保存する。
from pathlib import Path  # stage04へ保存する。
import numpy as np  # イベント集合と積算を扱う。
from settings import load_config, repo_path  # 外部JSONの入出力を使う。
from stage_coincidence_xy import scan_run  # 同じraw走査でXZ/YZとXYを同時に作る。
from track_runs import draw  # turboの4種類の積算図を共用する。

BRANCHES = [('single_hit','hg',1),('single_hit','lg',1),('two_hits','hg',2),('two_hits','lg',2)]  # 4条件を同時に処理する。


def accepted_events(config):  # ピクセル図と完全に同じ採用イベントを使う。
    events = {run:[[] for _ in BRANCHES] for run in config['runs']}  # ラン別の集合を用意する。
    for index,(branch,gain,multiplicity) in enumerate(BRANCHES):  # raw1/2とgainを分ける。
        key = 'single_hit_pixel_prefix' if multiplicity==1 else 'pixel_polarization_prefix'  # 採用CSVの出力名を読む。
        suffix = 'events' if multiplicity==1 else 'pairs'  # 1イベント1行のCSVを使う。
        path = Path(str(repo_path(config['outputs'][key]))+f'_{gain}_{suffix}.csv')  # 既に位置品質を選別した一覧を指定する。
        with path.open() as stream:  # rawを再較正せず採用IDだけ読む。
            for row in csv.DictReader(stream):  # イベント番号を復元する。
                if int(row['raw_hits'])!=multiplicity:  # raw数の異なる解析結果を流用しない。
                    raise ValueError(f'Unexpected multiplicity: {path}')  # 条件不一致を示す。
                events[row['run']][index].append((int(row['epoch'])<<32)+int(row['trig']))  # TPCと同じキーにする。
    report = json.loads((repo_path(config['result_dir'])/'stage_flag_quality_cutflow.json').read_text())  # シンチstage04の件数を読む。
    for index,(_,gain,multiplicity) in enumerate(BRANCHES):  # 採用イベント重複と件数を確認する。
        total = 0  # 全ラン合計を数える。
        for run in events:  # ランごとに一意なイベント集合へ直す。
            values = np.array(sorted(events[run][index]),np.int64)  # フレーム検索用の昇順にする。
            if len(np.unique(values))!=len(values):  # 同じイベントを複数ペアとして数える設定を除く。
                raise ValueError(f'Duplicate event: {run}, {gain}, {multiplicity}')  # raw2の1イベント1ペアを要求する。
            events[run][index] = values  # 検証した集合を保存する。
            total += len(values)  # 採用数を加算する。
        assert total==report['branches'][str(multiplicity)]['gains'][gain]['after_quality_events']  # シンチ図の採用数と一致させる。
    return events  # 4条件のイベント集合を返す。


def main():  # 品質選別後だけのTPC積算を作る。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行条件を用意する。
    parser.add_argument('--config',required=True,type=Path)  # 入出力を共通JSONで指定する。
    parser.add_argument('--plot-only',action='store_true')  # 保存配列から再描画する。
    args = parser.parse_args()  # 引数を読む。
    config = load_config(args.config)  # 共通の条件を読む。
    root = repo_path(config['result_dir'])  # 対象データセットの保存先を指定する。
    cache_dir = root/'shared/stage04_tracks'  # ラン別の4条件積算を置く。
    cache_dir.mkdir(parents=True,exist_ok=True)  # 保存先を作る。
    selections = accepted_events(config)  # HG/LGの採用イベントを読む。
    total_raw,total_xy = np.zeros((4,2,768,1024),np.int64),np.zeros((4,96,96),np.int64)  # 4条件を確保する。
    summaries = []  # ランごとの分母を残す。
    for run in config['runs']:  # ランを混ぜずに処理する。
        signature = hashlib.sha256(b''.join(values.tobytes() for values in selections[run])).hexdigest()  # 選別集合の変更を検出する。
        target = cache_dir/f'{run}.npz'  # 正常完了したランだけ再利用する。
        if not target.exists() and not args.plot_only:  # 未完了のランだけSSDを読む。
            with np.load(root/'shared/stage_tracks'/f'{run}.npz') as previous:  # 元の積算と入力が一致することを確認する。
                inputs = json.loads(str(previous['summary']))['inputs']  # 以前のrawファイル一覧を読む。
            for entry in inputs:  # 入力更新を混ぜない。
                stat = Path(entry['path']).stat()  # 現在のサイズと時刻を確認する。
                if [stat.st_size,stat.st_mtime_ns]!=[entry['bytes'],entry['mtime_ns']]:  # 同じrawであることを検証する。
                    raise ValueError(f'Raw changed: {entry["path"]}')  # 不一致なら停止する。
            histogram = np.zeros_like(total_raw)  # このランのXZ/YZを確保する。
            xy,coincident = scan_run(run,config,selections[run],raw_histogram=histogram)  # 採用フレームだけをデコードし、XYも同時に作る。
            with np.load(root/'shared/stage_tracks'/f'{run}.npz') as previous:  # stage03の上限を使って包含関係を確認する。
                for index,(_,_,multiplicity) in enumerate(BRANCHES):  # 品質選別でヒットが増えないことを検証する。
                    assert np.all(histogram[index]<=previous['histogram'][3 if multiplicity==1 else 4])  # raw1/2の部分集合になっているか確認する。
            clock_limit = int(np.floor((config['z_veto']-config['z_gem'])/config['velocity']+config['clock_offset']))  # 同じZしきい値へ直す。
            assert histogram[:,:,:,:clock_limit+1].sum()==0  # 品質選別後も全TPCのZ側イベント除外を保っているか確認する。
            summary = dict(run=run,selected_events=[len(values) for values in selections[run]],events_with_any_selected_coincidence=coincident,strip_clock_hits=[int(h.sum()) for h in histogram],xy_combinations=[int(h.sum()) for h in xy])  # ヒット数とイベント数を分ける。
            np.savez_compressed(target,raw=histogram,xy=xy,summary=json.dumps(summary),config=json.dumps(config),signature=signature)  # 再描画と再開用に保存する。
        with np.load(target) as saved:  # 正常完了したランを合算する。
            if json.loads(str(saved['config']))!=config or str(saved['signature'])!=signature:  # 異なる選別を混ぜない。
                raise ValueError(f'Selection/config mismatch: {target}')  # 再作成が必要と示す。
            total_raw += saved['raw']  # XZ/YZを足す。
            total_xy += saved['xy']  # XYを足す。
            summaries.append(json.loads(str(saved['summary'])))  # 採用数を記録する。
        print('COMPLETE STAGE04',run,flush=True)  # ランごとの完了を残す。
    for index,(branch,gain,multiplicity) in enumerate(BRANCHES):  # stage04のgain別フォルダへ保存する。
        directory = root/'with_z0_cut/stages'/f'04_{branch}_quality'/f'tpc_{gain}'  # シンチ図の直下へ分類する。
        directory.mkdir(parents=True,exist_ok=True)  # 保存先を作る。
        selected = sum(row['selected_events'][index] for row in summaries)  # 元シンチ図と同じイベント数を示す。
        title = f"{config.get('data_label', ', '.join(config['runs']))}; stage04 {gain.upper()}; raw{multiplicity}; N={selected:,} events"  # 条件を図中へ明記する。
        draw(total_raw[index],dict(config,color_percentiles=[10,99.5]),directory/'tpc_tracks',title,selection_label='Beam-axis THEN Z-side veto THEN raw multiplicity THEN gain quality; XZ/YZ raw strips; XY same-clock X-Y combinations',xy_histogram=total_xy[index])  # 同じ4種類の3パネルを作る。
        (directory/'selection.json').write_text(json.dumps(dict(stage='04',gain=gain,raw_multiplicity=multiplicity,selected_events=selected,source='exact accepted scintillator CSV event IDs; no additional TPC hit cut',runs=summaries),indent=2)+'\n')  # 条件と件数を残す。
    np.savez_compressed(cache_dir/'all_runs.npz',raw=total_raw,xy=total_xy,config=json.dumps(config))  # SSDなしの再描画を可能にする。
    write_cutflow(root/'stage04_track_cutflow.json', dict(branches=BRANCHES,runs=summaries), 'quality_tracks')  # HG/LGごとの最終TPC採用数を保存する。
    print('ALL STAGE04 TRACKS COMPLETE',root,flush=True)  # 全処理の完了を知らせる。


if __name__=='__main__':  # importだけではrawを読まない。
    main()  # コマンドから呼ばれたときだけ実行する。
