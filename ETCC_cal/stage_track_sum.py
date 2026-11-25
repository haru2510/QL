#!/usr/bin/env python3
"""ビーム軸→Z側→raw1/2の累積条件でTPC全ストリップ・クロックヒットを積算する。"""
from cutflow_report import write_cutflow  # 全run合計と日本語の説明を先頭に保存する。
import argparse  # 外部JSONと再描画を指定する。
import json  # 条件とイベント数を保存する。
from pathlib import Path  # 結果を段階別に保存する。
import numpy as np  # イベント対応と投影を計算する。
from numba import njit  # rawのビット走査を高速化する。
from settings import load_config, repo_path  # 共通のパス指定を使う。
from track_sum import decoder  # TPCの読み取りと配置を再利用する。
from track_runs import draw  # カラー調整版と全レンジ版を同じ方法で作る。
from pixel_polarization import vertex_table, fitted_axis  # 既存の頂点と保存軸を使う。

STAGES = ['00_before_cuts', '01_beam_axis', '02_after_z0', '03_raw1', '03_raw2']  # raw1とraw2は第2段から分岐する。
GROUPS = ['without_z0_cut/stages/00_before_cuts', 'without_z0_cut/stages/01_beam_axis', 'with_z0_cut/stages/02_beam_axis_then_z0', 'with_z0_cut/stages/03_raw_single_hit', 'with_z0_cut/stages/03_raw_two_hits']  # シンチ図と同じ保存先を使う。


@njit  # 全段階を1回のrawビット走査で処理する。
def accumulate_stages(buf, offsets, records, keep, masks, axis, positions, histogram):  # フレームごとの選別を全rawヒットへ適用する。
    hits = np.zeros(len(offsets), np.int64)  # ヒットなしを含めたイベント数を確認する。
    for frame in range(len(offsets)):  # ファイル中のフレームを順に読む。
        if not keep[frame]:  # 同じボードの重複フレームを除く。
            continue  # 初出だけを数える。
        first = offsets[frame] + decoder.HEAD_SIZE + decoder.FADC_SIZE + 4  # rawレコードの開始位置を求める。
        for record in range(records[frame]):  # 早い16クロックに限定せず全ヒットを使う。
            start = first + decoder.RECORD_SIZE * record  # 次のレコードへ進む。
            clock = decoder.read_word(buf, start) & 0x7ff  # クロックを取り出す。
            if clock >= decoder.N_CLOCKS:  # 不正クロックを除く。
                continue  # 有効範囲だけ積算する。
            for channel in range(decoder.N_FADC_CH):  # 4ワードのビットを読む。
                word = decoder.read_word(buf, start + 4 * (4-channel))  # raw仕様の順番で読む。
                for bit in range(32):  # 各ストリップのヒットを調べる。
                    if (word >> bit) & 1:  # ヒットのあるストリップを採用する。
                        coordinate = positions[channel * 32 + bit]  # 配置を共通座標ビンに直す。
                        hits[frame] += 1  # ToTを含むストリップ・クロックヒット数を数える。
                        for stage in range(histogram.shape[0]):  # 同じrawヒットを通過した段階へ積算する。
                            if masks[stage, frame]:  # イベント単位の条件を適用する。
                                histogram[stage, axis, coordinate, clock] += 1  # 段階別の投影へ加算する。
    return hits  # 分母とヒットありイベントの区別に使う。


def selection_sets(run, config, vertices, point, direction):  # rawを読む前にイベント集合を確定する。
    v = vertices.get(run, np.zeros((0, 4)))  # このランの全頂点を読む。
    distances = np.linalg.norm(np.cross(v[:, 1:4]-point, direction), axis=1)  # 無限直線との最短3D距離を計算する。
    beam = v[distances <= config['polarization']['beam_max_distance_mm'], 0].astype(np.int64)  # 第1段の通過イベントを得る。
    with np.load(repo_path(config['outputs']['veto_dir']) / f'z0_track_veto_{run}.npz') as saved:  # 全TPCによる除外リストを読む。
        after_z = np.setdiff1d(beam, saved['event'])  # ビーム軸条件の後にZ側除外をかける。
    with np.load(repo_path(config['outputs']['polarization_hit_cache_dir']) / f'{run}.npz') as saved:  # 全raw HAヒットを読む。
        events, counts = np.unique(saved['event'], return_counts=True)  # ピクセル有効性より前のraw数を確定する。
    return [None, beam, after_z, np.intersect1d(after_z, events[counts == 1]), np.intersect1d(after_z, events[counts == 2])]  # 最後の2条件は並列に分岐する。


def main():  # 各ランのrawを1回ずつ読み、段階別に合算する。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行条件を説明する。
    parser.add_argument('--config', required=True, type=Path)  # 全条件を外部JSONで指定する。
    parser.add_argument('--plot-only', action='store_true')  # 保存配列だけで再描画する。
    args = parser.parse_args()  # 指定を読む。
    config = load_config(args.config)  # 入出力と較正を読む。
    root = repo_path(config['result_dir'])  # 全ランの結果保存先を求める。
    cache_dir = root / 'shared/stage_tracks'  # 再開用のラン別配列を置く。
    cache_dir.mkdir(parents=True, exist_ok=True)  # 保存先を作る。
    vertices = vertex_table(repo_path(config['inputs']['vertices_all']))  # Z側選別前の頂点を読む。
    point, direction = fitted_axis(repo_path(config['inputs']['polarization_beam_axis']))  # 全ランへ同じ保存軸を使う。
    combined = np.zeros((5, 2, 768, 1024), np.int64)  # 5段階の配列を確保する。
    summaries = []  # 分母と通過数を保存する。
    for run in config['runs']:  # ランをまたいでイベント番号を混ぜない。
        target = cache_dir / f'{run}.npz'  # 正常完了したランだけ再利用する。
        if not target.exists() and not args.plot_only:  # 未完了のランだけrawを走査する。
            selected = selection_sets(run, config, vertices, point, direction)  # 各段階の通過集合を求める。
            histogram = np.zeros_like(combined)  # このランだけの配列を作る。
            all_events, with_hits, manifests = set(), set(), []  # 全TPCトリガーを分母にする。
            boards = decoder.collect_boards([str(repo_path(config['data_root'])/storage/run) for storage in config['storages']], 0)  # 完成rawを全数集める。
            if not boards:  # 空の入力はエラーにする。
                raise FileNotFoundError(run)  # ランを明示する。
            for (side, board), paths in boards.items():  # 同じボード内で重複とepochを管理する。
                seen, epoch, previous = set(), 0, None  # 巻き戻りの状態を初期化する。
                upic, local_board = divmod(board, decoder.N_BOARDS_PER_UPIC)  # uPIC内の位置を得る。
                axis, sign, const = decoder.STRIP_AXES[upic][0 if side == 'a' else 1]  # 向きを含めた配置を読む。
                mm = config['pitch_strip'] * (sign * (np.arange(decoder.N_STRIPS)+local_board*decoder.N_STRIPS)+const)  # mmへ換算する。
                positions = np.rint((mm+153.4)/0.4).astype(np.int64)  # 既存と同じ座標ビンを使う。
                for path in paths:  # ファイル単位で処理する。
                    stat = Path(path).stat()  # 入力を記録する。
                    manifests.append(dict(path=path, bytes=stat.st_size, mtime_ns=stat.st_mtime_ns))  # 検証用一覧を残す。
                    buf, trig, offsets, records = decoder.scan_file(path)  # フレーム位置を読む。
                    if not len(trig):  # 空ファイルを除く。
                        continue  # 次へ進む。
                    back = np.r_[previous is not None and trig[0] < previous, trig[1:] < trig[:-1]]  # トリガー巻き戻りを調べる。
                    epochs = epoch + np.cumsum(back)  # 同じボードのepochを保つ。
                    epoch, previous = int(epochs[-1]), int(trig[-1])  # 次のファイルへ引き継ぐ。
                    events = epochs * (1 << 32) + trig  # 照合キーを作る。
                    unique, first = np.unique(events, return_index=True)  # 同一ファイル内の重複も除く。
                    keep = np.zeros(len(events), bool)  # 初出フレームを示す。
                    keep[first] = ~np.isin(unique, np.fromiter(seen, dtype=np.int64))  # 既読の番号を除く。
                    seen.update(map(int, events[keep]))  # 重複除去用の集合を更新する。
                    masks = np.vstack([np.ones(len(events), bool)] + [np.isin(events, values) for values in selected[1:]])  # フレームへ累積条件を対応づける。
                    hits = accumulate_stages(buf, offsets, records, keep, masks, int(axis == 'y'), positions, histogram)  # rawは全段階まとめて1回走査する。
                    all_events.update(map(int, events[keep]))  # ヒットなしも含む全TPCイベントを数える。
                    with_hits.update(map(int, events[keep & (hits > 0)]))  # ヒットありの分母も別に記録する。
                print(f'{run} {side}{board:02d}: stage tracks accumulated', flush=True)  # ボード単位の進捗を示す。
            all_array = np.array(sorted(all_events), np.int64)  # 通過集合と実際のTPC記録を照合する。
            counts = [len(all_events)] + [int(np.isin(all_array, values).sum()) for values in selected[1:]]  # 同じ分母で各段階のイベント数を出す。
            summary = dict(run=run, all_tpc_events=len(all_events), tpc_events_with_hits=len(with_hits), events=dict(zip(STAGES, counts)), strip_clock_hits=dict(zip(STAGES, map(int, histogram.sum(axis=(1,2,3))))), inputs=manifests)  # 分母とヒット数を分ける。
            assert np.all(histogram[2] <= histogram[1]) and np.all(histogram[1] <= histogram[0])  # 累積条件の包含関係を検証する。
            assert np.all(histogram[3]+histogram[4] <= histogram[2])  # raw1とraw2が重ならないことを確認する。
            np.savez_compressed(target, histogram=histogram, summary=json.dumps(summary), config=json.dumps(config))  # 正常完了してから保存する。
        with np.load(target) as saved:  # 完了済みのランを合算する。
            if json.loads(str(saved['config'])) != config:  # 条件違いの再利用を防ぐ。
                raise ValueError(f'Config changed: {target}')  # 同じ名前で混ぜない。
            combined += saved['histogram']  # ランの投影だけを合算する。
            summaries.append(json.loads(str(saved['summary'])))  # 分母を記録する。
        print(f'COMPLETE stage tracks {run}', flush=True)  # ラン完了を表示する。
    plot_config = dict(config, color_percentiles=[10, 99.5])  # 空ビンを除く百分位でカラーを調整する。
    for index, group in enumerate(GROUPS):  # 段階ごとに独立したフォルダへ描く。
        directory = root / group  # シンチ図と同じ段階のフォルダを使う。
        directory.mkdir(parents=True, exist_ok=True)  # 保存先を作る。
        colors = draw(combined[index], plot_config, directory/'tpc_tracks', f"{config.get('data_label', ', '.join(config['runs']))} {STAGES[index]}", selection_label=['Completed raw; no event cuts; no ToT correction', 'Beam-axis distance cut; no Z veto; no ToT correction', 'Beam-axis THEN Z-side veto; no HA multiplicity cut; no ToT correction', 'Beam-axis THEN Z-side veto THEN exactly 1 raw HA hit; no ToT correction', 'Beam-axis THEN Z-side veto THEN exactly 2 raw HA hits; no ToT correction'][index])  # 共通関数で4種類の図を作る。
        (directory/'tpc_tracks.json').write_text(json.dumps(dict(stage=STAGES[index], selection='Beam-axis THEN Z-side veto THEN raw multiplicity; no gain-specific pixel validity cut; completed raw only', colors=colors, runs=summaries), indent=2)+'\n')  # 図の条件を保存する。
    np.savez_compressed(cache_dir/'all_runs.npz', histogram=combined, stages=np.array(STAGES), config=json.dumps(config))  # SSDなしの再描画配列を残す。
    write_cutflow(root/'stage_track_cutflow.json', dict(stages=STAGES,runs=summaries), 'tracks')  # TPCの分母を明示して保存する。


if __name__ == '__main__':  # importだけでrawを読まない。
    main()  # コマンドから呼ばれた場合だけ処理する。
