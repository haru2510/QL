#!/usr/bin/env python3
"""完成rawから段階別解析の全頂点・Z側除外・HG/LG重心キャッシュを作る。"""
import argparse  # 外部設定と実行範囲を指定する。
import csv  # コメント付き頂点CSVを保存する。
import json  # 入力一覧と条件を記録する。
import shutil  # 既存run06・07の計算結果を再利用する。
import time  # 実測の処理時間を表示する。
from pathlib import Path  # 入出力を扱う。
import numpy as np  # ファイル単位で行をまとめる。
from settings import load_config, repo_path  # 共通JSONのパス規則を使う。
from track_sum import decoder  # 検証済みのTPC読み取りを使う。
from pixel_polarization import legacy  # HG/LG重心を同時に読み取る。

FIELDS = ['run', 'epoch', 'trig', 'upic', 'clock', 'x', 'y', 'z', 'n_hits', 'fiducial', 'dist_to_line']  # 既存CSVと同じ列にする。


def reconstruct(run, config, destination):  # rawを一度読み、局所ファイルへ分割してメモリを抑える。
    directory = destination / run  # ラン別の中間保存先を決める。
    directory.mkdir(parents=True, exist_ok=True)  # 初回にフォルダを作る。
    final = directory / 'vertices.csv'  # 完了済み頂点の保存名を決める。
    marker = directory / 'complete.json'  # 全raw走査の完了を区別する。
    if marker.exists():  # 正常完了したランだけ再利用する。
        saved = json.loads(marker.read_text())  # 完了時の入力一覧を読む。
        for item in saved['inputs']:  # 入力が変わっていないか照合する。
            stat = Path(item['path']).stat()  # 現在のファイル情報を取得する。
            if [stat.st_size, stat.st_mtime_ns] != [item['bytes'], item['mtime_ns']]:  # 不一致を黙って流用しない。
                raise ValueError(f'Input changed: {item["path"]}')  # 手動で作り直す必要があることを示す。
        return final  # 保存済み頂点を使う。
    boards = decoder.collect_boards([str(repo_path(config['data_root']) / storage / run) for storage in config['storages']], 0)  # 完成rawを全て収集する。
    if not boards:  # 空入力を成功扱いしない。
        raise FileNotFoundError(run)  # 対象ランを表示する。
    size = config['preparation']['event_chunk_size']  # 結合するイベント番号幅を読む。
    keep_clocks = config['preparation']['keep_clocks']  # 既存頂点と同じ早い16レコードを保持する。
    clock_limit = int(np.floor((config['z_veto'] - config['z_gem']) / config['velocity'] + config['clock_offset']))  # Z側除外を整数クロックへ変換する。
    chunks, veto, all_events, manifests = {}, set(), set(), []  # 結合用ファイルと全ボード共通の判定を集める。
    serial = 0  # 中間ファイルの重複名を防ぐ。
    for (side, board), paths in boards.items():  # 同じボード内でepochを引き継ぐ。
        epoch, previous, seen = 0, None, set()  # トリガー巻き戻りと重複を管理する。
        upic, offset = divmod(board, decoder.N_BOARDS_PER_UPIC)  # uPIC内の配置を求める。
        axis, sign, const = decoder.STRIP_AXES[upic][0 if side == 'a' else 1]  # 配置と向きを読む。
        position = config['pitch_strip'] * (sign * (np.arange(decoder.N_STRIPS) + offset * decoder.N_STRIPS) + const)  # 各ストリップのmm座標を作る。
        for path in paths:  # SSD上のrawは1ファイルずつ読む。
            stat = Path(path).stat()  # 入力の同一性を記録する。
            manifests.append(dict(path=path, bytes=stat.st_size, mtime_ns=stat.st_mtime_ns))  # 将来の再利用に備える。
            buf, triggers, offsets, records = decoder.scan_file(path)  # 有効フレームの位置を読む。
            if not len(triggers):  # 空ファイルは記録だけ残す。
                continue  # 次へ進む。
            back = np.r_[previous is not None and triggers[0] < previous, triggers[1:] < triggers[:-1]]  # トリガー巻き戻りを判定する。
            epochs = epoch + np.cumsum(back)  # ボード内でepochを増やす。
            epoch, previous = int(epochs[-1]), int(triggers[-1])  # ファイル境界へ状態を引き継ぐ。
            events = epochs * (1 << 32) + triggers  # イベントキーを作る。
            unique, first = np.unique(events, return_index=True)  # 同一ファイル内の重複も除く。
            keep = np.zeros(len(events), bool)  # 初出フレームだけを残す。
            keep[first] = ~np.isin(unique, np.fromiter(seen, dtype=np.int64))  # 既読フレームを除く。
            seen.update(map(int, events[keep]))  # 次のファイルへ記録を引き継ぐ。
            all_events.update(map(int, events[keep]))  # 全TPCイベントの和集合を作る。
            capacity = int(np.minimum(records, keep_clocks).sum())  # 保持するレコード数だけ確保する。
            frame = np.empty(capacity, np.int64)  # 行の元フレームを保存する。
            clock = np.empty(capacity, np.int64)  # 行のクロックを保存する。
            count = np.empty(capacity, np.int64)  # ヒット数を保存する。
            total = np.empty(capacity, np.float64)  # 座標和を保存する。
            hits = np.zeros(len(events), np.int64)  # 全クロックのヒット数を保存する。
            length = decoder.fill_rows(buf, offsets, records, keep_clocks, position, frame, clock, count, total, hits)  # 従来と同じ頂点入力を作る。
            good = keep[frame[:length]]  # 初出フレームの行だけ残す。
            row_event = events[frame[:length][good]]  # 行のイベント番号を対応づける。
            row_clock, row_count, row_sum = clock[:length][good], count[:length][good], total[:length][good]  # 同じ行を切り出す。
            veto.update(map(int, row_event[row_clock <= clock_limit]))  # 最小クロックは必ず保持されるため全ヒットのZ側判定と同じになる。
            bucket_rows = row_event // size  # メモリ上の結合サイズを制限する。
            bucket_frames = events[keep] // size  # ヒットなしイベントも分母へ含める。
            for bucket in np.unique(bucket_frames):  # このファイルに含まれるイベント番号範囲へ分割する。
                use, frame_use = bucket_rows == bucket, bucket_frames == bucket  # 同じ範囲の行とフレームを選ぶ。
                output = directory / f'part_{serial:06d}.npz'  # ローカル中間ファイルの名を作る。
                np.savez(output, event=row_event[use], clock=row_clock[use], count=row_count[use], total=row_sum[use], upic=upic, side=int(side == 'c'), axis=int(axis == 'y'), frame_event=events[keep][frame_use], frame_hits=hits[keep][frame_use])  # 高速な非圧縮保存を使う。
                chunks.setdefault(int(bucket), []).append(output)  # 後で同じイベント範囲だけ結合する。
                serial += 1  # 次の保存名へ進む。
        print(f'{run} {side}{board:02d}: {len(seen):,} frames prepared', flush=True)  # ボード単位の進捗を知らせる。
    vertex_count = 0  # 頂点を持つイベント数を記録する。
    with final.open('w') as stream:  # ラン別CSVを作る。
        stream.write('# vertex = earliest coincident clock in the same uPIC; coordinates [mm]\n')  # 座標の意味を明記する。
        writer = csv.writer(stream)  # 表の書き出しを用意する。
        writer.writerow(FIELDS)  # 既存の列名を保存する。
        for bucket, files in sorted(chunks.items()):  # 最大10万イベントずつ結合する。
            columns = {name: [] for name in ('event', 'clock', 'count', 'total', 'upic', 'side', 'axis', 'frame_event', 'frame_hits')}  # 分割配列を集める。
            for path in files:  # ローカル中間ファイルを読む。
                with np.load(path) as part:  # 読み終えたファイルを閉じる。
                    for name in columns:  # 行情報とフレーム情報を集める。
                        value = part[name]  # 値を取り出す。
                        columns[name].append(np.full(len(part['event']), int(value), dtype=np.int64) if value.ndim == 0 else value)  # 配置は各行へ展開する。
            h = {name: np.concatenate(parts) for name, parts in columns.items()}  # この範囲だけメモリに載せる。
            key = h['event'] * 4096 + h['upic'] * 1024 + h['clock']  # (event,uPIC,clock)を整数キーにする。
            keys, inverse = np.unique(key, return_inverse=True)  # 同じ位置の行をまとめる。
            present = np.zeros((len(keys), 2), bool)  # アノード・カソードの有無を調べる。
            present[inverse, h['side']] = True  # 片側だけの点を除く準備をする。
            sums, numbers = np.zeros((len(keys), 2)), np.zeros((len(keys), 2))  # X/Yの座標和と数を作る。
            np.add.at(sums, (inverse, h['axis']), h['total'])  # 同じクロックの座標を加算する。
            np.add.at(numbers, (inverse, h['axis']), h['count'])  # 座標平均の分母を加算する。
            both = np.flatnonzero(present.all(axis=1))  # 両側のある候補を得る。
            candidate = keys[both]  # 候補のイベントキーを得る。
            order = np.lexsort((candidate % 1024, candidate // 4096))  # イベントごとに最も早いクロックを選ぶ。
            both, candidate = both[order], candidate[order]  # 候補と集計位置を揃える。
            first = np.r_[True, np.diff(candidate // 4096) != 0] if len(candidate) else np.zeros(0, bool)  # 各イベントの先頭を選ぶ。
            both, candidate = both[first], candidate[first]  # 頂点1点に絞る。
            frame_keys, frame_inv = np.unique(h['frame_event'], return_inverse=True)  # 全ボードのイベントを合わせる。
            frame_hits = np.bincount(frame_inv, weights=h['frame_hits']).astype(np.int64)  # 全クロックのヒット数を合計する。
            for index, value in zip(both, candidate):  # 頂点をCSVへ逐次保存する。
                event, remainder = divmod(int(value), 4096)  # イベントキーを復元する。
                upic, clock = divmod(remainder, 1024)  # uPICとクロックを復元する。
                x, y = sums[index] / np.maximum(numbers[index], 1)  # 座標を平均する。
                z = config['velocity'] * (clock - config['clock_offset']) + config['z_gem']  # Zへ換算する。
                writer.writerow([run, event >> 32, event & 0xffffffff, upic, clock, x, y, z, int(frame_hits[np.searchsorted(frame_keys, event)]), int(event not in veto), 'nan'])  # XYや高Z側のfiducial条件を入れない。
            vertex_count += len(candidate)  # 件数を加算する。
            print(f'{run} chunk {bucket}: {len(candidate):,} vertices', flush=True)  # 結合の進捗を表示する。
    veto_dir = repo_path(config['outputs']['veto_dir'])  # Z側判定の共有保存先を求める。
    veto_dir.mkdir(parents=True, exist_ok=True)  # 保存先を作る。
    np.savez_compressed(veto_dir / f'z0_track_veto_{run}.npz', event=np.array(sorted(veto), dtype=np.int64))  # 全TPCヒットに対する除外イベントを保存する。
    np.savez_compressed(directory / 'events.npz', event=np.array(sorted(all_events), dtype=np.int64))  # 全TPCイベントの分母を保存する。
    marker.write_text(json.dumps(dict(inputs=manifests, vertices=vertex_count, all_events=len(all_events), veto_events=len(veto), config=config), ensure_ascii=False, indent=2))  # 正常完了を最後に記録する。
    return final  # 完成したCSVを返す。


def main():  # 設定に従い全ランの入力キャッシュを準備する。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行条件を用意する。
    parser.add_argument('--config', type=Path, required=True)  # 入出力を外部JSONで指定する。
    args = parser.parse_args()  # 引数を読む。
    config = load_config(args.config)  # 共通設定を読む。
    root = repo_path(config['result_dir'])  # 結果保存先を求める。
    shared = root / 'shared'  # 中間配列を保存する。
    sources = []  # 合算する頂点CSVを集める。
    old = repo_path('ETCC_cal/output/20261124_run06_07')  # 既存の計算結果を再利用する。
    for run in config['runs']:  # ランを混ぜずに準備する。
        start = time.monotonic()  # 実測時間を記録する。
        if run in config['preparation']['runs_to_reconstruct']:  # run00〜05だけ新規に再構成する。
            sources.append(reconstruct(run, config, shared / 'vertex_chunks'))  # ファイル単位でTPCを読む。
        for folder, filename in [('selection', f'z0_track_veto_{run}.npz'), ('ha_hit_cache', f'{run}.npz')]:  # run06・07の保存キャッシュを流用する。
            target = shared / folder / filename  # 全ラン用の保存先を作る。
            target.parent.mkdir(parents=True, exist_ok=True)  # 保存先を確保する。
            source = old / 'shared' / folder / filename  # 既存ファイルを指定する。
            if not target.exists() and source.exists():  # 新規計算したものは上書きしない。
                shutil.copy2(source, target)  # 既存の完全な結果を再利用する。
        cache = shared / 'ha_hit_cache' / f'{run}.npz'  # HG/LG重心の保存名を指定する。
        if not cache.exists():  # 未計算のランだけHA rawを読む。
            directories = [str(repo_path(config['data_root']) / storage / run) for storage in config['storages']]  # 両storageを読む。
            boards = legacy.collect_files(directories, False, 0)  # 完成rawを全数使う。
            if not boards:  # 入力欠損を明示する。
                raise FileNotFoundError(f'No HA raw: {run}')  # 空キャッシュを保存しない。
            hits, frames, resets = legacy.read_hits(boards, include_centroids=True)  # HG/LG重心を同時に求める。
            np.savez_compressed(cache, **hits, metadata=json.dumps(dict(run=run, max_files_per_board=0, frames=frames, resets=resets, files=boards)))  # 既存解析で再利用できる形式にする。
        print(f'COMPLETE {run}: {time.monotonic()-start:.1f} s', flush=True)  # 実測所要時間を知らせる。
    if any(run in ('20261124_run06', '20261124_run07') and run not in config['preparation']['runs_to_reconstruct'] for run in config['runs']):  # 設定で要求した既存ランだけを再利用する。
        sources.append(old / 'without_z0_cut/vertices/all.csv')  # 円偏光run08〜12へrun06・07を混ぜない。
    destination = repo_path(config['inputs']['vertices_all'])  # 合算CSVの保存先を読む。
    destination.parent.mkdir(parents=True, exist_ok=True)  # 保存先を作る。
    with destination.open('w') as stream:  # 全ランのCSVを作る。
        stream.write('# vertex = earliest coincident clock in the same uPIC; coordinates [mm]\n')  # 既存形式を維持する。
        stream.write(','.join(FIELDS)+'\n')  # 共通の列名を書く。
        for source in sources:  # ラン別の結果を合算する。
            with source.open() as part:  # 巨大CSVを逐次コピーする。
                next(part); next(part)  # コメントと列名を除く。
                for row in part:  # 共用CSVでも要求したランだけを足す。
                    if row.split(',', 1)[0] in config['runs']:  # ラン集合の混入を防ぐ。
                        stream.write(row)  # データ行をコピーする。
    print(f'ALL PREPARED: {destination}', flush=True)  # 正常完了を知らせる。


if __name__ == '__main__':  # importだけではrawを読まない。
    main()  # コマンド指定で実行する。
