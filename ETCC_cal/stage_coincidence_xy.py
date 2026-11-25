#!/usr/bin/env python3
"""同イベント・同uPIC・同クロックのX×Y全組合せを積算し、TPCのXZ/YZの横に並べる。"""
import argparse  # 共通JSONと再描画を指定する。
import heapq  # ボードごとのフレームをイベント順に結合する。
import json  # 条件と件数を保存する。
from pathlib import Path  # rawと積算配列を扱う。
import numpy as np  # 3.2 mmのXYビンへ積算する。
from numba import njit  # rawビット走査とコインシデンスを高速化する。
from settings import load_config, repo_path  # 既存と同じ設定を使う。
from track_sum import decoder  # 検証済みTPCフレーム読み取りを使う。
from track_runs import draw  # 既存XZ/YZへXYパネルを追加する。
from stage_track_sum import STAGES, GROUPS, selection_sets, accumulate_stages  # 同じカットと保存先を使う。
from pixel_polarization import vertex_table, fitted_axis  # 保存された散乱点とビーム軸を読む。


@njit  # 1フレームのrawストリップを、クロック・XY表示ビンごとの8bitへ保存する。
def fill_occupancy(buf, offset, records, axis, strip_bins, occupancy, occupied_bins, bin_counts):  # 重複した同一ストリップはORで1回にまとめる。
    first = offset + decoder.HEAD_SIZE + decoder.FADC_SIZE + 4  # 最初のレコード位置を求める。
    for record in range(records):  # 全クロックを使い、早い16個には限定しない。
        start = first + decoder.RECORD_SIZE * record  # 現在のレコード位置を求める。
        clock = decoder.read_word(buf,start) & 0x7ff  # クロックを読む。
        if clock >= 1024:  # 不正クロックを除く。
            continue  # 有効範囲だけを使う。
        for channel in range(4):  # 128ストリップの4ワードを読む。
            word = decoder.read_word(buf,start+4*(4-channel))  # raw仕様に沿う順番を使う。
            for bit in range(32):  # ヒットビットを調べる。
                if (word >> bit) & 1:  # 有効ストリップを採用する。
                    coordinate = strip_bins[channel*32+bit]  # 0.4 mmの全体座標ビンを得る。
                    coarse, fine = coordinate//8, coordinate%8  # 3.2 mm表示ビンと内部の位置へ分ける。
                    if occupancy[axis,clock,coarse] == 0:  # 初めて使う表示ビンを記録する。
                        occupied_bins[axis,clock,bin_counts[axis,clock]] = coarse  # 後で空ビンを走査しないようにする。
                        bin_counts[axis,clock] += 1  # このクロックの使用ビン数を増やす。
                    occupancy[axis,clock,coarse] |= np.uint8(1 << fine)  # 同じストリップの重複レコードを二重計数しない。


@njit  # 同一イベントの全ボードが揃った時点でXYへ足す。
def commit_event(occupancy, occupied_bins, bin_counts, popcount, mask, histogram):  # 同じuPIC内だけでX/Yを組にする。
    total = 0  # このuPIC・イベントの同時ヒット組合せ数を記録する。
    for clock in range(1024):  # クロック差0の条件を厳密に使う。
        nx, ny = bin_counts[0,clock], bin_counts[1,clock]  # X/Yにヒットがあるか調べる。
        for ix in range(nx):  # X側の使用表示ビンだけを調べる。
            x = occupied_bins[0,clock,ix]  # 表示ビンを得る。
            for iy in range(ny):  # 同じクロックのY側だけを組み合わせる。
                y = occupied_bins[1,clock,iy]  # Y表示ビンを得る。
                combinations = int(popcount[occupancy[0,clock,x]]) * int(popcount[occupancy[1,clock,y]])  # 表示ビンへまとめる前の全strip組合せ数を数える。
                total += combinations  # 3D候補点数を加算する。
                for stage in range(histogram.shape[0]):  # 同じイベント条件を各段階へ適用する。
                    if mask[stage]:  # 通過した段階へだけ積算する。
                        histogram[stage,x,y] += combinations  # XYへ投影して足す。
        for axis in range(2):  # 使用したビンだけ初期化してメモリの書込みを抑える。
            for index in range(bin_counts[axis,clock]):  # 空ビンを触らない。
                occupancy[axis,clock,occupied_bins[axis,clock,index]] = 0  # 次のイベントへ持ち越さない。
            bin_counts[axis,clock] = 0  # ビン数も初期化する。
    return total  # 同時ヒットありイベントを数えるために返す。


def frame_stream(paths, tag, allowed_events=None):  # 同じボードのフレームをepoch付きイベント順に供給する。
    previous, epoch, seen = None, 0, set()  # 巻き戻りと重複を管理する。
    for path in paths:  # 大量のrawを1ファイルずつ読む。
        buf, triggers, offsets, records = decoder.scan_file(path)  # バッファと有効フレーム位置を読む。
        if not len(triggers):  # 空ファイルは飛ばす。
            continue  # 次のファイルへ進む。
        back = np.r_[previous is not None and triggers[0] < previous, triggers[1:] < triggers[:-1]]  # 従来と同じ巻き戻り条件を使う。
        epochs = epoch + np.cumsum(back)  # epochをフレームごとに付ける。
        epoch, previous = int(epochs[-1]), int(triggers[-1])  # 次ファイルへ状態を渡す。
        events = epochs*(1<<32)+triggers  # 他ボードと照合できるキーを作る。
        indices = np.flatnonzero(np.isin(events,allowed_events)) if allowed_events is not None else range(len(events))  # stage04では採用イベントだけを供給し、他フレームのビット走査を省く。
        for i in indices:  # 同じイベントを持つ他ボードとストリーム結合する。
            event = int(events[i])  # 集合のキーをPython整数にする。
            if event not in seen:  # 初出フレームだけを使う。
                seen.add(event)  # 次の重複を除くため記録する。
                yield event, buf, int(offsets[i]), int(records[i])  # バッファはこのファイルを使い切るまで保持する。
        print(tag,Path(path).name,'coincidence frames read',flush=True)  # 長いuPIC処理の進捗を残す。


def scan_run(run, config, selected, raw_histogram=None):  # 同じuPICの最大6ボードをイベントごとに結合する。
    directories = [str(repo_path(config['data_root'])/storage/run) for storage in config['storages']]  # 両storageを読む。
    boards = decoder.collect_boards(directories,0)  # 完成rawだけを使う。
    if not boards:  # 空入力を正常終了にしない。
        raise FileNotFoundError(run)  # ラン名を表示する。
    histogram = np.zeros((len(selected),96,96),np.int64)  # 各段階のXYを確保する。
    stage_sets = [None if events is None else set(map(int,events)) for events in selected]  # イベント条件を高速検索できるようにする。
    allowed_events = None if any(events is None for events in selected) else np.unique(np.concatenate(selected))  # 全イベント段階がなければ採用フレームだけ読み出す。
    coincident_events = set()  # 複数uPICに現れる同じイベントを1回だけ数える。
    popcount = np.array([bin(i).count('1') for i in range(256)],np.uint8)  # 各8bitのstrip数を事前計算する。
    for upic in sorted({board//decoder.N_BOARDS_PER_UPIC for _,board in boards}):  # 異なるuPIC間の偽組合せを作らない。
        occupancy = np.zeros((2,1024,96),np.uint8)  # 1イベント分のstrip占有だけ保持する。
        occupied_bins = np.empty((2,1024,96),np.int16)  # 非空表示ビンのリストを保持する。
        bin_counts = np.zeros((2,1024),np.int16)  # クロックごとの非空ビン数を保持する。
        streams, layouts, heap = [], [], []  # 各ボードの先頭フレームだけを保持する。
        for (side,board), paths in boards.items():  # 同じuPICのボードを集める。
            if board//decoder.N_BOARDS_PER_UPIC != upic:  # 別uPICは次で処理する。
                continue  # 混ぜずに飛ばす。
            axis, sign, const = decoder.STRIP_AXES[upic][0 if side=='a' else 1]  # ストリップの向きを読む。
            mm = config['pitch_strip']*(sign*(np.arange(128)+(board%decoder.N_BOARDS_PER_UPIC)*128)+const)  # 中心座標へ変換する。
            positions = np.rint((mm+153.4)/0.4).astype(np.int64)  # 既存XZ/YZと同じ0.4 mmビンへ変換する。
            stream = frame_stream(paths,f'{run} {side}{board:02d}',allowed_events)  # ボードのストリームを作る。
            index = len(streams)  # heap内でボードを区別する。
            streams.append(stream)  # generatorを保持する。
            layouts.append((int(axis=='y'),positions))  # X/Yの向きと位置を保存する。
            frame = next(stream,None)  # 最初のフレームを読む。
            if frame is not None:  # 空ボードはheapへ入れない。
                heapq.heappush(heap,(frame[0],index,frame))  # 最も早いイベントを選べるようにする。
        while heap:  # 同じイベントのフレームを全ボードから集める。
            event = heap[0][0]  # 現在の最小イベント番号を得る。
            mask = np.array([values is None or event in values for values in stage_sets],bool)  # 5段階のイベント条件を確定する。
            while heap and heap[0][0] == event:  # 同じイベントだけをoccupancyへまとめる。
                _, index, frame = heapq.heappop(heap)  # 次のボードのフレームを読む。
                axis, positions = layouts[index]  # そのボードの配置を読む。
                fill_occupancy(frame[1],frame[2],frame[3],axis,positions,occupancy,occupied_bins,bin_counts)  # rawストリップを共通の座標へ置く。
                if raw_histogram is not None:  # stage04では同じフレームからXZ/YZも同時に積算する。
                    from_stage = mask.reshape(len(mask),1)  # 1フレームのgain別通過条件へ直す。
                    accumulate_stages(frame[1],np.array([frame[2]],np.int64),np.array([frame[3]],np.int32),np.ones(1,np.bool_),from_stage,axis,positions,raw_histogram)  # rawとXYのためにSSDを二重走査しない。
                following = next(streams[index],None)  # 同じボードの次イベントを準備する。
                if following is not None:  # ストリーム終端ならheapへ戻さない。
                    heapq.heappush(heap,(following[0],index,following))  # 次のイベント順に並べる。
            if commit_event(occupancy,occupied_bins,bin_counts,popcount,mask,histogram):  # 3D同時ヒットを全組合せで積算する。
                coincident_events.add(event)  # 同時ヒットありのイベントを記録する。
        print(f'COMPLETE {run} uPIC{upic} XY coincidence',flush=True)  # uPICごとの完了を残す。
    if selected[0] is None and len(selected)==5:  # stage00〜03の累積条件にだけ包含関係を適用する。
        assert np.all(histogram[2]<=histogram[1]) and np.all(histogram[1]<=histogram[0])  # 累積カットを検証する。
        assert np.all(histogram[3]+histogram[4]<=histogram[2])  # raw1とraw2が交わらないことを検証する。
    return histogram, len(coincident_events)  # XY配列と同時ヒットありイベント数を返す。


def main():  # 同じカット条件のXZ/YZ/XYを作る。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行条件を用意する。
    parser.add_argument('--config',required=True,type=Path)  # 共通外部JSONを指定する。
    parser.add_argument('--plot-only',action='store_true')  # XY保存済みならSSDを使わず描く。
    args = parser.parse_args()  # 引数を読む。
    config = load_config(args.config)  # 入出力と条件を読む。
    root = repo_path(config['result_dir'])  # 結果保存先を求める。
    cache_dir = root/'shared/stage_coincidence_xy'  # XY専用の配列を保存する。
    cache_dir.mkdir(parents=True,exist_ok=True)  # 保存先を作る。
    vertices = vertex_table(repo_path(config['inputs']['vertices_all']))  # 保存全頂点を読む。
    point,direction = fitted_axis(repo_path(config['inputs']['polarization_beam_axis']))  # 同じビーム軸を使う。
    combined, summaries = np.zeros((5,96,96),np.int64), []  # 全ラン合計を用意する。
    for run in config['runs']:  # ランごとに結合する。
        target = cache_dir/f'{run}.npz'  # ラン別の再開用保存先を決める。
        if not target.exists() and not args.plot_only:  # 未完了ランだけrawを読む。
            with np.load(root/'shared/stage_tracks'/f'{run}.npz') as source:  # XZ/YZと同じ入力ファイルか確認する。
                original = json.loads(str(source['summary']))['inputs']  # 元のraw走査時の入力を読む。
            for entry in original:  # 再マウントやコピーで変わったデータを混ぜない。
                stat = Path(entry['path']).stat()  # 現在の入力情報を読む。
                if [stat.st_size,stat.st_mtime_ns] != [entry['bytes'],entry['mtime_ns']]:  # 同一性を確認する。
                    raise ValueError(f'Raw changed since XZ/YZ: {entry["path"]}')  # 不一致で停止する。
            xy,events = scan_run(run,config,selection_sets(run,config,vertices,point,direction))  # 全クロックのコインシデンスを取る。
            summary = dict(run=run,events_with_3d_coincidence=events,xy_combinations=[int(h.sum()) for h in xy])  # 点数とイベント数を分ける。
            np.savez_compressed(target,xy=xy,summary=json.dumps(summary),config=json.dumps(config))  # 完了したランだけ保存する。
        with np.load(target) as saved:  # 保存済みを合算する。
            if json.loads(str(saved['config'])) != config:  # 別条件を混ぜない。
                raise ValueError(f'Config mismatch: {target}')  # 条件違いを明示する。
            combined += saved['xy']  # ラン間では積算配列だけ加算する。
            summaries.append(json.loads(str(saved['summary'])))  # 件数を保存する。
    with np.load(root/'shared/stage_tracks/all_runs.npz') as raw:  # XZ/YZは既存の全rawヒット積算を維持する。
        for index,group in enumerate(GROUPS):  # 同じstageへ3パネルを保存する。
            prefix = root/group/'tpc_tracks'  # 4種類とも既存名を更新する。
            draw(raw['histogram'][index],dict(config,color_percentiles=[10,99.5]),prefix,f"{config.get('data_label', ', '.join(config['runs']))} {STAGES[index]}",selection_label='XZ/YZ: raw strip hits; XY: same-event, same-uPIC, same-clock X-Y combinations; no ToT correction',xy_histogram=combined[index])  # XYだけが3Dコインシデンス後であることを明記する。
    definition = dict(clock_tolerance=0,unit='all unique strip X-Y combinations at identical event, uPIC and clock; may include ambiguous pairs',runs=summaries)  # 3D候補点の定義を保存する。
    np.savez_compressed(cache_dir/'all_runs.npz',xy=combined,config=json.dumps(config))  # SSDなしの再描画用に保存する。
    (root/'stage_coincidence_xy.json').write_text(json.dumps(definition,indent=2)+'\n')  # 件数と定義を保存する。
    print('ALL XY COINCIDENCE COMPLETE',root,flush=True)  # 全処理の完了を知らせる。


if __name__=='__main__':  # importだけではrawを読まない。
    main()  # コマンドから呼ばれた場合だけ実行する。
