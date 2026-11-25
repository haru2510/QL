#!/usr/bin/env python3
"""最終raw1/2段階を、rawフラグのみの位置図（未決定はch中心）とHG/LG品質選別後のピクセル図に分ける。"""
from cutflow_report import write_cutflow  # 全run合計と日本語の説明を先頭に保存する。
import argparse  # 全解析と同じJSONを指定する。
import json  # 分母と品質選別後の件数を記録する。
from pathlib import Path  # 段階別フォルダを作る。
import numpy as np  # rawフラグのイベント数と位置を対応づける。
from settings import load_config, repo_path  # 共通の入出力を使う。
from pixel_polarization import rows, vertex_table, fitted_axis, select_events, pixel_positions  # 保存ビーム軸とイベント選別を共用する。
from plot_hit_pixel_maps import draw_pixel_map  # turboと半透明ビーム印を共用する。


def main():  # raw数とgain別の品質選別を独立した最終段階として描く。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行時の指定を用意する。
    parser.add_argument('--config', type=Path, required=True)  # 入力と出力を外部JSONで指定する。
    args = parser.parse_args()  # 指定を読む。
    config = load_config(args.config)  # 共通条件を読む。
    root = repo_path(config['outputs']['stage_pixel_map_root'])  # 段階別図の親を指定する。
    vertices = vertex_table(repo_path(config['inputs']['vertices_all']))  # Z側選別前の全頂点を使う。
    point, direction = fitted_axis(repo_path(config['inputs']['polarization_beam_axis']))  # 保存済みビーム軸を読む。
    geometry = {(int(row['board']), int(row['ch'])): row for row in rows(repo_path(config['inputs']['scintillator_geometry']))}  # ピクセル配置と位置未決定時のch中心を読む。
    calibration = repo_path(config['inputs']['ha_position_calibration_dir'])  # HGとLGの位置推定に使うHG較正を読む。
    options = dict(config['polarization'], beam_vertex_cut_enabled=True, z0_track_veto_enabled=True, require_exactly_two_raw_hits=False)  # ビーム軸→Z側の順で選別する。
    collected, summaries = {n: {g: [] for g in ('hg','lg')} for n in (1,2)}, []  # raw1/2のブロック位置とラン別件数を集める。
    for run in config['runs']:  # ランを混ぜずにイベントを照合する。
        with np.load(repo_path(config['outputs']['polarization_hit_cache_dir']) / f'{run}.npz') as cache:  # SSDを読み直さずrawフラグを読む。
            events, inverse, counts = np.unique(cache['event'], return_inverse=True, return_counts=True)  # HG/LGの品質選別前にrawヒット数を確定する。
            with np.load(repo_path(config['outputs']['veto_dir']) / f'z0_track_veto_{run}.npz') as veto:  # 全TPCヒットのZ側判定を読む。
                allowed, _ = select_events(events, counts, vertices.get(run, np.zeros((0,4))), point, direction, options, veto['event'])  # ビーム軸条件の後でZ側除外する。
            board, channel = cache['board'], cache['ch']  # rawフラグから分かるブロック番号を読む。
            record = {'run': run, 'raw_selected_events': {}, 'position_assignment': {}}  # ピクセル品質条件より前の分母を保存する。
            for multiplicity in (1,2):  # raw1とraw2は同じ条件から分岐する。
                selected = allowed & (counts == multiplicity)  # rawフラグ数だけで選ぶ。
                indices = np.flatnonzero(selected[inverse])  # 選択イベントの全rawヒットを取り出す。
                hits = {key: cache[key][indices] for key in cache.files if key != 'metadata'}  # raw1/2段階の全フラグを位置推定へ渡す。
                record['position_assignment'][str(multiplicity)] = {}  # gainごとにピクセルとch中心の内訳を残す。
                for gain in ('hg','lg'):  # HG/LGは異なる重心を使う。
                    product = pixel_positions(hits, gain, geometry, calibration, reject_saturation=False)  # 飽和フラグでは落とさず位置を試算する。
                    positions = product['position'].copy()  # 求められたピクセル位置を使う。
                    fallback = ~product['valid']  # ピクセル位置を決められなかったヒットを特定する。
                    for index in np.flatnonzero(fallback):  # 未決定のヒットもraw段階では落とさない。
                        row = geometry[(int(hits['board'][index]),int(hits['ch'][index]))]  # そのchの既知の配置を読む。
                        positions[index] = [float(row[key]) for key in ('x0','y0','z0')]  # ch中心を代替位置に使う。
                    assert len(positions) == int(selected.sum()) * multiplicity and np.isfinite(positions).all()  # 全rawフラグが有限な位置で残っていることを確認する。
                    collected[multiplicity][gain].append(positions)  # rawフラグ図もHG/LGで分ける。
                    record['position_assignment'][str(multiplicity)][gain] = {'pixel_hits': int(product['valid'].sum()), 'channel_centre_hits': int(fallback.sum()), 'saturated_raw_hits': int(np.count_nonzero(hits['satu_h' if gain == 'hg' else 'satu_l']))}  # 未決定・飽和を含む内訳を保存する。
                record['raw_selected_events'][str(multiplicity)] = int(selected.sum())  # 元のイベント数を保存する。
            summaries.append(record)  # ラン別の分母を残す。
    report = {'cut_order': ['beam_axis','z0_track_veto','raw_multiplicity','gain_specific_quality'], 'raw_position': 'attempt gain-specific pixels WITHOUT saturation rejection; unresolved hits at channel centres; ALL raw flags retained', 'quality_selection': 'gain-specific saturation, positive integral, finite centroid, calibration and geometry present, pixel inside 8x8; raw2 additionally requires both valid pixels and defined pair azimuth', 'runs': summaries, 'stage03': 'raw multiplicity ONLY, no gain quality rejection', 'stage04': 'gain-specific quality AFTER raw multiplicity', 'branches': {}}  # 図の違いを明記する。
    for multiplicity, folder in ((1,'03_raw_single_hit'),(2,'03_raw_two_hits')):  # 最終段階の直下を二つに分ける。
        parent = root / 'with_z0_cut/stages' / folder  # 元の最終図がある場所を指定する。
        before, after = parent, parent.parent / (folder.replace('03_raw_', '04_') + '_quality')  # フラグのみと品質選別後を別フォルダにする。
        before.mkdir(parents=True, exist_ok=True)  # rawフラグ段階を作る。
        after.mkdir(parents=True, exist_ok=True)  # 品質選別段階を作る。
        events = sum(row['raw_selected_events'][str(multiplicity)] for row in summaries)  # raw段階のイベント数を求める。
        branch = {'raw_events': events, 'raw_hit_flags': events*multiplicity, 'gains': {}}  # rawはHG/LG共通の分母にする。
        for gain in ('hg','lg'):  # 両gainで前後を確認できるようにする。
            points = np.vstack(collected[multiplicity][gain])  # このgainのピクセルとch中心代替点を合算する。
            fallback_count = sum(row['position_assignment'][str(multiplicity)][gain]['channel_centre_hits'] for row in summaries)  # ch中心へ代替した点数を合計する。
            title = f'{gain.upper()}: exactly {multiplicity} raw HA hit flags; {events:,} events / {len(points):,} flags\nBeam-axis THEN Z-side veto; BEFORE gain quality cuts; pixels where assignable, otherwise channel centres ({fallback_count:,} fallbacks)'  # 未決定位置だけch中心に置き、rawフラグを全て残したことを明記する。
            raw_plot = draw_pixel_map(points, config, title, 'raw hit flags / projected position', before/f'{gain}.png')  # 推定ピクセルと未決定のch中心で全rawフラグを描く。
            with np.load(after/f'{gain}.points.npz') as saved:  # 既に作成した最終段階の採用ピクセルを読む。
                valid_points = saved['positions']  # 同じ品質選別結果をそのまま使う。
            count = len(valid_points)//multiplicity  # 品質条件後の採用イベント数を求める。
            title = f'{gain.upper()}: exactly {multiplicity} raw HA hit flags THEN gain quality cuts; {count:,} events / {len(valid_points):,} hit uses\nBeam-axis THEN Z-side veto; assigned pixels; HG calibration (LG provisional)'  # 品質選別後はピクセル図であることを示す。
            draw_pixel_map(valid_points, config, title, 'accepted hit uses / projected pixel', after/f'{gain}.png', color_max=raw_plot['color_max'])  # 飽和などを除いた位置図を別に保存する。
            assert count <= events  # 品質選別後にイベントが増えていないことを確認する。
            branch['gains'][gain] = {'raw_pixel_hits': len(points)-fallback_count, 'raw_channel_centre_hits': fallback_count, 'after_quality_events': count, 'after_quality_hit_uses': len(valid_points), 'rejected_events': events-count}  # 飽和だけでなく全品質条件の減少数を保存する。
        report['branches'][str(multiplicity)] = branch  # raw1/2を別に記録する。
        for directory in (before,after):  # それぞれのフォルダ単独でも条件を確認できるようにする。
            (directory/'selection.json').write_text(json.dumps(report,indent=2)+'\n')
            (directory/'cutflow.json').write_text(json.dumps(dict(stage=directory.name, selection=report['stage03'] if directory == before else report['stage04'], counts=branch, position_assignment=report['raw_position'] if directory == before else 'assigned pixels only'),indent=2)+'\n')  # 各stageの条件と件数を独立して残す。  # 条件と件数を保存する。
    write_cutflow(root/'stage_flag_quality_cutflow.json', report, 'quality')  # 最終品質条件の合計と詳細を保存する。
    print(json.dumps(report['branches'],indent=2),flush=True)  # rawフラグと品質選別の違いを確認する。


if __name__ == '__main__':  # importだけでは図を作らない。
    main()  # コマンドからの実行時だけ処理する。
