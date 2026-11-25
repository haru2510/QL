#!/usr/bin/env python3
"""イベント残存数JSONの先頭へ、日本語の読み方と全run合計を追加する。"""
import json  # 機械向けの既存キーを維持して説明を追加する。
from pathlib import Path  # 保存先から共通の一覧を作る。

LABELS = {'00_before_cuts':'stage00：カット前', '01_beam_axis':'stage01：ビーム軸の距離条件後', '02_after_z0':'stage02：さらにZ側イベント除外後', '03_raw1':'stage03：raw 1ヒットの枝', '03_raw2':'stage03：raw 2ヒットの枝'}  # 1/2ヒットが分岐であることを示す。


def totals(data, kind):  # 元の数値から集計し、表示用の数値を手入力しない。
    if kind in ('scint','tracks'):  # stage00〜03の同じ分母での残存数を集計する。
        return {LABELS[key]:sum(run['events'][key] for run in data['runs']) for key in LABELS}  # ランごとのイベント数を足す。
    if kind=='quality':  # シンチの最終品質条件の数を取り出す。
        return {f'stage04：raw {n}ヒット → {gain.upper()}品質選別後':value['after_quality_events'] for n,branch in data['branches'].items() for gain,value in branch['gains'].items()}  # HGとLGを別々の枝として示す。
    if kind=='quality_tracks':  # stage04のTPC選別数も同じ順で示す。
        return {f'stage04：raw {multiplicity}ヒット → {gain.upper()}品質選別後':sum(run['selected_events'][index] for run in data['runs']) for index,(_,gain,multiplicity) in enumerate(data['branches'])}  # 4条件を独立に合計する。
    raise ValueError(kind)  # 未対応の集計を黙って書かない。


def make_overview(root):  # シンチのstage00〜04を1ファイルで読める入口にする。
    first,quality = root/'stage_cutflow.json', root/'stage_flag_quality_cutflow.json'  # 元の2つの集計を指定する。
    if not first.exists() or not quality.exists():  # まだ前後の集計が揃わない間は作らない。
        return  # 次の生成時に更新する。
    initial,final = json.loads(first.read_text()),json.loads(quality.read_text())  # 保存済み数値から一覧を作る。
    counts = dict(totals(initial,'scint'),**totals(final,'quality'))  # シンチの分母だけをつなぐ。
    content = {'データセット':root.name,'数え方':'イベント数。ヒット点数やペア内の2点の合計ではない。','分母':'stage00はrawシンチヒットが1個以上あるイベント。全TPCトリガーとは異なる。','カットの順番':['stage00 → stage01（ビーム軸）→ stage02（Z側）','stage02からstage03のraw 1ヒットとraw 2ヒットへ分岐','各raw枝からstage04のHG版とLG版へ分岐。HG→LGの連続カットではない。'],'全run合計_読みやすい表示':{key:f'{value:,} イベント' for key,value in counts.items()},'全run合計_数値':counts,'ラン別の詳細':{'stage00〜03':first.name,'stage04':quality.name}}  # 先に読む内容だけに絞る。
    track_path = root/'stage_track_cutflow.json'  # 全TPCを分母にした別集計を明示する。
    if track_path.exists():  # TPC集計がある場合だけ別枠へ載せる。
        track_counts = totals(json.loads(track_path.read_text()),'tracks')  # TPCの分母をシンチの列へ混ぜない。
        content['参考_TPCの別分母'] = {'分母':'全ボードの記録TPCイベントの和集合。ヒットなしを含む。','全run合計':{key:f'{value:,} イベント' for key,value in track_counts.items()},'詳細':track_path.name}  # 読み違いを防ぐ。
    (root/'event_cutflow_readable.json').write_text(json.dumps(content,ensure_ascii=False,indent=2)+'\n')  # 日本語と桁区切りで保存する。


def write_cutflow(path, data, kind):  # 既存形式を壊さず先頭に読みやすい要約を置く。
    path = Path(path)  # 保存先を扱いやすくする。
    prefix = ('_読み方','_集計単位','_全run合計')  # 再整形時に自分の説明だけを置き換える。
    original = {key:value for key,value in data.items() if key not in prefix}  # 元の数値キーは維持する。
    unit = 'stage00の分母はrawシンチヒットが1個以上あるイベント。' if kind=='scint' else 'stage00の分母は全記録TPCイベント（ヒットなしを含む）。' if kind=='tracks' else 'raw 1/2ヒットの枝ごとに、HGとLGの品質条件を独立に適用した採用イベント数。'  # 分母を明示する。
    formatted = {'_読み方':'まず全run合計を見る。stage03のraw1/raw2、stage04のHG/LGは別々の枝。runs以降はラン別の詳細。','_集計単位':'イベント数。'+unit,'_全run合計':{key:f'{value:,} イベント' for key,value in totals(original,kind).items()},**original}  # 説明と要約を先に配置する。
    path.write_text(json.dumps(formatted,ensure_ascii=False,indent=2)+'\n')  # 日本語をUnicodeエスケープせず保存する。
    make_overview(path.parent)  # 関連集計が揃っていれば共通入口も更新する。
