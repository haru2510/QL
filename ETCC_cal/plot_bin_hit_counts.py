#!/usr/bin/env python3
"""ビーム円柱内の投影ビンごとの元ヒット数分布を水平・円で比較する。"""
import argparse  # 比較用の外部JSONを受け取る。
import csv  # 比マップと同じ元件数を読む。
import json  # 分布の集計と条件を保存する。
from pathlib import Path  # 保存先を指定する。
import numpy as np  # ヒストグラムと分位点を計算する。
import matplotlib  # 保存専用で描く。
matplotlib.use('Agg')  # GUIを開かない。
import matplotlib.pyplot as plt  # HG/LGの分布を並べる。
from settings import repo_path  # リポジトリ基準のパスを使う。


def main():  # 元の計数を規格化せずに比較する。
    parser = argparse.ArgumentParser(description=__doc__)  # 設定の指定を用意する。
    parser.add_argument('--config',required=True,type=Path)  # 入力元と出力先をJSONから読む。
    args = parser.parse_args()  # 引数を読む。
    config = json.loads(args.config.read_text())  # 既存の比マップ設定を読む。
    root = repo_path(config['output_dir'])  # 円柱内CSVの保存先を読む。
    output = repo_path(config.get('hit_count_distribution_dir',str(root/'hit_count_distribution')))  # 分布図を独立したフォルダへ置く。
    output.mkdir(parents=True,exist_ok=True)  # 保存先を作る。
    plane = config.get('hit_count_projection','xy').lower()  # 今回はXYの投影ビンを使う。
    fig,axes = plt.subplots(2,2,figsize=(13,9),constrained_layout=True)  # HG/LGと低計数部・裾を分ける。
    summary = {'region':'3D distance to saved beam line <= zoom_radius_mm; hit-level selection', 'zoom_radius_mm':config.get('zoom_radius_mm',100),'projection':plane.upper(),'bin_population':'union of projected bins with >=1 hit in either sample; bins zero in BOTH samples are absent','stage':'04_two_hits_quality','normalization':'raw hit uses per projected bin, NOT unit-area fractions','gains':{}}  # 分母と0件の扱いを明記する。
    for row,gain in enumerate(('hg','lg')):  # gainによる採用差を分ける。
        source = root/f'{gain}_{plane}_zoom_ratio.csv'  # 比マップと同じ円柱選別後のビンを読む。
        with source.open() as stream:  # コメントなしCSVを読む。
            records = list(csv.DictReader(stream))  # 各ビンの水平・円の件数を取り出す。
        record = {'projected_bins':len(records),'datasets':{}}  # 両方正のビンに限定しない。
        for dataset,label,color in (('horizontal','Horizontal run00..07','tab:blue'),('circular','Circular run08..12','tab:orange')):  # 同じビン集合で2データを比較する。
            counts = np.array([int(item[dataset+'_count']) for item in records],np.int64)  # 0件を含む元の計数を読む。
            median = float(np.median(counts))  # 各ビンの計数の中央値を得る。
            statistics = {'sum_hit_uses':int(counts.sum()),'zero_bins':int((counts==0).sum()),'mean':float(counts.mean()),'median':median,'q25':float(np.quantile(counts,.25)),'q75':float(np.quantile(counts,.75)),'p90':float(np.quantile(counts,.90)),'p95':float(np.quantile(counts,.95)),'p99':float(np.quantile(counts,.99)),'max':int(counts.max()),'bins_ge_50':int((counts>=50).sum())}  # 裾を含む要約を残す。
            record['datasets'][dataset]=statistics  # データごとの集計を保存する。
            axes[row,0].hist(counts,bins=np.arange(-.5,51.5,1),histtype='step',color=color,lw=1.7,label=f'{label}; median={median:g}')  # 多数を占める0〜50件を1件刻みで描く。
            axes[row,0].axvline(median,color=color,ls=':',lw=1)  # 中央値を示す。
            frequency = np.bincount(counts)  # 全計数範囲の正確な頻度を得る。
            thresholds = np.flatnonzero(frequency)  # 実際に観測した計数値を使う。
            thresholds = thresholds[thresholds>0]  # 対数横軸には正値だけを使う。
            survival = np.cumsum(frequency[::-1])[::-1]  # N件以上ある投影ビン数を求める。
            axes[row,1].step(thresholds,survival[thresholds],where='post',color=color,lw=1.6,label=label)  # 数千件まで伸びる裾を累積分布で示す。
            with (output/f'{gain}_{dataset}_frequency.csv').open('w') as stream:  # 0件を含む全頻度を保存する。
                writer=csv.writer(stream);writer.writerow(['hits_per_projected_bin','number_of_projected_bins','number_of_projected_bins_ge_threshold'])  # 単位を列名へ入れる。
                writer.writerows((i,int(n),int(survival[i])) for i,n in enumerate(frequency))  # 空の計数階級も省略しない。
        axes[row,0].set(xlim=(-.5,50.5),xlabel='Hit uses per projected bin',ylabel='Number of projected bins',title=f'{gain.upper()}: low-count distribution (0..50 hits)')  # 横軸はイベント数でなくヒット使用数と示す。
        axes[row,1].set(xscale='log',yscale='log',xlabel='N: hit uses per projected bin',ylabel='Number of projected bins with >= N hits',title=f'{gain.upper()}: full-range cumulative tail')  # 全範囲は対数軸で示す。
        for ax in axes[row]:  # 図の補助情報を揃える。
            ax.legend(fontsize=8)  # 2データの色を明示する。
            ax.grid(alpha=.2)  # 件数を読みやすくする。
        summary['gains'][gain]=record  # gainの要約を保存する。
    fig.suptitle(f"Final double-hit quality; {plane.upper()} hit-count distribution inside beam cylinder R={config.get('zoom_radius_mm',100):g} mm\nSame union of bins for H and C; includes one-sided zeros, excludes bins zero in both; raw counts (exposures differ)")  # 集計対象を図中に明示する。
    fig.savefig(output/'hit_counts.png',dpi=160)  # 分布の全体像を保存する。
    plt.close(fig)  # 描画メモリを解放する。
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')  # 正確な数値を保存する。
    print(json.dumps(summary['gains'],indent=2),flush=True)  # 中央値・裾・最大値を確認できるようにする。


if __name__=='__main__':  # importだけでは図を作らない。
    main()  # コマンド実行時だけ処理する。
