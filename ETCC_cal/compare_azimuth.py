#!/usr/bin/env python3
"""stage04のシンチ2ヒット方位角について、規格化した水平／円の比を描く。"""
import argparse  # 比較専用JSONを受け取る。
import csv  # 方位角ビンごとの比と誤差を保存する。
import json  # 保存済みヒストグラムと条件を読む。
from pathlib import Path  # 入出力を扱う。
import numpy as np  # 規格化と統計共分散を計算する。
import matplotlib  # 保存用描画を使う。
matplotlib.use('Agg')  # GUIを開かない。
import matplotlib.pyplot as plt  # 分布と比を並べる。
from settings import repo_path  # リポジトリ基準のパスを使う。


def fraction_covariance(counts,sumw2):  # 総数規格化によるビン間相関も保存する。
    counts,sumw2 = np.asarray(counts,float),np.asarray(sumw2,float)  # 元の重み付き計数と重み二乗和を読む。
    total = counts.sum()  # 方位角分布の総重みを得る。
    if total<=0:  # 空の分布を規格化しない。
        raise ValueError('Empty azimuth distribution')  # 入力不正を示す。
    fraction = counts/total  # 各ビンの割合を求める。
    jacobian = (np.eye(len(counts))-fraction[:,None])/total  # n_i / sum(n) の偏微分を使う。
    covariance = (jacobian*sumw2[None,:])@jacobian.T  # 計数誤差を規格化後へ伝播する。
    return fraction,covariance  # 1イベント1ペアでは多項分布の共分散になる。


def compare(horizontal,circular,hvar,cvar):  # 2分布は異なるrunなので独立とする。
    ph,ch = fraction_covariance(horizontal,hvar)  # 水平を総数で規格化する。
    pc,cc = fraction_covariance(circular,cvar)  # 円も同じ方法で規格化する。
    valid = pc>0  # 分母0のビンを割らない。
    ratio,covariance = np.full(len(ph),np.nan),np.full((len(ph),len(ph)),np.nan)  # 未定義はNaNのまま残す。
    indices = np.flatnonzero(valid)  # 比を計算できるビンを得る。
    ratio[valid] = ph[valid]/pc[valid]  # 規格化した水平／円を求める。
    jh,jc = 1/pc[valid],-ph[valid]/pc[valid]**2  # 比の偏微分を求める。
    covariance[np.ix_(indices,indices)] = ch[np.ix_(indices,indices)]*np.outer(jh,jh)+cc[np.ix_(indices,indices)]*np.outer(jc,jc)  # 規格化の相関を含むデルタ法誤差を求める。
    return ph,pc,ratio,covariance,valid  # 統計誤差だけを返す。


def main():  # 既存の最終品質選別ヒストグラムだけから比を作る。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行指定を用意する。
    parser.add_argument('--config',required=True,type=Path)  # 入出力を外部JSONで指定する。
    args = parser.parse_args()  # 引数を読む。
    config = json.loads(args.config.read_text())  # 比較用の条件を読む。
    output = repo_path(config['output_dir'])  # 方位角比を専用フォルダへ置く。
    output.mkdir(parents=True,exist_ok=True)  # 保存先を作る。
    fig,axes = plt.subplots(2,2,figsize=(12,8),constrained_layout=True,sharex=True)  # HG/LGを横、元分布と比を縦に並べる。
    summary,arrays = {'stage':'04 double-hit quality','definition':'(H_phi/H_total)/(C_phi/C_total)','angle':'line joining two scintillator pixel centres, projected perpendicular to saved beam direction; modulo 180 degrees','errors':'statistical delta-method errors with normalization covariance; detector/calibration systematics excluded','gains':{}},{}  # 比と誤差の定義を残す。
    for column,gain in enumerate(('hg','lg')):  # gain別の採用イベントを比較する。
        hp = Path(str(repo_path(config['horizontal_prefix']))+f'_{gain}.json')  # 水平の最終分布を指定する。
        cp = Path(str(repo_path(config['circular_prefix']))+f'_{gain}.json')  # 円の最終分布を指定する。
        h,c = json.loads(hp.read_text()),json.loads(cp.read_text())  # rawを読み直さず保存計数を使う。
        if h['phi_edges_deg']!=c['phi_edges_deg'] or not np.allclose(h['beam_direction'],c['beam_direction']):  # 同じ角度定義とビンを使う。
            raise ValueError('Angular bins or beam direction differ')  # 条件違いなら停止する。
        for item in (h,c):  # 1イベント1ペアの最終条件を確認する。
            if not all(item['options'][key] for key in ('require_exactly_two_raw_hits','beam_vertex_cut_enabled','z0_track_veto_enabled')):  # stage04相当の前提を確認する。
                raise ValueError('Final raw2 beam/Z selection required')  # 別条件を流用しない。
        edges = np.asarray(h['phi_edges_deg'],float)  # 保存ビンの境界を使う。
        centers = (edges[:-1]+edges[1:])/2  # 10度ビンの中心を求める。
        nh,nc = np.asarray(h['weighted_counts'],float),np.asarray(c['weighted_counts'],float)  # 方位角計数を読む。
        ph,pc,ratio,covariance,valid = compare(nh,nc,h['sumw2'],c['sumw2'])  # 総数を揃えた比と誤差を得る。
        errors = np.sqrt(np.maximum(np.diag(covariance),0))  # 各ビンの統計誤差を描く。
        for p,source,label,color in ((ph,h,'Horizontal','tab:blue'),(pc,c,'Circular','tab:orange')):  # 元の規格化分布を上段へ示す。
            _,cov = fraction_covariance(source['weighted_counts'],source['sumw2'])  # 元分布にも規格化後の誤差を使う。
            axes[0,column].errorbar(centers,p,yerr=np.sqrt(np.maximum(np.diag(cov),0)),fmt='o-',ms=4,color=color,label=f'{label}: N={sum(source["weighted_counts"]):,.0f}')  # 統計量と形状を同時に確認する。
        axes[0,column].set(title=f'{gain.upper()}: normalized azimuth distributions',ylabel='Fraction of accepted events / bin')  # ヒット点数でなく採用ペア数と示す。
        axes[0,column].legend(fontsize=9)  # 水平・円の違いを示す。
        axes[1,column].errorbar(centers[valid],ratio[valid],yerr=errors[valid],xerr=np.diff(edges)[valid]/2,fmt='o',color='black',ms=4,capsize=2)  # 下段へ規格化した比を描く。
        axes[1,column].axhline(1,color='gray',ls='--')  # 同じ分布なら1になる基準線を示す。
        axes[1,column].set(title=f'{gain.upper()}: normalized horizontal / circular',xlabel='Azimuth [deg]',ylabel='(H / H total) / (C / C total)',ylim=config['ratio_y_limits'])  # 比の方向を明示する。
        for ax in axes[:,column]:  # 共通の横軸を使う。
            ax.set_xlim(0,180)  # 180度周期の範囲を表示する。
            ax.grid(alpha=.2)  # ビンの差を読みやすくする。
        with (output/f'{gain}_azimuth_ratio.csv').open('w') as stream:  # 元計数と比を保存する。
            writer=csv.writer(stream);writer.writerow(['phi_low_deg','phi_high_deg','phi_center_deg','horizontal_count','circular_count','horizontal_fraction','circular_fraction','normalized_H_over_C','stat_error'])  # 列名に定義を含める。
            writer.writerows(zip(edges[:-1],edges[1:],centers,nh,nc,ph,pc,ratio,errors))  # 全方位角ビンを残す。
        arrays.update({gain+'_ratio':ratio,gain+'_covariance':covariance,gain+'_edges_deg':edges})  # ビン間相関を後の解析に残す。
        summary['gains'][gain]={'horizontal_events':float(nh.sum()),'circular_events':float(nc.sum()),'ratio':ratio.tolist(),'stat_error':errors.tolist(),'horizontal_source':str(hp),'circular_source':str(cp)}  # 読みやすい集計を保存する。
    fig.suptitle('Scintillator azimuth: stage04 double-hit quality\nHorizontal run00..07 / circular run08..12; beam-axis then Z-side veto; HG position calibration (LG provisional)')  # 全体の条件を示す。
    fig.savefig(output/'azimuth_ratio.png',dpi=160)  # HG/LGの分布と比を1枚へ保存する。
    plt.close(fig)  # 描画メモリを解放する。
    np.savez_compressed(output/'azimuth_ratio.npz',**arrays)  # 共分散を数値として残す。
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')  # 集計を保存する。
    print('Saved',output/'azimuth_ratio.png',flush=True)  # 保存先を知らせる。


if __name__=='__main__':  # importだけでは描画しない。
    main()  # コマンドから呼ばれたときだけ処理する。
