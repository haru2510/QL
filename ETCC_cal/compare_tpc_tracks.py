#!/usr/bin/env python3
"""stage04のTPC積算について、同じ表示ビンの規格化した水平／円の比を作る。"""
import argparse  # 比較専用の外部JSONを受け取る。
import json  # 元の選別と規格化を保存する。
from pathlib import Path  # 保存配列と比較図を扱う。
import numpy as np  # ヒット数の合算と比を計算する。
import matplotlib  # 保存専用の描画を使う。
matplotlib.use('Agg')  # GUIを開かない。
import matplotlib.pyplot as plt  # 3投影を並べる。
from matplotlib.colors import Normalize, LogNorm, ListedColormap  # 比と元分布の色を分ける。
from settings import repo_path  # 外部パスをリポジトリ基準で解決する。
from plot_hit_pixel_maps import beam_reference  # run06・07の共通ビーム軸を表示する。


def binned_planes(raw, xy):  # 既存のTPC積算図と同じ表示ビンへ合算する。
    return {'XZ':raw[0].reshape(96,8,256,4).sum(axis=(1,3)).T,'YZ':raw[1].reshape(96,8,256,4).sum(axis=(1,3)).T,'XY':xy.T}  # XZ/YZは8strip×4clock、XYは8strip×8stripを使う。


def normalized_ratio(horizontal,circular,minimum=1):  # 各投影の総数で規格化してから割る。
    htotal,ctotal = int(horizontal.sum()),int(circular.sum())  # XZ/YZ/XYごとに独立した分母を使う。
    if htotal<=0 or ctotal<=0:  # 空の図を規格化しない。
        raise ValueError('Empty TPC projection')  # 元データを確認できるように止める。
    valid = circular>=minimum  # 分母0または指定未満は比を計算しない。
    ratio = np.full(horizontal.shape,np.nan)  # 擬似カウントで値を作らない。
    np.divide(horizontal,circular,out=ratio,where=valid)  # 同一表示ビンの水平／円を計算する。
    ratio *= ctotal/htotal  # 各マップの単位面積規格化へ変換する。
    return ratio,valid,htotal,ctotal  # 元の総数も保存する。


def decorate(ax,plane,edge,zedges,beam,alpha):  # TPC幾何とビーム軸を描く。
    line,marker,_,_ = beam  # 保存直線と底面シンチ平面での通過点を得る。
    if plane=='XY':  # XYはTPCの横断面を表示する。
        ax.set(xlabel='X [mm]',ylabel='Y [mm]',xlim=(edge[0],edge[-1]),ylim=(edge[-1],edge[0]))  # Y増加方向を下にする。
        ax.scatter(marker[0],marker[1],marker='*',c='magenta',edgecolors='black',s=140,alpha=alpha,zorder=5)  # 共通のビーム印を重ねる。
        ax.set_aspect('equal',adjustable='box')  # XYの距離縮尺を揃える。
    else:  # XZ/YZは全クロックの較正範囲を保つ。
        axis = 0 if plane=='XZ' else 1  # どちらの横座標か指定する。
        ax.set(xlabel=('X' if axis==0 else 'Y')+' [mm]',ylabel='Z [mm]',xlim=(edge[0],edge[-1]),ylim=(zedges[-1],zedges[0]))  # Z増加方向を下へ向ける。
        ax.plot(line[:,axis],line[:,2],color='magenta',ls='--',lw=1,alpha=.65)  # 軸位置を共通表示する。


def main():  # 水平・円のstage04完成配列だけを使い、SSDを読み直さず比較する。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行指定を用意する。
    parser.add_argument('--config',required=True,type=Path)  # 入力・出力・色範囲をJSONから読む。
    args = parser.parse_args()  # 引数を読む。
    config = json.loads(args.config.read_text())  # 比較専用設定を読む。
    output = repo_path(config['output_dir'])  # シンチ比較とは別フォルダへ置く。
    output.mkdir(parents=True,exist_ok=True)  # 保存先を作る。
    hflow = json.loads(repo_path(config['horizontal_cutflow']).read_text())  # 水平の採用イベント数を読む。
    cflow = json.loads(repo_path(config['circular_cutflow']).read_text())  # 円の採用イベント数を読む。
    if hflow['branches']!=cflow['branches']:  # 4条件の順序を取り違えない。
        raise ValueError('Branch order mismatch')  # 不一致なら停止する。
    with np.load(repo_path(config['horizontal_histogram'])) as saved:  # 元の完成済み積算を読む。
        hraw,hxy,hconfig = saved['raw'],saved['xy'],json.loads(str(saved['config']))  # 水平の全条件を読む。
    with np.load(repo_path(config['circular_histogram'])) as saved:  # 円も同じ形式で読む。
        craw,cxy,cconfig = saved['raw'],saved['xy'],json.loads(str(saved['config']))  # 円の全条件を読む。
    for key in ('velocity','clock_offset','z_gem','pitch_strip','z_veto'):  # 同じ座標・Z条件であることを確認する。
        if hconfig[key]!=cconfig[key]:  # 違う較正を混ぜない。
            raise ValueError(f'Calibration mismatch: {key}')  # 不一致項目を示す。
    assert hraw.shape==craw.shape==(4,2,768,1024) and hxy.shape==cxy.shape==(4,96,96)  # 元配列の定義を確認する。
    edge = np.linspace(-153.6,153.6,97)  # 横方向3.2 mmのビン境界を作る。
    clock_edges = np.arange(-.5,1024,4)  # 4clock表示ビンを使う。
    zedges = hconfig['velocity']*(clock_edges-hconfig['clock_offset'])+hconfig['z_gem']  # 同じ較正でmmへ変換する。
    beam = beam_reference(hconfig)  # 保存済みrun06・07の共通軸を使う。
    norm_ratio = Normalize(*config['ratio_color_limits'])  # 0〜2の中心1をjetの緑へ置く。
    summary = {'definition':'(H_bin/H_projection_total)/(C_bin/C_projection_total)','stage':'04 gain-specific scintillator quality event selection','normalization':'each gain, multiplicity and projection independently; shape comparison, not exposure/beam-flux normalization','units':{'XZ':'raw strip-clock hits','YZ':'raw strip-clock hits','XY':'same-event/same-uPIC/same-clock X-Y combinations'},'binning':{'transverse_mm':3.2,'clock':4},'colormap':config['ratio_colormap'],'ratio_color_limits':config['ratio_color_limits'],'minimum_circular_counts':config['minimum_circular_counts'],'branches':{}}  # 比の物理的な単位を明示する。
    for index,(branch,gain,multiplicity) in enumerate(hflow['branches']):  # raw1/2×HG/LGを別に描く。
        directory = output/f'04_{branch}_quality'  # シンチと同じstage名を使う。
        directory.mkdir(parents=True,exist_ok=True)  # 保存先を作る。
        hplanes,cplanes = binned_planes(hraw[index],hxy[index]),binned_planes(craw[index],cxy[index])  # 同じ表示ビンで比較する。
        events = [sum(run['selected_events'][index] for run in flow['runs']) for flow in (hflow,cflow)]  # 採用イベント数を明記する。
        record = {'horizontal_events':events[0],'circular_events':events[1],'projections':{}}  # 投影の総ヒット数とイベント数を区別する。
        fig,axes = plt.subplots(1,3,figsize=(18,6),constrained_layout=True)  # XZ・YZ・XYの比を並べる。
        comparison,grid = plt.subplots(3,3,figsize=(17,15),constrained_layout=True)  # 元の水平・円分布も並べて確認できるようにする。
        arrays = {}  # 元件数と比をNPZへ保存する。
        for column,plane in enumerate(('XZ','YZ','XY')):  # 3投影を処理する。
            h,c = hplanes[plane],cplanes[plane]  # 元の表示ビン件数を得る。
            ratio,valid,htotal,ctotal = normalized_ratio(h,c,config['minimum_circular_counts'])  # 同じビンどうしを割る。
            positive = np.concatenate((h[h>0]/htotal,c[c>0]/ctotal))  # 元分布は確率として比較する。
            count_norm = LogNorm(vmin=float(positive.min()),vmax=max(float(positive.max()),float(positive.min())*1.01))  # 水平と円で元分布の色範囲を共通にする。
            yedges = edge if plane=='XY' else zedges  # XYだけ縦軸がY座標になる。
            for row,(counts,total,label) in enumerate(((h,htotal,'Horizontal'),(c,ctotal,'Circular'))):  # 元の規格化済み積算も保存する。
                image=grid[row,column].pcolormesh(edge,yedges,np.ma.masked_where(counts==0,counts/total),cmap='turbo',norm=count_norm)  # 空ビンは白くする。
                comparison.colorbar(image,ax=grid[row,column],label='fraction of projection total / bin')  # 各投影の分母を示す。
                decorate(grid[row,column],plane,edge,zedges,beam,config['beam_marker_alpha'])  # 共通の幾何を使う。
                grid[row,column].set_title(f'{label}: {plane}')  # 元データと投影面を示す。
            for ax,canvas in ((axes[column],fig),(grid[2,column],comparison)):  # 比のみと一覧の両方へ描く。
                undefined=(h>0)&~valid  # 水平にヒットがあるのに円側の分母が足りない場所を特定する。
                ax.pcolormesh(edge,yedges,np.ma.masked_where(~undefined,np.ones_like(h)),cmap=ListedColormap(['lightgray']),vmin=0,vmax=1)  # 両方0は白、分母不足は灰色で区別する。
                image=ax.pcolormesh(edge,yedges,np.ma.masked_where(~valid,ratio),cmap=config['ratio_colormap'],norm=norm_ratio)  # 擬似カウントなしの比を描く。
                canvas.colorbar(image,ax=ax,label='Normalized H/C (green ~ 1)',extend='both')  # 採用配色の1の意味を示す。
                decorate(ax,plane,edge,zedges,beam,config['beam_marker_alpha'])  # 元分布と同じ範囲にする。
                ax.set_title(f'{plane}: normalized horizontal / circular')  # 比の向きを明示する。
            record['projections'][plane]={'horizontal_total':htotal,'circular_total':ctotal,'displayed_bins':int(valid.sum()),'both_zero_bins':int(((h==0)&(c==0)).sum()),'undefined_populated_bins':int(undefined.sum())}  # ゼロ分母の扱いを記録する。
            arrays.update({plane+'_horizontal_count':h,plane+'_circular_count':c,plane+'_normalized_ratio':ratio})  # 表示前の元数値を残す。
        title=f'Stage04 TPC: {gain.upper()}, exactly {multiplicity} raw HA hit(s) THEN gain quality\nHorizontal run00..07 / circular run08..12; H={events[0]:,}, C={events[1]:,} events; each projection normalized separately'  # イベント選別と規格化を示す。
        fig.suptitle(title)  # 比図の見出しを付ける。
        comparison.suptitle(title)  # 元図と比の一覧にも付ける。
        fig.savefig(directory/f'{gain}_tpc_ratio_z.png',dpi=150)  # ミリメートル表示の比を保存する。
        comparison.savefig(directory/f'{gain}_tpc_comparison_z.png',dpi=150)  # 元の2条件と比を保存する。
        plt.close(fig);plt.close(comparison)  # 描画メモリを解放する。
        np.savez_compressed(directory/f'{gain}_tpc_ratio.npz',**arrays,x_edges_mm=edge,z_edges_mm=zedges,clock_edges=clock_edges)  # 後の再表示と数値解析に使えるようにする。
        summary['branches'][f'raw{multiplicity}_{gain}']=record  # 4条件の集計を残す。
        print('COMPLETE TPC RATIO',branch,gain,flush=True)  # 条件ごとの完了を示す。
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')  # 定義と元の総数を保存する。
    print('ALL STAGE04 TPC RATIOS COMPLETE',output,flush=True)  # 全処理完了を記録する。


if __name__=='__main__':  # importだけでは実行しない。
    main()  # コマンドから呼ばれたときだけ実行する。
