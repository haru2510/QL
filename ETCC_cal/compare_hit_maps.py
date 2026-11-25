#!/usr/bin/env python3
"""最終raw2品質選別後の水平／円偏光シンチ位置分布を、総ヒット数で規格化して割る。"""
import argparse  # 外部JSONを指定する。
import csv  # 各ビンの元件数と比を保存する。
import json  # 規格化と除外ビン数を記録する。
from pathlib import Path  # 入出力を解決する。
import numpy as np  # 共通座標のヒット数を集計する。
import matplotlib  # 保存専用の描画を使う。
matplotlib.use('Agg')  # GUIを開かない。
import matplotlib.pyplot as plt  # 2条件と比を並べる。
from matplotlib.colors import TwoSlopeNorm, LogNorm  # 比は線形、元分布は対数で描く。
from settings import repo_path, load_config  # リポジトリ基準のパスと元解析設定を使う。
from plot_hit_pixel_maps import beam_reference  # 同じビーム星印を描く。


def ratio_bins(horizontal, circular, axes, minimum=1, bin_mm=0, normalization_totals=None):  # 投影してから同じ位置の件数を割る。
    combined = np.vstack((horizontal[:,axes],circular[:,axes]))  # 同じ投影座標へ集める。
    if bin_mm > 0:  # 見やすい図では共通の格子へ件数を合算する。
        combined = (np.floor(combined/bin_mm)+0.5)*bin_mm  # 両データに同じ幅と原点の表示ビンを使う。
    coordinates, inverse = np.unique(combined,axis=0,return_inverse=True)  # 両データの共通座標を作る。
    nh = np.bincount(inverse[:len(horizontal)],minlength=len(coordinates))  # 水平側の件数を数える。
    nc = np.bincount(inverse[len(horizontal):],minlength=len(coordinates))  # 円側の件数を数える。
    valid = nc >= minimum  # 円側0や指定した低統計ビンを割らない。
    raw = np.full(len(coordinates),np.nan)  # 分母0を未定義として保存する。
    np.divide(nh,nc,out=raw,where=valid)  # 擬似カウントを加えず元の比を計算する。
    htotal, ctotal = normalization_totals if normalization_totals is not None else (len(horizontal),len(circular))  # 拡大図でも全体と同じ規格化分母を保つ。
    normalized = raw*ctotal/htotal  # 各マップの総ヒット数で規格化した形状比へ直す。
    assert nh.sum()==len(horizontal) and nc.sum()==len(circular)  # 投影で件数が増減していないことを確認する。
    return coordinates,nh,nc,valid,raw,normalized  # 元件数を残して比の評価を可能にする。


def decorate(ax, axes, beam, config):  # 元マップと同じ幾何と半透明の星を使う。
    i,j = axes  # 表示する座標を得る。
    line,marker,fit_range,_ = beam  # 保存軸の投影を得る。
    if j==2:  # XZ/YZにはビーム直線を重ねる。
        ax.plot(line[:,i],line[:,j],color='magenta',ls='--',lw=1)  # fit範囲外も外挿として示す。
    ax.scatter(marker[i],marker[j],marker='*',s=150,c='magenta',edgecolors='black',linewidths=.6,alpha=config['beam_marker_alpha'],zorder=5)  # 下のピクセルが見える星印を置く。
    ax.set(xlabel='XYZ'[i]+' [mm]',ylabel='XYZ'[j]+' [mm]',xlim=(-215,215),ylim=(-215,215) if j==1 else (-5,375))  # 元図と同じ範囲にする。
    ax.invert_yaxis()  # YとZの増加方向を下へ向ける。
    ax.set_aspect('equal',adjustable='box')  # 距離の縮尺を揃える。


def main():  # HG/LGごとに水平分布・円分布・比を保存する。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行条件を用意する。
    parser.add_argument('--config',required=True,type=Path)  # 2つの入力と比較出力を外部JSONで指定する。
    args = parser.parse_args()  # 指定を読む。
    config = json.loads(args.config.read_text())  # 比較専用の条件を読む。
    equal_color = config.get('ratio_equal_color','white')  # 比が1の色を設定に合わせて説明する。
    output = repo_path(config['output_dir'])  # 比較結果を別フォルダに保存する。
    output.mkdir(parents=True,exist_ok=True)  # 初回に保存先を作る。
    source_config = load_config(repo_path(config['horizontal_config']))  # ビーム軸と選別条件を読む。
    circular_config = load_config(repo_path(config['circular_config']))  # 同じカットになっているか確認する。
    for key in ('velocity','clock_offset','z_gem','z_veto'):  # 座標とZ条件の一致を検証する。
        if source_config[key]!=circular_config[key]:  # 同条件と見なせない比較を止める。
            raise ValueError(f'Analysis condition mismatch: {key}')  # 不一致項目を示す。
    for key in ('beam_max_distance_mm','min_pair_distance_mm','common_pairs_only'):  # ペア選別も一致させる。
        if source_config['polarization'][key]!=circular_config['polarization'][key]:  # 形状比較の前提を確認する。
            raise ValueError(f'Pair condition mismatch: {key}')  # 不一致なら停止する。
    beam = beam_reference(source_config)  # 共通のrun06・07ビーム軸を使う。
    summary = {'definition':'(H_bin / H_total) / (C_bin / C_total); projection THEN division', 'stage':'04_two_hits_quality', 'normalization':'per-gain total accepted hit uses, NOT exposure or beam flux', 'uncertainty':'No significance or uncertainty assigned; two hits per event are correlated; low-count ratios fluctuate', 'projection_bin_mm':config.get('projection_bin_mm',0), 'ratio_colormap':config.get('ratio_colormap','RdBu_r'), 'ratio_color_limits':config['ratio_color_limits'], 'gains':{}}  # 形状比と絶対率を区別する。
    for gain in ('hg','lg'):  # gain別の品質選別結果を比較する。
        hpath = repo_path(config['horizontal_dir'])/f'{gain}.points.npz'  # 水平の採用ピクセル位置を読む。
        cpath = repo_path(config['circular_dir'])/f'{gain}.points.npz'  # 円の採用ピクセル位置を読む。
        with np.load(hpath) as saved:  # 読んだらファイルを閉じる。
            horizontal = saved['positions']  # 水平の全採用ヒットを得る。
        with np.load(cpath) as saved:  # 円側も同じ形式で読む。
            circular = saved['positions']  # 円の全採用ヒットを得る。
        if not len(horizontal) or not len(circular):  # 空データを規格化しない。
            raise ValueError(f'Empty map: {gain}')  # 不正入力を示す。
        record = {'horizontal_hit_uses':len(horizontal),'circular_hit_uses':len(circular),'horizontal_events':len(horizontal)//2,'circular_events':len(circular)//2,'raw_to_normalized_factor':len(circular)/len(horizontal),'projections':{}}  # 規格化係数を保存する。
        products = [ratio_bins(horizontal,circular,axes,config['minimum_circular_counts'],config.get('projection_bin_mm',0)) for axes in ((0,1),(0,2),(1,2))]  # XY・XZ・YZの比を求める。
        radius = config.get('zoom_radius_mm',100)  # 下段で拡大するビーム円柱の半径を読む。
        direction = beam[0][-1]-beam[0][0]  # 保存直線から方向を得る。
        direction = direction/np.linalg.norm(direction)  # 直線との距離を計算するため単位化する。
        hmask = np.linalg.norm(np.cross(horizontal-beam[0][0],direction),axis=1)<=radius  # 投影前の真の3D座標で円柱内ヒットを選ぶ。
        cmask = np.linalg.norm(np.cross(circular-beam[0][0],direction),axis=1)<=radius  # 円側にも同じ円柱を使う。
        zoom_products = [ratio_bins(horizontal[hmask],circular[cmask],axes,config['minimum_circular_counts'],config.get('projection_bin_mm',0),(len(horizontal),len(circular))) for axes in ((0,1),(0,2),(1,2))]  # 拡大内だけで再規格化しない。
        record['zoom'] = {'radius_mm':radius,'horizontal_hit_uses':int(hmask.sum()),'circular_hit_uses':int(cmask.sum()),'selection':'hit-level 3D distance to saved beam line; full-map normalization retained'}  # 切り出しがイベント再選別ではないことを残す。
        positive = np.concatenate([counts[counts>0]/total for product in products for counts,total in ((product[1],len(horizontal)),(product[2],len(circular)))])  # 元分布は共通カラー範囲にする。
        norm_counts = LogNorm(vmin=float(positive.min()),vmax=float(positive.max()))  # 確率密度の対数色を使う。
        norm_ratio = TwoSlopeNorm(vmin=config['ratio_color_limits'][0],vcenter=1.0,vmax=config['ratio_color_limits'][1])  # 比は0を表示できる線形色にする。
        fig,axes_grid = plt.subplots(3,3,figsize=(16,14),constrained_layout=True)  # 水平・円・比を縦に並べる。
        ratio_fig,ratio_axes = plt.subplots(2,3,figsize=(17,10),constrained_layout=True)  # 比だけの図も用意する。
        for column,(axes,product) in enumerate(zip(((0,1),(0,2),(1,2)),products)):  # 3つの投影を描く。
            coordinates,nh,nc,valid,raw,ratio = product  # 同一ビンの元件数と比を得る。
            plane = ''.join('XYZ'[i] for i in axes)  # 平面名を作る。
            with (output/f'{gain}_{plane.lower()}_ratio.csv').open('w') as stream:  # 低統計や分母0も監査できるCSVを保存する。
                writer = csv.writer(stream)  # 表を用意する。
                writer.writerow([plane[0]+'_mm',plane[1]+'_mm','horizontal_count','circular_count','horizontal_fraction','circular_fraction','raw_H_over_C','normalized_H_over_C','displayed'])  # 比の定義を列に残す。
                for xy,h,c,r,q,show in zip(coordinates,nh,nc,raw,ratio,valid):  # 全位置を保存する。
                    writer.writerow([*xy,int(h),int(c),h/len(horizontal),c/len(circular),r,q,int(show)])  # 元件数と規格化比を一緒に出す。
            for row,(counts,total,label) in enumerate(((nh,len(horizontal),'Horizontal'),(nc,len(circular),'Circular'))):  # 元の規格化分布も併記する。
                keep = counts>0  # 空ビンは描かない。
                image = axes_grid[row,column].scatter(coordinates[keep,0],coordinates[keep,1],c=counts[keep]/total,s=12,marker='s',cmap=config['colormap'],norm=norm_counts)  # 総数で割った各ビンの割合を描く。
                fig.colorbar(image,ax=axes_grid[row,column],label='accepted hit fraction / projected pixel')  # 元分布も規格化済みと示す。
                axes_grid[row,column].set_title(f'{label}: {plane}')  # データと平面を示す。
                decorate(axes_grid[row,column],axes,beam,config)  # 共通の範囲とビーム印を使う。
            for ax,canvas in ((axes_grid[2,column],fig),(ratio_axes[0,column],ratio_fig)):  # 比だけと比較一覧の両方を描く。
                ax.scatter(coordinates[~valid,0],coordinates[~valid,1],c='lightgray',s=config.get('ratio_marker_size',12),marker='s')  # 円が0または指定未満の場所は灰色にする。
                image = ax.scatter(coordinates[valid,0],coordinates[valid,1],c=ratio[valid],s=config.get('ratio_marker_size',12),marker='s',cmap=config.get('ratio_colormap','RdBu_r'),norm=norm_ratio)  # 規格化した水平／円を描く。
                canvas.colorbar(image,ax=ax,label=f'Normalized H/C ({equal_color} = 1; red = H excess)',extend='both')  # 1が同じ形状であると分かるラベルにする。
                ax.set_title(f'Normalized H/C: {plane}')  # 比の向きを明記する。
                decorate(ax,axes,beam,config)  # 半透明のビーム印を付ける。
            zoom_coordinates,zh,zc,zvalid,zraw,zratio = zoom_products[column]  # 円柱内だけを集計した比を得る。
            zoom_ax = ratio_axes[1,column]  # 全体図の真下へ拡大図を置く。
            zoom_ax.scatter(zoom_coordinates[~zvalid,0],zoom_coordinates[~zvalid,1],c='lightgray',s=config.get('zoom_marker_size',48),marker='s')  # 分母不足は拡大図でも灰色にする。
            zoom_image = zoom_ax.scatter(zoom_coordinates[zvalid,0],zoom_coordinates[zvalid,1],c=zratio[zvalid],s=config.get('zoom_marker_size',48),marker='s',cmap=config.get('ratio_colormap','RdBu_r'),norm=norm_ratio)  # 元ピクセル位置で同じ色範囲を使う。
            ratio_fig.colorbar(zoom_image,ax=zoom_ax,label='Normalized H/C (same scale)',extend='both')  # 上下で色の意味を変えない。
            decorate(zoom_ax,axes,beam,config)  # 同じビーム印と距離スケールを使う。
            zoom_ax.set_xlim(float(beam[0][:,axes[0]].min()-radius-5),float(beam[0][:,axes[0]].max()+radius+5))  # ビーム円柱の投影に横幅を合わせる。
            if axes[1]==1:  # XYは縦方向も円柱付近へ拡大する。
                zoom_ax.set_ylim(float(beam[0][:,1].max()+radius+5),float(beam[0][:,1].min()-radius-5))  # Y正方向を下向きに保つ。
            zoom_ax.set_title(f'{plane}: within beam cylinder R={radius:g} mm')  # 半径100 mmの切り出しと明示する。
            with (output/f'{gain}_{plane.lower()}_zoom_ratio.csv').open('w') as stream:  # 拡大図の元件数も保存する。
                writer=csv.writer(stream);writer.writerow([plane[0]+'_mm',plane[1]+'_mm','horizontal_count','circular_count','normalized_H_over_C','displayed'])  # 全体規格化の比を保存する。
                writer.writerows([*xy,int(h),int(c),r,int(v)] for xy,h,c,r,v in zip(zoom_coordinates,zh,zc,zratio,zvalid))  # 円柱内の全ビンを残す。
            record['projections'][plane] = {'displayed_bins':int(valid.sum()),'masked_bins':int((~valid).sum()),'horizontal_hits_in_masked_bins':int(nh[~valid].sum()),'ratio_above_color_max_bins':int(np.count_nonzero(ratio[valid]>config['ratio_color_limits'][1]))}  # 分母不足と色上限超過を記録する。
        title = f'{gain.upper()}: final double-hit quality; horizontal run00..07 / circular run08..12\nUnit-area normalization; spatial bin={config.get("projection_bin_mm",0):g} mm (0=native); {equal_color}=1; H={len(horizontal)//2:,}, C={len(circular)//2:,} events; gray: C count < {config["minimum_circular_counts"]}; HG calibration (LG provisional)'  # 条件と分母を図から読めるようにする。
        fig.suptitle(title)  # 比較一覧の見出しを付ける。
        ratio_fig.suptitle(title + f'\nTop: full map; bottom: 3D beam-cylinder R={radius:g} mm; same full-map normalization')  # 比だけの図にも同じ説明を付ける。
        fig.savefig(output/f'{gain}_comparison.png',dpi=150)  # 水平・円・比の一覧を保存する。
        ratio_fig.savefig(output/f'{gain}_ratio.png',dpi=150)  # XY・XZ・YZの比を保存する。
        plt.close(fig);plt.close(ratio_fig)  # 描画メモリを解放する。
        coords,nh,nc,valid,raw,ratio = ratio_bins(horizontal,circular,(0,1,2),config['minimum_circular_counts'])  # 投影前の実ピクセル中心ごとの比も保存する。
        with (output/f'{gain}_3d_pixel_ratio.csv').open('w') as stream:  # 真の3D位置での比を監査できるようにする。
            writer=csv.writer(stream);writer.writerow(['x_mm','y_mm','z_mm','horizontal_count','circular_count','raw_H_over_C','normalized_H_over_C','displayed'])  # 位置と件数の列を示す。
            writer.writerows([*xyz,int(h),int(c),r,q,int(show)] for xyz,h,c,r,q,show in zip(coords,nh,nc,raw,ratio,valid))  # 全位置を残す。
        summary['gains'][gain]=record  # gainごとの集計を保存する。
        print(gain,record,flush=True)  # 件数とマスク数を表示する。
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')  # 規格化・比の定義・件数を保存する。


if __name__=='__main__':  # importだけでは図を作らない。
    main()  # コマンドから呼ばれた場合だけ処理する。
