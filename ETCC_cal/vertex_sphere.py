#!/usr/bin/env python3
"""散乱点近傍のraw占有を保存し、球内積算とX/Y両分割の最下流代表点を求める。"""
import numpy as np  # 局所strip占有と投影を扱う。
from numba import njit  # 多重ヒット候補の球選別を高速化する。


@njit  # 3D候補を展開する前のコンパクトなraw占有を取り出す。
def take_local_cells(occupancy,occupied_bins,bin_counts,vertex,radius,velocity,offset,z_gem):  # 後の再計算用には球を囲む立方体を残す。
    cells=np.empty((int(bin_counts.sum()),4),np.int16)  # axis,clock,3.2mm-bin,8bit占有を保存する。
    used=0  # 保存したセル数を数える。
    for clock in range(1024):  # 全クロックの占有を最後に必ず初期化する。
        dz=velocity*(clock-offset)+z_gem-vertex[2]  # 散乱点からのZ差を得る。
        for axis in range(2):  # X/Yのraw stripを別々に残す。
            for i in range(bin_counts[axis,clock]):  # 非空ビンだけ走査する。
                coarse=occupied_bins[axis,clock,i]  # 全体3.2mmビン番号を得る。
                if abs(dz)<=radius and abs(-152.+3.2*coarse-vertex[axis])<=radius+2.:  # 球を覆う範囲だけ保存する。
                    cells[used,0]=axis  # 軸を保存する。
                    cells[used,1]=clock  # 元クロックを保存する。
                    cells[used,2]=coarse  # 元の位置ビンを保存する。
                    cells[used,3]=occupancy[axis,clock,coarse]  # 8本の実strip占有を保持する。
                    used+=1  # 保存数を増やす。
                occupancy[axis,clock,coarse]=0  # 次イベントと混ざらないように消す。
            bin_counts[axis,clock]=0  # 非空ビン数も初期化する。
    return cells[:used]  # 原点合わせや分割の変更時に再利用できる情報を返す。


@njit  # 1イベント・1uPICの局所rawから球内全X-Y候補を処理する。
def process_cells(cells,vertex,radius,velocity,offset,z_gem,mask,absolute_z,absolute_xy,relative,best_clock,sums,npoints):  # 最大Zの同率候補は重心にする。
    low=max(0,int(np.ceil((vertex[2]-radius-z_gem)/velocity+offset)))  # 球が含み得る最小クロックを得る。
    high=min(1023,int(np.floor((vertex[2]+radius-z_gem)/velocity+offset)))  # 最大クロックを得る。
    if high<low:  # 読み出し範囲と交わらない球を除く。
        return 0  # 球内点なしとする。
    occupancy=np.zeros((2,high-low+1,96),np.uint8)  # 球のZ範囲だけを確保する。
    indices=np.empty((2,high-low+1,96),np.int16)  # 非空位置ビンだけを使う。
    numbers=np.zeros((2,high-low+1),np.int16)  # クロックごとの非空ビン数を得る。
    for row in range(len(cells)):  # キャッシュしたrawセルを展開する。
        axis,clock,coarse,bits=cells[row,0],cells[row,1],cells[row,2],cells[row,3]  # 元情報を取り出す。
        if low<=clock<=high:  # 小さい半径へ再選別する場合にも対応する。
            k=clock-low  # 局所配列のクロックに直す。
            if occupancy[axis,k,coarse]==0:  # 同じセルの二重登録を避ける。
                indices[axis,k,numbers[axis,k]]=coarse  # ヒット位置を記録する。
                numbers[axis,k]+=1  # ビン数を増やす。
            occupancy[axis,k,coarse]|=np.uint8(bits)  # 重複stripはORで一回だけ数える。
    accepted=0  # 球内の3D候補数を得る。
    bins=relative.shape[-1]  # 散乱点原点の表示ビン数を得る。
    width=2*radius/bins  # 相対表示ビンの幅を得る。
    for k in range(high-low+1):  # 球のZ範囲だけ走査する。
        clock=k+low  # 元のclockへ戻す。
        z=velocity*(clock-offset)+z_gem  # Z座標を得る。
        dz=z-vertex[2]  # 球中心からのZ差を得る。
        for ix in range(numbers[0,k]):  # 同一uPIC・同クロックのXを得る。
            xb=indices[0,k,ix]  # Xの表示ビンを得る。
            for xbit in range(8):  # 実ストリップ中心で球を判定する。
                if not ((int(occupancy[0,k,xb])>>xbit)&1):  # 非ヒットstripを飛ばす。
                    continue  # 次のstripへ進む。
                x=-153.4+.4*(xb*8+xbit)  # X中心のmm座標を得る。
                dx=x-vertex[0]  # 散乱点からのX差を得る。
                if dx*dx+dz*dz>radius*radius:  # Yによらず外側なら組合せを省く。
                    continue  # 使えるXだけを残す。
                for iy in range(numbers[1,k]):  # 同じclockのYだけを組み合わせる。
                    yb=indices[1,k,iy]  # Yの表示ビンを得る。
                    for ybit in range(8):  # 実Ystripを走査する。
                        if not ((int(occupancy[1,k,yb])>>ybit)&1):  # 非ヒットを除く。
                            continue  # 次へ進む。
                        y=-153.4+.4*(yb*8+ybit)  # Y中心をmmへ直す。
                        dy=y-vertex[1]  # 球中心からのY差を得る。
                        if dx*dx+dy*dy+dz*dz>radius*radius:  # 真の3D球外を除く。
                            continue  # 球内候補だけを使う。
                        accepted+=1  # 同じ3D候補を全投影で数える。
                        rx=min(bins-1,max(0,int((dx+radius)/width)))  # 相対Xビンを得る。
                        ry=min(bins-1,max(0,int((dy+radius)/width)))  # 相対Yビンを得る。
                        rz=min(bins-1,max(0,int((dz+radius)/width)))  # 相対Zビンを得る。
                        for group in range(len(mask)):  # HG/LG・raw1/2の同じsource選別を使う。
                            if mask[group]:  # 採用イベントにだけ足す。
                                absolute_z[group,0,xb,clock//4]+=1  # 検出器座標のXZを積算する。
                                absolute_z[group,1,yb,clock//4]+=1  # YZも同じ候補数を積算する。
                                absolute_xy[group,xb,yb]+=1  # 検出器XYへ積算する。
                                relative[group,0,rx,ry]+=1  # 散乱点を原点にそろえたXYへ積算する。
                                relative[group,1,rx,rz]+=1  # 相対XZへ積算する。
                                relative[group,2,ry,rz]+=1  # 相対YZへ積算する。
                        for split in range(2):  # X分割・Y分割を両方計算する。
                            distance=dx if split==0 else dy  # 散乱点を通る分割面からの符号を見る。
                            if abs(distance)<1e-10:  # 面上の候補はどちらの側にも割り当てない。
                                continue  # 積算には残すが端点には使わない。
                            side=0 if distance<0 else 1  # 負側・正側へ分ける。
                            if clock>best_clock[split,side]:  # 正のドリフト速度なので最大clockが最大Zになる。
                                best_clock[split,side]=clock  # その側の最下流面を更新する。
                                sums[split,side,0]=x  # 新しい面のX和を始める。
                                sums[split,side,1]=y  # Y和を始める。
                                sums[split,side,2]=z  # Z和を始める。
                                npoints[split,side]=1  # 同率点数を初期化する。
                            elif clock==best_clock[split,side]:  # 同じ最大Zが複数なら重心を取る。
                                sums[split,side,0]+=x  # Xを加算する。
                                sums[split,side,1]+=y  # Yを加算する。
                                sums[split,side,2]+=z  # Zを加算する。
                                npoints[split,side]+=1  # 重心の分母を増やす。
    return accepted  # 球内点なしのイベントも区別できるようにする。
