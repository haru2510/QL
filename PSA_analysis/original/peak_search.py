#!/usr/bin/env python3

import argparse
import os

import h5py
import numpy as np

from natsort import natsorted


#
# Default analysis window
#
EMIN=30e3
EMAX=500e3

#
# 2D histogram binning
#
NBINS_X=250
NBINS_Y=300

X_MIN,X_MAX=-0.32,0.32
Y_MIN,Y_MAX=-0.16,0.16

#
# Peak finder settings
#
NPEAKS=64

SMOOTH_XY=(3,3)

THR_REL=0.01

MIN_DIST_BINS=3

ALLOW_EDGE_PEAKS=True


def draw_summary(
    counts_per_ch,
    peaks_per_ch,
    x_edges,
    y_edges,
    fig_path,
    draw
):

    try:

        import matplotlib.pyplot as plt

    except Exception as e:

        print(f"[WARN] matplotlib import failed : {e}")

        return

    fig,axes=plt.subplots(
        3,
        2,
        figsize=(8,12),
        constrained_layout=True
    )

    axes=axes.ravel()

    for ich in range(6):

        ax=axes[ich]

        c=counts_per_ch[ich]

        mesh=ax.pcolormesh(
            x_edges,
            y_edges,
            c.T,
            shading="auto",
            cmap="viridis"
        )

        fig.colorbar(
            mesh,
            ax=ax,
            fraction=0.046,
            pad=0.04
        )

        ax.set_title(f"ch{ich}")

        ax.set_xlabel("x_hg")
        ax.set_ylabel("y_hg")

        peaks=peaks_per_ch.get(ich,[])

        if peaks:

            px=[p[0] for p in peaks]
            py=[p[1] for p in peaks]

            ax.scatter(
                px,
                py,
                s=18,
                c="red",
                marker="x"
            )

            for ip,(xv,yv,_zv) in enumerate(peaks,start=1):

                ax.text(
                    xv,
                    yv,
                    str(ip),
                    fontsize=6,
                    color="white"
                )

    fig.suptitle("Peak map per channel")

    fig.savefig(
        fig_path,
        dpi=150
    )

    print(f"[SAVE] {fig_path}")

    if draw:

        plt.show()

    else:

        plt.close(fig)


def _box_smooth_1d(v,w):

    if w<=0:
        return v

    k=2*w+1

    vp=np.pad(
        v,
        (w,w),
        mode="reflect"
    )

    cs=np.concatenate((
        [0.0],
        np.cumsum(vp,dtype=np.float64)
    ))

    return (cs[k:]-cs[:-k])/float(k)


def box_smooth_2d(a,wx,wy):

    out=a

    if wy>0:

        tmp=np.empty_like(out)

        for ix in range(out.shape[0]):

            tmp[ix,:]=_box_smooth_1d(
                out[ix,:],
                wy
            )

        out=tmp

    if wx>0:

        tmp=np.empty_like(out)

        for iy in range(out.shape[1]):

            tmp[:,iy]=_box_smooth_1d(
                out[:,iy],
                wx
            )

        out=tmp

    return out


def find_peaks_localmax(
    counts,
    x_edges,
    y_edges,
    npeaks,
    thr_rel,
    smooth_xy,
    min_dist_bins,
    allow_edge_peaks
):

    if counts.size==0:
        return []

    wx,wy=int(smooth_xy[0]),int(smooth_xy[1])

    if wx>0 or wy>0:

        smoothed=box_smooth_2d(
            counts,
            wx,
            wy
        )

    else:

        smoothed=counts

    vmax=float(smoothed.max())

    if vmax<=0:
        return []

    thr=float(thr_rel)*vmax

    pad=np.pad(
        smoothed,
        ((1,1),(1,1)),
        mode="constant",
        constant_values=-np.inf
    )

    cen=pad[1:-1,1:-1]

    n1=pad[0:-2,1:-1]
    n2=pad[2:,1:-1]

    n3=pad[1:-1,0:-2]
    n4=pad[1:-1,2:]

    n5=pad[0:-2,0:-2]
    n6=pad[0:-2,2:]

    n7=pad[2:,0:-2]
    n8=pad[2:,2:]

    m=cen>=thr

    m&=cen>=n1
    m&=cen>=n2
    m&=cen>=n3
    m&=cen>=n4
    m&=cen>=n5
    m&=cen>=n6
    m&=cen>=n7
    m&=cen>=n8

    if (
        not allow_edge_peaks
        and cen.shape[0]>2
        and cen.shape[1]>2
    ):

        m[0,: ]=False
        m[-1,:]=False

        m[:,0 ]=False
        m[:,-1]=False

    idx=np.argwhere(m)

    if idx.size==0:
        return []

    cand_z=smoothed[
        idx[:,0],
        idx[:,1]
    ]

    order=np.argsort(cand_z)[::-1]

    x_centers=0.5*(x_edges[:-1]+x_edges[1:])
    y_centers=0.5*(y_edges[:-1]+y_edges[1:])

    min_d =int(max(0,min_dist_bins))
    min_d2=float(min_d*min_d)

    kept_bins=[]

    peaks=[]

    for k in order:

        ix,iy=int(idx[k,0]),int(idx[k,1])

        if min_d>0 and kept_bins:

            ok_far=True

            for jx,jy in kept_bins:

                dx=float(ix-jx)
                dy=float(iy-jy)

                if (dx*dx+dy*dy)<min_d2:

                    ok_far=False

                    break

            if not ok_far:
                continue

        peaks.append((
            float(x_centers[ix]),
            float(y_centers[iy]),
            float(counts[ix,iy])
        ))

        kept_bins.append((ix,iy))

        if len(peaks)>=int(npeaks):
            break

    return peaks


def process_single_directory(in_dir,args):

    fnames=natsorted([

        os.path.join(in_dir,fn)

        for fn in os.listdir(in_dir)

        if fn.endswith(".h5")
    ])

    if len(fnames)==0:

        print(f"[WARN] no h5 files : {in_dir}")

        return

    tag=args.tag.strip() if args.tag else os.path.basename(in_dir)

    out_dir_base=args.out_dir or os.path.join(in_dir,"peak_output")

    out_dir=f"{out_dir_base}_{tag}" if tag else out_dir_base

    counts_per_ch=[

        np.zeros(
            (NBINS_X,NBINS_Y),
            dtype=np.float64
        )

        for _ in range(6)
    ]

    x_edges=np.linspace(
        X_MIN,
        X_MAX,
        NBINS_X+1
    )

    y_edges=np.linspace(
        Y_MIN,
        Y_MAX,
        NBINS_Y+1
    )

    for name in fnames:

        print(name)

        with h5py.File(name,"r") as f:

            ene=f["hg"][:].astype(np.float64)
            x  =f["x_hg"][:].astype(np.float64)
            y  =f["y_hg"][:].astype(np.float64)
            ch =f["ch"][:].astype(np.int32)

        m=(ene>=args.emin)&(ene<=args.emax)

        if not np.any(m):
            continue

        x_m =x[m]
        y_m =y[m]
        ch_m=ch[m]

        for ich in range(6):

            mm=(ch_m==ich)

            if not np.any(mm):
                continue

            hist2d,_,_=np.histogram2d(
                x_m[mm],
                y_m[mm],
                bins =[NBINS_X,NBINS_Y],
                range=[[X_MIN,X_MAX],[Y_MIN,Y_MAX]]
            )

            counts_per_ch[ich]+=hist2d

    print("\n[DIAG] 2D histogram stats per channel")

    for ich,c in enumerate(counts_per_ch):

        print(
            f"  ch{ich}: "
            f"entries={int(c.sum())}, "
            f"max={float(c.max()):.1f}"
        )

    os.makedirs(out_dir,exist_ok=True)

    peaks_per_ch={}

    for ich,counts in enumerate(counts_per_ch):

        peaks=find_peaks_localmax(
            counts           =counts,
            x_edges          =x_edges,
            y_edges          =y_edges,
            npeaks           =args.npeaks,
            thr_rel          =THR_REL,
            smooth_xy        =SMOOTH_XY,
            min_dist_bins    =MIN_DIST_BINS,
            allow_edge_peaks =ALLOW_EDGE_PEAKS
        )

        peaks_per_ch[ich]=peaks

        print(f"\n[CH {ich}] peaks={len(peaks)}")

        for ip,(px,py,pz) in enumerate(peaks,start=1):

            print(
                f"  {ip:02d}: "
                f"x={px:+.6f}, "
                f"y={py:+.6f}, "
                f"z={pz:.1f}"
            )

    csv_name=f"position_{tag}.csv"

    csv_path=os.path.join(
        out_dir,
        csv_name
    )

    with open(csv_path,"w",encoding="utf-8") as w:

        w.write("ch,No,x,y,z\n")

        for ich,peaks in peaks_per_ch.items():

            for ip,(px,py,pz) in enumerate(peaks,start=1):

                w.write(
                    f"{ich},"
                    f"{ip},"
                    f"{px:.10g},"
                    f"{py:.10g},"
                    f"{pz:.10g}\n"
                )

    print(f"\n[SAVE] {csv_path}")

    fig_name=f"peaks_map_{tag}.png"

    fig_path=os.path.join(
        out_dir,
        fig_name
    )

    draw_summary(
        counts_per_ch=counts_per_ch,
        peaks_per_ch =peaks_per_ch,
        x_edges      =x_edges,
        y_edges      =y_edges,
        fig_path     =fig_path,
        draw         =args.show
    )


def main():

    ap=argparse.ArgumentParser(
        description="Find XY peaks per channel from decoded .h5 files"
    )

    ap.add_argument(
        "-i",
        "--input-dir",
        default="./",
        help="directory containing .h5 files or SN subdirs"
    )

    ap.add_argument(
        "-o",
        "--out-dir",
        default=None
    )

    ap.add_argument(
        "--tag",
        default=None
    )

    ap.add_argument(
        "--emin",
        type=float,
        default=EMIN
    )

    ap.add_argument(
        "--emax",
        type=float,
        default=EMAX
    )

    ap.add_argument(
        "--npeaks",
        type=int,
        default=NPEAKS
    )

    ap.add_argument(
        "--show",
        action="store_true"
    )

    args=ap.parse_args()

    root_dir=args.input_dir

    if not os.path.isdir(root_dir):

        raise RuntimeError(
            f"Input directory not found : {root_dir}"
        )

    direct_h5=[

        fn for fn in os.listdir(root_dir)

        if fn.endswith(".h5")
    ]

    #
    # normal mode
    #
    if len(direct_h5)>0:

        process_single_directory(
            root_dir,
            args
        )

        return

    #
    # SN subdir mode
    #
    subdirs=natsorted([

        os.path.join(root_dir,d)

        for d in os.listdir(root_dir)

        if os.path.isdir(os.path.join(root_dir,d))
    ])

    if len(subdirs)==0:

        raise RuntimeError(
            f"No h5 files or subdirs found : {root_dir}"
        )

    for subdir in subdirs:

        print("\n"+"="*80)

        print(f"[PROCESS] {subdir}")

        print("="*80)

        process_single_directory(
            subdir,
            args
        )


if __name__=="__main__":

    main()