#!/usr/bin/env python3
#
# 偏光ビーム試験の簡易解析: 対生成のイベントでシンチ（HA）のヒットが 2 つのものを選び、
# 2 点を結んだ線の方位角（ビーム軸のまわり）の分布から偏光の向きと変調を見る。
#
# raw ファイル（ha<NN>_<番号>.dat）を直接読む。中身の読み取りは psa_decode.py の read_file・decode_events。
# 1 つのランのボードは 2 つの storage（calltt、killua）に分かれるので、両方のディレクトリを -i に渡す。
# ボード同士はトリガーカウンタ（trig）でそろえる。
#
# ヒットの位置はブロック（8x8 の結晶）の中心。対生成の e+ e- は HG が飽和しやすく、
# 飽和したヒットは重心がブロックの中央に寄るので、ピクセルまでは求めない。
# 位置は param/geometry_sm3.csv（EricaNext の calib/psa/geometry.csv の写し）。
#
# ch の番号: ADC 番号 1, 2, 4, 8, 16, 32 → ch 0-5（EricaNext と同じ。2026-11-20 から）。
# psa_decode.py が .h5 に書く ch は 5 - ch（古い番号）なので、ここでは .h5 を使わず raw から直接読む。
#
# 方位角 phi: ビームの向き b に垂直な面に 2 点を結ぶ線を射影し、基準軸（既定 X = 東）の射影から測る。
# 線には向きが無いので 0-180 deg に畳む。分布を N(phi) = A [1 + mu cos(2 (phi - phi0))] で fit する。
# 検出器の形（HA の並び）だけでも分布は変わるので、偏光なし（または偏光の向きを回した）ランを
# --reference に渡して比（N / N_ref）でも fit する。
#
# 例:
#   python3 pair_polarization.py \
#     -i /mnt/nas_eureka/SMILE-3FM_data/20261201/calltt/20261201_run05 \
#        /mnt/nas_eureka/SMILE-3FM_data/20261201/killua/20261201_run05 \
#     --beam-dir 0 0 1 -o pair_run05
#

import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from psa_decode import read_file, decode_events, DEFAULT_TEMP_A, DEFAULT_TEMP_B  # noqa: E402

# このリポジトリ (analysis_sm3_FM) の場所
REPO_DIR = os.path.dirname(os.path.dirname(SRC_DIR))
DEFAULT_GEOMETRY = os.path.join(REPO_DIR, "PSA_analysis", "param", "geometry_sm3.csv")

N_CH = 6
FILE_RE = re.compile(r"^ha(\d+)_(\d+)\.dat$", re.IGNORECASE)


def parse_args():

    p = argparse.ArgumentParser(description="対生成のシンチの 2 ヒットを結んだ線の方位角で偏光を見る簡易解析（raw を直接読む）")
    p.add_argument("-i", "--input", nargs="+", required=True,
                   help="1 つのランの raw のディレクトリ（calltt と killua の両方）か ha<NN>_<番号>.dat。サブディレクトリも探す")
    p.add_argument("--reference", nargs="+", default=None,
                   help="比べるラン（偏光なし、または偏光の向きを回したもの）の raw。N / N_ref でも fit する")
    p.add_argument("-o", "--out", default="pair_polarization", help="出力の名前（<out>.png、<out>_pairs.csv、<out>.json）")
    p.add_argument("--geometry", default=DEFAULT_GEOMETRY, help="ブロックの位置（board,ch,x0,y0,z0,...）[mm]")
    p.add_argument("--beam-dir", type=float, nargs=3, default=[0., 0., 1.],
                   help="ビームの向き（検出器の座標。X: 東、Y: 南、Z: 下）。既定は上から下へ (0, 0, 1)")
    p.add_argument("--ref-axis", type=float, nargs=3, default=[1., 0., 0.],
                   help="方位角の基準の軸（ビームに垂直な面に射影して使う）。既定 X（東）")
    p.add_argument("--hg-min", type=float, default=1000.,
                   help="ヒットとして数える HG の下限 [ADC]（ノイズを落とす）")
    p.add_argument("--peak-id-window", type=int, default=3,
                   help="トリガー時刻のヒットだけ使う: ボードごとの peak_id の最頻値 ± この値。負なら使わない")
    p.add_argument("--n-hits", type=int, default=2, help="選ぶイベントのヒットの数（対生成なら 2）")
    p.add_argument("--min-distance", type=float, default=0.,
                   help="2 点の距離の下限 [mm]（隣のブロックへの漏れを落とすときに使う）")
    p.add_argument("--bins", type=int, default=18, help="方位角のビンの数（0-180 deg）")
    p.add_argument("--max-files", type=int, default=0,
                   help="ボードごとに読むファイルの数の上限（0 = 全部。試しに速く流すとき）")
    p.add_argument("--include-tmp", action="store_true", help="書きかけの .dat.tmp も読む")

    return p.parse_args()


def collect_files(inputs, include_tmp=False, max_files=0):
    """{board: [file, ...]}（番号順。同じ名前のファイルが 2 か所にあれば最初だけ）"""

    found = {}

    for item in inputs:
        if os.path.isdir(item):
            paths = glob.glob(os.path.join(item, "**", "ha*_*.dat*"), recursive=True)
        else:
            paths = glob.glob(item)
        for path in paths:
            name = os.path.basename(path)
            if name.endswith(".tmp"):
                if not include_tmp:
                    continue
                name = name[:-4]
            m = FILE_RE.match(name)
            if not m:
                continue
            key = (int(m.group(1)), int(m.group(2)))
            if key in found:
                print(f"[WARN] 同じファイルが 2 か所にある。最初の方を使う: {found[key]} / {path}")
                continue
            found[key] = path

    boards = {}
    for (board, num), path in sorted(found.items()):
        boards.setdefault(board, []).append(path)
    if max_files > 0:
        boards = {b: f[:max_files] for b, f in boards.items()}

    return boards


def read_hits(boards):
    """全ボードのヒット（DataCount の立った ch）: board, ch（ADC 番号）, trig, hg, lg, peak_id, satu_h, satu_l"""

    cols = {k: [] for k in ("board", "ch", "trig", "hg", "lg", "peak_id", "satu_h", "satu_l")}
    n_frames = {}

    for board, files in boards.items():
        seen = set()
        n_frames[board] = 0
        for path in files:
            t0 = time.time()
            buf, sizes, n = read_file(path)
            if buf is None or n == 0:
                print(f"[WARN] イベントが無い: {path}")
                continue
            (trig, clk, hit, temp_raw, temp_corr, ped, wave_sum, hg, lg, x_hg, y_hg, x_lg, y_lg,
             peak_id, hitnum, satu_h, satu_l, adc) = decode_events(buf, sizes, False, DEFAULT_TEMP_A, DEFAULT_TEMP_B, False)
            # 同じボードの同じ trig は 1 回だけ（ファイルの重なり）
            keep = np.array([t not in seen for t in trig.tolist()], dtype=bool)
            seen.update(trig[keep].tolist())
            n_frames[board] += int(keep.sum())
            ev, ch = np.nonzero((hit == 1) & keep[:, None])   # ch: decode_events の中の番号 = ADC 番号の順
            cols["board"].append(np.full(len(ev), board, dtype=np.int16))
            cols["ch"].append(ch.astype(np.int16))
            cols["trig"].append(trig[ev].astype(np.int64))
            cols["hg"].append(hg[ev, ch])
            cols["lg"].append(lg[ev, ch])
            cols["peak_id"].append(peak_id[ev, ch])
            cols["satu_h"].append(satu_h[ev, ch])
            cols["satu_l"].append(satu_l[ev, ch])
            print(f"ha{board:02d} {os.path.basename(path)}: {n} frames, {len(ev)} hits ({time.time() - t0:.1f} s)", flush=True)

    hits = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in cols.items()}

    return hits, n_frames


def select_on_time(hits, hg_min, window):
    """HG の下限と、ボードごとの peak_id の最頻値 ± window"""

    m = hits["hg"] >= hg_min
    windows = {}

    if window >= 0:
        on = np.zeros(len(m), dtype=bool)
        for b in np.unique(hits["board"][m]):
            mb = m & (hits["board"] == b)
            mode = int(np.bincount(hits["peak_id"][mb]).argmax())
            windows[int(b)] = (mode - window, mode + window)
            on |= mb & (np.abs(hits["peak_id"] - mode) <= window)
        m = on

    return m, windows


def load_geometry(path):

    rows = [l for l in open(path) if l.strip() and not l.startswith("#")]
    head = rows[0].strip().split(",")
    geo = {}
    for line in rows[1:]:
        v = dict(zip(head, line.strip().split(",")))
        geo[(int(v["board"]), int(v["ch"]))] = np.array([float(v["x0"]), float(v["y0"]), float(v["z0"])])

    return geo


def build_pairs(hits, mask, geo, n_hits, min_distance):
    """trig ごとにまとめ、ヒットが n_hits 個（geometry にあるもの）のイベントの 2 点"""

    idx = np.flatnonzero(mask)
    in_geo = np.array([(int(hits["board"][i]), int(hits["ch"][i])) in geo for i in idx], dtype=bool)
    missing = sorted(set((int(hits["board"][i]), int(hits["ch"][i])) for i in idx[~in_geo]))
    if missing:
        print(f"[WARN] geometry に無いブロックのヒットは使わない: {missing}")
    idx = idx[in_geo]

    order = np.argsort(hits["trig"][idx], kind="stable")
    idx = idx[order]
    trig = hits["trig"][idx]
    starts = np.flatnonzero(np.r_[True, trig[1:] != trig[:-1]])
    counts = np.diff(np.r_[starts, len(trig)])
    multiplicity = np.bincount(counts, minlength=8)

    pairs = []
    for s in starts[counts == n_hits]:
        group = idx[s:s + n_hits]
        if n_hits != 2:
            continue
        i, j = group
        p1 = geo[(int(hits["board"][i]), int(hits["ch"][i]))]
        p2 = geo[(int(hits["board"][j]), int(hits["ch"][j]))]
        if np.linalg.norm(p2 - p1) < min_distance:
            continue
        pairs.append((int(hits["trig"][i]), i, j, p1, p2))

    return pairs, multiplicity


def azimuth(pairs, beam, ref):
    """ビームに垂直な面での 2 点を結ぶ線の方位角 [deg]（0-180）"""

    b = np.asarray(beam, float)
    b /= np.linalg.norm(b)
    e1 = np.asarray(ref, float) - np.dot(ref, b) * b
    if np.linalg.norm(e1) < 1e-9:
        raise SystemExit("--ref-axis がビームと平行")
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(b, e1)
    d = np.array([p2 - p1 for _, _, _, p1, p2 in pairs]).reshape(-1, 3)
    phi = np.degrees(np.arctan2(d @ e2, d @ e1)) % 180.

    return phi, e1, e2


def fit_modulation(centres, counts, errors):
    """N(phi) = A [1 + mu cos(2 (phi - phi0))] を重み付き最小 2 乗で（線形: a + c cos 2phi + s sin 2phi）"""

    ok = errors > 0
    if ok.sum() < 3 or np.sum(counts[ok]) <= 0:
        return None
    t = np.radians(2 * centres[ok])
    X = np.column_stack([np.ones(ok.sum()), np.cos(t), np.sin(t)])
    w = 1 / errors[ok] ** 2
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    a, c, s = cov @ (X.T @ (w * counts[ok]))
    chi2 = float(np.sum(w * (counts[ok] - X @ np.array([a, c, s])) ** 2))
    amp = np.hypot(c, s)
    mu = amp / a
    phi0 = (np.degrees(np.arctan2(s, c)) / 2) % 180.
    # 誤差の伝播
    J_mu = np.array([-amp / a ** 2, c / (a * amp), s / (a * amp)])
    J_phi = np.degrees(np.array([0., -s / amp ** 2, c / amp ** 2])) / 2

    return {"A": float(a), "mu": float(mu), "mu_err": float(np.sqrt(J_mu @ cov @ J_mu)),
            "phi0_deg": float(phi0), "phi0_err_deg": float(np.sqrt(J_phi @ cov @ J_phi)),
            "chi2": chi2, "ndf": int(ok.sum() - 3)}


def analyse(inputs, args, geo, label):

    boards = collect_files(inputs, args.include_tmp, args.max_files)
    if not boards:
        raise SystemExit(f"raw ファイル（ha<NN>_<番号>.dat）が無い: {inputs}")
    print(f"[{label}] boards: " + ", ".join(f"ha{b:02d}({len(f)})" for b, f in boards.items()))
    hits, n_frames = read_hits(boards)
    mask, windows = select_on_time(hits, args.hg_min, args.peak_id_window)
    pairs, multiplicity = build_pairs(hits, mask, geo, args.n_hits, args.min_distance)
    phi, e1, e2 = azimuth(pairs, args.beam_dir, args.ref_axis)
    print(f"[{label}] hits {len(hits['trig'])}, after HG / peak_id {int(mask.sum())}, events with {args.n_hits} hits: {len(pairs)}")

    return dict(boards=boards, hits=hits, n_frames=n_frames, mask=mask, windows=windows, pairs=pairs,
                multiplicity=multiplicity, phi=phi, e1=e1, e2=e2)


def main():

    args = parse_args()
    geo = load_geometry(args.geometry)
    run = analyse(args.input, args, geo, "run")
    ref = analyse(args.reference, args, geo, "reference") if args.reference else None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    edges = np.linspace(0, 180, args.bins + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    n, _ = np.histogram(run["phi"], edges)
    fit = fit_modulation(centres, n.astype(float), np.sqrt(np.maximum(n, 1.)))
    result = {"input": args.input, "reference": args.reference, "geometry": args.geometry,
              "beam_dir": args.beam_dir, "ref_axis_projected": run["e1"].tolist(), "second_axis": run["e2"].tolist(),
              "hg_min": args.hg_min, "peak_id_windows": run["windows"], "frames_per_board": run["n_frames"],
              "multiplicity_after_cuts": run["multiplicity"].tolist(), "pairs": len(run["pairs"]),
              "fit_counts": fit}

    if ref:
        nr, _ = np.histogram(ref["phi"], edges)
        good = (n > 0) & (nr > 0)
        ratio = np.zeros(len(n))
        err = np.zeros(len(n))
        scale = nr.sum() / max(n.sum(), 1)
        ratio[good] = n[good] / nr[good] * scale
        err[good] = ratio[good] * np.sqrt(1 / n[good] + 1 / nr[good])
        result["fit_ratio"] = fit_modulation(centres, ratio, err)
        result["reference_pairs"] = len(ref["pairs"])

    # 図: 方位角の分布と fit、（あれば）比、選んだペアのヒットの位置、ヒットの数の分布
    ncol = 4 if ref else 3
    fig, axs = plt.subplots(1, ncol, figsize=(5.5 * ncol, 4.8))
    ax = axs[0]
    ax.errorbar(centres, n, yerr=np.sqrt(n), fmt="o", ms=4, color="k", label=f"run: {len(run['phi'])} pairs")
    xx = np.linspace(0, 180, 361)
    if fit:
        ax.plot(xx, fit["A"] * (1 + fit["mu"] * np.cos(np.radians(2 * (xx - fit["phi0_deg"])))), "r-",
                label=f"mu = {fit['mu']:.3f} +- {fit['mu_err']:.3f}\nphi0 = {fit['phi0_deg']:.1f} +- {fit['phi0_err_deg']:.1f} deg"
                      f"\nchi2/ndf = {fit['chi2']:.1f}/{fit['ndf']}")
    if ref:
        ax.step(centres, nr * n.sum() / max(nr.sum(), 1), where="mid", color="0.5", label="reference (scaled)")
    ax.set(xlabel="azimuth of the line between the 2 hits [deg]", ylabel="pairs / bin", xlim=(0, 180), ylim=(0, None),
           title="azimuth (0 = reference axis, around the beam)")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    col = 1
    if ref:
        ax = axs[1]
        fr = result["fit_ratio"]
        ax.errorbar(centres, ratio, yerr=err, fmt="o", ms=4, color="k")
        if fr:
            ax.plot(xx, fr["A"] * (1 + fr["mu"] * np.cos(np.radians(2 * (xx - fr["phi0_deg"])))), "r-",
                    label=f"mu = {fr['mu']:.3f} +- {fr['mu_err']:.3f}\nphi0 = {fr['phi0_deg']:.1f} +- {fr['phi0_err_deg']:.1f} deg")
            ax.legend(fontsize=8)
        ax.set(xlabel="azimuth [deg]", ylabel="N / N_ref (scaled to the same total)", xlim=(0, 180), title="run / reference")
        ax.grid(alpha=.3)
        col = 2
    ax = axs[col]
    if run["pairs"]:
        P = np.array([[*p1, *p2] for _, _, _, p1, p2 in run["pairs"]])
        allp = np.vstack([P[:, :3], P[:, 3:]])
        blocks, cnt = np.unique(np.round(allp, 1), axis=0, return_counts=True)
        sc = ax.scatter(blocks[:, 0], blocks[:, 1], c=cnt, s=40 + 160 * (blocks[:, 2] > 330), cmap="viridis",
                        marker="s", edgecolor="k", lw=.3)
        plt.colorbar(sc, ax=ax, label="hits in selected pairs")
    ax.set(xlabel="X [mm] (E)", ylabel="Y [mm] (S)", title="blocks of the selected pairs (large: bottom)")
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="datalim")
    ax = axs[col + 1]
    mult = run["multiplicity"]
    ax.bar(np.arange(len(mult)), mult, color="0.4")
    ax.set(xlabel="hits per trigger (after HG / peak_id)", ylabel="triggers", yscale="log", title="multiplicity")
    ax.grid(alpha=.3)
    fig.suptitle(f"{args.out}: beam {tuple(args.beam_dir)}, reference axis -> {np.round(run['e1'], 3).tolist()}, "
                 f"HG >= {args.hg_min:g}, peak_id window {args.peak_id_window}")
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=120)

    # ペアの一覧
    h = run["hits"]
    with open(f"{args.out}_pairs.csv", "w") as f:
        f.write("# ch = ADC number (1, 2, 4, 8, 16, 32 -> 0-5); positions are block centres [mm]; phi [deg]\n")
        f.write("trig,board1,ch1,hg1,satu_h1,x1,y1,z1,board2,ch2,hg2,satu_h2,x2,y2,z2,phi\n")
        for (trig, i, j, p1, p2), ph in zip(run["pairs"], run["phi"]):
            f.write(f"{trig},{h['board'][i]},{h['ch'][i]},{h['hg'][i]:.1f},{h['satu_h'][i]},{p1[0]:.2f},{p1[1]:.2f},{p1[2]:.2f},"
                    f"{h['board'][j]},{h['ch'][j]},{h['hg'][j]:.1f},{h['satu_h'][j]},{p2[0]:.2f},{p2[1]:.2f},{p2[2]:.2f},{ph:.2f}\n")
    with open(f"{args.out}.json", "w") as f:
        json.dump(result, f, indent=2, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(json.dumps({k: result[k] for k in ("pairs", "fit_counts") + (("fit_ratio",) if ref else ())}, indent=2))
    print(f"saved {args.out}.png, {args.out}_pairs.csv, {args.out}.json")


if __name__ == "__main__":
    main()
