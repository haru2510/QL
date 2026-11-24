#!/usr/bin/env python3
#
# 偏光ビーム試験の簡易解析: 対生成のイベントでシンチ（HA）のヒットが 2 つのものを選び、
# 2 点を結んだ線の方位角（ビーム軸のまわり）の分布から偏光の向きと変調を見る。
#
# raw ファイル（ha<NN>_<番号>.dat）を直接読む。中身の読み取りは PSA_analysis/src_onoda/psa_decode.py の read_file・decode_events。
# 1 つのランのボードは 2 つの storage（calltt、killua）に分かれるので、両方のディレクトリを読む。
# ボード同士はラン 1 つの中でトリガーカウンタ（trig）でそろえる。ボードの中で trig が戻ったら epoch を 1 つ進め、
# (epoch, trig) で組にする（EricaNext の PSA_ReadData と同じ）。複数のランは、ランごとに組にしてから方位角を足す。
#
# 入力と条件は、既定で ETCC_analysis/param/pair_polarization.json から読む（ラン、比べるラン、ビームの向き、カットなど）。
# コマンドラインで渡したものはそちらが優先する。-i / --runs を渡すと、JSON のラン（と比べるラン）は使わない。
#
# ヒットの位置はブロック（8x8 の結晶）の中心。対生成の e+ e- は HG が飽和しやすく、
# 飽和したヒットは重心がブロックの中央に寄るので、ピクセルまでは求めない。
# 位置は PSA_analysis/param/geometry_sm3.csv（EricaNext の calib/psa/geometry.csv の写し）。
#
# ch の番号: ADC 番号 1, 2, 4, 8, 16, 32 → ch 0-5（EricaNext と同じ。2026-11-20 から）。
# psa_decode.py が .h5 に書く ch は 5 - ch（古い番号）なので、ここでは .h5 を使わず raw から直接読む。
#
# 方位角 phi: ビームの向き b に垂直な面に 2 点を結ぶ線を射影し、基準軸（既定 X = 東）の射影から測る。
# 線には向きが無いので 0-180 deg に畳む。分布を N(phi) = A [1 + mu cos(2 (phi - phi0))] で fit する。
# 検出器の形（HA の並び）だけでも分布は変わるので、偏光なし（または偏光の向きを回した）ランを
# 比べるラン（reference）にして、比（N / N_ref）でも fit する。
#
# 例:
#   python3 ETCC_analysis/pair_polarization.py                                    # param/pair_polarization.json のとおり
#   python3 ETCC_analysis/pair_polarization.py --config my_beamtest.json
#   python3 ETCC_analysis/pair_polarization.py --runs 20261201_run05 20261201_run07 --reference-runs 20261201_run06
#   python3 ETCC_analysis/pair_polarization.py -i <calltt のラン> <killua のラン> -o pair_run05   # ディレクトリを直接（1 ラン）
#

import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

ETCC_DIR = os.path.dirname(os.path.abspath(__file__))
# このリポジトリ (analysis_sm3_FM) の場所。JSON の中の相対パスはここから
REPO_DIR = os.path.dirname(ETCC_DIR)
PSA_SRC_DIR = os.path.join(REPO_DIR, "PSA_analysis", "src_onoda")
if PSA_SRC_DIR not in sys.path:
    sys.path.insert(0, PSA_SRC_DIR)

from psa_decode import read_file, decode_events, DEFAULT_TEMP_A, DEFAULT_TEMP_B  # noqa: E402

DEFAULT_CONFIG = os.path.join(ETCC_DIR, "param", "pair_polarization.json")

N_CH = 6
FILE_RE = re.compile(r"^ha(\d+)_(\d+)\.dat$", re.IGNORECASE)
RUN_RE = re.compile(r"^(\d{8})_(run\d+)$")

# JSON にもコマンドラインにも無いときの値
DEFAULTS = {
    "nas": "/mnt/nas_eureka/SMILE-3FM_data",
    "storages": ["calltt", "killua"],
    "runs": [],
    "reference_runs": [],
    "out": "ETCC_analysis/output/pair_polarization",
    "geometry": "PSA_analysis/param/geometry_sm3.csv",
    "beam_dir": [0., 0., 1.],
    "ref_axis": [1., 0., 0.],
    "hg_min": 1000.,
    "peak_id_window": 3,
    "n_hits": 2,
    "min_distance": 0.,
    "bins": 18,
    "max_files": 0,
    "include_tmp": False,
}


def parse_args():

    p = argparse.ArgumentParser(description="対生成のシンチの 2 ヒットを結んだ線の方位角で偏光を見る簡易解析（raw を直接読む）")
    p.add_argument("--config", default=None,
                   help=f"入力と条件の JSON（既定 {os.path.relpath(DEFAULT_CONFIG, REPO_DIR)}。無ければ使わない）")
    p.add_argument("-i", "--input", nargs="+", default=None,
                   help="1 つのランの raw のディレクトリ（calltt と killua の両方）か ha<NN>_<番号>.dat。JSON のランの代わりに使う")
    p.add_argument("--reference", nargs="+", default=None,
                   help="-i と一緒に: 比べるラン（1 つ）の raw のディレクトリ")
    p.add_argument("--runs", nargs="+", default=None,
                   help="<日付>_runNN（例 20261201_run05）。NAS の calltt と killua を読む。複数なら足す。JSON のランの代わりに使う")
    p.add_argument("--reference-runs", nargs="+", default=None, help="--runs と同じ形で、比べるラン")
    p.add_argument("-o", "--out", default=None, help="出力の名前（<out>.png、<out>_pairs.csv、<out>.json）")
    p.add_argument("--geometry", default=None, help="ブロックの位置（board,ch,x0,y0,z0,...）[mm]")
    p.add_argument("--beam-dir", type=float, nargs=3, default=None,
                   help="ビームの向き（検出器の座標。X: 東、Y: 南、Z: 下）。既定は上から下へ (0, 0, 1)")
    p.add_argument("--ref-axis", type=float, nargs=3, default=None,
                   help="方位角の基準の軸（ビームに垂直な面に射影して使う）。既定 X（東）")
    p.add_argument("--hg-min", type=float, default=None, help="ヒットとして数える HG の下限 [ADC]（既定 1000）")
    p.add_argument("--peak-id-window", type=int, default=None,
                   help="トリガー時刻のヒットだけ使う: ボードごとの peak_id の最頻値 ± この値（既定 3）。負なら使わない")
    p.add_argument("--n-hits", type=int, default=None, help="選ぶイベントのヒットの数（既定 2）")
    p.add_argument("--min-distance", type=float, default=None, help="2 点の距離の下限 [mm]（既定 0）")
    p.add_argument("--bins", type=int, default=None, help="方位角のビンの数（0-180 deg、既定 18）")
    p.add_argument("--max-files", type=int, default=None,
                   help="ボードごとに読むファイルの数の上限（0 = 全部。試しに速く流すとき）")
    p.add_argument("--include-tmp", action="store_true", default=None, help="書きかけの .dat.tmp も読む")

    return p.parse_args()


def in_repo(path):
    """JSON の中の相対パスはリポジトリから"""
    return path if os.path.isabs(path) else os.path.join(REPO_DIR, path)


def load_settings(args):
    """DEFAULTS < JSON < コマンドライン"""

    s = dict(DEFAULTS)
    s["geometry"] = in_repo(s["geometry"])
    s["out"] = in_repo(s["out"])
    config = args.config or (DEFAULT_CONFIG if os.path.isfile(DEFAULT_CONFIG) else None)
    if config:
        if not os.path.isfile(config):
            raise SystemExit(f"設定ファイルが無い: {config}")
        with open(config) as f:
            js = json.load(f)
        unknown = set(js) - set(DEFAULTS) - {"_doc", "_units"}
        if unknown:
            raise SystemExit(f"{config}: 知らないキー {sorted(unknown)}")
        for k, v in js.items():
            if k.startswith("_"):
                continue
            s[k] = in_repo(v) if k in ("out", "geometry") and v else v
        print(f"config: {config}")
    s["config"] = config

    for k in ("out", "geometry", "beam_dir", "ref_axis", "hg_min", "peak_id_window", "n_hits",
              "min_distance", "bins", "max_files", "include_tmp"):
        v = getattr(args, k)
        if v is not None:
            s[k] = v
    if args.input is not None:
        s["runs"] = [{"name": "input", "dirs": args.input}]
        s["reference_runs"] = [{"name": "reference", "dirs": args.reference}] if args.reference else []
    elif args.runs is not None:
        s["runs"] = list(args.runs)
        s["reference_runs"] = list(args.reference_runs or [])
    elif args.reference_runs is not None:
        s["reference_runs"] = list(args.reference_runs)
    if not s["runs"]:
        raise SystemExit("ランが無い: ETCC_analysis/param/pair_polarization.json の runs、--runs か -i で渡す")

    return s


def resolve_run(spec, s):
    """ラン 1 つ: "20261201_run05"、{"date": ..., "run": ...}、{"dirs": [...]} → (名前, [ディレクトリ])"""

    if isinstance(spec, str):
        m = RUN_RE.match(spec)
        if not m:
            raise SystemExit(f"ランは <日付>_runNN の形で書く: {spec}")
        spec = {"date": m.group(1), "run": m.group(2)}
    if "dirs" in spec:
        dirs = [in_repo(d) if not os.path.isabs(d) else d for d in spec["dirs"]]
        name = spec.get("name") or os.path.basename(os.path.normpath(dirs[0]))
    else:
        date, run = str(spec["date"]), str(spec["run"])
        dirs = [os.path.join(s["nas"], date, st, f"{date}_{run}") for st in s["storages"]]
        name = spec.get("name", f"{date}_{run}")

    return name, dirs


def collect_files(inputs, include_tmp=False, max_files=0):
    """{board: [file, ...]}（番号順。同じ名前のファイルが 2 か所にあれば最初だけ）"""

    found = {}

    for item in inputs:
        if os.path.isdir(item):
            paths = glob.glob(os.path.join(item, "**", "ha*_*.dat*"), recursive=True)
        else:
            if not glob.glob(item):
                print(f"[WARN] 無い: {item}")
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
    """全ボードのヒット（DataCount の立った ch）: board, ch（ADC 番号）, trig, epoch, event, hg, lg, peak_id, satu_h, satu_l
    event = epoch * 2^32 + trig。trig がボードの中で戻ったら epoch を 1 つ進める"""

    cols = {k: [] for k in ("board", "ch", "trig", "epoch", "hg", "lg", "peak_id", "satu_h", "satu_l")}
    n_frames, n_resets = {}, {}

    for board, files in boards.items():
        seen = set()
        n_frames[board] = 0
        epoch, prev = 0, None
        n_resets[board] = 0
        for path in files:
            t0 = time.time()
            buf, sizes, n = read_file(path)
            if buf is None or n == 0:
                print(f"[WARN] イベントが無い: {path}")
                continue
            (trig, clk, hit, temp_raw, temp_corr, ped, wave_sum, hg, lg, x_hg, y_hg, x_lg, y_lg,
             peak_id, hitnum, satu_h, satu_l, adc) = decode_events(buf, sizes, False, DEFAULT_TEMP_A, DEFAULT_TEMP_B, False)
            t = trig.astype(np.int64)
            back = np.r_[prev is not None and t[0] < prev, t[1:] < t[:-1]]
            ep = epoch + np.cumsum(back)
            n_resets[board] += int(back.sum())
            epoch, prev = int(ep[-1]), int(t[-1])
            key = ep * (1 << 32) + t
            # 同じボードの同じ (epoch, trig) は 1 回だけ（ファイルの重なり）
            keep = np.array([k not in seen for k in key.tolist()], dtype=bool)
            seen.update(key[keep].tolist())
            n_frames[board] += int(keep.sum())
            ev, ch = np.nonzero((hit == 1) & keep[:, None])   # ch: decode_events の中の番号 = ADC 番号の順
            cols["board"].append(np.full(len(ev), board, dtype=np.int16))
            cols["ch"].append(ch.astype(np.int16))
            cols["trig"].append(t[ev])
            cols["epoch"].append(ep[ev].astype(np.int64))
            cols["hg"].append(hg[ev, ch])
            cols["lg"].append(lg[ev, ch])
            cols["peak_id"].append(peak_id[ev, ch])
            cols["satu_h"].append(satu_h[ev, ch])
            cols["satu_l"].append(satu_l[ev, ch])
            print(f"ha{board:02d} {os.path.basename(path)}: {n} frames, {len(ev)} hits ({time.time() - t0:.1f} s)", flush=True)
        if n_resets[board]:
            print(f"[WARN] ha{board:02d}: trig が {n_resets[board]} 回戻った（epoch で分けて組にする）")

    hits = {k: (np.concatenate(v) if v else np.zeros(0, dtype=np.int64)) for k, v in cols.items()}
    hits["event"] = hits["epoch"] * (1 << 32) + hits["trig"]

    return hits, n_frames, n_resets


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
    """(epoch, trig) ごとにまとめ、ヒットが n_hits 個（geometry にあるもの）のイベントの 2 点"""

    idx = np.flatnonzero(mask)
    in_geo = np.array([(int(hits["board"][i]), int(hits["ch"][i])) in geo for i in idx], dtype=bool)
    missing = sorted(set((int(hits["board"][i]), int(hits["ch"][i])) for i in idx[~in_geo]))
    if missing:
        print(f"[WARN] geometry に無いブロックのヒットは使わない: {missing}")
    idx = idx[in_geo]

    order = np.argsort(hits["event"][idx], kind="stable")
    idx = idx[order]
    event = hits["event"][idx]
    starts = np.flatnonzero(np.r_[True, event[1:] != event[:-1]]) if len(event) else np.zeros(0, dtype=int)
    counts = np.diff(np.r_[starts, len(event)])
    multiplicity = np.bincount(counts, minlength=8) if len(counts) else np.zeros(8, dtype=int)

    pairs = []
    if n_hits == 2:
        for s in starts[counts == 2]:
            i, j = idx[s:s + 2]
            p1 = geo[(int(hits["board"][i]), int(hits["ch"][i]))]
            p2 = geo[(int(hits["board"][j]), int(hits["ch"][j]))]
            if np.linalg.norm(p2 - p1) < min_distance:
                continue
            pairs.append((int(hits["epoch"][i]), int(hits["trig"][i]), i, j, p1, p2))

    return pairs, multiplicity


def axes(beam, ref):
    """ビームに垂直な面の 2 本の軸 e1（ref の射影、phi = 0）、e2 = b x e1"""

    b = np.asarray(beam, float)
    b /= np.linalg.norm(b)
    e1 = np.asarray(ref, float) - np.dot(ref, b) * b
    if np.linalg.norm(e1) < 1e-9:
        raise SystemExit("ref_axis がビームと平行")
    e1 /= np.linalg.norm(e1)

    return e1, np.cross(b, e1)


def azimuth(pairs, e1, e2):
    """ビームに垂直な面での 2 点を結ぶ線の方位角 [deg]（0-180）"""

    d = np.array([p2 - p1 for *_, p1, p2 in pairs]).reshape(-1, 3)

    return np.degrees(np.arctan2(d @ e2, d @ e1)) % 180.


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
    if a <= 0 or amp == 0:
        return None
    mu = amp / a
    phi0 = (np.degrees(np.arctan2(s, c)) / 2) % 180.
    # 誤差の伝播
    J_mu = np.array([-amp / a ** 2, c / (a * amp), s / (a * amp)])
    J_phi = np.degrees(np.array([0., -s / amp ** 2, c / amp ** 2])) / 2

    return {"A": float(a), "mu": float(mu), "mu_err": float(np.sqrt(J_mu @ cov @ J_mu)),
            "phi0_deg": float(phi0), "phi0_err_deg": float(np.sqrt(J_phi @ cov @ J_phi)),
            "chi2": chi2, "ndf": int(ok.sum() - 3)}


def analyse_run(name, dirs, s, geo, e1, e2):
    """ラン 1 つ: ヒットを読み、(epoch, trig) で組にし、方位角を出す"""

    boards = collect_files(dirs, s["include_tmp"], s["max_files"])
    if not boards:
        raise SystemExit(f"[{name}] raw ファイル（ha<NN>_<番号>.dat）が無い: {dirs}")
    print(f"[{name}] boards: " + ", ".join(f"ha{b:02d}({len(f)})" for b, f in boards.items()))
    hits, n_frames, n_resets = read_hits(boards)
    mask, windows = select_on_time(hits, s["hg_min"], s["peak_id_window"])
    pairs, multiplicity = build_pairs(hits, mask, geo, s["n_hits"], s["min_distance"])
    phi = azimuth(pairs, e1, e2)
    print(f"[{name}] hits {len(hits['trig'])}, after HG / peak_id {int(mask.sum())}, events with {s['n_hits']} hits: {len(pairs)}")

    return dict(name=name, dirs=dirs, hits=hits, n_frames=n_frames, n_resets=n_resets, windows=windows,
                pairs=pairs, multiplicity=multiplicity, phi=phi)


def main():

    args = parse_args()
    s = load_settings(args)
    geo = load_geometry(s["geometry"])
    e1, e2 = axes(s["beam_dir"], s["ref_axis"])
    runs = [analyse_run(*resolve_run(r, s), s, geo, e1, e2) for r in s["runs"]]
    refs = [analyse_run(*resolve_run(r, s), s, geo, e1, e2) for r in s["reference_runs"]]
    out = s["out"]
    if os.path.dirname(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    edges = np.linspace(0, 180, s["bins"] + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    hist = lambda phis: np.histogram(np.concatenate(phis) if phis else np.zeros(0), edges)[0]
    n = hist([r["phi"] for r in runs])
    fit = fit_modulation(centres, n.astype(float), np.sqrt(np.maximum(n, 1.)))

    def run_summary(r):
        h = hist([r["phi"]])
        return {"name": r["name"], "dirs": r["dirs"], "pairs": len(r["pairs"]), "frames_per_board": r["n_frames"],
                "trig_resets_per_board": r["n_resets"], "peak_id_windows": r["windows"],
                "multiplicity_after_cuts": r["multiplicity"].tolist(),
                "fit_counts": fit_modulation(centres, h.astype(float), np.sqrt(np.maximum(h, 1.)))}

    result = {"config": s["config"], "settings": {k: v for k, v in s.items() if k not in ("config",)},
              "ref_axis_projected": e1.tolist(), "second_axis": e2.tolist(),
              "runs": [run_summary(r) for r in runs], "reference_runs": [run_summary(r) for r in refs],
              "pairs": int(n.sum()), "fit_counts": fit}

    nr = None
    if refs:
        nr = hist([r["phi"] for r in refs])
        good = (n > 0) & (nr > 0)
        ratio = np.zeros(len(n))
        err = np.zeros(len(n))
        scale = nr.sum() / max(n.sum(), 1)
        ratio[good] = n[good] / nr[good] * scale
        err[good] = ratio[good] * np.sqrt(1 / n[good] + 1 / nr[good])
        result["reference_pairs"] = int(nr.sum())
        result["fit_ratio"] = fit_modulation(centres, ratio, err)

    # 図: 方位角の分布と fit、（あれば）比、選んだペアのヒットの位置、ヒットの数の分布
    ncol = 4 if refs else 3
    fig, axs = plt.subplots(1, ncol, figsize=(5.5 * ncol, 4.8))
    ax = axs[0]
    ax.errorbar(centres, n, yerr=np.sqrt(n), fmt="o", ms=4, color="k", label=f"runs: {int(n.sum())} pairs")
    xx = np.linspace(0, 180, 361)
    if fit:
        ax.plot(xx, fit["A"] * (1 + fit["mu"] * np.cos(np.radians(2 * (xx - fit["phi0_deg"])))), "r-",
                label=f"mu = {fit['mu']:.3f} +- {fit['mu_err']:.3f}\nphi0 = {fit['phi0_deg']:.1f} +- {fit['phi0_err_deg']:.1f} deg"
                      f"\nchi2/ndf = {fit['chi2']:.1f}/{fit['ndf']}")
    if refs:
        ax.step(centres, nr * n.sum() / max(nr.sum(), 1), where="mid", color="0.5", label="reference (scaled)")
    ax.set(xlabel="azimuth of the line between the 2 hits [deg]", ylabel="pairs / bin", xlim=(0, 180), ylim=(0, None),
           title="azimuth (0 = reference axis, around the beam)")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    col = 1
    if refs:
        ax = axs[1]
        fr = result["fit_ratio"]
        ax.errorbar(centres, ratio, yerr=err, fmt="o", ms=4, color="k")
        if fr:
            ax.plot(xx, fr["A"] * (1 + fr["mu"] * np.cos(np.radians(2 * (xx - fr["phi0_deg"])))), "r-",
                    label=f"mu = {fr['mu']:.3f} +- {fr['mu_err']:.3f}\nphi0 = {fr['phi0_deg']:.1f} +- {fr['phi0_err_deg']:.1f} deg")
            ax.legend(fontsize=8)
        ax.set(xlabel="azimuth [deg]", ylabel="N / N_ref (scaled to the same total)", xlim=(0, 180), title="runs / reference")
        ax.grid(alpha=.3)
        col = 2
    ax = axs[col]
    allpairs = [p for r in runs for p in r["pairs"]]
    if allpairs:
        P = np.array([[*p1, *p2] for *_, p1, p2 in allpairs])
        allp = np.vstack([P[:, :3], P[:, 3:]])
        blocks, cnt = np.unique(np.round(allp, 1), axis=0, return_counts=True)
        sc = ax.scatter(blocks[:, 0], blocks[:, 1], c=cnt, s=40 + 160 * (blocks[:, 2] > 330), cmap="viridis",
                        marker="s", edgecolor="k", lw=.3)
        plt.colorbar(sc, ax=ax, label="hits in selected pairs")
    ax.set(xlabel="X [mm] (E)", ylabel="Y [mm] (S)", title="blocks of the selected pairs (large: bottom)")
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="datalim")
    ax = axs[col + 1]
    mult = np.sum([r["multiplicity"][:8] for r in runs], axis=0)
    ax.bar(np.arange(len(mult)), mult, color="0.4")
    ax.set(xlabel="hits per trigger (after HG / peak_id)", ylabel="triggers", yscale="log", title="multiplicity")
    ax.grid(alpha=.3)
    names = ", ".join(r["name"] for r in runs) + ((" vs " + ", ".join(r["name"] for r in refs)) if refs else "")
    fig.suptitle(f"{os.path.basename(out)}: {names}\nbeam {tuple(s['beam_dir'])}, reference axis -> {np.round(e1, 3).tolist()}, "
                 f"HG >= {s['hg_min']:g}, peak_id window {s['peak_id_window']}")
    fig.tight_layout()
    fig.savefig(f"{out}.png", dpi=120)

    # ペアの一覧
    with open(f"{out}_pairs.csv", "w") as f:
        f.write("# ch = ADC number (1, 2, 4, 8, 16, 32 -> 0-5); positions are block centres [mm]; phi [deg]\n")
        f.write("run,epoch,trig,board1,ch1,hg1,satu_h1,x1,y1,z1,board2,ch2,hg2,satu_h2,x2,y2,z2,phi\n")
        for r, is_ref in [(r, False) for r in runs] + [(r, True) for r in refs]:
            h = r["hits"]
            for (ep, trig, i, j, p1, p2), ph in zip(r["pairs"], r["phi"]):
                f.write(f"{r['name']}{' (reference)' if is_ref else ''},{ep},{trig},"
                        f"{h['board'][i]},{h['ch'][i]},{h['hg'][i]:.1f},{h['satu_h'][i]},{p1[0]:.2f},{p1[1]:.2f},{p1[2]:.2f},"
                        f"{h['board'][j]},{h['ch'][j]},{h['hg'][j]:.1f},{h['satu_h'][j]},{p2[0]:.2f},{p2[1]:.2f},{p2[2]:.2f},{ph:.2f}\n")
    with open(f"{out}.json", "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(json.dumps({"pairs": result["pairs"], "per_run": {r["name"]: r["pairs"] for r in result["runs"]},
                      "fit_counts": fit, **({"fit_ratio": result["fit_ratio"]} if refs else {})}, indent=2))
    print(f"saved {out}.png, {out}_pairs.csv, {out}.json")


if __name__ == "__main__":
    main()
