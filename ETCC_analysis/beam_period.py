#!/usr/bin/env python3
#
# ビーム試験: イベントの時刻の周期を調べ、ビームのタイミング（パルスの周期と位相）に同期をとる。
#
# 時刻は各ボードのカウンタ（ボードごとに原点が違うので、ボードごとに解析する）:
#   HA  : ha<NN>_<番号>.dat の clk（1 MHz、32 bit）
#   TPC : tpc_enc_[ac]<NN>_<番号>.raw の N_time（100 us、32 bit）
# カウンタが 2^32 で戻ったら足し直す（つながった時刻にする）。ランが複数なら、ランごとに解析する。
#
# 手順（ボードごと）:
#   1. 時間間隔 Δt の分布、レートの時間変化
#   2. パワースペクトル: 時刻を bin で数え、segment 秒ずつ FFT して平均する（Leahy 規格化。雑音の平均は 2）。
#      spectra に (bin, segment) を並べると、速い周期（パルス）と遅い周期（スピル）を別々に見られる。
#      いちばん高いピークを取り、その 1/k（k = 2..max_subharmonic）にも有意なピークがあれば基本周波数はそちら。
#      frequency を与えればそれを使う（探さない）。周期が min_period_ticks カウントより短いものは探さない
#      （カウンタの刻みや読み出しの時間構造が偽のピークになる）。
#   3. 周波数を全イベントで詰める: Z^2_n（n = n_harmonics）を、使う時間を 5 倍ずつ伸ばしながら最大にする。
#      Z^2_n が n_sigma に届かなければ（畳んだ分布が平ら）周期は無かったことにする。
#   4. その周期で畳んだ位相の分布（profile）。最大の位相を phase_peak、最大値の on_level 倍を超える
#      連続した範囲を on-beam の窓にする。位相と時刻の 2 次元分布で、周期のずれ（カウンタとビームの周波数の差）を見る。
#
# 出力: <out>_<ラン>.png、<out>_<ラン>_events.csv（基準のボードのイベントごとの時刻・位相・on-beam。
#       HA どうし / TPC どうしは trig で他のボードと組にできる）、<out>.json（ボードごとの周期と位相、窓）。
#
# 入力と条件は既定で ETCC_analysis/param/beam_period.json から読む。コマンドラインで渡したものが優先する。
#
# 例:
#   python3 ETCC_analysis/beam_period.py
#   python3 ETCC_analysis/beam_period.py --runs 20261201_run05 --sources ha --max-files 2
#   python3 ETCC_analysis/beam_period.py -i <calltt のラン> <killua のラン> --frequency 1000
#

import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np
from natsort import natsorted
from numba import njit

ETCC_DIR = os.path.dirname(os.path.abspath(__file__))
# このリポジトリ (analysis_sm3_FM) の場所。JSON の中の相対パスはここから
REPO_DIR = os.path.dirname(ETCC_DIR)
for d in (os.path.join(REPO_DIR, "PSA_analysis", "src_onoda"), os.path.join(REPO_DIR, "TPC_analysis")):
    if d not in sys.path:
        sys.path.insert(0, d)

from psa_decode import find_events_jit, HEADER as HA_HEADER, FOOTER as HA_FOOTER  # noqa: E402
from check_tpc import (HEADER, FOOTER, HEAD_SIZE, FADC_SIZE, RECORD_SIZE, MIN_FRAME_SIZE,  # noqa: E402
                       N_CLOCKS, read_word)

DEFAULT_CONFIG = os.path.join(ETCC_DIR, "param", "beam_period.json")
HA_RE = re.compile(r"^ha(\d+)_(\d+)\.dat$", re.IGNORECASE)
TPC_RE = re.compile(r"^tpc_enc_([ac])(\d+)_(\d+)\.raw$")
RUN_RE = re.compile(r"^(\d{8})_(run\d+)$")
WRAP = 1 << 32

DEFAULTS = {
    "nas": "/mnt/nas_eureka/SMILE-3FM_data",
    "storages": ["calltt", "killua"],
    "runs": [],
    "out": "ETCC_analysis/output/beam_period",
    "sources": ["ha", "tpc"],
    "ha_tick": 1e-6,
    "tpc_tick": 1e-4,
    "spectra": [{"bin": 1e-6, "segment": 1.0}, {"bin": 1e-3, "segment": 100.0}],
    "max_segments": 300,
    "f_min": 0.5,
    "f_max": 1e5,
    "frequency": None,
    "n_sigma": 5.0,
    "max_subharmonic": 8,
    "min_period_ticks": 10,
    "n_harmonics": 4,
    "phase_bins": 100,
    "on_level": 0.5,
    "rate_bin": 1.0,
    "max_files": 0,
}


def parse_args():

    p = argparse.ArgumentParser(description="イベントの時刻の周期を調べ、ビームのパルスの周期と位相を出す（raw を直接読む）")
    p.add_argument("--config", default=None,
                   help=f"入力と条件の JSON（既定 {os.path.relpath(DEFAULT_CONFIG, REPO_DIR)}。無ければ使わない）")
    p.add_argument("-i", "--input", nargs="+", default=None, help="1 つのランの raw のディレクトリ（calltt と killua の両方）")
    p.add_argument("--runs", nargs="+", default=None, help="<日付>_runNN。NAS の calltt と killua を読む。ランごとに解析する")
    p.add_argument("-o", "--out", default=None, help="出力の名前（<out>_<ラン>.png、<out>_<ラン>_events.csv、<out>.json）")
    p.add_argument("--sources", nargs="+", choices=["ha", "tpc"], default=None, help="読むもの（既定 ha tpc）")
    p.add_argument("--frequency", type=float, default=None, help="ビームの周波数 [Hz]（与えると探さずに、そのまわりで詰める）")
    for k, t, h in (("f_min", float, "探す周波数の下限 [Hz]"), ("f_max", float, "探す周波数の上限 [Hz]"),
                    ("n_harmonics", int, "Z^2 の高調波の数"), ("phase_bins", int, "位相のビンの数"),
                    ("on_level", float, "on-beam の窓: profile の最大値のこの倍を超える範囲"),
                    ("max_files", int, "ボードごとに読むファイルの数の上限（0 = 全部）")):
        p.add_argument("--" + k.replace("_", "-"), type=t, default=None, help=h)

    return p.parse_args()


def in_repo(path):
    return path if os.path.isabs(path) else os.path.join(REPO_DIR, path)


def load_settings(args):
    """DEFAULTS < JSON < コマンドライン"""

    s = dict(DEFAULTS)
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
            if not k.startswith("_"):
                s[k] = in_repo(v) if k == "out" else v
        print(f"config: {config}")
    s["config"] = config
    for k in ("out", "sources", "frequency", "f_min", "f_max", "n_harmonics", "phase_bins", "on_level", "max_files"):
        v = getattr(args, k)
        if v is not None:
            s[k] = in_repo(v) if k == "out" else v
    if args.input is not None:
        s["runs"] = [{"name": "input", "dirs": args.input}]
    elif args.runs is not None:
        s["runs"] = list(args.runs)
    if not s["runs"]:
        raise SystemExit("ランが無い: ETCC_analysis/param/beam_period.json の runs、--runs か -i で渡す")

    return s


def resolve_run(spec, s):
    """"20261123_run03"、{"date": ..., "run": ...}、{"dirs": [...]} → (名前, [ディレクトリ])"""

    if isinstance(spec, str):
        m = RUN_RE.match(spec)
        if not m:
            raise SystemExit(f"ランは <日付>_runNN の形で書く: {spec}")
        spec = {"date": m.group(1), "run": m.group(2)}
    if "dirs" in spec:
        dirs = [in_repo(d) for d in spec["dirs"]]
        name = spec.get("name") or os.path.basename(os.path.normpath(dirs[0]))
    else:
        date, run = str(spec["date"]), str(spec["run"])
        dirs = [os.path.join(s["nas"], date, st, f"{date}_{run}") for st in s["storages"]]
        name = spec.get("name", f"{date}_{run}")

    return name, dirs


def collect_boards(dirs, sources, max_files):
    """{"ha03": [file, ...], "a00": [...], ...}（番号順。同じ名前が 2 か所にあれば最初だけ）"""

    found = {}
    for d in dirs:
        if not os.path.isdir(d):
            print(f"[WARN] 無い: {d}")
            continue
        for path in glob.glob(os.path.join(d, "**", "*"), recursive=True):
            name = os.path.basename(path)
            m = HA_RE.match(name)
            if m and "ha" in sources:
                key = (f"ha{int(m.group(1)):02d}", int(m.group(2)))
            else:
                m = TPC_RE.match(name)
                if not (m and "tpc" in sources):
                    continue
                key = (f"{m.group(1)}{int(m.group(2)):02d}", int(m.group(3)))
            if key in found:
                print(f"[WARN] 同じファイルが 2 か所にある。最初の方を使う: {found[key]} / {path}")
                continue
            if os.path.getsize(path) == 0:
                continue
            found[key] = path
    boards = {}
    for (tag, num), path in sorted(found.items()):
        boards.setdefault(tag, []).append(path)
    return {k: (natsorted(v)[:max_files] if max_files > 0 else natsorted(v)) for k, v in boards.items()}


@njit(cache=True, nogil=True)
def scan_tpc(buf, trig, tim):
    """TPC のフレームの N_trg と N_time（check_tpc.decode_file と同じ読み方）"""

    n_bytes = len(buf)
    pos = 0
    n = 0
    while pos + MIN_FRAME_SIZE <= n_bytes:
        if read_word(buf, pos) != HEADER:
            pos += 1
            continue
        t_trg = read_word(buf, pos + 4)
        t_time = read_word(buf, pos + 8)
        body = pos + HEAD_SIZE
        if read_word(buf, body) == FOOTER:
            trig[n] = t_trg
            tim[n] = t_time
            n += 1
            pos = body + 4
            continue
        record = body + FADC_SIZE + 4
        count = 0
        found = False
        while count < N_CLOCKS and record + 4 <= n_bytes:
            if read_word(buf, record) == FOOTER:
                found = True
                break
            if record + RECORD_SIZE > n_bytes:
                break
            record += RECORD_SIZE
            count += 1
        if not found and count < N_CLOCKS:
            pos += 4
            continue
        trig[n] = t_trg
        tim[n] = t_time
        n += 1
        pos = record + 4 if found else record
    return n


def read_counter(path):
    """1 ファイルの (trig, カウンタ)"""

    buf = np.fromfile(path, dtype=np.uint8)
    if os.path.basename(path).lower().startswith("ha"):
        starts, ends = find_events_jit(buf, *HA_HEADER, *HA_FOOTER)
        st = np.asarray(starts, dtype=np.int64)
        st = st[st + 12 <= len(buf)]
        word = lambda o: ((buf[st + o].astype(np.int64) << 24) | (buf[st + o + 1].astype(np.int64) << 16)
                          | (buf[st + o + 2].astype(np.int64) << 8) | buf[st + o + 3].astype(np.int64))
        return word(4), word(8)
    cap = len(buf) // MIN_FRAME_SIZE + 1
    trig = np.zeros(cap, np.uint32)
    tim = np.zeros(cap, np.uint32)
    n = scan_tpc(buf, trig, tim)
    return trig[:n].astype(np.int64), tim[:n].astype(np.int64)


def read_board(paths):
    """1 ボードの全ファイル: (epoch, trig, つながったカウンタ)。カウンタが 2^32 で戻ったら足し直す。
    同じ (trig, カウンタ) は 1 回だけ"""

    parts = [read_counter(p) for p in paths]
    trig = np.concatenate([a for a, _ in parts]) if parts else np.zeros(0, np.int64)
    cnt = np.concatenate([b for _, b in parts]) if parts else np.zeros(0, np.int64)
    if len(cnt) == 0:
        return trig, trig, cnt, 0
    d = np.diff(cnt) % WRAP
    d = np.where(d >= WRAP // 2, d - WRAP, d)      # 戻ったもの（順番の入れ替わり）は負のまま
    n_back = int(np.sum(d < 0))
    ticks = cnt[0] + np.r_[0, np.cumsum(d)]
    back = np.r_[False, trig[1:] < trig[:-1]]
    epoch = np.cumsum(back)
    key = np.unique(np.column_stack([epoch, trig, ticks]), axis=0, return_index=True)[1]
    key = np.sort(key)
    epoch, trig, ticks = epoch[key], trig[key], ticks[key]
    order = np.argsort(ticks, kind="stable")
    return epoch[order], trig[order], ticks[order], n_back


def power_spectrum(t, bin_s, segment_s, max_segments):
    """segment ごとの FFT の平均（Leahy: 2|F|^2 / N。雑音なら平均 2）。(周波数, パワー, segment の数)"""

    n_bins = int(round(segment_s / bin_s))
    if n_bins < 4:
        return None
    t0, span = t[0], t[-1] - t[0]
    n_seg = min(int(span // segment_s), max_segments)
    if n_seg < 1:
        return None
    edges = np.searchsorted(t, t0 + segment_s * np.arange(n_seg + 1))
    acc = np.zeros(n_bins // 2 + 1)
    used = 0
    for k in range(n_seg):
        tt = t[edges[k]:edges[k + 1]]
        if len(tt) < 2:
            continue
        idx = np.minimum(((tt - (t0 + k * segment_s)) / bin_s).astype(np.int64), n_bins - 1)
        c = np.bincount(idx, minlength=n_bins).astype(np.float64)
        f = np.fft.rfft(c)
        acc += 2 * np.abs(f) ** 2 / len(tt)
        used += 1
    if used == 0:
        return None
    freq = np.fft.rfftfreq(n_bins, bin_s)
    return freq[1:], acc[1:] / used, used


def find_fundamental(spectra, s):
    """いちばん有意なピークと、その 1/k に有意なピークがあればそちら"""

    best = None
    for sp in spectra:
        if sp is None:
            continue
        f, p, m = sp["freq"], sp["power"], sp["segments"]
        thr = 2 + s["n_sigma"] * 2 / np.sqrt(m)
        sel = (f >= s["f_min"]) & (f <= s["f_max"])
        if not sel.any():
            continue
        i = np.flatnonzero(sel)[np.argmax(p[sel])]
        signif = (p[i] - 2) / (2 / np.sqrt(m))
        if p[i] < thr or (best is not None and signif <= best["sigma"]):
            continue
        fp, df = f[i], f[1] - f[0]
        f0 = fp
        for k in range(s["max_subharmonic"], 1, -1):
            fk = fp / k
            if fk < max(s["f_min"], f[0]):
                continue
            j = np.flatnonzero(np.abs(f - fk) <= 1.5 * df)
            if len(j) and p[j].max() > thr:
                f0 = f[j[np.argmax(p[j])]]
                break
        best = {"peak_hz": float(fp), "fundamental_hz": float(f0), "sigma": float(signif), "resolution_hz": float(df),
                "spectrum": sp["label"]}
    return best


@njit(cache=True)
def z2_grid(t, freqs, n_harm):
    out = np.zeros(len(freqs))
    n = len(t)
    for i in range(len(freqs)):
        tot = 0.0
        for k in range(1, n_harm + 1):
            c = 0.0
            sn = 0.0
            w = 2 * np.pi * k * freqs[i]
            for j in range(n):
                ph = w * t[j]
                c += np.cos(ph)
                sn += np.sin(ph)
            tot += c * c + sn * sn
        out[i] = 2.0 * tot / max(n, 1)
    return out


def refine(t, f0, df0, n_harm, max_events=2_000_000):
    """Z^2_n を最大にする周波数。使う時間 L を 5 倍ずつ伸ばし、そのたびに f ± 1/L_前 を 121 点で探す"""

    t = t - t[0]
    span = t[-1]
    f, width, L_prev = f0, 1.5 * df0, 1 / df0
    while True:
        L = min(span, 5 * L_prev)
        u = t[t <= L]
        if len(u) > max_events:
            u = u[np.linspace(0, len(u) - 1, max_events).astype(np.int64)]
        grid = f + np.linspace(-width, width, 121)
        grid = grid[grid > 0]
        z = z2_grid(u, grid, n_harm)
        f = float(grid[np.argmax(z)])
        z_best = float(z.max())
        if L >= span:
            break
        width, L_prev = 1 / L, L
    return f, z_best, float(1 / max(span, 1e-12))


def profile(t, f, s):
    """位相の分布と on-beam の窓（最大値の on_level 倍を超える、最大の位相を含む連続した範囲）"""

    nb = s["phase_bins"]
    ph = ((t - t[0]) * f) % 1.0
    h = np.bincount(np.minimum((ph * nb).astype(np.int64), nb - 1), minlength=nb).astype(float)
    ipk = int(np.argmax(h))
    base = np.median(np.sort(h)[: max(nb // 4, 1)])
    level = base + s["on_level"] * (h[ipk] - base)
    # 最大のビンから左右に、level を超えるビンが続くところまで（位相は輪なので回り込む）
    lo, n_on = ipk, 1
    while n_on < nb and h[(lo - 1) % nb] > level:
        lo, n_on = lo - 1, n_on + 1
    while n_on < nb and h[(lo + n_on) % nb] > level:
        n_on += 1
    lo %= nb
    width = n_on / nb
    win = (lo / nb, ((lo + n_on) % nb) / nb)
    on = ((ph - win[0]) % 1.0) < width
    frac_on = float(on.mean())
    pulsed = float(1 - base * nb / max(h.sum(), 1))
    return dict(phase=ph, hist=h, peak_phase=(ipk + 0.5) / nb, window=win, width=width, on=on,
                on_fraction=frac_on, pulsed_fraction=pulsed, baseline=float(base))


def analyse_board(tag, ticks, tick, s):
    t = (ticks - ticks[0]) * tick
    span = float(t[-1]) if len(t) else 0.0
    res = {"board": tag, "events": int(len(t)), "span_s": span, "tick_s": tick}
    if len(t) < 10 or span <= 0:
        res["note"] = "too few events"
        return res, None
    res["mean_rate_hz"] = len(t) / span
    spectra = []
    for sp in s["spectra"]:
        b = max(sp["bin"], tick)
        r = power_spectrum(t, b, sp["segment"], s["max_segments"])
        if r:
            spectra.append({"freq": r[0], "power": r[1], "segments": r[2],
                            "label": f"bin {b:g} s, segment {sp['segment']:g} s"})
    # 周期が数カウントしかないと位相が決まらない（カウンタの刻みや読み出しの時間構造が偽のピークになる）
    s = dict(s, f_max=min(s["f_max"], 1 / (s["min_period_ticks"] * tick)))
    if s["frequency"]:
        f0 = float(s["frequency"])
        df0 = min([sp["freq"][1] - sp["freq"][0] for sp in spectra], default=1 / span)
        res["search"] = {"given_hz": f0}
    else:
        cand = find_fundamental(spectra, s)
        res["search"] = cand
        if cand is None:
            res["note"] = "no significant period"
            return res, dict(t=t, spectra=spectra)
        f0, df0 = cand["fundamental_hz"], cand["resolution_hz"]
    f, z2, df = refine(t, f0, df0, s["n_harmonics"])
    # 畳んだ分布が平らでないか: Z^2_n は周期が無ければ自由度 2n の χ^2
    from scipy.stats import chi2, norm
    z2_sigma = float(norm.isf(max(chi2.sf(z2, 2 * s["n_harmonics"]), 1e-300)))
    if not s["frequency"] and z2_sigma < s["n_sigma"]:
        res["search"]["refined_hz"] = f
        res["search"]["z2_sigma"] = z2_sigma
        res["note"] = "peak not confirmed by folding"
        return res, dict(t=t, spectra=spectra)
    pr = profile(t, f, s)
    res.update({"frequency_hz": f, "period_s": 1 / f, "period_ticks": 1 / f / tick, "frequency_resolution_hz": df,
                f"z2_{s['n_harmonics']}": z2, "z2_sigma": z2_sigma,
                "phase_origin_ticks": int(ticks[0]),
                "peak_phase": pr["peak_phase"], "on_window_phase": list(pr["window"]), "on_width": pr["width"],
                "on_fraction": pr["on_fraction"], "pulsed_fraction": pr["pulsed_fraction"]})
    return res, dict(t=t, spectra=spectra, f=f, prof=pr)


def plot_run(name, out, per_board, ref, s):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res, d = per_board[ref]
    t = d["t"]
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))
    ax = axs[0, 0]
    nb = max(int(np.ceil(t[-1] / s["rate_bin"])), 1)
    h, e = np.histogram(t, bins=nb, range=(0, nb * s["rate_bin"]))
    ax.step(e[:-1], h / s["rate_bin"], where="post")
    ax.set(xlabel="time [s]", ylabel="rate [Hz]", title=f"{ref}: rate ({s['rate_bin']:g} s bins)")

    ax = axs[0, 1]
    dt = np.diff(t)
    dt = dt[dt > 0]
    if len(dt):
        lo = max(dt.min(), res["tick_s"])
        ax.hist(dt, bins=np.logspace(np.log10(lo), np.log10(dt.max() * 1.01), 200), histtype="step")
        ax.set_xscale("log")
        ax.set_yscale("log")
    if "frequency_hz" in res:
        for k in range(1, 6):
            ax.axvline(k / res["frequency_hz"], color="m", ls=":", lw=0.8)
    ax.set(xlabel="Δt [s]", ylabel="events", title="time between events (magenta: n × period)")

    for ax, sp in zip((axs[0, 2], axs[1, 0]), (d["spectra"] + [None, None])[:2]):
        if sp is None:
            ax.axis("off")
            continue
        ax.plot(sp["freq"], sp["power"], lw=0.6)
        ax.axhline(2 + s["n_sigma"] * 2 / np.sqrt(sp["segments"]), color="r", ls="--", lw=0.8,
                   label=f"{s['n_sigma']:g} σ")
        if "frequency_hz" in res:
            ax.axvline(res["frequency_hz"], color="m", ls=":", lw=0.8, label="f0")
        ax.set_xscale("log")
        ax.set(xlabel="frequency [Hz]", ylabel="Leahy power", title=f"power spectrum ({sp['label']}, {sp['segments']} seg.)")
        ax.legend(fontsize=8)

    ax = axs[1, 1]
    if "prof" in d:
        pr = d["prof"]
        nbp = len(pr["hist"])
        x = (np.arange(2 * nbp) + 0.5) / nbp
        ax.step(x, np.tile(pr["hist"], 2), where="mid")
        lo, w = pr["window"][0], pr["width"]
        for c in (0, 1, 2):
            ax.axvspan(lo + c - 1, lo + c - 1 + w, color="orange", alpha=0.25)
        ax.set_xlim(0, 2)
        ax.set(xlabel="phase", ylabel="events",
               title=f"folded at {res['frequency_hz']:.6g} Hz (period {res['period_ticks']:.3f} ticks)\n"
                     f"on window (orange): {pr['window'][0]:.3f}-{pr['window'][1]:.3f}, "
                     f"{100 * pr['on_fraction']:.1f} % of events, pulsed {100 * pr['pulsed_fraction']:.1f} %")
        ax = axs[1, 2]
        ax.hist2d(t, pr["phase"], bins=[min(200, max(nb, 10)), 50], cmap="turbo")
        ax.set(xlabel="time [s]", ylabel="phase", title="phase vs time (a slope = period mismatch)")
    else:
        axs[1, 1].text(0.5, 0.5, res.get("note", ""), ha="center")
        axs[1, 2].axis("off")

    found = [f"{b}: {r['frequency_hz']:.6g} Hz, peak phase {r['peak_phase']:.3f}, on {100 * r['on_fraction']:.0f} %"
             for b, (r, _) in per_board.items() if "frequency_hz" in r]
    n_none = len(per_board) - len(found)
    lines = "; ".join(found[:6]) + (f"; ... ({len(found) - 6} more)" if len(found) > 6 else "")
    lines += (f"   [{n_none} boards without a period: see the JSON]" if n_none else "")
    fig.suptitle(f"{os.path.basename(out)}: {name} (reference board {ref}; every board has its own counter, so the phases differ)\n"
                 + lines, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, .93))
    fig.savefig(f"{out}_{name}.png", dpi=100)
    plt.close(fig)


def main():

    args = parse_args()
    s = load_settings(args)
    out = s["out"]
    if os.path.dirname(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)
    summary = {"config": s["config"], "settings": {k: v for k, v in s.items() if k != "config"}, "runs": []}

    for spec in s["runs"]:
        name, dirs = resolve_run(spec, s)
        boards = collect_boards(dirs, s["sources"], s["max_files"])
        if not boards:
            print(f"[WARN] [{name}] raw が無い: {dirs}")
            continue
        per_board, raw = {}, {}
        for tag, paths in boards.items():
            t0 = time.time()
            epoch, trig, ticks, n_back = read_board(paths)
            tick = s["ha_tick"] if tag.startswith("ha") else s["tpc_tick"]
            res, d = analyse_board(tag, ticks, tick, s)
            res["files"] = len(paths)
            res["counter_backsteps"] = n_back
            per_board[tag] = (res, d)
            raw[tag] = (epoch, trig, ticks)
            msg = (f"{res['frequency_hz']:.6g} Hz, peak phase {res['peak_phase']:.3f}, on {100 * res['on_fraction']:.1f} %"
                   if "frequency_hz" in res else res.get("note", ""))
            print(f"[{name}] {tag}: {len(paths)} files, {res['events']} events, {res['span_s']:.1f} s: {msg} "
                  f"({time.time() - t0:.1f} s)", flush=True)
            if n_back:
                print(f"[WARN] {tag}: カウンタが {n_back} 回戻った（順番の入れ替わり。並べ直して使う）")
        ok = [b for b, (r, d) in per_board.items() if d is not None]
        if not ok:
            continue
        # 基準のボード: 周期が出たものの中で、HA を先に、イベントの多いもの
        ref = max(ok, key=lambda b: ("frequency_hz" in per_board[b][0], b.startswith("ha"), per_board[b][0]["events"]))
        plot_run(name, out, per_board, ref, s)
        res, d = per_board[ref]
        if "prof" in d:
            epoch, trig, ticks = raw[ref]
            with open(f"{out}_{name}_events.csv", "w") as f:
                f.write(f"# board {ref}; frequency {res['frequency_hz']:.9g} Hz; phase = ((ticks - {res['phase_origin_ticks']}) "
                        f"* tick * f) mod 1; on window {res['on_window_phase'][0]:.4f}-{res['on_window_phase'][1]:.4f}\n")
                f.write("board,epoch,trig,ticks,t_s,phase,on_beam\n")
                for e, tr, tk, ts, ph, on in zip(epoch, trig, ticks, d["t"], d["prof"]["phase"], d["prof"]["on"]):
                    f.write(f"{ref},{e},{tr},{tk},{ts:.6f},{ph:.5f},{int(on)}\n")
        summary["runs"].append({"name": name, "dirs": dirs, "reference_board": ref,
                                "boards": {b: r for b, (r, _) in per_board.items()}})
        print(f"saved {out}_{name}.png" + (f", {out}_{name}_events.csv" if "prof" in d else ""))

    with open(f"{out}.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(f"saved {out}.json")


if __name__ == "__main__":
    main()
