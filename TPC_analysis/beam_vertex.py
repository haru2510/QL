#!/usr/bin/env python3
#
# ビーム試験の簡易解析: TPC の散乱点（対生成の頂点）の 3 次元の分布から、γ線ビームの通った道を出す。
#
# raw（tpc_enc_a<NN>_<番号>.raw、tpc_enc_c<NN>_<番号>.raw）を直接読む。フレームの読み方は check_tpc.py と同じ。
# 1 つのランのボードは 2 つの storage（calltt: アノード、killua: カソード）に分かれるので両方読み、
# ラン 1 つの中でトリガーカウンタ（N_trg）でそろえる（ボードの中で N_trg が戻ったら epoch を 1 つ進める）。
# 複数のランは、ランごとに頂点を出してから足す。
#
# 頂点: イベントの中で、同じ uPIC のアノードとカソードの両方にヒットがある一番早いクロック（z が一番小さい点）。
# x と y はそのクロックのアノード・カソードのストリップの位置の平均（uPIC ごとの向きは check_tpc.py の STRIP_AXES、
# EricaNext の TOT_Skewness::Get_TPC_XYZ と同じ）。ToT 補正と 3 次元の組み合わせ（偽の点が出る）はしない。
# 対生成では e+ e- が頂点から前に飛ぶので、ビームが上から下（+Z）なら、一番上流の点 = z が一番小さい点。
# z = velocity × (clock − clock_offset) + z_gem（EricaNext と同じ。既定値は EricaNext の calib/tpc/drift_velocity.dat）。
#
# fiducial: 頂点が壁（±153.6 mm）とドリフト領域の両端（z_gem、317 mm）から margin だけ内側。
# 外で対生成して入ってきた e± の事象（頂点が壁や端に張り付く）を落とす。
# ノイズ: ボードの最初のファイルで、フレームの hot_fraction 以上で鳴るストリップは使わない。
# ヒットが min_hits より少ない事象は使わない。メモリを抑えるため、フレームごとにクロックでまとめ（数、位置の和）、
# 早いクロックから keep_clocks 個だけ持つ（頂点には一番早いところしか要らない）。
#
# ビームの軸: fiducial を通った頂点に 3 次元の直線を fit する（主成分。軸から遠い点を外して fit し直す）。
#
# 入力と条件は既定で TPC_analysis/param/beam_vertex.json から読む。コマンドラインで渡したものが優先する。
#
# 例:
#   python3 TPC_analysis/beam_vertex.py                                   # JSON のとおり
#   python3 TPC_analysis/beam_vertex.py --runs 20261123_run03 20261123_run04 --max-files 2
#   python3 TPC_analysis/beam_vertex.py -i <calltt のラン> <killua のラン> -o TPC_analysis/beam_output/run03
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

TPC_DIR = os.path.dirname(os.path.abspath(__file__))
if TPC_DIR not in sys.path:
    sys.path.insert(0, TPC_DIR)

from check_tpc import (HEADER, FOOTER, HEAD_SIZE, FADC_SIZE, RECORD_SIZE, MIN_FRAME_SIZE,  # noqa: E402
                       N_CLOCKS, N_FADC_CH, N_STRIPS, N_BOARDS_PER_UPIC, STRIP_AXES, read_word)

# このリポジトリ (analysis_sm3_FM) の場所。JSON の中の相対パスはここから
REPO_DIR = os.path.dirname(TPC_DIR)
DEFAULT_CONFIG = os.path.join(TPC_DIR, "param", "beam_vertex.json")
FILE_RE = re.compile(r"^tpc_enc_([ac])(\d+)_(\d+)\.raw$")
RUN_RE = re.compile(r"^(\d{8})_(run\d+)$")
HALF_XY = 153.6        # ドリフト領域の半分の幅 [mm]（384 ストリップ × 0.4 mm）
Z_TOP = 317.0          # μ-PIC から天板まで [mm]

DEFAULTS = {
    "nas": "/mnt/nas_eureka/SMILE-3FM_data",
    "storages": ["calltt", "killua"],
    "runs": [],
    "out": "TPC_analysis/beam_output/beam_vertex",
    "pitch_strip": 0.4,
    "clock_offset": 133.69,
    "velocity": 0.39358,
    "z_gem": 5.0,
    "margin": 10.0,
    "min_hits": 10,
    "hot_fraction": 0.05,
    "fit_cut_sigma": 3.0,
    "fit_iterations": 5,
    "keep_clocks": 16,
    "max_files": 0,
}


def parse_args():

    p = argparse.ArgumentParser(description="TPC の散乱点の 3 次元の分布からビームの道を出す簡易解析（raw を直接読む）")
    p.add_argument("--config", default=None,
                   help=f"入力と条件の JSON（既定 {os.path.relpath(DEFAULT_CONFIG, REPO_DIR)}。無ければ使わない）")
    p.add_argument("-i", "--input", nargs="+", default=None,
                   help="1 つのランの raw のディレクトリ（calltt と killua の両方）。JSON のランの代わりに使う")
    p.add_argument("--runs", nargs="+", default=None,
                   help="<日付>_runNN（例 20261123_run03）。NAS の calltt と killua を読む。複数なら足す")
    p.add_argument("-o", "--out", default=None, help="出力の名前（<out>.png、<out>_vertices.csv、<out>.json）")
    for k, t, h in (("clock_offset", float, "GEM 面のクロック [clock]"), ("velocity", float, "ドリフト速度 [mm/clock]"),
                    ("z_gem", float, "GEM 面の z [mm]"), ("margin", float, "fiducial の余白 [mm]"),
                    ("min_hits", int, "事象の点の数の下限"), ("hot_fraction", float, "これ以上の割合の事象で鳴るストリップを除く"),
                    ("max_files", int, "ボードごとに読むファイルの数の上限（0 = 全部）"),
                    ("keep_clocks", int, "フレームごとに残す早いクロックの数")):
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
    for k in ("out", "clock_offset", "velocity", "z_gem", "margin", "min_hits", "hot_fraction", "max_files", "keep_clocks"):
        v = getattr(args, k)
        if v is not None:
            s[k] = v
    if args.input is not None:
        s["runs"] = [{"name": "input", "dirs": args.input}]
    elif args.runs is not None:
        s["runs"] = list(args.runs)
    if not s["runs"]:
        raise SystemExit("ランが無い: TPC_analysis/param/beam_vertex.json の runs、--runs か -i で渡す")

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


def collect_boards(dirs, max_files):
    """{(side, board): [file, ...]}（番号順。同じ名前が 2 か所にあれば最初だけ）"""

    found = {}
    for d in dirs:
        if not os.path.isdir(d):
            print(f"[WARN] 無い: {d}")
        for path in glob.glob(os.path.join(d, "**", "tpc_enc_*.raw"), recursive=True):
            m = FILE_RE.match(os.path.basename(path))
            if not m:
                continue
            key = (m.group(1), int(m.group(2)), int(m.group(3)))
            if key in found:
                print(f"[WARN] 同じファイルが 2 か所にある。最初の方を使う: {found[key]} / {path}")
                continue
            if os.path.getsize(path) == 0:
                print(f"[WARN] 空のファイル: {path}")
                continue
            found[key] = path
    boards = {}
    for (side, board, num), path in sorted(found.items()):
        boards.setdefault((side, board), []).append(path)
    return {k: (natsorted(v)[:max_files] if max_files > 0 else natsorted(v)) for k, v in boards.items()}


@njit(cache=True, nogil=True)
def scan_frames(buf, trig, n_hits_out, offsets, n_records):
    """フレームの位置と N_trg、ヒットの数。check_tpc.decode_file と同じ読み方（FADC は読まない）"""

    n_bytes = len(buf)
    pos = 0
    n = 0
    while pos + MIN_FRAME_SIZE <= n_bytes:
        if read_word(buf, pos) != HEADER:
            pos += 1
            continue
        event_trig = read_word(buf, pos + 4)
        body = pos + HEAD_SIZE
        if read_word(buf, body) == FOOTER:
            trig[n] = event_trig
            n_records[n] = 0
            offsets[n] = pos
            n_hits_out[n] = 0
            n += 1
            pos = body + 4
            continue
        first = body + FADC_SIZE + 4
        record = first
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
        hits = 0
        for index in range(count):
            start = first + RECORD_SIZE * index
            if (read_word(buf, start) & 0x7ff) >= N_CLOCKS:
                continue
            for ch in range(N_FADC_CH):
                w = read_word(buf, start + 4 * (4 - ch))
                for bit in range(32):
                    if (w >> bit) & 1:
                        hits += 1
        trig[n] = event_trig
        n_records[n] = count
        offsets[n] = pos
        n_hits_out[n] = hits
        n += 1
        pos = record + 4 if found else record
    return n


@njit(cache=True, nogil=True)
def count_strips(buf, offsets, n_records, counts):
    """ストリップごとに、鳴ったフレームの数（鳴りっぱなしのストリップを探す）"""

    seen = np.zeros(N_STRIPS, np.bool_)
    for f in range(len(offsets)):
        seen[:] = False
        first = int(offsets[f]) + HEAD_SIZE + FADC_SIZE + 4
        for index in range(n_records[f]):
            start = first + RECORD_SIZE * index
            if (read_word(buf, start) & 0x7ff) >= N_CLOCKS:
                continue
            for ch in range(N_FADC_CH):
                w = read_word(buf, start + 4 * (4 - ch))
                for bit in range(32):
                    if (w >> bit) & 1:
                        seen[32 * ch + bit] = True
        for k in range(N_STRIPS):
            if seen[k]:
                counts[k] += 1


@njit(cache=True, nogil=True)
def fill_rows(buf, offsets, n_records, hot, keep_clocks, strip_pos,
              row_frame, row_clock, row_count, row_sum, frame_hits):
    """フレームごとに、鳴りっぱなしでないストリップのヒットをクロックごとにまとめ（数、位置の和）、
    早いクロックから keep_clocks 個だけ行にする。frame_hits はフレームのヒットの数。返り値は行の数"""

    n_rows = 0
    clk = np.zeros(N_CLOCKS, np.int64)
    cnt = np.zeros(N_CLOCKS, np.int64)
    sm = np.zeros(N_CLOCKS, np.float64)
    for f in range(len(offsets)):
        first = int(offsets[f]) + HEAD_SIZE + FADC_SIZE + 4
        m = 0
        total = 0
        for index in range(n_records[f]):
            start = first + RECORD_SIZE * index
            clock = read_word(buf, start) & 0x7ff
            if clock >= N_CLOCKS:
                continue
            c = 0
            ssum = 0.0
            for ch in range(N_FADC_CH):
                w = read_word(buf, start + 4 * (4 - ch))
                for bit in range(32):
                    if (w >> bit) & 1:
                        k = 32 * ch + bit
                        if not hot[k]:
                            c += 1
                            ssum += strip_pos[k]
            if c > 0 and m < N_CLOCKS:
                clk[m] = clock
                cnt[m] = c
                sm[m] = ssum
                m += 1
                total += c
        frame_hits[f] = total
        order = np.argsort(clk[:m])
        for q in range(min(m, keep_clocks)):
            i = order[q]
            row_frame[n_rows] = f
            row_clock[n_rows] = clk[i]
            row_count[n_rows] = cnt[i]
            row_sum[n_rows] = sm[i]
            n_rows += 1
    return n_rows


def scan_file(path):
    buf = np.fromfile(path, dtype=np.uint8)
    cap = len(buf) // MIN_FRAME_SIZE + 1
    trig = np.zeros(cap, np.uint32)
    nh = np.zeros(cap, np.int64)
    offsets = np.zeros(cap, np.int64)
    nrec = np.zeros(cap, np.int32)
    n = scan_frames(buf, trig, nh, offsets, nrec)
    return buf, trig[:n].astype(np.int64), offsets[:n], nrec[:n]


def find_hot(paths, fraction):
    """最初のファイルで、フレームの fraction 以上で鳴るストリップ"""

    buf, trig, offsets, nrec = scan_file(paths[0])
    counts = np.zeros(N_STRIPS, np.int64)
    count_strips(buf, offsets, nrec, counts)
    return counts >= fraction * max(len(trig), 1)


def read_board(paths, hot, strip_pos, keep_clocks):
    """1 ボードの全ファイル: 行（event, clock, 数, 位置の和）と、フレームごとの (event, ヒットの数)。
    event = epoch * 2^32 + N_trg（ボードの中で N_trg が戻ったら epoch を進める）"""

    rows = {k: [] for k in ("event", "clock", "count", "sum")}
    fev, fhits = [], []
    epoch, prev, n_frames, n_back = 0, None, 0, 0
    seen = set()
    for path in paths:
        buf, trig, offsets, nrec = scan_file(path)
        n = len(trig)
        if n == 0:
            continue
        back = np.r_[prev is not None and trig[0] < prev, trig[1:] < trig[:-1]]
        ep = epoch + np.cumsum(back)
        n_back += int(back.sum())
        epoch, prev = int(ep[-1]), int(trig[-1])
        key = ep * (1 << 32) + trig
        keep = np.array([k not in seen for k in key.tolist()], dtype=bool)
        seen.update(key[keep].tolist())
        n_frames += int(keep.sum())
        cap = int(np.minimum(nrec, keep_clocks).sum())
        rf = np.zeros(cap, np.int64)
        rc = np.zeros(cap, np.int64)
        rn = np.zeros(cap, np.int64)
        rs = np.zeros(cap, np.float64)
        fh = np.zeros(n, np.int64)
        k = fill_rows(buf, offsets, nrec, hot, keep_clocks, strip_pos, rf, rc, rn, rs, fh)
        rf, rc, rn, rs = rf[:k], rc[:k], rn[:k], rs[:k]
        good = keep[rf]
        rows["event"].append(key[rf[good]])
        rows["clock"].append(rc[good])
        rows["count"].append(rn[good])
        rows["sum"].append(rs[good])
        fev.append(key[keep])
        fhits.append(fh[keep])
    cat = lambda a, d: np.concatenate(a) if a else np.zeros(0, d)
    return ({k: cat(v, np.float64 if k == "sum" else np.int64) for k, v in rows.items()},
            cat(fev, np.int64), cat(fhits, np.int64), n_frames, n_back)


def run_vertices(name, dirs, s):
    """ラン 1 つ: 全ボードを読み、事象ごとの頂点（z が一番小さい、同じ uPIC のアノードとカソードが同じクロックにある点）"""

    boards = collect_boards(dirs, s["max_files"])
    if not boards:
        raise SystemExit(f"[{name}] raw（tpc_enc_*.raw）が無い: {dirs}")
    print(f"[{name}] boards: " + ", ".join(f"{sd}{b:02d}({len(f)})" for (sd, b), f in boards.items()))
    cols = {k: [] for k in ("event", "upic", "side", "axis", "clock", "count", "sum")}
    ev_parts, hit_parts = [], []
    frames, resets, hot_list = {}, {}, {}
    for (side, board), paths in boards.items():
        t0 = time.time()
        upic, offset = board // N_BOARDS_PER_UPIC, board % N_BOARDS_PER_UPIC
        axis_name, sign, const = STRIP_AXES[upic][0 if side == "a" else 1]
        strip_pos = s["pitch_strip"] * (sign * (np.arange(N_STRIPS) + offset * N_STRIPS) + const)
        hot = find_hot(paths, s["hot_fraction"])
        rows, fev, fh, nf, nb = read_board(paths, hot, strip_pos, s["keep_clocks"])
        tag = f"{side}{board:02d}"
        frames[tag], resets[tag], hot_list[tag] = nf, nb, np.flatnonzero(hot).tolist()
        if nb:
            print(f"[WARN] {tag}: N_trg が {nb} 回戻った（epoch で分けて組にする）")
        n = len(rows["event"])
        cols["event"].append(rows["event"])
        cols["upic"].append(np.full(n, upic, np.int64))
        cols["side"].append(np.full(n, 0 if side == "a" else 1, np.int64))
        cols["axis"].append(np.full(n, 0 if axis_name == "x" else 1, np.int64))
        cols["clock"].append(rows["clock"])
        cols["count"].append(rows["count"])
        cols["sum"].append(rows["sum"])
        ev_parts.append(fev)
        hit_parts.append(fh)
        print(f"  {tag}: {len(paths)} files, {nf} frames, hot strips {int(hot.sum())} ({time.time() - t0:.1f} s)", flush=True)
    h = {k: np.concatenate(v) for k, v in cols.items()}

    # 事象のヒットの数（全ボードの和）
    fev = np.concatenate(ev_parts)
    fh = np.concatenate(hit_parts)
    ue, inv = np.unique(fev, return_inverse=True)
    n_hits_ev = np.bincount(inv.ravel(), weights=fh).astype(np.int64)
    n_events = int(np.sum(n_hits_ev > 0))
    good_ev = ue[n_hits_ev >= s["min_hits"]]
    keep = np.isin(h["event"], good_ev)
    h = {k: v[keep] for k, v in h.items()}

    # (event, upic, clock) ごとに、アノードとカソードの両方がある一番早いもの
    key3 = np.column_stack([h["event"], h["upic"], h["clock"]])
    u3, inv3 = np.unique(key3, axis=0, return_inverse=True)
    inv3 = inv3.ravel()
    has = np.zeros((len(u3), 2), bool)
    has[inv3, h["side"]] = True
    sums = np.zeros((len(u3), 2))
    nums = np.zeros((len(u3), 2))
    np.add.at(sums, (inv3, h["axis"]), h["sum"])
    np.add.at(nums, (inv3, h["axis"]), h["count"])
    both = np.flatnonzero(has.all(1))
    cand = u3[both]
    order = np.lexsort((cand[:, 2], cand[:, 0]))
    both, cand = both[order], cand[order]
    first = np.r_[True, cand[1:, 0] != cand[:-1, 0]] if len(cand) else np.zeros(0, bool)
    vi, vtx = both[first], cand[first]
    xy = sums[vi] / np.maximum(nums[vi], 1)
    z = s["velocity"] * (vtx[:, 2] - s["clock_offset"]) + s["z_gem"]
    nh = dict(zip(ue.tolist(), n_hits_ev.tolist()))
    print(f"[{name}] events with hits {n_events}, with >= {s['min_hits']} hits {len(good_ev)}, with a vertex {len(vtx)}")

    return dict(name=name, dirs=dirs, frames=frames, resets=resets, hot_strips=hot_list, events=n_events,
                vtx=np.column_stack([xy[:, 0], xy[:, 1], z]), vtx_event=vtx[:, 0], vtx_upic=vtx[:, 1], vtx_clock=vtx[:, 2],
                vtx_nhits=np.array([nh[e] for e in vtx[:, 0].tolist()], dtype=np.int64))


def fiducial(v, s):
    m = s["margin"]
    return ((np.abs(v[:, 0]) < HALF_XY - m) & (np.abs(v[:, 1]) < HALF_XY - m)
            & (v[:, 2] > s["z_gem"] + m) & (v[:, 2] < Z_TOP - m))


def fit_line(p, cut_sigma, iterations):
    """3 次元の直線（主成分）。軸からの距離が cut_sigma × RMS を超える点を外して fit し直す"""

    use = np.ones(len(p), bool)
    for _ in range(iterations):
        if use.sum() < 3:
            return None
        c = p[use].mean(0)
        w, v = np.linalg.eigh(np.cov((p[use] - c).T))
        d = v[:, -1]
        if d[2] < 0:
            d = -d
        r = np.linalg.norm(np.cross(p - c, d), axis=1)
        rms = np.sqrt(np.mean(r[use] ** 2))
        new = r < cut_sigma * rms
        if np.array_equal(new, use):
            break
        use = new
    theta = np.degrees(np.arccos(abs(d[2])))
    phi = np.degrees(np.arctan2(d[1], d[0]))
    return {"point_mm": c.tolist(), "direction": d.tolist(), "angle_from_z_deg": float(theta), "azimuth_deg": float(phi),
            "used": int(use.sum()), "distance_rms_mm": float(rms),
            "distance_r50_r68_r90_mm": [float(x) for x in np.percentile(r[use], [50, 68, 90])],
            "use": use, "r": r}


def main():

    args = parse_args()
    s = load_settings(args)
    runs = [run_vertices(*resolve_run(r, s), s) for r in s["runs"]]
    out = s["out"]
    if os.path.dirname(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)

    V = np.vstack([r["vtx"] for r in runs])
    fid = fiducial(V, s)
    line = fit_line(V[fid], s["fit_cut_sigma"], s["fit_iterations"]) if fid.sum() >= 3 else None
    result = {"config": s["config"], "settings": {k: v for k, v in s.items() if k != "config"},
              "runs": [{"name": r["name"], "dirs": r["dirs"], "events": r["events"], "vertices": int(len(r["vtx"])),
                        "fiducial": int(fiducial(r["vtx"], s).sum()), "hot_strips": r["hot_strips"],
                        "frames_per_board": r["frames"], "trig_resets_per_board": r["resets"]} for r in runs],
              "vertices": int(len(V)), "fiducial": int(fid.sum()),
              "fiducial_mm": {"x": [-HALF_XY + s["margin"], HALF_XY - s["margin"]], "y": [-HALF_XY + s["margin"], HALF_XY - s["margin"]],
                              "z": [s["z_gem"] + s["margin"], Z_TOP - s["margin"]]},
              "beam_line": {k: v for k, v in line.items() if k not in ("use", "r")} if line else None}

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig, axs = plt.subplots(2, 3, figsize=(18, 11))
    lim = HALF_XY + 10
    zr = (-10, Z_TOP + 10)
    F = V[fid]
    for row, (P, tag) in enumerate(((V, f"all vertices (N={len(V):,})"), (F, f"fiducial (N={len(F):,})"))):
        for col, (i, j, xr, yr, lab) in enumerate(((0, 1, (-lim, lim), (-lim, lim), ("X [mm] (E)", "Y [mm] (S)")),
                                                   (0, 2, (-lim, lim), zr, ("X [mm] (E)", "Z [mm] (down)")),
                                                   (1, 2, (-lim, lim), zr, ("Y [mm] (S)", "Z [mm] (down)")))):
            ax = axs[row, col]
            if len(P):
                hh = ax.hist2d(P[:, i], P[:, j], bins=[np.linspace(*xr, 82), np.linspace(*yr, 82)], cmap="turbo",
                               norm=LogNorm(vmin=1))
                plt.colorbar(hh[3], ax=ax, label="vertices / bin")
            m = s["margin"]
            box = {0: (-HALF_XY + m, HALF_XY - m), 1: (-HALF_XY + m, HALF_XY - m), 2: (s["z_gem"] + m, Z_TOP - m)}
            ax.add_patch(plt.Rectangle((box[i][0], box[j][0]), box[i][1] - box[i][0], box[j][1] - box[j][0],
                                       fill=False, ec="k", ls="--", lw=1.2))
            if line and row == 1:
                c, d = np.array(line["point_mm"]), np.array(line["direction"])
                t = np.linspace(-500, 500, 3)
                ax.plot(c[i] + t * d[i], c[j] + t * d[j], "m-", lw=1.5)
            ax.set(xlim=xr, ylim=yr, xlabel=lab[0], ylabel=lab[1], title=tag)
            if j == 1 or j == 2:
                ax.invert_yaxis()
            ax.set_aspect("equal", adjustable="box")
    if line:
        txt = (f"beam line (fiducial vertices, {line['used']:,} used): angle from Z {line['angle_from_z_deg']:.1f} deg, "
               f"azimuth {line['azimuth_deg']:.1f} deg; point ({', '.join(f'{x:.1f}' for x in line['point_mm'])}) mm; "
               f"distance to the line R50/R68/R90 = {', '.join(f'{x:.1f}' for x in line['distance_r50_r68_r90_mm'])} mm")
    else:
        txt = "beam line: not fitted (fewer than 3 fiducial vertices)"
    names = ", ".join(r["name"] for r in runs)
    fig.suptitle(f"{os.path.basename(out)}: {names}\nvertex = earliest clock with anode and cathode hits in the same uPIC; "
                 f"fiducial margin {s['margin']:g} mm (black dashed box); beam line: magenta; v {s['velocity']} mm/clock, offset {s['clock_offset']}, z_gem {s['z_gem']}\n{txt}",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, .93))
    fig.savefig(f"{out}.png", dpi=100)

    with open(f"{out}_vertices.csv", "w") as f:
        f.write("# vertex = earliest clock with anode and cathode hits in the same uPIC; x, y, z [mm] (X: E, Y: S, Z: down)\n")
        f.write("run,epoch,trig,upic,clock,x,y,z,n_hits,fiducial,dist_to_line\n")
        k = 0
        for r in runs:
            fr = fiducial(r["vtx"], s)
            for e, u, c, v, nh, fd in zip(r["vtx_event"], r["vtx_upic"], r["vtx_clock"], r["vtx"], r["vtx_nhits"], fr):
                dist = ""
                if line:
                    dist = f"{np.linalg.norm(np.cross(v - np.array(line['point_mm']), np.array(line['direction']))):.2f}"
                f.write(f"{r['name']},{e >> 32},{e & 0xffffffff},{u},{c},{v[0]:.2f},{v[1]:.2f},{v[2]:.2f},{nh},{int(fd)},{dist}\n")
                k += 1
    with open(f"{out}.json", "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(json.dumps({"vertices": result["vertices"], "fiducial": result["fiducial"], "beam_line": result["beam_line"]}, indent=2))
    print(f"saved {out}.png, {out}_vertices.csv, {out}.json")


if __name__ == "__main__":
    main()
