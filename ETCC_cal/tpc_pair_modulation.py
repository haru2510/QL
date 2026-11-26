#!/usr/bin/env python3
"""TPCで散乱点から2本の枝が出るイベント（対生成の候補）の e+e- 面の向きを、偏光の違うデータで比べる。

1. 散乱点: tpc_electron_modulation.py と同じ（Z <= z_veto のrawヒットなし、ビーム軸から --beam-radius 以内）。
2. 枝: 散乱点から --radius 以内、3 mm以上下のrawストリップを、XZとYZで別々に見る。各投影で
   - 1本の直線（tpc_electron_modulation.one_slope）
   - 散乱点を通る2本の直線（角度の投票で2本を初期化し、近い方への割当てと直線fitを繰り返す）
   の両方をfitして保存する。どちらを採るかは解析のときに決める（--rescan なしで条件を変えられる）。
3. 分類（検出器のX・Y、ビームはほぼZに平行）:
   - DD: XZもYZも2本 → e+e- の横方向の差は (Δx', ±Δy')。XとYの枝の組み合わせが2通りあり、どちらが正しいか
     わからないので、X軸に対して折り返した角 ψ = atan(|Δy'| / |Δx'|)（0°: 面がX、90°: 面がY）を使う。
   - DS: XZは2本、YZは1本 → 面はほぼXZ（ψ ≈ 0°）。 SD: その逆で ψ ≈ 90°。
   - SS: 両方1本（対生成として使わない）。
4. 比べる: 面がXに近いもの N_X = DS + DD(ψ<45°) と Yに近いもの N_Y = SD + DD(ψ>=45°) の非対称
   Q = (N_X − N_Y)/(N_X + N_Y) をデータごとに出し、差を取る。DDのψ分布は割り算して A[1 + μ cos 2ψ] でfitする。
   折り返しているので、偏光がX・Yに対して45°のときの情報（sin 2φ）は失われる。

ランごとの結果は <result_dir>/with_z0_cut/tpc_pair_branches_r<radius>/<run>.npz に保存し、あれば読み直さない。
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
from scipy.optimize import curve_fit

from pixel_polarization import fitted_axis
from settings import load_config, repo_path
from stage_coincidence_xy import scan_run
from tpc_electron_modulation import INLIER_MM, one_slope, read_vertices
from tpc_projected_pair import projected_hits
from vertex_sphere import take_local_cells

SCAN_MIN_SEPARATION_DEG = 3.0  # 保存する2本の枝の角度差の下限（解析では --min-separation で絞る）


@njit
def two_branch(hits, axis):
    """1つの投影で散乱点を通る2本の直線。[傾き0, 傾き1, 点数0, 点数1, Z幅0, Z幅1, Zの最大0, Zの最大1]（傾きの小さい順）。"""
    out = np.full(8, np.nan)
    n_all = 0
    for k in range(hits.shape[0]):
        if hits[k, 0] == axis:
            n_all += 1
    if n_all < 12:
        return out
    z = np.empty(n_all)
    t = np.empty(n_all)
    j = 0
    for k in range(hits.shape[0]):
        if hits[k, 0] == axis:
            z[j] = hits[k, 1]
            t[j] = hits[k, 2]
            j += 1
    directions = np.empty(2)
    available = np.ones(n_all, np.bool_)
    for branch in range(2):  # 1本目を投票で決め、その点を除いてもう1本
        best_count, best_angle = -1, 0.0
        for grid in range(-170, 171):  # 0.5°刻み
            angle = grid * np.pi / 360.0
            count = np.sum(available & (np.abs(t * np.cos(angle) - z * np.sin(angle)) <= INLIER_MM))
            if count > best_count:
                best_count, best_angle = count, angle
        directions[branch] = best_angle
        available &= np.abs(t * np.cos(best_angle) - z * np.sin(best_angle)) > INLIER_MM
    masks = np.zeros((2, n_all), np.bool_)
    for _ in range(6):  # 近い方の枝に割り当てて、それぞれ直線fit
        d0 = np.abs(t * np.cos(directions[0]) - z * np.sin(directions[0]))
        d1 = np.abs(t * np.cos(directions[1]) - z * np.sin(directions[1]))
        masks[0] = (d0 <= d1) & (d0 <= INLIER_MM)
        masks[1] = (d1 < d0) & (d1 <= INLIER_MM)
        for branch in range(2):
            m = masks[branch]
            if m.sum() < 6 or np.ptp(z[m]) < 5.0:
                return out
            zz, tt = z[m], t[m]
            directions[branch] = np.arctan(np.sum(zz * tt) / np.sum(zz * zz))
    if abs(directions[0] - directions[1]) < SCAN_MIN_SEPARATION_DEG * np.pi / 180:
        return out
    order = np.argsort(directions)
    for k in range(2):
        b = order[k]
        m = masks[b]
        out[k] = np.tan(directions[b])
        out[2 + k] = m.sum()
        out[4 + k] = np.ptp(z[m])
        out[6 + k] = z[m].max()
    return out


def scan(config, run, args, beam_point, beam, target):
    """1ランの散乱点ごとに、XZ・YZの1本fitと2本fitを保存する。"""
    vertices = read_vertices(config, run, beam_point, beam, args.beam_radius)
    ids = np.array(sorted(vertices), dtype=np.int64)
    lookup = {int(e): i for i, e in enumerate(ids)}
    vertex = np.array([vertices[int(e)] for e in ids]).reshape(-1, 3)
    calibration = tuple(config[k] for k in ("velocity", "clock_offset", "z_gem"))
    parts = {}

    def process(event, upic, occupancy, occupied, numbers):
        i = lookup.get(event)
        if i is None:
            take_local_cells(occupancy, occupied, numbers, np.zeros(3), 0.0, *calibration)
            return 0
        cells = take_local_cells(occupancy, occupied, numbers, vertex[i], args.radius, *calibration)
        hits = projected_hits(cells, vertex[i], args.radius, *calibration)
        if hits:
            parts.setdefault(i, []).append(np.asarray(hits))
        return len(hits)

    print(f"{run}: {len(ids):,} vertices within {args.beam_radius:g} mm of the beam axis after the Z veto", flush=True)
    scan_run(run, config, [ids], event_processor=process)
    single = np.full((len(ids), 2, 5), np.nan)  # 軸ごとに slope, n_inlier, n_all, z_span, rms
    double = np.full((len(ids), 2, 8), np.nan)  # 軸ごとに two_branch の8個
    for i, chunks in parts.items():
        hits = np.unique(np.concatenate(chunks), axis=0)
        for axis in range(2):
            single[i, axis] = one_slope(hits, axis)
            double[i, axis] = two_branch(hits, axis)
    temp = target.with_suffix(".tmp.npz")
    np.savez_compressed(temp, ids=ids, vertex=vertex, single=single, double=double,
                        info=json.dumps({"run": run, "beam_radius_mm": args.beam_radius, "radius_mm": args.radius,
                                         "inlier_mm": INLIER_MM, "scan_min_separation_deg": SCAN_MIN_SEPARATION_DEG}))
    temp.replace(target)
    print("SAVED", target, flush=True)


def classify(single, double, args):
    """各イベントの分類（0: なし, 1: SS, 2: DS, 3: SD, 4: DD）と、DDの折り返し角 ψ [deg]。"""
    n_all = single[:, :, 2]
    one = (np.isfinite(single[:, :, 0]) & (single[:, :, 3] >= args.min_span)
           & (single[:, :, 1] / np.maximum(n_all, 1) >= args.min_inlier_fraction))
    sep = np.degrees(np.abs(np.arctan(double[:, :, 1]) - np.arctan(double[:, :, 0])))
    two = (np.isfinite(double[:, :, 0]) & (sep >= args.min_separation)
           & (np.minimum(double[:, :, 4], double[:, :, 5]) >= args.min_branch_span)
           & ((double[:, :, 2] + double[:, :, 3]) / np.maximum(n_all, 1) >= args.min_inlier_fraction))
    two = np.nan_to_num(two, nan=False).astype(bool)
    cls = np.zeros(len(single), np.int8)
    cls[one[:, 0] & one[:, 1] & ~two[:, 0] & ~two[:, 1]] = 1
    cls[two[:, 0] & ~two[:, 1] & one[:, 1]] = 2
    cls[~two[:, 0] & one[:, 0] & two[:, 1]] = 3
    cls[two[:, 0] & two[:, 1]] = 4
    dx = double[:, 0, 1] - double[:, 0, 0]
    dy = double[:, 1, 1] - double[:, 1, 0]
    psi = np.where(cls == 4, np.degrees(np.arctan2(np.abs(dy), np.abs(dx))), np.nan)
    return cls, psi


def load_dataset(config_path, runs, args, beam_point, beam):
    config = load_config(config_path)
    if runs:
        missing = set(runs) - set(config["runs"])
        if missing:
            raise SystemExit(f"{config_path}: runs {sorted(missing)} are not in the config")
        config["runs"] = list(runs)
    folder = repo_path(config["result_dir"]) / f"with_z0_cut/tpc_pair_branches_r{args.radius:g}"
    folder.mkdir(parents=True, exist_ok=True)
    cls, psi = [], []
    for run in config["runs"]:
        target = folder / f"{run}.npz"
        if args.rescan or not target.exists():
            scan(config, run, args, beam_point, beam, target)
        if args.scan_only:
            continue
        with np.load(target) as saved:
            c, p = classify(saved["single"], saved["double"], args)
        cls.append(c)
        psi.append(p)
    if args.scan_only:
        return None
    return np.concatenate(cls), np.concatenate(psi)


def cos2(psi, a, mu):
    return a * (1 + mu * np.cos(2 * np.radians(psi)))


def asymmetry(cls, psi):
    nx = int(np.sum(cls == 2) + np.sum((cls == 4) & (psi < 45)))
    ny = int(np.sum(cls == 3) + np.sum((cls == 4) & (psi >= 45)))
    n = nx + ny
    q = (nx - ny) / n
    return {"N_X": nx, "N_Y": ny, "Q": q, "Q_err": float(np.sqrt((1 - q * q) / n))}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--a", type=Path, required=True, help="分子のデータの設定")
    parser.add_argument("--b", type=Path, help="分母のデータの設定（--scan-only なら省略可）")
    parser.add_argument("--runs-a", nargs="+", help="--a の設定のうち使うラン（既定: 全ラン）")
    parser.add_argument("--runs-b", nargs="+", help="--b の設定のうち使うラン（既定: 全ラン）")
    parser.add_argument("--label-a", default="A")
    parser.add_argument("--label-b", default="B")
    parser.add_argument("--beam-radius", type=float, default=5.0, help="散乱点とビーム軸の距離の上限 [mm]")
    parser.add_argument("--radius", type=float, default=40.0, help="枝を見る散乱点からの距離 [mm]")
    parser.add_argument("--min-span", type=float, default=10.0, help="1本fitのZ幅の下限 [mm]")
    parser.add_argument("--min-branch-span", type=float, default=25.0, help="2本fitの各枝のZ幅の下限 [mm]（短いと1本の曲がった飛跡を2本に割った偽物が入る）")
    parser.add_argument("--min-separation", type=float, default=6.0, help="1つの投影での2本の角度差の下限 [deg]")
    parser.add_argument("--min-inlier-fraction", type=float, default=0.7, help="直線（1本または2本）に乗ったストリップの割合の下限")
    parser.add_argument("--bins", type=int, default=9, help="ψ（0-90°）のビン数")
    parser.add_argument("--out", type=Path, help="保存先")
    parser.add_argument("--rescan", action="store_true")
    parser.add_argument("--scan-only", action="store_true", help="--a のランをスキャンして保存するだけ")
    args = parser.parse_args()

    beam_point, beam = fitted_axis(repo_path(load_config(args.a)["inputs"]["polarization_beam_axis"]))
    if args.scan_only:
        load_dataset(args.a, args.runs_a, args, beam_point, beam)
        return
    labels = [args.label_a, args.label_b]
    data = [load_dataset(args.a, args.runs_a, args, beam_point, beam), load_dataset(args.b, args.runs_b, args, beam_point, beam)]
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    colors = ["#2a78d6", "#eb6834"]
    names = {0: "none", 1: "SS (1+1)", 2: "DS (2 in XZ, 1 in YZ)", 3: "SD (1 in XZ, 2 in YZ)", 4: "DD (2+2)"}
    result = {"labels": labels, "conditions": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
              "classes": [{names[k]: int(np.sum(c == k)) for k in names} for c, _ in data],
              "asymmetry": [asymmetry(c, p) for c, p in data]}
    qa, qb = result["asymmetry"]
    result["Q_difference"] = {"value": qa["Q"] - qb["Q"], "err": float(np.hypot(qa["Q_err"], qb["Q_err"]))}

    edges = np.linspace(0, 90, args.bins + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    hists = [np.histogram(p[c == 4], bins=edges)[0].astype(float) for c, p in data]
    na, nb = hists[0].sum(), hists[1].sum()
    ratio = (hists[0] / na) / np.where(hists[1] > 0, hists[1] / nb, np.nan)
    err = ratio * np.sqrt(1 / np.maximum(hists[0], 1) + 1 / np.maximum(hists[1], 1))
    p, cov = curve_fit(cos2, centres, ratio, p0=[1, 0], sigma=err, absolute_sigma=True)
    chi2 = float(np.sum(((ratio - cos2(centres, *p)) / err) ** 2))
    result["DD_ratio_fit"] = {"A": p[0], "mu": p[1], "mu_err": float(np.sqrt(cov[1, 1])), "chi2": chi2, "ndf": args.bins - 2,
                              "psi_edges_deg": edges.tolist(), "counts": [h.astype(int).tolist() for h in hists]}

    fig, axes = plt.subplots(1, 4, figsize=(22, 5.2), constrained_layout=True)
    ax = axes[0]
    keys = [1, 2, 3, 4]
    width = 0.38
    for k, ((c, _), label, color) in enumerate(zip(data, labels, colors)):
        frac = [100 * np.sum(c == key) / len(c) for key in keys]
        bars = ax.bar(np.arange(4) + (k - 0.5) * width, frac, width=width - 0.03, color=color, label=f"{label}: {len(c):,} vertices")
        for bar, f in zip(bars, frac):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{f:.1f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(np.arange(4), [names[k].replace(" (", "\n(") for k in keys], fontsize=8)
    ax.set(ylabel="% of beam-axis vertices", title="Event classes")
    ax.legend(fontsize=8, frameon=False)
    ax = axes[1]
    for h, label, color in zip(hists, labels, colors):
        n = h.sum()
        ax.errorbar(centres, h / n * args.bins, yerr=np.sqrt(h) / n * args.bins, fmt="o-", ms=4, color=color, label=f"{label}: {int(n):,}")
    ax.set(xlabel="ψ = atan(|Δy'| / |Δx'|) [deg]  (0: e+e- plane ∥ X, 90: ∥ Y)", ylabel="normalized (mean = 1)",
           title="DD: folded pair-plane angle", xlim=(0, 90))
    ax.legend(fontsize=8, frameon=False)
    ax = axes[2]
    ax.errorbar(centres, ratio, yerr=err, fmt="o", color="black", ms=4)
    xx = np.linspace(0, 90, 181)
    ax.plot(xx, cos2(xx, *p), color="#e34948", label=f"μ = {p[1]:.3f} ± {np.sqrt(cov[1, 1]):.3f}, χ²/ndf = {chi2:.1f}/{args.bins - 2}")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set(xlabel="ψ [deg]", ylabel=f"{labels[0]} / {labels[1]} (each normalized)", title="DD ratio, fit A[1 + μ cos 2ψ]", xlim=(0, 90))
    ax.legend(fontsize=8, frameon=False)
    ax = axes[3]
    for k, (q, label, color) in enumerate(zip(result["asymmetry"], labels, colors)):
        ax.errorbar([k], [q["Q"]], yerr=[q["Q_err"]], fmt="o", color=color, ms=7, capsize=4)
        ax.text(k + 0.08, q["Q"], f"{q['Q']:+.4f} ± {q['Q_err']:.4f}\nN_X {q['N_X']:,} / N_Y {q['N_Y']:,}", va="center", fontsize=8)
    d = result["Q_difference"]
    ax.set_xticks([0, 1], labels, fontsize=8)
    ax.set_xlim(-0.5, 1.9)
    ax.set(ylabel="Q = (N_X − N_Y) / (N_X + N_Y)", title=f"Plane ∥X vs ∥Y (DS + DD ψ<45 vs SD + DD ψ≥45)\nΔQ = {d['value']:+.4f} ± {d['err']:.4f}")
    fig.suptitle(f"TPC two-branch (pair candidate) events at the beam axis: {labels[0]} vs {labels[1]}   "
                 f"(r ≤ {args.radius:g} mm, branch sep ≥ {args.min_separation:g}°, branch Z span ≥ {args.min_branch_span:g} mm, "
                 f"inlier ≥ {args.min_inlier_fraction:g})")
    fig.savefig(out / "pair_modulation.png", dpi=120)
    (out / "pair_modulation.json").write_text(json.dumps(result, indent=1, default=float) + "\n")
    print("Saved", out / "pair_modulation.png")
    for label, cl, q in zip(labels, result["classes"], result["asymmetry"]):
        print(label, cl, f"Q = {q['Q']:+.4f} ± {q['Q_err']:.4f}")
    print(f"ΔQ = {d['value']:+.4f} ± {d['err']:.4f};  DD ratio μ = {p[1]:.4f} ± {np.sqrt(cov[1, 1]):.4f}, χ²/ndf = {chi2:.1f}/{args.bins - 2}")


if __name__ == "__main__":
    main()
