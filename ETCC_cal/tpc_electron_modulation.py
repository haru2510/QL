#!/usr/bin/env python3
"""TPCの反跳電子の初めの向きの方位角分布（モジュレーションカーブ）を、偏光の違う2つのデータで比べる。

1. 散乱点: vertices_all のうち、Z <= z_veto のrawヒットを持つイベントを除き（veto_dirのz0_track_veto_<run>.npz）、
   フィットしたビーム軸（inputs.polarization_beam_axis）から --beam-radius 以内のもの。
2. 向き: 散乱点から --radius 以内、散乱点より3 mm以上下のrawストリップを、XZとYZで別々に、散乱点を通る1本の直線で
   fitする（角度の投票で初期値を決め、1.6 mm以内のストリップで最小二乗を繰り返す）。XとYを組にしないので、
   同じクロックの組合せの曖昧さが入らない。dx/dz と dy/dz から、ビーム軸に対する天頂角θと方位角φ（0-360°、
   ビーム軸に垂直な面に射影した検出器のX軸が0°、Yが90°）を出す。
3. 比べる: θ >= θmin の φ 分布を2つのデータで割り、A [1 + μ cos 2(φ − φ0)] でfitする。偏光の向きが90°違う2つなら、
   比のμは1つのデータのモジュレーションの約2倍。検出器の形（uPICの隙間、死んだストリップ、X/Yの差）は割り算で消える。

ランごとの結果は <result_dir>/with_z0_cut/tpc_electron_direction/<run>.npz に保存し、あれば読み直さない（--rescan で作り直す）。
"""
import argparse
import csv
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
from tpc_projected_pair import projected_hits
from vertex_sphere import take_local_cells

INLIER_MM = 1.6  # 直線からこの距離以内のストリップを飛跡とみなす（tpc_projected_pair.two_slopesと同じ）


@njit
def one_slope(hits, axis):
    """1つの投影 (dz, 横の座標) の、散乱点を通る1本の直線の傾き・使った点の数・全点数・Zの幅・残差のRMS。"""
    n_all = 0
    for k in range(hits.shape[0]):
        if hits[k, 0] == axis:
            n_all += 1
    if n_all < 6:
        return np.nan, 0, n_all, 0.0, np.nan
    z = np.empty(n_all)
    t = np.empty(n_all)
    j = 0
    for k in range(hits.shape[0]):
        if hits[k, 0] == axis:
            z[j] = hits[k, 1]
            t[j] = hits[k, 2]
            j += 1
    best_count, best_angle = -1, 0.0
    for grid in range(-85, 86):  # 1°刻みの投票
        angle = grid * np.pi / 180.0
        count = np.sum(np.abs(t * np.cos(angle) - z * np.sin(angle)) <= INLIER_MM)
        if count > best_count:
            best_count, best_angle = count, angle
    slope = np.tan(best_angle)
    mask = np.abs(t - slope * z) * np.cos(best_angle) <= INLIER_MM
    for _ in range(5):
        if mask.sum() < 6:
            return np.nan, int(mask.sum()), n_all, 0.0, np.nan
        zz, tt = z[mask], t[mask]
        slope = np.sum(zz * tt) / np.sum(zz * zz)
        angle = np.arctan(slope)
        mask = np.abs(t - slope * z) * np.cos(angle) <= INLIER_MM
    zz, tt = z[mask], t[mask]
    if len(zz) < 6:
        return np.nan, len(zz), n_all, 0.0, np.nan
    rms = np.sqrt(np.mean(((tt - slope * zz) * np.cos(np.arctan(slope))) ** 2))
    return slope, len(zz), n_all, zz.max() - zz.min(), rms


def read_vertices(config, run, beam_point, beam, beam_radius):
    """Z側除外を通り、ビーム軸から beam_radius 以内の散乱点 {event: (x, y, z)}。"""
    veto = np.load(repo_path(config["outputs"]["veto_dir"]) / f"z0_track_veto_{run}.npz")["event"]
    out = {}
    with repo_path(config["inputs"]["vertices_all"]).open() as stream:
        for row in csv.DictReader(line for line in stream if not line.startswith("#")):
            if row["run"] != run:
                continue
            p = np.array([float(row[a]) for a in "xyz"])
            if p[2] <= config["z_veto"] or np.linalg.norm(np.cross(p - beam_point, beam)) > beam_radius:
                continue
            out.setdefault((int(row["epoch"]) << 32) | int(row["trig"]), p)
    events = np.array(sorted(out), dtype=np.int64)
    keep = ~np.isin(events, veto)
    return {int(e): out[int(e)] for e in events[keep]}


def scan(config, run, args, beam_point, beam, target):
    """1ランの散乱点ごとに、XZ・YZの傾きを出して保存する。"""
    vertices = read_vertices(config, run, beam_point, beam, args.beam_radius)
    ids = np.array(sorted(vertices), dtype=np.int64)
    lookup = {int(e): i for i, e in enumerate(ids)}
    vertex = np.array([vertices[int(e)] for e in ids]).reshape(-1, 3)
    calibration = tuple(config[k] for k in ("velocity", "clock_offset", "z_gem"))
    parts = {}

    def process(event, upic, occupancy, occupied, numbers):
        i = lookup.get(event)
        if i is None:  # 選んでいないイベント（scan_runは選んだものだけ渡すが念のため）。占有は消す。
            take_local_cells(occupancy, occupied, numbers, np.zeros(3), 0.0, *calibration)
            return 0
        cells = take_local_cells(occupancy, occupied, numbers, vertex[i], args.radius, *calibration)
        hits = projected_hits(cells, vertex[i], args.radius, *calibration)
        if hits:
            parts.setdefault(i, []).append(np.asarray(hits))
        return len(hits)

    print(f"{run}: {len(ids):,} vertices within {args.beam_radius:g} mm of the beam axis after the Z veto", flush=True)
    scan_run(run, config, [ids], event_processor=process)
    fit = np.full((len(ids), 2, 5), np.nan)  # 軸ごとに slope, n_inlier, n_all, z_span, rms
    for i, chunks in parts.items():
        hits = np.unique(np.concatenate(chunks), axis=0)  # 複数uPICの重なりを1回に
        for axis in range(2):
            fit[i, axis] = one_slope(hits, axis)
    temp = target.with_suffix(".tmp.npz")
    np.savez_compressed(temp, ids=ids, vertex=vertex, fit=fit,
                        info=json.dumps({"run": run, "beam_radius_mm": args.beam_radius, "radius_mm": args.radius,
                                         "inlier_mm": INLIER_MM}))
    temp.replace(target)
    print("SAVED", target, flush=True)


def angles(fit, beam, x_ref, y_ref):
    """傾きからビーム軸に対する θ [deg] と φ [deg, 0-360]。"""
    d = np.stack([fit[:, 0, 0], fit[:, 1, 0], np.ones(len(fit))], axis=1)
    d /= np.linalg.norm(d, axis=1)[:, None]
    theta = np.degrees(np.arccos(np.clip(d @ beam, -1, 1)))
    phi = np.degrees(np.arctan2(d @ y_ref, d @ x_ref)) % 360
    return theta, phi


def load_dataset(config_path, args, beam_point, beam, x_ref, y_ref, runs=None):
    """1つのデータ（設定）のラン（runs 指定時はその一部）を読み、品質の条件を通った θ, φ を返す。"""
    config = load_config(config_path)
    if runs:  # 同じ偏光の中でランを分ける比較（系統誤差の確認）に使う。
        missing = set(runs) - set(config["runs"])
        if missing:
            raise SystemExit(f"{config_path}: runs {sorted(missing)} are not in the config")
        config["runs"] = list(runs)
    folder = repo_path(config["result_dir"]) / "with_z0_cut/tpc_electron_direction"
    folder.mkdir(parents=True, exist_ok=True)
    theta, phi, counts = [], [], {"vertices": 0, "both_projections_fitted": 0, "single_track_quality": 0}
    for run in config["runs"]:
        target = folder / f"{run}.npz"
        if args.rescan or not target.exists():
            scan(config, run, args, beam_point, beam, target)
        with np.load(target) as saved:
            fit = saved["fit"]
        ok = np.isfinite(fit[:, :, 0]).all(axis=1) & (fit[:, :, 3] >= args.min_span).all(axis=1)
        frac = fit[:, :, 1] / np.maximum(fit[:, :, 2], 1)
        good = ok & (frac >= args.min_inlier_fraction).all(axis=1)
        counts["vertices"] += len(fit)
        counts["both_projections_fitted"] += int(ok.sum())
        counts["single_track_quality"] += int(good.sum())
        t, p = angles(fit[good], beam, x_ref, y_ref)
        theta.append(t)
        phi.append(p)
    return config, np.concatenate(theta), np.concatenate(phi), counts


def modulation(phi, a, mu, phi0):
    return a * (1 + mu * np.cos(2 * np.radians(phi - phi0)))


def fit_curve(centres, values, errors):
    ok = errors > 0
    p0 = [values[ok].mean(), 0.05, 0.0]
    best = None
    for start in (0.0, 45.0, 90.0, 135.0):  # φ0の初期値を変えて一番よいもの
        try:
            p, cov = curve_fit(modulation, centres[ok], values[ok], p0=[p0[0], 0.05, start], sigma=errors[ok], absolute_sigma=True, maxfev=20000)
        except RuntimeError:
            continue
        chi2 = float(np.sum(((values[ok] - modulation(centres[ok], *p)) / errors[ok]) ** 2))
        if best is None or chi2 < best[2]:
            best = (p, cov, chi2)
    p, cov, chi2 = best
    if p[1] < 0:  # μ >= 0 に、φ0 を 90°回して合わせる
        p[1], p[2] = -p[1], p[2] + 90
    p[2] %= 180
    err = np.sqrt(np.diag(cov))
    return {"A": p[0], "mu": p[1], "mu_err": err[1], "phi0_deg": p[2], "phi0_err_deg": err[2], "chi2": chi2, "ndf": int(ok.sum()) - 3}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=Path, default=Path(__file__).parent / "config/20261124_run06_07.json", help="分子のデータ（既定: 水平 1124 run06・07）")
    parser.add_argument("--b", type=Path, default=Path(__file__).parent / "config/20261125_vertical.json", help="分母のデータ（既定: 垂直 1125）")
    parser.add_argument("--runs-a", nargs="+", help="--a の設定のうち使うラン（既定: 全ラン）")
    parser.add_argument("--runs-b", nargs="+", help="--b の設定のうち使うラン（既定: 全ラン）")
    parser.add_argument("--label-a", default="horizontal (1124 run06-07)")
    parser.add_argument("--label-b", default="vertical (1125 run00-05)")
    parser.add_argument("--beam-radius", type=float, default=5.0, help="散乱点とビーム軸の距離の上限 [mm]")
    parser.add_argument("--radius", type=float, default=30.0, help="向きを出すのに使う散乱点からの距離 [mm]")
    parser.add_argument("--min-span", type=float, default=10.0, help="各投影で直線に乗ったストリップのZの幅の下限 [mm]")
    parser.add_argument("--min-inlier-fraction", type=float, default=0.7, help="各投影で直線に乗ったストリップの割合の下限（1本の飛跡らしさ）")
    parser.add_argument("--theta-min", type=float, nargs="+", default=[10.0, 20.0, 30.0], help="使う天頂角の下限 [deg]（複数なら行を分けて描く）")
    parser.add_argument("--bins", type=int, default=18, help="φのビン数（0-360°）")
    parser.add_argument("--out", type=Path, help="保存先（既定: output/20261125_vertical_vs_horizontal/tpc_electron_modulation）")
    parser.add_argument("--rescan", action="store_true", help="ランごとの保存を作り直す")
    args = parser.parse_args()

    config_a = load_config(args.a)
    beam_point, beam = fitted_axis(repo_path(config_a["inputs"]["polarization_beam_axis"]))
    x_ref = np.array([1.0, 0.0, 0.0]) - beam[0] * beam
    x_ref /= np.linalg.norm(x_ref)
    y_ref = np.cross(beam, x_ref)
    data = [load_dataset(args.a, args, beam_point, beam, x_ref, y_ref, args.runs_a), load_dataset(args.b, args, beam_point, beam, x_ref, y_ref, args.runs_b)]
    out = args.out or repo_path("ETCC_cal/output/20261125_vertical_vs_horizontal/tpc_electron_modulation")
    out.mkdir(parents=True, exist_ok=True)

    edges = np.linspace(0, 360, args.bins + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    labels = [args.label_a, args.label_b]
    colors = ["#2a78d6", "#eb6834"]
    rows = len(args.theta_min)
    fig, axes = plt.subplots(rows + 1, 3, figsize=(18, 4.2 * (rows + 1)), constrained_layout=True)
    result = {"labels": labels, "conditions": {k: v for k, v in vars(args).items() if k not in ("out",)}, "counts": [d[3] for d in data], "theta_min": {}}
    result["conditions"] = {k: (str(v) if isinstance(v, Path) else v) for k, v in result["conditions"].items()}

    ax = axes[0, 0]
    tb = np.linspace(0, 90, 46)
    for (config, theta, phi, counts), label, color in zip(data, labels, colors):
        ax.hist(theta, bins=tb, histtype="step", color=color, lw=1.5, density=True, label=f"{label}: {len(theta):,}")
    for t in args.theta_min:
        ax.axvline(t, color="0.5", ls=":")
    ax.set(xlabel="θ from the beam axis [deg]", ylabel="density", title="Initial electron direction: polar angle")
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    for (config, theta, phi, counts), label, color in zip(data, labels, colors):
        ax.hist(phi, bins=edges, histtype="step", color=color, lw=1.5, density=True, label=label)
    ax.set(xlabel="φ [deg] (0 = detector X, 90 = Y)", ylabel="density", title="φ, all θ (detector acceptance included)", xlim=(0, 360))
    ax.legend(fontsize=8)
    axes[0, 2].axis("off")
    text = ["Selection (per dataset):"]
    for label, d in zip(labels, data):
        c = d[3]
        text.append(f"{label}: vertices {c['vertices']:,} → both projections fitted {c['both_projections_fitted']:,} → single-track {c['single_track_quality']:,}")
    text += ["", f"vertex: Z veto (no raw hit at Z <= 15 mm), within {args.beam_radius:g} mm of the fitted beam axis",
             f"direction: strips within {args.radius:g} mm of the vertex, one line through the vertex per projection",
             f"quality: Z span >= {args.min_span:g} mm and inlier fraction >= {args.min_inlier_fraction:g} in both XZ and YZ"]
    axes[0, 2].text(0, 1, "\n".join(text), va="top", fontsize=8, wrap=True)

    for r, tmin in enumerate(args.theta_min, start=1):
        hists = []
        for (config, theta, phi, counts) in data:
            h, _ = np.histogram(phi[theta >= tmin], bins=edges)
            hists.append(h.astype(float))
        ax = axes[r, 0]
        for h, label, color in zip(hists, labels, colors):
            n = h.sum()
            ax.errorbar(centres, h / n * args.bins, yerr=np.sqrt(h) / n * args.bins, fmt="o-", ms=3, color=color, label=f"{label}: {int(n):,}")
        ax.set(xlabel="φ [deg]", ylabel="normalized counts (mean = 1)", title=f"θ >= {tmin:g}°: φ distributions", xlim=(0, 360))
        ax.legend(fontsize=8)
        fits = {}
        for k, (h, label) in enumerate(zip(hists, labels)):
            fits[label] = fit_curve(centres, h, np.sqrt(np.maximum(h, 1)))
        na, nb = hists[0].sum(), hists[1].sum()
        ratio = (hists[0] / na) / np.where(hists[1] > 0, hists[1] / nb, np.nan)
        err = ratio * np.sqrt(1 / np.maximum(hists[0], 1) + 1 / np.maximum(hists[1], 1))
        fr = fit_curve(centres, ratio, err)
        ax = axes[r, 1]
        ax.errorbar(centres, ratio, yerr=err, fmt="o", color="black", ms=4, label="data")
        xx = np.linspace(0, 360, 361)
        ax.plot(xx, modulation(xx, fr["A"], fr["mu"], fr["phi0_deg"]), color="#e34948",
                label=f"μ = {fr['mu']:.3f} ± {fr['mu_err']:.3f}, φ0 = {fr['phi0_deg']:.0f} ± {fr['phi0_err_deg']:.0f}°\nχ²/ndf = {fr['chi2']:.1f}/{fr['ndf']}")
        ax.axhline(1, color="0.6", lw=0.8)
        ax.set(xlabel="φ [deg]", ylabel=f"N({labels[0].split()[0]}) / N({labels[1].split()[0]}) (each normalized)",
               title=f"θ >= {tmin:g}°: ratio, fit A[1 + μ cos2(φ − φ0)]", xlim=(0, 360))
        ax.legend(fontsize=8)
        ax = axes[r, 2]
        ax.axis("off")
        lines = [f"θ >= {tmin:g}°", f"ratio: μ = {fr['mu']:.4f} ± {fr['mu_err']:.4f} (one dataset ≈ μ/2 = {fr['mu']/2:.4f})",
                 f"       φ0 = {fr['phi0_deg']:.1f} ± {fr['phi0_err_deg']:.1f}°, χ²/ndf = {fr['chi2']:.1f}/{fr['ndf']}", "",
                 "single dataset fits (include detector acceptance):"]
        for label, f in fits.items():
            lines.append(f"  {label}: μ = {f['mu']:.3f} ± {f['mu_err']:.3f}, φ0 = {f['phi0_deg']:.0f}°")
        ax.text(0, 1, "\n".join(lines), va="top", fontsize=9, family="monospace")
        result["theta_min"][f"{tmin:g}"] = {"events": [int(na), int(nb)], "ratio_fit": fr, "single_fits": fits,
                                             "phi_edges_deg": edges.tolist(), "counts": [h.astype(int).tolist() for h in hists]}
    fig.suptitle(f"TPC recoil-electron azimuth: {labels[0]} / {labels[1]}  (initial direction within {args.radius:g} mm of the vertex)")
    fig.savefig(out / "modulation.png", dpi=120)
    (out / "modulation.json").write_text(json.dumps(result, indent=1, default=float) + "\n")
    print("Saved", out / "modulation.png")
    for t, r in result["theta_min"].items():
        f = r["ratio_fit"]
        print(f"θ >= {t}°: N = {r['events']}, ratio μ = {f['mu']:.4f} ± {f['mu_err']:.4f}, φ0 = {f['phi0_deg']:.1f}°, χ²/ndf = {f['chi2']:.1f}/{f['ndf']}")


if __name__ == "__main__":
    main()
