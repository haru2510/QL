"""stage4: 比較（分子のサンプル ÷ 分母のサンプル）ごとに、モジュレーションと有意度を出す。

対生成（analysis.pair）:
  - Q = (N_X − N_Y)/(N_X + N_Y)、N_X = DS + DD(ψ < 45°)、N_Y = SD + DD(ψ >= 45°)。サンプルごとに出して差 ΔQ を取る。
  - DD の ψ（0〜90°）のヒストグラムを各サンプルで規格化して割り、A [1 + μ cos 2ψ] で fit する。
    折り返しているので、偏光が X・Y に対して φ0 傾いていると見える μ は μ_true · cos 2φ0（45° では 0）。
  - DS・DD・SD の図（pair_curve）: DS（ψ ≈ 0° 側）、DD の ψ の各ビン、SD（ψ ≈ 90° 側）の 11 点を、対生成の候補の総数で
    規格化して割る。fit は DD の点だけで A [1 + μ c_k]（c_k = ビン内の cos 2ψ の平均）。DS・SD は参照点で fit に使わない。
    DS・SD は ψ が範囲 [psi_lo, psi_hi] でしか決まらない。範囲内で ψ が一様と仮定して 11 点全部を fit した値も、
    参考として JSON（fit_all_points_reference）と summary に残す（DS・SD が DD より小さいモジュレーションを示し、合わない）。
  - 有意度: 平らな線との Δχ²（1自由度）と ΔQ / 誤差。
1本の飛跡（analysis.single_track、SS のみ）:
  - θ >= θmin ごとに φ（0〜360°）のヒストグラムの比を A [1 + μ cos 2(φ − φ0)] で fit（Δχ² は2自由度）。

アクセプタンス（検出器の形）は割り算で消す。分母には偏光の無い（円偏光の）データを使うのが基本。
全体図（analysis.overview）: 複数のサンプルを1枚に重ね、基準のサンプルに対する比を描く（DD だけで fit）。

出力: stage4_modulation/<cut_set>/overview.png, overview.json（全体図）
      stage4_modulation/<cut_set>/<comparison>/{pair.png, pair.json, pair_curve.png, pair_curve.json, single_track.png, single_track.json}, summary.md, summary.json
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .fitting import asymmetry, fit_cos2_fixed_phase, fit_cos2_free_phase, fit_cos2_mean, normalized_ratio
from .stage3_classify import CLASS_NAMES, DD, DS, SD, SS

COLORS = ("#2a78d6", "#eb6834")  # 分子, 分母


def load_classified(pipeline, cut_name, sample_name):
    with np.load(pipeline.stage_dir(3, cut_name) / f"{sample_name}.npz") as saved:
        keys = ("cls", "psi", "psi_lo", "psi_hi", "theta", "phi", "pairing", "phi_pair")
        return {key: saved[key] for key in keys if key in saved.files}


# ---------------------------------------------------------------- 対生成
def pair_result(data, labels, bins):
    """Q と DD の ψ の比の fit。"""
    edges = np.linspace(0, 90, bins + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    asym, hists, classes = [], [], []
    for d in data:
        cls, psi = d["cls"], d["psi"]
        dd = cls == DD
        n_x = np.sum(cls == DS) + np.sum(dd & (psi < 45))
        n_y = np.sum(cls == SD) + np.sum(dd & (psi >= 45))
        asym.append(asymmetry(n_x, n_y))
        hists.append(np.histogram(psi[dd], bins=edges)[0])
        classes.append({CLASS_NAMES[k]: int(np.sum(cls == k)) for k in CLASS_NAMES} | {"beam_axis_vertices": int(len(cls))})
    ratio, error, ok = normalized_ratio(hists[0], hists[1])
    fit, model = fit_cos2_fixed_phase(centres[ok], ratio[ok], error[ok])
    dq = asym[0]["Q"] - asym[1]["Q"]
    dq_err = float(np.hypot(asym[0]["Q_err"], asym[1]["Q_err"]))
    result = {"labels": labels, "classes": classes, "asymmetry": asym,
              "delta_Q": {"value": dq, "err": dq_err, "significance_sigma": abs(dq) / dq_err},
              "dd_psi_fit": fit, "psi_edges_deg": edges.tolist(), "dd_psi_counts": [h.tolist() for h in hists]}
    return result, (centres, hists, ratio, error, ok, model)


def plot_pair(result, plot_data, title, target):
    centres, hists, ratio, error, ok, model = plot_data
    labels = result["labels"]
    bins = len(centres)
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.2), constrained_layout=True)

    ax = axes[0]  # 分類の割合
    keys = ["SS", "DS", "SD", "DD"]
    for k, (classes, label, color) in enumerate(zip(result["classes"], labels, COLORS)):
        total = classes["beam_axis_vertices"]
        frac = [100 * classes[key] / total for key in keys]
        bars = ax.bar(np.arange(4) + (k - 0.5) * 0.38, frac, width=0.35, color=color, label=f"{label}: {total:,} vertices")
        for bar, f in zip(bars, frac):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{f:.1f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(np.arange(4), ["SS\n1+1", "DS\nXZ 2, YZ 1", "SD\nXZ 1, YZ 2", "DD\n2+2"], fontsize=8)
    ax.set(ylabel="% of beam-axis vertices", title="Event classes")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]  # DD の ψ 分布
    for h, label, color in zip(hists, labels, COLORS):
        n = h.sum()
        ax.errorbar(centres, h / n * bins, yerr=np.sqrt(h) / n * bins, fmt="o-", ms=4, color=color, label=f"{label}: {int(n):,}")
    ax.set(xlabel="ψ [deg]  (0: e+e- plane ∥ X, 90: ∥ Y)", ylabel="normalized (mean = 1)", title="DD: folded pair-plane angle", xlim=(0, 90))
    ax.legend(fontsize=8, frameon=False)

    ax = axes[2]  # DD の比と fit
    fit = result["dd_psi_fit"]
    half = 45.0 / bins  # 横棒 = ψ のビンの幅
    ax.errorbar(centres[ok], ratio[ok], yerr=error[ok], xerr=half, fmt="o", color="black", ms=4)
    xx = np.linspace(0, 90, 181)
    ax.plot(xx, model(xx, fit["A"], fit["mu"]), color="#e34948",
            label=f"μ = {fit['mu']:+.3f} ± {fit['mu_err']:.3f}  ({fit['significance_sigma']:.1f}σ vs flat)\n"
                  f"χ²/ndf = {fit['chi2']:.1f}/{fit['ndf']} (flat: {fit['chi2_flat']:.1f}/{fit['ndf_flat']})")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set(xlabel="ψ [deg]", ylabel="ratio (each normalized)", title=f"DD ratio {labels[0]} / {labels[1]}: A[1 + μ cos 2ψ]", xlim=(0, 90))
    ax.legend(fontsize=8, frameon=False)

    ax = axes[3]  # Q
    for k, (q, label, color) in enumerate(zip(result["asymmetry"], labels, COLORS)):
        ax.errorbar([k], [q["Q"]], yerr=[q["Q_err"]], fmt="o", color=color, ms=7, capsize=4)
        ax.text(k + 0.08, q["Q"], f"{q['Q']:+.4f} ± {q['Q_err']:.4f}\nN_X {q['N_X']:,} / N_Y {q['N_Y']:,}", va="center", fontsize=8)
    ax.axhline(0, color="0.6", lw=0.8)
    d = result["delta_Q"]
    ax.set_xticks([0, 1], labels, fontsize=8)
    ax.set_xlim(-0.5, 1.9)
    ax.set(ylabel="Q = (N_X − N_Y)/(N_X + N_Y)",
           title=f"plane ∥X (DS + DD ψ<45) vs ∥Y (SD + DD ψ≥45)\nΔQ = {d['value']:+.4f} ± {d['err']:.4f} ({d['significance_sigma']:.1f}σ)")
    fig.suptitle(title)
    fig.savefig(target, dpi=120)
    plt.close(fig)


def mean_cos2(lo_deg, hi_deg):
    """ψ が [lo, hi] で一様なときの cos 2ψ の平均（lo = hi ならその値）。"""
    lo, hi = np.radians(lo_deg), np.radians(hi_deg)
    width = hi - lo
    with np.errstate(divide="ignore", invalid="ignore"):
        averaged = (np.sin(2 * hi) - np.sin(2 * lo)) / (2 * width)
    return np.where(width > 1e-9, averaged, np.cos(2 * lo))


def pair_points(d, bins):
    """1サンプルの対生成の候補を点に分ける: 0 = DS, 1..bins = DD の ψ のビン, bins+1 = SD。

    戻り値: 各候補の点番号、各候補の cos 2ψ の平均（DS・SD は範囲内で一様と仮定）、ψ の範囲 [lo, hi]。
    """
    edges = np.linspace(0, 90, bins + 1)
    cls = d["cls"]
    point = np.full(len(cls), -1)
    point[cls == DS] = 0
    dd = cls == DD
    point[dd] = 1 + np.clip(np.digitize(d["psi"][dd], edges) - 1, 0, bins - 1)
    point[cls == SD] = bins + 1
    used = point >= 0
    lo, hi = d["psi_lo"][used], d["psi_hi"][used]
    return point[used], mean_cos2(lo, hi), np.column_stack([lo, hi])


def pair_curve_result(data, labels, bins):
    """DS・DD（ψ のビン）・SD の点を、対生成の候補の総数で規格化して割り、A [1 + μ c_k] で fit する。"""
    edges = np.linspace(0, 90, bins + 1)
    names = ["DS"] + [f"DD {edges[k]:.0f}-{edges[k + 1]:.0f}" for k in range(bins)] + ["SD"]
    counts, cos_sums, ranges = [], np.zeros(bins + 2), [[] for _ in range(bins + 2)]
    for d in data:
        point, c, psi_range = pair_points(d, bins)
        counts.append(np.bincount(point, minlength=bins + 2))
        cos_sums += np.bincount(point, weights=c, minlength=bins + 2)
        for k in (0, bins + 1):  # DS・SD の ψ の範囲（図の横棒用）
            ranges[k].append(psi_range[point == k])
    total = counts[0] + counts[1]
    c_k = cos_sums / np.maximum(total, 1)  # 2つのサンプルを合わせた点ごとの cos 2ψ の平均
    ratio, error, ok = normalized_ratio(counts[0], counts[1])
    dd = np.zeros(bins + 2, bool)
    dd[1:bins + 1] = True
    fit = fit_cos2_mean(c_k[ok & dd], ratio[ok & dd], error[ok & dd])          # 主: DD の点だけ
    fit_all = fit_cos2_mean(c_k[ok], ratio[ok], error[ok])                  # 参考: DS・SD も入れた 11 点
    psi_eff = np.degrees(0.5 * np.arccos(np.clip(c_k, -1, 1)))  # cos 2ψ_eff = c_k となる ψ（図の横軸）
    spans = {}
    for k, name in ((0, "DS"), (bins + 1, "SD")):
        both = np.concatenate(ranges[k])
        spans[name] = [float(np.median(both[:, 0])), float(np.median(both[:, 1]))]
    result = {"labels": labels, "points": names, "counts": [c.tolist() for c in counts], "mean_cos2": c_k.tolist(),
              "psi_effective_deg": psi_eff.tolist(), "ratio": ratio.tolist(), "ratio_err": error.tolist(),
              "median_psi_range_deg": spans, "fit": fit, "fit_all_points_reference": fit_all}
    return result, (psi_eff, counts, ratio, error, ok)


def plot_pair_curve(result, plot_data, title, target):
    psi_eff, counts, ratio, error, ok = plot_data
    labels, fit = result["labels"], result["fit"]
    n_points = len(psi_eff)
    edge_points = {0: "DS", n_points - 1: "SD"}
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.4), constrained_layout=True)

    ax = axes[0]  # 各点の割合（対生成の候補の総数で規格化）
    for k, (c, label, color) in enumerate(zip(counts, labels, COLORS)):
        frac = 100 * c / c.sum()
        ax.errorbar(psi_eff + (k - 0.5) * 0.8, frac, yerr=100 * np.sqrt(c) / c.sum(), fmt="o", ms=5, color=color,
                    label=f"{label}: {int(c.sum()):,} pair candidates")
    for k, name in edge_points.items():
        ax.annotate(name, (psi_eff[k], 100 * counts[0][k] / counts[0].sum()), textcoords="offset points", xytext=(6, 6), fontsize=9)
    ax.set(xlabel="ψ [deg] (DS, SD: effective ψ with cos 2ψ = mean over the allowed range)", ylabel="% of pair candidates",
           title="Pair candidates: DS, DD (ψ bins), SD", xlim=(-3, 93))
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]  # 比と fit
    dd = np.array([k not in edge_points for k in range(n_points)])
    edges = np.linspace(0, 90, n_points - 1)  # DD のビンの端（横棒 = ビンの幅）
    use = ok & dd
    bin_index = np.flatnonzero(use) - 1
    xerr = np.clip([psi_eff[use] - edges[bin_index], edges[bin_index + 1] - psi_eff[use]], 0, None)
    ax.errorbar(psi_eff[use], ratio[use], yerr=error[use], xerr=xerr, fmt="o", color="black", ms=5, label="DD (ψ bins; bar = bin width)")
    for k, name in edge_points.items():
        lo, hi = result["median_psi_range_deg"][name]
        ax.errorbar([psi_eff[k]], [ratio[k]], yerr=[error[k]], xerr=[[psi_eff[k] - lo], [hi - psi_eff[k]]], fmt="s", ms=7,
                    mfc="white", mec="#7b4fd0", ecolor="#7b4fd0", capsize=3, label="DS / SD (reference, not fitted)" if k == 0 else None)
        ax.annotate(name, (psi_eff[k], ratio[k]), textcoords="offset points", xytext=(6, 8), fontsize=9, color="#7b4fd0")
    xx = np.linspace(0, 90, 181)
    ax.plot(xx, fit["A"] * (1 + fit["mu"] * np.cos(2 * np.radians(xx))), color="#e34948",
            label=f"fit to DD only, A[1 + μ cos 2ψ]: μ = {fit['mu']:+.3f} ± {fit['mu_err']:.3f} ({fit['significance_sigma']:.1f}σ vs flat)\n"
                  f"χ²/ndf = {fit['chi2']:.1f}/{fit['ndf']} (flat: {fit['chi2_flat']:.1f}/{fit['ndf_flat']})")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set(xlabel="ψ [deg]  (0: e+e- plane ∥ X, 90: ∥ Y)", ylabel="ratio (each normalized to all pair candidates)",
           title=f"{labels[0]} / {labels[1]}", xlim=(-3, 93))
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.text(0.01, 0.02, "DD: bar = bin width, used in the fit.  DS / SD: reference only; point at the effective ψ, bar = median allowed ψ range",
            transform=ax.transAxes, fontsize=7, color="0.35")
    fig.suptitle(title)
    fig.savefig(target, dpi=120)
    plt.close(fig)


def plot_overview(pipeline, cut_name, settings, bins, target):
    """複数のサンプル（例: 水平・垂直・円偏光）を1枚に重ねる。

    左: 各サンプルの DS・DD（ψ のビン）・SD の割合（対生成の候補の総数に対する %）。
    右: 基準のサンプル（例: 円偏光）に対する比。DD の点だけで A [1 + μ c_k] を fit した線を引き、DS・SD は参考の点（fit に使わない）。
    """
    names = settings["samples"]
    reference = settings["reference"]
    data = {n: load_classified(pipeline, cut_name, n) for n in names}
    points = {n: pair_points(data[n], bins) for n in names}
    counts = {n: np.bincount(points[n][0], minlength=bins + 2) for n in names}
    cos_sum = sum(np.bincount(points[n][0], weights=points[n][1], minlength=bins + 2) for n in names)
    c_k = cos_sum / np.maximum(sum(counts.values()), 1)            # 全サンプルを合わせた点ごとの cos 2ψ の平均
    psi_eff = np.degrees(0.5 * np.arccos(np.clip(c_k, -1, 1)))     # 図の横軸
    edges = np.linspace(0, 90, bins + 1)
    dd = np.zeros(bins + 2, bool)
    dd[1:bins + 1] = True
    palette = ["#2a78d6", "#eb6834", "#1baf7a", "#7b4fd0"]
    others = [n for n in names if n != reference]
    color = {n: palette[i] for i, n in enumerate(others)} | {reference: "#7a7a7a"}
    label = {n: pipeline.samples[n]["label"] for n in names}

    fig, axes = plt.subplots(1, 2, figsize=(16, 5.6), constrained_layout=True)
    ax = axes[0]
    for i, n in enumerate(names):
        share = 100 * counts[n] / counts[n].sum()
        error = 100 * np.sqrt(counts[n]) / counts[n].sum()
        x = psi_eff + (i - (len(names) - 1) / 2) * 0.9  # 重ならないように少しずらす
        ax.errorbar(x[dd], share[dd], yerr=error[dd], fmt="o", ms=5, color=color[n], label=f"{label[n]}: {int(counts[n].sum()):,} pair candidates")
        ax.errorbar(x[~dd], share[~dd], yerr=error[~dd], fmt="s", ms=7, color=color[n])
    for k, text in ((0, "DS"), (bins + 1, "SD")):
        ax.annotate(text, (psi_eff[k], max(100 * counts[n][k] / counts[n].sum() for n in names)), textcoords="offset points",
                    xytext=(8, 4), fontsize=10)
    ax.set(xlabel="ψ [deg]  (0: e+e- plane ∥ X, 90: ∥ Y;  DS, SD: effective ψ)", ylabel="% of pair candidates",
           title="Pair candidates per point (circle: DD ψ bins, square: DS / SD)", xlim=(-3, 93))
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    xx = np.linspace(0, 90, 181)
    fits = {}
    for i, n in enumerate(others):
        ratio, error, ok = normalized_ratio(counts[n], counts[reference])
        use = ok & dd
        fit = fit_cos2_mean(c_k[use], ratio[use], error[use])
        fits[n] = fit
        x = psi_eff + (i - (len(others) - 1) / 2) * 0.9
        bin_index = np.flatnonzero(use) - 1
        xerr = np.clip([x[use] - edges[bin_index], edges[bin_index + 1] - x[use]], 0, None)  # ずらした点でも負にしない
        ax.errorbar(x[use], ratio[use], yerr=error[use], xerr=xerr, fmt="o", ms=5, color=color[n],
                    label=f"{label[n]} / {label[reference]}: DD fit μ = {fit['mu']:+.3f} ± {fit['mu_err']:.3f} ({fit['significance_sigma']:.1f}σ)")
        ax.errorbar(x[~dd], ratio[~dd], yerr=error[~dd], fmt="s", ms=7, mfc="white", mec=color[n], ecolor=color[n])
        ax.plot(xx, fit["A"] * (1 + fit["mu"] * np.cos(2 * np.radians(xx))), color=color[n], lw=1.5, alpha=0.8)
    for k, text in ((0, "DS"), (bins + 1, "SD")):
        ax.annotate(text, (psi_eff[k], 1.0), textcoords="offset points", xytext=(-6, -16), fontsize=10, color="0.3")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set(xlabel="ψ [deg]", ylabel=f"ratio to {label[reference]} (each normalized to all pair candidates)",
           title=f"Modulation relative to {label[reference]} (line: fit to DD points only)", xlim=(-3, 93))
    ax.legend(fontsize=8, frameon=False)
    ax.text(0.01, 0.02, "filled circles: DD ψ bins (bar = bin width), used in the fit;  open squares: DS / SD (ψ only bounded), not fitted",
            transform=ax.transAxes, fontsize=7, color="0.35")
    fig.suptitle(f"[{cut_name}] pair candidates: " + ", ".join(label[n] for n in names))
    fig.savefig(target, dpi=120)
    plt.close(fig)
    return {"samples": names, "reference": reference, "counts": {n: counts[n].tolist() for n in names},
            "mean_cos2": c_k.tolist(), "dd_fit": fits}


# ---------------------------------------------------------------- 組み方を決めた DD の φ（0〜180°）
def pair_phi_result(data, labels, bins, edges=None):
    """組み方を決めた DD の面の方位角 φ_pair の比を A [1 + μ cos 2(φ − φ0)] で fit する（φ0 も自由）。

    edges: ビンの境目 [deg]（増える順で、最後 − 最初 = 180）。None なら 0〜180° を bins 等分。
    φ_pair は 180° 周期なので、[edges[0], edges[0] + 180) に移してから数える。幅の違うビンは、ビンの中で平均したモデルで fit する。
    """
    edges = np.linspace(0, 180, bins + 1) if edges is None else np.asarray(edges, float)
    if not np.isclose(edges[-1] - edges[0], 180) or np.any(np.diff(edges) <= 0):
        raise ValueError(f"phi_edges_deg must increase and span 180 deg: {edges.tolist()}")
    centres = 0.5 * (edges[1:] + edges[:-1])
    widths = np.diff(edges)
    hists = []
    for d in data:
        phi = d["phi_pair"][np.isfinite(d["phi_pair"])]
        hists.append(np.histogram((phi - edges[0]) % 180 + edges[0], bins=edges)[0])
    ratio, error, ok = normalized_ratio(hists[0], hists[1])
    fit, model = fit_cos2_free_phase(centres[ok], ratio[ok], error[ok], widths_deg=widths[ok])
    decided = [{"DD": int(np.sum(d["cls"] == DD)), "decided": int(np.isfinite(d["phi_pair"]).sum())} for d in data]
    result = {"labels": labels, "phi_edges_deg": edges.tolist(), "counts": [h.tolist() for h in hists], "events": decided, "fit": fit}
    return result, (centres, hists, ratio, error, ok, model)


def plot_pair_phi(result, plot_data, title, target):
    centres, hists, ratio, error, ok, model = plot_data
    labels, fit = result["labels"], result["fit"]
    edges = np.array(result["phi_edges_deg"])
    half = 0.5 * np.diff(edges)
    lo, hi = edges[0], edges[-1]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.2), constrained_layout=True)
    ax = axes[0]
    for h, label, color, ev in zip(hists, labels, COLORS, result["events"]):
        n = h.sum()
        density = h / n / (2 * half) * 180  # 平均 = 1（幅の違うビンでも、度あたりの割合で比べる）
        ax.errorbar(centres, density, yerr=density / np.sqrt(np.maximum(h, 1)), xerr=half, fmt="o", ms=4, color=color,
                    label=f"{label}: {ev['decided']:,} of {ev['DD']:,} DD paired")
    ax.set(xlabel="φ_pair [deg]  (pair plane; 0: ∥ X, 90: ∥ Y)", ylabel="normalized (mean = 1)",
           title="DD with resolved XZ/YZ pairing", xlim=(lo, hi))
    ax.legend(fontsize=8, frameon=False)
    ax = axes[1]
    ax.errorbar(centres[ok], ratio[ok], yerr=error[ok], xerr=half[ok], fmt="o", color="black", ms=4)
    xx = np.linspace(lo, hi, 361)
    ax.plot(xx, model(xx, fit["A"], fit["mu"], fit["phi0_deg"]), color="#e34948",
            label=f"μ = {fit['mu']:.3f} ± {fit['mu_err']:.3f}, φ0 = {fit['phi0_deg']:.0f} ± {fit['phi0_err_deg']:.0f}° "
                  f"({fit['significance_sigma']:.1f}σ vs flat)\nχ²/ndf = {fit['chi2']:.1f}/{fit['ndf']}")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set(xlabel="φ_pair [deg]", ylabel="ratio (each normalized)", title=f"{labels[0]} / {labels[1]}: A[1 + μ cos 2(φ − φ0)]", xlim=(lo, hi))
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(title)
    fig.savefig(target, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- 1本の飛跡
def single_track_result(data, labels, theta_mins, bins):
    edges = np.linspace(0, 360, bins + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    rows, plots = {}, []
    for theta_min in theta_mins:
        hists = [np.histogram(d["phi"][(d["cls"] == SS) & (d["theta"] >= theta_min)], bins=edges)[0] for d in data]
        ratio, error, ok = normalized_ratio(hists[0], hists[1])
        fit, model = fit_cos2_free_phase(centres[ok], ratio[ok], error[ok])
        rows[f"{theta_min:g}"] = {"events": [int(h.sum()) for h in hists], "fit": fit, "counts": [h.tolist() for h in hists]}
        plots.append((theta_min, hists, ratio, error, ok, model, fit))
    return {"labels": labels, "phi_edges_deg": edges.tolist(), "theta_min": rows}, (centres, plots)


def plot_single_track(result, plot_data, title, target):
    centres, plots = plot_data
    labels = result["labels"]
    bins = len(centres)
    fig, axes = plt.subplots(len(plots), 2, figsize=(13, 4 * len(plots)), squeeze=False, constrained_layout=True)
    xx = np.linspace(0, 360, 361)
    for row, (theta_min, hists, ratio, error, ok, model, fit) in enumerate(plots):
        ax = axes[row, 0]
        for h, label, color in zip(hists, labels, COLORS):
            n = h.sum()
            ax.errorbar(centres, h / n * bins, yerr=np.sqrt(h) / n * bins, fmt="o-", ms=3, color=color, label=f"{label}: {int(n):,}")
        ax.set(xlabel="φ [deg] (0 = detector X)", ylabel="normalized (mean = 1)", title=f"SS, θ >= {theta_min:g}°", xlim=(0, 360))
        ax.legend(fontsize=8, frameon=False)
        ax = axes[row, 1]
        ax.errorbar(centres[ok], ratio[ok], yerr=error[ok], xerr=180.0 / bins, fmt="o", color="black", ms=4)  # 横棒 = φ のビンの幅
        ax.plot(xx, model(xx, fit["A"], fit["mu"], fit["phi0_deg"]), color="#e34948",
                label=f"μ = {fit['mu']:.3f} ± {fit['mu_err']:.3f}, φ0 = {fit['phi0_deg']:.0f} ± {fit['phi0_err_deg']:.0f}° "
                      f"({fit['significance_sigma']:.1f}σ)\nχ²/ndf = {fit['chi2']:.1f}/{fit['ndf']}")
        ax.axhline(1, color="0.6", lw=0.8)
        ax.set(xlabel="φ [deg]", ylabel="ratio (each normalized)", title=f"{labels[0]} / {labels[1]}: A[1 + μ cos 2(φ − φ0)]", xlim=(0, 360))
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle(title)
    fig.savefig(target, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------- まとめ
def run(pipeline, cut_set_names=None):
    print("== stage4: modulation", flush=True)
    pair_settings = pipeline.analysis["pair"]
    single_settings = pipeline.analysis["single_track"]
    for cut_name in pipeline.selected_cut_sets(cut_set_names):
        summary = []
        for comparison in pipeline.comparisons:
            names = (comparison["numerator"], comparison["denominator"])
            labels = [pipeline.samples[n]["label"] for n in names]
            data = [load_classified(pipeline, cut_name, n) for n in names]
            folder = pipeline.stage_dir(4, cut_name, comparison["name"])
            title = f"[{cut_name}] {comparison['name']}: {labels[0]} / {labels[1]}"
            entry = {"comparison": comparison["name"], "numerator": labels[0], "denominator": labels[1],
                     "note": comparison.get("note", "")}
            if pair_settings.get("enabled", True):
                result, plot_data = pair_result(data, labels, pair_settings["psi_bins"])
                (folder / "pair.json").write_text(json.dumps(result, indent=1) + "\n")
                plot_pair(result, plot_data, title + "  (pair candidates)", folder / "pair.png")
                curve, curve_plot = pair_curve_result(data, labels, pair_settings["psi_bins"])
                (folder / "pair_curve.json").write_text(json.dumps(curve, indent=1) + "\n")
                plot_pair_curve(curve, curve_plot, title + "  (pair candidates: DS + DD + SD)", folder / "pair_curve.png")
                entry["pair"] = {"events": [a["N_X"] + a["N_Y"] for a in result["asymmetry"]],
                                 "all_points_mu": curve["fit_all_points_reference"]["mu"],
                                 "all_points_mu_err": curve["fit_all_points_reference"]["mu_err"],
                                 "all_points_chi2": curve["fit_all_points_reference"]["chi2"],
                                 "all_points_ndf": curve["fit_all_points_reference"]["ndf"],
                                 "Q": [a["Q"] for a in result["asymmetry"]], "delta_Q": result["delta_Q"],
                                 "dd_mu": result["dd_psi_fit"]["mu"], "dd_mu_err": result["dd_psi_fit"]["mu_err"],
                                 "dd_sigma": result["dd_psi_fit"]["significance_sigma"]}
            paired = [int(np.isfinite(d["phi_pair"]).sum()) if "phi_pair" in d else 0 for d in data]
            if pipeline.pairing.get("enabled", False) and min(paired) >= 50:  # 組み方の材料が無い（足りない）サンプルは飛ばす
                result, plot_data = pair_phi_result(data, labels, pair_settings.get("phi_bins", 12), pair_settings.get("phi_edges_deg"))
                (folder / "pair_phi.json").write_text(json.dumps(result, indent=1) + "\n")
                plot_pair_phi(result, plot_data, title + "  (DD, resolved pairing)", folder / "pair_phi.png")
                f = result["fit"]
                entry["pair_phi"] = {"events": [e["decided"] for e in result["events"]], "mu": f["mu"], "mu_err": f["mu_err"],
                                     "phi0_deg": f["phi0_deg"], "phi0_err_deg": f["phi0_err_deg"], "sigma": f["significance_sigma"]}
            if single_settings.get("enabled", True):
                result, plot_data = single_track_result(data, labels, single_settings["theta_min_deg"], single_settings["phi_bins"])
                (folder / "single_track.json").write_text(json.dumps(result, indent=1) + "\n")
                plot_single_track(result, plot_data, title + "  (single tracks, SS)", folder / "single_track.png")
                first = result["theta_min"][f"{single_settings['theta_min_deg'][0]:g}"]
                entry["single_track"] = {"theta_min_deg": single_settings["theta_min_deg"][0], "events": first["events"],
                                         "mu": first["fit"]["mu"], "mu_err": first["fit"]["mu_err"],
                                         "phi0_deg": first["fit"]["phi0_deg"], "sigma": first["fit"]["significance_sigma"]}
            summary.append(entry)
            print_entry(cut_name, entry)
        folder = pipeline.stage_dir(4, cut_name)
        overview = pipeline.analysis.get("overview", {"enabled": False})
        if overview.get("enabled", False):
            result = plot_overview(pipeline, cut_name, overview, pair_settings["psi_bins"], folder / "overview.png")
            (folder / "overview.json").write_text(json.dumps(result, indent=1) + "\n")
            print(f"  saved {folder / 'overview.png'}", flush=True)
        (folder / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
        (folder / "summary.md").write_text(summary_markdown(cut_name, pipeline.cut_sets[cut_name], summary))
        print(f"  saved {folder / 'summary.md'}", flush=True)


def print_entry(cut_name, entry):
    text = f"  [{cut_name}] {entry['comparison']}:"
    if "pair" in entry:
        p = entry["pair"]
        text += (f" pair ΔQ = {p['delta_Q']['value']:+.4f} ± {p['delta_Q']['err']:.4f} ({p['delta_Q']['significance_sigma']:.1f}σ),"
                 f" DD μ = {p['dd_mu']:+.3f} ± {p['dd_mu_err']:.3f} ({p['dd_sigma']:.1f}σ),"
                 f" (reference DS+DD+SD μ = {p['all_points_mu']:+.3f}, χ²/ndf = {p['all_points_chi2']:.1f}/{p['all_points_ndf']})")
    if "pair_phi" in entry:
        q = entry["pair_phi"]
        text += f" | DD φ: μ = {q['mu']:.3f} ± {q['mu_err']:.3f} at φ0 = {q['phi0_deg']:.0f} ± {q['phi0_err_deg']:.0f}° ({q['sigma']:.1f}σ)"
    if "single_track" in entry:
        s = entry["single_track"]
        text += f" | SS θ>={s['theta_min_deg']:g}: μ = {s['mu']:.3f} ± {s['mu_err']:.3f} at {s['phi0_deg']:.0f}° ({s['sigma']:.1f}σ)"
    print(text, flush=True)


def summary_markdown(cut_name, cuts, summary):
    lines = [f"# Modulation summary: {cut_name}", "", "```json", json.dumps(cuts, indent=1), "```", "",
             "| comparison | pair events | Q (num / den) | ΔQ | DD μ (cos 2ψ) | reference: DS+DD+SD μ (χ²/ndf) | DD φ (resolved pairing): μ, φ0 | single-track μ, φ0 |", "|---|---|---|---|---|---|---|---|"]
    for e in summary:
        p, s = e.get("pair"), e.get("single_track")
        pair_cells = ["-"] * 5
        if p:
            d = p["delta_Q"]
            pair_cells = [f"{p['events'][0]:,} / {p['events'][1]:,}", f"{p['Q'][0]:+.3f} / {p['Q'][1]:+.3f}",
                          f"{d['value']:+.3f} ± {d['err']:.3f} ({d['significance_sigma']:.1f}σ)",
                          f"{p['dd_mu']:+.3f} ± {p['dd_mu_err']:.3f} ({p['dd_sigma']:.1f}σ)",
                          f"{p['all_points_mu']:+.3f} ± {p['all_points_mu_err']:.3f} ({p['all_points_chi2']:.1f}/{p['all_points_ndf']})"]
        q = e.get("pair_phi")
        phi_cell = (f"{q['mu']:.3f} ± {q['mu_err']:.3f}, {q['phi0_deg']:.0f} ± {q['phi0_err_deg']:.0f}° ({q['sigma']:.1f}σ)" if q else "-")
        single_cell = (f"θ≥{s['theta_min_deg']:g}°: {s['mu']:.3f} ± {s['mu_err']:.3f}, {s['phi0_deg']:.0f}° ({s['sigma']:.1f}σ)"
                       if s else "-")
        lines.append(f"| {e['comparison']} ({e['numerator']} / {e['denominator']}) | " + " | ".join(pair_cells) + f" | {phi_cell} | {single_cell} |")
    lines += ["", "σ: statistical only. DD μ: fit to DD ψ bins, Δχ² vs flat with 1 dof (main result). Reference DS+DD+SD μ assumes ψ of DS/SD uniform in its allowed range and does not fit well (DS/SD show smaller modulation). Single-track uses 2 dof (μ, φ0)."]
    return "\n".join(lines) + "\n"
