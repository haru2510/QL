"""stage4: 比較（分子のサンプル ÷ 分母のサンプル）ごとに、モジュレーションと有意度を出す。

対生成（analysis.pair）:
  - Q = (N_X − N_Y)/(N_X + N_Y)、N_X = DS + DD(ψ < 45°)、N_Y = SD + DD(ψ >= 45°)。サンプルごとに出して差 ΔQ を取る。
  - DD の ψ（0〜90°）のヒストグラムを各サンプルで規格化して割り、A [1 + μ cos 2ψ] で fit する。
    折り返しているので、偏光が X・Y に対して φ0 傾いていると見える μ は μ_true · cos 2φ0（45° では 0）。
  - DS・DD・SD をまとめたカーブ: DS（ψ ≈ 0° 側）、DD の ψ の各ビン、SD（ψ ≈ 90° 側）の 11 点を、対生成の候補の総数で
    規格化して割る。DS・SD は ψ が範囲 [psi_lo, psi_hi] でしか決まらないので、範囲内で ψ が一様と仮定して
    cos 2ψ を平均した c_k を点ごとに出し、A [1 + μ c_k] で fit する（DD のビンは c_k = ビン内の cos 2ψ の平均）。
  - 有意度: 平らな線との Δχ²（1自由度）と ΔQ / 誤差。
1本の飛跡（analysis.single_track、SS のみ）:
  - θ >= θmin ごとに φ（0〜360°）のヒストグラムの比を A [1 + μ cos 2(φ − φ0)] で fit（Δχ² は2自由度）。

アクセプタンス（検出器の形）は割り算で消す。分母には偏光の無い（円偏光の）データを使うのが基本。
出力: stage4_modulation/<cut_set>/<comparison>/{pair.png, pair.json, pair_curve.png, pair_curve.json, single_track.png, single_track.json}, summary.md, summary.json
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
        return {key: saved[key] for key in ("cls", "psi", "psi_lo", "psi_hi", "theta", "phi")}


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
    ax.errorbar(centres[ok], ratio[ok], yerr=error[ok], fmt="o", color="black", ms=4)
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


def pair_curve_result(data, labels, bins):
    """DS・DD（ψ のビン）・SD の点を、対生成の候補の総数で規格化して割り、A [1 + μ c_k] で fit する。"""
    edges = np.linspace(0, 90, bins + 1)
    names = ["DS"] + [f"DD {edges[k]:.0f}-{edges[k + 1]:.0f}" for k in range(bins)] + ["SD"]
    counts, cos_sums, ranges = [], np.zeros(bins + 2), [[] for _ in range(bins + 2)]
    for d in data:
        cls, lo, hi = d["cls"], d["psi_lo"], d["psi_hi"]
        point = np.full(len(cls), -1)
        point[cls == DS] = 0
        dd = cls == DD
        point[dd] = 1 + np.clip(np.digitize(d["psi"][dd], edges) - 1, 0, bins - 1)
        point[cls == SD] = bins + 1
        used = point >= 0
        counts.append(np.bincount(point[used], minlength=bins + 2))
        c = mean_cos2(lo[used], hi[used])
        cos_sums += np.bincount(point[used], weights=c, minlength=bins + 2)
        for k in (0, bins + 1):  # DS・SD の ψ の範囲（図の横棒用）
            ranges[k].append(np.column_stack([lo[point == k], hi[point == k]]))
    total = counts[0] + counts[1]
    c_k = cos_sums / np.maximum(total, 1)  # 2つのサンプルを合わせた点ごとの cos 2ψ の平均
    ratio, error, ok = normalized_ratio(counts[0], counts[1])
    fit = fit_cos2_mean(c_k[ok], ratio[ok], error[ok])
    psi_eff = np.degrees(0.5 * np.arccos(np.clip(c_k, -1, 1)))  # cos 2ψ_eff = c_k となる ψ（図の横軸）
    spans = {}
    for k, name in ((0, "DS"), (bins + 1, "SD")):
        both = np.concatenate(ranges[k])
        spans[name] = [float(np.median(both[:, 0])), float(np.median(both[:, 1]))]
    result = {"labels": labels, "points": names, "counts": [c.tolist() for c in counts], "mean_cos2": c_k.tolist(),
              "psi_effective_deg": psi_eff.tolist(), "ratio": ratio.tolist(), "ratio_err": error.tolist(),
              "median_psi_range_deg": spans, "fit": fit}
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
    ax.errorbar(psi_eff[ok & dd], ratio[ok & dd], yerr=error[ok & dd], fmt="o", color="black", ms=5, label="DD (ψ bins)")
    for k, name in edge_points.items():
        lo, hi = result["median_psi_range_deg"][name]
        ax.errorbar([psi_eff[k]], [ratio[k]], yerr=[error[k]], xerr=[[psi_eff[k] - lo], [hi - psi_eff[k]]], fmt="s", ms=7,
                    color="#7b4fd0", capsize=3)
        ax.annotate(name, (psi_eff[k], ratio[k]), textcoords="offset points", xytext=(6, 8), fontsize=9, color="#7b4fd0")
    xx = np.linspace(0, 90, 181)
    ax.plot(xx, fit["A"] * (1 + fit["mu"] * np.cos(2 * np.radians(xx))), color="#e34948",
            label=f"A[1 + μ cos 2ψ]: μ = {fit['mu']:+.3f} ± {fit['mu_err']:.3f} ({fit['significance_sigma']:.1f}σ vs flat)\n"
                  f"χ²/ndf = {fit['chi2']:.1f}/{fit['ndf']} (flat: {fit['chi2_flat']:.1f}/{fit['ndf_flat']})")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set(xlabel="ψ [deg]  (0: e+e- plane ∥ X, 90: ∥ Y)", ylabel="ratio (each normalized to all pair candidates)",
           title=f"{labels[0]} / {labels[1]}", xlim=(-3, 93))
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.text(0.01, 0.02, "DS / SD: point at the effective ψ, horizontal bar = median allowed ψ range (ψ assumed uniform in it)",
            transform=ax.transAxes, fontsize=7, color="0.35")
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
        ax.errorbar(centres[ok], ratio[ok], yerr=error[ok], fmt="o", color="black", ms=4)
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
                                 "curve_mu": curve["fit"]["mu"], "curve_mu_err": curve["fit"]["mu_err"],
                                 "curve_sigma": curve["fit"]["significance_sigma"],
                                 "Q": [a["Q"] for a in result["asymmetry"]], "delta_Q": result["delta_Q"],
                                 "dd_mu": result["dd_psi_fit"]["mu"], "dd_mu_err": result["dd_psi_fit"]["mu_err"],
                                 "dd_sigma": result["dd_psi_fit"]["significance_sigma"]}
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
        (folder / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
        (folder / "summary.md").write_text(summary_markdown(cut_name, pipeline.cut_sets[cut_name], summary))
        print(f"  saved {folder / 'summary.md'}", flush=True)


def print_entry(cut_name, entry):
    text = f"  [{cut_name}] {entry['comparison']}:"
    if "pair" in entry:
        p = entry["pair"]
        text += (f" pair ΔQ = {p['delta_Q']['value']:+.4f} ± {p['delta_Q']['err']:.4f} ({p['delta_Q']['significance_sigma']:.1f}σ),"
                 f" DD μ = {p['dd_mu']:+.3f} ± {p['dd_mu_err']:.3f} ({p['dd_sigma']:.1f}σ),"
                 f" DS+DD+SD μ = {p['curve_mu']:+.3f} ± {p['curve_mu_err']:.3f} ({p['curve_sigma']:.1f}σ)")
    if "single_track" in entry:
        s = entry["single_track"]
        text += f" | SS θ>={s['theta_min_deg']:g}: μ = {s['mu']:.3f} ± {s['mu_err']:.3f} at {s['phi0_deg']:.0f}° ({s['sigma']:.1f}σ)"
    print(text, flush=True)


def summary_markdown(cut_name, cuts, summary):
    lines = [f"# Modulation summary: {cut_name}", "", "```json", json.dumps(cuts, indent=1), "```", "",
             "| comparison | pair events | Q (num / den) | ΔQ | DS+DD+SD curve μ | DD-only μ | single-track μ, φ0 |", "|---|---|---|---|---|---|---|"]
    for e in summary:
        p, s = e.get("pair"), e.get("single_track")
        pair_cells = ["-"] * 5
        if p:
            d = p["delta_Q"]
            pair_cells = [f"{p['events'][0]:,} / {p['events'][1]:,}", f"{p['Q'][0]:+.3f} / {p['Q'][1]:+.3f}",
                          f"{d['value']:+.3f} ± {d['err']:.3f} ({d['significance_sigma']:.1f}σ)",
                          f"{p['curve_mu']:+.3f} ± {p['curve_mu_err']:.3f} ({p['curve_sigma']:.1f}σ)",
                          f"{p['dd_mu']:+.3f} ± {p['dd_mu_err']:.3f} ({p['dd_sigma']:.1f}σ)"]
        single_cell = (f"θ≥{s['theta_min_deg']:g}°: {s['mu']:.3f} ± {s['mu_err']:.3f}, {s['phi0_deg']:.0f}° ({s['sigma']:.1f}σ)"
                       if s else "-")
        lines.append(f"| {e['comparison']} ({e['numerator']} / {e['denominator']}) | " + " | ".join(pair_cells) + f" | {single_cell} |")
    lines += ["", "σ: statistical only. Curve μ (DS + DD ψ bins + SD, ψ of DS/SD assumed uniform in its allowed range) and DD-only μ use Δχ² vs flat with 1 dof; single-track uses 2 dof (μ, φ0)."]
    return "\n".join(lines) + "\n"
