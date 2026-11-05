#!/usr/bin/env python3
#
# Cs-137 と Co-60 のスペクトルをガウスで fit し、ch ごとに ADC (hg) からエネルギーへの換算係数を求める。
#
#   Cs : 662 keV の光電ピークを ガウス 1 本 + 直線の背景 で fit
#   Co : 1173 keV と 1332 keV を ガウス 2 本 + 直線の背景 で fit (2 本の中心は別々に動かす)
#   換算 : 3 点に hg = a * E + b を当てはめる (b はオフセット。0 に固定しない)
#
# ゲインは温度で 1 ℃ あたり約 8 % 変わるので、両方の線源で同じ温度の区間 (既定 2 ℃ 幅) だけを使う。
# 区間は HA ごとに「2 つの線源のヒット数の少ない方」がいちばん多くなるものを選ぶ (--temp-low で固定もできる)。
# ヒットはトリガの時刻に当たる peak_id のものだけを使い、飽和したヒットは除く (temp_spectra.py と同じ)。
#
#   出力: <out-dir>/<HA>_energy_calib.png   スペクトルと fit、換算の直線
#         <out-dir>/energy_calib.csv        HA, ch ごとの a, b とピークの位置・幅 (HA ごとに追記・上書き)
#
# 例: python3 energy_calib.py --ha ha02
#     python3 energy_calib.py --ha ha00 ha01 ha02 --cs-dir .../h5_runs --co-dir .../h5_runs_run12
#

import argparse
import csv
import glob
import os
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from temp_spectra import PEAK_ID_MIN, on_time_range  # noqa: E402

# このリポジトリ (analysis_sm3_FM) の場所。既定のパスはここからの相対で決める
REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


DEFAULT_CS_DIR = os.path.join(REPO_DIR, "PSA_analysis", "h5_runs")
DEFAULT_CO_DIR = os.path.join(REPO_DIR, "PSA_analysis", "h5_runs_run12")
DEFAULT_OUTPUT_DIR = os.path.join(REPO_DIR, "check_output", "energy_calib")

N_CHANNELS = 6

E_CS = 661.66
E_CO1 = 1173.23
E_CO2 = 1332.49

HG_MAX = 120e3
HG_BINS = 480

# 換算係数を使ってよいかの判定 (quality 列)。どれかに当たれば "bad:<理由>"
CS_MIN_OK = 24e3            # Cs のピークがこれより下だと hit threshold の段差 (約 2 万) に重なり、位置がずれる
CO_RATIO_OK = (1.55, 1.95)  # 1173 keV / 662 keV のピークの比。外れたら Co の fit が別の山を拾っている
A_RELERR_OK = 0.03          # a の相対誤差


def parse_args():

    parser = argparse.ArgumentParser(description="ADC -> energy calibration from Cs-137 and Co-60 peaks")
    parser.add_argument("--ha", nargs="+", required=True, help="HA to calibrate, e.g. ha02 ha03")
    parser.add_argument("--cs-dir", default=DEFAULT_CS_DIR, help="parent of HA directories of the Cs runs")
    parser.add_argument("--co-dir", default=DEFAULT_CO_DIR, help="parent of HA directories of the Co runs")
    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--temp-width", type=float, default=2.0, help="width of the temperature range [degC]")
    parser.add_argument("--temp-low", type=float, default=None,
                        help="lower edge of the temperature range (default: chosen per HA)")
    parser.add_argument("--cs-min", type=float, default=15e3,
                        help="search the Cs photopeak above this [ADC x sample]")
    parser.add_argument("--linear", action="store_true", help="draw the spectra in linear scale (default: log)")

    return parser.parse_args()


def load(directory):
    """1 つの HA の .h5 を読み、トリガの時刻に当たる飽和していないヒットの hg, ch, 温度 を返す。"""

    parts = {"hg": [], "ch": [], "temp": [], "peak_id": []}

    for path in natsorted(glob.glob(os.path.join(directory, "*.h5"))):

        with h5py.File(path, "r") as h5:

            selected = (h5["peak_id"][:] > PEAK_ID_MIN) & (h5["satu_h"][:] == 0) & (h5["satu_l"][:] == 0)

            parts["hg"].append(h5["hg"][:][selected])
            parts["ch"].append(h5["ch"][:][selected])
            parts["temp"].append(h5["temp"][:][selected].mean(axis=1))
            parts["peak_id"].append(h5["peak_id"][:][selected])

    if not parts["hg"]:
        raise FileNotFoundError(f"No .h5 files in {directory}")

    data = {name: np.concatenate(value) for name, value in parts.items()}
    low, high = on_time_range(data["peak_id"], None)
    on_time = (data["peak_id"] >= low) & (data["peak_id"] <= high)

    return {name: value[on_time] for name, value in data.items()}


def choose_temp_range(cs, co, width):
    """2 つの線源のヒット数の少ない方が最大になる温度の区間 [low, low + width) を返す。"""

    lows = np.arange(np.floor(min(cs["temp"].min(), co["temp"].min()) / width) * width,
                     max(cs["temp"].max(), co["temp"].max()), width)
    counts = [min(((d["temp"] >= t) & (d["temp"] < t + width)).sum() for d in (cs, co)) for t in lows]

    return float(lows[int(np.argmax(counts))])


def gauss(x, amplitude, mean, sigma):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def model_cs(x, amplitude, mean, sigma, slope, offset):
    return gauss(x, amplitude, mean, sigma) + slope * (x - mean) + offset


def model_co(x, amplitude1, mean1, sigma1, amplitude_ratio, mean_ratio, slope, offset):
    """1173 keV と 1332 keV の 2 本。1332 keV の中心は mean1 * mean_ratio、幅は sigma1 * sqrt(mean_ratio)。"""
    mean2, sigma2 = mean1 * mean_ratio, sigma1 * np.sqrt(mean_ratio)
    return (gauss(x, amplitude1, mean1, sigma1) + gauss(x, amplitude1 * amplitude_ratio, mean2, sigma2)
            + slope * (x - mean1) + offset)


def co_lines(p, cov):
    """model_co のパラメータから 2 本の (高さ, 中心, 幅) と、中心・幅の誤差を返す。

    1332 keV の中心 mean1 * r と幅 sigma1 * sqrt(r) は 2 つのパラメータの関数なので、共分散も入れて誤差を伝える。
    戻り値: (line1, line2), (err_mean1, err_sigma1), (err_mean2, err_sigma2)
    """
    amplitude1, mean1, sigma1, amplitude_ratio, mean_ratio = p[:5]
    mean2, sigma2 = mean1 * mean_ratio, sigma1 * np.sqrt(mean_ratio)

    # d(mean2)/d(mean1, r), d(sigma2)/d(sigma1, r)
    j_mean2 = np.zeros(len(p)); j_mean2[1], j_mean2[4] = mean_ratio, mean1
    j_sigma2 = np.zeros(len(p)); j_sigma2[2], j_sigma2[4] = np.sqrt(mean_ratio), 0.5 * sigma1 / np.sqrt(mean_ratio)

    err = np.sqrt(np.diag(cov))
    return (((amplitude1, mean1, sigma1), (amplitude1 * amplitude_ratio, mean2, sigma2)),
            (err[1], err[2]),
            (np.sqrt(j_mean2 @ cov @ j_mean2), np.sqrt(j_sigma2 @ cov @ j_sigma2)))


def fit(model, centers, counts, window, p0, bounds):
    """window の中だけを Poisson の誤差で fit する。戻り値は (パラメータ, 誤差, 共分散, chi2, ndf) か None。"""

    inside = (centers >= window[0]) & (centers <= window[1])
    if inside.sum() < len(p0) + 3:
        return None

    try:
        popt, pcov = curve_fit(model, centers[inside], counts[inside], p0=p0, bounds=bounds,
                               sigma=np.sqrt(np.maximum(counts[inside], 1)), absolute_sigma=True, maxfev=20000)
    except (RuntimeError, ValueError):
        return None

    x, y = centers[inside], counts[inside]
    chi2 = float(np.sum((y - model(x, *popt)) ** 2 / np.maximum(y, 1)))
    ndf = int(inside.sum() - len(popt))

    return popt, np.sqrt(np.diag(pcov)), pcov, chi2, ndf


def smoothed_argmax(centers, counts, low, high):

    smooth = np.convolve(counts, np.ones(5) / 5, mode="same")
    inside = (centers >= low) & (centers <= high)

    return centers[inside][np.argmax(smooth[inside])]


def fit_cs(centers, counts, search_min):
    """Cs の光電ピーク。低い側は hit threshold の肩があるので、窓は上側を広く取り、3 回まわして窓を合わせる。"""

    mean = smoothed_argmax(centers, counts, search_min, HG_MAX)
    sigma = 0.08 * mean
    result = None

    for _ in range(3):

        window = (mean - 1.0 * sigma, mean + 2.0 * sigma)
        peak = counts[np.argmin(abs(centers - mean))]
        p0 = [peak, mean, sigma, 0.0, 0.0]
        bounds = ([0, window[0], 0.01 * mean, -np.inf, -np.inf], [np.inf, window[1], 0.5 * mean, np.inf, np.inf])
        out = fit(model_cs, centers, counts, window, p0, bounds)
        if out is None:
            break
        result = out
        mean, sigma = abs(out[0][1]), abs(out[0][2])

    return result


def fit_co(centers, counts, cs_mean):
    """Co の 2 本。Cs のピークの 1.4-2.2 倍の範囲で山を探し、その位置を基準に窓を決めて 1 回 fit する。

    1332/1173 = 1.136 で、オフセット b > 0 ならこれより小さくなるので、中心の比は 1.06-1.16 に制限する。
    2 本の高さの比 (1332 / 1173) は 0.5-1.5 に制限する。窓を fit の結果で広げると、
    1173 keV のガウスが低い側のコンプトンの連続部分に流れるので、窓は山の位置で固定する。
    """

    peak = smoothed_argmax(centers, counts, 1.4 * cs_mean, 2.2 * cs_mean)
    window = (0.85 * peak, 1.30 * peak)
    height = counts[np.argmin(abs(centers - peak))]
    mean1 = 0.96 * peak

    p0 = [0.7 * height, mean1, 0.05 * peak, 0.9, 1.11, 0.0, 0.0]
    bounds = ([0, 0.80 * peak, 0.01 * peak, 0.5, 1.06, -np.inf, -np.inf],
              [np.inf, 1.10 * peak, 0.15 * peak, 1.5, 1.16, np.inf, np.inf])

    return fit(model_co, centers, counts, window, p0, bounds)


def linear_fit(energy, mean, error):
    """hg = a * E + b を重み付きで当てはめる。

    戻り値: a, b, a の誤差, b の誤差, chi2, ndf
    a, b の誤差は、chi2/ndf > 1 のとき sqrt(chi2/ndf) 倍に広げる (点のばらつきが fit の誤差より大きいため)。
    """

    weight = 1.0 / np.asarray(error) ** 2
    design = np.vstack([energy, np.ones_like(energy)]).T
    cov = np.linalg.inv(design.T @ (design * weight[:, None]))
    a, b = cov @ design.T @ (weight * mean)

    chi2 = float(np.sum(weight * (mean - (a * energy + b)) ** 2))
    ndf = len(energy) - 2
    scale = np.sqrt(max(chi2 / ndf, 1.0)) if ndf > 0 else 1.0

    return a, b, scale * np.sqrt(cov[0, 0]), scale * np.sqrt(cov[1, 1]), chi2, ndf


def calibrate_ha(ha, args):

    cs = load(os.path.join(args.cs_dir, ha))
    co = load(os.path.join(args.co_dir, ha))

    temp_low = args.temp_low if args.temp_low is not None else choose_temp_range(cs, co, args.temp_width)
    temp_high = temp_low + args.temp_width
    for d in (cs, co):
        keep = (d["temp"] >= temp_low) & (d["temp"] < temp_high)
        for name in list(d):
            d[name] = d[name][keep]

    print(f"[READ] {ha}: {temp_low:.0f}-{temp_high:.0f} degC, Cs {len(cs['hg'])} hits, Co {len(co['hg'])} hits")

    edges = np.linspace(0, HG_MAX, HG_BINS + 1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    energy = np.array([E_CS, E_CO1, E_CO2])

    fig, axes = plt.subplots(4, N_CHANNELS, figsize=(4.2 * N_CHANNELS, 14),
                             gridspec_kw={"height_ratios": [1, 1, 1, 0.55]})
    rows = []

    for ch in range(N_CHANNELS):

        counts_cs = np.histogram(cs["hg"][cs["ch"] == ch], bins=edges)[0].astype(float)
        counts_co = np.histogram(co["hg"][co["ch"] == ch], bins=edges)[0].astype(float)

        result_cs = fit_cs(centers, counts_cs, args.cs_min) if counts_cs.sum() > 1000 else None
        result_co = fit_co(centers, counts_co, result_cs[0][1]) if (result_cs and counts_co.sum() > 1000) else None

        ax_cs, ax_co, ax_line, ax_res = axes[0, ch], axes[1, ch], axes[2, ch], axes[3, ch]
        ax_cs.step(centers / 1e3, counts_cs, where="mid", color="tab:blue", lw=0.8)
        ax_co.step(centers / 1e3, counts_co, where="mid", color="tab:blue", lw=0.8)
        ax_cs.set_title(f"{ha} ch{ch}  Cs-137 ({int(counts_cs.sum())} hits)", fontsize=9)
        ax_co.set_title(f"{ha} ch{ch}  Co-60 ({int(counts_co.sum())} hits)", fontsize=9)

        if result_cs is None or result_co is None:
            for ax in (ax_cs, ax_co, ax_line, ax_res):
                ax.text(0.5, 0.5, "fit failed", transform=ax.transAxes, ha="center", color="red")
            print(f"  ch{ch}: fit failed")
            continue

        (amp, mean_cs, sig_cs, slope, offset), err_cs, _, chi2_cs, ndf_cs = result_cs
        err_mean_cs, err_sig_cs = err_cs[1], err_cs[2]
        xs = np.linspace(mean_cs - 1.0 * sig_cs, mean_cs + 2.0 * sig_cs, 200)
        ax_cs.plot(xs / 1e3, model_cs(xs, *result_cs[0]), "r-", lw=1.2)
        ax_cs.plot(xs / 1e3, gauss(xs, amp, mean_cs, sig_cs), "r:", lw=1)
        ax_cs.axvline(mean_cs / 1e3, color="r", lw=0.6, alpha=0.6)
        fwhm_cs = 2.355 * sig_cs / mean_cs * 100
        err_fwhm_cs = fwhm_cs * np.hypot(err_sig_cs / sig_cs, err_mean_cs / mean_cs)
        ax_cs.text(0.97, 0.95,
                   f"662 keV\n mean {mean_cs/1e3:.3f} $\\pm$ {err_mean_cs/1e3:.3f}k\n"
                   f" sigma {sig_cs/1e3:.3f} $\\pm$ {err_sig_cs/1e3:.3f}k\n"
                   f" FWHM {fwhm_cs:.1f} $\\pm$ {err_fwhm_cs:.1f} %\n"
                   f" $\\chi^2$/ndf {chi2_cs:.0f}/{ndf_cs}",
                   transform=ax_cs.transAxes, ha="right", va="top", fontsize=7)

        p_co, _, cov_co, chi2_co, ndf_co = result_co
        ((amp1, mean1, sig1), (amp2, mean2, sig2)), (err_mean1, err_sig1), (err_mean2, err_sig2) = \
            co_lines(p_co, cov_co)
        xs = np.linspace(mean1 - 2.0 * sig1, mean2 + 2.0 * sig2, 300)
        ax_co.plot(xs / 1e3, model_co(xs, *p_co), "r-", lw=1.2)
        ax_co.plot(xs / 1e3, gauss(xs, amp1, mean1, sig1), "m:", lw=1)
        ax_co.plot(xs / 1e3, gauss(xs, amp2, mean2, sig2), "g:", lw=1)
        for m in (mean1, mean2):
            ax_co.axvline(m / 1e3, color="r", lw=0.6, alpha=0.6)
        ax_co.text(0.97, 0.95,
                   f"1173 keV\n mean {mean1/1e3:.3f} $\\pm$ {err_mean1/1e3:.3f}k\n"
                   f" sigma {sig1/1e3:.3f} $\\pm$ {err_sig1/1e3:.3f}k\n"
                   f"1332 keV\n mean {mean2/1e3:.3f} $\\pm$ {err_mean2/1e3:.3f}k\n"
                   f" sigma {sig2/1e3:.3f} $\\pm$ {err_sig2/1e3:.3f}k\n"
                   f" $\\chi^2$/ndf {chi2_co:.0f}/{ndf_co}",
                   transform=ax_co.transAxes, ha="right", va="top", fontsize=7)

        means = np.array([mean_cs, mean1, mean2])
        errors = np.array([err_mean_cs, err_mean1, err_mean2])  # ガウス fit の中心の誤差
        a, b, err_a, err_b, chi2_line, ndf_line = linear_fit(energy, means, errors)
        a0 = np.sum(energy * means) / np.sum(energy ** 2)  # b = 0 のとき

        ex = np.linspace(0, 1500, 50)
        ax_line.errorbar(energy, means / 1e3, yerr=errors / 1e3, fmt="o", color="k", ms=4)
        ax_line.plot(ex, (a * ex + b) / 1e3, "r-", lw=1,
                     label=f"a = {a:.3f} $\\pm$ {err_a:.3f} /keV\nb = {b/1e3:+.3f} $\\pm$ {err_b/1e3:.3f}k\n"
                           f"$\\chi^2$/ndf = {chi2_line:.1f}/{ndf_line}")
        ax_line.plot(ex, a0 * ex / 1e3, "b--", lw=0.8, label=f"b = 0: hg = {a0:.2f} E")
        residual = means - (a * energy + b)
        reasons = []
        if mean_cs < CS_MIN_OK:
            reasons.append("cs_near_threshold")
        if not CO_RATIO_OK[0] <= mean1 / mean_cs <= CO_RATIO_OK[1]:
            reasons.append("co_ratio")
        if err_a / a > A_RELERR_OK:
            reasons.append("a_err")
        quality = "ok" if not reasons else "bad:" + "+".join(reasons)
        ax_line.set_title(f"ch{ch} [{quality}]: residual " + ", ".join(f"{r/1e3:+.2f}k" for r in residual),
                          fontsize=8, color="k" if quality == "ok" else "red")
        ax_line.legend(fontsize=7, loc="upper left")

        ax_res.errorbar(energy, residual, yerr=errors, fmt="o", color="k", ms=4, capsize=3)
        ax_res.axhline(0, color="r", lw=0.8)
        ax_res.set_title(f"ch{ch}: data - fit (error = Gaussian fit error of the mean)", fontsize=7)

        print(f"  ch{ch}: Cs {mean_cs/1e3:.2f}k, Co {mean1/1e3:.2f}k {mean2/1e3:.2f}k -> "
              f"a = {a:.3f} +- {err_a:.3f} /keV, b = {b/1e3:+.2f} +- {err_b/1e3:.2f}k  [{quality}]")

        rows.append({
            "HA": ha, "ch": ch, "temp_low": temp_low, "temp_high": temp_high,
            "a_per_keV": f"{a:.6g}", "a_err": f"{err_a:.3g}", "b": f"{b:.6g}", "b_err": f"{err_b:.3g}",
            "chi2_line": f"{chi2_line:.3g}", "ndf_line": ndf_line,
            "a_b0_per_keV": f"{a0:.6g}",
            "mean_662": f"{mean_cs:.6g}", "mean_662_err": f"{err_mean_cs:.3g}",
            "sigma_662": f"{sig_cs:.6g}", "sigma_662_err": f"{err_sig_cs:.3g}",
            "mean_1173": f"{mean1:.6g}", "mean_1173_err": f"{err_mean1:.3g}",
            "sigma_1173": f"{sig1:.6g}", "sigma_1173_err": f"{err_sig1:.3g}",
            "mean_1332": f"{mean2:.6g}", "mean_1332_err": f"{err_mean2:.3g}",
            "sigma_1332": f"{sig2:.6g}", "sigma_1332_err": f"{err_sig2:.3g}",
            "chi2_cs": f"{chi2_cs:.4g}", "ndf_cs": ndf_cs, "chi2_co": f"{chi2_co:.4g}", "ndf_co": ndf_co,
            "hits_cs": int(counts_cs.sum()), "hits_co": int(counts_co.sum()),
            "quality": quality,
        })

    for ax in axes[:2].ravel():
        ax.set_xlabel("hg [k ADC x sample]", fontsize=8)
        ax.set_xlim(0, HG_MAX / 1e3)
        if not args.linear:
            top = max((np.max(line.get_ydata()) for line in ax.get_lines() if len(line.get_ydata())), default=10)
            ax.set_yscale("log")
            ax.set_ylim(1, 2 * top)
        ax.grid(alpha=0.3)
    for ax in axes[2]:
        ax.set_xlabel("energy [keV]", fontsize=8)
        ax.set_ylabel("peak hg [k]", fontsize=8)
        ax.set_xlim(0, 1500)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3)
    for ax in axes[3]:
        ax.set_xlabel("energy [keV]", fontsize=8)
        ax.set_ylabel("residual [ADC x sample]", fontsize=8)
        ax.set_xlim(0, 1500)
        ax.grid(alpha=0.3)

    fig.suptitle(f"{ha}: energy calibration, {temp_low:.0f}-{temp_high:.0f} degC (mean of 6 sensors), "
                 f"trigger-time hits, saturated removed")
    fig.tight_layout()

    os.makedirs(args.out_dir, exist_ok=True)
    fig_path = os.path.join(args.out_dir, f"{ha}_energy_calib.png")
    fig.savefig(fig_path, dpi=80)
    plt.close(fig)
    print(f"[SAVE] {fig_path}")

    return rows


def write_csv(rows, out_dir):
    """HA ごとに行を入れ替えて energy_calib.csv に書く (ほかの HA の行は残す)。"""

    path = os.path.join(out_dir, "energy_calib.csv")
    old = []
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            old = [row for row in csv.DictReader(f) if row["HA"] not in {r["HA"] for r in rows}]

    all_rows = natsorted(old + [{k: str(v) for k, v in r.items()} for r in rows], key=lambda r: (r["HA"], r["ch"]))
    if not all_rows:
        return

    with open(path, "w", newline="", encoding="utf-8") as f:
        # 列は今回の行に合わせる (古い版の行で足りない列は空にする)
        fieldnames = list(rows[0].keys()) if rows else list(all_rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames, restval="", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"[SAVE] {path}")


def main():

    args = parse_args()
    rows = []

    for ha in args.ha:
        rows += calibrate_ha(ha, args)

    write_csv(rows, args.out_dir)


if __name__ == "__main__":
    main()
