"""モジュレーションの fit と有意度。"""
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit


def normalized_ratio(numerator, denominator):
    """2つのヒストグラムをそれぞれ合計で規格化した比と、その統計誤差（ポアソン）。"""
    a = numerator.astype(float)
    b = denominator.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = (a / a.sum()) / (b / b.sum())
        error = ratio * np.sqrt(1 / a + 1 / b)
    ok = (a > 0) & (b > 0)
    return ratio, error, ok


def flat_chi2(values, errors):
    """重み付き平均（平らな線）での χ²。"""
    w = 1 / errors ** 2
    mean = np.sum(w * values) / w.sum()
    return float(np.sum(((values - mean) / errors) ** 2))


def significance(delta_chi2, extra_parameters):
    """平らな線に対する Δχ²（自由度 extra_parameters）の p 値を、両側ガウスの σ に直す（1自由度なら √Δχ²）。"""
    p = stats.chi2.sf(max(delta_chi2, 0.0), extra_parameters)
    return float(stats.norm.isf(p / 2)) if p > 0 else float("inf"), float(p)


def fit_cos2_fixed_phase(psi_deg, ratio, error):
    """A [1 + μ cos 2ψ]（位相は 0 に固定。折り返した ψ 用）。μ の符号: 正なら ψ = 0（X）側が多い。"""
    model = lambda x, a, mu: a * (1 + mu * np.cos(2 * np.radians(x)))
    p, cov = curve_fit(model, psi_deg, ratio, p0=[1.0, 0.0], sigma=error, absolute_sigma=True)
    chi2 = float(np.sum(((ratio - model(psi_deg, *p)) / error) ** 2))
    chi2_flat = flat_chi2(ratio, error)
    sigma, p_value = significance(chi2_flat - chi2, 1)
    return {"A": float(p[0]), "mu": float(p[1]), "mu_err": float(np.sqrt(cov[1, 1])), "chi2": chi2, "ndf": len(ratio) - 2,
            "chi2_flat": chi2_flat, "ndf_flat": len(ratio) - 1, "significance_sigma": sigma, "p_value": p_value}, model


def fit_cos2_free_phase(phi_deg, ratio, error, widths_deg=None):
    """A [1 + μ cos 2(φ − φ0)]（μ >= 0 にそろえ、φ0 は 0〜180°）。

    widths_deg を渡すと、各ビンの中で平均したモデル A [1 + μ cos 2(c − φ0) sin w / w]（c = ビンの中心、w = ビン幅 [rad]）で fit する
    （幅の違うビンを混ぜても μ が偏らない）。返す model は平均する前の曲線。
    """
    model = lambda x, a, mu, phi0: a * (1 + mu * np.cos(2 * np.radians(x - phi0)))
    if widths_deg is not None:
        w = np.radians(np.asarray(widths_deg, float))
        damping = np.sin(w) / w
        fit_model = lambda x, a, mu, phi0: a * (1 + mu * damping * np.cos(2 * np.radians(x - phi0)))
    else:
        fit_model = model
    best = None
    for start in (0.0, 45.0, 90.0, 135.0):  # 位相の初期値を変えて一番よいものを使う
        try:
            p, cov = curve_fit(fit_model, phi_deg, ratio, p0=[1.0, 0.05, start], sigma=error, absolute_sigma=True, maxfev=20000)
        except RuntimeError:
            continue
        chi2 = float(np.sum(((ratio - fit_model(phi_deg, *p)) / error) ** 2))
        if best is None or chi2 < best[2]:
            best = (p, cov, chi2)
    p, cov, chi2 = best
    if p[1] < 0:
        p[1], p[2] = -p[1], p[2] + 90
    p[2] %= 180
    err = np.sqrt(np.diag(cov))
    chi2_flat = flat_chi2(ratio, error)
    sigma, p_value = significance(chi2_flat - chi2, 2)
    return {"A": float(p[0]), "mu": float(p[1]), "mu_err": float(err[1]), "phi0_deg": float(p[2]), "phi0_err_deg": float(err[2]),
            "chi2": chi2, "ndf": len(ratio) - 3, "chi2_flat": chi2_flat, "ndf_flat": len(ratio) - 1,
            "significance_sigma": sigma, "p_value": p_value}, model


def asymmetry(n_x, n_y):
    """Q = (N_X − N_Y)/(N_X + N_Y) と二項分布の誤差。"""
    n = n_x + n_y
    q = (n_x - n_y) / n
    return {"N_X": int(n_x), "N_Y": int(n_y), "Q": float(q), "Q_err": float(np.sqrt((1 - q * q) / n))}


def fit_cos2_mean(mean_cos2, ratio, error):
    """点ごとに cos 2ψ の平均 c_k が決まっているときの A [1 + μ c_k]（DS・SD のように ψ に幅がある点も入れられる）。"""
    model = lambda c, a, mu: a * (1 + mu * c)
    p, cov = curve_fit(model, mean_cos2, ratio, p0=[1.0, 0.0], sigma=error, absolute_sigma=True)
    chi2 = float(np.sum(((ratio - model(mean_cos2, *p)) / error) ** 2))
    chi2_flat = flat_chi2(ratio, error)
    sigma, p_value = significance(chi2_flat - chi2, 1)
    return {"A": float(p[0]), "mu": float(p[1]), "mu_err": float(np.sqrt(cov[1, 1])), "chi2": chi2, "ndf": len(ratio) - 2,
            "chi2_flat": chi2_flat, "ndf_flat": len(ratio) - 1, "significance_sigma": sigma, "p_value": p_value}
