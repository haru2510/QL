"""1つの投影（XZ または YZ）のストリップに、散乱点を通る直線を1本・2本でfitする。

入力は散乱点を原点にした (dz, t)。dz は散乱点からの Z 差（下流が正）、t は X または Y の差 [mm]。
直線は t = slope * dz（散乱点を通る）。直線から INLIER_MM 以内のストリップをその直線に乗った点とみなす。
"""
import numpy as np
from numba import njit

INLIER_MM = 1.6               # 直線に乗ったとみなす距離
MIN_POINTS = 6                # 1本の直線に必要な点の数
TWO_LINE_MIN_TOTAL = 12       # 2本fitを試す最低の点の数
TWO_LINE_MIN_SPAN_MM = 5.0    # fitの途中で各枝に要求する Z の幅（これ未満なら2本fitは失敗）
TWO_LINE_MIN_SEP_DEG = 3.0    # 保存する2本の角度差の下限（解析ではカットでさらに絞る）
ITERATIONS = 6                # 割当てと fit の繰り返し回数

SINGLE_FIELDS = ("slope", "n_inlier", "n_all", "z_span", "rms")
DOUBLE_FIELDS = ("slope0", "slope1", "n0", "n1", "z_span0", "z_span1", "z_max0", "z_max1")


@njit(cache=True)
def _distance(dz, t, angle):
    """角度 angle（Z軸から）の、散乱点を通る直線までの距離。"""
    return np.abs(t * np.cos(angle) - dz * np.sin(angle))


@njit(cache=True)
def _vote(dz, t, usable, step_deg, limit_deg):
    """直線の角度を step_deg 刻みで振り、usable な点が一番多く乗る角度を返す。"""
    best_count, best_angle = -1, 0.0
    n_steps = int(limit_deg / step_deg)
    for k in range(-n_steps, n_steps + 1):
        angle = np.radians(k * step_deg)
        count = np.sum(usable & (_distance(dz, t, angle) <= INLIER_MM))
        if count > best_count:
            best_count, best_angle = count, angle
    return best_angle


@njit(cache=True)
def one_line(dz, t):
    """1本の直線。戻り値は SINGLE_FIELDS の順 (slope, n_inlier, n_all, z_span, rms)。失敗は slope = NaN。"""
    out = np.full(5, np.nan)
    n_all = len(dz)
    out[2] = n_all
    if n_all < MIN_POINTS:
        out[1] = 0
        return out
    angle = _vote(dz, t, np.ones(n_all, np.bool_), 1.0, 85.0)
    mask = _distance(dz, t, angle) <= INLIER_MM
    for _ in range(5):  # 乗った点で最小二乗し、乗る点を選び直す
        if mask.sum() < MIN_POINTS:
            out[1] = mask.sum()
            return out
        angle = np.arctan(np.sum(dz[mask] * t[mask]) / np.sum(dz[mask] ** 2))
        mask = _distance(dz, t, angle) <= INLIER_MM
    if mask.sum() < MIN_POINTS:
        out[1] = mask.sum()
        return out
    out[0] = np.tan(angle)
    out[1] = mask.sum()
    out[3] = dz[mask].max() - dz[mask].min()
    out[4] = np.sqrt(np.mean(_distance(dz[mask], t[mask], angle) ** 2))
    return out


@njit(cache=True)
def two_lines(dz, t):
    """散乱点を通る2本の直線。戻り値は DOUBLE_FIELDS の順（傾きの小さい枝が0）。失敗は全部 NaN。

    1. 投票で1本目の角度を決め、その直線に乗った点を除いてもう一度投票して2本目を決める。
    2. 各点を近い方の直線に割り当て（INLIER_MM より遠い点はどちらにも入れない）、それぞれ直線fitし直す。
       これを ITERATIONS 回繰り返す。途中で枝の点が MIN_POINTS 未満か Z 幅が TWO_LINE_MIN_SPAN_MM 未満なら失敗。
    """
    out = np.full(8, np.nan)
    n_all = len(dz)
    if n_all < TWO_LINE_MIN_TOTAL:
        return out
    angles = np.empty(2)
    usable = np.ones(n_all, np.bool_)
    for branch in range(2):
        angles[branch] = _vote(dz, t, usable, 0.5, 85.0)
        usable &= _distance(dz, t, angles[branch]) > INLIER_MM
    masks = np.zeros((2, n_all), np.bool_)
    for _ in range(ITERATIONS):
        d0 = _distance(dz, t, angles[0])
        d1 = _distance(dz, t, angles[1])
        masks[0] = (d0 <= d1) & (d0 <= INLIER_MM)
        masks[1] = (d1 < d0) & (d1 <= INLIER_MM)
        for branch in range(2):
            m = masks[branch]
            if m.sum() < MIN_POINTS or np.ptp(dz[m]) < TWO_LINE_MIN_SPAN_MM:
                return out
            angles[branch] = np.arctan(np.sum(dz[m] * t[m]) / np.sum(dz[m] ** 2))
    if abs(angles[0] - angles[1]) < np.radians(TWO_LINE_MIN_SEP_DEG):
        return out
    order = np.argsort(angles)
    for k in range(2):
        m = masks[order[k]]
        out[k] = np.tan(angles[order[k]])
        out[2 + k] = m.sum()
        out[4 + k] = np.ptp(dz[m])
        out[6 + k] = dz[m].max()
    return out


def fit_event(hits):
    """1イベントの (axis, dz, t) の配列から、XZ・YZ それぞれの1本fit [2,5] と2本fit [2,8]。"""
    single = np.full((2, 5), np.nan)
    double = np.full((2, 8), np.nan)
    for axis in range(2):
        points = hits[hits[:, 0] == axis]
        dz = np.ascontiguousarray(points[:, 1])
        t = np.ascontiguousarray(points[:, 2])
        single[axis] = one_line(dz, t)
        double[axis] = two_lines(dz, t)
    return single, double


def fit_constants():
    """保存ファイルに残す fit の定数。"""
    return {"inlier_mm": INLIER_MM, "min_points": MIN_POINTS, "two_line_min_total": TWO_LINE_MIN_TOTAL,
            "two_line_min_span_mm": TWO_LINE_MIN_SPAN_MM, "two_line_min_separation_deg": TWO_LINE_MIN_SEP_DEG,
            "iterations": ITERATIONS}
