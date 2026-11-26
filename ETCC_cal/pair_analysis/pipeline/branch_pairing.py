"""XZ と YZ の2本の枝を下流まで追い、どの XZ の枝とどの YZ の枝が同じ粒子かを決める材料を作る。

XZ（アノード）と YZ（カソード）は別々に読むので、2本の飛跡があると、XZ の枝 (a0, a1) と YZ の枝 (b0, b1) の
組み方が2通りある（A: a0-b0 と a1-b1、B: a0-b1 と a1-b0）。同じ粒子なら同じクロックで同じ電子の雲を見ているので、
次の2つで組み方を決められる。

1. クロックごとのストリップ数の相関: 電離の揺らぎで電荷が多いクロックでは、両方の投影で鳴るストリップが増える。
   XZ の枝と YZ の枝のストリップ数の相関（Pearson）を取り、相関の大きい組み方を選ぶ。
2. 枝の端の Z: e+ と e- は止まる・TPC から出る位置が違う。XZ と YZ で終わる Z がそろう組み方を選ぶ。

枝の追い方（投影ごと）:
- 出発点は stage2 の2本 fit の直線（散乱点から scan.radius_mm 以内で決めたもの）。
- 下流へクロックごとに、そのクロックのストリップを「隣り合う塊」（cluster_gap_mm より離れたら別）に分け、
  各枝の予測位置から window_mm 以内の塊を割り当てる。1つの塊は1本の枝にだけ割り当てる（近い組から順に）。
- 予測位置は、直近 refit_length_mm の塊の重心に直線を当てて延ばす（点が少ないうちは出発点の直線）。
- 2本の予測位置の差が shared_mm より小さいクロックは、どちらの枝の塊か分からないので、位置を更新せず、
  ストリップ数も記録しない（出発点の直線のまま延ばす）。
- 枝の塊が gap_mm 続けて無ければ、その枝は終わり（端の Z = 最後に塊があった Z）。
"""
import numpy as np


def follow_branches(points, slopes, settings, mm_per_clock, return_paths=False, trace=None):
    """1投影の2本を追う。points は (dz, t) [N, 2]、slopes は出発点の2本の傾き。

    戻り値: counts（各枝の {クロック番号: ストリップ数}）、ends（各枝の端の dz）、shared（区別できないクロック番号の集合）。
    クロック番号は dz / mm_per_clock を丸めたもの（XZ と YZ で同じ番号になる）。
    return_paths=True なら、各枝に割り当てた塊の重心の列 [(dz, t), ...] も返す（stage5 の絵用）。
    trace にリストを渡すと、クロックごとの処理（塊・予測位置・割り当て・重なり）を辞書で追記する（README の図用）。
    """
    window, gap = settings["window_mm"], settings["gap_mm"]
    state = [{"dz": [0.0], "t": [0.0], "slope": float(s), "last": 0.0, "alive": True} for s in slopes]
    counts = [dict(), dict()]
    shared = set()
    for dz in np.unique(points[:, 0]):
        clock = int(round(dz / mm_per_clock))
        t = np.sort(points[points[:, 0] == dz, 1])
        groups = np.split(t, np.flatnonzero(np.diff(t) > settings["cluster_gap_mm"]) + 1)
        centre = np.array([g.mean() for g in groups])
        size = np.array([len(g) for g in groups])
        predicted = [predict(b, dz, settings["refit_length_mm"]) for b in state]
        record = None
        if trace is not None:
            record = {"dz": float(dz), "strips": t, "groups": groups, "predicted": list(predicted),
                      "alive": [b["alive"] for b in state], "shared": False, "assigned": {}}
            trace.append(record)
        if abs(predicted[0] - predicted[1]) < settings["shared_mm"]:
            shared.add(clock)
            if record is not None:
                record["shared"] = True
            for b in state:  # 重なりの中では終わらせない
                b["last"] = dz
            continue
        options = sorted((abs(centre[g] - predicted[k]), k, g) for k in range(2) for g in range(len(groups))
                         if state[k]["alive"] and abs(centre[g] - predicted[k]) <= window)
        assigned, taken = {}, set()
        for _, k, g in options:
            if k not in assigned and g not in taken:
                assigned[k] = g
                taken.add(g)
        if record is not None:
            record["assigned"] = dict(assigned)
        for k, b in enumerate(state):
            if not b["alive"]:
                continue
            if k in assigned:
                g = assigned[k]
                b["dz"].append(dz)
                b["t"].append(float(centre[g]))
                b["last"] = dz
                counts[k][clock] = int(size[g])
            elif dz - b["last"] > gap:
                b["alive"] = False
    if return_paths:
        paths = [np.column_stack([b["dz"][1:], b["t"][1:]]) for b in state]
        return counts, [b["last"] for b in state], shared, paths
    return counts, [b["last"] for b in state], shared


def predict(branch, dz, refit_length):
    """直近 refit_length の塊の重心に直線を当てて dz での位置を予測する（5点未満なら出発点の直線）。"""
    z, t = np.array(branch["dz"]), np.array(branch["t"])
    near = z >= z[-1] - refit_length
    if near.sum() >= 5:
        slope, intercept = np.polyfit(z[near], t[near], 1)
        return slope * dz + intercept
    return branch["slope"] * dz


def correlation(count_x, count_y, excluded, min_clocks):
    """2つの枝のクロックごとのストリップ数の相関と、使ったクロックの数。"""
    clocks = sorted((set(count_x) & set(count_y)) - excluded)
    if len(clocks) < min_clocks:
        return np.nan, len(clocks)
    a = np.array([count_x[c] for c in clocks], float)
    b = np.array([count_y[c] for c in clocks], float)
    if a.std() == 0 or b.std() == 0:
        return np.nan, len(clocks)
    return float(np.corrcoef(a, b)[0, 1]), len(clocks)


def pairing_features(hits, double, settings, mm_per_clock):
    """1イベントの組み方の材料。

    戻り値: r[2,2]（XZ の枝 i と YZ の枝 j の相関）, n[2,2]（使ったクロック数）, end[2,2]（[XZ/YZ, 枝] の端の dz）。
    """
    follows = []
    for axis in range(2):
        points = hits[hits[:, 0] == axis][:, 1:]
        follows.append(follow_branches(points, double[axis, :2], settings, mm_per_clock))
    (count_x, end_x, shared_x), (count_y, end_y, shared_y) = follows
    excluded = shared_x | shared_y
    r = np.full((2, 2), np.nan)
    n = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            r[i, j], n[i, j] = correlation(count_x[i], count_y[j], excluded, settings["min_clocks"])
    return r, n, np.array([end_x, end_y])


def decide(r, end, rule):
    """組み方を決める: +1 = A（a0-b0, a1-b1）、-1 = B（a0-b1, a1-b0）、0 = 決められない。どの方法で決めたかも返す。

    相関の差が min_r_difference 以上ならそれで決め、そうでなく端の Z の差が min_end_difference_mm 以上なら端で決める。
    """
    r_a = np.nanmean([r[0, 0], r[1, 1]]) if np.isfinite([r[0, 0], r[1, 1]]).any() else np.nan
    r_b = np.nanmean([r[0, 1], r[1, 0]]) if np.isfinite([r[0, 1], r[1, 0]]).any() else np.nan
    if np.isfinite(r_a) and np.isfinite(r_b) and abs(r_a - r_b) >= rule["min_r_difference"]:
        return (1 if r_a > r_b else -1), "correlation"
    end_a = abs(end[0, 0] - end[1, 0]) + abs(end[0, 1] - end[1, 1])
    end_b = abs(end[0, 0] - end[1, 1]) + abs(end[0, 1] - end[1, 0])
    if abs(end_a - end_b) >= rule["min_end_difference_mm"]:
        return (1 if end_a < end_b else -1), "end_z"
    return 0, "undecided"
