"""stage5: 分類ごとに、散乱点まわりの XZ・YZ の raw ストリップと、fit した直線を描く（目で確かめる用）。

event_display の設定: cut_sets（描く cut set）, dataset, run, classes, count（各分類の先頭から何件）, skip, context_mm,
pairing_count（枝の組み方の絵にする DD の件数。0 なら描かない）。
黒: fit に使った範囲（scan.radius_mm 以内）のストリップ、灰色: context_mm 以内の全ストリップ。
出力: stage5_event_display/<cut_set>/<run>_<class>.png
      stage5_event_display/<cut_set>/<run>_DD_pairing_<skip>.png（枝を追った様子と、クロックごとのストリップ数）
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from stage_coincidence_xy import scan_run
from tpc_projected_pair import projected_hits
from vertex_sphere import take_local_cells

from .branch_pairing import correlation, decide, follow_branches
from .plot_style import use_japanese_font
from .stage3_classify import CLASS_NAMES, DD, DS, SD, classify

BRANCH_COLORS = ("#2a78d6", "#eb6834")
LINE_COLOR = "#1baf7a"


def read_hits(pipeline, dataset, run_name, events, vertices, context):
    """選んだイベントの、散乱点から context mm 以内の (axis, dz, t)。"""
    config = pipeline.etcc_config(dataset)
    calibration = (config["velocity"], config["clock_offset"], config["z_gem"])
    row_of = {int(e): i for i, e in enumerate(events)}
    hits_of = {}

    def collect(event, upic, occupancy, occupied, counts):
        row = row_of.get(event)
        if row is None:
            take_local_cells(occupancy, occupied, counts, np.zeros(3), 0.0, *calibration)
            return 0
        cells = take_local_cells(occupancy, occupied, counts, vertices[row], context, *calibration)
        hits = projected_hits(cells, vertices[row], context, *calibration)
        if hits:
            hits_of.setdefault(row, []).append(np.asarray(hits))
        return len(hits)

    scan_run(run_name, config, [np.sort(events)], event_processor=collect)
    return {row: np.unique(np.concatenate(chunks), axis=0) for row, chunks in hits_of.items()}


def draw_projection(ax, hits, axis, single, double, is_two, radius, context):
    name = "XY"[axis]
    points = hits[hits[:, 0] == axis] if len(hits) else np.zeros((0, 3))
    near = np.hypot(points[:, 1], points[:, 2]) <= radius
    ax.scatter(points[~near, 1], points[~near, 2], s=3, color="#b8b7b2", lw=0)
    ax.scatter(points[near, 1], points[near, 2], s=4, color="black", lw=0)
    z = np.array([0.0, radius])
    if is_two:
        for branch in range(2):
            ax.plot(z, double[branch] * z, color=BRANCH_COLORS[branch], lw=1.4, alpha=0.85)
        a0, a1 = np.degrees(np.arctan(double[:2]))
        label = f"2本: {a0:.1f}°, {a1:.1f}°（差 {a1 - a0:.1f}°）、Z 幅 {double[4]:.0f}/{double[5]:.0f} mm"
    else:
        if np.isfinite(single[0]):
            ax.plot(z, single[0] * z, color=LINE_COLOR, lw=1.4, alpha=0.85)
        label = f"1本: {np.degrees(np.arctan(single[0])):.1f}°、乗った点 {int(single[1])}/{int(single[2])}"
    ax.add_patch(plt.Circle((0, 0), radius, fill=False, ls=":", color="0.5"))
    ax.plot(0, 0, "r*", ms=9)
    ax.set(xlim=(-5, context), ylim=(-context / 2, context / 2), xlabel="Z - Z_散乱点 [mm]", ylabel=f"{name} - {name}_散乱点 [mm]")
    ax.set_aspect("equal")
    ax.set_title(f"{name}Z: {label}", fontsize=8)


def run(pipeline, cut_set_names=None):
    settings = pipeline.event_display
    if not settings.get("enabled", False):
        print("== stage5: event display disabled in the config", flush=True)
        return
    print("== stage5: event display", flush=True)
    use_japanese_font()  # このステージの図は README の説明用なので日本語にする
    dataset, run_name = settings["dataset"], settings["run"]
    radius, context = pipeline.scan["radius_mm"], settings["context_mm"]
    saved = np.load(pipeline.stage_dir(2, dataset) / f"{run_name}.npz")
    names = cut_set_names or settings.get("cut_sets") or list(pipeline.cut_sets)
    for cut_name, cuts in pipeline.selected_cut_sets(names).items():
        base = saved["beam_distance"] <= cuts["beam_radius_mm"]
        if cuts["z_veto"]:
            base &= ~saved["z_vetoed"]
        cls, psi = classify(saved["single"], saved["double"], cuts)
        wanted = [code for code, name in CLASS_NAMES.items() if name in settings["classes"]]
        picks = {code: np.flatnonzero(base & (cls == code))[settings["skip"]:settings["skip"] + settings["count"]] for code in wanted}
        rows = np.concatenate(list(picks.values())).astype(int)
        hits = read_hits(pipeline, dataset, run_name, saved["event"][rows], saved["vertex"][rows], context)
        hits = {int(rows[k]): h for k, h in hits.items()}
        folder = pipeline.stage_dir(5, cut_name)
        for code, chosen in picks.items():
            if len(chosen) == 0:
                continue
            fig, axes = plt.subplots(len(chosen), 2, figsize=(11, 3.3 * len(chosen)), squeeze=False, constrained_layout=True)
            for r, i in enumerate(chosen):
                event_hits = hits.get(int(i), np.zeros((0, 3)))
                for axis in range(2):
                    is_two = code == DD or (code == DS and axis == 0) or (code == SD and axis == 1)
                    draw_projection(axes[r, axis], event_hits, axis, saved["single"][i, axis], saved["double"][i, axis], is_two, radius, context)
                event = int(saved["event"][i])
                x, y, z = saved["vertex"][i]
                extra = f", ψ = {psi[i]:.0f}°" if code == DD else ""
                axes[r, 0].text(0.01, 0.97, f"{run_name} event {event >> 32}:{event & 0xFFFFFFFF}\n散乱点 ({x:.1f}, {y:.1f}, {z:.1f}) mm{extra}",
                                transform=axes[r, 0].transAxes, va="top", fontsize=7)
            fig.suptitle(f"[{cut_name}] {CLASS_NAMES[code]} のイベント、{run_name}: 黒 = fit に使ったストリップ（{radius:g} mm 以内）、灰色 = {context:g} mm 以内")
            target = folder / f"{run_name}_{CLASS_NAMES[code]}_{settings['skip']}.png"
            fig.savefig(target, dpi=110)
            plt.close(fig)
            print(f"  saved {target}", flush=True)
        draw_pairing(pipeline, cut_name, cuts, saved, run_name, dataset, folder, settings)


def draw_pairing(pipeline, cut_name, cuts, saved, run_name, dataset, folder, settings):
    """DD の枝を追った様子と、組み方の決め方（ストリップ数の相関）を描く。

    1行 = 1イベント。左: XZ、中: YZ（灰色の点 = ストリップ、色の線 = 枝に割り当てた塊の中心。決めた組み方で
    同じ粒子の XZ と YZ の枝を同じ色にする）。右: クロックごとのストリップ数（実線 = XZ の枝、点線 = YZ の枝）。
    """
    count = settings.get("pairing_count", 0)
    if count <= 0 or not pipeline.pairing.get("enabled", False):
        return
    base = saved["beam_distance"] <= cuts["beam_radius_mm"]
    if cuts["z_veto"]:
        base &= ~saved["z_vetoed"]
    cls, _ = classify(saved["single"], saved["double"], cuts)
    rows = np.flatnonzero(base & (cls == DD))[settings["skip"]:settings["skip"] + 4 * count]
    config = pipeline.etcc_config(dataset)
    mm_per_clock = config["velocity"]
    context = pipeline.pairing["context_mm"]
    hits = read_hits(pipeline, dataset, run_name, saved["event"][rows], saved["vertex"][rows], context)
    drawn = []
    for k, i in enumerate(rows):  # 相関で組み方が決まったイベントを count 件
        if k not in hits:
            continue
        h = hits[k]
        follows = [follow_branches(h[h[:, 0] == axis][:, 1:], saved["double"][i, axis, :2], pipeline.pairing, mm_per_clock,
                                   return_paths=True) for axis in range(2)]
        excluded = follows[0][2] | follows[1][2]
        r = np.full((2, 2), np.nan)
        for a in range(2):
            for b in range(2):
                r[a, b] = correlation(follows[0][0][a], follows[1][0][b], excluded, pipeline.pairing["min_clocks"])[0]
        end = np.array([follows[0][1], follows[1][1]])
        pairing, method = decide(r, end, pipeline.pairing)
        if method == "correlation":
            drawn.append((i, h, follows, r, pairing))
        if len(drawn) == count:
            break
    if not drawn:
        return
    colors = ("#2a78d6", "#eb6834")
    fig, axes = plt.subplots(len(drawn), 3, figsize=(17, 3.6 * len(drawn)), squeeze=False, constrained_layout=True)
    for row, (i, h, follows, r, pairing) in enumerate(drawn):
        partner = [0, 1] if pairing == 1 else [1, 0]  # XZ の枝 a に対応する YZ の枝
        color_y = [colors[partner.index(b)] for b in range(2)]
        for axis, name in enumerate("XY"):
            ax = axes[row, axis]
            p = h[h[:, 0] == axis]
            ax.scatter(p[:, 1], p[:, 2], s=2, color="#b8b7b2", lw=0)
            paths = follows[axis][3]
            for b in range(2):
                color = colors[b] if axis == 0 else color_y[b]
                label = f"{'a' if axis == 0 else 'b'}{b}"
                if len(paths[b]):
                    ax.plot(paths[b][:, 0], paths[b][:, 1], color=color, lw=1.4, label=label)
            ax.plot(0, 0, "r*", ms=9)
            ax.set(xlabel="Z - Z_散乱点 [mm]", ylabel=f"{name} - {name}_散乱点 [mm]", xlim=(-5, context * 0.8))
            ax.legend(fontsize=7, frameon=False, loc="upper left")
            ax.set_title(f"{name}Z: 追った枝（同じ色 = 同じ粒子と判定）", fontsize=9)
        ax = axes[row, 2]
        for a in range(2):
            for axis, style in ((0, "-"), (1, "--")):
                branch = a if axis == 0 else partner[a]
                series = follows[axis][0][branch]
                clocks = np.array(sorted(series))
                if len(clocks):
                    ax.plot(clocks * mm_per_clock, [series[c] for c in clocks], style, color=colors[a], lw=1,
                            label=f"{'a' if axis == 0 else 'b'}{branch}")
        r_a, r_b = np.nanmean([r[0, 0], r[1, 1]]), np.nanmean([r[0, 1], r[1, 0]])
        chosen = "組み方 A（a0-b0, a1-b1）" if pairing == 1 else "組み方 B（a0-b1, a1-b0）"
        ax.set(xlabel="Z - Z_散乱点 [mm]", ylabel="その枝の1クロックのストリップ数",
               title=f"r00 {r[0, 0]:.2f}, r11 {r[1, 1]:.2f}, r01 {r[0, 1]:.2f}, r10 {r[1, 0]:.2f}  →  r_A {r_a:.2f}, r_B {r_b:.2f}: {chosen}")
        ax.title.set_fontsize(8)
        ax.legend(fontsize=7, frameon=False, ncol=2)
        event = int(saved["event"][i])
        axes[row, 0].text(0.01, 0.03, f"{run_name} event {event >> 32}:{event & 0xFFFFFFFF}", transform=axes[row, 0].transAxes, fontsize=7)
    fig.suptitle(f"[{cut_name}] DD の枝の組み方、{run_name}: 実線 = XZ の枝、点線 = それと組んだ YZ の枝（同じ色 = 同じ粒子）")
    target = folder / f"{run_name}_DD_pairing_{settings['skip']}.png"
    fig.savefig(target, dpi=110)
    plt.close(fig)
    print(f"  saved {target}", flush=True)
