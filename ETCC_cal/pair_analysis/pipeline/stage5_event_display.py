"""stage5: 分類ごとに、散乱点まわりの XZ・YZ の raw ストリップと、fit した直線を描く（目で確かめる用）。

event_display の設定: cut_sets（描く cut set）, dataset, run, classes, count（各分類の先頭から何件）, skip, context_mm。
黒: fit に使った範囲（scan.radius_mm 以内）のストリップ、灰色: context_mm 以内の全ストリップ。
出力: stage5_event_display/<cut_set>/<run>_<class>.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from stage_coincidence_xy import scan_run
from tpc_projected_pair import projected_hits
from vertex_sphere import take_local_cells

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
        label = f"2 lines: {a0:.1f}°, {a1:.1f}° (sep {a1 - a0:.1f}°), Z span {double[4]:.0f}/{double[5]:.0f} mm"
    else:
        if np.isfinite(single[0]):
            ax.plot(z, single[0] * z, color=LINE_COLOR, lw=1.4, alpha=0.85)
        label = f"1 line: {np.degrees(np.arctan(single[0])):.1f}°, on line {int(single[1])}/{int(single[2])}"
    ax.add_patch(plt.Circle((0, 0), radius, fill=False, ls=":", color="0.5"))
    ax.plot(0, 0, "r*", ms=9)
    ax.set(xlim=(-5, context), ylim=(-context / 2, context / 2), xlabel="Z − Z_vertex [mm]", ylabel=f"{name} − {name}_vertex [mm]")
    ax.set_aspect("equal")
    ax.set_title(f"{name}Z: {label}", fontsize=8)


def run(pipeline, cut_set_names=None):
    settings = pipeline.event_display
    if not settings.get("enabled", False):
        print("== stage5: event display disabled in the config", flush=True)
        return
    print("== stage5: event display", flush=True)
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
                axes[r, 0].text(0.01, 0.97, f"{run_name} event {event >> 32}:{event & 0xFFFFFFFF}\nvertex ({x:.1f}, {y:.1f}, {z:.1f}) mm{extra}",
                                transform=axes[r, 0].transAxes, va="top", fontsize=7)
            fig.suptitle(f"[{cut_name}] {CLASS_NAMES[code]} events, {run_name}: black = strips used (≤ {radius:g} mm), grey = up to {context:g} mm")
            target = folder / f"{run_name}_{CLASS_NAMES[code]}_{settings['skip']}.png"
            fig.savefig(target, dpi=110)
            plt.close(fig)
            print(f"  saved {target}", flush=True)
