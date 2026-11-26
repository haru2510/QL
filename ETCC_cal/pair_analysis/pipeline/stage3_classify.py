"""stage3: cut set ごとにカットをかけ、イベントを分類し、カットの流れ（各段階で残る数）を表にする。

分類（投影ごとに「2本」「きれいな1本」を判定してから組み合わせる）:
  投影が「2本」: 2本fitが成立し、角度差 >= two_branch.min_separation_deg、各枝の Z 幅 >= two_branch.min_branch_span_mm、
                 2本に乗ったストリップの割合 >= two_branch.min_inlier_fraction
  投影が「1本」: 1本fitが成立し、Z 幅 >= single_line.min_span_mm、乗った割合 >= single_line.min_inlier_fraction、かつ「2本」でない
  SS = 1本+1本（1本の飛跡: コンプトン電子の候補）
  DS = XZ 2本 + YZ 1本（対生成の候補、e+e- 面がほぼ X）
  SD = XZ 1本 + YZ 2本（対生成の候補、面がほぼ Y）
  DD = 2本+2本（対生成の候補、面の角 ψ を計算する）

DD の面の角: XZ・YZ それぞれで2本の傾きの差 Δx' = (dx/dz)1 − (dx/dz)0, Δy' = (dy/dz)1 − (dy/dz)0。
XZ と YZ の枝の組み方が2通りあり区別できないので、X 軸で折り返した ψ = atan(|Δy'|/|Δx'|)（0〜90°）を使う。

1本の飛跡（SS）の向き: 傾き (dx/dz, dy/dz, 1) から、ビーム軸に対する θ と方位角 φ（0〜360°、0 = 検出器の X）。

DS・SD の ψ は1つに決まらないので、とりうる範囲 [psi_lo, psi_hi] を保存する（psi_range の説明を参照）。

出力: stage3_classify/<cut_set>/<sample>.npz（run, event, cls, psi, psi_lo, psi_hi, theta, phi）, cutflow.json, cutflow.md
"""
import json

import numpy as np

from pixel_polarization import fitted_axis

CLASS_NAMES = {0: "other", 1: "SS", 2: "DS", 3: "SD", 4: "DD"}
SS, DS, SD, DD = 1, 2, 3, 4


def load_sample(pipeline, sample):
    """サンプルのランの stage2 結果をつなげる。run は sample["runs"] の何番目か。"""
    parts = {key: [] for key in ("event", "vertex", "beam_distance", "z_vetoed", "single", "double", "run")}
    for index, run_name in enumerate(sample["runs"]):
        with np.load(pipeline.stage_dir(2, sample["dataset"]) / f"{run_name}.npz") as saved:
            for key in ("event", "vertex", "beam_distance", "z_vetoed", "single", "double"):
                parts[key].append(saved[key])
            parts["run"].append(np.full(len(saved["event"]), index, np.int16))
    return {key: np.concatenate(values) for key, values in parts.items()}


def projection_flags(single, double, cuts):
    """投影ごとの判定 [N, 2]。カットの流れ用に、2本の条件を1つずつ足した途中の結果も返す。"""
    n_all = np.maximum(single[:, :, 2], 1)
    line = cuts["single_line"]
    branch = cuts["two_branch"]
    with np.errstate(invalid="ignore"):
        good_line = (np.isfinite(single[:, :, 0]) & (single[:, :, 3] >= line["min_span_mm"])
                     & (single[:, :, 1] / n_all >= line["min_inlier_fraction"]))
        separation = np.degrees(np.abs(np.arctan(double[:, :, 1]) - np.arctan(double[:, :, 0])))
        steps = {"two_line_fit": np.isfinite(double[:, :, 0])}
        steps["separation"] = steps["two_line_fit"] & (separation >= branch["min_separation_deg"])
        steps["branch_span"] = steps["separation"] & (np.minimum(double[:, :, 4], double[:, :, 5]) >= branch["min_branch_span_mm"])
        steps["inlier_fraction"] = steps["branch_span"] & ((double[:, :, 2] + double[:, :, 3]) / n_all >= branch["min_inlier_fraction"])
    two = steps["inlier_fraction"]
    one = good_line & ~two
    return one, two, steps


def classify(single, double, cuts):
    """各イベントの分類（CLASS_NAMES）と DD の ψ [deg]。"""
    one, two, _ = projection_flags(single, double, cuts)
    cls = np.zeros(len(single), np.int8)
    cls[one[:, 0] & one[:, 1]] = SS
    cls[two[:, 0] & one[:, 1]] = DS
    cls[one[:, 0] & two[:, 1]] = SD
    cls[two[:, 0] & two[:, 1]] = DD
    dx = double[:, 0, 1] - double[:, 0, 0]
    dy = double[:, 1, 1] - double[:, 1, 0]
    psi = np.where(cls == DD, np.degrees(np.arctan2(np.abs(dy), np.abs(dx))), np.nan)
    return cls, psi


def psi_range(cls, psi, double, cuts):
    """各イベントの ψ がとりうる範囲 [ψ_lo, ψ_hi] [deg]。

    DD: ψ_lo = ψ_hi = ψ。
    DS: YZ は1本に見えたので、YZ での2本の傾きの差は tan(min_separation_deg) より小さいとみなす。
        ψ は 0 から ψ_max = atan(tan(min_separation_deg) / |Δx'|)（45° を上限）の間。
    SD: 同様に ψ は 90° − atan(tan(min_separation_deg) / |Δy'|) から 90° の間。
    それ以外: NaN。
    """
    limit = np.tan(np.radians(cuts["two_branch"]["min_separation_deg"]))
    dx = np.abs(double[:, 0, 1] - double[:, 0, 0])
    dy = np.abs(double[:, 1, 1] - double[:, 1, 0])
    lo = np.full(len(cls), np.nan)
    hi = np.full(len(cls), np.nan)
    dd, ds, sd = cls == DD, cls == DS, cls == SD
    lo[dd] = hi[dd] = psi[dd]
    lo[ds] = 0.0
    hi[ds] = np.minimum(np.degrees(np.arctan(limit / dx[ds])), 45.0)
    lo[sd] = 90.0 - np.minimum(np.degrees(np.arctan(limit / dy[sd])), 45.0)
    hi[sd] = 90.0
    return lo, hi


def track_angles(single, beam_direction):
    """1本fitの傾きから、ビーム軸に対する θ [deg] と方位角 φ [deg, 0-360, 0 = 検出器の X]。"""
    direction = np.stack([single[:, 0, 0], single[:, 1, 0], np.ones(len(single))], axis=1)
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    x_ref = np.array([1.0, 0.0, 0.0]) - beam_direction[0] * beam_direction  # X をビームに垂直な面へ射影
    x_ref /= np.linalg.norm(x_ref)
    y_ref = np.cross(beam_direction, x_ref)
    theta = np.degrees(np.arccos(np.clip(direction @ beam_direction, -1, 1)))
    phi = np.degrees(np.arctan2(direction @ y_ref, direction @ x_ref)) % 360
    return theta, phi


def cut_flow(pipeline, sample, data, cuts, base, cls):
    """各段階で残る数（上から順にカットを足す）。"""
    stage1 = json.loads((pipeline.stage_dir(1) / "summary.json").read_text())[sample["dataset"]]
    runs = [stage1[r] for r in sample["runs"]]
    triggers = sum(r["triggers_max_trig"] for r in runs)
    flow = [("all triggers (max trig, approx.)", triggers)]
    if cuts["z_veto"]:
        flow.append(("no raw TPC hit at Z <= z_veto", triggers - sum(r["z_vetoed_events"] for r in runs)))
        flow.append(("has a vertex (Z > z_veto)", sum(r["vertices_after_z_veto"] for r in runs)))
    else:
        flow.append(("has a vertex", sum(r["vertices"] for r in runs)))
    flow.append((f"vertex within {cuts['beam_radius_mm']:g} mm of the beam axis", int(base.sum())))
    single, double = data["single"][base], data["double"][base]
    flow.append(("strips in both XZ and YZ", int((single[:, :, 2] > 0).all(axis=1).sum())))
    _, _, steps = projection_flags(single, double, cuts)
    branch = cuts["two_branch"]
    labels = {"two_line_fit": "[pair] 2-line fit in XZ or YZ",
              "separation": f"[pair]  + branch separation >= {branch['min_separation_deg']:g} deg",
              "branch_span": f"[pair]  + each branch Z span >= {branch['min_branch_span_mm']:g} mm",
              "inlier_fraction": f"[pair]  + strips on the 2 lines >= {branch['min_inlier_fraction']:g}"}
    for key, label in labels.items():
        flow.append((label, int(steps[key].any(axis=1).sum())))
    c = cls[base]
    flow.append(("[pair]  + other projection is 2 lines or a clean line (DS+SD+DD)", int(np.isin(c, (DS, SD, DD)).sum())))
    for code in (DD, DS, SD):
        flow.append((f"[pair]      {CLASS_NAMES[code]}", int((c == code).sum())))
    flow.append(("[single track] SS (clean single line in XZ and YZ)", int((c == SS).sum())))
    return [{"step": step, "events": int(n), "percent_of_triggers": 100.0 * n / triggers if triggers else None} for step, n in flow]


def run(pipeline, cut_set_names=None):
    print("== stage3: classify", flush=True)
    _, beam_direction = fitted_axis(pipeline.beam_axis_file)
    for cut_name, cuts in pipeline.selected_cut_sets(cut_set_names).items():
        folder = pipeline.stage_dir(3, cut_name)
        flows = {}
        for sample_name, sample in pipeline.samples.items():
            data = load_sample(pipeline, sample)
            base = data["beam_distance"] <= cuts["beam_radius_mm"]
            if cuts["z_veto"]:
                base &= ~data["z_vetoed"]
            cls, psi = classify(data["single"], data["double"], cuts)
            psi_lo, psi_hi = psi_range(cls, psi, data["double"], cuts)
            theta, phi = track_angles(data["single"], beam_direction)
            np.savez_compressed(folder / f"{sample_name}.npz", run=data["run"][base], event=data["event"][base],
                                cls=cls[base], psi=psi[base], psi_lo=psi_lo[base], psi_hi=psi_hi[base], theta=theta[base], phi=phi[base],
                                runs=np.array(sample["runs"]))
            flows[sample_name] = cut_flow(pipeline, sample, data, cuts, base, cls)
            counts = {CLASS_NAMES[k]: int((cls[base] == k).sum()) for k in CLASS_NAMES}
            print(f"  [{cut_name}] {sample_name}: {counts}", flush=True)
        (folder / "cutflow.json").write_text(json.dumps({"cuts": cuts, "samples": flows}, indent=1) + "\n")
        (folder / "cutflow.md").write_text(cutflow_markdown(cut_name, cuts, flows, pipeline))
        print(f"  saved {folder / 'cutflow.md'}", flush=True)


def cutflow_markdown(cut_name, cuts, flows, pipeline):
    """サンプルを列に並べたカットの流れの表。"""
    names = list(flows)
    lines = [f"# Cut flow: {cut_name}", "", "```json", json.dumps(cuts, indent=1), "```", "",
             "| step | " + " | ".join(pipeline.samples[n]["label"] for n in names) + " |",
             "|---|" + "---|" * len(names)]
    for i, row in enumerate(flows[names[0]]):
        cells = []
        for n in names:
            r = flows[n][i]
            cells.append(f"{r['events']:,} ({r['percent_of_triggers']:.3g}%)" if r["percent_of_triggers"] is not None else f"{r['events']:,}")
        lines.append(f"| {row['step']} | " + " | ".join(cells) + " |")
    lines += ["", "Percent is relative to all triggers. Trigger counts use the largest TPC trig number (TPC events are trig 1..N)."]
    return "\n".join(lines) + "\n"
