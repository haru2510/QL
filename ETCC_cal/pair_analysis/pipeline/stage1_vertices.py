"""stage1: 散乱点を読み、Z 側除外の印とビーム軸からの距離を付けて、ランごとに保存する。

- 1イベントに散乱点が複数（複数の uPIC）あるときは、CSV で最初に出てくるものだけを使う。
- 保存するのは、ビーム軸から scan.beam_radius_mm 以内の散乱点だけ（stage2 で raw を読む候補）。
  Z 側除外にかかったイベントも印を付けて残す（stage2 で読むのは、z_veto=false の cut set があるときだけ）。
- 件数（全トリガー・Z 側除外・散乱点の数など）は summary.json に残し、stage3 のカットの流れの表で使う。

出力: stage1_vertices/<dataset>/<run>.npz（event, vertex[N,3], beam_distance, z_vetoed）と summary.json
"""
import csv
import json

import numpy as np

from pixel_polarization import fitted_axis
from settings import repo_path


def beam_distance(points, beam_point, beam_direction):
    """点からビーム軸（直線）までの距離 [mm]。"""
    return np.linalg.norm(np.cross(points - beam_point, beam_direction), axis=1)


def read_first_vertices(csv_path, runs):
    """{run: {event: (x, y, z)}}。event = (epoch << 32) | trig。各イベントの最初の行だけを使う。"""
    table = {run: {} for run in runs}
    with csv_path.open() as stream:
        for row in csv.DictReader(line for line in stream if not line.startswith("#")):
            events = table.get(row["run"])
            if events is None:
                continue
            event = (int(row["epoch"]) << 32) | int(row["trig"])
            if event not in events:
                events[event] = (float(row["x"]), float(row["y"]), float(row["z"]))
    return table


def run(pipeline, rebuild=False):
    print("== stage1: vertices", flush=True)
    beam_point, beam_direction = fitted_axis(pipeline.beam_axis_file)
    scan_radius = pipeline.scan["beam_radius_mm"]
    summary_file = pipeline.stage_dir(1) / "summary.json"
    summary = json.loads(summary_file.read_text()) if summary_file.exists() and not rebuild else {}
    for name, dataset in pipeline.datasets.items():
        folder = pipeline.stage_dir(1, name)
        todo = [r for r in dataset["runs"] if rebuild or not (folder / f"{r}.npz").exists() or r not in summary.get(name, {})]
        if not todo:
            print(f"  {name}: all runs already done", flush=True)
            continue
        config = pipeline.etcc_config(name)
        z_veto = config["z_veto"]
        print(f"  {name}: reading {config['inputs']['vertices_all']}", flush=True)
        table = read_first_vertices(repo_path(config["inputs"]["vertices_all"]), todo)
        for run_name in todo:
            vetoed_events = np.load(repo_path(config["outputs"]["veto_dir"]) / f"z0_track_veto_{run_name}.npz")["event"]
            events = np.array(sorted(table[run_name]), dtype=np.int64)
            points = np.array([table[run_name][int(e)] for e in events]).reshape(-1, 3)
            # Z 側除外: イベント全体の印に加え、散乱点そのものが Z <= z_veto のものも除外扱いにする
            z_vetoed = np.isin(events, vetoed_events) | (points[:, 2] <= z_veto)
            distance = beam_distance(points, beam_point, beam_direction)
            keep = distance <= scan_radius
            np.savez_compressed(folder / f"{run_name}.npz", event=events[keep], vertex=points[keep],
                                beam_distance=distance[keep], z_vetoed=z_vetoed[keep])
            max_event = int(max(vetoed_events.max(initial=0), (events & 0xFFFFFFFF).max(initial=0)))
            summary.setdefault(name, {})[run_name] = {
                "triggers_max_trig": max_event,
                "z_vetoed_events": int(len(vetoed_events)),
                "vertices": int(len(events)),
                "vertices_after_z_veto": int((~z_vetoed).sum()),
                "within_scan_radius": int(keep.sum()),
                "within_scan_radius_after_z_veto": int((keep & ~z_vetoed).sum()),
            }
            print(f"    {run_name}: {len(events):,} vertices, {(~z_vetoed).sum():,} after Z veto, "
                  f"{(keep & ~z_vetoed).sum():,} within {scan_radius:g} mm of the beam axis", flush=True)
            summary_file.write_text(json.dumps(summary, indent=1) + "\n")
    summary_file.write_text(json.dumps(summary, indent=1) + "\n")
    print(f"  saved {summary_file}", flush=True)
