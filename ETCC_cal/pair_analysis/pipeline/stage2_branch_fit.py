"""stage2: raw TPC を読み、散乱点ごとに XZ・YZ の1本fitと2本fitを保存する（時間がかかるのはここだけ）。

- 対象は stage1 の散乱点のうち、ビーム軸から scan.beam_radius_mm 以内のもの
  （Z 側除外のイベントは、z_veto=false の cut set があるときだけ読む）。
- 散乱点から scan.radius_mm 以内、散乱点より 3 mm 以上下流のストリップを、XZ と YZ で別々に使う。
  X と Y を組にしないので、同じクロックの X・Y の組合せの曖昧さは入らない。
- どちらの fit もカットをかけずに保存する。カットは stage3 で cut set ごとにかける（raw を読み直さずに条件を変えられる）。
- ランごとに別プロセスで並列に回す（scan.workers）。

出力: stage2_branch_fit/<dataset>/<run>.npz
  event, vertex, beam_distance, z_vetoed（stage1 と同じ並び）, single[N,2,5], double[N,2,8], info(JSON)
  single/double の中身は track_fit.SINGLE_FIELDS / DOUBLE_FIELDS（2つ目の添字 0 = XZ, 1 = YZ）。
"""
import json
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from stage_coincidence_xy import scan_run
from tpc_projected_pair import projected_hits
from vertex_sphere import take_local_cells

from .config import Pipeline
from .track_fit import fit_constants, fit_event


def output_file(pipeline, dataset, run_name):
    return pipeline.stage_dir(2, dataset) / f"{run_name}.npz"


def scan_one_run(config_path, dataset, run_name):
    """1ランの raw を読み、選んだ散乱点ごとに fit する（別プロセスで呼ばれる）。"""
    pipeline = Pipeline(config_path)
    config = pipeline.etcc_config(dataset)
    selected = np.load(pipeline.stage_dir(1, dataset) / f"{run_name}.npz")
    use = selected["beam_distance"] <= pipeline.scan["beam_radius_mm"]
    if not pipeline.scan_includes_z_vetoed:
        use &= ~selected["z_vetoed"]
    events = selected["event"][use]
    vertices = selected["vertex"][use]
    row_of = {int(e): i for i, e in enumerate(events)}
    radius = pipeline.scan["radius_mm"]
    calibration = (config["velocity"], config["clock_offset"], config["z_gem"])
    hits_of = {}

    def collect(event, upic, occupancy, occupied, counts):
        """scan_run が uPIC ごと・イベントごとに呼ぶ。散乱点のまわりのストリップを取り出す。"""
        row = row_of.get(event)
        if row is None:  # 選んでいないイベント: 占有を消すだけ
            take_local_cells(occupancy, occupied, counts, np.zeros(3), 0.0, *calibration)
            return 0
        cells = take_local_cells(occupancy, occupied, counts, vertices[row], radius, *calibration)
        hits = projected_hits(cells, vertices[row], radius, *calibration)  # (axis, dz, t)、dz >= 3 mm
        if hits:
            hits_of.setdefault(row, []).append(np.asarray(hits))
        return len(hits)

    print(f"  {dataset} {run_name}: scanning {len(events):,} vertices", flush=True)
    scan_run(run_name, config, [events], event_processor=collect)
    single = np.full((len(events), 2, 5), np.nan)
    double = np.full((len(events), 2, 8), np.nan)
    for row, chunks in hits_of.items():
        hits = np.unique(np.concatenate(chunks), axis=0)  # 複数 uPIC に同じストリップが出たら1回に
        single[row], double[row] = fit_event(hits)
    info = {"dataset": dataset, "run": run_name, "radius_mm": radius, "scan_beam_radius_mm": pipeline.scan["beam_radius_mm"],
            "include_z_vetoed": pipeline.scan_includes_z_vetoed, "fit": fit_constants()}
    target = output_file(pipeline, dataset, run_name)
    temporary = target.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, event=events, vertex=vertices, beam_distance=selected["beam_distance"][use],
                        z_vetoed=selected["z_vetoed"][use], single=single, double=double, info=json.dumps(info))
    temporary.replace(target)
    return f"{dataset} {run_name}: saved {target}"


def run(pipeline, rescan=False):
    print("== stage2: branch fit (raw TPC scan)", flush=True)
    tasks = []
    for dataset, settings in pipeline.datasets.items():
        for run_name in settings["runs"]:
            target = output_file(pipeline, dataset, run_name)
            if rescan or not target.exists() or not same_scan(target, pipeline):
                tasks.append((dataset, run_name))
    if not tasks:
        print("  all runs already scanned", flush=True)
        return
    print(f"  {len(tasks)} runs to scan with {pipeline.scan['workers']} workers", flush=True)
    with ProcessPoolExecutor(max_workers=pipeline.scan["workers"]) as pool:
        futures = [pool.submit(scan_one_run, str(pipeline.path), dataset, run_name) for dataset, run_name in tasks]
        for future in as_completed(futures):
            print("  " + future.result(), flush=True)


def same_scan(target, pipeline):
    """保存済みファイルが今のスキャン条件で作られたものか。"""
    with np.load(target) as saved:
        info = json.loads(str(saved["info"]))
    return (info["radius_mm"] == pipeline.scan["radius_mm"]
            and info["scan_beam_radius_mm"] == pipeline.scan["beam_radius_mm"]
            and info["include_z_vetoed"] == pipeline.scan_includes_z_vetoed
            and info["fit"] == fit_constants())
