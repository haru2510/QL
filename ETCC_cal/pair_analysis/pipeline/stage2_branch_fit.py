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

枝の組み方（設定の pairing.enabled が true のとき、fit の後にもう一度 raw を読む）:
- 対象は、XZ と YZ の両方で2本 fit が成立したイベント（どの cut set の DD にもなりうるもの）。
- 散乱点から pairing.context_mm 以内のストリップで2本の枝を下流まで追い、組み方の材料を保存する（branch_pairing.py）。
出力: stage2_branch_fit/<dataset>/<run>_pairing.npz
  event, r[N,2,2]（XZ の枝 i と YZ の枝 j のストリップ数の相関）, n_clocks[N,2,2], end_dz[N,2,2]（[XZ/YZ, 枝] の端）, info(JSON)
"""
import json
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from stage_coincidence_xy import scan_run
from tpc_projected_pair import projected_hits
from vertex_sphere import take_local_cells

from .branch_pairing import pairing_features
from .config import Pipeline
from .track_fit import fit_constants, fit_event

PAIRING_SCAN_KEYS = ("context_mm", "window_mm", "gap_mm", "shared_mm", "cluster_gap_mm", "refit_length_mm", "min_clocks")


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


def pairing_file(pipeline, dataset, run_name):
    return pipeline.stage_dir(2, dataset) / f"{run_name}_pairing.npz"


def pair_one_run(config_path, dataset, run_name):
    """1ランの raw をもう一度読み、両投影で2本 fit が成立したイベントの枝を追う（別プロセスで呼ばれる）。"""
    pipeline = Pipeline(config_path)
    settings = pipeline.pairing
    config = pipeline.etcc_config(dataset)
    fits = np.load(output_file(pipeline, dataset, run_name))
    rows = np.flatnonzero(np.isfinite(fits["double"][:, :, 0]).all(axis=1))
    events, vertices = fits["event"][rows], fits["vertex"][rows]
    row_of = {int(e): k for k, e in enumerate(events)}
    calibration = (config["velocity"], config["clock_offset"], config["z_gem"])
    context = settings["context_mm"]
    hits_of = {}

    def collect(event, upic, occupancy, occupied, counts):
        k = row_of.get(event)
        if k is None:
            take_local_cells(occupancy, occupied, counts, np.zeros(3), 0.0, *calibration)
            return 0
        cells = take_local_cells(occupancy, occupied, counts, vertices[k], context, *calibration)
        hits = projected_hits(cells, vertices[k], context, *calibration)
        if hits:
            hits_of.setdefault(k, []).append(np.asarray(hits))
        return len(hits)

    print(f"  {dataset} {run_name}: following branches of {len(events):,} two-line events", flush=True)
    scan_run(run_name, config, [np.sort(events)], event_processor=collect)
    r = np.full((len(events), 2, 2), np.nan)
    n_clocks = np.zeros((len(events), 2, 2))
    end_dz = np.full((len(events), 2, 2), np.nan)
    for k, chunks in hits_of.items():
        hits = np.unique(np.concatenate(chunks), axis=0)
        r[k], n_clocks[k], end_dz[k] = pairing_features(hits, fits["double"][rows[k]], settings, config["velocity"])
    info = {"dataset": dataset, "run": run_name, **{key: settings[key] for key in PAIRING_SCAN_KEYS}}
    target = pairing_file(pipeline, dataset, run_name)
    temporary = target.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, event=events, r=r, n_clocks=n_clocks, end_dz=end_dz, info=json.dumps(info))
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
    run_parallel(pipeline, scan_one_run, tasks, "fit")
    if not pipeline.pairing.get("enabled", False):
        return
    tasks = []
    for dataset, settings in pipeline.datasets.items():
        for run_name in settings["runs"]:
            target = pairing_file(pipeline, dataset, run_name)
            if rescan or not target.exists() or not same_pairing(target, pipeline):
                tasks.append((dataset, run_name))
    run_parallel(pipeline, pair_one_run, tasks, "branch pairing")


def run_parallel(pipeline, function, tasks, what):
    if not tasks:
        print(f"  {what}: all runs already done", flush=True)
        return
    print(f"  {what}: {len(tasks)} runs with {pipeline.scan['workers']} workers", flush=True)
    with ProcessPoolExecutor(max_workers=pipeline.scan["workers"]) as pool:
        futures = [pool.submit(function, str(pipeline.path), dataset, run_name) for dataset, run_name in tasks]
        for future in as_completed(futures):
            print("  " + future.result(), flush=True)


def same_pairing(target, pipeline):
    """保存済みの組み方の材料が、今の枝の追い方の条件で作られたものか。"""
    with np.load(target) as saved:
        info = json.loads(str(saved["info"]))
    return all(info.get(key) == pipeline.pairing[key] for key in PAIRING_SCAN_KEYS)


def same_scan(target, pipeline):
    """保存済みファイルが今のスキャン条件で作られたものか。"""
    with np.load(target) as saved:
        info = json.loads(str(saved["info"]))
    return (info["radius_mm"] == pipeline.scan["radius_mm"]
            and info["scan_beam_radius_mm"] == pipeline.scan["beam_radius_mm"]
            and info["include_z_vetoed"] == pipeline.scan_includes_z_vetoed
            and info["fit"] == fit_constants())
