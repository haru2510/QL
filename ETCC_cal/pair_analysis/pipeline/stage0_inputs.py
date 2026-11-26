"""stage0: 入力がそろっているかを確かめ、ランごとのトリガー数と Z 側除外の数を数える。

このパイプラインの入力（既存の ETCC_cal の出力。作り直すときは prepare_stages.py）:
  - 散乱点 CSV（ETCC 設定の inputs.vertices_all）: 同じ uPIC でアノードとカソードがそろう一番早いクロックの点。
  - Z 側除外（ETCC 設定の outputs.veto_dir/z0_track_veto_<run>.npz）: Z <= z_veto の raw ヒットを持つイベント番号。
  - raw TPC（ETCC 設定の data_root/storages/<run>）: stage2 と stage5 で読む。

トリガー数: TPC のイベントは epoch 0、trig 1..N なので、Z 側除外と散乱点に現れる最大の trig を N とみなす（近似）。
"""
import json

import numpy as np

from settings import repo_path


def run(pipeline):
    print("== stage0: inputs", flush=True)
    beam_ok = pipeline.beam_axis_file.exists()
    summary = {"beam_axis_fit": str(pipeline.beam_axis_file), "beam_axis_fit_exists": beam_ok, "datasets": {}}
    problems = [] if beam_ok else [f"missing beam axis fit: {pipeline.beam_axis_file}"]
    for name in pipeline.datasets:
        config = pipeline.etcc_config(name)
        vertices_csv = repo_path(config["inputs"]["vertices_all"])
        veto_dir = repo_path(config["outputs"]["veto_dir"])
        if not vertices_csv.exists():
            problems.append(f"{name}: missing vertices CSV {vertices_csv}")
        runs = {}
        for run_name in pipeline.datasets[name]["runs"]:
            veto_file = veto_dir / f"z0_track_veto_{run_name}.npz"
            raw_dirs = [repo_path(config["data_root"]) / storage / run_name for storage in config["storages"]]
            info = {"veto_file": str(veto_file), "raw_dirs_found": [str(d) for d in raw_dirs if d.exists()]}
            if not veto_file.exists():
                problems.append(f"{name} {run_name}: missing {veto_file}")
            else:
                vetoed = np.load(veto_file)["event"]
                info["z_vetoed_events"] = int(len(vetoed))
                info["max_vetoed_event"] = int(vetoed.max()) if len(vetoed) else 0
            if not info["raw_dirs_found"]:
                problems.append(f"{name} {run_name}: raw directory not found under {config['data_root']}")
            runs[run_name] = info
        summary["datasets"][name] = {"etcc_config": pipeline.datasets[name]["etcc_config"], "vertices_csv": str(vertices_csv),
                                     "z_veto_mm": config["z_veto"], "runs": runs}
    summary["problems"] = problems
    out = pipeline.stage_dir(0) / "inputs.json"
    out.write_text(json.dumps(summary, indent=1) + "\n")
    for problem in problems:
        print("  PROBLEM:", problem)
    print(f"  saved {out}", flush=True)
    if problems:
        raise SystemExit("stage0: inputs are missing (see above)")
