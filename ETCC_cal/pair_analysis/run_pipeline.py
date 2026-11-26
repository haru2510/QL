#!/usr/bin/env python3
"""TPC の飛跡による偏光解析を、設定JSONに従って stage0 から stage5 まで順に実行する。

  stage0  入力の確認（散乱点 CSV・Z 側除外・raw の場所・ビーム軸）
  stage1  散乱点に Z 側除外の印とビーム軸からの距離を付ける
  stage2  raw TPC を読み、散乱点ごとに XZ・YZ の1本fit・2本fit（時間がかかるのはここだけ。結果は再利用）
  stage3  cut set ごとにカットをかけて分類（SS / DS / SD / DD）し、カットの流れの表を作る
  stage4  比較ごとにモジュレーション（対生成: Q と DD の ψ、1本の飛跡: φ）と有意度
  stage5  分類ごとのイベントの絵

例:
  python3 ETCC_cal/pair_analysis/run_pipeline.py --config ETCC_cal/pair_analysis/config/20261124_25_polarization.json
  python3 ETCC_cal/pair_analysis/run_pipeline.py --config ... --stages 3 4 --cut-sets default span15
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # ETCC_cal の共通モジュール（settings など）を使う

from pipeline import stage0_inputs, stage1_vertices, stage2_branch_fit, stage3_classify, stage4_modulation, stage5_event_display
from pipeline.config import Pipeline


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, required=True, help="パイプラインの設定JSON")
    parser.add_argument("--stages", type=int, nargs="+", help="実行するステージ（既定: 設定の run_stages、無ければ 0〜5）")
    parser.add_argument("--cut-sets", nargs="+", help="stage3〜5 で使う cut set（既定: 設定の全部）")
    parser.add_argument("--rebuild", action="store_true", help="stage1・stage2 の保存済み結果を使わずに作り直す")
    args = parser.parse_args()

    pipeline = Pipeline(args.config)
    stages = args.stages or pipeline.raw.get("run_stages", [0, 1, 2, 3, 4, 5])
    print(f"pipeline '{pipeline.name}' -> {pipeline.output_dir}")
    print(f"stages {stages}, cut sets {args.cut_sets or list(pipeline.cut_sets)}", flush=True)
    pipeline.output_dir.mkdir(parents=True, exist_ok=True)
    (pipeline.output_dir / "config_used.json").write_text(pipeline.path.read_text())  # 使った設定を結果の隣に残す

    steps = {
        0: lambda: stage0_inputs.run(pipeline),
        1: lambda: stage1_vertices.run(pipeline, rebuild=args.rebuild),
        2: lambda: stage2_branch_fit.run(pipeline, rescan=args.rebuild),
        3: lambda: stage3_classify.run(pipeline, args.cut_sets),
        4: lambda: stage4_modulation.run(pipeline, args.cut_sets),
        5: lambda: stage5_event_display.run(pipeline, args.cut_sets),
    }
    for stage in sorted(stages):
        start = time.time()
        steps[stage]()
        print(f"   stage{stage} done in {time.time() - start:.0f} s", flush=True)


if __name__ == "__main__":
    main()
