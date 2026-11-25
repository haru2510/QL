#!/usr/bin/env python3
"""共通JSONで準備→シンチ図→TPC段階別積算を順に実行する。"""
import argparse  # 設定と準備中プロセスを指定する。
import os  # 既に動いている準備処理の状態を確認する。
import subprocess  # 各既存スクリプトを独立したプロセスで実行する。
import sys  # 使用するPythonを統一する。
import time  # 実行中プロセスの完了を待つ。
from pathlib import Path  # リポジトリとログを指定する。
from settings import load_config, repo_path  # 全処理の入出力を共通にする。


def main():  # 各処理の成功を確認してから次へ進む。
    parser = argparse.ArgumentParser(description=__doc__)  # 実行条件を用意する。
    parser.add_argument('--config', required=True, type=Path)  # 全解析の外部JSONを指定する。
    parser.add_argument('--wait-preparation-pid', type=int)  # 今回既に開始した準備処理を重複起動しない。
    parser.add_argument('--after-workflow-log', type=Path)  # 前の全工程が正常完了した場合だけ次データを開始する。
    args = parser.parse_args()  # 引数を読む。
    config_path = args.config.resolve()  # 子プロセスへ絶対パスを渡す。
    config = load_config(config_path)  # 入出力とランを読む。
    root = repo_path(config['result_dir'])  # 結果の保存先を決める。
    logs = root / 'shared/workflow_logs'  # 処理ごとのログを残す。
    logs.mkdir(parents=True, exist_ok=True)  # 保存先を作る。
    if args.wait_preparation_pid:  # 既に動作している準備処理がある場合だけ待つ。
        while True:  # 初回処理の終了を確認する。
            try:  # プロセスの存在を調べる。
                os.kill(args.wait_preparation_pid, 0)  # シグナルを送らず生存確認する。
            except ProcessLookupError:  # 正常終了・異常終了のどちらでも待機を終える。
                break  # 再開時の検証へ進む。
            time.sleep(5)  # CPUを使わず待つ。
    if args.after_workflow_log and 'ALL STAGES COMPLETE' not in args.after_workflow_log.read_text():  # 待機した解析の正常終了を確認する。
        raise RuntimeError(f'Previous workflow did not complete: {args.after_workflow_log}')  # 失敗した前段を無視しない。
    commands = [('prepare_stages.py', []), ('stage_pixel_maps.py', []), ('pixel_polarization.py', ['--from-cache']), ('single_hit_pixels.py', []), ('stage_flag_quality_maps.py', []), ('stage_track_sum.py', []), ('stage_coincidence_xy.py', []), ('stage_quality_tracks.py', [])]  # キャッシュ確認を含めて全処理を順番に実行する。
    for script, options in commands:  # 前段が成功した場合だけ次を実行する。
        command = [sys.executable, '-u', str(Path(__file__).parent/script), '--config', str(config_path), *options]  # shell展開なしで引数を渡す。
        print('START', script, flush=True)  # 全体ログへ処理名を出す。
        with (logs/f'{Path(script).stem}.log').open('w') as stream:  # 詳細を専用ログへ出す。
            subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)  # 異常時はその場で停止する。
        print('COMPLETE', script, flush=True)  # 正常完了を記録する。
    print('ALL STAGES COMPLETE', root, flush=True)  # 全処理の完了を示す。


if __name__ == '__main__':  # importだけでは実行しない。
    main()  # コマンドから呼ばれた場合に実行する。
