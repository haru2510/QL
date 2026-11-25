#!/usr/bin/env python3
"""水平・円のstage04積算が正常完了した後、保存配列からTPC比を作る。"""
import argparse  # 待機対象と比較JSONを指定する。
import json  # 出力先を読む。
import os  # 既存解析の終了だけを確認する。
import subprocess  # 比較コードを実行する。
import sys  # 同じPythonを使う。
import time  # 待機中のCPU負荷を避ける。
from pathlib import Path  # スクリプトとログを扱う。
from settings import repo_path  # リポジトリ基準のパスを使う。


def main():  # stage04の正常終了を待ってから処理する。
    parser = argparse.ArgumentParser(description=__doc__)  # 引数を用意する。
    parser.add_argument('--wait-pid',required=True,type=int)  # stage04ドライバの終了を待つ。
    parser.add_argument('--workflow-log',required=True,type=Path)  # 両データの完了マーカーを確認する。
    parser.add_argument('--config',required=True,type=Path)  # 比較の入力・出力・表示条件を指定する。
    args = parser.parse_args()  # 引数を読む。
    print('WAITING FOR STAGE04',args.wait_pid,flush=True)  # 起動済みの待機処理と示す。
    while True:  # 既存処理を止めずに待つ。
        try:  # プロセスの存在を調べる。
            os.kill(args.wait_pid,0)  # シグナルは送らない。
        except ProcessLookupError:  # 終了したら次へ進む。
            break  # 待機を終える。
        time.sleep(5)  # CPUを使わず待つ。
    if 'ALL STAGE04 FOLLOWUPS COMPLETE' not in args.workflow_log.read_text():  # 前段失敗を無視しない。
        raise RuntimeError(f'Stage04 did not finish: {args.workflow_log}')  # 部分配列で比較しない。
    config = json.loads(args.config.read_text())  # 比較の保存先を読む。
    output = repo_path(config['output_dir'])  # シンチ比較と別フォルダにする。
    output.mkdir(parents=True,exist_ok=True)  # 初回に保存先を作る。
    print('START STAGE04 TPC RATIO',flush=True)  # 実処理の開始を記録する。
    with (output/'workflow.log').open('w') as stream:  # 描画と集計の詳細を保存する。
        subprocess.run([sys.executable,'-u',str(Path(__file__).parent/'compare_tpc_tracks.py'),'--config',str(args.config.resolve())],check=True,stdout=stream,stderr=subprocess.STDOUT)  # SSDのraw再走査はしない。
    print('ALL STAGE04 TPC RATIOS COMPLETE',output,flush=True)  # 比較の完了を記録する。


if __name__=='__main__':  # importだけでは実行しない。
    main()  # コマンドから呼ばれた場合だけ処理する。
