#!/usr/bin/env python3
"""実行中の解析が正常終了した後、複数データセットへstage04の品質選別後TPC積算を追加する。"""
import argparse  # 待機対象と設定を指定する。
import os  # 実行中解析の生存確認をする。
import subprocess  # XY追加スクリプトを順に実行する。
import sys  # 同じPythonを使う。
import time  # CPUを使わず前段を待つ。
from pathlib import Path  # スクリプトとログを扱う。
from settings import load_config, repo_path  # 共通の出力先を使う。


def main():  # 前段成功後にだけ追加解析を開始する。
    parser = argparse.ArgumentParser(description=__doc__)  # 指定内容を用意する。
    parser.add_argument('--wait-pid',required=True,type=int)  # 既存の円偏光解析を待つ。
    parser.add_argument('--workflow-log',required=True,type=Path)  # 正常終了マーカーを確認するログを指定する。
    parser.add_argument('--configs',required=True,nargs='+',type=Path)  # 順に処理する外部JSONを指定する。
    args = parser.parse_args()  # 引数を読む。
    print('WAITING FOR WORKFLOW',args.wait_pid,flush=True)  # 待機中も起動済みであることを記録する。
    while True:  # 前段の終了まで待つ。
        try:  # プロセスの有無だけを確認する。
            os.kill(args.wait_pid,0)  # シグナルは送らない。
        except ProcessLookupError:  # 終了したら次の確認へ進む。
            break  # 待機を終える。
        time.sleep(5)  # CPU負荷を避ける。
    if 'ALL XY FOLLOWUPS COMPLETE' not in args.workflow_log.read_text():  # 前段が失敗した場合は続行しない。
        raise RuntimeError(f'Previous workflow failed: {args.workflow_log}')  # エラーを明示する。
    for path in args.configs:  # データセットを混ぜずに順に処理する。
        config = load_config(path)  # 入出力を読む。
        log = repo_path(config['result_dir'])/'shared/workflow_logs/stage_quality_tracks.log'  # 各データセットへ詳細を残す。
        log.parent.mkdir(parents=True,exist_ok=True)  # 保存先を作る。
        print('START STAGE04',path,flush=True)  # 全体ログに実行状況を残す。
        with log.open('w') as stream:  # rawごとの詳細を分けて保存する。
            subprocess.run([sys.executable,'-u',str(Path(__file__).parent/'stage_quality_tracks.py'),'--config',str(path.resolve())],check=True,stdout=stream,stderr=subprocess.STDOUT)  # 失敗時は次へ進まない。
        print('COMPLETE STAGE04',path,flush=True)  # 正常完了を示す。
    print('ALL STAGE04 FOLLOWUPS COMPLETE',flush=True)  # 全体完了を記録する。


if __name__=='__main__':  # importだけでは実行しない。
    main()  # コマンドから呼ばれた場合だけ実行する。
