#!/usr/bin/env python3
"""各解析で共通に使うJSON設定と入出力パス。相対パスはリポジトリ直下が基準。"""
import json  # 外部JSONを読み込む。
from pathlib import Path  # 入出力パスを安全に扱う。

REPO = Path(__file__).resolve().parents[1]  # コードの場所からリポジトリを求める。
DEFAULT_CONFIG = Path(__file__).parent / "config/20261124_run06_07.json"  # 引数省略時にも外部JSONを使う。


def repo_path(value):  # 設定内の相対パスをリポジトリ基準に解決する。
    path = Path(value)  # 文字列をパスに変換する。
    return path if path.is_absolute() else REPO / path  # 絶対パスならそのまま使う。


def load_config(path):  # 全解析で同じ形式のJSONを読み込む。
    config = json.loads(Path(path).read_text())  # 外部JSONを辞書にする。
    if config["velocity"] <= 0:  # クロックからZへの変換が単調増加か確認する。
        raise ValueError("velocity must be positive [mm/clock]")  # 不正な較正値で解析しない。
    if config["pitch_strip"] != 0.4:  # 現在のストリップ配置の前提を確認する。
        raise ValueError("This strip geometry requires pitch_strip = 0.4 mm")  # 配置と積算ビンの不一致を防ぐ。
    return config  # 読み取った設定を呼び出し側へ返す。


def input_file(config, key):  # 入力ファイルをJSONから選ぶ。
    return repo_path(config["inputs"][key])  # 入力欄のパスをリポジトリ基準に変換する。


def output_file(config, key):  # 出力ファイルをJSONから選ぶ。
    path = repo_path(config["outputs"][key])  # 出力欄のパスをリポジトリ基準に変換する。
    path.parent.mkdir(parents=True, exist_ok=True)  # ファイルの親フォルダを必要なら作る。
    return path  # 解決済みの保存先を返す。
