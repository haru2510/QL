#!/usr/bin/env python3
"""Z側除外の前後で、1イベントのシンチ（HA）ヒットフラグ数の分布を円グラフとヒストグラムで描く。

入力は select_vertices.py が保存した selection/ha_hit_flags_<run>.npz と z0_track_veto_<run>.npz、
全イベント数は track_sum.json。HA rawは読まない。ヒットはヘッダのDataCountのビットで、HG・時刻のカットはしない。
"""
import argparse  # 設定ファイルのパスを受け取る。
import json  # 全イベント数を読む。
from pathlib import Path  # 入出力パスを扱う。
import numpy as np  # イベントごとのヒット数を数える。
import matplotlib  # 保存用の描画環境を用意する。
matplotlib.use("Agg")  # GUIを使わずPNGを保存する。
import matplotlib.pyplot as plt  # 円グラフとヒストグラムを描く。
from settings import DEFAULT_CONFIG, load_config, repo_path, output_file  # 入出力は他のrun06・07コードと同じJSONを使う。

COLORS = ["#b8b7b2", "#2a78d6", "#eb6834", "#1baf7a"]  # 0、1、2、3以上（0は灰色、残りは固定順のカテゴリ色）。
LABELS = ["0 hits", "1 hit", "2 hits", "3+ hits"]  # 円グラフの区分。
MAX_HITS = 12  # ヒストグラムの最後のビン（これ以上はまとめる）。


def hit_counts(name, config, all_events):  # 1ランの全イベントと除外後イベントのヒット数の配列を作る。
    selection = repo_path(config["outputs"]["veto_dir"])  # キャッシュの場所。
    data = np.load(selection / f"ha_hit_flags_{name}.npz")  # HAのイベント番号とヒット数。
    veto = np.load(selection / f"z0_track_veto_{name}.npz")["event"]  # Z側除外のイベント番号。
    counts = np.zeros(all_events + 1, np.int64)  # trig 1..N（epoch 0）ごとのヒット数。HAに記録が無ければ0。
    event, n_hits = data["event"], data["n_hits"]  # 保存済みの値を取り出す。
    inside = (event >= 1) & (event <= all_events)  # TPCのイベントに対応するものだけを使う。
    counts[event[inside]] = n_hits[inside]  # ヒット数を入れる。
    keep = np.ones(all_events + 1, bool)  # 除外を通るイベント。
    keep[0] = False  # trig 0 は存在しない。
    keep[veto] = False  # Z<=15 mmのヒットがあるイベントを除く。
    return counts[1:], counts[keep]  # 全イベントと除外後イベントのヒット数。


def main():  # 保存済みの集計から図を作る。
    parser = argparse.ArgumentParser(description=__doc__)  # 引数の説明を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 共通の解析条件を読む。
    args = parser.parse_args()  # 指定を読む。
    config = load_config(args.config)  # 共通JSONを読む。
    scan = json.loads(output_file(config, "track_metadata").read_text())  # 全イベント数を読む。
    totals = {run["run"]: run["all_events"] for run in scan["runs"]}  # ランごとの全イベント数。
    before, after = [], []  # 全ランを足した配列を作る。
    for name in config["runs"]:  # ランを混ぜずに対応づけてから足す。
        a, b = hit_counts(name, config, totals[name])  # このランの前後のヒット数。
        before.append(a)  # 除外前に加える。
        after.append(b)  # 除外後に加える。
    sets = [("All events", np.concatenate(before)), (f"After Z <= {config['z_veto']:g} mm veto", np.concatenate(after))]  # 比べる2つ。

    fig = plt.figure(figsize=(13, 9.5), constrained_layout=True)  # 上段に円グラフ2つ、下段にヒストグラム。
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])  # 2行2列。
    for col, (title, n) in enumerate(sets):  # 前後それぞれの円グラフ。
        ax = fig.add_subplot(grid[0, col])  # 上段の枠。
        parts = [np.sum(n == 0), np.sum(n == 1), np.sum(n == 2), np.sum(n >= 3)]  # 0、1、2、3以上のイベント数。
        wedges, _ = ax.pie(parts, colors=COLORS, startangle=90, counterclock=False, wedgeprops=dict(edgecolor="white", linewidth=2))  # 白い隙間で区分を分ける。
        for wedge, count in zip(wedges, parts):  # 各区分に割合と件数を書く。
            if count / len(n) < 0.03:  # 細すぎる区分は凡例に任せる。
                continue  # 文字を重ねない。
            angle = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)  # 区分の中央の角度。
            ax.text(0.68 * np.cos(angle), 0.68 * np.sin(angle), f"{100 * count / len(n):.1f}%\n({count:,})", ha="center", va="center", fontsize=11, color="white", fontweight="bold")  # 割合と件数。
        ax.set_title(f"{title}\n{len(n):,} events, mean {n.mean():.2f} hits/event", fontsize=12)  # イベント数と平均。
        ax.legend(wedges, [f"{label}: {100 * c / len(n):.1f}%" for label, c in zip(LABELS, parts)], loc="center left", bbox_to_anchor=(0.92, 0.5), frameon=False, fontsize=10)  # 0ヒットなど細い区分も凡例で読めるようにする。
        ax.set_aspect("equal")  # 円を正円にする。

    ax = fig.add_subplot(grid[1, :])  # 下段のヒストグラム。
    bins = np.arange(MAX_HITS + 1)  # 0..MAX_HITS（最後は以上）。
    width = 0.4  # 2本の棒の幅。
    for k, ((title, n), color) in enumerate(zip(sets, ["#86b6ef", "#1c5cab"])):  # 前を淡い青、後を濃い青。
        h = np.array([np.sum(np.minimum(n, MAX_HITS) == b) for b in bins]) / len(n) * 100  # イベントの割合 [%]。
        ax.bar(bins + (k - 0.5) * width, h, width=width - 0.04, color=color, label=title)  # 前後を並べる。
    ax.set_xticks(bins, [str(b) for b in bins[:-1]] + [f"{MAX_HITS}+"])  # 最後のビンは以上。
    ax.set(xlabel="HA hit flags per event (all boards)", ylabel="Fraction of events [%]")  # 軸。
    ax.grid(axis="y", color="#e4e3df", lw=0.8)  # 控えめな目盛り線。
    ax.set_axisbelow(True)  # 線を棒の後ろに置く。
    for side in ("top", "right"):  # 不要な枠を消す。
        ax.spines[side].set_visible(False)  # 枠線を隠す。
    ax.legend(frameon=False)  # 前後の凡例。

    fig.suptitle(", ".join(config["runs"]) + ": HA hit multiplicity before / after the TPC Z-side veto\n"
                 "Hit = DataCount bit in the HA frame header (no HG / peak_id / timing cut); selection = no raw TPC hit at Z <= "
                 f"{config['z_veto']:g} mm only", fontsize=12)  # 何を数えたかと選別の条件。
    destination = repo_path(config["outputs"]["veto_dir"]) / "scint_hit_multiplicity.png"  # summary.jsonの隣に保存する。
    fig.savefig(destination, dpi=150)  # 図を保存する。
    plt.close(fig)  # 描画メモリを解放する。
    print(f"Saved {destination}")  # 保存先を知らせる。


if __name__ == "__main__":  # importだけでは描かない。
    main()  # コマンド実行時にだけ処理する。
