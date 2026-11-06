#!/usr/bin/env python3
#
# zeno から届く間引きイベントのうち HA の枠を受け、check_data.py と同じ確認図を積算しながら更新する。
#
# 受ける形式は ql_tpc.py と同じで、HA の枠は側 'h'、番号は HA 番号 (ha03 なら 3)。
# TPC の枠 (側 'a', 'c') は読み飛ばすので、ql_tpc.py と同じ送信元につないでよい。
# zeno が無いときは ql_replay.py が同じ形式で送る。起動するとブラウザが開き、図は自動で更新される。
#
# ブラウザには次の 2 枚を、モニターに収まるように並べる。
#   バタフライマップ : 全 HA を、底面 (ha00-05) を中央に、側面 (N: ha06, 07 / E: ha15, 16 /
#                      S: ha12, 13 / W: ha09, 10) をその上下左右に置いた十字形に並べる。
#                      HA ごとに 6 ch のバタフライマップを描き、ヒットの無い ch は赤く塗る。
#                      HA ごとの直近のトリガレートを書き、レートが高いほどタイルを濃いオレンジにする
#   ペデスタル       : HA と ch ごとの HG0-3 / LG0-3 のペデスタル (中央値と標準偏差)。
#                      面ごとに色帯を付け、バタフライマップと同じ色で並べる
# 別のページ (/temp) には MPPC の温度を出す。
#   上の段 : HA ごとの温度 (6 センサの平均) の時間変化。run が替わっても消さず、直近 12 時間を持つ
#   下の段 : HA ごとの今の温度と、6 センサの幅 (最小-最大)
# 温度は psa_param の HA ごとの係数で較正する。係数が ha02 の写し (CAL_NOTE が placeholder) の HA には
# "*" を付ける。psa_param に無い HA は psa_decode の既定の係数を使い、"?" を付ける。
#
# 解析は psa_decode.py の decode_events をそのまま使う。
#

import argparse
import asyncio
import csv
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

HERE = Path(__file__).resolve().parent

# psa_decode.py は analysis_sm3_FM. から始まる名前で sub_scr を読むので、その親も入れる
sys.path[:0] = [str(HERE), str(HERE.parents[2]), str(HERE.parents[1] / "TPC_analysis")]

from check_data import (  # noqa: E402
    N_BINS_X,
    N_BINS_Y,
    N_CHANNELS,
    X_MAX,
    X_MIN,
    Y_MAX,
    Y_MIN,
)
from psa_decode import (  # noqa: E402
    DEFAULT_TEMP_A,
    DEFAULT_TEMP_B,
    decode_events,
    get_sn_from_filename,
    get_unit_by_sn,
)
from sub_scr.get_info import get_temp_calibration, load_psa_config  # noqa: E402
from ql_tpc import (  # noqa: E402
    RUN_RESET_TRIGGERS,
    prepare_view,
    receive_loop,
    save_figure,
)

# このリポジトリ (analysis_sm3_FM) の場所。既定のパスはここからの相対で決める
REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


DEFAULT_OUTPUT_DIR = os.path.join(REPO_DIR, "check_output")
DEFAULT_PSA_PARAM = os.path.join(REPO_DIR, "PSA_analysis", "param", "psa_param_20261029.csv")

# HA の番号の範囲 (ha00-ha17)
N_HA = 18

# 温度の時間変化は、この秒数ごとの平均にまとめて、直近 TEMP_KEEP 秒を持つ
TEMP_BIN = 10.0
TEMP_KEEP = 12 * 3600.0

# HA の clk は 1 MHz の 32 bit カウンタ
HA_CLOCK = 1e6
HA_CLOCK_WRAP = 2 ** 32


SIGNALS = ["HG0", "HG1", "HG2", "HG3", "LG0", "LG1", "LG2", "LG3"]

# HA の配置。一覧の図は、底面を中央に、側面をその上下左右に置いた十字形に並べる。
# 同じ面の中は番号順。(面の名前, 色, [HA 番号])
HA_GROUPS = [
    ("Bottom", "#eef3fb", [0, 1, 2, 3, 4, 5]),
    ("N side", "#fdf0e6", [6, 7]),
    ("E side", "#eaf6ea", [15, 16]),
    ("S side", "#fbeaea", [12, 13]),
    ("W side", "#f2ecf8", [9, 10]),
]
GROUP_PLACES = {           # 十字形の 4 行 x 5 列のうち、面を置く場所 (行, 列) と、面の中の並び (行数, 列数)
    "Bottom": ((slice(1, 3), slice(1, 4)), (2, 3)),
    "N side": ((0, slice(1, 4)), (1, 2)),
    "S side": ((3, slice(1, 4)), (1, 2)),
    "W side": ((slice(1, 3), 0), (2, 1)),
    "E side": ((slice(1, 3), 4), (2, 1)),
}
# HA_GROUPS に無い HA (未接続の ha08 などがつながったとき) を置く、十字形の四隅
SPARE_PLACES = [(0, 0), (0, 4), (3, 0), (3, 4)]

# トリガレートに応じてタイルの背景を塗る色
RATE_COLORS = colormaps["Oranges"]

# ch あたり平均でこれだけヒットがたまってから、ヒットの無い ch を赤く塗る。少ないうちは偶然 0 になる
DEAD_CH_HITS = 10

# 一覧の図のトリガレートは、直近のこの秒数で求める
RATE_WINDOW = 60.0

# 一覧の図は 1 枚が小さいので、ビンをさらに粗くする。細かいままだと 1 ビンが 1 画素より小さくなる
OVERVIEW_BINS = (N_BINS_X // 5, N_BINS_Y // 5)

# 左にバタフライマップ、右にペデスタルを並べる
PAGE = """<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<title>PSA QL</title>
<style>
html, body { margin: 0; height: 100%; overflow: hidden; background: #fff; }
body { display: grid; height: 100vh; grid-template-columns: 62fr 38fr; }
img { width: 100%; height: 100%; min-height: 0; object-fit: contain; }
#status { position: fixed; right: 8px; bottom: 4px; font: 12px sans-serif; color: #555; }
</style>
</head>
<body>
<img id="overview" alt="">
<img id="pedestal" alt="">
<div id="status">図を待っています</div>
<script>
const statusLine = document.getElementById("status");
function watch(name, interval) {
  const figure = document.getElementById(name);
  function refresh() {
    const next = new Image();
    next.onload = () => {
      figure.src = next.src;
      statusLine.textContent = "読み込み " + new Date().toLocaleTimeString();
    };
    next.src = "/" + name + ".png?t=" + Date.now();
  }
  refresh();
  setInterval(refresh, interval);
}
watch("overview", SUMMARY_MS);
watch("pedestal", SUMMARY_MS);
</script>
</body>
</html>
"""

# MPPC の温度は別のページ (/temp) に出す
TEMP_PAGE = """<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<title>PSA QL temperature</title>
<style>
html, body { margin: 0; height: 100%; overflow: hidden; background: #fff; }
img { display: block; width: 100%; height: 100vh; object-fit: contain; }
#status { position: fixed; right: 8px; bottom: 4px; font: 12px sans-serif; color: #555; }
</style>
</head>
<body>
<img id="temperature" alt="">
<div id="status">図を待っています</div>
<script>
const statusLine = document.getElementById("status");
const figure = document.getElementById("temperature");
function refresh() {
  const next = new Image();
  next.onload = () => {
    figure.src = next.src;
    statusLine.textContent = "読み込み " + new Date().toLocaleTimeString();
  };
  next.src = "/temperature.png?t=" + Date.now();
}
refresh();
setInterval(refresh, SUMMARY_MS);
</script>
</body>
</html>
"""


def parse_args():

    parser = argparse.ArgumentParser(description="PSA (HA) quick look for thinned events from zeno")
    parser.add_argument("--host", default="zeno", help="host sending the thinned events")
    parser.add_argument("--port", type=int, default=10012, help="TCP port of the thinned events")
    parser.add_argument("--ha", type=int, nargs="+", default=None,
                        help="HA numbers to show in the accumulated figure (default: all received)")
    parser.add_argument("--update", type=float, default=10.0,
                        help="update interval of the accumulated figures [s]")
    parser.add_argument("--reconnect", type=float, default=3.0, help="reconnect interval [s]")
    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--psa-param", default=DEFAULT_PSA_PARAM,
                        help="psa_param csv with the temperature calibration of each HA")
    parser.add_argument("--name", default="ql_psa", help="name of the output PNG")
    parser.add_argument("--http-bind", default="127.0.0.1", help="address of the viewer")
    parser.add_argument("--http-port", type=int, default=8091, help="TCP port of the viewer")
    parser.add_argument("--no-open", action="store_true", help="do not open a browser")

    return parser.parse_args()


def load_temp_calibrations(path):
    """psa_param から、HA 番号ごとの温度の係数 (temp_a, temp_b) と印を返す。

    印は "" (較正済み)、"*" (係数がほかの HA の写し)、"?" (psa_param に無く、既定の係数を使う)。
    """

    calibrations = {}

    try:
        psa_config = load_psa_config({"psa_param_file": path})
        with open(path, encoding="utf-8") as f:
            notes = {row.get("HA", "").strip().lower(): row.get("CAL_NOTE", "") for row in csv.DictReader(f)}
    except (OSError, KeyError, ValueError) as exc:
        print(f"[WARN] cannot read {path}: {exc}. default temperature calibration is used", flush=True)
        psa_config, notes = [], {}

    for index in range(N_HA):

        name = f"ha{index:02d}"
        unit = None
        if any(unit_config.get("HA", "").lower() == name for unit_config in psa_config):
            unit = get_unit_by_sn(get_sn_from_filename(f"{name}_00000.dat", psa_config), psa_config)

        if unit is None:
            calibrations[index] = (DEFAULT_TEMP_A, DEFAULT_TEMP_B, "?")
            continue

        temp_a, temp_b = get_temp_calibration(unit)
        mark = "*" if "placeholder" in notes.get(name, "").lower() else ""
        calibrations[index] = (np.ascontiguousarray(temp_a, dtype=np.float32),
                               np.ascontiguousarray(temp_b, dtype=np.float32), mark)

    return calibrations


def decode_frame(frame, temp_a=DEFAULT_TEMP_A, temp_b=DEFAULT_TEMP_B):
    """1 枠を psa_decode の decode_events で読み、ヒットした ch ごとの値を返す。"""

    buf = np.frombuffer(frame, dtype=np.uint8).reshape(1, -1)
    sizes = np.array([len(frame)], dtype=np.int32)

    (trig, clock, hit, _, temp, ped, _, hg, _, x_hg, y_hg, _, _,
     peak_id, _, satu_h, _, adc) = decode_events(buf, sizes, False, temp_a, temp_b, True)

    hits = []

    # psa_decode.py と同じく、ch 番号は 5 - (データ中の並び)
    for index in np.flatnonzero(hit[0] == 1):
        hits.append({
            "ch": N_CHANNELS - 1 - int(index),
            "hg": float(hg[0, index]),
            "x_hg": float(x_hg[0, index]),
            "y_hg": float(y_hg[0, index]),
            "peak_id": int(peak_id[0, index]),
            "satu_h": int(satu_h[0, index]),
            "ped": ped[0, index].copy(),
            "wave": adc[0, index, :, 2:].copy(),
        })

    return int(trig[0]), int(clock[0]), temp[0].copy(), hits


class HaState:
    """1 台の HA の積算。"""

    def __init__(self):

        self.trig = []
        self.seconds = []
        self.temp = []
        self.hits = {name: [] for name in ("ch", "hg", "x_hg", "y_hg", "peak_id", "satu_h")}
        self.ped = []
        self.previous_clock = None
        self.elapsed = 0

    def add(self, trig, clock, temp, hits):

        # clk の 1 周 (約 72 分) をまたいでも時刻が戻らないようにする
        if self.previous_clock is not None:
            self.elapsed += (clock - self.previous_clock) % HA_CLOCK_WRAP
        self.previous_clock = clock

        self.trig.append(trig)
        self.seconds.append(self.elapsed / HA_CLOCK)
        self.temp.append(temp)

        for hit in hits:
            for name in self.hits:
                self.hits[name].append(hit[name])
            self.ped.append(hit["ped"])

    def snapshot(self):

        state = {name: np.array(values) for name, values in self.hits.items()}
        state["ped"] = np.array(self.ped, dtype=np.float32).reshape(-1, len(SIGNALS))
        state["trig"] = np.array(self.trig, dtype=np.int64)
        state["seconds"] = np.array(self.seconds, dtype=np.float64)
        state["temp"] = np.array(self.temp, dtype=np.float32).reshape(-1, 6)

        return state


class TempHistory:
    """1 台の HA の温度の時間変化。TEMP_BIN 秒ごとに、6 センサそれぞれの和とイベント数をまとめる。"""

    def __init__(self):

        self.bins = {}

    def add(self, now, temp):

        key = int(now // TEMP_BIN)
        total, count = self.bins.get(key, (np.zeros(6), 0))
        self.bins[key] = (total + temp, count + 1)

        # 新しいビンができたときだけ、TEMP_KEEP より古いビンを捨てる
        if count == 0:
            oldest = int((now - TEMP_KEEP) // TEMP_BIN)
            for old in [k for k in self.bins if k < oldest]:
                del self.bins[old]

    def snapshot(self):
        """(時刻 [unix s], 6 センサの平均, センサごとの平均 [n, 6]) を時刻の順に返す。"""

        keys = sorted(self.bins)
        times = np.array([(k + 0.5) * TEMP_BIN for k in keys])
        sensors = np.array([self.bins[k][0] / self.bins[k][1] for k in keys]).reshape(-1, 6)

        return times, sensors.mean(axis=1), sensors


def recent_temp(state):
    """直近 RATE_WINDOW 秒の温度 (6 センサの平均) [℃]。"""

    if state is None or len(state["temp"]) == 0:
        return None

    selected = state["seconds"] >= state["seconds"][-1] - RATE_WINDOW

    return float(state["temp"][selected].mean())


def butterfly_counts(state, ch, bins):
    """1 ch のバタフライマップのヒストグラムと、その ch のヒット数を返す。"""

    selected = state["ch"] == ch if len(state["ch"]) else np.zeros(0, dtype=bool)
    n_hits = int(selected.sum())
    counts, x_edges, y_edges = np.histogram2d(
        state["x_hg"][selected] if n_hits else [], state["y_hg"][selected] if n_hits else [],
        bins=bins, range=[[X_MIN, X_MAX], [Y_MIN, Y_MAX]])

    return counts, x_edges, y_edges, n_hits


def recent_rate(state):
    """直近 RATE_WINDOW 秒のトリガレート [Hz]。間引く前の trig の増え方を HA の clk の時間で割る。"""

    seconds, trig = state["seconds"], state["trig"]

    if len(trig) < 2:
        return None

    selected = seconds >= seconds[-1] - RATE_WINDOW

    if selected.sum() < 2:
        selected = slice(-2, None)

    span = seconds[selected][-1] - seconds[selected][0]

    return (trig[selected][-1] - trig[selected][0]) / span if span > 0 else None


def draw_ha_tile(tile, index, state, max_rate, temp_mark=""):
    """一覧の図の 1 台分。6 ch のバタフライマップを 2 行 3 列で小さく描き、ヒットの無い ch は赤く塗る。

    どの HA に集中しているかが分かるよう、トリガレートを題に書き、全 HA で最も高いレートに対する
    割合でタイルの背景を濃く塗る。
    """

    n_hits = 0 if state is None else len(state["ch"])
    rate = None if state is None else recent_rate(state)
    temp = recent_temp(state)
    temp_text = "" if temp is None else f"  {temp:.1f}\u00b0C{temp_mark}"

    if state is None:
        label = f"ha{index:02d}  no data"
    elif rate is None:
        label = f"ha{index:02d}  -- Hz{temp_text}  ({n_hits} hits)"
    else:
        label = f"ha{index:02d}  {rate:.3g} Hz{temp_text}  ({n_hits} hits)"
        if max_rate:
            tile.set_facecolor(RATE_COLORS(0.05 + 0.6 * rate / max_rate))

    tile.suptitle(label, fontsize=9, fontweight="bold")
    axes = tile.subplots(2, 3).ravel()

    for ch, ax in enumerate(axes[:N_CHANNELS]):

        ax.set_xticks([])
        ax.set_yticks([])

        if state is None:
            ax.set_facecolor("#e6e6e6")
            continue

        counts, _, _, n_ch = butterfly_counts(state, ch, OVERVIEW_BINS)

        if counts.max() > 0:
            ax.imshow(np.ma.masked_equal(counts.T, 0), origin="lower", aspect="auto",
                      extent=(X_MIN, X_MAX, Y_MIN, Y_MAX), norm=LogNorm(vmin=1),
                      cmap="viridis", interpolation="nearest")
        elif n_hits >= DEAD_CH_HITS * N_CHANNELS:
            # ほかの ch には十分ヒットがあるのに、この ch だけ無い
            ax.set_facecolor("#f4b6b6")

        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.set_title(f"ch{ch}  {n_ch}", fontsize=7, pad=1)


def draw_overview_figure(states, title, out_path, marks=None):
    """全 HA を、底面を中央に、側面をその上下左右に置いた十字形に並べる。

    配置にある HA は、データがまだ無くても灰色の枠を出す。配置に無い HA は四隅に置く。
    """

    rates = [rate for rate in (recent_rate(state) for state in states.values()) if rate is not None]
    max_rate = max(rates) if rates else None

    fig = Figure(figsize=(14, 10.5), constrained_layout=True)
    grid = fig.add_gridspec(4, 5)

    for label, color, members in HA_GROUPS:

        place, (n_rows, n_cols) = GROUP_PLACES[label]
        group = fig.add_subfigure(grid[place], facecolor=color, edgecolor="0.6", linewidth=1)
        group.suptitle(label, fontsize=10, fontweight="bold", x=0.01, ha="left")
        tiles = np.atleast_1d(group.subfigures(n_rows, n_cols, wspace=0.03, hspace=0.03)).ravel()

        for tile, index in zip(tiles, members):
            draw_ha_tile(tile, index, states.get(index), max_rate, (marks or {}).get(index, ""))

    grouped = {index for _, _, members in HA_GROUPS for index in members}
    spares = [index for index in states if index not in grouped]

    for place, index in zip(SPARE_PLACES, spares):
        tile = fig.add_subfigure(grid[place], facecolor="#f2f2f2", edgecolor="0.6", linewidth=1)
        draw_ha_tile(tile, index, states[index], max_rate, (marks or {}).get(index, ""))

    fig.suptitle(title, fontsize=11)

    save_figure(fig, out_path)


def draw_pedestal_figure(states, title, out_path):
    """HA と ch ごとのペデスタルを、HG0-3 と LG0-3 に分けて描く。

    並びはバタフライマップの面の順 (Bottom, N, E, S, W)。配置にある HA は、データが無くても
    場所を空けておく。配置に無い HA は最後に足す。
    """

    grouped = [(label, color, members) for label, color, members in HA_GROUPS]
    spares = [index for index in states if all(index not in members for _, _, members in HA_GROUPS)]
    if spares:
        grouped.append(("other", "#f2f2f2", spares))

    fig = Figure(figsize=(9, 10.5), constrained_layout=True)
    ax_hg, ax_lg = fig.subplots(2, 1)

    positions = {}
    for label, color, members in grouped:
        for index in members:
            for ch in range(N_CHANNELS):
                positions[(index, ch)] = len(positions)

    for ax, signals, name in ((ax_hg, range(4), "HG"), (ax_lg, range(4, 8), "LG")):

        # 面ごとの色帯と、HA の境目の線
        for label, color, members in grouped:
            first = positions[(members[0], 0)] - 0.5
            last = positions[(members[-1], N_CHANNELS - 1)] + 0.5
            ax.axvspan(first, last, color=color, zorder=0)
            ax.text((first + last) / 2, 1.01, label, transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=9, fontweight="bold")
            for index in members[1:]:
                ax.axvline(positions[(index, 0)] - 0.5, color="0.75", linewidth=0.6, zorder=1)

        for offset, signal in enumerate(signals):

            xs, medians, spreads = [], [], []

            for index, state in states.items():
                for ch in range(N_CHANNELS):
                    selected = state["ch"] == ch if len(state["ch"]) else np.zeros(0, dtype=bool)
                    if (index, ch) not in positions or not selected.any():
                        continue
                    ped = state["ped"][selected, signal]
                    xs.append(positions[(index, ch)] + (offset - 1.5) * 0.15)
                    medians.append(np.median(ped))
                    spreads.append(ped.std())

            if xs:
                ax.errorbar(xs, medians, yerr=spreads, fmt="o", markersize=2.5, capsize=1.5,
                            linewidth=0.8, color=f"C{offset}", label=SIGNALS[signal], zorder=2)

        ticks = [positions[(index, 0)] + (N_CHANNELS - 1) / 2 for _, _, members in grouped for index in members]
        names = [f"ha{index:02d}" for _, _, members in grouped for index in members]
        ax.set_xticks(ticks, names, fontsize=8, rotation=90)
        ax.set_xlim(-0.5, len(positions) - 0.5)
        ax.set_ylabel(f"{name} pedestal [ADC]")
        ax.set_title(f"{name} pedestal (median, std; ch0-5 left to right in each HA)", fontsize=10,
                     pad=16)
        if any(len(state["ch"]) for state in states.values()):
            ax.legend(fontsize=7, ncol=4, loc="lower right")

    fig.suptitle(title, fontsize=10)

    save_figure(fig, out_path)


def draw_temperature_figure(histories, marks, title, out_path):
    """MPPC の温度。上の段は HA ごとの時間変化 (6 センサの平均)、下の段は今の温度と 6 センサの幅。

    並びと色の帯はバタフライマップの面の順 (Bottom, N, E, S, W)。
    """

    grouped = [(label, color, members) for label, color, members in HA_GROUPS]
    spares = [index for index in histories if all(index not in members for _, _, members in HA_GROUPS)]
    if spares:
        grouped.append(("other", "#f2f2f2", spares))
    order = [index for _, _, members in grouped for index in members]

    fig = Figure(figsize=(14, 9), constrained_layout=True)
    ax_time, ax_now = fig.subplots(2, 1, gridspec_kw={"height_ratios": [1.4, 1]})
    now = time.time()

    # 上の段: 時間変化。横軸は今からの時間 [分]
    for number, index in enumerate(order):
        if index not in histories:
            continue
        times, mean, _ = histories[index].snapshot()
        if len(times) == 0:
            continue
        style = "--" if marks.get(index) else "-"
        ax_time.plot((times - now) / 60, mean, style, marker=".", markersize=3,
                     color=colormaps["tab20"](number % 20), linewidth=1,
                     label=f"ha{index:02d}{marks.get(index, '')}")

    ax_time.set_xlabel("time from now [min]")
    ax_time.set_ylabel("MPPC temperature [\u00b0C] (mean of 6 sensors)")
    ax_time.grid(alpha=0.3)
    if ax_time.lines:
        ax_time.legend(fontsize=7, ncol=9, loc="upper left")
    ax_time.set_title("temperature history (10 s average, kept for 12 h, not reset by a new run)", fontsize=10)

    # 下の段: 今の温度 (直近 RATE_WINDOW 秒) と 6 センサの最小-最大
    positions = {index: number for number, index in enumerate(order)}
    for label, color, members in grouped:
        first, last = positions[members[0]] - 0.5, positions[members[-1]] + 0.5
        ax_now.axvspan(first, last, color=color, zorder=0)
        ax_now.text((first + last) / 2, 1.01, label, transform=ax_now.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

    for index, history in histories.items():
        times, mean, sensors = history.snapshot()
        recent = times >= now - RATE_WINDOW
        if not recent.any():
            recent = slice(-1, None)
        value = float(mean[recent].mean())
        sensor_now = sensors[recent].mean(axis=0)
        x = positions[index]
        ax_now.vlines(x, sensor_now.min(), sensor_now.max(), color="0.4", linewidth=3, zorder=2)
        ax_now.plot(x, value, "o", color="C3" if not marks.get(index) else "0.5", zorder=3)
        ax_now.text(x, sensor_now.max(), f"{value:.1f}{marks.get(index, '')}", ha="center", va="bottom",
                    fontsize=8, zorder=4)

    ax_now.set_xticks(range(len(order)), [f"ha{index:02d}" for index in order], fontsize=8, rotation=90)
    ax_now.set_xlim(-0.5, len(order) - 0.5)
    ax_now.set_ylabel("MPPC temperature [\u00b0C]")
    ax_now.grid(alpha=0.3, axis="y")
    ax_now.set_title(f"now (last {RATE_WINDOW:.0f} s): point = mean of 6 sensors, bar = min-max of the sensors. "
                     "* calibration copied from another HA, ? not in psa_param", fontsize=9, pad=16)

    fig.suptitle(title, fontsize=10)

    save_figure(fig, out_path)


class QuickLook:

    def __init__(self, args):

        self.args = args
        self.paths = {name: os.path.join(args.out_dir, f"{args.name}_{name}.png")
                      for name in ("overview", "pedestal", "temperature")}
        self.bad_packets = 0
        self.calibrations = load_temp_calibrations(args.psa_param)
        self.marks = {index: mark for index, (_, _, mark) in self.calibrations.items()}
        # 温度の時間変化は run が替わっても消さない
        self.temp_histories = {}
        self.reset()

    def reset(self):

        self.states = {}
        self.latest = {}
        self.n_events = 0

    def add(self, _side, index, frame):

        temp_a, temp_b, _ = self.calibrations.get(index, (DEFAULT_TEMP_A, DEFAULT_TEMP_B, "?"))

        try:
            trig, clock, temp, hits = decode_frame(frame, temp_a, temp_b)
        except Exception:
            self.bad_packets += 1
            return

        # どれかの HA のトリガ番号が大きく戻ったら新しい run とする
        if trig + RUN_RESET_TRIGGERS < self.latest.get(index, -1):
            print(f"[RUN] ha{index:02d} trig went back ({self.latest[index]} -> {trig}). "
                  "restart the accumulation", flush=True)
            self.reset()

        self.latest[index] = max(self.latest.get(index, -1), trig)
        self.states.setdefault(index, HaState()).add(trig, clock, temp, hits)
        self.temp_histories.setdefault(index, TempHistory()).add(time.time(), temp)
        self.n_events += 1

    def count_bad(self):

        self.bad_packets += 1

    def snapshot(self):

        shown = sorted(self.states) if self.args.ha is None else [i for i in self.args.ha if i in self.states]
        states = {index: self.states[index].snapshot() for index in shown}

        now = time.strftime("%H:%M:%S")
        title = (f"PSA QL {now}: {len(self.states)} HA, {self.n_events} thinned events, "
                 f"broken packets {self.bad_packets}")

        return states, title

    async def update_loop(self):
        """図を 1 枚ずつ順に描く。matplotlib を同時に使わないようにするため。"""

        while True:

            await asyncio.sleep(self.args.update)

            if len(self.states) == 0:
                continue

            # 描画に失敗しても受信は止めない
            try:
                states, title = self.snapshot()
                await asyncio.to_thread(draw_overview_figure, states, title, self.paths["overview"], self.marks)
                await asyncio.to_thread(draw_pedestal_figure, states, title, self.paths["pedestal"])
                await asyncio.to_thread(draw_temperature_figure, dict(self.temp_histories), self.marks, title,
                                        self.paths["temperature"])
                print(f"[DRAW] {title}", flush=True)
            except Exception:
                traceback.print_exc()


async def run(args):

    quick_look = QuickLook(args)

    # start_viewer はページの EVENT_MS も埋める。このページには 1 イベントの図が無いので、同じ間隔にする
    args.event_update = args.update
    prepare_view(args, quick_look.paths, PAGE, {"/temp": TEMP_PAGE})

    await asyncio.gather(receive_loop(args, "h", quick_look.add, quick_look.count_bad),
                         quick_look.update_loop())


def main():

    args = parse_args()

    try:
        asyncio.run(run(args))
    except KeyboardInterrupt:
        print("[STOP]")


if __name__ == "__main__":

    main()
