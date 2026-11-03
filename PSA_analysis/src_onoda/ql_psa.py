#!/usr/bin/env python3
#
# zeno から届く間引きイベントのうち HA の枠を受け、check_data.py と同じ確認図を積算しながら更新する。
#
# 受ける形式は ql_tpc.py と同じで、HA の枠は側 'h'、番号は HA 番号 (ha03 なら 3)。
# TPC の枠 (側 'a', 'c') は読み飛ばすので、ql_tpc.py と同じ送信元につないでよい。
# zeno が無いときは ql_replay.py が同じ形式で送る。起動するとブラウザが開き、図は自動で更新される。
#
# ブラウザには次の 4 枚を、モニターに収まるように並べる。
#   1 イベントの図 : 最新イベントの HG0-3 と LG0-3 の波形 (ペデスタルを引いた値)
#   レートの図     : HA ごとのトリガレート、ch ごとのペデスタル
#   一覧の図       : 全 HA を、底面 (ha00-05) を中央に、側面 (N: ha06, 07 / E: ha15, 16 /
#                    S: ha12, 13 / W: ha09, 10) をその上下左右に置いた十字形に並べる。
#                    HA ごとに 6 ch のバタフライマップを小さく描き、ヒットの無い ch は赤く塗る
#   詳細の図       : ページのボタン (または ← →) で選んだ 1 台の、6 ch のバタフライマップ、
#                    hg と peak_id の分布、最新イベントの HG / LG の波形、ペデスタル
#
# 解析は psa_decode.py の decode_events をそのまま使う。
#

import argparse
import asyncio
import json
import math
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

HERE = Path(__file__).resolve().parent

# psa_decode.py は analysis_sm3_FM. から始まる名前で sub_scr を読むので、その親も入れる
sys.path[:0] = [str(HERE), str(HERE.parents[2]), str(HERE.parents[1] / "TPC_analysis")]

from check_data import (  # noqa: E402
    HG_BINS,
    HG_MAX,
    HG_MIN,
    N_BINS_X,
    N_BINS_Y,
    N_CHANNELS,
    N_SAMPLES,
    X_MAX,
    X_MIN,
    Y_MAX,
    Y_MIN,
)
from psa_decode import DEFAULT_TEMP_A, DEFAULT_TEMP_B, decode_events  # noqa: E402
from ql_tpc import (  # noqa: E402
    RUN_RESET_TRIGGERS,
    prepare_view,
    receive_loop,
    save_figure,
)


DEFAULT_OUTPUT_DIR = "/home/onoda/analysis_sm3_FM/check_output"

# HA の clk は 1 MHz の 32 bit カウンタ
HA_CLOCK = 1e6
HA_CLOCK_WRAP = 2 ** 32

# 間引くとイベントが少ないので、バタフライマップのビンは check_data.py の半分にする
BUTTERFLY_BINS = (N_BINS_X // 2, N_BINS_Y // 2)

SIGNALS = ["HG0", "HG1", "HG2", "HG3", "LG0", "LG1", "LG2", "LG3"]

# 詳細の図を用意する HA 番号の範囲 (HA_0 から HA_17)
N_HA = 18

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

# 一覧の図は 1 枚が小さいので、ビンをさらに粗くする。細かいままだと 1 ビンが 1 画素より小さくなる
OVERVIEW_BINS = (N_BINS_X // 5, N_BINS_Y // 5)

# 上段に 1 イベントの図とレートの図、HA を選ぶボタン、下段の左に一覧、右に選んだ HA の詳細を並べる
PAGE = """<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<title>PSA QL</title>
<style>
html, body { margin: 0; height: 100%; overflow: hidden; background: #fff; font: 13px sans-serif; }
body { display: grid; height: 100vh; grid-template-columns: 56fr 44fr;
       grid-template-rows: 25vh 3.6vh 71.4vh; }
img { width: 100%; height: 100%; min-height: 0; object-fit: contain; }
#bar { grid-column: 1 / 3; display: flex; align-items: center; gap: 10px; padding: 0 6px; }
#bar .group { display: flex; align-items: center; gap: 2px; padding: 1px 4px; border-radius: 4px; }
#bar button { padding: 1px 5px; font: 12px sans-serif; }
#bar button.selected { background: #1f77b4; color: #fff; }
#bar button:disabled { opacity: 0.35; }
#status { margin-left: auto; color: #555; }
</style>
</head>
<body>
<img id="event" alt="">
<img id="stats" alt="">
<div id="bar"><span id="buttons" style="display: contents"></span><span id="status">図を待っています</span></div>
<img id="overview" alt="">
<img id="detail" alt="">
<script>
const statusLine = document.getElementById("status");
const buttons = document.getElementById("buttons");
const GROUPS = HA_GROUPS_JSON;
const N_HA = N_HA_COUNT;
const order = [];
let selected = null;
let available = [];

function name(index) { return "ha" + String(index).padStart(2, "0"); }

function load(id, url) {
  const next = new Image();
  next.onload = () => {
    document.getElementById(id).src = next.src;
    statusLine.textContent = "読み込み " + new Date().toLocaleTimeString();
  };
  next.src = url + "?t=" + Date.now();
}

function select(index) {
  selected = index;
  for (const button of buttons.querySelectorAll("button")) {
    button.classList.toggle("selected", Number(button.dataset.index) === index);
  }
  refreshDetail();
}

function refreshDetail() {
  if (selected !== null) load("detail", "/" + name(selected) + ".png");
}

// 詳細の図がある HA だけボタンを押せるようにする。最初は図のある最初の HA を選ぶ
function refreshButtons() {
  const found = [];
  let pending = N_HA;
  for (let index = 0; index < N_HA; index++) {
    const probe = new Image();
    probe.onload = () => { found.push(index); done(); };
    probe.onerror = done;
    probe.src = "/" + name(index) + ".png?t=" + Date.now();
  }
  function done() {
    if (--pending > 0) return;
    available = order.filter((index) => found.includes(index));
    for (const button of buttons.querySelectorAll("button")) {
      button.disabled = !available.includes(Number(button.dataset.index));
    }
    if ((selected === null || !available.includes(selected)) && available.length > 0) {
      select(available[0]);
    }
  }
}

// 面ごとにボタンをまとめる。面に入っていない HA は「other」にまとめる
const grouped = GROUPS.flatMap((group) => group[2]);
const others = [...Array(N_HA).keys()].filter((index) => !grouped.includes(index));
for (const [label, color, members] of GROUPS.concat([["other", "#f2f2f2", others]])) {
  const box = document.createElement("span");
  box.className = "group";
  box.style.background = color;
  box.append(label + ":");
  for (const index of members) {
    const button = document.createElement("button");
    button.textContent = name(index);
    button.dataset.index = index;
    button.disabled = true;
    button.onclick = () => select(index);
    box.appendChild(button);
    order.push(index);
  }
  buttons.appendChild(box);
}

// ← → で、図のある HA をボタンの順に選ぶ
document.addEventListener("keydown", (event) => {
  if (available.length === 0 || !["ArrowLeft", "ArrowRight"].includes(event.key)) return;
  const step = event.key === "ArrowRight" ? 1 : -1;
  const position = Math.max(available.indexOf(selected), 0);
  select(available[(position + step + available.length) % available.length]);
});

function every(interval, refresh) { refresh(); setInterval(refresh, interval); }
every(EVENT_MS, () => load("event", "/event.png"));
every(SUMMARY_MS, () => load("stats", "/stats.png"));
every(SUMMARY_MS, () => load("overview", "/overview.png"));
every(SUMMARY_MS, () => { refreshButtons(); refreshDetail(); });
</script>
</body>
</html>
""".replace("N_HA_COUNT", str(N_HA)).replace("HA_GROUPS_JSON", json.dumps(HA_GROUPS))


def parse_args():

    parser = argparse.ArgumentParser(description="PSA (HA) quick look for thinned events from zeno")
    parser.add_argument("--host", default="zeno", help="host sending the thinned events")
    parser.add_argument("--port", type=int, default=10012, help="TCP port of the thinned events")
    parser.add_argument("--ha", type=int, nargs="+", default=None,
                        help="HA numbers to show in the accumulated figure (default: all received)")
    parser.add_argument("--update", type=float, default=10.0,
                        help="update interval of the accumulated figures [s]")
    parser.add_argument("--event-update", type=float, default=1.0,
                        help="update interval of the one-event figure [s]")
    parser.add_argument("--rate-bin", type=float, default=60.0, help="bin width of the trigger rate [s]")
    parser.add_argument("--reconnect", type=float, default=3.0, help="reconnect interval [s]")
    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--name", default="ql_psa", help="name of the output PNG")
    parser.add_argument("--http-bind", default="127.0.0.1", help="address of the viewer")
    parser.add_argument("--http-port", type=int, default=8091, help="TCP port of the viewer")
    parser.add_argument("--no-open", action="store_true", help="do not open a browser")

    return parser.parse_args()


def decode_frame(frame):
    """1 枠を psa_decode の decode_events で読み、ヒットした ch ごとの値を返す。"""

    buf = np.frombuffer(frame, dtype=np.uint8).reshape(1, -1)
    sizes = np.array([len(frame)], dtype=np.int32)

    (trig, clock, hit, _, temp, ped, _, hg, _, x_hg, y_hg, _, _,
     peak_id, _, satu_h, _, adc) = decode_events(buf, sizes, False, DEFAULT_TEMP_A, DEFAULT_TEMP_B, True)

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
        self.latest_event = None    # ヒットのあった最新イベント (trig, ヒット)。詳細の図の波形に使う
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

        if len(hits) > 0:
            self.latest_event = (trig, hits)

    def snapshot(self):

        state = {name: np.array(values) for name, values in self.hits.items()}
        state["ped"] = np.array(self.ped, dtype=np.float32).reshape(-1, len(SIGNALS))
        state["trig"] = np.array(self.trig, dtype=np.int64)
        state["seconds"] = np.array(self.seconds, dtype=np.float64)
        state["latest_event"] = self.latest_event

        return state


def trigger_rate(state, rate_bin):
    """区間ごとの トリガ番号の増え方 / 時間 [Hz] を返す。"""

    seconds, trig = state["seconds"], state["trig"]

    if len(trig) < 2:
        return None

    n_bins = int(seconds.max() // rate_bin) + 1
    edges = np.arange(n_bins + 1) * rate_bin
    rate = np.full(n_bins, np.nan)

    for index in range(n_bins):

        selected = (seconds >= edges[index]) & (seconds < edges[index + 1])

        if selected.sum() < 2:
            continue

        span = seconds[selected].max() - seconds[selected].min()

        if span > 0:
            rate[index] = (trig[selected].max() - trig[selected].min()) / span

    return edges, rate


def draw_waveforms(ax_hg, ax_lg, hits, note="pedestal subtracted", legend_size=7):
    """1 イベントの HG0-3 と LG0-3 の波形 (ペデスタルを引いた値) を描き、ヒットごとの説明を返す。"""

    samples = np.arange(N_SAMPLES)
    labels = []

    for hit in hits:

        wave = hit["wave"] - hit["ped"]
        style = "-" if len(labels) == 0 else "--"

        for signal in range(4):
            ax_hg.plot(samples, wave[:, signal], style, color=f"C{signal}", label=f"ch{hit['ch']} {SIGNALS[signal]}")
            ax_lg.plot(samples, wave[:, 4 + signal], style, color=f"C{signal}",
                       label=f"ch{hit['ch']} {SIGNALS[4 + signal]}")

        labels.append(f"ch{hit['ch']}: hg {hit['hg']:.3g}, x_hg {hit['x_hg']:+.3f}, "
                      f"y_hg {hit['y_hg']:+.3f}, peak_id {hit['peak_id']}"
                      + (", saturated" if hit["satu_h"] else ""))

    for ax, name in ((ax_hg, "HG"), (ax_lg, "LG")):
        ax.set_title(f"{name} waveform ({note})")
        ax.set_xlabel("sample")
        ax.set_ylabel("ADC")
        ax.set_xlim(0, N_SAMPLES - 1)
        if len(hits) > 0:
            ax.legend(fontsize=legend_size, ncol=2)

    return labels


def draw_event_figure(event, out_path):
    """最新イベントの HG0-3 と LG0-3 の波形を描く。"""

    fig = Figure(figsize=(12, 3.6), constrained_layout=True)
    ax_hg, ax_lg = fig.subplots(1, 2)

    index, trig, hits = event
    labels = draw_waveforms(ax_hg, ax_lg, hits)

    fig.suptitle(f"ha{index:02d}   trig {trig}   " + ("   ".join(labels) if labels else "no hit"),
                 fontsize=10)

    save_figure(fig, out_path)


def draw_stats_figure(states, title, out_path, args):
    """HA ごとのトリガレートと、ch ごとのペデスタルを描く。"""

    fig = Figure(figsize=(12, 3.6), constrained_layout=True)
    ax_rate, ax_ped = fig.subplots(1, 2)

    for index, state in states.items():

        curve = trigger_rate(state, args.rate_bin)

        if curve is not None:
            edges, rate = curve
            ax_rate.stairs(rate, edges, label=f"ha{index:02d}")

    ax_rate.set_title(f"trigger rate ({args.rate_bin:g} s bins)")
    ax_rate.set_xlabel("time from the first event [s]")
    ax_rate.set_ylabel("trigger / s")
    if len(states) > 0:
        ax_rate.legend(fontsize=6, ncol=math.ceil(len(states) / 6))

    names = []
    medians = []

    for index, state in states.items():
        for ch in range(N_CHANNELS):
            selected = state["ch"] == ch if len(state["ch"]) else np.zeros(0, dtype=bool)
            names.append(f"ha{index:02d}")
            medians.append(np.median(state["ped"][selected], axis=0) if selected.any()
                           else np.full(len(SIGNALS), np.nan))

    if len(medians) > 0:
        medians = np.array(medians)
        positions = np.arange(len(names))
        for signal, name in enumerate(SIGNALS):
            ax_ped.plot(positions, medians[:, signal], marker="o", markersize=3,
                        linestyle="-" if signal < 4 else "--", label=name)
        # 目盛りは HA ごとに 1 つ (ch0 の位置) にし、HA の境目に線を引く
        ax_ped.set_xticks(positions[::N_CHANNELS], names[::N_CHANNELS], fontsize=7, rotation=90)
        ax_ped.set_xticks(positions, minor=True)
        for edge in positions[::N_CHANNELS][1:]:
            ax_ped.axvline(edge - 0.5, color="0.85", linewidth=0.8)
        ax_ped.legend(fontsize=6, ncol=4)

    ax_ped.set_title("pedestal median per ch")
    ax_ped.set_ylabel("pedestal [ADC]")

    fig.suptitle(title, fontsize=10)

    save_figure(fig, out_path)


def butterfly_counts(state, ch, bins):
    """1 ch のバタフライマップのヒストグラムと、その ch のヒット数を返す。"""

    selected = state["ch"] == ch if len(state["ch"]) else np.zeros(0, dtype=bool)
    n_hits = int(selected.sum())
    counts, x_edges, y_edges = np.histogram2d(
        state["x_hg"][selected] if n_hits else [], state["y_hg"][selected] if n_hits else [],
        bins=bins, range=[[X_MIN, X_MAX], [Y_MIN, Y_MAX]])

    return counts, x_edges, y_edges, n_hits


def draw_ha_tile(tile, index, state):
    """一覧の図の 1 台分。6 ch のバタフライマップを 2 行 3 列で小さく描き、ヒットの無い ch は赤く塗る。"""

    n_hits = 0 if state is None else len(state["ch"])
    tile.suptitle(f"ha{index:02d}  {n_hits} hits" if state is not None else f"ha{index:02d}  no data",
                  fontsize=9, fontweight="bold")
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
        elif n_hits > 0:
            # ほかの ch にはヒットがあるのに、この ch だけ無い
            ax.set_facecolor("#f4b6b6")

        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.set_title(f"ch{ch}  {n_ch}", fontsize=7, pad=1)


def draw_overview_figure(states, title, out_path):
    """全 HA を、底面を中央に、側面をその上下左右に置いた十字形に並べる。

    配置にある HA は、データがまだ無くても灰色の枠を出す。配置に無い HA は四隅に置く。
    """

    fig = Figure(figsize=(14, 10.5), constrained_layout=True)
    grid = fig.add_gridspec(4, 5)

    for label, color, members in HA_GROUPS:

        place, (n_rows, n_cols) = GROUP_PLACES[label]
        group = fig.add_subfigure(grid[place], facecolor=color, edgecolor="0.6", linewidth=1)
        group.suptitle(label, fontsize=10, fontweight="bold", x=0.01, ha="left")
        tiles = np.atleast_1d(group.subfigures(n_rows, n_cols, wspace=0.03, hspace=0.03)).ravel()

        for tile, index in zip(tiles, members):
            draw_ha_tile(tile, index, states.get(index))

    grouped = {index for _, _, members in HA_GROUPS for index in members}
    spares = [index for index in states if index not in grouped]

    for place, index in zip(SPARE_PLACES, spares):
        tile = fig.add_subfigure(grid[place], facecolor="#f2f2f2", edgecolor="0.6", linewidth=1)
        draw_ha_tile(tile, index, states[index])

    fig.suptitle(title, fontsize=11)

    save_figure(fig, out_path)


def draw_ha_pedestal(ax, state, signals, name):
    """1 台分の、signals (HG0-3 か LG0-3) のペデスタル。ch ごとに中央値と標準偏差を描く。"""

    positions = np.arange(len(signals))

    for ch in range(N_CHANNELS):

        selected = state["ch"] == ch if len(state["ch"]) else np.zeros(0, dtype=bool)

        if not selected.any():
            continue

        ped = state["ped"][selected][:, list(signals)]
        ax.errorbar(positions + (ch - 2.5) * 0.12, np.median(ped, axis=0), yerr=ped.std(axis=0),
                    fmt="o", markersize=3, capsize=2, label=f"ch{ch}")

    ax.set_xticks(positions, [SIGNALS[signal] for signal in signals], fontsize=7)
    ax.set_xlim(-0.5, len(signals) - 0.5)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_title(f"{name} pedestal", fontsize=10)
    if name == "HG":
        ax.set_ylabel("median, std [ADC]")
        if len(state["ch"]) > 0:
            ax.legend(fontsize=6, ncol=2)


def draw_detail_figure(index, state, out_path):
    """1 台分の詳細。6 ch のバタフライマップ (2 行 3 列)、hg と peak_id の分布、最新イベントの
    HG / LG の波形、ペデスタルを描く。"""

    fig = Figure(figsize=(11, 13.5), constrained_layout=True)
    grid = fig.add_gridspec(4, 3)
    ax_hg = fig.add_subplot(grid[2, 0:2])
    ax_peak = fig.add_subplot(grid[2, 2])
    ax_wave_hg = fig.add_subplot(grid[3, 0])
    ax_wave_lg = fig.add_subplot(grid[3, 1])
    # HG と LG はペデスタルの高さが大きく違うので、縦軸を分ける
    ax_ped_hg, ax_ped_lg = grid[3, 2].subgridspec(1, 2, wspace=0.05).subplots()

    for ch in range(N_CHANNELS):

        ax = fig.add_subplot(grid[ch // 3, ch % 3])
        counts, x_edges, y_edges, n_hits = butterfly_counts(state, ch, BUTTERFLY_BINS)

        if counts.max() > 0:
            ax.pcolormesh(x_edges, y_edges, counts.T, norm=LogNorm(vmin=1), cmap="viridis")

        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.set_title(f"ch{ch} ({n_hits} hits)", fontsize=10)
        ax.set_xlabel("x_hg")
        ax.set_ylabel("y_hg")

        if n_hits > 0:
            selected = state["ch"] == ch
            ax_hg.hist(state["hg"][selected], bins=HG_BINS // 4, range=(HG_MIN, HG_MAX),
                       histtype="step", label=f"ch{ch}")
            ax_peak.hist(state["peak_id"][selected], bins=N_SAMPLES, range=(-0.5, N_SAMPLES - 0.5),
                         histtype="step", label=f"ch{ch}")

    ax_hg.set_yscale("log")
    ax_hg.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    ax_hg.set_title("hg", fontsize=10)
    ax_hg.set_xlabel("hg [ADC x sample]")
    ax_peak.set_title("peak_id", fontsize=10)
    ax_peak.set_xlabel("peak_id [sample]")

    for ax in (ax_hg, ax_peak):
        if len(state["ch"]) > 0:
            ax.legend(fontsize=7)

    # 最新イベントの波形。ヒットの説明は図が混むので、どの trig かだけを題に足す
    if state["latest_event"] is not None:
        trig, hits = state["latest_event"]
        draw_waveforms(ax_wave_hg, ax_wave_lg, hits, note=f"trig {trig}", legend_size=6)
    else:
        for ax, name in ((ax_wave_hg, "HG"), (ax_wave_lg, "LG")):
            ax.set_title(f"{name} waveform (no hit yet)", fontsize=10)

    draw_ha_pedestal(ax_ped_hg, state, range(4), "HG")
    draw_ha_pedestal(ax_ped_lg, state, range(4, 8), "LG")

    group = next((label for label, _, members in HA_GROUPS if index in members), "other")
    fig.suptitle(f"ha{index:02d} ({group})  {len(state['ch'])} hits", fontsize=12, fontweight="bold")

    save_figure(fig, out_path)


class QuickLook:

    def __init__(self, args):

        self.args = args
        names = ["event", "stats", "overview"] + [f"ha{index:02d}" for index in range(N_HA)]
        self.paths = {name: os.path.join(args.out_dir, f"{args.name}_{name}.png") for name in names}
        self.bad_packets = 0
        self.reset()

    def reset(self):

        self.states = {}
        self.latest = {}
        self.latest_event = None
        self.n_events = 0

    def add(self, _side, index, frame):

        try:
            trig, clock, temp, hits = decode_frame(frame)
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
        self.n_events += 1

        if len(hits) > 0:
            self.latest_event = (index, trig, hits)

    def count_bad(self):

        self.bad_packets += 1

    def snapshot(self):

        shown = sorted(self.states) if self.args.ha is None else [i for i in self.args.ha if i in self.states]
        states = {index: self.states[index].snapshot() for index in shown}

        now = time.strftime("%H:%M:%S")
        title = (f"PSA QL {now}: {len(self.states)} HA, {self.n_events} thinned events, "
                 f"broken packets {self.bad_packets}")

        return states, title

    def draw_details(self, states):
        """HA ごとに詳細の図を描く。今の積算に無い HA の図は、前の run のものなので消す。"""

        for index in range(N_HA):
            path = self.paths[f"ha{index:02d}"]
            if index in states:
                draw_detail_figure(index, states[index], path)
            elif os.path.exists(path):
                os.remove(path)

    async def update_loop(self):
        """図を 1 枚ずつ順に描く。matplotlib を同時に使わないようにするため。"""

        drawn_event = None
        next_accumulated = time.monotonic() + self.args.update

        while True:

            await asyncio.sleep(self.args.event_update)

            # 描画に失敗しても受信は止めない
            try:
                event = self.latest_event
                if event is not None and event is not drawn_event:
                    drawn_event = event
                    await asyncio.to_thread(draw_event_figure, event, self.paths["event"])

                if len(self.states) > 0 and time.monotonic() >= next_accumulated:
                    next_accumulated = time.monotonic() + self.args.update
                    states, title = self.snapshot()
                    await asyncio.to_thread(draw_stats_figure, states, title, self.paths["stats"], self.args)
                    await asyncio.to_thread(draw_overview_figure, states, title, self.paths["overview"])
                    await asyncio.to_thread(self.draw_details, states)
                    print(f"[DRAW] {title}", flush=True)
            except Exception:
                traceback.print_exc()


async def run(args):

    quick_look = QuickLook(args)

    prepare_view(args, quick_look.paths, PAGE)

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
