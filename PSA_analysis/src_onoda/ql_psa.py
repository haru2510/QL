#!/usr/bin/env python3
#
# zeno から届く間引きイベントのうち HA の枠を受け、check_data.py と同じ確認図を積算しながら更新する。
#
# 受ける形式は ql_tpc.py と同じで、HA の枠は側 'h'、番号は HA 番号 (ha03 なら 3)。
# TPC の枠 (側 'a', 'c') は読み飛ばすので、ql_tpc.py と同じ送信元につないでよい。
# zeno が無いときは ql_replay.py が同じ形式で送る。起動するとブラウザが開き、図は自動で更新される。
#
# ブラウザには次の 3 枚を、モニターに収まるように並べる。
#   1 イベントの図 : 最新イベントの HG0-3 と LG0-3 の波形 (ペデスタルを引いた値)
#   レートの図     : HA ごとのトリガレート、ch ごとのペデスタル
#   積算の図       : HA ごとに 1 行。6 ch のバタフライマップ、hg の分布、peak_id の分布
#
# 解析は psa_decode.py の decode_events をそのまま使う。
#

import argparse
import asyncio
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
    VIEW_PAGE,
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

PAGE = VIEW_PAGE.replace("TPC QL", "PSA QL")


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


def draw_event_figure(event, out_path):
    """最新イベントの HG0-3 と LG0-3 の波形を描く。"""

    fig = Figure(figsize=(12, 3.6), constrained_layout=True)
    ax_hg, ax_lg = fig.subplots(1, 2)

    index, trig, hits = event
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
        ax.set_title(f"{name} waveform (pedestal subtracted)")
        ax.set_xlabel("sample")
        ax.set_ylabel("ADC")
        ax.set_xlim(0, N_SAMPLES - 1)
        if len(hits) > 0:
            ax.legend(fontsize=7, ncol=2)

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
        ax_rate.legend(fontsize=7)

    names = []
    medians = []

    for index, state in states.items():
        for ch in range(N_CHANNELS):
            selected = state["ch"] == ch if len(state["ch"]) else np.zeros(0, dtype=bool)
            names.append(f"ha{index:02d}\nch{ch}")
            medians.append(np.median(state["ped"][selected], axis=0) if selected.any()
                           else np.full(len(SIGNALS), np.nan))

    if len(medians) > 0:
        medians = np.array(medians)
        positions = np.arange(len(names))
        for signal, name in enumerate(SIGNALS):
            ax_ped.plot(positions, medians[:, signal], marker="o", markersize=3,
                        linestyle="-" if signal < 4 else "--", label=name)
        ax_ped.set_xticks(positions, names, fontsize=6)
        ax_ped.legend(fontsize=6, ncol=2)

    ax_ped.set_title("pedestal median per ch")
    ax_ped.set_ylabel("pedestal [ADC]")

    fig.suptitle(title, fontsize=10)

    save_figure(fig, out_path)


def draw_summary_figure(states, out_path):
    """HA ごとに 1 行、6 ch のバタフライマップ、hg の分布、peak_id の分布を描く。"""

    n_rows = max(len(states), 1)
    fig = Figure(figsize=(24, 3.2 * n_rows), constrained_layout=True)
    grid = fig.add_gridspec(n_rows, N_CHANNELS + 2)

    for row, (index, state) in enumerate(states.items()):

        ax_hg = fig.add_subplot(grid[row, N_CHANNELS])
        ax_peak = fig.add_subplot(grid[row, N_CHANNELS + 1])

        for ch in range(N_CHANNELS):

            selected = state["ch"] == ch if len(state["ch"]) else np.zeros(0, dtype=bool)
            n_hits = int(selected.sum())

            ax = fig.add_subplot(grid[row, ch])
            counts, x_edges, y_edges = np.histogram2d(
                state["x_hg"][selected] if n_hits else [], state["y_hg"][selected] if n_hits else [],
                bins=BUTTERFLY_BINS, range=[[X_MIN, X_MAX], [Y_MIN, Y_MAX]])

            if counts.max() > 0:
                ax.pcolormesh(x_edges, y_edges, counts.T, norm=LogNorm(vmin=1), cmap="viridis")

            ax.set_xlim(X_MIN, X_MAX)
            ax.set_ylim(Y_MIN, Y_MAX)
            ax.set_title(f"ha{index:02d} ch{ch} ({n_hits} hits)", fontsize=9)
            ax.set_xlabel("x_hg")
            ax.set_ylabel("y_hg")

            if n_hits > 0:
                ax_hg.hist(state["hg"][selected], bins=HG_BINS // 4, range=(HG_MIN, HG_MAX),
                           histtype="step", label=f"ch{ch}")
                ax_peak.hist(state["peak_id"][selected], bins=N_SAMPLES, range=(-0.5, N_SAMPLES - 0.5),
                             histtype="step", label=f"ch{ch}")

        ax_hg.set_yscale("log")
        ax_hg.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        ax_hg.set_title(f"ha{index:02d} hg", fontsize=9)
        ax_hg.set_xlabel("hg [ADC x sample]")
        ax_peak.set_title(f"ha{index:02d} peak_id", fontsize=9)
        ax_peak.set_xlabel("peak_id [sample]")

        for ax in (ax_hg, ax_peak):
            if len(state["ch"]) > 0:
                ax.legend(fontsize=6)

    save_figure(fig, out_path)


class QuickLook:

    def __init__(self, args):

        self.args = args
        self.paths = {name: os.path.join(args.out_dir, f"{args.name}_{name}.png")
                      for name in ("event", "stats", "summary")}
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
                    await asyncio.to_thread(draw_summary_figure, states, self.paths["summary"])
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
