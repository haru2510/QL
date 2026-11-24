#!/usr/bin/env python3
"""TPCの全ストリップ・クロックヒットをXZ/YZに積算する。ToT補正・XY選別はしない。"""
import argparse  # 設定と再描画モードをコマンドラインで選ぶ。
import json  # 条件・件数・入力ファイル情報を保存する。
import sys  # 既存のrawデコーダを読み込む検索先を追加する。
from pathlib import Path  # 入出力パスを扱う。
import numpy as np  # クロック別のヒストグラムとイベント番号を扱う。
from numba import njit  # rawのビット走査を機械語にして処理時間を抑える。
import matplotlib  # 保存専用の描画環境を設定する。
matplotlib.use("Agg")  # GUIを開かずPNGを出力する。
import matplotlib.pyplot as plt  # カット前後の投影図を描く。
from matplotlib.colors import LogNorm  # ヒット数の色を対数で表示する。

REPO = Path(__file__).resolve().parents[1]  # このスクリプトの親からリポジトリの位置を求める。
sys.path.insert(0, str(REPO / "ETCC_analysis"))  # 既存デコーダが置かれた場所を指定する。
import beam_vertex as decoder  # 検証済みのTPCフレーム読み取りを再使用する。
from settings import DEFAULT_CONFIG, load_config, repo_path, output_file  # 外部JSONと入出力パスを全コードで共通に扱う。


@njit  # 各ビットを読む内側のループだけを高速化する。
def accumulate(buf, offsets, n_records, keep, reject, axis, positions, clock_limit, all_hist, kept_hist):  # 1ファイルの投影と低Zヒット判定を行う。
    bad = np.zeros(len(offsets), np.bool_)  # 各フレームに低Zヒットがあるか記録する。
    hit_counts = np.zeros(len(offsets), np.int64)  # フレームごとの全rawヒット数を記録する。
    for frame in range(len(offsets)):  # ファイル内のフレームを順に読む。
        if not keep[frame]:  # 重複して記録されたトリガーは再計数しない。
            continue  # このフレームを飛ばす。
        first = offsets[frame] + decoder.HEAD_SIZE + decoder.FADC_SIZE + 4  # 最初のストリップレコードの位置を求める。
        for index in range(n_records[frame]):  # 保存された全レコードを使い、早い16クロックには限定しない。
            start = first + decoder.RECORD_SIZE * index  # 現在のレコードのバイト位置を求める。
            clock = decoder.read_word(buf, start) & 0x7ff  # 下位11ビットからクロックを得る。
            if clock >= decoder.N_CLOCKS:  # 有効範囲外のクロックを除く。
                continue  # 不正クロックをヒット数に含めない。
            for channel in range(decoder.N_FADC_CH):  # 32ストリップ分のワードを4個読む。
                word = decoder.read_word(buf, start + 4 * (4 - channel))  # raw仕様の逆順に格納されたビット列を読む。
                for bit in range(32):  # ワード内の各ストリップを調べる。
                    if ((word >> bit) & 1) == 0:  # ヒットフラグが立っていない場合を除く。
                        continue  # 非ヒットのストリップを飛ばす。
                    strip = channel * 32 + bit  # ボード内のストリップ番号0〜127に直す。
                    coordinate = positions[strip]  # 配置と向きを反映した座標ビン番号を得る。
                    all_hist[axis, coordinate, clock] += 1  # 全イベント側に1ヒットを加える。
                    hit_counts[frame] += 1  # 分母確認用にrawヒット総数も数える。
                    if not reject[frame]:  # イベント全体が低Zカットを通過しているか判定する。
                        kept_hist[axis, coordinate, clock] += 1  # 通過したイベントのヒットだけを加える。
                    if clock <= clock_limit:  # Z<=しきい値に対応する整数クロックを判定する。
                        bad[frame] = True  # アノード・カソードの片側でも除外対象にする。
    return bad, hit_counts  # 新規に求めた除外判定とヒット数を返す。


def scan_run(name, config, cached_veto):  # 1ランをファイル単位で走査してメモリ使用量を抑える。
    directories = [str(repo_path(config["data_root"]) / storage / name) for storage in config["storages"]]  # 両storageの入力先を作る。
    boards = decoder.collect_boards(directories, 0)  # 完成したrawファイルを全ボード分集める。
    if not boards:  # 入力が空の結果を正常終了として扱わない。
        raise FileNotFoundError(f"TPC raw files not found: {directories}")  # 見つからない入力先を明示する。
    shape = (2, 768, decoder.N_CLOCKS)  # X/Y、0.4 mmの768ストリップ、全クロックを保持する。
    all_hist = np.zeros(shape, np.int64)  # カット前のヒット数を確保する。
    kept_hist = np.zeros(shape, np.int64)  # カット後のヒット数を確保する。
    all_events = set()  # ボードをまたいだ全記録トリガーの和集合を作る。
    hit_events = set()  # rawヒットが1個以上あるトリガーも別に数える。
    fresh_veto = set()  # 今回読み直したrawから除外トリガーを独立に求める。
    cached_sorted = np.array(sorted(cached_veto or set()), dtype=np.int64)  # 保存済み除外集合を検索用の配列にする。
    manifests = []  # 読んだファイルのパス・サイズ・更新時刻を残す。
    frame_counts = {}  # 各ボードで読み取れたフレーム数を残す。
    clock_limit = int(np.floor((config["z_veto"] - config["z_gem"]) / config["velocity"] + config["clock_offset"]))  # mmの条件を整数クロックに換算する。
    for (side, board), paths in boards.items():  # アノード・カソードの各ボードを処理する。
        upic, local_board = divmod(board, decoder.N_BOARDS_PER_UPIC)  # uPIC番号とuPIC内ボード番号を求める。
        axis_name, sign, constant = decoder.STRIP_AXES[upic][0 if side == "a" else 1]  # 配置に従ってX/Yと正負の向きを選ぶ。
        axis = 0 if axis_name == "x" else 1  # ヒストグラムのX=0、Y=1に対応させる。
        strips = np.arange(decoder.N_STRIPS) + local_board * decoder.N_STRIPS  # uPIC内のストリップ番号を求める。
        positions_mm = config["pitch_strip"] * (sign * strips + constant)  # ストリップ中心の座標[mm]を計算する。
        positions = np.rint((positions_mm + 153.4) / 0.4).astype(np.int64)  # 共通の0.4 mm座標ビンに変換する。
        if np.any((positions < 0) | (positions >= 768)):  # 想定外の配置を配列の範囲外参照にしない。
            raise ValueError("TPC strip geometry is outside the configured 768 bins")  # 配置不一致を明示する。
        seen = set()  # 同じボードの重複フレームを1回だけ読む。
        epoch, previous, frames = 0, None, 0  # トリガーカウンタ巻き戻りの状態を初期化する。
        for path in paths:  # 大量のrawを1ファイルずつ処理する。
            stat = Path(path).stat()  # 入力のサイズと更新時刻を確認する。
            manifests.append({"path": path, "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})  # 将来の入力照合に必要な情報を保存する。
            buf, trig, offsets, records = decoder.scan_file(path)  # rawからフレームの位置とトリガー番号を得る。
            if len(trig) == 0:  # 空ファイルをヒットなしイベントと混同しない。
                print(f"[WARN] no frames: {path}", flush=True)  # 入力の問題をログに残す。
                continue  # 次のファイルへ進む。
            back = np.r_[previous is not None and trig[0] < previous, trig[1:] < trig[:-1]]  # カウンタが減った場所を調べる。
            epochs = epoch + np.cumsum(back)  # 巻き戻り後を別のイベント番号空間に分ける。
            epoch, previous = int(epochs[-1]), int(trig[-1])  # 次のファイルへ状態を引き継ぐ。
            events = epochs * (1 << 32) + trig  # (epoch,trig)を1個の整数にまとめる。
            _, first = np.unique(events, return_index=True)  # 同一ファイル内の重複も除く。
            keep = np.zeros(len(events), dtype=bool)  # 最初の出現かどうかを格納する。
            keep[first] = [int(events[i]) not in seen for i in first]  # 前のファイルにもないフレームだけを使う。
            seen.update(map(int, events[keep]))  # このボードで既に読んだトリガーを更新する。
            reject = np.isin(events, cached_sorted)  # 保存済みのイベント単位の除外判定を対応づける。
            bad, counts = accumulate(buf, offsets, records, keep, reject, axis, positions, clock_limit, all_hist, kept_hist)  # 全ヒットを積算し、除外集合も再確認する。
            all_events.update(map(int, events[keep]))  # 全ボードでの記録トリガーの和集合を更新する。
            hit_events.update(map(int, events[keep & (counts > 0)]))  # 少なくとも1ヒットあるイベントを数える。
            fresh_veto.update(map(int, events[keep & bad]))  # 低Zヒットが1個でもあるイベントを集める。
            frames += int(keep.sum())  # ボードの有効フレーム数を加える。
        frame_counts[f"{side}{board:02d}"] = frames  # ボードごとの読み取り件数を保存する。
        print(f"{name} {side}{board:02d}: {frames:,} frames", flush=True)  # 読み込みの進捗を知らせる。
    if cached_veto is not None and fresh_veto != cached_veto:  # 再接続やコピーで入力が変わっていないか確認する。
        raise RuntimeError(f"{name}: cached veto differs from current raw; rerun the veto scan before publishing")  # 古い選別集合で新しい図を公開しない。
    assert np.all(kept_hist <= all_hist)  # カット後のビン数がカット前を超えないことを確認する。
    negative_limit = int(np.ceil(config["clock_offset"] - config["z_gem"] / config["velocity"]))  # Z<0となるクロック範囲の上端を求める。
    metadata = {  # イベント数とヒット数を別々の項目にする。
        "run": name,  # ランをまたいだトリガー番号の混同を防ぐ。
        "all_events": len(all_events),  # 記録された全トリガーの和集合。
        "events_with_hits": len(hit_events),  # TPCが1ヒット以上あるイベント数。
        "veto_events": len(fresh_veto),  # 低Z側ヒットを含むイベント数。
        "passing_events": len(all_events - fresh_veto),  # ヒットなしも含む通過イベント数。
        "passing_events_with_hits": len(hit_events - fresh_veto),  # TPCヒットを持つ通過イベント数。
        "raw_strip_clock_hits": int(all_hist.sum()),  # 全rawストリップ・クロックヒット数。
        "passing_strip_clock_hits": int(kept_hist.sum()),  # 通過イベントでのrawヒット数。
        "negative_z_hits": int(all_hist[:, :, :negative_limit].sum()),  # 負のZに変換されるヒット数。
        "frames_per_board": frame_counts,  # ボードごとの読み取り確認用件数。
        "input_manifest": manifests,  # 入力ファイルのパス・サイズ・時刻。
    }  # 検証に使う情報をまとめる。
    return all_hist, kept_hist, metadata, fresh_veto  # 投影と再確認した除外イベントを返す。


def draw(histogram_path, out):  # 保存したヒストグラムだけから図を作り直す。
    data = np.load(histogram_path)  # SSDを使わず圧縮済み積算を読む。
    config = json.loads(str(data["config"]))  # 積算時の較正を読み戻す。
    x_edges = np.linspace(-153.6, 153.6, 97)  # 8ストリップをまとめた3.2 mmビンを作る。
    clock_edges = np.arange(-0.5, decoder.N_CLOCKS, 4)  # 4クロックごとのビン境界を作る。
    z_edges = config["velocity"] * (clock_edges - config["clock_offset"]) + config["z_gem"]  # ビン境界も同じ較正式でmmに直す。
    arrays = [data["all_hits"], data["kept_hits"]]  # カット前後の同じ定義のヒット数を取り出す。
    fig, axes = plt.subplots(2, 2, figsize=(13, 11), constrained_layout=True)  # XZ/YZをカット前後で並べる。
    for axis in range(2):  # XとYの直接観測された投影を別々に描く。
        images = [h[axis].reshape(96, 8, 256, 4).sum(axis=(1, 3)).T for h in arrays]  # 表示ビンへ単純加算し、イベントで重み付けしない。
        norm = LogNorm(vmin=1, vmax=max(int(images[0].max()), 2))  # カット前後で色のスケールを共通にする。
        for row in range(2):  # 同じ投影のカット前後を描く。
            ax = axes[row, axis]  # 対応するパネルを選ぶ。
            counts = np.ma.masked_less(images[row], 1)  # 空ビンを白く表示する。
            image = ax.pcolormesh(x_edges, z_edges, counts, cmap="turbo", norm=norm)  # rawストリップ・クロックヒットの積算を描く。
            fig.colorbar(image, ax=ax, label="raw strip-clock hits / bin")  # 色がイベント数ではないことを示す。
            ax.axhline(config["z_veto"], color="black", ls="--", lw=1)  # Z側カットの位置を描く。
            ax.axhline(0, color="gray", ls=":", lw=1)  # 負のZが分かるよう原点も描く。
            tag = "All recorded events" if row == 0 else f"After any-hit Z <= {config['z_veto']:g} mm event veto"  # 使用したイベント集合を説明する。
            ax.set(xlabel=("X" if axis == 0 else "Y") + " [mm]", ylabel="Z [mm] (down)", title=tag)  # 軸の単位を明示する。
            ax.invert_yaxis()  # Zが増える方向を下向きにする。
    fig.suptitle(", ".join(config["runs"]) + ": raw TPC track accumulation\nNo vertex selection, XY selection, ToT correction, or HA multiplicity requirement")  # 頂点図との違いを明記する。
    fig.savefig(out, dpi=150)  # カット前後の単純積算図を保存する。
    plt.close(fig)  # 描画メモリを解放する。


def main():  # コマンドからの実行をまとめる。
    parser = argparse.ArgumentParser(description=__doc__)  # コマンドの説明を作る。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 既定の条件を指定する。
    parser.add_argument("--rebuild-veto", action="store_true")  # 較正や入力が変わった場合は除外集合も作り直す。
    parser.add_argument("--plot-only", action="store_true")  # SSDなしで図だけを再作成できるようにする。
    args = parser.parse_args()  # 利用者の指定を読み取る。
    config = load_config(args.config)  # 共通の入力・較正・保存先を読む。
    out_dir = repo_path(config["result_dir"])  # 結果の保存先をリポジトリ基準で解決する。
    out_dir.mkdir(parents=True, exist_ok=True)  # 初回実行なら保存先を作る。
    histogram_path = output_file(config, "track_histogram")  # 再描画用の積算データを保存する名前を決める。
    histogram_path.parent.mkdir(parents=True, exist_ok=True)  # 積算専用の保存先を作る。
    if not args.plot_only:  # rawを読む必要がある場合だけ走査する。
        total_all = np.zeros((2, 768, decoder.N_CLOCKS), np.int64)  # 全ラン合計のカット前ヒストグラムを確保する。
        total_kept = np.zeros_like(total_all)  # 全ラン合計のカット後ヒストグラムを確保する。
        summaries = []  # 各ランの分母と入力情報を集める。
        for name in config["runs"]:  # ランを混ぜずに個別に積算する。
            cache_path = repo_path(config["outputs"]["veto_dir"]) / f"z0_track_veto_{name}.npz"  # ランごとのイベント除外リストを指定する。
            if args.rebuild_veto or not cache_path.exists():  # 初回か入力変更時はrawから除外集合を先に作る。
                _, _, _, fresh = scan_run(name, config, None)  # 最初の走査では全ボードの除外判定だけを確定する。
                cache_path.parent.mkdir(parents=True, exist_ok=True)  # 除外リストの保存先を作る。
                np.savez_compressed(cache_path, event=np.array(sorted(fresh), dtype=np.int64))  # 全ボードの判定を統合して保存する。
            cache = np.load(cache_path)  # 積算に使う確定した除外集合を読む。
            veto = set(map(int, cache["event"]))  # 全ボード共通のイベント除外集合に直す。
            all_hist, kept_hist, metadata, _ = scan_run(name, config, veto)  # rawを走査して除外集合の一致も確認する。
            total_all += all_hist  # カット前の投影をラン間で足す。
            total_kept += kept_hist  # カット後の投影をラン間で足す。
            summaries.append(metadata)  # イベント数の分母を保存する。
        np.savez_compressed(histogram_path, all_hits=total_all, kept_hits=total_kept, config=json.dumps(config))  # 全ビンを保存してSSDなしの再描画を可能にする。
        result = {"count_unit": "one raw strip hit at one clock; no ToT correction", "config": config, "runs": summaries}  # ヒット数の定義と条件を添える。
        output_file(config, "track_metadata").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")  # 入力と件数を読みやすいJSONへ保存する。
    draw(histogram_path, output_file(config, "track_figure"))  # 保存済みの配列から積算図を作る。
    print(f"Saved {output_file(config, 'track_figure')}", flush=True)  # 最終図の場所を報告する。


if __name__ == "__main__":  # importだけではSSDへアクセスしない。
    main()  # コマンドで呼ばれたときにだけ処理する。
