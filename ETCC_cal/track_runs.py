#!/usr/bin/env python3
"""run00〜07の全raw TPCヒットを、選別なしで積算する。カラー調整版と全レンジ版を保存する。"""
import argparse  # 外部JSONと再描画だけのモードを指定する。
import json  # 設定・入力一覧・カラースケールを保存する。
from pathlib import Path  # キャッシュと出力パスを扱う。
import numpy as np  # ヒット積算配列を計算する。
import matplotlib  # 画面を開かない描画環境を指定する。
matplotlib.use("Agg")  # 結果をPNGへ直接出力する。
import matplotlib.pyplot as plt  # XZとYZの積算を並べる。
from matplotlib.colors import LogNorm  # ヒット数を対数色で表示する。
from settings import load_config, repo_path, output_file  # 全解析共通のJSON入出力を使う。
from track_sum import scan_run, accumulate, decoder  # 同じrawフレーム読み取り・重複除去・ストリップ配置を使う。

DEFAULT_CONFIG = Path(__file__).parent / "config/20261124_run00_07.json"  # 全ラン積算用の外部設定を既定にする。


def read_temporary(name, config, completed):  # 完成rawを読み直さず、末尾tmpの有効な新規フレームだけを追加する。
    histogram = np.zeros((2, 768, 1024), dtype=np.int64)  # tmpから増えたヒットだけを保持する。
    files = []  # tmpの回収状況を記録する。
    for storage in config["storages"]:  # 両storageの末尾ファイルを調べる。
        directory = repo_path(config["data_root"]) / storage / name  # 対象ランのフォルダを求める。
        for path in sorted(directory.rglob("tpc_enc_*.raw.tmp")):  # 書きかけのファイルも全て確認する。
            if path.with_suffix("").exists():  # 同名の完成rawがあるときは完成側を優先する。
                continue  # 同じファイルを二重に数えない。
            match = decoder.FILE_RE.fullmatch(path.name[:-4])  # tmp拡張子だけを外し、ボード番号を読む。
            if match is None:  # TPC rawの形式に一致しないファイルを除く。
                continue  # 別形式のデータを誤読しない。
            side, board, sequence = match.groups()  # 読み出し側・ボード・連番を取り出す。
            board = int(board)  # 配置計算に使う整数へ変換する。
            previous = []  # 同じボードの完成rawのうち、このtmpより前のものを集める。
            for entry in completed["input_manifest"]:  # キャッシュを作った完成rawの一覧を使う。
                raw_match = decoder.FILE_RE.fullmatch(Path(entry["path"]).name)  # 完成rawの番号を読む。
                if raw_match and raw_match[1] == side and int(raw_match[2]) == board and int(raw_match[3]) < int(sequence):  # 同じボードの前方の連番だけを選ぶ。
                    previous.append((int(raw_match[3]), entry["path"]))  # 一番後ろの完成rawを選べるようにする。
            last_trigger = -1  # 完成rawがなければtmpの先頭から使う。
            if previous:  # 完成rawがある場合は境界のトリガーを確認する。
                _, triggers, _, _ = decoder.scan_file(max(previous)[1])  # 最後の完成ファイルだけを読み直す。
                if len(triggers):  # 空でない完成ファイルなら末尾を取得する。
                    last_trigger = int(triggers[-1])  # これ以前は既に完成raw側に含まれている。
            buf, trig, offsets, records = decoder.scan_file(str(path))  # tmp内から読み取れる有効フレームを回収する。
            if len(trig) and (np.any(trig[1:] < trig[:-1]) or int(trig[-1]) < last_trigger):  # 末尾でカウンタが巻き戻った場合を検出する。
                raise RuntimeError(f"Trigger reset in tmp: {path}; full event-aware rescan required")  # 別epochを誤って重複扱いしない。
            _, first = np.unique(trig, return_index=True)  # tmp内の重複トリガーも一度だけ数える。
            keep = np.zeros(len(trig), dtype=bool)  # 使用するフレームを示す配列を用意する。
            keep[first] = trig[first] > last_trigger  # 完成rawの末尾より新しいトリガーだけを追加する。
            upic, local_board = divmod(board, decoder.N_BOARDS_PER_UPIC)  # uPIC番号とボード内の配置を求める。
            axis_name, sign, constant = decoder.STRIP_AXES[upic][0 if side == "a" else 1]  # 既存積算と同じX/Yと向きを使う。
            strips = np.arange(decoder.N_STRIPS) + local_board * decoder.N_STRIPS  # uPIC内のストリップ番号を求める。
            positions_mm = config["pitch_strip"] * (sign * strips + constant)  # ストリップ中心をmmに直す。
            positions = np.rint((positions_mm + 153.4) / 0.4).astype(np.int64)  # 保存済み積算と同じビン番号へ変換する。
            reject = np.zeros(len(trig), dtype=bool)  # 全データ積算なので低Zイベントも除外しない。
            unused = np.zeros_like(histogram)  # 共通関数のカット後出力は使わない。
            _, counts = accumulate(buf, offsets, records, keep, reject, 0 if axis_name == "x" else 1, positions, -1, histogram, unused)  # 全有効クロックのヒットをそのまま足す。
            files.append({"path": str(path), "bytes": path.stat().st_size, "parsed_frames": len(trig), "added_frames": int(keep.sum()), "added_hits": int(counts.sum())})  # 読めた範囲と追加件数を明示する。
    return histogram, files  # tmpで増えた積算と入力一覧を返す。


def draw(hist, config, prefix, label):  # 保存した全ヒットからmm版・クロック版の図を作る。
    projections = [hist[i].reshape(96, 8, 256, 4).sum(axis=(1, 3)).T for i in range(2)]  # 8ストリップ・4クロックごとに単純加算する。
    positive = np.concatenate([p[p > 0] for p in projections])  # 空ビンを除いてカラー範囲を見積もる。
    if len(positive) == 0:  # ヒットが全くない入力は描かない。
        raise ValueError(f"No hits: {label}")  # 空の入力を明示する。
    percentiles = config["color_percentiles"]  # 利用者がJSONに指定した上下の百分位を使う。
    low, high = np.percentile(positive, percentiles)  # 大多数のヒットの濃淡が分かる範囲を求める。
    low = max(float(low), 1)  # 対数色の下限を正に保つ。
    high = max(float(high), low * 1.01)  # 一定値しかない場合でも上下限を分ける。
    x_edges = np.linspace(-153.6, 153.6, 97)  # 3.2 mm表示ビンの境界を使う。
    clocks = np.arange(-0.5, 1024, 4)  # 全有効クロック0〜1023を省略せず表示する。
    z_edges = config["velocity"] * (clocks - config["clock_offset"]) + config["z_gem"]  # 同じ全範囲を既存較正でmmにも換算する。
    scales = [("contrast", low, high), ("fullscale", 1.0, max(float(positive.max()), 2))]  # 調整版と下限1の全レンジ版を両方残す。
    for units, y_edges in (("z", z_edges), ("clock", clocks)):  # 較正に依存しないクロック表示も必ず作る。
        for mode, vmin, vmax in scales:  # カラースケールの違いだけを切り替える。
            fig, axes = plt.subplots(1, 2, figsize=(12, 7), constrained_layout=True)  # X投影とY投影を横に並べる。
            for i, ax in enumerate(axes):  # 同一スケールで両投影を描く。
                image = ax.pcolormesh(x_edges, y_edges, np.ma.masked_less(projections[i], 1), cmap="turbo", norm=LogNorm(vmin=vmin, vmax=vmax))  # 空ビンはマスクし、データ値は変更しない。
                fig.colorbar(image, ax=ax, label="raw strip-clock hits / bin", extend="both")  # 色範囲の外の値があることも示す。
                ax.set(xlabel=("X" if i == 0 else "Y") + " [mm]", ylabel="Z [mm] (down)" if units == "z" else "raw clock", title=("XZ" if i == 0 else "YZ") if units == "z" else ("X-clock" if i == 0 else "Y-clock"))  # 表示量と単位を明記する。
                ax.invert_yaxis()  # Zまたはクロックが増える方向を下向きにする。
            fig.suptitle(f"{label}: all recorded raw TPC hits\nNo event veto, vertex / XY / HA selection or ToT correction; {mode}: {vmin:.3g} .. {vmax:.3g}")  # 何も選別していないことと色範囲を示す。
            destination = prefix.with_name(prefix.name + f"_{units}_{mode}.png")  # mm/クロックとカラー方式をファイル名で区別する。
            fig.savefig(destination, dpi=150)  # SSDなしでも再利用できる図を保存する。
            plt.close(fig)  # 次の描画へメモリを残さない。
    return {"percentiles": percentiles, "contrast_vmin": low, "contrast_vmax": high, "fullscale_vmin": 1, "fullscale_vmax": int(positive.max())}  # 再現可能なカラー設定を返す。


def main():  # 個別ランと全ランの合計を作る。
    parser = argparse.ArgumentParser(description=__doc__)  # 入力設定の説明を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # raw入力・結果パス・カラーを外部JSONから指定する。
    parser.add_argument("--resume", action="store_true")  # 完了済みランのキャッシュを検証して読み直しを省く。
    parser.add_argument("--plot-only", action="store_true")  # 保存済み配列だけで再描画できるようにする。
    args = parser.parse_args()  # 利用者の指定を読む。
    config = load_config(args.config)  # 全ラン積算の設定を読み取る。
    cache_dir = repo_path(config["outputs"]["run_cache_dir"])  # 個別ランの積算データの保存先を求める。
    cache_dir.mkdir(parents=True, exist_ok=True)  # 初回実行ならフォルダを作る。
    total = np.zeros((2, 768, 1024), dtype=np.int64)  # 合算用の全座標・全クロック配列を用意する。
    metadata = []  # 個別ランの入力と件数を残す。
    for name in config["runs_to_scan"]:  # 今回追加するrun00〜05を個別に読む。
        cache = cache_dir / f"{name}.npz"  # ラン別の積算データ保存先を決める。
        details = cache_dir / f"{name}.json"  # 入力と件数の保存先を決める。
        reuse = args.resume and cache.exists() and details.exists()  # 配列と入力記録の両方が揃ったランだけを再使用する。
        if reuse:  # 前回の完了済みデータを使ってよいか確認する。
            saved = np.load(cache)  # 保存時の条件を取り出す。
            saved_config = json.loads(str(saved["config"]))  # 比較する較正と入力先を読む。
            for key in ("velocity", "clock_offset", "z_gem", "pitch_strip", "data_root", "storages"):  # 入力と座標変換が変わっていないか調べる。
                if saved_config[key] != config[key]:  # 同名でも条件が違うキャッシュを流用しない。
                    reuse = False  # このランをrawから作り直す。
            manifests = json.loads(details.read_text())["input_manifest"]  # 前回の入力ファイル情報を取り出す。
            for entry in manifests:  # ファイルのサイズと更新時刻も確認する。
                stat = Path(entry["path"]).stat()  # 記録された入力の現在の状態を読む。
                if stat.st_size != entry["bytes"] or stat.st_mtime_ns != entry["mtime_ns"]:  # コピーや更新が入っていないか調べる。
                    reuse = False  # 入力が変わっていたら再使用をやめる。
            saved.close()  # 条件確認用に開いた圧縮ファイルを閉じる。
        if not args.plot_only and not reuse:  # 完了済みキャッシュがない場合だけSSDを走査する。
            hist, _, information, veto = scan_run(name, config, None)  # 全ヒットを読む。除外集合は計算するが積算には適用しない。
            information.pop("passing_strip_clock_hits", None)  # 未選別の走査でカット後ヒット数を誤表示しない。
            information["selection"] = "None: all recorded raw TPC hits"  # 生の全ヒット集計であることを明記する。
            np.savez_compressed(cache, all_hits=hist, config=json.dumps(config))  # 正規化せず全ビンを保存する。
            details.write_text(json.dumps(information, indent=2) + "\n")  # raw入力の一覧と件数を保存する。
        data = np.load(cache)  # 新規または保存済みの全ヒット積算を読む。
        hist = data["all_hits"]  # カットなしの配列だけを使用する。

        information = json.loads(details.read_text())  # このランの件数を読み戻す。
        if config.get("include_tmp", False):  # 全データ版では末尾のtmpも含める。
            tmp_cache = cache_dir / f"{name}_tmp.npz"  # 完成rawのキャッシュとtmpの追加分を分ける。
            tmp_details = cache_dir / f"{name}_tmp.json"  # tmpの回収情報を別に保存する。
            if not args.plot_only:  # 再描画時はSSDを読み直さない。
                added, sources = read_temporary(name, config, information)  # 末尾tmpの新規ヒットを回収する。
                np.savez_compressed(tmp_cache, all_hits=added)  # tmpの追加分だけを保存する。
                tmp_details.write_text(json.dumps(sources, indent=2) + "\n")  # 有効フレーム数と追加ヒット数を保存する。
            hist = hist + np.load(tmp_cache)["all_hits"]  # 完成rawとtmpを二重計数せず合算する。
            information["temporary_files"] = json.loads(tmp_details.read_text())  # tmpの回収情報を合算の記録に添える。
            information["completed_raw_strip_clock_hits"] = information["raw_strip_clock_hits"]  # 完成rawだけの件数も残す。
            information["raw_strip_clock_hits"] = int(hist.sum())  # 全ヒット数をtmp込みへ更新する。
            information["all_events_completed_raw"] = information.pop("all_events")  # tmp追加後の全イベント数と誤解しない名前にする。
            information["events_with_hits_completed_raw"] = information.pop("events_with_hits")  # 完成rawのイベント数を明確に区別する。
        total += hist  # 異なるランのイベント番号を混ぜずヒストグラムだけを合算する。
        information["colors"] = draw(hist, config, cache_dir / name, name)  # 個別ランの両座標・両カラー方式を保存する。
        metadata.append(information)  # 合算の構成要素を記録する。
        print(f"Saved all-hit accumulation: {name}", flush=True)  # ランごとの完了を報告する。
    reference_path = repo_path(config["inputs"]["run06_07_track_histogram"])  # 既に走査済みのrun06・07の全ヒット積算を指定する。
    reference = np.load(reference_path)  # run06・07のrawを読み直さず積算を使う。
    old_config = json.loads(str(reference["config"]))  # 保存時の較正値を読む。
    for key in ("velocity", "clock_offset", "z_gem", "pitch_strip"):  # mm変換と座標ビンの整合性を確認する。
        if old_config[key] != config[key]:  # 較正値が異なる積算を黙って合算しない。
            raise ValueError(f"Reference calibration mismatch: {key}")  # 不一致の項目を表示する。
    total += reference["all_hits"]  # カット後でなくカット前のrun06・07を合算する。
    result_path = output_file(config, "combined_histogram")  # 全8ランの保存先をJSONから読む。
    np.savez_compressed(result_path, all_hits=total, config=json.dumps(config))  # 全データの積算を保存する。
    prefix = output_file(config, "combined_figure_prefix")  # 合算図の名前をJSONから読む。
    colors = draw(total, config, prefix, "20261124 run00 .. run07")  # 全ランのmm・クロック版を作る。
    result = {"selection": "None; all recorded raw TPC hits; raw clocks 0..1023", "config": config, "new_runs": metadata, "run06_07_source": str(reference_path), "total_strip_clock_hits": int(total.sum()), "colors": colors}  # 条件と保存したヒット総数を明示する。
    output_file(config, "combined_metadata").write_text(json.dumps(result, indent=2) + "\n")  # 再描画と検証に必要な情報を保存する。
    print(f"Saved all-run accumulation: {result_path}", flush=True)  # 全データの保存完了を知らせる。


if __name__ == "__main__":  # importだけでrawを読み始めない。
    main()  # コマンド実行時にだけ処理する。
