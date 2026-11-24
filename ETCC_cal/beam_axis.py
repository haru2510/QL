#!/usr/bin/env python3
"""Z断面の局所ピーク中心を直線fitする。座標・幅はmm、傾きはmm/mm。"""
import argparse  # 頂点CSVと結果の保存先を指定する。
import json  # fitの数値と条件を保存する。
from pathlib import Path  # 実行ディレクトリに依存しないパスを扱う。
import numpy as np  # 配列・ヒストグラム・重み付き直線fitを計算する。
from scipy.optimize import least_squares  # Poisson偏差を最小化する。
from scipy.special import erf  # Gaussianを各ビンの境界間で積分する。
import matplotlib  # 画像保存用の描画環境を選ぶ。
matplotlib.use("Agg")  # GUIを開かず結果を保存する。
import matplotlib.pyplot as plt  # 全体図・局所図・ピーク位置を描く。
from matplotlib.colors import LogNorm  # ヒット数の広い範囲を対数色で表す。
from settings import DEFAULT_CONFIG, load_config, output_file  # 入力CSVと結果の保存先を外部JSONから読む。

REFERENCE_Z = 155.0  # 軸位置を表す基準Z [mm]。


def peak_fit(points, roi=(40, 90, 8, 44), step=2.0):  # 1断面のGaussian信号と滑らかな背景をfitする。
    x_edges = np.arange(roi[0], roi[1] + 0.01, step)  # Xの2 mmビン境界を作る。
    y_edges = np.arange(roi[2], roi[3] + 0.01, step)  # Yの2 mmビン境界を作る。
    counts = np.histogram2d(points[:, 0], points[:, 1], bins=(x_edges, y_edges))[0]  # ROI内の頂点を数える。
    x = (x_edges[1:] + x_edges[:-1]) / 2  # 背景計算用のXビン中心を求める。
    y = (y_edges[1:] + y_edges[:-1]) / 2  # 背景計算用のYビン中心を求める。
    if counts.sum() == 0:  # 空の断面を無意味なfit結果にしない。
        raise ValueError("No vertices in the peak-fit ROI")  # 入力範囲の問題を明示する。

    def expected(parameters):  # 8パラメータから各ビンの期待頂点数を求める。
        signal, mx, my, sx, sy, background, gx, gy = parameters  # 信号数・中心・幅・背景・背景勾配を分ける。
        px = np.diff(0.5 * (1 + erf((x_edges - mx) / (np.sqrt(2) * sx))))  # X方向Gaussianを各ビンで積分する。
        py = np.diff(0.5 * (1 + erf((y_edges - my) / (np.sqrt(2) * sy))))  # Y方向Gaussianを各ビンで積分する。
        core = signal * px[:, None] * py[None, :]  # 無相関の2次元Gaussian信号を作る。
        exponent = gx * (x[:, None] - 65) + gy * (y[None, :] - 25)  # 背景の緩いX/Y勾配を表す。
        return core + background * np.exp(exponent)  # 常に正の背景を信号に足す。

    def residual(parameters):  # Poisson偏差を最小二乗用の残差に変換する。
        mean = np.maximum(expected(parameters), 1e-10)  # log(0)と0除算を防ぐ。
        log_term = counts * np.log(np.maximum(counts, 1e-10) / mean)  # 0件のビンも扱える対数項を計算する。
        deviance = 2 * (mean - counts + log_term)  # 各ビンのPoisson偏差を求める。
        signed = np.sign(mean - counts) * np.sqrt(np.maximum(deviance, 0))  # 数値丸めの負値を0にして符号付き残差を作る。
        return signed.ravel()  # 最適化関数が受け取る1次元配列にする。

    background = max(float(np.median(counts)), 0.1)  # 背景の初期値を中央値から決める。
    signal = max(float(counts.sum() - background * counts.size), 10)  # 背景を超える成分を信号の初期値にする。
    initial = [signal, 63, 22, 5, 4, background, 0, 0]  # 観測されたピーク付近から最適化を始める。
    lower = [0, 45, 14, 0.6, 0.6, 0.001, -0.07, -0.07]  # 信号数・幅・背景の下限と中心探索範囲を指定する。
    upper = [1e7, 80, 34, 18, 18, 1e5, 0.07, 0.07]  # 背景への吸い込みを避けるため探索上限も指定する。
    fit = least_squares(residual, initial, bounds=(lower, upper), max_nfev=1000)  # ビン積分したモデルで中心を求める。
    if not fit.success:  # 収束しない断面を有効な中心として採用しない。
        raise RuntimeError(f"Peak fit failed: {fit.message}")  # 最適化の失敗理由を示す。
    dof = counts.size - len(initial)  # fitの自由度を求める。
    deviance = float(np.sum(fit.fun**2))  # 最終的なPoisson偏差を合計する。
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac)  # 残差の局所曲率からパラメータの共分散を見積もる。
    covariance *= max(1, deviance / dof)  # モデルからのずれが大きい場合は誤差を膨らませる。
    errors = np.sqrt(np.maximum(np.diag(covariance), 0))  # 各パラメータの誤差を取り出す。
    return {  # 描画と直線fitに必要な数値を返す。
        "x_mm": float(fit.x[1]),  # 断面のXピーク中心 [mm]。
        "y_mm": float(fit.x[2]),  # 断面のYピーク中心 [mm]。
        "sx_mm": float(fit.x[3]),  # GaussianコアのX幅 [mm]。
        "sy_mm": float(fit.x[4]),  # GaussianコアのY幅 [mm]。
        "x_error_mm": float(errors[1]),  # X中心のモデル上の誤差 [mm]。
        "y_error_mm": float(errors[2]),  # Y中心のモデル上の誤差 [mm]。
        "signal_counts": float(fit.x[0]),  # Gaussian信号の積分頂点数。
        "deviance": deviance,  # 断面fitの適合度を確認する値。
        "ndf": dof,  # 偏差を評価する自由度。
        "parameters": fit.x.tolist(),  # 全8パラメータを後から検証できるよう残す。
    }  # 断面ごとの結果を完成させる。


def fit_slices(vertices, roi=(40, 90, 8, 44), shift=0):  # Zを20 mmごとに区切って各XYピークを求める。
    rows = []  # 各断面の結果を保存する。
    for low in range(15 + shift, 295, 20):  # 通常は15〜295 mm、感度確認では境界を10 mmずらす。
        selected = vertices[(vertices[:, 2] >= low) & (vertices[:, 2] < low + 20)]  # このZ範囲に入る全頂点を選ぶ。
        result = peak_fit(selected, roi)  # 局所ピークの信号と背景を分離する。
        result.update(z_mm=low + 10, z_low_mm=low, z_high_mm=low + 20, vertices=len(selected))  # 断面の中心・範囲・入力数を添える。
        rows.append(result)  # Z順に結果を蓄える。
    return rows  # 断面のピーク中心列を返す。


def fit_line(rows, axis):  # X(Z)またはY(Z)を重み付き直線でfitする。
    z = np.array([row["z_mm"] for row in rows])  # 各断面の中心Z [mm]を並べる。
    values = np.array([row[f"{axis}_mm"] for row in rows])  # その座標のピーク中心 [mm]を並べる。
    errors = np.array([row[f"{axis}_error_mm"] for row in rows])  # 断面fitの中心誤差 [mm]を並べる。
    errors = np.maximum(errors, 0.05)  # 極小誤差の断面だけで軸が決まらないよう0.05 mmを下限にする。
    design = np.column_stack((np.ones(len(z)), z - REFERENCE_Z))  # 位置と傾きの2パラメータの設計行列を作る。
    weights = 1 / errors**2  # 中心誤差の逆分散を重みにする。
    covariance = np.linalg.inv(design.T @ (weights[:, None] * design))  # 重み付き直線fitの共分散を求める。
    parameters = covariance @ (design.T @ (weights * values))  # 基準Zでの位置と傾きを求める。
    chi2 = float(np.sum(((values - design @ parameters) / errors)**2))  # 直線からの残差を評価する。
    covariance *= max(1, chi2 / (len(z) - 2))  # 直線で説明できない変動を誤差に反映する。
    return {  # 数値の意味と単位が分かるキー名で保存する。
        "at_z155_mm": float(parameters[0]),  # Z=155 mmでの軸位置 [mm]。
        "slope_mm_per_mm": float(parameters[1]),  # Zが1 mm増えたときの座標変化 [mm/mm]。
        "errors": np.sqrt(np.diag(covariance)).tolist(),  # 位置[mm]と傾き[mm/mm]の誤差。
        "chi2": chi2,  # 直線としての適合度。
        "ndf": len(z) - 2,  # 位置と傾きの2パラメータを引いた自由度。
    }  # 直線fitの結果を返す。


def angle(x_line, y_line):  # 3次元軸とZ軸の角度を求める。
    transverse = np.hypot(x_line["slope_mm_per_mm"], y_line["slope_mm_per_mm"])  # 横方向の傾きを合成する。
    return float(np.degrees(np.arctan(transverse)))  # ラジアンを度に直して返す。


def draw(vertices, result, output, local=False):  # 全体図と局所表示を同じfit結果から描く。
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)  # 投影3枚とfit診断3枚を用意する。
    z_line = np.linspace(15, 295, 200)  # 通常のfit範囲内でだけ軸を描く。
    x_line, y_line = result["x_line"], result["y_line"]  # 求めた2方向の直線を取り出す。
    x = x_line["at_z155_mm"] + x_line["slope_mm_per_mm"] * (z_line - REFERENCE_Z)  # 軸のX座標 [mm]を計算する。
    y = y_line["at_z155_mm"] + y_line["slope_mm_per_mm"] * (z_line - REFERENCE_Z)  # 軸のY座標 [mm]を計算する。
    projections = [(0, 1, "XY"), (0, 2, "XZ"), (1, 2, "YZ")]  # 見る座標平面を指定する。
    for column, (i, j, title) in enumerate(projections):  # 各投影を同じ順番で描く。
        selected = vertices  # 標準表示ではXYの追加選別をしない。
        if local and title == "XZ":  # 局所版のXZだけ帯状の表示範囲を使う。
            selected = vertices[np.abs(vertices[:, 1] - 22) < 15]  # 表示用にY=7〜37 mmを選ぶ。
            title += " display only: |Y - 22| < 15 mm"  # イベント除外やfit範囲との違いをラベルに書く。
        if local and title == "YZ":  # 局所版のYZだけ別の帯状範囲を使う。
            selected = vertices[np.abs(vertices[:, 0] - 63) < 20]  # 表示用にX=43〜83 mmを選ぶ。
            title += " display only: |X - 63| < 20 mm"  # 条件を図から読めるようにする。
        edges_x = np.linspace(-164, 164, 83)  # 検出器全体を含む横方向の表示範囲を使う。
        edges_y = np.linspace(-164, 164, 83) if j == 1 else np.linspace(-50, 365, 84)  # 負Zと高Z側も表示範囲に含める。
        counts = np.histogram2d(selected[:, i], selected[:, j], bins=(edges_x, edges_y))[0]  # 表示する頂点数を数える。
        ax = axes[0, column]  # 上段の対応するパネルを選ぶ。
        image = ax.pcolormesh(edges_x, edges_y, np.ma.masked_less(counts.T, 1), norm=LogNorm(vmin=1), cmap="turbo")  # 空ビンをマスクし対数色で描く。
        fig.colorbar(image, ax=ax, label="vertices / bin")  # 色の意味を表示する。
        ax.plot(x if i == 0 else y, y if j == 1 else z_line, "m-", lw=1.5)  # 各平面に同じ3次元軸を射影する。
        ax.set(xlabel=("X", "Y", "Z")[i] + " [mm]", ylabel=("X", "Y", "Z")[j] + " [mm]", title=title)  # 座標の単位を明記する。
        ax.invert_yaxis()  # YとZが増える方向を下にする。
    rows = result["slices"]  # 各Z断面のfit結果を取り出す。
    z = np.array([row["z_mm"] for row in rows])  # 断面中心のZ座標を並べる。
    for column, (coordinate, fitted) in enumerate((("x", x_line), ("y", y_line))):  # XとYのピーク中心変化を別々に描く。
        ax = axes[1, column]  # 下段の中心位置パネルを選ぶ。
        centers = [row[f"{coordinate}_mm"] for row in rows]  # 各断面の中心位置を取り出す。
        errors = [row[f"{coordinate}_error_mm"] for row in rows]  # 断面fitの中心誤差を取り出す。
        ax.errorbar(z, centers, xerr=10, yerr=errors, fmt="o", ms=4, color="black", label="Local Gaussian centre")  # 断面幅と中心誤差を区別して描く。
        prediction = fitted["at_z155_mm"] + fitted["slope_mm_per_mm"] * (z_line - REFERENCE_Z)  # 直線fitの予測位置を求める。
        ax.plot(z_line, prediction, "m-", label="Line fit")  # 中心列への直線を描く。
        ax.set(xlabel="Z [mm]", ylabel=coordinate.upper() + " centre [mm]", title=f"{coordinate.upper()}(Z): local peak centres")  # 全体のPCAと区別する。
        ax.grid(alpha=0.25)  # 小さな位置変化を読みやすくする。
        ax.legend(fontsize=8)  # 点と線の意味を表示する。
    ax = axes[1, 2]  # 右下をピーク幅の診断に使う。
    ax.plot(z, [row["sx_mm"] for row in rows], "o-", label="sigma X")  # XのGaussianコア幅を描く。
    ax.plot(z, [row["sy_mm"] for row in rows], "s-", label="sigma Y")  # YのGaussianコア幅を描く。
    ax.set(xlabel="Z [mm]", ylabel="Gaussian core width [mm]", title="Vertex peak width")  # ビーム自体の幅と断定しない。
    ax.legend()  # XとYの幅を区別する。
    mode = "Local display (cuts labelled)" if local else "Full display (no extra XY selection)"  # 表示版の違いを明記する。
    fig.suptitle(f"{result['run_label']}: local-peak axis; angle from Z = {result['angle_from_z_deg']:.3f} deg\n{mode}; event selection is inherited from input CSV")  # 不要なPCA表記をタイトルから外す。
    fig.savefig(output, dpi=140)  # fit図を指定の保存先へ出力する。
    plt.close(fig)  # 図のメモリを解放する。


def main():  # 入力から数値結果と2種類の図を作る。
    parser = argparse.ArgumentParser(description=__doc__)  # コマンドの説明を用意する。
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)  # 共通の外部JSONを指定する。
    parser.add_argument("--input", type=Path)  # 選別済み頂点CSVを指定する。
    parser.add_argument("--out", type=Path)  # 出力の共通名を指定する。
    args = parser.parse_args()  # 実行時の指定を読む。
    config = load_config(args.config)  # 入出力欄を持つ外部設定を読む。
    args.input = args.input or output_file(config, "vertices_selected")  # 選別済み頂点CSVを入力にする。
    args.out = args.out or output_file(config, "beam_axis_prefix")  # JSONに記した保存名を使う。
    vertices = np.loadtxt(args.input, delimiter=",", skiprows=2, usecols=(5, 6, 7))  # X/Y/Zだけをメモリに読み込む。
    args.out.parent.mkdir(parents=True, exist_ok=True)  # 保存先を用意する。
    rows = fit_slices(vertices)  # 標準のZ断面とROIで局所中心を求める。
    x_line, y_line = fit_line(rows, "x"), fit_line(rows, "y")  # 断面中心に2方向の直線をfitする。
    checks = []  # 解析範囲を変えた感度確認を残す。
    for roi, shift in [((36, 94, 6, 48), 0), ((44, 86, 10, 42), 0), ((40, 90, 8, 44), 10)]:  # XY範囲とZの区切りを少し変える。
        alternative = fit_slices(vertices, roi, shift)  # 同じデータに別の断面条件を適用する。
        lx, ly = fit_line(alternative, "x"), fit_line(alternative, "y")  # 別条件の軸を求める。
        checks.append({"roi_mm": roi, "z_shift_mm": shift, "x": lx, "y": ly, "angle_from_z_deg": angle(lx, ly)})  # 位置と傾きの変動を保存する。
    result = {  # 結果の数値と適用条件をまとめる。
        "method": "Local XY peak in Z slices; bin-integrated 2D Gaussian + exponential background; weighted lines",  # 軸の求め方を明記する。
        "source": str(args.input.resolve()),  # どの選別済みCSVを使ったか記録する。
        "run_label": ", ".join(config["runs"]),  # タイトルのラン名も外部JSONに従う。
        "vertices": len(vertices),  # 入力イベント数を保存する。
        "roi_mm": [40, 90, 8, 44],  # fitに使用する局所XY範囲。
        "xy_bin_mm": 2,  # ピークfitのXYビン幅。
        "slice_thickness_mm": 20,  # Z断面の厚さ。
        "line_reference_z_mm": REFERENCE_Z,  # 軸位置の基準Z。
        "line_center_error_floor_mm": 0.05,  # 直線fitの重みに使う誤差下限。
        "slices": rows,  # 各断面の中心・幅・適合度。
        "x_line": x_line,  # X(Z)の係数と誤差。
        "y_line": y_line,  # Y(Z)の係数と誤差。
        "angle_from_z_deg": angle(x_line, y_line),  # 3次元軸とZ軸の角度。
        "sensitivity_checks": checks,  # 条件変更による変動。
        "display_cuts_only": {"local_XZ_y_mm": [7, 37], "local_YZ_x_mm": [43, 83]},  # 局所図だけの表示条件をfit条件と区別する。
        "caveat": "Statistical/model errors do not include drift calibration, alignment, or vertex-estimator bias.",  # 絶対精度とは区別する。
    }  # 結果辞書を完成させる。
    args.out.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")  # 数値結果をSSDなしでも読める形で保存する。
    draw(vertices, result, args.out.with_suffix(".png"))  # XY追加選別なしの標準図を保存する。
    draw(vertices, result, args.out.with_name(args.out.name + "_local").with_suffix(".png"), local=True)  # 条件を明記した局所表示も保存する。
    print(json.dumps({"x_line": x_line, "y_line": y_line, "angle_from_z_deg": result["angle_from_z_deg"]}, indent=2))  # レビュー用の主要結果を報告する。


if __name__ == "__main__":  # import時には解析を始めない。
    main()  # コマンド実行時だけ全処理を呼び出す。
