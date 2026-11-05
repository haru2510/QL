#!/usr/bin/env python3
#
# 複数の run の HA の .dat をまとめて .h5 にデコードし、HA ごとのディレクトリに置く。
#
#   出力: <out-dir>/<HA>/<run_id>_<HA>_<番号>.h5   (例: h5_runs/ha06/20261103_run02_ha06_00012.h5)
#
# run が違っても .dat の名前は同じ (ha06_00000.dat) なので、run_id を付けて名前がぶつからないようにする。
# HA ごとにまとめるので、peak_search.py -i <out-dir> で HA ごとに全 run を合わせたピーク探索ができ、
# check_data.py や grid_spect.py にも HA のディレクトリをそのまま渡せる。
#
# デコードは psa_decode.py の process_dat_file をそのまま使う。温度較正は psa_param の HA 列から引く。
# 書き途中で止まった run の .dat.tmp も、読める所まで読む (--no-tmp で除く)。
# 出力が既にあるファイルは飛ばすので、途中で止めても続きから再開できる。
#

import argparse
import glob
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

from natsort import natsorted

# このリポジトリ (analysis_sm3_FM) の場所。既定のパスはここからの相対で決める
REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.insert(0, os.path.dirname(REPO_DIR))
from analysis_sm3_FM.PSA_analysis.src_onoda.psa_decode import (  # noqa: E402
    get_sn_from_filename,
    get_unit_by_sn,
    process_dat_file,
)
from analysis_sm3_FM.PSA_analysis.src_onoda.sub_scr.get_info import (  # noqa: E402
    get_temp_calibration,
    load_psa_config,
)


DEFAULT_OUTPUT_DIR = os.path.join(REPO_DIR, "PSA_analysis", "h5_runs")
DEFAULT_PSA_PARAM = os.path.join(REPO_DIR, "PSA_analysis", "param", "psa_param_20261029.csv")

FILE_PATTERN = re.compile(r"^(ha\d+)_(\d+)\.dat(\.tmp)?$")


def parse_args():

    parser = argparse.ArgumentParser(description="Decode HA .dat files of several runs into per-HA .h5 directories")
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="run directories (glob patterns are allowed), e.g. /mnt/d/UVSOR2026_late/*/20261103_run0[1-3]")
    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--ha", nargs="+", default=None, help="HA to decode, e.g. ha06 ha07 (default: all)")
    parser.add_argument("--psa-param", default=DEFAULT_PSA_PARAM, help="psa_param csv (HA, FPGA_SN, temperature calibration)")
    parser.add_argument("--no-tmp", action="store_true", help="skip .dat.tmp left by runs that did not stop normally")
    parser.add_argument("--workers", type=int, default=4, help="files decoded at the same time")

    return parser.parse_args()


def collect_jobs(args):
    """(入力 .dat, 出力 .h5, HA) の一覧を返す。"""

    directories = natsorted({path for item in args.input for path in glob.glob(item) if os.path.isdir(path)})

    if len(directories) == 0:
        raise FileNotFoundError(f"No run directories found: {args.input}")

    jobs = []

    for directory in directories:

        run_id = os.path.basename(os.path.normpath(directory))

        # 長い run は <run_id>/part0, part1, ... に分かれているので、親の名前を run_id にする
        if re.fullmatch(r"part\d+", run_id):
            run_id = os.path.basename(os.path.dirname(os.path.normpath(directory)))

        for name in natsorted(os.listdir(directory)):

            match = FILE_PATTERN.match(name)

            if match is None or (match.group(3) and args.no_tmp):
                continue

            ha, number = match.group(1), match.group(2)

            if args.ha and ha not in args.ha:
                continue

            out_path = os.path.join(args.out_dir, ha, f"{run_id}_{ha}_{number}.h5")
            jobs.append((os.path.join(directory, name), out_path, ha))

    return directories, jobs


def temperature_calibrations(jobs, psa_param):
    """HA ごとの温度較正 (temp_a, temp_b) を返す。見つからない HA は None (psa_decode の既定値)。"""

    psa_config = load_psa_config({"psa_param_file": psa_param})
    calibrations = {}

    for ha in sorted({ha for _, _, ha in jobs}):

        unit = get_unit_by_sn(get_sn_from_filename(f"{ha}_00000.dat", psa_config), psa_config)

        if unit is None:
            print(f"[WARN] {ha}: not in {psa_param}, default temperature calibration is used")
            calibrations[ha] = (None, None)
        else:
            calibrations[ha] = get_temp_calibration(unit)

    return calibrations


def decode(job):

    src, dst, temp_a, temp_b = job

    # 途中で止めたときに書きかけの .h5 が残って次回飛ばされないよう、書き終えてから置き換える
    partial = f"{dst}.part"

    try:
        process_dat_file(src, partial, temp_a=temp_a, temp_b=temp_b)
    except ValueError as exc:
        # イベントの無い .dat (止まった run の空の .tmp など)
        return f"[SKIP] {src}: {exc}"

    os.replace(partial, dst)

    return None


def main():

    args = parse_args()

    directories, jobs = collect_jobs(args)
    calibrations = temperature_calibrations(jobs, args.psa_param)

    todo = [(src, dst, *calibrations[ha]) for src, dst, ha in jobs if not os.path.exists(dst)]

    print(f"[READ] {len(directories)} runs, {len(jobs)} files, {len(jobs) - len(todo)} already decoded")

    for ha in sorted({ha for _, _, ha in jobs}):
        os.makedirs(os.path.join(args.out_dir, ha), exist_ok=True)

    with ProcessPoolExecutor(args.workers) as pool:
        for index, message in enumerate(pool.map(decode, todo, chunksize=4), start=1):
            if message:
                print(message)
            if index % 200 == 0 or index == len(todo):
                print(f"[DECODE] {index}/{len(todo)}", flush=True)

    for ha in sorted({ha for _, _, ha in jobs}):
        n_files = len(glob.glob(os.path.join(args.out_dir, ha, "*.h5")))
        print(f"[HA {ha}] {n_files} files in {os.path.join(args.out_dir, ha)}")


if __name__ == "__main__":

    main()
