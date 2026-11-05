#!/usr/bin/env python3

import argparse
import datetime
import multiprocessing
from multiprocessing import RLock
import os
import socket
import subprocess
import time

from tqdm import tqdm

from sub_scr import send_rbcp_cmd
from sub_scr.get_info import load_run_config, load_psa_config, get_unit_list


RUN_FILE="./run_file.json"
SITCP_PORT=24
RECV_SIZE=2**14
DEFAULT_MB=100


def parse_args():

    parser=argparse.ArgumentParser()

    parser.add_argument("--size", type=float, default=DEFAULT_MB)
    parser.add_argument("--mins", type=float, default=None)
    parser.add_argument("--filenum", type=int, default=None)
    parser.add_argument("--link", type=str, default="")
    parser.add_argument("-q", "--quiet", action="store_true")

    return parser.parse_args()


def log(message):

    with tqdm.get_lock():

        print(message, flush=True)


def update_link(work_dir, output_file):

    os.makedirs(work_dir, exist_ok=True)

    link_path=os.path.join(work_dir, os.path.basename(output_file))

    if os.path.islink(link_path) or os.path.exists(link_path):
        os.remove(link_path)

    subprocess.run(["ln", "-s", output_file, link_path], check=False)


def build_output_file(save_dir, unit_label, file_index):

    timestamp=datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    return os.path.join(
        save_dir,
        f"{timestamp}_{unit_label}_{file_index}.dat"
    )


def connect_sitcp(ip):

    sock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    sock.settimeout(10)

    sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 2**20)

    sock.connect((ip, SITCP_PORT))

    sock.settimeout(60)

    return sock


def get_data(
    sock,
    save_dir,
    work_dir,
    unit_label,
    max_file_size,
    daq_end_with,
    max_file_num,
    quiet,
    position,
    stop_event
):

    file_index=0

    os.makedirs(save_dir, exist_ok=True)

    output_file=build_output_file(save_dir, unit_label, file_index)

    file=open(output_file, "wb")

    pbar=tqdm(
        total=max_file_size,
        unit="B",
        unit_scale=True,
        desc=unit_label,
        disable=quiet,
        position=position,
        leave=True,
        dynamic_ncols=False
    )

    prev_size=0

    try:

        while not stop_event.is_set():

            if datetime.datetime.now()>=daq_end_with:
                break

            try:

                response=sock.recv(RECV_SIZE)

            except socket.timeout:

                continue

            if not response:
                continue

            file.write(response)

            current_size=file.tell()

            pbar.update(current_size-prev_size)

            prev_size=current_size

            if current_size>=max_file_size:

                file.flush()

                os.fsync(file.fileno())

                update_link(work_dir, output_file)

                file.close()

                pbar.close()

                file_index+=1

                if max_file_num is not None and file_index>=max_file_num:

                    log(f"[{unit_label}] reached max file count : {max_file_num}")

                    break

                output_file=build_output_file(save_dir, unit_label, file_index)

                file=open(output_file, "wb")

                pbar=tqdm(
                    total=max_file_size,
                    unit="B",
                    unit_scale=True,
                    desc=unit_label,
                    disable=quiet,
                    position=position,
                    leave=True,
                    dynamic_ncols=False
                )

                prev_size=0

    except KeyboardInterrupt:

        stop_event.set()

    finally:

        if not file.closed:

            file.flush()

            os.fsync(file.fileno())

            update_link(work_dir, output_file)

            file.close()

        pbar.close()


def run_unit(unit, args, save_dir, work_base_dir, position, stop_event):

    psa_sn=unit["FPGA_SN"]

    sitcp_ip=unit["ip_address"]

    max_file_size=args.size*1024*1024

    if args.mins is None:

        daq_time_limit=datetime.timedelta(days=5)

    else:

        daq_time_limit=datetime.timedelta(minutes=args.mins)

    daq_end_with=datetime.datetime.now()+daq_time_limit

    if args.link=="":

        work_dir=os.path.join(work_base_dir, psa_sn)

    else:

        work_dir=os.path.join(work_base_dir, args.link, psa_sn)

    log(f"[{psa_sn}] reset FIFO")

    send_rbcp_cmd.parse_and_execute(
        "wrs 0x0086 0x00ff",
        ip=sitcp_ip,
        verbose=False
    )

    time.sleep(0.6)

    send_rbcp_cmd.parse_and_execute(
        "wrs 0x0086 0x0000",
        ip=sitcp_ip,
        verbose=False
    )

    sock=None

    while sock is None and not stop_event.is_set():

        try:

            sock=connect_sitcp(sitcp_ip)

            log(f"[{psa_sn}] connected")

        except Exception as exc:

            log(f"[{psa_sn}] {exc}")

            time.sleep(1)

    if stop_event.is_set():
        return

    try:

        get_data(
            sock=sock,
            save_dir=save_dir,
            work_dir=work_dir,
            unit_label=psa_sn,
            max_file_size=max_file_size,
            daq_end_with=daq_end_with,
            max_file_num=args.filenum,
            quiet=args.quiet,
            position=position,
            stop_event=stop_event
        )

    finally:

        if sock is not None:
            sock.close()


def main():

    args=parse_args()

    tqdm.set_lock(RLock())

    run_config=load_run_config(RUN_FILE)

    psa_config=load_psa_config(run_config)

    units=get_unit_list(run_config, psa_config)

    save_dir=run_config["data_dir"]

    work_base_dir=run_config["workdir"]

    stop_event=multiprocessing.Event()

    processes=[]

    try:

        for position, unit in enumerate(units):

            process=multiprocessing.Process(
                target=run_unit,
                args=(unit, args, save_dir, work_base_dir, position, stop_event)
            )

            process.start()

            processes.append(process)

        for process in processes:
            process.join()

    except KeyboardInterrupt:

        print("\nKeyboardInterrupt")

        stop_event.set()

        for process in processes:
            process.join()

    print("------------------ DAQ END -------------------")


if __name__=="__main__":

    main()
