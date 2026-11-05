#!/usr/bin/env python3

import multiprocessing
import time
from datetime import datetime

import numpy as np

from plot_temp import run_temperature_plot
from sub_scr import send_rbcp_cmd
from sub_scr.get_info import get_unit_list, load_psa_config, load_run_config
from sub_scr.read_temp import temperature_daq

RUN_FILE = "./run_file.json"
TEMP_INTERVAL_SEC = 3
HV_INTERVAL_SEC = 3
BAR_LENGTH = 50
DEFAULT_TEMP_A = 3.8896501447771166
DEFAULT_TEMP_B = 1232.017507557546
DEFAULT_HV_A = 0.01564
DEFAULT_HV_B = -25.95474


def return_HV(val_0x, param_a, param_b):
    return param_a * int(val_0x, 16) + param_b


def unit_param(unit, group, key, default):
    return unit.get(group, {}).get(key, default)


def unit_serial(unit):
    return unit.get("PSA_SN") or unit.get("FPGA_SN") or "None"


def make_bar(current_val):
    current_pos = int(current_val / 0x3FFF * BAR_LENGTH)
    current_pos = max(0, min(current_pos, BAR_LENGTH))

    bar = [" "] * (BAR_LENGTH + 1)
    for i in range(current_pos + 1):
        bar[i] = "="
    bar[current_pos] = "V"
    return "".join(bar)


def start_temperature_processes(units, latest_temp, temp_output_dir, delay_sec):
    processes = []

    for unit in units:
        psa_sn = unit_serial(unit)
        ip = unit["ip_address"]
        temp_a = unit_param(unit, "temperature_parameter", "a", DEFAULT_TEMP_A)
        temp_b = unit_param(unit, "temperature_parameter", "b", DEFAULT_TEMP_B)
        hv_a = unit_param(unit, "hv_parameter", "a", DEFAULT_HV_A)
        hv_b = unit_param(unit, "hv_parameter", "b", DEFAULT_HV_B)

        latest_temp[psa_sn] = np.nan
        process = multiprocessing.Process(
            target=temperature_daq,
            args=(
                ip,
                psa_sn,
                temp_a,
                temp_b,
                hv_a,
                hv_b,
                latest_temp,
                TEMP_INTERVAL_SEC,
                temp_output_dir,
            ),
        )
        process.start()
        processes.append(process)
        time.sleep(delay_sec)

    return processes


def start_plot_process(units, temp_output_dir, start_at):
    serial_to_ip = {str(unit_serial(unit)): str(unit["ip_address"]) for unit in units}
    process = multiprocessing.Process(
        target=run_temperature_plot,
        args=(temp_output_dir, serial_to_ip, start_at),
    )
    process.start()
    return process


def stop_processes(processes):
    for process in processes:
        if process.is_alive():
            process.terminate()

    for process in processes:
        process.join()


def main():
    monitor_start_at = datetime.now()

    run_config = load_run_config(RUN_FILE)
    psa_config = load_psa_config(run_config)
    units = sorted(
        get_unit_list(run_config, psa_config), key=lambda x: x["position"]["index"]
    )

    temp_output_dir = run_config["temp_data_dir"]
    delay_sec = HV_INTERVAL_SEC / len(units)

    manager = multiprocessing.Manager()
    latest_temp = manager.dict()
    temp_history = {unit_serial(unit): [] for unit in units}

    process_list = start_temperature_processes(
        units,
        latest_temp,
        temp_output_dir,
        delay_sec,
    )
    process_list.append(start_plot_process(units, temp_output_dir, monitor_start_at))

    try:
        while True:
            lines = []

            for unit in units:
                now = datetime.now().strftime("%H:%M:%S")
                psa_sn = unit_serial(unit)
                ip = unit["ip_address"]
                index = unit["position"]["index"]
                hv_a = unit_param(unit, "hv_parameter", "a", DEFAULT_HV_A)
                hv_b = unit_param(unit, "hv_parameter", "b", DEFAULT_HV_B)

                response = send_rbcp_cmd.parse_and_execute(
                    "rd 0x0080 2", ip=ip, verbose=False
                )

                if response is None:
                    lines.append(f"[{index:02d}] {ip}  {now}")
                    lines.append("[ NO RESPONSE ]")
                    lines.append("")
                    time.sleep(delay_sec)
                    continue

                current_hex = hex(int.from_bytes(response[8:], "big"))
                current_val = int(current_hex, 16)
                current_hv = return_HV(current_hex, hv_a, hv_b)
                bar_str = make_bar(current_val)

                temp = latest_temp[psa_sn]
                if np.isnan(temp):
                    temp_str = "TEMP : NaN"
                else:
                    temp_history[psa_sn].append(temp)
                    if len(temp_history[psa_sn]) > 50:
                        temp_history[psa_sn].pop(0)

                    diff = temp - np.mean(temp_history[psa_sn])
                    trend = "(=)"
                    if diff > 0.25:
                        trend = "(+)"
                    elif diff < -0.25:
                        trend = "(-)"
                    temp_str = f"TEMP : {temp:6.2f} degC {trend}"

                lines.append(f"[{index:02d}] {ip}  {now}")
                lines.append(
                    f"[{bar_str}] "
                    f"HV : {current_hv:6.1f} V "
                    f"({current_hex})   "
                    f"{temp_str}"
                )
                lines.append("")
                time.sleep(delay_sec)

            print("\033[H\033[J", end="")
            print("\n".join(lines), end="", flush=True)

    except KeyboardInterrupt:
        pass
    finally:
        stop_processes(process_list)


if __name__ == "__main__":
    main()
