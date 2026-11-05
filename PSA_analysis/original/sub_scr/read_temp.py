#!/usr/bin/env python3

import socket
import struct
import numpy as np
import time
import datetime
import os


RBCP_VER = 0xFF
RBCP_CMD_WR = 0x80
RBCP_CMD_RD = 0xC0
UDP_BUF_SIZE = 2048
DEFAULT_PORT = 4660
FILE_ROTATE_SEC = 3600


def send_rbcp_command(ip, port, command, address, length, data=b''):

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    sock.settimeout(1.0)

    header = struct.pack(
        "!BBBBI",
        RBCP_VER,
        command,
        0x01,
        length,
        address
    )

    packet = header + data

    sock.sendto(packet, (ip, port))

    try:

        response = sock.recv(UDP_BUF_SIZE)

        return response

    except socket.timeout:

        return None

    finally:

        sock.close()


def read_temp_adc(ip, port=DEFAULT_PORT):

    response = send_rbcp_command(
        ip,
        port,
        RBCP_CMD_RD,
        0x0300,
        2
    )

    if response is None:
        return None

    return int.from_bytes(response[8:], "big")


def read_hv_raw(ip, port=DEFAULT_PORT):

    response = send_rbcp_command(
        ip,
        port,
        RBCP_CMD_RD,
        0x0080,
        2
    )

    if response is None:
        return None

    return int.from_bytes(response[8:], "big")


def convert_temp(adc, a, b):

    return (-b + adc) / a


def convert_hv(raw, a, b):

    return a * raw + b


def temperature_daq(
    ip,
    psa_sn,
    temp_a,
    temp_b,
    hv_a,
    hv_b,
    latest_temp,
    interval_sec=5,
    output_dir="/home/smile/data/tempdata/"
):

    os.makedirs(output_dir, exist_ok=True)

    file_start = datetime.datetime.now()

    filename = os.path.join(
        output_dir,
        f"{datetime.datetime.strftime(file_start,'%Y%m%d-%H%M%S')}_{psa_sn}.dat"
    )

    temp_data = np.full(60, np.nan)

    latest_temp[psa_sn] = np.nan

    while True:

        now = datetime.datetime.now()

        if (now - file_start).total_seconds() > FILE_ROTATE_SEC:

            file_start = now

            filename = os.path.join(
                output_dir,
                f"{datetime.datetime.strftime(file_start,'%Y%m%d-%H%M%S')}_{psa_sn}.dat"
            )

        temp_raw = read_temp_adc(ip)

        hv_raw = read_hv_raw(ip)

        if temp_raw is not None:

            temp = convert_temp(
                temp_raw,
                temp_a,
                temp_b
            )

            latest_temp[psa_sn] = temp

            temp_data[:-1] = temp_data[1:]

            temp_data[-1] = temp

        else:

            temp = np.nan

            latest_temp[psa_sn] = np.nan

        if hv_raw is not None:

            hv = convert_hv(
                hv_raw,
                hv_a,
                hv_b
            )

        else:

            hv = np.nan

        with open(filename, "a") as f:

            f.write(
                f"{now},"
                f"{temp:.4f},"
                f"{temp_raw},"
                f"{hv:.4f},"
                f"{hv_raw}\n"
            )

        time.sleep(interval_sec)


