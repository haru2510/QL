#!/usr/bin/env python3

import sys,time,argparse
from sub_scr import send_rbcp_cmd
from sub_scr.get_info import load_run_config,load_psa_config,get_unit_list


RUN_FILE="./run_file.json"


parser=argparse.ArgumentParser(description="RBCP Command Sender")

# parser.add_argument("--hit",type=lambda x:int(x,0),default=0x5000,help="Hit threshold value")
parser.add_argument("--hit",type=lambda x:int(x,0),default=0x4ccc,help="Hit threshold value")
parser.add_argument("--data",type=lambda x:int(x,0),default=0x0a28,help="Data threshold value")
# parser.add_argument("--data",type=lambda x:int(x,0),default=0x0a28,help="Data threshold value")

args=parser.parse_args()

hit_th=args.hit
data_th=args.data


if not (0x0<=hit_th<=0xffff):
    sys.exit(f'err. wrong hit_th {hex(hit_th)}')

if not (0x0<=data_th<=0xffff):
    sys.exit(f'err. wrong data_th {hex(data_th)}')


run_config=load_run_config(RUN_FILE)
psa_config=load_psa_config(run_config)

units=get_unit_list(run_config,psa_config)


default_cmd=[
    f"wrs 0x0100 {hex(hit_th)}",
    f"wrs 0x0102 {hex(hit_th)}",
    f"wrs 0x0104 0xffff",
    f"wrs 0x0106 0xffff",
    f"wrs 0x0108 0xffff",
    f"wrs 0x010a 0xffff",
    f"wrs 0x0200 {hex(data_th)}",
    f"wrs 0x0202 {hex(data_th)}",
    f"wrs 0x0204 {hex(data_th)}",
    f"wrs 0x0206 {hex(data_th)}",
    f"wrs 0x0208 {hex(data_th)}",
    f"wrs 0x020a {hex(data_th)}",
]

for unit in units:

    psa_sn = unit["FPGA_SN"]
    ip = unit["ip_address"]

    print(f"\n=== {psa_sn} ({ip}) ===")

    for cmd_i in default_cmd:

        print(f"input    : {cmd_i}")

        response = send_rbcp_cmd.parse_and_execute(
            cmd_i,
            ip=ip,
            verbose=True
        )

        time.sleep(0.2)
