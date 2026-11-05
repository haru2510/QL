#!/usr/bin/env python3
import sys,time
from datetime import datetime

from sub_scr import send_rbcp_cmd
from sub_scr.get_info import load_run_config,load_psa_config,get_unit_list


RUN_FILE="./run_file.json"

HV_SAFE_STEP_THRESHOLD=0x1000
HV_MAX_HEX=0x3fff


def convert_HV(hv_val,a,b):

    return int((hv_val-b)/a)


def return_HV(val,a,b):

    return a*val+b


if __name__=="__main__":

    if len(sys.argv)!=2:

        print("err. python HV_cal.py [HV val]")

        sys.exit()

    input_val=sys.argv[1]

    run_config=load_run_config(RUN_FILE)
    psa_config=load_psa_config(run_config)

    units=get_unit_list(run_config,psa_config)

    log_file=f'{run_config["data_dir"]}/logHV.log'

    with open(log_file,"a") as f:

        now=datetime.now().strftime("%Y%m%d-%H%M%S")

        for unit in units:

            psa_sn=unit["FPGA_SN"]

            ip=unit["ip_address"]

            a=unit["hv_parameter"]["a"]
            b=unit["hv_parameter"]["b"]

            hv_max_v=return_HV(HV_MAX_HEX,a,b)

            print(f"\n=== {psa_sn} ({ip}) ===")

            response=send_rbcp_cmd.parse_and_execute("rd 0x0080 2", ip=ip, verbose=False)

            current_val=int.from_bytes(response[8:],"big")

            current_hex=f"0x{current_val:04x}"

            current_hv=return_HV(current_val,a,b)

            #
            # direct hex input
            #
            if input_val.startswith("0x"):

                target_val=int(input_val,16)

                if target_val>HV_MAX_HEX:

                    print(f"[WARN] target exceeds 0x{HV_MAX_HEX:04x}")

                    print(f"max HV = {hv_max_v:.1f} V")

                    continue

            #
            # voltage input
            #
            else:

                hv_val=float(input_val)

                if not (0<=hv_val<=300):

                    print("err. 0 <= val <= 300")

                    continue

                target_val=convert_HV(hv_val,a,b)

                if target_val>HV_MAX_HEX:

                    print(f"[WARN] converted value exceeds 0x{HV_MAX_HEX:04x}")

                    print(f"requested HV : {hv_val:.1f} V")

                    print(f"max HV       : {hv_max_v:.1f} V")

                    continue

            target_hex=f"0x{target_val:04x}"

            target_hv=return_HV(target_val,a,b)

            diff=abs(target_val-current_val)

            diff_flag="  X" if diff>=HV_SAFE_STEP_THRESHOLD else ""

            print(f"{current_hv:.1f}V ({current_hex}) --> {target_hv:.1f}V ({target_hex}) {diff_flag}")

            bar_length=40

            #
            # gradual ramp
            #
            if diff>=HV_SAFE_STEP_THRESHOLD:

                for i in range(1,6):

                    step_val=int(
                        current_val+
                        (target_val-current_val)*i/5
                    )

                    step_hex=f"0x{step_val:04x}"

                    step_hv=return_HV(step_val,a,b)

                    step_pos=int(step_val/HV_MAX_HEX*bar_length)

                    target_pos=int(target_val/HV_MAX_HEX*bar_length)

                    bar=[" "]*(bar_length+1)

                    for j in range(step_pos+1):

                        bar[j]="="

                    bar[step_pos]="C"
                    bar[target_pos]="T"

                    print(f'\r[{"".join(bar)}] {step_hv:6.1f}V ({step_hex})', end="", flush=True)

                    send_rbcp_cmd.parse_and_execute(f"wrs 0x0080 {step_hex}", ip=ip, verbose=False)

                    time.sleep(1)

                print("")

            #
            # direct write
            #
            else:

                step_pos=int(target_val/HV_MAX_HEX*bar_length)

                bar=[" "]*(bar_length+1)

                for j in range(step_pos+1):

                    bar[j]="="

                bar[step_pos]="C"

                print(f'[{"".join(bar)}] {target_hv:.1f}V ({target_hex})')

                send_rbcp_cmd.parse_and_execute(f"wrs 0x0080 {target_hex}", ip=ip, verbose=False)

            f.write(
                f"{now}, "
                f"{psa_sn}, "
                f"{ip}, "
                f"{current_hv:.1f}, "
                f"{target_hv:.1f}, "
                f"{target_hex}\n"
            )

            time.sleep(0.5)