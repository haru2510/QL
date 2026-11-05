#!/usr/bin/env python3

import os
import csv
import json
import numpy as np


DEFAULT_RUN_FILE="./run_file.json"

TEMP_KEYS=["TL","TR","ML","MR","BL","BR"]
DAC_KEYS =["DAC1","DAC2"]


DEFAULT_RUN_CONFIG={

    "workdir":"./",

    "data_dir":"/home/smile/data",

    "temp_data_dir":"/home/smile/data/tempdata",

    "psa":{
        "ip":[-1]
    }
}


DEFAULT_TEMP_PARAM={

    "sn":"None",

    "a":3.8896501447771166,
    "b":1232.017507557546
}


DEFAULT_SINGLE_UNIT={

    "FPGA_SN":"None",
    "ADC_SN" :"None",

    "ip_address":"192.168.10.16",

    "scinti_thin":"",

    "position":{
        "surface":"TOP",
        "index":-1
    },

    #
    # backward compatibility
    # currently TL is used
    #
    "temperature_parameter":{

        "a":DEFAULT_TEMP_PARAM["a"],
        "b":DEFAULT_TEMP_PARAM["b"]
    },

    #
    # full temperature parameters
    #
    "temperature_parameters":{

        k:DEFAULT_TEMP_PARAM.copy()

        for k in TEMP_KEYS
    },

    "hv_parameter":{
        "a":0.01564,
        "b":-25.95474
    },

    "dac_parameter":{

        k:{
            "sn":"None",
            "a":0.0,
            "b":0.0
        }

        for k in DAC_KEYS
    }
}


def safe_float(v,default=0.0):

    if v is None:
        return default

    v=str(v).strip()

    if v=="":
        return default

    return float(v)


def load_run_config(run_file=DEFAULT_RUN_FILE):

    if not os.path.exists(run_file):
        return DEFAULT_RUN_CONFIG

    with open(run_file) as f:
        return json.load(f)


def load_psa_config(run_config):

    psa_param_file=run_config.get(
        "psa_param_dir",
        run_config.get("psa_param_file",None)
    )

    if psa_param_file is None:

        print("[WARN] psa_param_file is not defined")

        return [DEFAULT_SINGLE_UNIT]

    if not os.path.exists(psa_param_file):

        print(f"[WARN] PSA param file not found : {psa_param_file}")

        return [DEFAULT_SINGLE_UNIT]

    print(f"[INFO] PSA param file : {psa_param_file}")

    psa_config=[]

    with open(psa_param_file) as f:

        for row in csv.DictReader(f):

            temp_params={

                k:{
                    "sn":row.get(f"{k}_SN","None"),

                    "a":safe_float(
                        row.get(f"{k}_A"),
                        DEFAULT_TEMP_PARAM["a"]
                    ),

                    "b":safe_float(
                        row.get(f"{k}_B"),
                        DEFAULT_TEMP_PARAM["b"]
                    )
                }

                for k in TEMP_KEYS
            }

            dac_params={

                k:{
                    "sn":row.get(f"{k}_SN","None"),

                    "a":safe_float(
                        row.get(f"{k}_A"),
                        0.0
                    ),

                    "b":safe_float(
                        row.get(f"{k}_B"),
                        0.0
                    )
                }

                for k in DAC_KEYS
            }

            psa_config.append({

                "FPGA_SN":row.get("FPGA_SN","None"),
                "ADC_SN" :row.get("ADC_SN","None"),

                "ip_address":row.get("IP",""),

                "scinti_thin":row.get(
                    "SCINTI_THIN",
                    ""
                ),

                "position":{

                    "surface":row.get(
                        "SURFACE",
                        ""
                    ),

                    "index":int(
                        row.get(
                            "SLOT",
                            -1
                        ) or -1
                    )
                },

                #
                # backward compatibility
                #
                "temperature_parameter":{

                    "a":temp_params["TL"]["a"],
                    "b":temp_params["TL"]["b"]
                },

                #
                # full temperature parameters
                #
                "temperature_parameters":temp_params,

                "hv_parameter":{

                    "a":safe_float(
                        row.get("HV_A"),
                        0.01564
                    ),

                    "b":safe_float(
                        row.get("HV_B"),
                        -25.95474
                    )
                },

                "dac_parameter":dac_params
            })

    print(f"[INFO] PSA units : {len(psa_config)}")

    return psa_config


def get_ip_list(run_config):

    return run_config["psa"]["ip"]


def get_unit_config(ip_last,psa_config):

    if int(ip_last)==-1:
        return DEFAULT_SINGLE_UNIT

    ip_addr=f"192.168.100.{ip_last}"

    for unit in psa_config:

        if unit["ip_address"]==ip_addr:
            return unit

    raise ValueError(
        f"unit not found : {ip_addr}"
    )


def get_unit_list(run_config,psa_config):

    return [

        get_unit_config(ip,psa_config)

        for ip in get_ip_list(run_config)
    ]


def get_temp_calibration(unit):

    temp_params=unit["temperature_parameters"]

    temp_a=np.array([

        temp_params["TL"]["a"],
        temp_params["TR"]["a"],
        temp_params["ML"]["a"],
        temp_params["MR"]["a"],
        temp_params["BL"]["a"],
        temp_params["BR"]["a"]

    ],dtype=np.float32)[::-1]

    temp_b=np.array([

        temp_params["TL"]["b"],
        temp_params["TR"]["b"],
        temp_params["ML"]["b"],
        temp_params["MR"]["b"],
        temp_params["BL"]["b"],
        temp_params["BR"]["b"]

    ],dtype=np.float32)[::-1]

    return np.ascontiguousarray(temp_a),np.ascontiguousarray(temp_b)