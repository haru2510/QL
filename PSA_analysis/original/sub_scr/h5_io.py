#!/usr/bin/env python3

import os,h5py
import numpy as np


def build_time_array(clk):

    time_arr=(clk.astype(np.int64)-int(clk[0]))/1e6

    wrap=((2**32)-1)/1e6

    time_arr[time_arr<0]+=wrap

    return time_arr


def build_hit_mask(hit):

    return np.where(hit==1)


def build_hit_data(decoded):

    trig     =decoded["trig"]
    clk      =decoded["clk"]
    hit      =decoded["hit"]

    temp_raw =decoded["temp_raw"]
    temp_corr=decoded["temp_corr"]

    ped      =decoded["ped"]
    wave_sum =decoded["wave_sum"]

    hg=decoded["hg"]
    lg=decoded["lg"]

    x_hg=decoded["x_hg"]
    y_hg=decoded["y_hg"]

    x_lg=decoded["x_lg"]
    y_lg=decoded["y_lg"]

    peak_id=decoded["peak_id"]

    hitnum=decoded["hitnum"]

    satu_h=decoded["satu_h"]
    satu_l=decoded["satu_l"]

    time_arr=build_time_array(clk)

    ev_idx,ch_idx=build_hit_mask(hit)

    return {
        "event_time":time_arr,

        "ev_idx":ev_idx,
        "ch_idx":ch_idx,

        "trig_sel" :trig[ev_idx],
        "time_sel" :time_arr[ev_idx],

        "ch_sel":5-ch_idx,

        "hg_sel":hg[ev_idx,ch_idx],
        "lg_sel":lg[ev_idx,ch_idx],

        "x_hg_sel":x_hg[ev_idx,ch_idx],
        "y_hg_sel":y_hg[ev_idx,ch_idx],

        "x_lg_sel":x_lg[ev_idx,ch_idx],
        "y_lg_sel":y_lg[ev_idx,ch_idx],

        "ped_sel" :ped[ev_idx,ch_idx],
        "wave_sel":wave_sum[ev_idx,ch_idx],

        "temp_sel"    :temp_corr[ev_idx],
        "temp_raw_sel":temp_raw[ev_idx],

        "peak_id_sel":peak_id[ev_idx,ch_idx],

        "satu_h_sel":satu_h[ev_idx,ch_idx],
        "satu_l_sel":satu_l[ev_idx,ch_idx]
    }


def write_main_h5(output_path,hit_data):

    with h5py.File(output_path,'w') as f:

        f.create_dataset('trig'      ,data=hit_data["trig_sel"])
        f.create_dataset('event_time',data=hit_data["time_sel"])

        f.create_dataset('ch'    ,data=hit_data["ch_sel"])
        f.create_dataset('hitnum',data=hit_data["hitnum"])

        f.create_dataset('wave_sum',data=hit_data["wave_sel"])

        f.create_dataset('temp'    ,data=hit_data["temp_sel"])
        f.create_dataset('temp_raw',data=hit_data["temp_raw_sel"])

        f.create_dataset('x_hg',data=hit_data["x_hg_sel"])
        f.create_dataset('y_hg',data=hit_data["y_hg_sel"])

        f.create_dataset('x_lg',data=hit_data["x_lg_sel"])
        f.create_dataset('y_lg',data=hit_data["y_lg_sel"])

        f.create_dataset('hg',data=hit_data["hg_sel"])
        f.create_dataset('lg',data=hit_data["lg_sel"])

        f.create_dataset('ped',data=hit_data["ped_sel"])

        f.create_dataset('peak_id',data=hit_data["peak_id_sel"])

        f.create_dataset('satu_h',data=hit_data["satu_h_sel"])
        f.create_dataset('satu_l',data=hit_data["satu_l_sel"])

        f.attrs["temp_order"]="TL,TR,ML,MR,BL,BR"


def write_adc_h5(adc_output_path,decoded,hit_data):

    with h5py.File(adc_output_path,'w') as f:

        f.create_dataset('trig'      ,data=decoded["trig"])
        f.create_dataset('event_time',data=hit_data["event_time"])

        f.create_dataset('hit'   ,data=decoded["hit"])
        f.create_dataset('hitnum',data=decoded["hitnum"])

        f.create_dataset('temp'    ,data=decoded["temp_corr"])
        f.create_dataset('temp_raw',data=decoded["temp_raw"])

        f.create_dataset('ped',data=decoded["ped"])

        f.create_dataset('peak_id',data=decoded["peak_id"])

        f.create_dataset('satu_h',data=decoded["satu_h"])
        f.create_dataset('satu_l',data=decoded["satu_l"])

        f.create_dataset('adc',data=decoded["adc_wave"])

        f.attrs["temp_order"]="TL,TR,ML,MR,BL,BR"


def save_h5(decoded,output_path,keep_adc=False):

    output_dir =os.path.dirname(output_path)
    output_stem=os.path.splitext(os.path.basename(output_path))[0]

    adc_output_path=os.path.join(
        output_dir,
        f"{output_stem}_adc.h5"
    )

    hit_data=build_hit_data(decoded)

    write_main_h5(output_path,hit_data)

    if keep_adc:

        write_adc_h5(
            adc_output_path,
            decoded,
            hit_data
        )

        return output_path,adc_output_path

    return output_path