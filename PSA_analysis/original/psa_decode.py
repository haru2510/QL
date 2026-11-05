#!/usr/bin/env python3

import os,sys,time,glob,argparse,h5py,numba
import numpy as np

from numba import njit,prange
from natsort import natsorted

from sub_scr.get_info import (
    load_run_config,
    load_psa_config,
    get_temp_calibration
)

numba.set_num_threads(4)

RUN_FILE="./run_file.json"

DEFAULT_TEMP_A=np.ascontiguousarray(np.full(6,3.925820834,dtype=np.float32))
DEFAULT_TEMP_B=np.ascontiguousarray(np.full(6,1230.29635,dtype=np.float32))

DEFAULT_H5_STORE_DIR="./h5_store"

HEADER=b'\xeb\x90\x19\x64'
FOOTER=b'\x75\x50\x49\x43'

PEDESTAL_RANGE=23

WAVE_FRONT=4
WAVE_BACK=9

WAVE_DATA_START=25
WAVE_DATA_STOP=41


def get_unit_by_sn(sn,psa_config):

    for unit in psa_config:

        if unit.get("FPGA_SN")==sn:
            return unit

    return None


@njit(cache=True)
def find_events_jit(buf,H0,H1,H2,H3,F0,F1,F2,F3):

    n=buf.size

    starts=[]
    ends=[]

    i=0

    while i<n-4:

        if buf[i]==H0 and buf[i+1]==H1 and buf[i+2]==H2 and buf[i+3]==H3:

            j=i+4

            while j<n-4:

                if buf[j]==F0 and buf[j+1]==F1 and buf[j+2]==F2 and buf[j+3]==F3:

                    starts.append(i)
                    ends.append(j+4)

                    i=j+4

                    break

                j+=1

        i+=1

    return starts,ends


def read_file(file_path):

    buf_list=[]

    with open(file_path,'rb') as f:

        while True:

            chunk=f.read(64*1024*1024)

            if not chunk:
                break

            buf_list.append(np.frombuffer(chunk,dtype=np.uint8))

    if len(buf_list)==0:
        return None,None,0

    buf=np.concatenate(buf_list)

    starts,ends=find_events_jit(
        buf,
        HEADER[0],HEADER[1],HEADER[2],HEADER[3],
        FOOTER[0],FOOTER[1],FOOTER[2],FOOTER[3]
    )

    n=len(starts)

    if n==0:
        return None,None,0

    lengths=[ends[i]-starts[i] for i in range(n)]

    max_len=max(lengths)

    events_buf=np.zeros((n,max_len),dtype=np.uint8)

    event_sizes=np.zeros(n,dtype=np.int32)

    for i in range(n):

        L=lengths[i]

        event_sizes[i]=L

        events_buf[i,:L]=buf[starts[i]:ends[i]]

    return events_buf,event_sizes,n


@njit(cache=True,fastmath=True,nogil=True)
def convert_adc_into(b_adc,out_row):

    r1=(b_adc[0]<<24)|(b_adc[1]<<16)|(b_adc[2]<<8)|b_adc[3]
    r2=(b_adc[4]<<24)|(b_adc[5]<<16)|(b_adc[6]<<8)|b_adc[7]
    r3=(b_adc[8]<<24)|(b_adc[9]<<16)|(b_adc[10]<<8)|b_adc[11]
    r4=(b_adc[12]<<24)|(b_adc[13]<<16)|(b_adc[14]<<8)|b_adc[15]

    out_row[0]=(r1>>24)&0b111111
    out_row[1]=(r1>>16)&0b1111111
    out_row[2]=(r1>>2)&0b11111111111111

    out_row[3]=((r1<<12)&0b11000000000000)+((r2>>20)&0b00111111111111)

    out_row[4]=(r2>>6)&0b11111111111111

    out_row[5]=((r2<<8)&0b11111100000000)+((r3>>24)&0b00000011111111)

    out_row[6]=(r3>>10)&0b11111111111111

    out_row[7]=((r3<<4)&0b11111111110000)+((r4>>26)&0b00000000001111)

    out_row[8]=(r4>>14)&0b11111111111111
    out_row[9]=r4&0b11111111111111


@njit(cache=True,fastmath=True)
def decode_events(events_buf,event_sizes,constpeak,temp_a,temp_b,keep_adc):

    n_events=events_buf.shape[0]

    trig=np.zeros(n_events,dtype=np.uint32)
    clk=np.zeros(n_events,dtype=np.uint32)

    hit=np.zeros((n_events,6),dtype=np.uint8)

    temp_raw=np.zeros((n_events,6),dtype=np.uint16)
    temp_corr=np.zeros((n_events,6),dtype=np.float32)

    ped=np.zeros((n_events,6,8),dtype=np.float32)
    wave_sum=np.zeros((n_events,6,8),dtype=np.float32)

    hg=np.zeros((n_events,6),dtype=np.float32)
    lg=np.zeros((n_events,6),dtype=np.float32)

    satu_h=np.zeros((n_events,6),dtype=np.uint8)
    satu_l=np.zeros((n_events,6),dtype=np.uint8)

    x_hg=np.zeros((n_events,6),dtype=np.float32)
    y_hg=np.zeros((n_events,6),dtype=np.float32)

    x_lg=np.zeros((n_events,6),dtype=np.float32)
    y_lg=np.zeros((n_events,6),dtype=np.float32)

    peak_id=np.zeros((n_events,6),dtype=np.int32)

    hitnum=np.zeros(n_events,dtype=np.uint8)

    if keep_adc:
        adc_wave=np.zeros((n_events,6,50,10),dtype=np.float32)
    else:
        adc_wave=np.zeros((0,0,0,0),dtype=np.float32)

    for i in prange(n_events):

        L=event_sizes[i]

        if L<0x28:
            continue

        ev=events_buf[i,:L]

        trig[i]=(ev[4]<<24)|(ev[5]<<16)|(ev[6]<<8)|ev[7]
        clk[i]=(ev[8]<<24)|(ev[9]<<16)|(ev[10]<<8)|ev[11]

        dc_bit=(ev[13]<<8)|ev[14]
        datacnt=(dc_bit>>6)&0b111111

        for j in range(6):
            hit[i,j]=(datacnt>>(5-j))&1

        b14=(ev[0x14]<<24)|(ev[0x15]<<16)|(ev[0x16]<<8)|ev[0x17]
        b58=(ev[0x18]<<24)|(ev[0x19]<<16)|(ev[0x1A]<<8)|ev[0x1B]
        b9C=(ev[0x1C]<<24)|(ev[0x1D]<<16)|(ev[0x1E]<<8)|ev[0x1F]

        raw_t=np.empty(6,dtype=np.uint16)

        raw_t[0]=(b14>>14)&0x3FFF
        raw_t[1]=(b14>>0)&0x3FFF

        raw_t[2]=(b58>>14)&0x3FFF
        raw_t[3]=(b58>>0)&0x3FFF

        raw_t[4]=(b9C>>14)&0x3FFF
        raw_t[5]=(b9C>>0)&0x3FFF

        for ch in range(6):

            temp_raw[i,ch]=raw_t[ch]

            temp_corr[i,ch]=(raw_t[ch]-temp_b[ch])/temp_a[ch]

        adc=np.zeros((6,50,10),dtype=np.float32)

        data_start=0x24
        data_end=L-4

        if data_end>data_start:

            nblocks=(data_end-data_start)//16

            row=np.empty(10,dtype=np.uint16)

            for bi in range(nblocks):

                bs=data_start+bi*16

                convert_adc_into(ev[bs:bs+16],row)

                adc_num=int(row[0])

                ch=-1

                if adc_num==1:
                    ch=0
                elif adc_num==2:
                    ch=1
                elif adc_num==4:
                    ch=2
                elif adc_num==8:
                    ch=3
                elif adc_num==16:
                    ch=4
                elif adc_num==32:
                    ch=5

                if ch!=-1:

                    rid=int(row[1])-1

                    if 0<=rid<50:

                        for k in range(10):
                            adc[ch,rid,k]=float(row[k])

                        for k in range(2,6):

                            if adc[ch,rid,k]>=16383:

                                satu_h[i,ch]=1

                                break

                        for k in range(6,10):

                            if adc[ch,rid,k]>=16383:

                                satu_l[i,ch]=1

                                break

        for ch in range(6):

            for bidx in range(8):

                acc=0.0

                for r in range(PEDESTAL_RANGE):
                    acc+=adc[ch,r,2+bidx]

                ped[i,ch,bidx]=acc/PEDESTAL_RANGE

        if keep_adc:

            for ch in range(6):

                for r in range(50):

                    for k in range(10):
                        adc_wave[i,ch,r,k]=adc[ch,r,k]

        for ch in range(6):

            for r in range(50):

                for bidx in range(8):
                    adc[ch,r,2+bidx]-=ped[i,ch,bidx]

        for ch in range(6):

            best_val=-1e30
            best_r=0

            for r in range(50):

                v=adc[ch,r,6]+adc[ch,r,7]+adc[ch,r,8]+adc[ch,r,9]

                if v>best_val:
                    best_val=v
                    best_r=r

            peak_id[i,ch]=best_r

            if constpeak:

                ws=WAVE_DATA_START
                we=WAVE_DATA_STOP

            else:

                ws=max(0,best_r-WAVE_FRONT)
                we=min(49,best_r+WAVE_BACK)

            for bidx in range(8):

                acc=0.0

                prev=adc[ch,ws,2+bidx]

                for r in range(ws+1,we+1):

                    v=adc[ch,r,2+bidx]

                    acc+=0.5*(prev+v)

                    prev=v

                wave_sum[i,ch,bidx]=acc

            I0,I1,I2,I3=wave_sum[i,ch,0],wave_sum[i,ch,1],wave_sum[i,ch,2],wave_sum[i,ch,3]
            J0,J1,J2,J3=wave_sum[i,ch,4],wave_sum[i,ch,5],wave_sum[i,ch,6],wave_sum[i,ch,7]

            SH=I0+I1+I2+I3
            SL=J0+J1+J2+J3

            hg[i,ch]=SH
            lg[i,ch]=SL

            if SH!=0:
                x_hg[i,ch]=(I0-I1-I2+I3)/SH
                y_hg[i,ch]=(-I0-I1+I2+I3)/SH

            if SL!=0:
                x_lg[i,ch]=(J0-J1-J2+J3)/SL
                y_lg[i,ch]=(-J0-J1+J2+J3)/SL

        for j in range(6):
            hitnum[i]+=hit[i,j]

    return trig,clk,hit,temp_raw,temp_corr,ped,wave_sum,hg,lg,x_hg,y_hg,x_lg,y_lg,peak_id,hitnum,satu_h,satu_l,adc_wave


def process_dat_file(input_path,output_path=None,idx=None,total_files=None,constpeak=False,keep_adc=False,temp_a=None,temp_b=None):

    events_buf,event_sizes,n_events=read_file(input_path)

    if events_buf is None or n_events==0:
        raise ValueError(f"No valid events found in {input_path}")

    if temp_a is None or temp_b is None:
        temp_a=DEFAULT_TEMP_A
        temp_b=DEFAULT_TEMP_B

    trig,clk,hit,temp_raw,temp_corr,ped,wave_sum,hg,lg,x_hg,y_hg,x_lg,y_lg,peak_id,hitnum,satu_h,satu_l,adc_wave = decode_events(events_buf, event_sizes, constpeak, temp_a, temp_b, keep_adc)

    time_arr=(clk.astype(np.int64)-int(clk[0]))/1e6

    wrap=((2**32)-1)/1e6

    time_arr[time_arr<0]+=wrap

    mask=(hit==1)

    ev_idx,ch_idx=np.where(mask)

    trig_sel=trig[ev_idx]
    time_sel=time_arr[ev_idx]

    ch_sel=5-ch_idx

    hg_sel=hg[ev_idx,ch_idx]
    lg_sel=lg[ev_idx,ch_idx]

    x_hg_sel=x_hg[ev_idx,ch_idx]
    y_hg_sel=y_hg[ev_idx,ch_idx]

    x_lg_sel=x_lg[ev_idx,ch_idx]
    y_lg_sel=y_lg[ev_idx,ch_idx]

    ped_sel=ped[ev_idx,ch_idx]
    wave_sel=wave_sum[ev_idx,ch_idx]

    temp_sel=temp_corr[ev_idx]

    peak_id_sel=peak_id[ev_idx,ch_idx]

    satu_h_sel=satu_h[ev_idx,ch_idx]
    satu_l_sel=satu_l[ev_idx,ch_idx]
    
    temp_raw_sel=temp_raw[ev_idx]

    if output_path is None or os.path.isdir(output_path):

        dat_dir=output_path if output_path else os.path.dirname(input_path)

        stem=os.path.splitext(os.path.basename(input_path))[0]

        output_path=os.path.join(dat_dir,f"{stem}.h5")

    output_dir=os.path.dirname(output_path)

    output_stem=os.path.splitext(os.path.basename(output_path))[0]

    adc_output_path=os.path.join(output_dir,f"{output_stem}_adc.h5")

    with h5py.File(output_path,'w') as f_out:

        f_out.create_dataset('trig',data=trig_sel)
        f_out.create_dataset('event_time',data=time_sel)

        f_out.create_dataset('ch',data=ch_sel)

        f_out.create_dataset('hitnum',data=hitnum[ev_idx])

        f_out.create_dataset('wave_sum',data=wave_sel)

        f_out.create_dataset('temp',data=temp_sel)

        f_out.create_dataset('x_hg',data=x_hg_sel)
        f_out.create_dataset('y_hg',data=y_hg_sel)

        f_out.create_dataset('x_lg',data=x_lg_sel)
        f_out.create_dataset('y_lg',data=y_lg_sel)

        f_out.create_dataset('hg',data=hg_sel)
        f_out.create_dataset('lg',data=lg_sel)

        f_out.create_dataset('ped',data=ped_sel)

        f_out.create_dataset('peak_id',data=peak_id_sel)

        f_out.create_dataset('satu_h',data=satu_h_sel)
        f_out.create_dataset('satu_l',data=satu_l_sel)
        
        f_out.create_dataset('temp_raw',data=temp_raw_sel)
        
        f_out.attrs["temp_order"]="TL,TR,ML,MR,BL,BR"

    if keep_adc:

        with h5py.File(adc_output_path,'w') as f_adc:

            f_adc.create_dataset('trig',data=trig)
            f_adc.create_dataset('event_time',data=time_arr)

            f_adc.create_dataset('hit',data=hit)

            f_adc.create_dataset('hitnum',data=hitnum)

            f_adc.create_dataset('temp',data=temp_corr)

            f_adc.create_dataset('ped',data=ped)

            f_adc.create_dataset('peak_id',data=peak_id)

            f_adc.create_dataset('satu_h',data=satu_h)
            f_adc.create_dataset('satu_l',data=satu_l)

            f_adc.create_dataset('adc',data=adc_wave)
            
            f_adc.create_dataset('temp_raw',data=temp_raw)
            
            f_adc.attrs["temp_order"]="TL,TR,ML,MR,BL,BR"
            

        return output_path,adc_output_path

    return output_path


def safe_symlink(target_path,link_path):

    os.makedirs(os.path.dirname(link_path) or ".",exist_ok=True)

    if os.path.lexists(link_path):

        if os.path.islink(link_path):

            try:

                if os.path.realpath(link_path)==os.path.realpath(target_path):
                    return

            except OSError:
                pass

            os.unlink(link_path)

        else:
            raise FileExistsError(f"Cannot create symlink because file exists: {link_path}")

    os.symlink(target_path,link_path)


if __name__=='__main__':

    parser=argparse.ArgumentParser(description='Fast JIT-optimized decoder (pro_decoder_v2): .dat -> .h5')
    parser.add_argument('-i','--input',nargs='+',default=['.'])
    parser.add_argument('-o','--output',default=None)
    parser.add_argument('--retry',action='store_true')
    parser.add_argument('-c','--constpeak',action='store_true')
    parser.add_argument('--adc',action='store_true')
    parser.add_argument('--link-dir',default=None)

    args=parser.parse_args()

    run_config=load_run_config(RUN_FILE)
    psa_config=load_psa_config(run_config)

    target_files=[]

    for path in args.input:

        if os.path.isfile(path):

            if path.lower().endswith('.dat'):

                stem=os.path.basename(path)

                sn=stem.split('_')[2]

                target_files.append({
                    "sn":sn,
                    "path":path
                })

            continue

        if os.path.isdir(path):

            entries=[os.path.join(path,x) for x in os.listdir(path)]

            direct_dat=[
                p for p in entries
                if os.path.isfile(p)
                and p.lower().endswith('.dat')
            ]

            if len(direct_dat)>0:

                for p in direct_dat:

                    stem=os.path.basename(p)

                    sn=stem.split('_')[2]

                    target_files.append({
                        "sn":sn,
                        "path":p
                    })

                continue

            for subdir in [p for p in entries if os.path.isdir(p)]:

                sn=os.path.basename(subdir)

                for fname in os.listdir(subdir):

                    if not fname.lower().endswith('.dat'):
                        continue

                    target_files.append({
                        "sn":sn,
                        "path":os.path.join(subdir,fname)
                    })

            continue

        matches=glob.glob(path,recursive=True)

        if len(matches)==0:

            print(f"[WARN] Input not found: {path}")

            continue

        for p in matches:

            if os.path.isfile(p) and p.lower().endswith('.dat'):

                stem=os.path.basename(p)

                sn=stem.split('_')[2]

                target_files.append({
                    "sn":sn,
                    "path":p
                })

    if len(target_files)==0:

        print("[ERROR] No target .dat files found.")

        sys.exit(1)

    target_files=natsorted(
        target_files,
        key=lambda x:(x["sn"],x["path"])
    )

    total_files=len(target_files)

    print(f"Total target files: {total_files}")

    os.makedirs(DEFAULT_H5_STORE_DIR,exist_ok=True)

    for idx,target in enumerate(target_files,start=1):

        sn=target["sn"]

        path=target["path"]

        unit=get_unit_by_sn(sn,psa_config)

        if unit is None:

            print(f"[WARN] calibration not found: {sn}")

            temp_a,temp_b=DEFAULT_TEMP_A,DEFAULT_TEMP_B

        else:

            temp_a,temp_b=get_temp_calibration(unit)

        stem=os.path.splitext(os.path.basename(path))[0]

        link_dir=args.link_dir if args.link_dir else (os.path.dirname(path) or ".")

        os.makedirs(link_dir,exist_ok=True)

        h5_path=os.path.join(DEFAULT_H5_STORE_DIR,f"{stem}.h5")

        adc_h5_path=os.path.splitext(h5_path)[0]+'_adc.h5'

        h5_link_path=os.path.join(link_dir,f"{stem}.h5")

        adc_h5_link_path=os.path.join(link_dir,f"{stem}_adc.h5")

        if (not args.retry) and os.path.exists(h5_path) and ((not args.adc) or os.path.exists(adc_h5_path)):

            safe_symlink(h5_path,h5_link_path)

            if args.adc:
                safe_symlink(adc_h5_path,adc_h5_link_path)

            continue

        try:

            out=process_dat_file(
                path,
                h5_path,
                idx,
                total_files,
                constpeak=args.constpeak,
                keep_adc=args.adc,
                temp_a=temp_a,
                temp_b=temp_b
            )

            if args.adc:

                h5_out,adc_out=out

                safe_symlink(h5_out,h5_link_path)

                safe_symlink(adc_out,adc_h5_link_path)

                print(f"[{idx}/{total_files}] Finished: {sn} : {h5_out}, {adc_out}")

            else:

                safe_symlink(out,h5_link_path)

                print(f"[{idx}/{total_files}] Finished: {sn} : {out}")

        except Exception as e:

            print(f"[ERROR] Failed to process {path}: {e}",file=sys.stderr)

    print('~~task is over~~')