#!/usr/bin/env python3

import numba
import numpy as np

from numba import njit,prange

numba.set_num_threads(4)

HEADER=b'\xeb\x90\x19\x64'
FOOTER=b'\x75\x50\x49\x43'

PEDESTAL_RANGE=23

WAVE_FRONT=4
WAVE_BACK =9

WAVE_DATA_START=25
WAVE_DATA_STOP =41


@njit(cache=True)
def find_events_jit(buf,H0,H1,H2,H3,F0,F1,F2,F3):

    n=buf.size

    starts=[]
    ends  =[]

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

    events_buf =np.zeros((n,max_len),dtype=np.uint8)
    event_sizes=np.zeros(n,dtype=np.int32)

    for i in range(n):

        L=lengths[i]

        event_sizes[i]=L
        events_buf[i,:L]=buf[starts[i]:ends[i]]

    return events_buf,event_sizes,n


@njit(cache=True,fastmath=True,nogil=True)
def convert_adc_into(b_adc,out_row):

    r1=(b_adc[0] <<24)|(b_adc[1] <<16)|(b_adc[2] <<8)|b_adc[3]
    r2=(b_adc[4] <<24)|(b_adc[5] <<16)|(b_adc[6] <<8)|b_adc[7]
    r3=(b_adc[8] <<24)|(b_adc[9] <<16)|(b_adc[10]<<8)|b_adc[11]
    r4=(b_adc[12]<<24)|(b_adc[13]<<16)|(b_adc[14]<<8)|b_adc[15]

    out_row[0]=(r1>>24)&0b111111
    out_row[1]=(r1>>16)&0b1111111
    out_row[2]=(r1>> 2)&0b11111111111111

    out_row[3]=((r1<<12)&0b11000000000000)+((r2>>20)&0b00111111111111)

    out_row[4]=(r2>> 6)&0b11111111111111

    out_row[5]=((r2<<8)&0b11111100000000)+((r3>>24)&0b00000011111111)

    out_row[6]=(r3>>10)&0b11111111111111

    out_row[7]=((r3<<4)&0b11111111110000)+((r4>>26)&0b00000000001111)

    out_row[8]=(r4>>14)&0b11111111111111
    out_row[9]= r4     &0b11111111111111


@njit(cache=True,fastmath=True)
def decode_events(events_buf,event_sizes,constpeak,temp_a,temp_b,keep_adc):

    n_events=events_buf.shape[0]

    trig=np.zeros(n_events,dtype=np.uint32)
    clk =np.zeros(n_events,dtype=np.uint32)

    hit=np.zeros((n_events,6),dtype=np.uint8)

    temp_raw =np.zeros((n_events,6),dtype=np.uint16)
    temp_corr=np.zeros((n_events,6),dtype=np.float32)

    ped     =np.zeros((n_events,6,8),dtype=np.float32)
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
        clk[i] =(ev[8]<<24)|(ev[9]<<16)|(ev[10]<<8)|ev[11]

        dc_bit =(ev[13]<<8)|ev[14]
        datacnt=(dc_bit>>6)&0b111111

        for j in range(6):
            hit[i,j]=(datacnt>>(5-j))&1

        b14=(ev[0x14]<<24)|(ev[0x15]<<16)|(ev[0x16]<<8)|ev[0x17]
        b58=(ev[0x18]<<24)|(ev[0x19]<<16)|(ev[0x1A]<<8)|ev[0x1B]
        b9C=(ev[0x1C]<<24)|(ev[0x1D]<<16)|(ev[0x1E]<<8)|ev[0x1F]

        raw_t=np.empty(6,dtype=np.uint16)

        raw_t[0]=(b14>>14)&0x3FFF
        raw_t[1]=(b14>> 0)&0x3FFF

        raw_t[2]=(b58>>14)&0x3FFF
        raw_t[3]=(b58>> 0)&0x3FFF

        raw_t[4]=(b9C>>14)&0x3FFF
        raw_t[5]=(b9C>> 0)&0x3FFF

        for ch in range(6):

            temp_raw[i,ch]=raw_t[ch]

            temp_corr[i,ch]=(raw_t[ch]-temp_b[ch])/temp_a[ch]

        adc=np.zeros((6,50,10),dtype=np.float32)

        data_start=0x24
        data_end  =L-4

        if data_end>data_start:

            nblocks=(data_end-data_start)//16

            row=np.empty(10,dtype=np.uint16)

            for bi in range(nblocks):

                bs=data_start+bi*16

                convert_adc_into(ev[bs:bs+16],row)

                adc_num=int(row[0])

                ch=-1

                if   adc_num==1 : ch=0
                elif adc_num==2 : ch=1
                elif adc_num==4 : ch=2
                elif adc_num==8 : ch=3
                elif adc_num==16: ch=4
                elif adc_num==32: ch=5

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
            best_r  =0

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

                ws=max(0 ,best_r-WAVE_FRONT)
                we=min(49,best_r+WAVE_BACK)

            for bidx in range(8):

                acc =0.0
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

                x_hg[i,ch]=( I0-I1-I2+I3)/SH
                y_hg[i,ch]=(-I0-I1+I2+I3)/SH

            if SL!=0:

                x_lg[i,ch]=( J0-J1-J2+J3)/SL
                y_lg[i,ch]=(-J0-J1+J2+J3)/SL

        for j in range(6):
            hitnum[i]+=hit[i,j]

    return (
        trig,clk,hit,
        temp_raw,temp_corr,
        ped,wave_sum,
        hg,lg,
        x_hg,y_hg,
        x_lg,y_lg,
        peak_id,
        hitnum,
        satu_h,satu_l,
        adc_wave
    )


def decode_dat_file(file_path,constpeak,temp_a,temp_b,keep_adc=False):

    events_buf,event_sizes,n_events=read_file(file_path)

    if events_buf is None or n_events==0:
        raise ValueError(f"No valid events found in {file_path}")

    out=decode_events(
        events_buf,
        event_sizes,
        constpeak,
        temp_a,
        temp_b,
        keep_adc
    )

    return {
        "trig"     :out[0],
        "clk"      :out[1],
        "hit"      :out[2],
        "temp_raw" :out[3],
        "temp_corr":out[4],
        "ped"      :out[5],
        "wave_sum" :out[6],
        "hg"       :out[7],
        "lg"       :out[8],
        "x_hg"     :out[9],
        "y_hg"     :out[10],
        "x_lg"     :out[11],
        "y_lg"     :out[12],
        "peak_id"  :out[13],
        "hitnum"   :out[14],
        "satu_h"   :out[15],
        "satu_l"   :out[16],
        "adc_wave" :out[17]
    }