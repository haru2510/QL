#!/usr/bin/env python3
#
# Storage Node で、HA の閉じた .dat を名前の順に読み、trig が 1000 以上戻った箇所を出す。
# QL で「trig went back」が出たときに、データそのものに戻りがあるかを確かめる。
# ノードに置かなくてよい。QL 用 PC から ssh の標準入力で渡す。標準ライブラリだけで動く。
#
#   ssh calltt python3 - /data/smile_ics/storage_g0 6 100 < PSA_analysis/src_onoda/check_ha_order.py
#   ssh calltt python3 - /data/smile_ics/storage_g0 15 9550-9556 < PSA_analysis/src_onoda/check_ha_order.py
#     引数: storage のディレクトリ, HA 番号, 最後から何ファイル読むか (または番号の範囲 9550-9556),
#           [run ID (省略時は最新)]
#
import os, re, sys
root, ha, which = sys.argv[1], int(sys.argv[2]), sys.argv[3]
H, F = bytes.fromhex("eb901964"), bytes.fromhex("75504943")
runs = sorted(d for d in os.listdir(root) if re.match(r"\d{8}_run\d+$", d))
run = sys.argv[4] if len(sys.argv) > 4 else runs[-1]
d = os.path.join(root, run)
pat = re.compile(rf"ha{ha:02d}_(\d+)\.dat$")
files = sorted((int(m.group(1)), n) for n in os.listdir(d) for m in [pat.match(n)] if m)
if "-" in which:
    first_seq, last_seq = map(int, which.split("-"))
    files = [(seq, n) for seq, n in files if first_seq <= seq <= last_seq]
else:
    files = files[-int(which):]
print(f"{run} ha{ha:02d}: {len(files)} files ({files[0][1]} .. {files[-1][1]})")
prev = None; n = 0; drops = 0
for seq, name in files:
    data = open(os.path.join(d, name), "rb").read(); pos = 0; first = None
    while pos + 40 <= len(data):
        if data[pos:pos + 4] != H:
            nxt = data.find(H, pos + 1)
            if nxt < 0: break
            pos = nxt; continue
        size = 40 if data[pos + 36:pos + 40] == F else 840
        trig = int.from_bytes(data[pos + 4:pos + 8], "big"); n += 1
        first = trig if first is None else first
        if prev is not None and trig + 1000 < prev:
            drops += 1
            if drops <= 10: print(f"  戻り: {name} の {pos} バイト目 {prev} -> {trig} (差 {prev - trig})")
        prev = trig; pos += size
    if "-" in which or seq in (files[0][0], files[-1][0]):
        print(f"  {name}: {len(data)} バイト, trig {first} .. {prev}")
print(f"枠 {n}, 1000 以上の戻り {drops} 回")
