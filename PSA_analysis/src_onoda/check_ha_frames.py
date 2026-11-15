#!/usr/bin/env python3
#
# Storage Node で、最新 run の HA の .dat (書き込み中の .tmp も) を読み、枠の長さ (ヘッダから次のヘッダ
# までのバイト数) と trig の範囲を HA ごとに出す。QL で「no data」のままの HA があるときに、
# ql_feed.py が枠と認める形 (40 + 800 x ヒットした ch の数のバイトで、フッタで終わる) かを確かめる。
# ノードに置かなくてよい。QL 用 PC から ssh の標準入力で渡す。標準ライブラリだけで動く。
#
#   ssh killua python3 - /data/smile_ics/storage_g1 0 1 < PSA_analysis/src_onoda/check_ha_frames.py
#     引数: storage のディレクトリ, HA 番号 (複数可), [--run run ID (省略時は最新)]
#
import os, re, sys
from collections import Counter
H, F = bytes.fromhex("eb901964"), bytes.fromhex("75504943")
THIN = 450
argv = sys.argv[1:]
run = None
if "--run" in argv:
    i = argv.index("--run"); run = argv[i + 1]; del argv[i:i + 2]
root, has = argv[0], [int(a) for a in argv[1:]]
runs = sorted(d for d in os.listdir(root) if re.match(r"\d{8}_run\d+$", d))
run = run or runs[-1]
d = os.path.join(root, run)
names = os.listdir(d)
print(f"{d}")
for ha in has:
    pat = re.compile(rf"ha{ha:02d}_(\d+)\.dat(\.tmp)?$")
    files = sorted((int(m.group(1)), n) for n in names for m in [pat.match(n)] if m)
    others = sorted(n for n in names if n.startswith(f"ha{ha:02d}") and not pat.match(n))[:5]
    if not files:
        print(f"ha{ha:02d}: ql_feed.py が拾う名前のファイルが無い。似た名前: {others}")
        continue
    lengths, ok, n_frames, thin, trigs = Counter(), 0, 0, 0, []
    for seq, name in files[-3:]:
        data = open(os.path.join(d, name), "rb").read()
        pos = data.find(H)
        while 0 <= pos < len(data):
            nxt = data.find(H, pos + 4)
            end = nxt if nxt >= 0 else len(data)
            lengths[end - pos] += 1
            n_frames += 1
            trig = int.from_bytes(data[pos + 4:pos + 8], "big"); trigs.append(trig)
            if any(data[pos + s - 4:pos + s] == F for s in range(40, 4841, 800) if pos + s <= len(data)):
                ok += 1
                thin += trig % THIN == 0
            pos = nxt
    print(f"ha{ha:02d}: {len(files)} files ({files[0][1]} .. {files[-1][1]}), 最後の 3 ファイルの枠 {n_frames}")
    print(f"  枠の長さ [バイト: 個数] {dict(lengths.most_common(5))}")
    print(f"  ql_feed.py が枠と認める {ok}, そのうち trig が {THIN} の倍数で QL へ送る {thin}")
    if trigs:
        print(f"  trig {trigs[0]} .. {trigs[-1]} (最小 {min(trigs)}, 最大 {max(trigs)})")
