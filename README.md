# SMILE-3 QL

Storage Node（calltt, killua）が書いている run から TPC と HA（PSA）の QL を出し、あわせて
全データを QL 用 PC へ同期する。

```
calltt / killua                                QL 用 PC
  ql_feed.py（間引いて流す）──ssh、magallan を ProxyJump──▶ ql_stream.py
                                                          │ localhost:10012
                                                          ▼
                                        ql_tpc.py（:8090） / ql_psa.py（:8091）
  /data/smile_ics/storage_gX ──ql_sync.sh（全データ）──▶ ~/smile_ics_ql/storage_gX/<run>/
```

QL はノード側で間引いた分（既定 450 分の 1）だけを受けるので、全データの同期を待たず、
数秒の遅れで更新される。全データの同期は裏で続く。

## 使い方

```bash
cd ~/analysis_sm3_FM
./ql_live.sh                    # QL と全データの同期をまとめて起動する。Ctrl-C で全部止まる
./ql_live.sh --since 20261101   # この日付より前の run は同期しない
QL_SYNC=0 ./ql_live.sh          # QL だけ出し、全データの同期はしない
QL_FOLLOW_ARGS=--from-start ./ql_live.sh   # 今の run を最初から積算する
./ql_sync.sh --watch            # 同期だけ
```

起動したときに走っている run は、起動した時点から積算する。そのあとに始まった run は
最初から積算する。

ブラウザで TPC は http://127.0.0.1:8090 、HA は http://127.0.0.1:8091 を開く。
図の PNG は `check_output/` に書く。

## ファイル

| ファイル | 役割 |
| --- | --- |
| `ql_live.sh` | 下をまとめて起動する |
| `TPC_analysis/ql_stream.py` | 各ノードで `ql_feed.py` を ssh で動かし、間引いたイベントを合わせて流す |
| `TPC_analysis/ql_feed.py` | ノード側で run のファイルを追いかけて間引く。ssh で渡すのでノードに置かなくてよい |
| `TPC_analysis/ql_tpc.py` | TPC の QL |
| `PSA_analysis/src_onoda/ql_psa.py` | HA の QL |
| `ql_sync.sh` | Storage Node の全データを QL 用 PC へ同期する |
| `TPC_analysis/ql_follow.py` | 同期したファイルから間引いたイベントを流す（同期済みのデータで QL を出し直すとき） |
| `TPC_analysis/ql_replay.py` | 保存済みの run を流す試験用の送信器 |

残りは QL が使う解析の関数。

## QL 用 PC の準備

1. このリポジトリを `~/analysis_sm3_FM` に clone する。`psa_decode.py` がこのディレクトリ名で
   モジュールを読むので、名前は変えない。

   ```bash
   git clone git@github.com:haru2510/QL.git ~/analysis_sm3_FM
   ```

   private リポジトリなので、QL 用 PC の公開鍵をリポジトリの Deploy keys（読み取りのみ）に
   登録しておく。

2. Python のパッケージと rsync を入れる。Storage Node には python3 と rsync が要る
   （`ql_feed.py` は標準ライブラリだけで動く）。

   ```bash
   pip install numpy matplotlib natsort numba h5py
   ```

3. `~/.ssh/config` に magallan と Storage Node を書き、公開鍵を 3 台に登録する
   （`ssh-copy-id magallan` などで登録できる）。

   ```
   Host magallan
     HostName 192.168.1.xx
     User msgc
   Host calltt killua
     User msgc
     ProxyJump magallan
   Host calltt
     HostName 192.168.100.12
   Host killua
     HostName 192.168.100.13
   ```

## 更新

```bash
cd ~/analysis_sm3_FM && git pull
```
