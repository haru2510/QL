# SMILE-3 QL

Storage Node（calltt, killua）が保存した run を QL 用 PC へ同期し、TPC と HA（PSA）の QL を出す。

```
calltt / killua ──(ssh、magallan を ProxyJump)──▶ QL 用 PC ~/smile_ics_ql/storage_gX/<run>/
                                                    │ ql_follow.py が新しいファイルを間引いて流す
                                                    ▼ localhost:10012
                                       ql_tpc.py（:8090） / ql_psa.py（:8091）
```

## 使い方

```bash
cd ~/analysis_sm3_FM
./ql_live.sh                    # 同期と QL をまとめて起動する。Ctrl-C で全部止まる
./ql_live.sh --since 20261101   # この日付より前の run は同期しない
./ql_sync.sh --watch            # 同期だけ
QL_FOLLOW_ARGS=--from-start ./ql_live.sh   # 今の run を最初から積算する
```

起動したときに走っている run は、それまでに届いている分を読み飛ばし、起動した時点から
積算する（長い run の途中で起動しても、すぐに表示が始まる）。そのあとに始まった run は
最初から積算する。

ブラウザで TPC は http://127.0.0.1:8090 、HA は http://127.0.0.1:8091 を開く。
図の PNG は `check_output/` に書く。

書き込み中のファイルも、書き足された分ずつ同期して読むので、QL は同期の間隔（既定 30 秒）
程度の遅れで更新される。

## ファイル

| ファイル | 役割 |
| --- | --- |
| `ql_live.sh` | 下の 4 つをまとめて起動する |
| `ql_sync.sh` | Storage Node の run を QL 用 PC へ同期する |
| `TPC_analysis/ql_follow.py` | 同期したファイルから間引いたイベントを流す |
| `TPC_analysis/ql_tpc.py` | TPC の QL |
| `PSA_analysis/src_onoda/ql_psa.py` | HA の QL |
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

2. Python のパッケージと rsync を入れる。

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
