# SMILE-3 QL

Storage Node（calltt, killua）が書いている run から TPC と HA（PSA）の QL を出す。全データを
QL 用 PC へ同期するコマンドは別にある。

```
calltt / killua                                QL 用 PC
  ql_feed.py（間引いて流す）──ssh、magallan を ProxyJump──▶ ql_stream.py
                                                          │ localhost:10012
                                                          ▼
                                        ql_tpc.py（:8090） / ql_psa.py（:8091）
  /data/smile_ics/storage_gX ──ql_sync.sh（全データ）──▶ ~/smile_ics_ql/storage_gX/<run>/
```

QL はノード側で間引いた分（既定 450 分の 1）だけを受けるので、数秒の遅れで更新される。
全データの同期とは独立に動く。

## 使い方

```bash
cd ~/analysis_sm3_FM
./ql_live.sh                    # QL を起動する。Ctrl-C で全部止まる
./ql_live.sh --from-start       # 今の run を最初から積算する
```

全データの同期は、別のターミナルで実行する。

```bash
./ql_sync.sh --watch                  # 30 秒ごとに同期し続ける
./ql_sync.sh --watch --since 20261101 # この日付より前の run は同期しない
./ql_sync.sh                          # 1 回だけ同期する
```

同期先は `QL_DEST`（既定 `~/smile_ics_ql`）で変える。全データの同期はノードと回線への負荷が
大きいので、RUN 中に要らなければ RUN の合間に 1 回だけ実行する。

起動したときに走っている run は、起動した時点から積算する。そのあとに始まった run は
最初から積算する。

ブラウザで TPC は http://127.0.0.1:8090 、HA は http://127.0.0.1:8091 を開く。TPC の FADC 波形と
その積算は http://127.0.0.1:8090/fadc 、HA の MPPC の温度は http://127.0.0.1:8091/temp に出る。
起動するとこの 4 つのページが開く。温度は `PSA_analysis/param/psa_param_20261029.csv` の HA ごとの係数で
較正する（係数がほかの HA の写しの HA には `*`、psa_param に無い HA には `?` を付ける）。
図の PNG は `check_output/` に書く。

## ファイル

| ファイル | 役割 |
| --- | --- |
| `ql_live.sh` | `ql_stream.py`、`ql_tpc.py`、`ql_psa.py` をまとめて起動する |
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

## 解析のコード（QL 以外）

QL のほかに、HA（PSA）のデータを後から解析するスクリプトも `PSA_analysis/src_onoda/` に置いている。
詳しい使い方と結果の見方は Notion の「analysis_sm3_FM のコードの解説」にまとめている。

```
analysis_sm3_FM/
├── ql_live.sh, ql_sync.sh          QL の起動と全データの同期
├── TPC_analysis/                   TPC の QL と確認図（check_tpc.py）
├── PSA_analysis/
│   ├── src_onoda/                  HA の解析と QL（ここのコードを使う）
│   │   └── sub_scr/                デコードの下回り、パラメータの読み込み、SiTCP の通信
│   ├── param/                      psa_param*.csv（HA の SN・HV・温度の較正）、run_file.json
│   ├── original/                   最初にもらった元のコード（比べるために残す。使わない）
│   └── （以下は git に入れない）h5_runs*/ など .h5、peak_*/ など結果、logs/
└── check_output/                   確認図（git に入れない）
```

HA の解析の流れ（Storage Node の run を .h5 にしてから解析する）:

| 順 | スクリプト | すること |
| --- | --- | --- |
| 1 | `decode_runs.py` | 複数の run の `.dat` を HA ごとの `.h5` にデコードする（`psa_decode.py` を使う） |
| 2 | `check_data.py` | バタフライマップ、hg、peak_id、ペデスタルの確認図 |
| 3 | `temp_spectra.py` | 温度で区切った hg スペクトルと光電ピークの位置 |
| 4 | `peak_search.py` | バタフライマップの 64 個のピークを探す |
| 5 | `grid_spect.py` | 位置補正とピクセルごとのスペクトル |
| 6 | `pixel_id.py` | 位置補正を新しいデータに当てはめ、ピクセル番号を付ける |
| 7 | `energy_calib.py` | Cs-137 と Co-60 のピークから ADC → エネルギーの換算係数を出す |
| - | `check_waveform.py` | `.dat` の生波形を重ねた図 |
| - | `check_ha_order.py` | Storage Node で HA の `.dat` の trig の戻りを調べる |

```bash
python3 PSA_analysis/src_onoda/decode_runs.py -i "/mnt/d/UVSOR2026_late/*/20261103_run12/part*" -o PSA_analysis/h5_runs_run12
python3 PSA_analysis/src_onoda/check_data.py -i PSA_analysis/h5_runs_run12/ha02 --tag run12_ha02
python3 PSA_analysis/src_onoda/energy_calib.py --ha ha02
```

既定の入力・出力先は、どれもこのリポジトリの中からの相対で決まる。画面に図を出さずに保存だけするときは
`MPLBACKEND=Agg` を付ける（XLaunch が無いと peak_search などが止まる）。
