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

`ql_sync.sh` は zeno の `/opt/smile-ics/data/scaler/tcu_scaler*.jsonl` も
`$QL_DEST/scaler/` へ同期する。生カウンター（`dead_time_10mhz`、`real_time_10mhz`）と
差分レート（`dead_time_ratio`）の最新確認用・時次ファイルを含む。`--since` は run のみに適用し、
Scaler は全期間を同期する。取得元は `QL_SCALER_NODE` と `QL_SCALER_DIR` で変更でき、
`QL_SCALER_NODE=none ./ql_sync.sh` で Scaler の同期を無効にできる。

起動したときに走っている run は、起動した時点から積算する。そのあとに始まった run は
最初から積算する。

ブラウザで TPC は http://127.0.0.1:8090 、HA は http://127.0.0.1:8091 を開く。TPC の FADC 波形と
その積算は http://127.0.0.1:8090/fadc 、HA の MPPC の温度は http://127.0.0.1:8091/temp に出る。
起動するとこの 4 つのページが開く（前に開いたタブが残っていれば、新しく開かずにそのタブで表示を再開する。QL が止まっている間は、ページの右下にそう出る）。温度は `PSA_analysis/param/psa_param_20261029.csv` の HA ごとの係数で
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
| `ETCC_analysis/` | ビーム試験の解析（偏光、TPC の頂点とビームの軸、ビームの周期）。「ETCC のビーム試験の解析」を参照。QL では使わない |

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

3. `~/.ssh/config` に magallan、Storage Node、zeno を書き、公開鍵を 4 台に登録する
   （`ssh-copy-id magallan` などで登録できる）。

   ```
   Host magallan
     HostName 192.168.1.xx
     User msgc
   Host calltt killua zeno
     User msgc
     ProxyJump magallan
   Host zeno
     HostName 192.168.100.11
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

## ETCC のビーム試験の解析（`ETCC_analysis/`）

HA と TPC の両方を使うビーム試験の解析は `ETCC_analysis/` に置く。入力と条件の JSON は `ETCC_analysis/param/`、出力は既定で `ETCC_analysis/output/`。
どれも raw を直接読み、ランは `"20261201_run05"` と書けば `<nas>/20261201/{calltt,killua}/20261201_run05` を読む（`{"date": ..., "run": ...}` や、ディレクトリを直接 `{"name": ..., "dirs": [...]}` でもよい）。
コマンドラインで渡したものは JSON より優先する。

| スクリプト | JSON | すること |
| --- | --- | --- |
| `pair_polarization.py` | `pair_polarization.json` | HA のヒットが 2 つのイベントで、2 点を結んだ線のビーム軸まわりの方位角から偏光の変調と向きを fit する |
| `beam_vertex.py` | `beam_vertex.json` | TPC の散乱点（対生成の頂点）の 3 次元の分布と、ビームの軸 |
| `beam_period.py` | `beam_period.json` | イベントの時刻の周期から、ビームのパルスの周期と位相（on-beam の窓）を出す |

### 偏光ビーム試験の簡易解析（`pair_polarization.py`）

対生成のイベントで、HA のヒットが 2 つ（`n_hits`）のものを選び、2 点を結んだ線をビームに垂直な面に射影した方位角 φ（0〜180°）の分布を `A [1 + μ cos 2(φ − φ0)]` で fit する。
`.h5` を作らずに raw（`ha<NN>_<番号>.dat`）を直接読む。1 つのランのボードは 2 つの storage（calltt、killua）に分かれるので両方読み、ラン 1 つの中でボード同士をトリガーカウンタ（`trig`）でそろえる（ボードの中で `trig` が戻ったら epoch を 1 つ進め、(epoch, `trig`) で組にする。EricaNext の `PSA_ReadData` と同じ）。
複数のランは、ランごとに組にしてから方位角を足す（ランをまたいで `trig` を混ぜない）。

**入力と条件は `ETCC_analysis/param/pair_polarization.json` に書く**（既定でこれを読む）。ラン（`runs`）、比べるラン（`reference_runs`）、ビームの向き、カットなど。各キーの説明はファイルの `_doc`。
今はテスト用に手元の ha03 のラン（ペアは無い）が入っている。
コマンドラインで渡したものは JSON より優先する。`--runs` や `-i` を渡すと、JSON の `runs` と `reference_runs` は使わない。

実行のしかた（WSL で。`/usr/bin/python3` を使う。numba、h5py、matplotlib が入っている方）:

```bash
cd ~/analysis_sm3_FM

# 1. ETCC_analysis/param/pair_polarization.json のとおりに流す（出力は既定で ETCC_analysis/output/pair_polarization.*）
/usr/bin/python3 ETCC_analysis/pair_polarization.py

# 2. 別の設定ファイル（ビーム試験用に写して書き換えたものなど）
cp ETCC_analysis/param/pair_polarization.json ETCC_analysis/param/pair_beamtest.json   # runs などを書き換える
/usr/bin/python3 ETCC_analysis/pair_polarization.py --config ETCC_analysis/param/pair_beamtest.json

# 3. ランをコマンドラインで（複数なら足す）。まず --max-files 2 で数ファイルだけ試すとよい
/usr/bin/python3 ETCC_analysis/pair_polarization.py \
  --runs 20261201_run05 20261201_run07 --reference-runs 20261201_run06 \
  --beam-dir 0 0 1 --max-files 2 -o ETCC_analysis/output/run05_07_vs_06_test

# 4. ディレクトリを直接（1 ラン。NAS 以外の場所にあるとき）
/usr/bin/python3 ETCC_analysis/pair_polarization.py \
  -i /path/to/calltt/20261201_run05 /path/to/killua/20261201_run05 -o ETCC_analysis/output/run05
```

- 日付・ラン番号（`20261201`、`run05` など）は例。実際のランに置き換える。
- ヒットの位置はブロック（8×8 の結晶）の中心で、`PSA_analysis/param/geometry_sm3.csv`（EricaNext の `calib/psa/geometry.csv` の写し）から取る。e± は HG が飽和しやすく、飽和すると重心がブロックの中央に寄るので、ピクセルまでは求めない。
- ch は ADC 番号（1, 2, 4, 8, 16, 32 → 0〜5、EricaNext と同じ）。`psa_decode.py` の `.h5` の ch（5 − ch）とは逆なので注意。
- 座標は検出器の座標（X: 東、Y: 南、Z: 下）。`beam_dir` はビームの向き（上から下へなら `[0, 0, 1]`、東から西へなら `[-1, 0, 0]`）、`ref_axis` は φ = 0 の軸（ビームに垂直な面に射影して使う）。試験の配置に合わせる。
- ヒットは HG ≥ `hg_min`（既定 1000 ADC）で、peak_id がボードごとの最頻値 ± `peak_id_window`（既定 3）のもの。
- 検出器の形（HA の並び）だけでも φ の分布は平らにならないので、偏光の大きさは比べるランとの比（`N / N_ref`）の fit で見る。偏光の向きを 90° 回したランと比べたときは、比の μ は偏光の効果の 2 倍（2 で割って読む）。偏光なしのランならそのまま。
- 出力: `<out>.png`（全ランを足した φ の分布と fit、比、選んだペアのブロック、1 トリガーのヒットの数）、`<out>_pairs.csv`（ラン、epoch、trig ごとのペアの一覧）、`<out>.json`（全体とランごとの fit、条件、ボードごとのフレーム数と trig の戻った回数）。

### ビーム試験: TPC の散乱点の 3 次元の分布からビームの道を出す（`beam_vertex.py`）

γ線ビームが TPC のガスの中で相互作用した点（対生成の頂点）を事象ごとに出し、その 3 次元の分布と、頂点に fit した直線（ビームの軸）を描く。
TPC の raw（`tpc_enc_a*`、`tpc_enc_c*`）を直接読む（読み方は `TPC_analysis/check_tpc.py` と同じ）。ランの書き方と JSON の扱いは `pair_polarization.py` と同じで、既定で `ETCC_analysis/param/beam_vertex.json` を読む。

- **頂点**: 同じ uPIC のアノードとカソードの両方にヒットがある一番早いクロック（z が一番小さい点）。x、y はそのクロックのストリップの位置の平均。対生成の e± は頂点から前に飛ぶので、ビームが上から下（+Z）なら一番上流の点になる。ToT 補正と 3 次元の組み合わせ（偽の点が出る）はしない。
  EricaNext の `TOT_Skewness` の散乱点（歪度で飛跡の端を選ぶ）は電子 1 本の前提なので、対生成には使わない。
- **z** = `velocity` × (クロック − `clock_offset`) + `z_gem`（EricaNext と同じ。既定値は EricaNext の `calib/tpc/drift_velocity.dat`）。**z の絶対値は、シンチ（HA）などの外のトリガーで時刻の原点が決まっているランでしか意味が無い**。TPC 自身のトリガーのラン（TPC_cal など）では、一番早いヒットがいつもトリガーの時刻に来るので、頂点が z の端に張り付く。
- **fiducial**: 頂点が壁（±153.6 mm）とドリフト領域の両端（`z_gem`、317 mm）から `margin`（既定 10 mm）だけ内側。外（容器、天板、HA）で対生成して入ってきた e± の事象を落とす。
- **ノイズのカットはしない**（鳴りっぱなしのストリップもヒットの少ない事象も使う）。事象のヒットの数は `_vertices.csv` の `n_hits` に書くので、要れば後で切る。
- **ビームの軸**: fiducial を通った頂点に 3 次元の直線を fit する（主成分。軸からの距離が `fit_cut_sigma` × RMS を超える頂点を外して fit し直す）。Z からの角度、方位角、軸からの距離の分布（R50、R68、R90）を出す。
- メモリを抑えるため、フレームごとにクロックでまとめ（ヒットの数と位置の和）、早いクロックから `keep_clocks`（既定 16）個だけ持つ。

```bash
cd ~/analysis_sm3_FM
# JSON のとおり（今は 20261116 run03、ETCC のモード。NAS にある）
/usr/bin/python3 ETCC_analysis/beam_vertex.py
# ビーム試験のラン（まず --max-files 2 で試す。複数なら足す）
/usr/bin/python3 ETCC_analysis/beam_vertex.py --runs 20261123_run03 20261123_run04 --max-files 2 \
  -o ETCC_analysis/output/run03_04_test
# ディレクトリを直接
/usr/bin/python3 ETCC_analysis/beam_vertex.py -i /path/to/calltt/20261123_run03 /path/to/killua/20261123_run03 \
  -o ETCC_analysis/output/run03
```

出力: `<out>.png`（全頂点と fiducial の頂点の X–Y、X–Z、Y–Z の分布、fiducial の箱（黒の破線）、fit した軸（マゼンタ））、`<out>_vertices.csv`（事象ごとの頂点、fiducial かどうか、軸からの距離）、`<out>.json`（軸の fit、条件、ボードごとのフレーム数）。

### ビームのタイミングに同期をとる: イベントの時刻の周期（`beam_period.py`）

パルスのビーム（レーザー逆コンプトンなど）に同期をとるため、イベントの時刻の周期と位相を出す。時刻は各ボードのカウンタ:
HA は `ha<NN>_*.dat` の clk（1 MHz）、TPC は `tpc_enc_*.raw` の N_time（100 us）。どちらも 32 bit で、2^32 で戻ったら足し直す。
**ボードごとにカウンタの原点が違うので、ボードごと・ランごとに解析する**（位相はボードの間で比べられない。同じボードの中、または trig で組にしたイベントで使う）。

- **パワースペクトル**: 時刻を `bin` で数え、`segment` 秒ずつ FFT して平均する（Leahy。雑音なら 2）。`spectra` に速い周期（パルス、既定 1 us / 1 s）と遅い周期（スピル、既定 1 ms / 100 s）を並べる。
  いちばん高いピークの 1/k にも有意なピークがあれば、基本周波数はそちら。周波数が分かっていれば `frequency`（`--frequency`）に書く。
- **偽のピーク**: 周期が `min_period_ticks`（既定 10 カウント。HA なら 10 us、TPC なら 1 ms）より短いものは探さない（カウンタの刻みや読み出しの時間構造が出る。線源のラン 20261029 run06 で 421 kHz に出た）。
  畳んだ分布の Z^2_n が `n_sigma` に届かなければ、周期は無かったことにする。
- **周期を詰める**: Z^2_n（`n_harmonics`）を、使う時間を 5 倍ずつ伸ばしながら最大にする。
- **on-beam の窓**: 畳んだ位相の分布で、底から最大値までの `on_level`（既定 0.5）を超える、最大の位相を含む連続した範囲。位相と時刻の 2 次元図に傾きがあれば、カウンタとビームの周波数がずれている。

```bash
cd ~/analysis_sm3_FM
# JSON のとおり
/usr/bin/python3 ETCC_analysis/beam_period.py
# ビーム試験のラン、HA だけ、まず数ファイル
/usr/bin/python3 ETCC_analysis/beam_period.py --runs 20261201_run05 --sources ha --max-files 2
# 周波数が分かっているとき（そのまわりで詰める）
/usr/bin/python3 ETCC_analysis/beam_period.py --runs 20261201_run05 --frequency 1000
```

出力: `<out>_<ラン>.png`（基準のボードのレート、Δt の分布、スペクトル、畳んだ位相の分布と on-beam の窓、位相と時刻）、
`<out>_<ラン>_events.csv`（基準のボード（周期が出た中で HA を先に、イベントの多いもの）のイベントごとの epoch、trig、カウンタ、位相、on-beam。trig で他の解析と組にする）、
`<out>.json`（ボードごとの周波数、周期（秒とカウント）、位相の原点、ピークの位相、窓）。

既定の入力・出力先は、どれもこのリポジトリの中からの相対で決まる。画面に図を出さずに保存だけするときは
`MPLBACKEND=Agg` を付ける（XLaunch が無いと peak_search などが止まる）。
