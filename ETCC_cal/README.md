# ETCCビーム試験解析 — 小野田チェック用

2026-11-24のrun06・07について、頂点分布、Z側のイベント除外、ビーム軸、TPC飛跡積算、シンチのヒット分布、イベント表示を整理した。
コードの各処理行には日本語の解説を付けた。Pythonは `/usr/bin/python3` を使う。

## 最初に見る結果

正規の保存先は [`output/20261124_run06_07/`](output/20261124_run06_07/)。

| 結果 | ファイル |
| --- | --- |
| 頂点XY・XZ・YZ、Z側除外あり・追加のXY選別なし | [vertices/no_xy_selection.png](output/20261124_run06_07/vertices/no_xy_selection.png) |
| Z側除外も含めイベント除外なしの全頂点 | [vertices/no_event_veto.png](output/20261124_run06_07/vertices/no_event_veto.png) |
| ビーム軸fit、全体表示 | [beam_axis/axis_z0veto.png](output/20261124_run06_07/beam_axis/axis_z0veto.png) |
| 同じfitのビーム付近の表示 | [beam_axis/axis_z0veto_local.png](output/20261124_run06_07/beam_axis/axis_z0veto_local.png) |
| fitの係数・各Z断面の中心・感度確認 | [beam_axis/axis_z0veto.json](output/20261124_run06_07/beam_axis/axis_z0veto.json) |
| 全TPCヒットのXZ・YZ積算、カット前後 | [tracks/track_sum.png](output/20261124_run06_07/tracks/track_sum.png) |
| SSDなしの飛跡積算再描画用ヒストグラム | [tracks/track_sum.npz](output/20261124_run06_07/tracks/track_sum.npz) |
| 全イベントの分母・残存数、選別後のシンチヒットの割合（`scintillator`） | [selection/summary.json](output/20261124_run06_07/selection/summary.json) |
| Z側除外の前後のシンチヒット数の分布（円グラフ・ヒストグラム） | [selection/scint_hit_multiplicity.png](output/20261124_run06_07/selection/scint_hit_multiplicity.png) |
| 全rawシンチヒットのブロック別分布 | [scintillator/all_raw_blocks.png](output/20261124_run06_07/scintillator/all_raw_blocks.png) |
| run06・trigger 22204の候補表示 | [events/run06_trig22204.png](output/20261124_run06_07/events/run06_trig22204.png) / [回転できるHTML](output/20261124_run06_07/events/run06_trig22204.html) |

## run00〜07の全データ・選別なし飛跡積算

完成rawに加え、run00〜05の末尾 `.raw.tmp` から読める有効フレームも含めた。Z側除外、頂点選別、XY選別、HAのヒット数条件、ToT補正はいずれも適用していない。rawの全有効クロック0〜1023を表示する。

- [全ラン合算・mm版・カラー調整](output/20261124_run00_07/tracks/all_runs_z_contrast.png)
- [全ラン合算・クロック版・カラー調整](output/20261124_run00_07/tracks/all_runs_clock_contrast.png)
- [全ラン合算・mm版・全レンジ](output/20261124_run00_07/tracks/all_runs_z_fullscale.png)
- [全ラン合算・クロック版・全レンジ](output/20261124_run00_07/tracks/all_runs_clock_fullscale.png)
- [積算配列](output/20261124_run00_07/tracks/all_runs.npz) / [入力・件数・カラー範囲](output/20261124_run00_07/tracks/all_runs.json)
- [run00〜05の個別図とキャッシュ](output/20261124_run00_07/runs/)

全ヒット数は **19,552,033,609ストリップ・クロックヒット**。イベント数ではない。
`.raw.tmp` から追加したのは、延べ181,289ボードフレーム・63,959,991ヒット。完成rawとの境界を確認して重複を除き、不完全な末尾フレームはデコーダで除外した。
完成rawのキャッシュとtmpの追加分は別々に保存し、再描画時に足す。カウンタ巻き戻りがtmpの境界で見つかった場合は誤集計せず停止する。

カラー調整版は、空ビンを除く表示ビンの10〜99.5百分位（今回12,228.9〜2,331,741.875）を対数色の範囲にした。範囲外の値をデータから削除してはいない。
全レンジ版も下限1〜最大4,235,423で保存。XZとYZは同じ色範囲。個別ランはそれぞれ色範囲を決めているため、色だけでラン間の強度を比較しない。
`config/20261124_run00_07.json` の `color_percentiles` で変更し、`track_runs.py --plot-only` でSSDなしに再描画できる。

run00はアノード側の完成rawが約8.5〜9万フレーム、カソード側が約125万フレームと偏っていた。今回は保存されている全ヒットを使い、両側が揃ったイベントへの限定はしない。利用可能な入力の偏りは、位置分布の非対称性にも関係し得る。

頂点・Z側除外・ビーム軸の結果は、引き続きrun06・07についての結果であり、全8ランの軸fitを行ったわけではない。

## フォルダとコード

```text
ETCC_cal/
  README.md                 この説明
  config/                   外部入力先・ランごとの設定（--configで指定）
  settings.py               JSON読み取りとパス解決を共通化
  param/20261124_run06_07.json  整理前の設定を残した互換用コピー
  track_sum.py              raw全ヒットの積算、除外集合の再確認、イベント数集計
  track_runs.py             run00〜07全ヒットのmm/clock積算とカラー調整
  select_vertices.py        保存済み全頂点へイベント単位のZ側除外を適用、選別前後のシンチヒットフラグを集計
  beam_axis.py              Zごとの局所ピーク中心と軸のfit、全体図と局所図
  pixel_polarization.py      HG/LGピクセル方位角、任意のビーム直線・raw2選別
  plot_vertices.py          追加のXY選別なしの頂点投影
  plot_saved_event.py        保存した候補点群からPNG・HTMLを再描画
  plot_scintillator.py       保存した全rawシンチ件数から再描画
  plot_scint_hits.py         Z側除外の前後のシンチヒット数の分布を描く（SSD不要）
  output/20261124_run00_07/    全8ランの選別なしraw積算・クロック表示
  output/20261124_run06_07/
    vertices/               全頂点・選別後頂点のCSVと図
    beam_axis/              採用する局所ピークfit
    tracks/                 raw飛跡積算の配列・入力一覧・図
    selection/              除外イベントのNPZと残存数
    scintillator/           全rawシンチヒットの集計
    events/                 個別候補のPNG・HTML・点群JSON
```

rawのフレーム仕様とストリップ配置は既存の [`ETCC_analysis/beam_vertex.py`](../ETCC_analysis/beam_vertex.py) と [`TPC_analysis/check_tpc.py`](../TPC_analysis/check_tpc.py) を使用する。
今回の整理コードでは、クロック読み取り・イベントの重複除去・座標変換・イベント除外・積算・軸fitを区別して読めるようにした。

`output/` は大きなデータのためGit対象外。コード・設定・このREADMEはGit追跡対象。
旧 `ETCC_analysis/output` の平置き61ファイルは削除せず [`archive/20261124_run06_07/`](../ETCC_analysis/output/archive/20261124_run06_07/) へ種類別に移動した。
旧パスと新パスの対応はその中の `file_moves.json`。`prototypes/` は探索時の実行履歴で、現在の正規コードではない。

## 外部入力の設定

[`config/`](config/README.md) に設定JSONを置く。run06・07用コードは引数省略時に `config/20261124_run06_07.json` を読む。`param/` は互換用コピー。

`data_root`・`storages`・`runs` がraw入力、`inputs.vertices_all` が全頂点CSV。`outputs` で選別後CSV、各PNG、NPZ、集計JSON、fitの保存名を指定する。選別後のCSVは下流の描画・軸fitが `outputs.vertices_selected` から読む。

例: `/usr/bin/python3 ETCC_cal/beam_axis.py --config ETCC_cal/config/20261124_run06_07.json`。コマンドの `--input`・`--out` は任意の上書き。

`track_runs.py` は全ランの選別なし積算専用で、既定の設定は `config/20261124_run00_07.json`。`runs_to_scan` は新規走査するrun00〜05、`inputs.run06_07_track_histogram` は保存済みのrun06・07のカット前積算。`include_tmp=true` で末尾の有効フレームも加える。

## 入力とトリガーの対応

- SSD: `/mnt/d/UVSOR2026_late/{calltt,killua}/20261124_runNN/`。
- `calltt`: アノード、`killua`: カソード。両側のrawを読む。
- 同じラン内で `(epoch, trig)` をイベント識別子にする。トリガーカウンタが戻ったらepochを進める。
- 同じボードの重複フレームは一度だけ数え、全イベント数は全ボードのイベント番号の和集合で求める。
- run06・07の再走査では、保存済みのZ側除外集合と今回のrawの判定が完全一致した。
- 入力ファイルのパス・サイズ・更新時刻は `tracks/track_sum.json` の `input_manifest` に保存。

以前の読み込みでは一部のrun07ファイルが空だったが、再接続後には読み取れた。古い偏光図と最新の全ヒット集計では入力状態が違うため、そのまま差分を物理量として比較しない。

## 座標・較正値・負のZ

座標系は X: 東、Y: 南、Z: 下。距離の単位はmm。

```text
Z [mm] = velocity [mm/clock] × (clock − clock_offset [clock]) + z_gem [mm]
       = 0.39358 × (clock − 133.69) + 5.0
```

| パラメータ | 値 | 単位 | 意味 |
| --- | ---: | --- | --- |
| pitch_strip | 0.4 | mm | ストリップ間隔 |
| velocity | 0.39358 | mm/clock | 既存のドリフト速度 |
| clock_offset | 133.69 | clock | GEM面に対応するクロック |
| z_gem | 5.0 | mm | GEM面の位置 |
| z_veto | 15.0 | mm | イベントを除外する低Z側のしきい値 |

この式では**整数クロック120以下のヒットが負のZ**になる。rawの有効クロックは0〜1023で、mm換算では約−47.62〜355.02 mmまである。
負のZを自動的に0へ丸めたり消したりはしていない。

負のZは「物理的な飛跡が必ず検出領域外にあった」という確定ではない。トリガーとの時間関係、ノイズ、時刻原点やドリフト較正のずれなどでも生じ得る。
今回のrunからドリフト速度を再較正したわけではない。GEM側と天板側の端が明瞭でないと、既知の長さをクロック差で割る方法では速度・時刻原点を安定に決められない。
端の確認には、Z側カット後に人工的に生じた端を使わず、全データのクロック表示を使う。読み出し窓の端（clock=0・1023）や図の表示端も、GEM・天板の物理的な端と同一視しない。

## 3種類の条件を区別する

### 1. イベント全体を除外するZ側条件

**全TPCボードのアノードまたはカソードにZ ≤ 15 mmのrawストリップヒットが1個でもあれば、そのイベント全体を除外する。**

- 現在の較正では整数クロック159以下に対応する。
- アノードとカソードの同時計数である必要はない。
- 頂点だけのZ判定ではない。頂点が内側でも、別の場所に低Zヒットがあれば除外する。
- XYの壁と高Z側にはfiducial除外を適用しない。
- 低Zヒットが実際の外来飛跡かノイズかは、この条件だけでは区別しない。

保存済み頂点はイベントごとに独立に計算されているため、除外対象のイベント番号をrawから求めてから対応する行を除けば、生存イベントの頂点を再計算する必要はない。
`vertices/z0veto.csv` の `z0_veto_pass` はこの条件の通過を意味する。旧全方向fiducialフラグは使わず、旧PCA軸からの距離も空欄にしている。

### 2. ビーム軸のfitに使う局所範囲

- 主解析のZ断面: 15〜295 mmを20 mmごと。上下端の余白をイベント除外として追加しているわけではない。
- 各断面のXY fit範囲: X=40〜90 mm、Y=8〜44 mm。
- XYのビン幅: 2 mm。各ビン内で積分した2次元Gaussian + 正の滑らかな背景 `b exp(gx(X−65)+gy(Y−25))` をfitする。
- GaussianのX/Y相関は導入していない。8パラメータは信号数、X/Y中心、X/Y幅、背景、背景のX/Y勾配。
- 中心探索の制限: X=45〜80 mm、Y=14〜34 mm。幅は0.6〜18 mm。これは今回観測された局所集中の周辺を解析する設定。
- 背景から分離した各断面の中心を重み付き直線にfitする。誤差の下限は0.05 mm。
- fitの偏差やχ²が大きい場合は共分散を膨らませる。系統誤差が全て評価されたという意味ではない。

全頂点へのPCA（主成分分析）は、背景と検出器配置による広がりを拾って約41°となったため不採用。
現在採用している軸は**局所ピーク中心への直線fit**。

### 3. 局所図だけの表示条件

標準の全体図は追加のXY選別なし。局所表示版だけは背景を減らして見るため次の帯を使う。

| 図 | 表示する範囲 |
| --- | --- |
| XZ局所表示 | `|Y − 22| < 15 mm`、すなわちY=7〜37 mm |
| YZ局所表示 | `|X − 63| < 20 mm`、すなわちX=43〜83 mm |

この帯幅は表示のために選んだもので、物理的な最適性を検証したカットではない。イベント除外条件にも軸fitのXY範囲にも使わない。
全体図と局所図は同じfit結果を表示する。

## 頂点の定義と制約

同じuPICのアノードとカソードに同時ヒットがある最も早いクロックを選び、そのクロックのストリップ位置の平均をX/Y、上記の較正式をZとする。
各ボードで早い16ヒットクロックを保持した従来の推定法を使用。ToT補正、電子と陽電子の分離、ノイズストリップ除去、対生成の確定識別はしていない。

`vertices/all.csv` は、旧処理で10万トリガーずつ分割して求めた全頂点を引き継いだもの。run06で分割前後の109,267点が一致した。
全run07を一括処理する旧 `beam_vertex.py` はメモリ不足で終了したため、大量データでそのまま一括実行しない。
rawから頂点を完全再構成するコードの整理は既存 `beam_vertex.py` 側に残る。現在の `select_vertices.py` は保存済み全頂点CSVを入力にする。

## run06・07の最新結果

| ラン | 全記録イベント | Z側条件の通過 | 通過率 | 通過かつ頂点あり |
| --- | ---: | ---: | ---: | ---: |
| run06 | 122,889 | 42,072 | 34.24% | 34,937 |
| run07 | 1,330,557 | 456,982 | 34.35% | 375,721 |
| 合計 | 1,453,446 | 499,054 | 34.34% | 410,658 |

分母はヒットなしも含む記録トリガー。頂点を求められたイベントだけを分母にしていない。
カット後の頂点は全記録イベントの28.25%。頂点のZだけを15 mmで切る場合より27,375イベントを追加で除外した。

採用する軸（座標はmm）:

```text
X = 62.26454 − 0.01727826 × (Z − 155)
Y = 22.13397 + 0.00123810 × (Z − 155)
```

- Z軸からの3次元角度: 0.99241°。
- ZX平面: Xが減る向きに約0.990°。ZY平面: Yが増える向きに約0.071°。
- XY傾きの相関を無視した誤差伝播では約0.034°。絶対精度の保証ではない。
- fit範囲・Z区切りを変えると0.963〜0.997°。
- Xの直線fitはχ²/自由度=63.45/12で、完全な直線とは一致していない。
- ドリフト較正、位置合わせ、頂点推定の偏りの系統誤差は未評価。

## 単純な飛跡積算の定義

`track_sum.py` は、全レコードの**1ストリップ・1クロックのヒットを1個**としてXZ/YZへ加える。
頂点だけの積算でも、イベントあたり1票でもない。ヒットの多いイベント、長いToTの信号はその分多く数える。

- rawの全有効クロック0〜1023を使う。頂点推定の「早い16クロック」の制限は使わない。
- XYの追加選別、HAのヒット数条件、ToT補正はしない。
- XZとYZは直接観測されたストリップ座標を使う。曖昧なXYの全組み合わせを作らない。
- 保存配列はX/Y × 768座標ビン × 1024クロック。図では8ストリップ・4クロックずつ足して表示する。
- run06・07のカット前後の図は同じカラースケール。空ビンはマスク、LogNormの下限は1。
- `track_sum.npz` があればSSDなしでカラースケールやビン表示を変更できる。

## Z≈310〜315 mmの頂点分布の段差

Z側除外を外した `vertices/no_event_veto.png` でも段差は残る。追加のXY表示条件やZ≤15 mmの除外が作った段差ではない。
頂点は最初の同時ヒットだけを1点として数え、飛跡積算は後続のストリップ・クロックヒットも全て数えるため、両者の端は一致するとは限らない。
全TPC積算には頂点を求められなかったイベントも含むため、イベント集合も完全には一致しない。
Z側除外後の例では、頂点数が305〜310 mmで7,298点、315〜320 mmで631点に減る一方、通過全イベントのrawヒット数は約3,460万→約3,327万で緩やか。
TPC上端の想定位置317 mmに近いが、物理端と断定してドリフト速度を較正してはいない。

## 今後のドリフト速度推定: 散乱点分布の端を使う

全ヒットの飛跡積算よりも、**散乱点（推定頂点）のクロック分布の端を使う方法を優先して検討する**。
全ヒット積算には信号の後続ヒットやToTの広がりも含まれる一方、散乱点は最初の同時ヒットを選ぶため、相互作用位置の端を追いやすい可能性がある。
実際に、散乱点分布ではZ≈310〜315 mm付近の急減が見えるが、全ヒット積算の同じ位置には明瞭な段差が出ていない。

推定には、Z≤15 mmのイベント除外を含むfiducialカットを外した全散乱点を使い、**既存の速度でmm換算したZではなく、元の頂点クロック**をfitする。
GEM側・天板側に対応する端を確認できた場合は、既知の位置差とクロック差から求める。

```text
velocity [mm/clock] = (Z_top − Z_gem) [mm] / (clock_top − clock_gem) [clock]
clock_offset [clock] = clock_gem
```

端のクロックは、遷移の幅と背景を含むモデルで推定し、fit範囲・ビーム付近の選別・頂点推定法を変えた変動も確認する。
現在見えている段差が天板などの物理端に対応するか、読み出し窓・ノイズ・頂点選別の影響ではないかを先に検証する。
一方の端しか確認できない場合は、独立に決めた時刻原点などの制約がなければ、ドリフト速度と時刻原点の両方を決められない。

これは今後の解析方針であり、**散乱点からの新しいドリフト速度はまだ算出していない**。現行結果は引き続き0.39358 mm/clockを使用している。

## シンチとイベント表示

シンチヒット数はrawのDataCountフラグを数える。HGやpeak_idによる追加カットはしない。
全rawヒットのブロック別分布と、rawでちょうど2ヒットのペア解析は別の集計。

chは現行EricaNextのADC番号順（ADC 1,2,4,8,16,32 → ch0〜5）。旧HDF5出力番号は `5 − ch`。
全raw分布で最多の2ブロックは、ha01・現行ch2（出力ch3）とha01・現行ch0（出力ch5）。

イベント図のシンチ位置はブロック中心。TPCの3D点は同一クロックのXY組み合わせであり、偽の点を含み得る。XZ・YZ投影も一緒に確認する。
run06・trigger 22204では二股の形が見え、ha01とha15の現行ch4（どちらも出力ch1）が同じトリガーでヒットしている。ただし対生成の確定ではない。
保存したtrigger 22204・55668はZ側イベント除外にも抵触しないことを確認した。

## LG重心によるシンチのピクセル位置（HG位置較正を暫定流用）

`reconstruct_ha_lg.py` は保存済み候補と同じトリガーのHA rawを読み、LG4信号のペデスタル差し引き後の積分から重心を求める。

```text
x_lg = (J0 − J1 − J2 + J3) / (J0 + J1 + J2 + J3)
y_lg = (−J0 − J1 + J2 + J3) / (J0 + J1 + J2 + J3)
```

LG用の位置較正は作っていないため、ユーザー指定によりHG用 `position_haNN.csv` の64ピークと(n,m)対応を流用する。
EricaNextと同じthin-plate spline（smoothing=0）をLG重心に適用し、`floor(corrected + 0.5)` で1〜8のピクセルを選ぶ。
位置は `geometry_sm3.csv` の向き・ピッチで検出器座標へ変換する。HG/LGの重心・積分・飽和フラグ、使用較正ファイルのハッシュ、元のブロック中心もJSONへ保存する。
LGも飽和、LG積分が非正、または補正後が8×8範囲外なら、ピクセル位置を採用せずブロック中心へ戻し、図とJSONに明記する。
HG用較正のLGへの流用は暫定であり、HGとLGの位置応答が同じであることを検証したわけではない。エネルギーのkeV較正をLGへ流用する処理は行っていない。

既存2候補の4ヒットはHG・LGとも飽和フラグなし。LG位置は全て8×8内だった。

- [trigger 22204・LG位置版PNG](output/20261124_run06_07/events/lg/run06_trig22204_lg.png) / [回転できるHTML](output/20261124_run06_07/events/lg/run06_trig22204_lg.html)
- [trigger 55668・LG位置版PNG](output/20261124_run06_07/events/lg/run06_trig55668_lg.png) / [回転できるHTML](output/20261124_run06_07/events/lg/run06_trig55668_lg.html)
- [位置・飽和フラグの一覧](output/20261124_run06_07/events/lg/lg_position_summary.json)

入力・出力は同じ `config/20261124_run06_07.json` に指定する。
`inputs.event_candidates`、`inputs.ha_position_calibration_dir`、`inputs.scintillator_geometry`、`outputs.lg_event_dir` を使う。

```bash
/usr/bin/python3 ETCC_cal/reconstruct_ha_lg.py --config ETCC_cal/config/20261124_run06_07.json
/usr/bin/python3 ETCC_cal/plot_saved_event.py \
  --input ETCC_cal/output/20261124_run06_07/events/lg/run06_trig22204_lg.json \
  --out ETCC_cal/output/20261124_run06_07/events/lg/run06_trig22204_lg
```

## HG/LGのピクセル位置を使う偏光解析

[`pixel_polarization.py`](pixel_polarization.py) は同じrawからHG重心版・LG重心版の2通りを作る。
ブロック中心ではなく、`position_haNN.csv` の64ピークをthin-plate splineで補正したピクセル(n,m)の中心を使用する。
ペアCSVにはrun・epoch・trig・rawヒット数・ビーム距離・両ヒットのboard/ch/n/m・XYZ・方位角・重みを保存する。
HGもLGも、対象gainがADC飽和、積分が非正、較正なし、配置なし、または8×8外なら位置を採用しない。この解析では無効ピクセルをブロック中心で代用しない。
LGにはHGの位置較正を暫定流用している。複数ピクセルが発光したときには光量重心が混ざるため、推定した1ピクセルが実際の全エネルギー付与位置を表すとは限らない。

入力・出力は同じ `config/20261124_run06_07.json` の `inputs`・`outputs`、切替条件は `polarization` に指定する。

| キー | 既定値 | 意味 |
| --- | --- | --- |
| `require_exactly_two_raw_hits` | true | ピクセルで除外する前のraw DataCountフラグがちょうど2個。falseなら有効ピクセルの全無順序ペアを作る |
| `beam_vertex_cut_enabled` | false | 同じrun・epoch・trigの散乱点が保存済みビーム直線から指定距離以内かを選別。ON時の頂点なしは除外 |
| `beam_max_distance_mm` | 10 | 直線との最短3D距離の上限[mm]。XYの差だけではない |
| `common_pairs_only` | false | trueならHG/LG両方で有効な同一rawヒットペアに揃える。方位角は各gainで別計算 |
| `pair_weighting` | per_event | 同一イベントの採用ペアに重みを等分し総和1。per_pairなら各ペア1票 |
| `beam_direction_source` | fitted_line | 方位角を定義する軸の向きを保存fitから得る。configuredなら`beam_direction`を使う |
| `min_pair_distance_mm` | 0 | ピクセル間の最小距離[mm]。ビームに垂直な投影長0の組は常に除外 |
| `phi_bins` | 18 | 0〜180度の分布のビン数 |
| `max_files_per_board` | 0 | 0は全完成HA raw、正値は試験用のボードごとのファイル制限 |
| `include_tmp` | false | trueならHAの`.dat.tmp`も対象にする |

`inputs.polarization_vertices` は既定で `vertices/all.csv`。Z≤15 mmのイベント除外を暗黙に追加していない。
Z側条件も必要なら、その入力を `vertices/z0veto.csv` へ変更する。
`inputs.polarization_beam_axis` は他コードで求めた `beam_axis/axis_z0veto.json`。
方位角用のビーム方向と距離条件は別で、configuredで方向を指定しても、距離選別には保存された直線を使う。
ビーム条件OFFかつconfiguredの場合は保存軸や散乱点ファイルを読む必要はない。
rawの3ヒットから無効ピクセルを落として2個になっても、raw2条件ONでは採用しない。

```bash
# 同じJSONでHG版とLG版を同時作成する（初回はHA rawを読む）。
/usr/bin/python3 ETCC_cal/pixel_polarization.py --config ETCC_cal/config/20261124_run06_07.json

# 重心付きrawヒットNPZを再利用して条件を切り替える（SSD不要）。
/usr/bin/python3 ETCC_cal/pixel_polarization.py --from-cache

# ビーム距離10 mmとraw2条件を有効にした別名出力の設定例。
/usr/bin/python3 ETCC_cal/pixel_polarization.py \
  --config ETCC_cal/config/20261124_run06_07_polarization_beam10.json --from-cache

# CLIでも独立に上書きできる。出力を残す場合はJSON側で別のprefixを指定する。
# --beam-cut / --no-beam-cut
# --require-two-hits / --no-require-two-hits
```

run06・07の完成HA rawによる実行結果:

| 条件 | HGピクセルペア | LGピクセルペア |
| --- | ---: | ---: |
| raw2=true、ビーム条件OFF | 145,233 | 253,355 |
| raw2=true、ビーム直線から10 mm以内 | 29,599 | 50,015 |

- [HG版](output/20261124_run06_07/polarization/pixel_hg.png) / [LG版](output/20261124_run06_07/polarization/pixel_lg.png)
- [ビーム条件ON・HG版](output/20261124_run06_07/polarization/pixel_beam10_hg.png) / [LG版](output/20261124_run06_07/polarization/pixel_beam10_lg.png)
- [HGピクセル付きCSV](output/20261124_run06_07/polarization/pixel_hg_pairs.csv) / [LGピクセル付きCSV](output/20261124_run06_07/polarization/pixel_lg_pairs.csv)
- [両gainの集計](output/20261124_run06_07/polarization/pixel_comparison.json)

件数差には対象gainの飽和・8×8範囲外・位置による投影長0の違いが含まれる。同じイベントでの比較には`common_pairs_only=true`を使う。
変調fitはアクセプタンス未補正の生の方位角分布へのfitで、得たmuを較正済み偏光度とは解釈しない。実際に単純なcos2phiに乗らない大きな構造もある。
多重ヒットの全ペアは互いに相関する。sumw2は表示用の誤差で、そのペア間共分散を含めた偏光測定の最終誤差ではない。

実装検証では、人工イベントのraw3→有効2の誤採用防止、傾いた直線との距離、頂点なし、切替の独立性、重み、投影長0を確認した。
実run06キャッシュでもraw2条件OFF・ビーム条件ON・共通ペアONを検証し、多重ヒットの採用とイベント総重み1を確認した。

## 再実行

以下はリポジトリ直下で実行する。新たな依存インストールは不要。

```bash
# 全rawヒットを再走査する。保存済み除外集合とrawの一致も確認する。
/usr/bin/python3 ETCC_cal/track_sum.py

# 入力・較正が変わったときは、除外集合をrawから再作成する（2回走査する）。
/usr/bin/python3 ETCC_cal/track_sum.py --rebuild-veto

# 保存した全頂点へイベント単位の除外を適用し、分母と残存率を集計する。
# シンチ（HA）のヒットフラグも数える。初回だけHA rawのヘッダを読み、selection/ha_hit_flags_<run>.npz に保存する。
# rawが変わったときは --rebuild-ha で読み直す。
/usr/bin/python3 ETCC_cal/select_vertices.py

# summary.json の scintillator:
#   ヒット = HAフレームのヘッダのDataCountのビット（データのあるch）。HG・peak_id・時刻のカットはしない。
#   TPCとは (epoch, trig) で対応づける。TPCのイベントは epoch 0、trig 1..all_events。
#   events_with_hit_*: HAヒットが1個以上あるイベント数（all: 全イベント、passing: Z側除外の通過、vertex: 通過して頂点もある）と、その割合 [%]。
#   hits_*: ヒットフラグの総数。multiplicity_*: 1イベントのヒット数の分布（HAに記録が無いイベントは0）。

# 除外の前後のヒット数の分布を円グラフとヒストグラムで描く（selection/のキャッシュだけを使う。SSD不要）。
/usr/bin/python3 ETCC_cal/plot_scint_hits.py

# 局所ピーク軸をfitし、全体図と局所表示を両方保存する（SSD不要）。
/usr/bin/python3 ETCC_cal/beam_axis.py

# 追加XY選別なしの頂点分布だけを描く（SSD不要）。
/usr/bin/python3 ETCC_cal/plot_vertices.py \
  --input ETCC_cal/output/20261124_run06_07/vertices/z0veto.csv \
  --out ETCC_cal/output/20261124_run06_07/vertices/no_xy_selection.png

# Z側除外も外した全頂点を描く（SSD不要）。
/usr/bin/python3 ETCC_cal/plot_vertices.py \
  --input ETCC_cal/output/20261124_run06_07/vertices/all.csv \
  --out ETCC_cal/output/20261124_run06_07/vertices/no_event_veto.png \
  --selection-label 'All reconstructed vertices; no Z-side, XY, or high-Z event veto'

# 保存済みの全ヒット積算を再描画する（SSD不要）。
/usr/bin/python3 ETCC_cal/track_sum.py --plot-only

# 同じJSONのinputs/outputsからイベント図・シンチ図を再描画する（SSD不要）。
/usr/bin/python3 ETCC_cal/plot_saved_event.py
/usr/bin/python3 ETCC_cal/plot_scintillator.py

# run00〜07の全完成rawヒットを選別なしで積算する。
/usr/bin/python3 ETCC_cal/track_runs.py --config ETCC_cal/config/20261124_run00_07.json

# 中断後は、完了済みランの条件・入力サイズ・更新時刻を確認して再利用する。
/usr/bin/python3 ETCC_cal/track_runs.py --resume

# 全データ積算のカラースケールだけを変更して再描画する（SSD不要）。
/usr/bin/python3 ETCC_cal/track_runs.py --plot-only
```

レビュー時は、最初に `track_sum.py` の `accumulate()`（何を1ヒットと数えるか）、`scan_run()`（イベント対応と除外集合の検証）、`select_vertices.py`（事象全体の除外）、`beam_axis.py` の `peak_fit()`・`fit_line()`（モデルと重み）を確認すると条件を追いやすい。
