# 外部入力の設定

SSD・NASなどの入力先、対象ラン、較正値、保存先を指定するJSONを置く。
[`20261124_run06_07.json`](20261124_run06_07.json) はrun06・07の頂点選別・積算・軸fitの既定設定。
各キーの説明と単位はJSONの `_doc` と `_units` に記載。
相対パスの基準は、このフォルダではなくリポジトリ直下 `analysis_sm3_FM/`。

```bash
# このフォルダの設定を明示して積算する。
/usr/bin/python3 ETCC_cal/track_sum.py --config ETCC_cal/config/20261124_run06_07.json

# 同じ入力・較正・保存先で頂点を選別する。
/usr/bin/python3 ETCC_cal/select_vertices.py --config ETCC_cal/config/20261124_run06_07.json
```

まず `data_root`（`calltt/`・`killua/` の親）、`runs`、`result_dir` を確認する。
較正や入力を変えた場合、除外リストを再作成するには `track_sum.py --rebuild-veto --config ...` を使う。
run06・07用コード（`track_sum.py`、`select_vertices.py`、`beam_axis.py`、`plot_vertices.py`、`plot_saved_event.py`、`plot_scintillator.py`、`reconstruct_ha_lg.py`、`pixel_polarization.py`）とも、`--config` を省略すると `config/20261124_run06_07.json` を読む。
`inputs.vertices_all` が全頂点CSV、`outputs.vertices_selected` が選別後の頂点CSV。`outputs` で図・集計・fitの出力先も設定する。
`beam_axis.py` と `plot_vertices.py` の `--input`・`--out` はJSONより優先する任意の上書き。

`20261124_run00_07.json` は、全ランの選別なし飛跡積算を行う `track_runs.py` 用。入力先と追加走査ラン、run06・07の保存済み積算、カラーの百分位、全出力先を設定できる。

`pixel_polarization.py` も同じrun06・07用JSONを使う。HG/LG両版を同時に作り、`polarization`でraw2条件・散乱点のビーム直線距離条件・共通ペア・重みを切り替える。`20261124_run06_07_polarization_beam10.json` は距離10 mm選別をONにして別名保存する例。

`single_hit_pixels.py` も同じrun06・07用JSONを使う。`inputs.single_hit_cache_dir` が重心付きrawヒットNPZの入力、`outputs.single_hit_pixel_prefix` が出力名。`single_hit_pixels` で任意のビーム条件とHG/LG共通イベントへの制限を切り替える。

`beam_axis_overlay` は1/2ヒットの位置分布図に重ねるビーム印だけの設定。イベント選別とは独立で、`enabled`、`z_range_mm`、`marker_z_mm`を指定する。保存済みCSVからの再描画は `plot_hit_pixel_maps.py --config ...`。

## run00〜07の段階別解析

`20261124_run00_07_stages.json` を `run_stage_workflow.py` と `notebooks/analysis_workflow.ipynb` に渡す。既存の `20261124_run00_07.json` はtmp込み・未選別の全積算用なので区別する。段階別解析は完成rawのみで、run06・07の保存ビーム軸を共通使用する。

`20261124_horizontal_over_circular.json` は最終2ヒット品質選別後の水平／円の形状比用。入力マップ、出力先、円側の最小件数、比のカラー範囲を指定する。
