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
6本（`track_sum.py`、`select_vertices.py`、`beam_axis.py`、`plot_vertices.py`、`plot_saved_event.py`、`plot_scintillator.py`）とも、`--config` を省略すると `config/20261124_run06_07.json` を読む。
`inputs.vertices_all` が全頂点CSV、`outputs.vertices_selected` が選別後の頂点CSV。`outputs` で図・集計・fitの出力先も設定する。
`beam_axis.py` と `plot_vertices.py` の `--input`・`--out` はJSONより優先する任意の上書き。

`20261124_run00_07.json` は、全ランの選別なし飛跡積算を行う `track_runs.py` 用。入力先と追加走査ラン、run06・07の保存済み積算、カラーの百分位、全出力先を設定できる。
