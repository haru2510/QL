#!/usr/bin/env bash
# QL 用 PC で、TPC / PSA の QL をまとめて起動する。Ctrl-C で全部止まる。
# 全データの同期はしない。同期するときは、別のターミナルで ./ql_sync.sh --watch を実行する。
#
#   ./ql_live.sh                     # いちばん新しい run の QL を、起動した時点から出す
#   ./ql_live.sh --from-start        # 今の run を最初から積算する
#   QL_VIEWER_ARGS="--no-open --http-bind 0.0.0.0" ./ql_live.sh   # ql_tpc.py / ql_psa.py に渡す
#   QL_NODES="calltt:storage_g0" ./ql_live.sh                     # ノードを選ぶ
# 引数は TPC_analysis/ql_stream.py に渡す（--thin など。ql_stream.py --help を参照）。
#
# 起動するもの
#   TPC_analysis/ql_stream.py: 各 Storage Node で ql_feed.py を ssh で動かし、ノード側で間引いた
#                              イベントを受けて localhost:10012 で流す
#   TPC_analysis/ql_tpc.py   : TPC の QL（ブラウザ http://127.0.0.1:8090、FADC は http://127.0.0.1:8090/fadc）
#   PSA_analysis/src_onoda/ql_psa.py: HA の QL（ブラウザ http://127.0.0.1:8091）
# 図の PNG は check_output/ に書く。
#
# ノードへの ssh の事前準備は ql_sync.sh の冒頭を参照する。ノードには python3 が要る。

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
QL_NODES=${QL_NODES:-"calltt:storage_g0 killua:storage_g1"}
PYTHON=${PYTHON:-python3}
OUT_DIR=${QL_OUT_DIR:-$HERE/check_output}
read -r -a VIEWER_ARGS <<<"${QL_VIEWER_ARGS:-}"

# Ctrl-C や 1 つが落ちたときに、起動した全部を止める
trap 'trap - EXIT; kill 0 2>/dev/null' INT TERM EXIT

"$PYTHON" "$HERE/TPC_analysis/ql_stream.py" --nodes "$QL_NODES" "$@" &
"$PYTHON" "$HERE/TPC_analysis/ql_tpc.py" --host 127.0.0.1 -o "$OUT_DIR" ${VIEWER_ARGS[@]+"${VIEWER_ARGS[@]}"} &
"$PYTHON" "$HERE/PSA_analysis/src_onoda/ql_psa.py" --host 127.0.0.1 -o "$OUT_DIR" ${VIEWER_ARGS[@]+"${VIEWER_ARGS[@]}"} &

wait -n
echo "ql_live: どれかが終了したので全部止める" >&2
