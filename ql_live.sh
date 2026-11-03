#!/usr/bin/env bash
# QL 用 PC で、Storage Node からの同期と TPC / PSA の QL をまとめて起動する。Ctrl-C で全部止まる。
#
#   ./ql_live.sh                     # いちばん新しい run の QL を出す
#   ./ql_live.sh --since 20261101    # 同期する run を日付で絞る（ql_sync.sh に渡す）
#   QL_DEST=/data/ql QL_INTERVAL=10 ./ql_live.sh
#   QL_VIEWER_ARGS="--no-open --http-bind 0.0.0.0" ./ql_live.sh   # ql_tpc.py / ql_psa.py に渡す
#   QL_FOLLOW_ARGS=--from-start ./ql_live.sh   # 今の run を最初から読む（既定は起動した時点から）
#
# 起動するもの
#   ql_sync.sh --watch       : Storage Node の run を $QL_DEST へ同期し続ける
#   TPC_analysis/ql_follow.py: 同期したファイルから間引いたイベントを localhost:10012 で流す
#   TPC_analysis/ql_tpc.py   : TPC の QL（ブラウザ http://127.0.0.1:8090）
#   PSA_analysis/src_onoda/ql_psa.py: HA の QL（ブラウザ http://127.0.0.1:8091）
# 図の PNG は check_output/ に書く。
#
# 同期の事前準備は ql_sync.sh の冒頭を参照する。

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
export QL_DEST=${QL_DEST:-$HOME/smile_ics_ql}
QL_INTERVAL=${QL_INTERVAL:-30}
PYTHON=${PYTHON:-python3}
OUT_DIR=${QL_OUT_DIR:-$HERE/check_output}
read -r -a VIEWER_ARGS <<<"${QL_VIEWER_ARGS:-}"
read -r -a FOLLOW_ARGS <<<"${QL_FOLLOW_ARGS:-}"

# Ctrl-C や 1 つが落ちたときに、起動した全部を止める
trap 'trap - EXIT; kill 0 2>/dev/null' INT TERM EXIT

"$HERE/ql_sync.sh" --watch "$QL_INTERVAL" "$@" &
"$PYTHON" "$HERE/TPC_analysis/ql_follow.py" --data-dir "$QL_DEST" ${FOLLOW_ARGS[@]+"${FOLLOW_ARGS[@]}"} &
"$PYTHON" "$HERE/TPC_analysis/ql_tpc.py" --host 127.0.0.1 -o "$OUT_DIR" ${VIEWER_ARGS[@]+"${VIEWER_ARGS[@]}"} &
"$PYTHON" "$HERE/PSA_analysis/src_onoda/ql_psa.py" --host 127.0.0.1 -o "$OUT_DIR" ${VIEWER_ARGS[@]+"${VIEWER_ARGS[@]}"} &

wait -n
echo "ql_live: どれかが終了したので全部止める" >&2
