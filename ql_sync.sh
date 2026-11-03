#!/usr/bin/env bash
# Storage Node（calltt, killua）の run を、この QL 用 PC へ同期する。QL 用 PC で実行する。
# smile-ics のプロセスとは独立に動き、RUN 中でも使える。
#
#   ./ql_sync.sh                         # 1 回だけ同期する
#   ./ql_sync.sh --watch                 # 30 秒ごとに同期し続ける（Ctrl-C で止める）
#   ./ql_sync.sh --watch 10              # 同期の間隔を秒で指定する
#   ./ql_sync.sh --since 20261101        # この日付より前の run は送らない
#   ./ql_sync.sh --dry-run               # 転送されるファイルだけ見る
#   QL_DEST=/data/ql ./ql_sync.sh
#
# 同期するもの: Storage Node の /data/smile_ics/storage_gX → $QL_DEST/storage_gX/
#   - 書き込み中のファイル（*.tmp）は送らない。閉じたファイルから順に届く
#   - Storage Node で消えたファイルは消さない（collect_runs.sh の削除に追従しない）
#   - --since（または QL_SINCE）を付けると、run ID（YYYYMMDD_runNN）の日付がそれより前の run は
#     送らない。metrics など run 以外のディレクトリは常に送る
#   - Storage Node 側の rsync は nice / ionice で優先度を下げ、DAQ の書き込みを妨げない
#
# QL 用 PC（192.168.1.x）は検出器ネットワーク（192.168.100.x）に直接届かないので、両方に
# つながる magallan を ssh の ProxyJump で経由する。データは magallan のディスクには書かれない。
#
# 事前準備（QL 用 PC で 1 回だけ）
#   1. ~/.ssh/config に次を書く（magallan の HostName は 192.168.1 側のアドレス）
#        Host magallan
#          HostName 192.168.1.xx
#          User msgc
#        Host calltt killua
#          User msgc
#          ProxyJump magallan
#        Host calltt
#          HostName 192.168.100.12
#        Host killua
#          HostName 192.168.100.13
#   2. QL 用 PC の公開鍵を 3 台に登録する
#        ssh-copy-id magallan; ssh-copy-id calltt; ssh-copy-id killua
#   3. QL 用 PC と Storage Node の両方に rsync を入れる

set -euo pipefail

# QL_DEST: 同期先のディレクトリ。ql_follow.py の --data-dir と合わせる
QL_DEST=${QL_DEST:-$HOME/smile_ics_ql}
# QL_NODES: Storage Node（~/.ssh/config の Host 名）と storage process id の組
QL_NODES=${QL_NODES:-"calltt:storage_g0 killua:storage_g1"}
# QL_SINCE: この日付（YYYYMMDD、JST）より前の run は送らない。空なら全 run を送る
QL_SINCE=${QL_SINCE:-}
SRC_ROOT=/data/smile_ics

INTERVAL=0
DRY_RUN=0

die() { printf 'ql_sync: %s\n' "$*" >&2; exit 1; }

while [ $# -gt 0 ]; do
  case $1 in
    --watch)
      INTERVAL=30
      if [ $# -gt 1 ] && [[ $2 =~ ^[0-9]+$ ]]; then
        INTERVAL=$2
        shift
      fi
      [ "$INTERVAL" -gt 0 ] || die "--watch の間隔は 1 秒以上にする"
      ;;
    --since)
      [ $# -gt 1 ] || die "--since には日付（YYYYMMDD）を付ける"
      QL_SINCE=$2
      shift
      ;;
    --dry-run) DRY_RUN=1 ;;
    -h|--help) sed -n '2,/^$/p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) die "不明なオプション: $1" ;;
  esac
  shift
done

# 2026-11-01 のような書き方も受け付け、run ID と同じ YYYYMMDD にそろえる
QL_SINCE=${QL_SINCE//-/}
[ -z "$QL_SINCE" ] || [[ $QL_SINCE =~ ^[0-9]{8}$ ]] || die "--since の日付は YYYYMMDD で書く: $QL_SINCE"

mkdir -p "$QL_DEST" || die "$QL_DEST を作れない"

# 二重に起動すると同じファイルを取り合うので、1 つに限る
exec 9>"$QL_DEST/.ql_sync.lock"
flock -n 9 || die "別の ql_sync.sh が動いている"

SSH_CMD="ssh -o BatchMode=yes -o ConnectTimeout=10"

# 1 ノード分を同期し、1 行で結果を出す
sync_node() {
  local node=$1 process_id=$2 listing name out status files size
  local opts=(-a --stats --exclude='*.tmp' -e "$SSH_CMD"
              --rsync-path='nice -n 19 ionice -c3 rsync')
  [ "$DRY_RUN" -eq 1 ] && opts+=(--dry-run --itemize-changes)

  # --since: ノードの run を一覧し、日付が古い run を除外する。一覧に失敗したら、全 run を
  # 送ってしまわないよう、このノードは飛ばす
  if [ -n "$QL_SINCE" ]; then
    if ! listing=$($SSH_CMD "$node" "ls -1 '$SRC_ROOT/$process_id'" 2>&1 </dev/null); then
      printf '  %-10s %s: 失敗 (run の一覧を取れない)\n' "$process_id" "$node"
      printf '%s\n' "$listing" | head -3 | sed 's/^/    /'
      return 1
    fi
    while read -r name; do
      [[ $name =~ ^([0-9]{8})_run ]] && [ "${BASH_REMATCH[1]}" \< "$QL_SINCE" ] \
        && opts+=(--exclude "/$name/")
    done <<<"$listing"
  fi

  mkdir -p "$QL_DEST/$process_id"
  out=$(rsync "${opts[@]}" "$node:$SRC_ROOT/$process_id/" "$QL_DEST/$process_id/" 2>&1 </dev/null) \
    && status=0 || status=$?
  # 24 は転送中に元のファイルが消えた場合（collect_runs.sh の削除など）。次回で揃うので失敗にしない
  if [ "$status" -ne 0 ] && [ "$status" -ne 24 ]; then
    printf '  %-10s %s: 失敗 (rsync %d)\n' "$process_id" "$node" "$status"
    printf '%s\n' "$out" | grep -E '^(rsync|ssh)[: ]|No such|denied|refused|timed out|resolve' \
      | head -3 | sed 's/^/    /'
    return 1
  fi
  [ "$DRY_RUN" -eq 1 ] && printf '%s\n' "$out" | grep '^>f' | sed 's/^/    /'
  files=$(printf '%s\n' "$out" | sed -n 's/^Number of regular files transferred: //p')
  size=$(printf '%s\n' "$out" | sed -n 's/^Total transferred file size: //p')
  printf '  %-10s %s: %s files, %s\n' "$process_id" "$node" "${files:-0}" "${size:-0 bytes}"
}

sync_once() {
  local entry rc=0
  printf '[%s] %s へ同期%s\n' "$(date '+%F %T')" "$QL_DEST" "${QL_SINCE:+（$QL_SINCE 以降の run）}"
  for entry in $QL_NODES; do
    sync_node "${entry%:*}" "${entry##*:}" || rc=1
  done
  df -h "$QL_DEST" | tail -1 | awk '{printf "  空き容量 %s / %s (%s 使用)\n", $4, $2, $5}'
  # 全ノードを 1 周したことを ql_follow.py に知らせる。dry-run では何も届いていないので触らない
  [ "$DRY_RUN" -eq 1 ] || touch "$QL_DEST/.ql_sync_done"
  return $rc
}

if [ "$INTERVAL" -eq 0 ]; then
  sync_once
  exit
fi

# --watch: 失敗しても止めずに次の周期で再試行する
while true; do
  sync_once || printf 'ql_sync: 同期に失敗した。%s 秒後に再試行する\n' "$INTERVAL" >&2
  # ロックの記述子を sleep に渡さない。渡すと止めたあとも sleep が終わるまでロックが残る
  sleep "$INTERVAL" 9>&-
done
