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
#   - 書き込み中のファイル（*.tmp）も送る。データのファイル（*.raw, *.dat, *.tmp）は後ろに
#     書き足されるだけなので、rsync --append で増えた分だけを送る。ノードで .tmp が閉じて名前が
#     変わったら、こちらでも名前を変えてから続きを送る（送り直さない）
#   - データ以外（run_meta.json など）は書き換わることがあるので、--append を使わずに送る
#   - Storage Node で消えたファイルは消さない（collect_runs.sh の削除に追従しない）
#   - --since（または QL_SINCE）を付けると、run ID（YYYYMMDD_runNN）の日付がそれより前の run は
#     送らない。metrics など run 以外のディレクトリは常に送る
#   - 毎回、各ノードのいちばん新しい run を先に送り、そのあとで過去の run と metrics を送る
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

# 1 周に何度もつなぐので、接続を使い回す（ControlMaster）。magallan を経由する分の手間が減る。
# ソケットのパスは 108 バイトまでなので、名前は短くする
SSH_CMD="ssh -o BatchMode=yes -o ConnectTimeout=10 -o ControlMaster=auto"
SSH_CMD+=" -o ControlPath=${XDG_RUNTIME_DIR:-/tmp}/qls-%r@%h:%p -o ControlPersist=60"

# データのファイルは後ろに書き足されるだけなので、増えた分だけ送る
DATA_FILTER=(--append --include='*/' --include='*.raw' --include='*.dat' --include='*.tmp'
             --exclude='*')
# データ以外。書き換わることがあるので --append は使わない
META_FILTER=(--exclude='*.raw' --exclude='*.dat' --exclude='*.tmp')

# rsync を 1 回走らせ、1 行で結果を出す
run_rsync() {
  local label=$1 node=$2 src=$3 dest=$4 out status files size
  shift 4
  local opts=(-a --stats -e "$SSH_CMD" --rsync-path='nice -n 19 ionice -c3 rsync' "$@")
  [ "$DRY_RUN" -eq 1 ] && opts+=(--dry-run --itemize-changes)

  mkdir -p "$dest"
  out=$(rsync "${opts[@]}" "$node:$src" "$dest" 2>&1 </dev/null) && status=0 || status=$?
  # 24 は転送中に元のファイルが消えた場合（collect_runs.sh の削除など）。次回で揃うので失敗にしない
  if [ "$status" -ne 0 ] && [ "$status" -ne 24 ]; then
    printf '  %-26s %s: 失敗 (rsync %d)\n' "$label" "$node" "$status"
    printf '%s\n' "$out" | grep -E '^(rsync|ssh)[: ]|No such|denied|refused|timed out|resolve' \
      | head -3 | sed 's/^/    /'
    return 1
  fi
  [ "$DRY_RUN" -eq 1 ] && printf '%s\n' "$out" | grep '^>f' | sed 's/^/    /'
  files=$(printf '%s\n' "$out" | sed -n 's/^Number of regular files transferred: //p')
  # 実際にネットワークで受け取った量。--append では、ファイル全体ではなく増えた分になる
  size=$(printf '%s\n' "$out" | sed -n 's/^Total bytes received: //p')
  printf '  %-26s %s: %s files, 受信 %s bytes\n' "$label" "$node" "${files:-0}" "${size:-0}"
}

# ノードの run を一覧し、LATEST（いちばん新しい run）と EXCLUDES（--since より古い run）、
# NODE_TMPS（ノードで書き込み中のファイル。process id からの相対パス）を決める
list_runs() {
  local node=$1 process_id=$2 listing name
  LATEST=
  EXCLUDES=()
  if ! listing=$($SSH_CMD "$node" "cd '$SRC_ROOT/$process_id' && ls -1 && echo '--tmp--' \
                 && find . -name '*.tmp' -type f" 2>&1 </dev/null); then
    printf '  %-26s %s: 失敗 (run の一覧を取れない)\n' "$process_id" "$node"
    printf '%s\n' "$listing" | head -3 | sed 's/^/    /'
    return 1
  fi
  NODE_TMPS=$(sed '1,/^--tmp--$/d' <<<"$listing")
  listing=$(sed '/^--tmp--$/,$d' <<<"$listing")
  while read -r name; do
    [[ $name =~ ^([0-9]{8})_run[0-9]+$ ]] || continue
    if [ -n "$QL_SINCE" ] && [ "${BASH_REMATCH[1]}" \< "$QL_SINCE" ]; then
      EXCLUDES+=(--exclude "/$name/")
    elif [ -z "$LATEST" ] || [ "$LATEST" \< "$name" ]; then
      LATEST=$name
    fi
  done <<<"$listing"
}

# ノードで閉じて名前が変わった .tmp を、こちらでも同じ名前に変える。変えずにおくと、rsync が
# 閉じたファイルを別のファイルとして最初から送り直してしまう
rename_closed() {
  local process_id=$1 path rel
  local -A writing=()
  [ "$DRY_RUN" -eq 1 ] && return 0
  while read -r rel; do
    [ -n "$rel" ] && writing[$rel]=1
  done <<<"$NODE_TMPS"
  while IFS= read -r -d '' path; do
    rel="./${path#"$QL_DEST/$process_id/"}"
    [ -n "${writing[$rel]:-}" ] && continue
    if [ -e "${path%.tmp}" ]; then
      rm -f "$path"
    else
      mv "$path" "${path%.tmp}"
    fi
  done < <(find "$QL_DEST/$process_id" -name '*.tmp' -type f -print0 2>/dev/null)
}

# 1 周分の同期。QL がすぐ出るよう、まず各ノードのいちばん新しい run だけを送って ql_follow.py に
# 知らせ、そのあとで残り（過去の run と metrics）を送る。初回に過去の run が多くても QL は待たない
sync_once() {
  local entry node process_id rc=0
  local -A latest=() excludes=() listed=()
  local started=$SECONDS
  printf '[%s] %s へ同期%s\n' "$(date '+%F %T')" "$QL_DEST" "${QL_SINCE:+（$QL_SINCE 以降の run）}"

  for entry in $QL_NODES; do
    node=${entry%:*}
    process_id=${entry##*:}
    # 一覧に失敗したノードは、古い run まで送ってしまわないよう、この周は飛ばす
    list_runs "$node" "$process_id" || { rc=1; continue; }
    listed[$entry]=1
    latest[$entry]=$LATEST
    # 除外が無いときに printf '%q ' が空の引数 '' を作らないよう、あるときだけ書き出す
    excludes[$entry]=
    [ ${#EXCLUDES[@]} -eq 0 ] || excludes[$entry]=$(printf '%q ' "${EXCLUDES[@]}")
    rename_closed "$process_id"
    [ -n "$LATEST" ] || continue
    run_rsync "$process_id/$LATEST" "$node" "$SRC_ROOT/$process_id/$LATEST/" \
      "$QL_DEST/$process_id/$LATEST/" "${DATA_FILTER[@]}" || rc=1
  done
  [ "$DRY_RUN" -eq 1 ] || touch "$QL_DEST/.ql_sync_done"
  # QL はここまで終わると更新される。かかった時間を出し、遅れの元を見分けられるようにする
  printf '  いちばん新しい run まで %d 秒\n' $((SECONDS - started))

  for entry in $QL_NODES; do
    [ -n "${listed[$entry]:-}" ] || continue
    node=${entry%:*}
    process_id=${entry##*:}
    eval "EXCLUDES=(${excludes[$entry]})"
    # 除外は include より先に書く。rsync は最初に当てはまった規則を使う
    run_rsync "$process_id（残り）" "$node" "$SRC_ROOT/$process_id/" "$QL_DEST/$process_id/" \
      --exclude "/${latest[$entry]:-.}/" ${EXCLUDES[@]+"${EXCLUDES[@]}"} "${DATA_FILTER[@]}" || rc=1
    run_rsync "$process_id（データ以外）" "$node" "$SRC_ROOT/$process_id/" "$QL_DEST/$process_id/" \
      ${EXCLUDES[@]+"${EXCLUDES[@]}"} "${META_FILTER[@]}" || rc=1
  done

  printf '  1 周 %d 秒\n' $((SECONDS - started))
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
