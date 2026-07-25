#!/bin/bash
# 自动持续提交 TD/SOC 任务的脚本
# 依赖：
#   - 当前目录：orcacal
#   - ./find_missing_td_ids.sh  子目录参数，如 pic16inp
#   - ./list_td_jobs.sh              输出正在队列中的 TD 编号
#   - ./subtdsbatch.sh        按编号提交 TD 任务
#
# 逻辑与 auto_submit_opt.sh 相同，只是针对 TD 任务。

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# 用法: ./auto_submit_td.sh <子目录> [MAX_JOBS] [SLEEP_SEC]
if [[ $# -lt 1 ]]; then
  echo "用法: $0 <子目录> [最大并行数] [检查间隔秒]" >&2
  exit 1
fi

SUBDIR="$1"; shift
MAX_JOBS="${1:-20}"; [[ $# -gt 0 ]] && shift
SLEEP_SEC="${1:-1800}"

# 自动推断前缀（如 pic16inp -> pic16）
PREFIX="$(basename "$SUBDIR")"
PREFIX="${PREFIX%inp}"
PREFIX="${PREFIX%/}"

echo "自动提交 TD 脚本启动："
echo "  子目录: $SUBDIR"
echo "  前缀: $PREFIX"
echo "  最大并行任务数: $MAX_JOBS"
echo "  检查间隔: $SLEEP_SEC 秒"
echo "按 Ctrl+C 终止。"

TMP_BASE="/tmp/${PREFIX}_auto_submit_td_$$"  # $$ 是当前脚本的 PID
MISSING="${TMP_BASE}_missing.txt"
RUNNING="${TMP_BASE}_running.txt"
PENDING="${TMP_BASE}_pending.txt"
# 退出时自动清理临时文件（可选，不想要可以删掉这三行）
cleanup() {
  rm -f "$MISSING" "$RUNNING" "$PENDING" || true
}
trap cleanup EXIT
while true; do
  # 1. 计算 TD 缺失编号（基于 TD log）
  "$SCRIPT_DIR/find_missing_td_ids.sh" "$SUBDIR" "$PREFIX" > "$MISSING" || true
  if [[ ! -s "$MISSING" ]]; then
    echo "$(date '+%F %T') 没有缺失 TD 编号（或无法获取），休眠 $SLEEP_SEC 秒..."
    sleep "$SLEEP_SEC"
    continue
  fi

  # 2. 获取正在运行/排队中的 TD 编号（基于 squeue 作业名）
  "$SCRIPT_DIR/list_td_jobs.sh" "$PREFIX" > "$RUNNING" || true

  # TD 运行/排队计数
  if [[ -s "$RUNNING" ]]; then
    td_running_count=$(wc -l < "$RUNNING")
  else
    td_running_count=0
  fi

  # 当前用户整体 squeue 任务数（不限于 TD），用于全局并发上限
  total_queue_count=$(squeue -u "$USER" -h 2>/dev/null | wc -l | tr -d ' ')
  total_queue_count=${total_queue_count:-0}

  # 4. 计算待提交编号 = 缺失编号 - 正在运行编号
  if [[ -s "$RUNNING" ]]; then
    comm -23 "$MISSING" "$RUNNING" > "$PENDING"
  else
    cp "$MISSING" "$PENDING"
  fi

  if [[ ! -s "$PENDING" ]]; then
    echo "$(date '+%F %T') 没有待提交 TD 编号，当前总任务: $total_queue_count，TD运行: $td_running_count，休眠 $SLEEP_SEC 秒..."
    sleep "$SLEEP_SEC"
    continue
  fi

  # 5. 如当前总任务数小于上限，则补提交
  if (( total_queue_count >= MAX_JOBS )); then
    echo "$(date '+%F %T') 当前总任务数($total_queue_count) 已达到上限($MAX_JOBS)，TD运行: $td_running_count，休眠 $SLEEP_SEC 秒..."
    sleep "$SLEEP_SEC"
    continue
  fi

  # 剩余可提交的配额 = MAX_JOBS - 当前总任务数
  need=$(( MAX_JOBS - total_queue_count ))
  pending_count=$(wc -l < "$PENDING")

  if (( pending_count <= 0 )); then
    echo "$(date '+%F %T') pending_count=0，休眠 $SLEEP_SEC 秒..."
    sleep "$SLEEP_SEC"
    continue
  fi

  if (( need > pending_count )); then
    to_submit=$pending_count
  else
    to_submit=$need
  fi

  ids=$(head -n "$to_submit" "$PENDING" | tr '\n' ' ')

  echo "$(date '+%F %T') 总任务: $total_queue_count，TD运行: $td_running_count，待提交: $pending_count，将提交: $to_submit 个编号：$ids"

  # 6. 调用 subtdsbatch.sh 提交这些 TD 任务
  "$SCRIPT_DIR/subtdsbatch.sh" "$SUBDIR" $ids
  

  echo "$(date '+%F %T') 本轮 TD 提交完成，休眠 $SLEEP_SEC 秒..."
  sleep "$SLEEP_SEC"
done
