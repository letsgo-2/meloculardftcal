#!/bin/bash
# 自动持续提交 Opt/Freq 任务的脚本
# 依赖：
#   - 当前目录：orcacal
#   - ./find_missing_opt_ids.sh  子目录参数，如 pic16inp
#   - ./list_opt_jobs.sh           输出正在队列中的编号
#   - ./suboptsbatch.sh           按编号提交任务
#
# 逻辑：
#   1. 每次循环：
#      - 计算缺失的编号（根据 log）；
#      - 获取当前正在运行/排队的编号（根据 squeue 作业名）；
#      - 待提交编号 = 缺失编号 - 正在运行编号；
#      - 若正在运行的任务数 < MAX_JOBS，则从待提交编号中取若干，用 suboptsbatch.sh 提交，补到 MAX_JOBS；
#   2. 每隔 SLEEP_SEC 秒重复一次。

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# 用法: ./auto_submit_opt.sh <子目录> [MAX_JOBS] [SLEEP_SEC]
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

echo "自动提交脚本启动："
echo "  子目录: $SUBDIR"
echo "  前缀: $PREFIX"
echo "  最大并行任务数: $MAX_JOBS"
echo "  检查间隔: $SLEEP_SEC 秒"
echo "按 Ctrl+C 终止。"


TMP_BASE="/tmp/${PREFIX}_auto_submit_$$"  # $$ 是当前脚本的 PID
MISSING="${TMP_BASE}_missing.txt"
RUNNING="${TMP_BASE}_running.txt"
PENDING="${TMP_BASE}_pending.txt"
# 退出时自动清理临时文件（可选，不想要可以删掉这三行）
cleanup() {
  rm -f "$MISSING" "$RUNNING" "$PENDING" || true
}
trap cleanup EXIT

while true; do
  # 1. 计算缺失编号（基于 log）
  "$SCRIPT_DIR/find_missing_opt_ids.sh" "$SUBDIR" "$PREFIX" > "$MISSING" || true
  if [[ ! -s "$MISSING" ]]; then
    echo "$(date '+%F %T') 没有缺失编号（或无法获取），休眠 $SLEEP_SEC 秒..."
    sleep "$SLEEP_SEC"
    continue
  fi

  # 2. 获取正在运行/排队中的编号（基于 squeue 作业名）
  "$SCRIPT_DIR/list_opt_jobs.sh" "$PREFIX" > "$RUNNING" || true
  # 如果没有正在运行的 pic16 作业，文件可能为空，这种情况也允许

  # 3. 统计当前正在运行的任务数
  if [[ -s "$RUNNING" ]]; then
    running_count=$(wc -l < "$RUNNING")
  else
    running_count=0
  fi

    # 当前用户整体 squeue 任务数，用于全局并发上限
  total_queue_count=$(squeue -u "$USER" -h 2>/dev/null | wc -l | tr -d ' ')
  total_queue_count=${total_queue_count:-0}
  # 4. 计算待提交编号 = 缺失编号 - 正在运行编号
  #    注意：两个文件都应是已排序的编号列表
  if [[ -s "$RUNNING" ]]; then
    comm -23 "$MISSING" "$RUNNING" > "$PENDING"
  else
    cp "$MISSING" "$PENDING"
  fi

  if [[ ! -s "$PENDING" ]]; then
    echo "$(date '+%F %T') 没有待提交编号，opt任务运行: $running_count，当前总任务数: $total_queue_count，休眠 $SLEEP_SEC 秒..."
    sleep "$SLEEP_SEC"
    continue
  fi

  # 5. 如当前运行任务数小于上限，则补提交
  if (( total_queue_count >= MAX_JOBS )); then
    echo "$(date '+%F %T') 当前总任务数($total_queue_count) 已达到上限($MAX_JOBS)，opt任务运行: $running_count，休眠 $SLEEP_SEC 秒..."
    sleep "$SLEEP_SEC"
    continue
  fi

  # 需要补充的任务数
  need=$(( MAX_JOBS - total_queue_count ))

  # 待提交任务数
  pending_count=$(wc -l < "$PENDING")

  if (( pending_count <= 0 )); then
    echo "$(date '+%F %T') pending_count=0，休眠 $SLEEP_SEC 秒..."
    sleep "$SLEEP_SEC"
    continue
  fi

  # 实际要提交的任务数
  if (( need > pending_count )); then
    to_submit=$pending_count
  else
    to_submit=$need
  fi

  ids=$(head -n "$to_submit" "$PENDING" | tr '\n' ' ')

  echo "$(date '+%F %T') 当前运行opt任务: $running_count，待提交: $pending_count，将提交: $to_submit 个编号：$ids"

  # 6. 调用 suboptsbatch.sh 提交这些任务（传入子目录，前缀自动推断）
  "$SCRIPT_DIR/suboptsbatch.sh" "$SUBDIR" $ids

  echo "$(date '+%F %T') 本轮提交完成，休眠 $SLEEP_SEC 秒..."
  sleep "$SLEEP_SEC"
done
