#!/bin/bash
# 用法（在 orcacal 目录下）:
#   ./list_td_jobs.sh <前缀> [用户名]
#   例如:
#     ./list_td_jobs.sh pic16
#     ./list_td_jobs.sh pic16 myuser
#
# 功能：
#   使用 squeue 提取当前用户在队列中的作业名称（NAME 列），
#   过滤“TD 任务”（作业名以 td 结尾的 pic 作业，例如 pic16_0001td 或 pic16_0001_td），
#   输出对应的四位编号列表（去重排序），便于与 TD 缺失编号比对。

PREFIX="$1"
USER_NAME="${2:-$USER}"

if [ -z "$PREFIX" ]; then
  echo "用法: $0 <前缀> [用户名]" >&2
  exit 1
fi

if ! command -v squeue >/dev/null 2>&1; then
  echo "错误: 未找到 squeue 命令，请确认在 Slurm 环境中运行。"
  exit 1
fi

NAMES_FILE=$(mktemp "/tmp/${PREFIX}_running_td_jobs.XXXXXX")
cleanup() {
  rm -f "$NAMES_FILE"
}
trap cleanup EXIT

squeue -u "$USER_NAME" -o "%j" 2>/dev/null \
  | sed '1d' \
  | grep -E "^${PREFIX}_[0-9]+.*td$" \
  | sort -u > "$NAMES_FILE"

if [ ! -s "$NAMES_FILE" ]; then
  exit 0
fi

# 提取编号：兼容 pic16_0001td 和 pic16_0001_td 两种命名
sed -n -E "s/.*${PREFIX}_([0-9]{4})_?td$/\\1/p" "$NAMES_FILE" | sort -u
