#!/bin/bash
# 用法（在 orcacal 目录下）:
#   ./list_opt_jobs.sh [前缀] [用户名]
#   示例: ./list_opt_jobs.sh pic16
#         ./list_opt_jobs.sh pic13 myuser
#
# 功能：
#   使用 squeue 提取当前用户在队列中的作业名称（NAME 列），
#   过滤名称与 pic 计算相关的作业（如 pic16_XXXX 形式），
#   输出唯一且排序后的作业名称和对应编号，方便与 log 编号比对。

PREFIX="${1:-pic16}"
USER_NAME="${2:-$USER}"

# 检查 squeue 是否可用
if ! command -v squeue >/dev/null 2>&1; then
  echo "错误: 未找到 squeue 命令，请确认在 Slurm 环境中运行。"
  exit 1
fi

# echo "当前用户: $USER_NAME"
# echo "队列中与 pic 相关的作业名称 (去重排序)："

NAMES_FILE=$(mktemp "/tmp/${PREFIX}_running_jobs.XXXXXX")
cleanup() {
  rm -f "$NAMES_FILE"
}
trap cleanup EXIT

# 提取作业名称列，过滤以 ${PREFIX}_ 开头且后接数字的 name，排除 td 结尾的任务，并去重排序
squeue -u "$USER_NAME" -o "%j" 2>/dev/null \
  | sed '1d' \
  | grep -E "^${PREFIX}_[0-9]+" \
  | grep -vE '(td|_td)$' \
  | sort -u > "$NAMES_FILE"

if [ ! -s "$NAMES_FILE" ]; then
  # echo "没有找到名称以 pic 开头的作业。"
  exit 0
fi

# cat /tmp/pic_jobs_names.txt

# echo
#echo "如果作业名为 ${PREFIX}_XXXX 形式，则对应的编号为（四位数字部分）："
sed -n -E "s/.*${PREFIX}_([0-9]{4}).*/\\1/p" "$NAMES_FILE" | sort -u
