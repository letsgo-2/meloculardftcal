#!/bin/bash
# 用法（在 orcacal 目录下）:
#   bash find_missing_opt_ids.sh <子目录> [前缀]
#
# 功能：
#   在给定子目录中查找 pic16_*.sh 文件（排除 *td.sh），
#   提取编号作为总列表；然后查找包含 "ORCA TERMINATED NORMALLY"
#   且不以 td.log 结尾的主日志，提取编号作为已完成列表，
#   计算缺失的编号并输出。

SUBDIR="$1"
PREFIX="${2:-}"

if [ -z "$SUBDIR" ]; then
  echo "用法: $0 子目录名 [前缀]" >&2
  exit 1
fi

if [ ! -d "$SUBDIR" ]; then
  exit 1
fi

# 自动推断前缀（如 pic16inp -> pic16）
if [ -z "$PREFIX" ]; then
  PREFIX="$(basename "$SUBDIR")"
  PREFIX="${PREFIX%inp}"
  PREFIX="${PREFIX%/}"
fi

cd "$SUBDIR" || { echo "无法进入目录: $SUBDIR"; exit 1; }

ALL_IDS=$(mktemp "/tmp/${PREFIX}_all_ids.XXXXXX")
OK_IDS=$(mktemp "/tmp/${PREFIX}_ok_ids.XXXXXX")
cleanup() {
  rm -f "$ALL_IDS" "$OK_IDS"
}
trap cleanup EXIT

# 1. 提取所有 ${PREFIX}_*.sh 文件的编号（排除 *td.sh）作为总列表
ls "${PREFIX}"_*.sh 2>/dev/null \
  | grep -v 'td.sh' \
  | sed -E "s/.*${PREFIX}_([0-9]{4})\\.sh/\\1/" \
  | sort -n > "$ALL_IDS"

if [ ! -s "$ALL_IDS" ]; then
  # echo "未在目录 $SUBDIR 中找到 pic16_*.sh 文件"
  exit 0
fi

# 2. 提取包含 ORCA TERMINATED NORMALLY 且不是 *td.log 的主日志编号
grep -l 'ORCA TERMINATED NORMALLY' "${PREFIX}"_*.log 2>/dev/null \
  | grep -v 'td.log' \
  | sed -E "s/.*${PREFIX}_([0-9]{4})\\.log/\\1/" \
  | sort -n > "$OK_IDS"

# 3. 计算缺失编号（总列表 - 已完成列表）
# echo "目录: $PWD"
# echo "缺失的编号为:"
comm -23 "$ALL_IDS" "$OK_IDS"
