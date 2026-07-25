#!/bin/bash
# 用法（在 orcacal 目录下）:
#   bash find_missing_td_ids.sh pic16inp
#
# 功能：
#   在给定子目录中查找 pic16_*td.log，筛选出包含
#   "ORCA TERMINATED NORMALLY" 的 TD 日志，提取四位编号。
#   仅对“基础态主日志已正常结束”的编号做 TD 缺失检测：
#   TD 待跑编号 = 基础态成功编号 - TD 成功编号。

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

BASE_OK_IDS=$(mktemp "/tmp/${PREFIX}_base_ok_ids.XXXXXX")
TD_OK_IDS=$(mktemp "/tmp/${PREFIX}_td_ok_ids.XXXXXX")
cleanup() {
  rm -f "$BASE_OK_IDS" "$TD_OK_IDS"
}
trap cleanup EXIT

# 基础态成功编号（主日志）
grep -l 'ORCA TERMINATED NORMALLY' "${PREFIX}"_*.log 2>/dev/null \
  | grep -v 'td.log' \
  | sed -E "s/.*${PREFIX}_([0-9]{4})\.log/\1/" \
  | sort -n > "$BASE_OK_IDS"

# TD 成功编号
grep -l 'ORCA TERMINATED NORMALLY' "${PREFIX}"_*td.log 2>/dev/null \
  | sed -E "s/.*${PREFIX}_([0-9]{4})td\.log/\1/" \
  | sort -n > "$TD_OK_IDS"

if [ ! -s "$BASE_OK_IDS" ]; then
  # 基础态没有成功的编号，则无需跑 TD
  exit 0
fi

# 待跑 TD 编号 = 基础态成功编号 - TD 成功编号
if [ -s "$TD_OK_IDS" ]; then
  comm -23 "$BASE_OK_IDS" "$TD_OK_IDS"
else
  cat "$BASE_OK_IDS"
fi
