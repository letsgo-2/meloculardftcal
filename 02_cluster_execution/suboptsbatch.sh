#!/usr/bin/env bash
set -euo pipefail

# 用法示例：
#   ./suboptsbatch.sh pic16inp 20-30
#   ./suboptsbatch.sh pic13inp 1-10
#   ./suboptsbatch.sh pic16inp "[20-30],[50-60],78"
#
# 第一个参数：脚本所在目录（必填，如 pic16inp 或 pic13inp）
# 前缀自动从目录名推断：去掉末尾的 "inp" 得到 pic16、pic13 等

if [[ $# -lt 2 ]]; then
    echo "用法: $0 <脚本目录> 序号/区间 [...]" >&2
    echo "示例: $0 pic16inp 1-10" >&2
    echo "     $0 pic13inp 20-30 50-60 78" >&2
    echo "     $0 pic16inp \"[20-30],[50-60],78\"" >&2
    exit 1
fi

SCRIPTDIR="$1"; shift
PREFIX="$(basename "$SCRIPTDIR")"
PREFIX="${PREFIX%inp}"

if [[ $# -lt 1 ]]; then
    echo "错误: 需要至少一个序号或区间参数" >&2
    exit 1
fi

if [[ ! -d "$SCRIPTDIR" ]]; then
    echo "错误: 找不到目录: $SCRIPTDIR" >&2
    exit 1
fi

cd "$SCRIPTDIR"

if [[ $# -lt 1 ]]; then
    echo "用法: PREFIX=<前缀> $0 序号/区间 [...]" >&2
    echo "示例: PREFIX=pic16 $0 1-10" >&2
    echo "     PREFIX=pic13 $0 20-30 50-60 78" >&2
    echo "     PREFIX=pic16 $0 \"[20-30],[50-60],78\"" >&2
    exit 1
fi


#sleep 21600


# 把所有参数合并成一个字符串，方便同时支持空格和逗号分隔
spec="$*"
# 逗号改成空格，这样  "[20-30],[50-60],78" → "[20-30] [50-60] 78"
spec=${spec//,/ }

submit_one() {
    local idx="$1"
    # 先转十进制，避免前导 0 被当作八进制导致 0008 等报错
    local idx_dec=$((10#$idx))
    # 补零成四位，例如 1 -> 0001
    local idx_padded
    idx_padded=$(printf "%04d" "$idx_dec")
    local script="${PREFIX}_${idx_padded}.sh"

    if [[ ! -f "$script" ]]; then
        echo "警告: 找不到脚本 $script，跳过。" >&2
        return
    fi

    echo "提交: sbatch $script"
    sbatch "$script"
}

# 逐个 token 解析
for token in $spec; do
    # 去掉可选的方括号前后缀，例如 [20-30] -> 20-30
    token="${token#[}"
    token="${token%]}"

    # 如果是区间 a-b
    if [[ "$token" == *-* ]]; then
        start=${token%-*}
        end=${token#*-}

        # 确保是整数
        if ! [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]]; then
            echo "跳过非法区间: $token" >&2
            continue
        fi

        # 转十进制，避免 0008 被当作八进制
        start=$((10#$start))
        end=$((10#$end))

        if (( start > end )); then
            echo "跳过区间(起始>结束): $start-$end" >&2
            continue
        fi

        for ((i=start; i<=end; i++)); do
            submit_one "$i"
        done
    else
        # 单个数字，例如 78
        if ! [[ "$token" =~ ^[0-9]+$ ]]; then
            echo "跳过非法序号: $token" >&2
            continue
        fi
        submit_one "$token"
    fi
done
