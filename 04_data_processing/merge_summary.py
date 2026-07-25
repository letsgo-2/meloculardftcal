#!/usr/bin/env python3
"""
将 data/pic16.csv 添加编号列（pic16_0001 ...），并把 data/summary.csv 里的新字段合并到末尾，输出 pic16index.csv。

示例用法：
    # 使用默认路径，输出到 data/pic16index.csv
    python merge_pic16_summary.py

    # 自定义输出路径
    python merge_pic16_summary.py --output data/pic16index_merged.csv

    # 保留 summary.csv 的 SMILES 列（默认会移除以避免重复）
    python merge_pic16_summary.py --keep-summary-smiles
"""

import argparse
from pathlib import Path
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(
        description="合并 原始csv文件（例如：pic16.csv） 与 生成的计算统计文件 summary.csv，生成含编号的 picxxindex.csv"
    )
    parser.add_argument(
        "--piccsvfile",
        default="",
        required=True,
        help="原始 pic16.csv 路径（例如: data/pic16.csv）",
    )
    parser.add_argument(
        "--summary",
        default="",
        required=True,
        help="生成的计算统计文件picxxsummary.csv 路径（例如: data/pic16inpout_summary.csv），需包含 filename 列",
    )
    parser.add_argument(
        "--output",
        default="",
        required=True,
        help="输出文件路径（默认: 基于 summary 文件名生成，例如 summary.csv -> summary_index.csv）",
    )
    parser.add_argument(
        "--keep-summary-smiles",
        action="store_true",
        help="保留 summary.csv 中的 SMILES 列（默认删除以避免重复列）",
    )

    args = parser.parse_args()

    pic16_path = Path(args.piccsvfile)
    summary_path = Path(args.summary)
    
    # 如果 output 未指定，基于 summary 文件名生成
    if args.output:
        output_path = Path(args.output)
    else:
        # 使用 summary 文件的目录和文件名（stem）加上 _index.csv
        output_path = summary_path.parent / f"{summary_path.stem}_index.csv"

    if not pic16_path.exists():
        raise FileNotFoundError(f"未找到 原始csv文件: {pic16_path}")
    if not summary_path.exists():
        raise FileNotFoundError(f"未找到 summary 文件: {summary_path}")

    # 读取 pic16；如果已有 ID 列则复用，否则按顺序生成
    df = pd.read_csv(pic16_path)
    if "ID" not in df.columns:
        df.insert(0, "ID", [f"{pic16_path.stem}_{i:04d}" for i in range(1, len(df) + 1)])

    # 读取 summary
    df_sum = pd.read_csv(summary_path)
    if "filename" not in df_sum.columns:
        raise ValueError("summary.csv 中缺少 filename 列，无法对应合并")

    # 重命名 filename -> ID
    df_sum = df_sum.rename(columns={"filename": "ID"})

    # 默认去掉 summary 的 SMILES 列，避免与 pic16 中重复
    if not args.keep_summary_smiles and "SMILES" in df_sum.columns:
        df_sum = df_sum.drop(columns=["SMILES"])

    # 合并：保持 pic16 行顺序，summary 的新增列附加在末尾
    existing_cols = set(df.columns)
    new_cols = [c for c in df_sum.columns if c not in existing_cols]
    df_out = df.merge(df_sum, on="ID", how="left")
    # 调整列顺序：原列在前，再追加 summary 新列
    df_out = df_out[list(df.columns) + new_cols]

    # 保存
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(output_path, index=False)
    print(f"完成，写入 {output_path}，共 {len(df_out)} 行，列数 {df_out.shape[1]}")


if __name__ == "__main__":
    main()

