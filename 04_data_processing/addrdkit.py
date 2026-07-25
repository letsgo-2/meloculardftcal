#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RDKit 药化指标计算脚本
=====================

功能说明：
    本脚本用于为包含 SMILES 字符串的 CSV 文件添加五个核心药化指标列：
    1. LogP - 脂水分配系数（脂溶性指标）
    2. TPSA - 拓扑极性表面积（渗透性指标）
    3. Rotatable_Bonds - 可旋转键数量（分子柔性指标）
    4. charge - 形式电荷 (Formal Charge)
    5. SA_Score - 合成可及性评分（合成难度指标）

使用方法：
    python addrdkit.py <输入文件> --output <输出文件>
    
    示例：
    python addrdkit.py data/pic25inpdft2.csv --output data/pic25inpdft3.csv
    python addrdkit.py data/pic3inpdft2.csv -o data/pic3inpdft3.csv

依赖要求：
    - pandas: 数据处理
    - rdkit: 分子描述符计算
    - sascorer.py: SA Score 计算（可选，位于脚本同目录）

输入要求：
    - CSV 文件必须包含 'SMILES' 列
    - SMILES 字符串应为有效的化学结构表示

输出说明：
    - 输出文件包含原文件的所有列
    - 新增五列：LogP, TPSA, Rotatable_Bonds, charge, SA_Score
    - 如果 SMILES 解析失败或缺少 sascorer.py，对应值显示为 "N/A"

作者：Auto
日期：2024
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os
import argparse

# 尝试导入 SA Score 计算脚本
try:
    import sascorer
except ImportError:
    print("【警告】未在当前目录找到 sascorer.py。SA_Score 将显示为 N/A。")
    print("请从 https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py 下载。")
    sascorer = None

def calculate_rdkit_metrics(smiles, missing_cols):
    """
    根据 SMILES 字符串计算五个核心药化指标 (仅计算缺失列)
    
    参数:
        smiles (str): 分子的 SMILES 字符串表示
        missing_cols (list): 需要计算的列名列表
        
    返回:
        dict: 包含计算结果的字典
    """
    res = {col: "N/A" for col in missing_cols}
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # 1. LogP
            if "LogP" in missing_cols:
                res["LogP"] = round(Descriptors.MolLogP(mol), 3)
            # 2. TPSA
            if "TPSA" in missing_cols:
                res["TPSA"] = round(Descriptors.TPSA(mol), 3)
            # 3. Rotatable Bonds
            if "Rotatable_Bonds" in missing_cols:
                res["Rotatable_Bonds"] = Descriptors.NumRotatableBonds(mol)
            # 4. charge
            if "charge" in missing_cols:
                res["charge"] = Chem.GetFormalCharge(mol)
            # 5. SA Score
            if "SA_Score" in missing_cols and sascorer:
                res["SA_Score"] = round(sascorer.calculateScore(mol), 3)
    except Exception as e:
        print(f"SMILES 解析错误: {smiles} -> {e}")
    return res

def main():
    """
    主函数：解析命令行参数，读取 CSV 文件，计算药化指标并保存结果
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description='为 CSV 文件添加 RDKit 药化指标列',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  %(prog)s data/pic25inpdft2.csv --output data/pic25inpdft3.csv
  %(prog)s data/pic3inpdft2.csv -o data/pic3inpdft3.csv

注意事项:
  - 输入文件必须包含 'SMILES' 列
  - 如果未找到 sascorer.py，SA_Score 列将显示为 N/A
  - 每处理 50 条记录会显示进度信息
        """
    )
    parser.add_argument('input', type=str, help='输入 CSV 文件路径（例如: data/pic3inpdft2.csv）')
    parser.add_argument('--output', '-o', type=str, required=True, 
                        help='输出 CSV 文件路径（例如: data/pic3inpdft3.csv）')
    
    args = parser.parse_args()
    
    INPUT_CSV = args.input
    OUTPUT_CSV = args.output
    
    if not os.path.exists(INPUT_CSV):
        print(f"错误: 找不到输入文件 {INPUT_CSV}")
        return

    print(f"正在读取主表: {INPUT_CSV} ...")
    df = pd.read_csv(INPUT_CSV)

    if 'SMILES' not in df.columns:
        print("错误: CSV 文件中未找到 'SMILES' 列。")
        return

    # --- 核心改进：检测哪些列需要计算 ---
    target_cols = ["LogP", "TPSA", "Rotatable_Bonds", "charge", "SA_Score"]
    missing_cols = [col for col in target_cols if col not in df.columns]

    if not missing_cols:
        print("【提示】所有药化指标列已存在，无需计算。直接复制文件到输出路径。")
        df.to_csv(OUTPUT_CSV, index=False)
        return

    print(f"检测到缺失列: {missing_cols}，开始批量计算...")
    
    # 批量计算
    metrics_list = []
    for idx, row in df.iterrows():
        metrics = calculate_rdkit_metrics(row['SMILES'], missing_cols)
        metrics_list.append(metrics)
        
        if (idx + 1) % 50 == 0:
            print(f"进度: {idx+1}/{len(df)} 已完成")

    # 转换为 DataFrame 并与原表拼接
    df_metrics = pd.DataFrame(metrics_list)
    df_final = pd.concat([df, df_metrics], axis=1)

    # 保存结果
    df_final.to_csv(OUTPUT_CSV, index=False)
    print(f"\n【处理成功】补全后的表格已保存至: {OUTPUT_CSV}")
    print(f"补全的列: {', '.join(missing_cols)}")

if __name__ == "__main__":
    main()