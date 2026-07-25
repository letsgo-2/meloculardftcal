#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
光敏剂特性评估工具 v2.5 (plot_summary2.py)
功能更新：
    1. 保持个体图与全汇总图不变。
    2. 新增：Generation 汇总图 (除去 pic0)。
    3. 新增：背景对比图 (pic0 灰色透明底 + Generation 彩色叠加)。
    4. 统一锁定：X[-0.1, 3.5], Y[0, 1.1], Color[0, 5], 增加 (0,1)-(1,0) 虚线。
    5. 新增：从 publishinpdftv22.csv 加载基准点，在个体图和 overlay 图中用红色五角星标记。
"""

import pandas as pd
import numpy as np
import argparse
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns

def load_pictos_labels(script_dir):
    paths_to_try = [
        os.path.join(script_dir, 'data', 'pictos.txt'),
        os.path.join(os.path.dirname(script_dir), 'data', 'pictos.txt'),
        'data/pictos.txt'
    ]
    mapping = {}
    for p in paths_to_try:
        if os.path.exists(p):
            with open(p, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'): continue
                    if ':' in line:
                        parts = line.split(':', 1)
                        mapping[parts[0].strip().lower()] = parts[1].strip().upper()
            break
    return mapping

def get_mapped_info(filename, mapping):
    basename = os.path.basename(filename)
    match = re.search(r'pic(\d+)', basename, re.IGNORECASE)
    if match:
        pic_id = match.group(1)
        pic_key = f"pic{pic_id}".lower()
        label = mapping.get(pic_key, pic_key.upper())
        return label, pic_id
    return os.path.splitext(basename)[0], "999"

def load_reference_points(script_dir):
    """从 publishinpdftv22.csv 加载基准点数据"""
    paths_to_try = [
        os.path.join(script_dir, 'data', 'publishinpdftv22.csv'),
        os.path.join(os.path.dirname(script_dir), 'data', 'publishinpdftv22.csv'),
        'data/publishinpdftv22.csv'
    ]
    
    for p in paths_to_try:
        if os.path.exists(p):
            try:
                df = pd.read_csv(p)
                ref_points = {}
                
                for _, row in df.iterrows():
                    index_val = str(row.get('index', '')).strip()
                    if not index_val or not index_val.startswith('pic'):
                        continue
                    
                    # 提取 pic 序号
                    match = re.search(r'pic(\d+)', index_val, re.IGNORECASE)
                    if not match:
                        continue
                    pic_id = match.group(1)
                    
                    # 跳过 pic0
                    if pic_id == '0':
                        continue
                    
                    s1 = row.get('S1_energy_eV', np.nan)
                    t1 = row.get('T1_energy_eV', np.nan)
                    
                    # 检查数据是否有效
                    try:
                        if pd.notna(s1) and pd.notna(t1):
                            s1_val = float(s1)
                            t1_val = float(t1)
                            
                            if s1_val > 0:  # 确保 S1 > 0，避免除零
                                st_gap = s1_val - t1_val
                                t1_s1_ratio = t1_val / s1_val
                                
                                ref_points[pic_id] = {
                                    'ST_Gap': st_gap,
                                    'T1_S1_Ratio': t1_s1_ratio,
                                    'S1': s1_val,
                                    'T1': t1_val
                                }
                                print(f"  📍 已加载基准点 pic{pic_id}: ST_Gap={st_gap:.4f}, T1/S1={t1_s1_ratio:.4f}")
                    except (ValueError, TypeError):
                        continue
                
                return ref_points
            except Exception as e:
                print(f"  ⚠️  加载基准点数据时出错: {e}")
                return {}
    
    print("  ⚠️  未找到 publishinpdftv22.csv 文件")
    return {}

def apply_plot_style(title):
    """统一的坐标轴、虚线和样式配置"""
    plt.xlim(-0.1, 3.5)
    plt.ylim(0.0, 1.1)
    plt.plot([0, 1], [1, 0], linestyle='--', color='gray', alpha=0.5, linewidth=1.2, zorder=1)
    plt.xlabel(r'$S_1 - T_1$ ($\Delta E_{ST}$, eV)', fontsize=13, fontweight='bold')
    plt.ylabel('$T_1 / S_1$ Ratio', fontsize=13, fontweight='bold')
    plt.title(title, fontsize=14, pad=15, fontweight='bold')
    plt.grid(True, linestyle=':', alpha=0.3, zorder=0)
    sns.despine()

def save_standard_scatter(df, path, title, ref_point=None):
    """标准彩色散点图生成逻辑，可选标记基准点"""
    if df.empty: return
    plt.figure(figsize=(10, 8), dpi=300)
    sns.set_style("ticks")
    
    scatter = plt.scatter(
        df['ST_Gap'], df['T1_S1_Ratio'], 
        c=df['S1'], cmap='viridis', vmin=0, vmax=5,
        s=55, alpha=0.8, edgecolor='white', linewidth=0.5, zorder=3
    )
    cbar = plt.colorbar(scatter)
    cbar.set_label('$S_1$ Energy (eV)', fontsize=11)
    
    # 标记基准点（红色五角星，黑色边缘）
    if ref_point is not None:
        plt.scatter(
            ref_point['ST_Gap'], ref_point['T1_S1_Ratio'],
            marker='*', s=100, color='red', edgecolor='black', linewidth=1.5,
            zorder=5, label='Reference Point'
        )
        plt.legend(loc='upper right', frameon=True, fontsize=10)
    
    apply_plot_style(title)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    print(f"  ✅ 已生成: {path}")

def save_overlay_scatter(gen_df, pic0_df, path, title, ref_points=None):
    """背景对比图：pic0 灰色透明 + Generation 彩色 + 基准点标记"""
    plt.figure(figsize=(10, 8), dpi=300)
    sns.set_style("ticks")
    
    # 1. 底色：pic0 (灰色透明)
    if not pic0_df.empty:
        plt.scatter(
            pic0_df['ST_Gap'], pic0_df['T1_S1_Ratio'], 
            color='lightgrey', s=40, alpha=0.3, label='RefDataset', zorder=2
        )
    
    # 2. 高亮：Generation (彩色)
    if not gen_df.empty:
        scatter = plt.scatter(
            gen_df['ST_Gap'], gen_df['T1_S1_Ratio'], 
            c=gen_df['S1'], cmap='viridis', vmin=0, vmax=5,
            s=60, alpha=0.8, edgecolor='white', linewidth=0.5, zorder=3, label='Generation'
        )
        cbar = plt.colorbar(scatter)
        cbar.set_label('$S_1$ Energy (eV)', fontsize=11)
    
    # 3. 标记所有基准点（红色五角星）
    if ref_points is not None and len(ref_points) > 0:
        ref_st_gaps = []
        ref_ratios = []
        for pic_id, ref_data in ref_points.items():
            ref_st_gaps.append(ref_data['ST_Gap'])
            ref_ratios.append(ref_data['T1_S1_Ratio'])
        
        if ref_st_gaps:
            plt.scatter(
                ref_st_gaps, ref_ratios,
                marker='*', s=400, color='red', edgecolor='black', linewidth=1.5,
                zorder=5, label='Reference Points'
            )
    
    apply_plot_style(title)
    plt.legend(loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    print(f"  ✅ 已生成对比图: {path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs='+')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = 'plots'
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    mapping = load_pictos_labels(script_dir)
    
    # 加载基准点数据
    print("\n--- Step 0: 正在加载基准点数据 ---")
    ref_points = load_reference_points(script_dir)
    
    all_data_list = []
    
    print("\n--- Step 1: 正在处理数据与个体图 ---")
    for f in args.input:
        if not os.path.exists(f): continue
        label, pic_id = get_mapped_info(f, mapping)
        df = pd.read_csv(f)
        
        s1_col = 'S1' if pic_id == '0' else 'S1_energy_eV'
        t1_col = 'T1' if pic_id == '0' else 'T1_energy_eV'
        
        if s1_col in df.columns and t1_col in df.columns:
            temp_df = df[[s1_col, t1_col]].copy()
            temp_df.columns = ['S1', 'T1']
            temp_df['ST_Gap'] = temp_df['S1'] - temp_df['T1']
            temp_df['T1_S1_Ratio'] = temp_df['T1'] / temp_df['S1']
            temp_df['pic_id'] = pic_id
            temp_df['Label'] = label
            temp_df = temp_df.replace([np.inf, -np.inf], np.nan).dropna()
            
            # 获取对应的基准点（仅对非 pic0 的数据）
            ref_point = None
            if pic_id != '0' and pic_id in ref_points:
                ref_point = ref_points[pic_id]
                print(f"  📍 将在个体图中标记基准点: pic{pic_id} ({label})")
            
            # 保存个体图（包含基准点标记）
            save_standard_scatter(
                temp_df, 
                os.path.join(output_dir, f"scatter_ST_{label}.png"), 
                f"Property Map: {label}",
                ref_point=ref_point
            )
            all_data_list.append(temp_df)

    if not all_data_list: return
    full_df = pd.concat(all_data_list, ignore_index=True)
    
    # 数据切片
    pic0_df = full_df[full_df['pic_id'] == '0']
    gen_df = full_df[full_df['pic_id'] != '0']

    print("\n--- Step 2: 正在生成汇总对比图 ---")
    
    # 1. 所有点的联合图 (含 pic0)
    save_standard_scatter(full_df, os.path.join(output_dir, "scatter_ST_combined_ALL.png"), "Combined Library (All Scaffolds)")
    
    # 2. 除去 pic0 的 Generation 联合图
    save_standard_scatter(gen_df, os.path.join(output_dir, "scatter_ST_generation_only.png"), "Generation Library (Excluding REF)")
    
    # 3. 背景对比图 (pic0 灰色透明 + Generation 彩色 + 基准点标记)
    save_overlay_scatter(
        gen_df, pic0_df, 
        os.path.join(output_dir, "scatter_ST_overlay_background.png"), 
        "Generation vs. Dataset",
        ref_points=ref_points
    )

    print(f"\n✨ 分析圆满完成！共有 {len(full_df)} 个分子进入统计。")

if __name__ == "__main__":
    main()