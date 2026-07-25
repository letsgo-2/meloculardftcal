#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
光敏剂骨架属性分布综合评估工具 (plot_summary.py v22.0)
===================================================
更新说明：
    1. 边际分布归一化：强制使用 common_norm=False，确保不同样本量下的分布形状伸缩一致。
    2. 字段映射：pic0 (S1) vs 其他 pic (S1_energy_eV) 自动适配。
    3. 坐标锁定：MW (0-1200), S1 (1.0-5.0 eV), SA (0-8), 波长 (200-1000 nm)。
    4. 增加图片输出路径的实时 print 提示。

用法示例（在仓库 plot/ 目录下执行，路径按实际调整）：
    # 多个 pic 输入 CSV（文件名需含 pic 数字，如 pic0、pic3inpout_m.csv），生成个体图 + plots/ 下汇总 KDE/散点
    python3 plot_summary.py ../data/xhx/pic0.csv ../data/xhx/pic3inpout_m.csv

    # 通配一批文件
    python3 plot_summary.py ../data/xhx/pic*.csv

    # 指定 SMILES 列名（若表头不是默认探测顺序时使用）
    python3 plot_summary.py ../data/xhx/pic3.csv --col SMILES

    # 指定图 1（三属性 KDE 对比）输出文件名（仍保存在 plots/ 下）
    python3 plot_summary.py ../data/xhx/pic0.csv ../data/xhx/pic9.csv --output my_scaffold_kde.png

依赖：需能读取 pictos.txt、publishinpdftv22.csv（与脚本相对路径见代码内 paths_to_try），输出目录为当前工作目录下的 plots/。
"""

import pandas as pd
import numpy as np
import argparse
import os
import re
import warnings
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from scipy import stats
from scipy.interpolate import interp1d

# 环境配置：适配服务器环境
warnings.filterwarnings("ignore")
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_pictos_labels(script_dir):
    """加载映射关系 picX:label"""
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
                        mapping[parts[0].strip().lower()] = parts[1].strip()
            break
    return mapping

def format_final_label(raw_label):
    """标准化标签显示"""
    if not raw_label: return "Unknown"
    if raw_label.lower() == 'other': return 'RefDataset'
    return raw_label.upper() if raw_label.lower().startswith('s') else raw_label

def load_reference_points(script_dir, pictos_labels):
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
                    
                    # 获取对应的标签
                    pic_key = f"pic{pic_id}".lower()
                    raw_label = pictos_labels.get(pic_key, None)
                    if not raw_label:
                        continue
                    label = format_final_label(raw_label)
                    
                    # 提取属性值
                    sa_score = row.get('SA_Score', np.nan)
                    logp = row.get('LogP', np.nan)
                    s1_energy = row.get('S1_energy_eV', np.nan)
                    
                    # 计算 Molecular_Weight（如果缺失）
                    molecular_weight = np.nan
                    smiles = row.get('SMILES', '')
                    if pd.notna(smiles) and smiles:
                        try:
                            mol = Chem.MolFromSmiles(str(smiles))
                            if mol:
                                molecular_weight = Descriptors.MolWt(mol)
                        except:
                            pass
                    
                    # 计算 S1_wavelength_nm（从 S1_energy_eV）
                    s1_wavelength = np.nan
                    if pd.notna(s1_energy) and s1_energy > 0:
                        s1_wavelength = 1240.0 / float(s1_energy)
                    
                    # 构建基准点数据
                    point_data = {
                        'label': label,
                        'pic_id': pic_id,
                        'SA_Score': float(sa_score) if pd.notna(sa_score) else None,
                        'LogP': float(logp) if pd.notna(logp) else None,
                        'Molecular_Weight': float(molecular_weight) if pd.notna(molecular_weight) else None,
                        'S1': float(s1_energy) if pd.notna(s1_energy) else None,
                        'S1_wavelength_nm': float(s1_wavelength) if pd.notna(s1_wavelength) else None
                    }
                    
                    # 只保存有有效数据的点
                    if any(v is not None for k, v in point_data.items() if k in ['SA_Score', 'LogP', 'Molecular_Weight', 'S1', 'S1_wavelength_nm']):
                        ref_points[pic_id] = point_data
                        sa_str = f"{sa_score:.2f}" if pd.notna(sa_score) else "N/A"
                        logp_str = f"{logp:.2f}" if pd.notna(logp) else "N/A"
                        mw_str = f"{molecular_weight:.1f}" if pd.notna(molecular_weight) else "N/A"
                        s1_str = f"{s1_energy:.3f}" if pd.notna(s1_energy) else "N/A"
                        wave_str = f"{s1_wavelength:.1f}" if pd.notna(s1_wavelength) else "N/A"
                        print(f"  📍 已加载基准点 pic{pic_id} ({label}): SA={sa_str}, LogP={logp_str}, MW={mw_str}, S1={s1_str}, Wave={wave_str}")
                
                return ref_points
            except Exception as e:
                print(f"  ⚠️  加载基准点数据时出错: {e}")
                return {}
    
    print("  ⚠️  未找到 publishinpdftv22.csv 文件")
    return {}

def get_mapped_label(filename, pictos_labels):
    """解析文件名并映射标签"""
    basename = os.path.basename(filename)
    match = re.search(r'pic(\d+)', basename, re.IGNORECASE)
    if match:
        pic_id = match.group(1)
        pic_key = f"pic{pic_id}".lower()
        raw_mapped = pictos_labels.get(pic_key, pic_key)
        return format_final_label(raw_mapped), pic_id
    return basename, "999"

def process_properties(df, smiles_col):
    """属性补全及 S1 -> 波长转换"""
    prop_tasks = {'Molecular_Weight': Descriptors.MolWt, 'LogP': Descriptors.MolLogP}
    missing = [c for c in prop_tasks if c not in df.columns]
    
    mols = []
    if missing or 'SA_Score' not in df.columns:
        mols = [Chem.MolFromSmiles(s) if isinstance(s, str) else None for s in df[smiles_col]]
        
    for col in missing:
        df[col] = [prop_tasks[col](m) if m else np.nan for m in mols]
        
    if 'SA_Score' not in df.columns:
        try:
            from rdkit.Contrib.SA_Score import sascorer
            df['SA_Score'] = [sascorer.calculateScore(m) if m else np.nan for m in mols]
        except: df['SA_Score'] = np.nan

    if 'S1' in df.columns:
        df['S1_wavelength_nm'] = df['S1'].apply(lambda x: 1240.0 / x if (pd.notnull(x) and x > 0) else np.nan)
    return df

def plot_joint_scatter(df, x_col, y_col, title, output_path, x_label, y_label, x_lim, y_lim, color='#3498db', hue_col=None, palette=None, ref_point=None, ref_points_list=None):
    """带归一化一致性边际分布的散点图，可选标记参考点"""
    try:
        plot_df = df.dropna(subset=[x_col, y_col]).copy()
        if plot_df.empty: return

        # 设置绘图风格
        sns.set_style("white")
        # 创建 JointGrid 确保比例一致 (ratio=5)
        g = sns.JointGrid(data=plot_df, x=x_col, y=y_col, hue=hue_col, palette=palette, height=10, ratio=5)
        
        if hue_col:
            # 核心：使用 common_norm=False 保证小样本组的分布曲线不会被压扁，实现伸缩一致
            g.plot_joint(sns.scatterplot, s=50, alpha=0.6, edgecolor='white', linewidth=0.4)
            g.plot_marginals(sns.kdeplot, fill=True, alpha=0.3, common_norm=False)
        else:
            g.plot_joint(sns.scatterplot, s=50, alpha=0.6, color=color, edgecolor='white', linewidth=0.4)
            g.plot_marginals(sns.kdeplot, color=color, fill=True, alpha=0.3)
        
        # 强制锁定坐标轴（先设置范围，确保参考点在可见范围内）
        g.ax_joint.set_xlim(x_lim)
        g.ax_joint.set_ylim(y_lim)
        g.ax_joint.set_xlabel(x_label, fontsize=14, fontweight='bold')
        g.ax_joint.set_ylabel(y_label, fontsize=14, fontweight='bold')
        g.ax_joint.grid(linestyle='--', alpha=0.5, zorder=0)
        
        # 标记参考点（单个或列表）- 在设置坐标轴范围之后标记
        ref_point_added = False
        if ref_point is not None:
            x_val = ref_point.get(x_col, None)
            y_val = ref_point.get(y_col, None)
            if x_val is not None and y_val is not None:
                # 检查参考点是否在可见范围内
                if x_lim[0] <= x_val <= x_lim[1] and y_lim[0] <= y_val <= y_lim[1]:
                    point_color = ref_point.get('color', 'red')
                    g.ax_joint.scatter(x_val, y_val, marker='*', s=300, color=point_color, edgecolors='black', linewidth=1.5, zorder=10, label='Reference Point')
                    ref_point_added = True
        
        if ref_points_list is not None and len(ref_points_list) > 0:
            first_ref_point = True
            for ref_pt in ref_points_list:
                x_val = ref_pt.get(x_col, None)
                y_val = ref_pt.get(y_col, None)
                if x_val is not None and y_val is not None:
                    # 检查参考点是否在可见范围内
                    if x_lim[0] <= x_val <= x_lim[1] and y_lim[0] <= y_val <= y_lim[1]:
                        point_color = ref_pt.get('color', 'red')
                        # 只为第一个参考点添加 label，用于图例
                        label = 'Reference Points' if first_ref_point else ''
                        g.ax_joint.scatter(x_val, y_val, marker='*', s=300, color=point_color, edgecolors='black', linewidth=1.5, zorder=10, label=label)
                        first_ref_point = False
                        ref_point_added = True
        
        # 如果有参考点，确保图例显示
        if ref_point_added:
            handles, labels = g.ax_joint.get_legend_handles_labels()
            if handles:
                g.ax_joint.legend(handles=handles, labels=labels, loc='upper right', fontsize=9, framealpha=0.9)
        
        plt.subplots_adjust(top=0.93)
        g.fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✅ 散点边际图已保存: {output_path}")
    except Exception as e:
        print(f"  ❌ 绘图失败 ({title}): {e}")

def main():
    parser = argparse.ArgumentParser(
        description="多骨架 CSV：SA/MW/LogP KDE 与 SA–波长、MW–S1 联合散点图（pic0=RefDataset，其余=Generation）。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
用法示例：
  %(prog)s ../data/xhx/pic0.csv ../data/xhx/pic3inpout_m.csv
  %(prog)s ../data/xhx/pic*.csv --col SMILES
  %(prog)s pic0.csv pic9.csv -o my_kde.png

说明：
  - 位置参数 input：一个或多个 CSV；文件名中需能匹配 pic<number> 以映射骨架标签。
  - pic0 使用列 S1；非 pic0 优先 S1_energy_eV（见脚本内字段映射）。
  - 图 1 文件名由 -o/--output 指定，默认 scaffold_comparison_kde.png；其余固定名见 plots/。
""",
    )
    parser.add_argument(
        "input",
        nargs="+",
        help="输入 CSV 路径（可多个）",
    )
    parser.add_argument(
        "--col",
        default="SMILES",
        help="SMILES 列名（默认 SMILES；不存在时会尝试 SMILES/smiles）",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="scaffold_comparison_kde.png",
        help="图 1（三列 KDE）保存文件名，仅文件名，目录固定为 plots/（默认 scaffold_comparison_kde.png）",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    pictos_labels = load_pictos_labels(script_dir)
    os.makedirs('plots', exist_ok=True)
    
    # 加载基准点数据（提前加载，用于后续散点图标记）
    print("\n--- Step 0: 正在加载基准点数据 ---")
    ref_points = load_reference_points(script_dir, pictos_labels)

    all_data_list = []
    print("\n--- Step 1: 读取数据与智能字段解析 ---")
    
    for f in args.input:
        if not os.path.exists(f): continue
        df_temp = pd.read_csv(f)
        label, pic_id = get_mapped_label(f, pictos_labels)
        
        # 不同 ID 的能量字段映射逻辑
        energy_src = 'S1' if pic_id == '0' else 'S1_energy_eV'
        if energy_src in df_temp.columns:
            df_temp['S1'] = df_temp[energy_src]
        else:
            # 自动降级寻找
            for alt in ['S1', 'S1_energy_eV', 'energy']:
                if alt in df_temp.columns:
                    df_temp['S1'] = df_temp[alt]
                    break
        
        smiles_col = next((c for c in [args.col, 'SMILES', 'smiles'] if c in df_temp.columns), None)
        if smiles_col:
            df_temp = process_properties(df_temp, smiles_col)
        
        df_temp['Scaffold_Label'] = label
        df_temp['pic_id'] = pic_id
        df_temp['Category'] = 'RefDataset' if pic_id == '0' else 'Generation'
        all_data_list.append(df_temp)
        
        # 生成个体散点图 (锁定范围)
        idx_color = sns.color_palette("husl", 11)[int(pic_id)%11 if pic_id.isdigit() else 0]
        
        # 获取对应的参考点（如果存在）- 使用红色标记
        ref_point_sa_wave = None
        ref_point_mw_s1 = None
        if pic_id in ref_points:
            ref_pt = ref_points[pic_id].copy()
            ref_pt['color'] = 'red'  # 所有散点图的参考点都使用红色
            # SA-Wave 参考点
            if ref_pt.get('SA_Score') is not None and ref_pt.get('S1_wavelength_nm') is not None:
                ref_point_sa_wave = {k: ref_pt.get(k) for k in ['SA_Score', 'S1_wavelength_nm', 'color']}
            # MW-S1 参考点
            if ref_pt.get('Molecular_Weight') is not None and ref_pt.get('S1') is not None:
                ref_point_mw_s1 = {k: ref_pt.get(k) for k in ['Molecular_Weight', 'S1', 'color']}
        
        # SA-Wave
        plot_joint_scatter(df_temp, 'SA_Score', 'S1_wavelength_nm', f"SA vs Wave: {label}", 
                           os.path.join('plots', f"individual_scatter_sa_wave_{label}.png"),
                           'SA Score', 'S1 Wavelength (nm)', (0, 8), (200, 1000), color=idx_color,
                           ref_point=ref_point_sa_wave)
        # MW-S1
        plot_joint_scatter(df_temp, 'Molecular_Weight', 'S1', f"MW vs S1: {label}", 
                           os.path.join('plots', f"individual_scatter_mw_s1_{label}.png"),
                           'Molecular Weight', 'S1 Energy (eV)', (0, 1600), (1.0, 5.0), color=idx_color,
                           ref_point=ref_point_mw_s1)

    if not all_data_list: return
    full_df = pd.concat(all_data_list, ignore_index=True)
    if 'logp' in full_df.columns: full_df.rename(columns={'logp': 'LogP'}, inplace=True)
    gen_df = full_df[full_df['Category'] == 'Generation']
    other_df = full_df[full_df['Category'] == 'RefDataset']

    # --- KDE 汇总图部分 ---
    print("\n--- Step 2: 生成汇总 KDE 曲线对比图 ---")
    props = ['SA_Score', 'Molecular_Weight', 'LogP']
    
    fig1, axes1 = plt.subplots(1, 3, figsize=(22, 7), dpi=300)
    for i, prop in enumerate(props):
        sns.kdeplot(data=full_df, x=prop, hue='Scaffold_Label', ax=axes1[i], palette='husl', linewidth=2.0, common_norm=False)
    plt.tight_layout(); plt.savefig(os.path.join('plots', args.output)); plt.close(); print(f"  ✅ 图 1 已保存: {args.output}")

    fig2, axes2 = plt.subplots(1, 3, figsize=(22, 7), dpi=300)
    sns.set_palette(['#e74c3c', '#95a5a6'])
    for i, prop in enumerate(props):
        sns.kdeplot(data=full_df, x=prop, hue='Category', hue_order=['Generation', 'RefDataset'], ax=axes2[i], linewidth=3.5, fill=True, alpha=0.2, common_norm=False)
    plt.tight_layout(); plt.savefig(os.path.join('plots', 'library_total_distribution.png')); plt.close(); print(f"  ✅ 图 2 已保存: library_total_distribution.png")

    # 加载基准点数据
    print("\n--- Step 2.1: 正在加载基准点数据 ---")
    ref_points = load_reference_points(script_dir, pictos_labels)
    
    fig3, axes3 = plt.subplots(1, 3, figsize=(22, 7), dpi=300)
    
    # 准备图例句柄（统一的图例）
    legend_handles = [
        Line2D([0],[0],color='#e74c3c',marker='o',label='Generation'), 
        Line2D([0],[0],color='#95a5a6',marker='s',label='RefDaset'), 
        Line2D([0],[0],color='grey',alpha=0.4,label='Indiv.')
    ]
    
    # 创建标签到颜色的映射（与 seaborn kdeplot 使用的 palette 一致）
    # 使用明确的 palette 字典确保颜色映射一致
    unique_labels = sorted(gen_df['Scaffold_Label'].unique())
    palette_colors = sns.color_palette("husl", len(unique_labels))
    palette_dict = dict(zip(unique_labels, palette_colors))
    
    for i, prop in enumerate(props):
        ax = axes3[i]
        # 绘制每个标签的 KDE 曲线，直接使用原始标签
        for label in unique_labels:
            label_data = gen_df[gen_df['Scaffold_Label'] == label][prop].dropna()
            if len(label_data) > 0:
                # 将 Series 转换为 DataFrame
                label_df = pd.DataFrame({prop: label_data})
                sns.kdeplot(data=label_df, x=prop, label=label, ax=ax, 
                           color=palette_dict[label], linewidth=1.6, alpha=0.8, common_norm=False)
        
        if not other_df.empty: sns.kdeplot(data=other_df, x=prop, color='#95a5a6', marker='s', markevery=0.1, label='RefDataset', ax=ax, common_norm=False)
        sns.kdeplot(data=gen_df, x=prop, color='#e74c3c', marker='o', markevery=0.08, label='Generation', ax=ax, common_norm=False)
        
        # 标记基准点（五角星，画在对应的 KDE 曲线上）
        if ref_points:
            for pic_id, point_data in ref_points.items():
                prop_value = point_data.get(prop, None)
                label = point_data.get('label', None)
                
                if prop_value is not None and label and label in palette_dict:
                    color = palette_dict[label]
                    
                    # 获取该标签对应的数据，计算 KDE 找到基准点在曲线上的 y 值
                    label_data = gen_df[gen_df['Scaffold_Label'] == label][prop].dropna()
                    if len(label_data) > 0:
                        try:
                            # 计算 KDE
                            kde = stats.gaussian_kde(label_data)
                            # 创建 x 范围（覆盖基准点位置）
                            x_min = min(label_data.min(), prop_value)
                            x_max = max(label_data.max(), prop_value)
                            x_range = np.linspace(x_min, x_max, 500)
                            y_range = kde(x_range)
                            
                            # 使用插值找到基准点 x 位置对应的 y 值
                            interp_func = interp1d(x_range, y_range, kind='linear', 
                                                  bounds_error=False, fill_value='extrapolate')
                            y_value = interp_func(prop_value)
                            
                            # 绘制五角星在曲线上
                            ax.scatter(prop_value, y_value, marker='*', s=250, color=color, edgecolors='none', zorder=11)
                        except Exception:
                            # 如果计算失败，使用备用位置
                            y_min, y_max = ax.get_ylim()
                            ax.scatter(prop_value, y_max * 0.95, marker='*', s=250, color=color, edgecolors='none', zorder=11)
                    else:
                        # 如果没有数据，使用备用位置
                        y_min, y_max = ax.get_ylim()
                        ax.scatter(prop_value, y_max * 0.95, marker='*', s=250, color=color, edgecolors='none', zorder=11)
        
        # 显示图例（包含所有 KDE 曲线标签，使用原始标签）
        # 获取当前图例中的标签，按原始标签排序
        handles, labels = ax.get_legend_handles_labels()
        
        # 将标签分为 KDE 曲线标签和其他标签
        kde_labels = [l for l in labels if l in unique_labels]
        other_labels = [l for l in labels if l not in unique_labels]
        
        # 按照 unique_labels 的顺序排序 KDE 标签（已排序）
        kde_labels_sorted = [l for l in unique_labels if l in kde_labels]
        
        # 重新排列 handles 和 labels
        new_handles = []
        new_labels = []
        label_to_handle = dict(zip(labels, handles))
        
        # 先添加 KDE 曲线标签（按排序顺序）
        for label in kde_labels_sorted:
            if label in label_to_handle:
                new_handles.append(label_to_handle[label])
                new_labels.append(label)
        
        # 再添加其他标签（Generation, RefDaset等）
        for label in other_labels:
            if label in label_to_handle:
                new_handles.append(label_to_handle[label])
                new_labels.append(label)
        
        # 更新图例
        ax.legend(new_handles, new_labels, loc='upper right', fontsize=8, framealpha=0.9, ncol=1)
    
    plt.tight_layout(); plt.savefig(os.path.join('plots', 'generation_integrated_kde.png')); plt.close(); print(f"  ✅ 图 3 已保存: generation_integrated_kde.png")

    # --- 汇总散点图部分 ---
    print("\n--- Step 3: 生成汇总散点分布图 (固定伸缩比例) ---")
    
    # 准备参考点列表（用于汇总图）
    # SA-Wave 参考点列表
    ref_points_sa_wave = []
    for pic_id, ref_pt in ref_points.items():
        if ref_pt.get('SA_Score') is not None and ref_pt.get('S1_wavelength_nm') is not None:
            ref_points_sa_wave.append({
                'SA_Score': ref_pt['SA_Score'],
                'S1_wavelength_nm': ref_pt['S1_wavelength_nm'],
                'color': 'red'
            })
    
    # MW-S1 参考点列表
    ref_points_mw_s1 = []
    for pic_id, ref_pt in ref_points.items():
        if ref_pt.get('Molecular_Weight') is not None and ref_pt.get('S1') is not None:
            ref_points_mw_s1.append({
                'Molecular_Weight': ref_pt['Molecular_Weight'],
                'S1': ref_pt['S1'],
                'color': 'red'
            })
    
    # SA-Wave 汇总图
    plot_joint_scatter(gen_df, 'SA_Score', 'S1_wavelength_nm', "Generation Library: SA vs Wavelength", 
                       os.path.join('plots', 'generation_combined_scatter_sa_wave.png'), 
                       'SA Score', 'S1 Wavelength (nm)', (0, 8), (200, 1000), color='#e74c3c',
                       ref_points_list=ref_points_sa_wave)
    plot_joint_scatter(full_df, 'SA_Score', 'S1_wavelength_nm', "Full Library Comparison: SA vs Wavelength", 
                       os.path.join('plots', 'library_total_scatter_sa_wave_comparison.png'), 
                       'SA Score', 'S1 Wavelength (nm)', (0, 8), (200, 1000), hue_col='Category', 
                       palette={'Generation': '#e74c3c', 'RefDataset': '#95a5a6'},
                       ref_points_list=ref_points_sa_wave)

    # MW-S1 汇总图
    plot_joint_scatter(gen_df, 'Molecular_Weight', 'S1', "Generation Library: MW vs S1 Energy", 
                       os.path.join('plots', 'generation_combined_scatter_mw_s1.png'), 
                       'Molecular Weight', 'S1 Energy (eV)', (0, 1600), (1.0, 5.0), color='#e74c3c',
                       ref_points_list=ref_points_mw_s1)
    plot_joint_scatter(full_df, 'Molecular_Weight', 'S1', "Full Library Comparison: MW vs S1 Energy", 
                       os.path.join('plots', 'library_total_scatter_mw_s1_comparison.png'), 
                       'Molecular Weight', 'S1 Energy (eV)', (0, 1600), (1.0, 5.0), hue_col='Category', 
                       palette={'Generation': '#3498db', 'RefDataset': '#95a5a6'},
                       ref_points_list=ref_points_mw_s1)

    print(f"\n📊 [Done] 所有分析任务圆满完成。已确保边际分布伸缩一致。")

if __name__ == "__main__":
    main()