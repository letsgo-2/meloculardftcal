#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frontier Orbital Energies & Band Gap 绘图脚本 (增强版 v5)

功能说明：
    1. 生成小提琴图、箱线图。
    2. HOMO vs LUMO 散点图（带 1eV 投影轴）。
    3. 新增：Gap (LUMO-HOMO) vs. Dipole (D) 散点图。
"""
from __future__ import annotations

from typing import List, Dict, Optional, Tuple
import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines
from scipy.stats import gaussian_kde 


def load_pictos_labels(pictos_file: Path) -> Dict[str, str]:
    labels = {}
    if not pictos_file.exists():
        print(f"警告: pictos文件不存在: {pictos_file}")
        return labels
    try:
        with open(pictos_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or ':' not in line:
                    continue
                pic_name, label = line.split(':', 1)
                labels[pic_name.strip()] = label.strip()
    except Exception as e:
        print(f"警告: 无法读取pictos文件: {e}")
    return labels


def get_label_from_filename(filename: str, pictos_labels: Dict[str, str]) -> str:
    pic_name = None
    for pic in pictos_labels.keys():
        if filename.startswith(pic):
            pic_name = pic
            break
    if pic_name and pic_name in pictos_labels:
        label = pictos_labels[pic_name]
        if label.lower().startswith('s'):
            return f"S{label[1:].upper()}"
        return label
    return Path(filename).stem


def sort_labels_by_s_order(labels: List[str]) -> List[str]:
    def get_s_number(label):
        label_upper = label.upper()
        if label_upper.startswith('S'):
            try:
                return int(label_upper[1:])
            except:
                return 999
        return 999
    sorted_labels = sorted(labels, key=get_s_number)
    formatted_labels = []
    for label in sorted_labels:
        if label.lower().startswith('s'):
            formatted_labels.append(f"S{label[1:].upper()}")
        else:
            formatted_labels.append(label)
    return formatted_labels


def extract_data_from_csv(csv_file: Path, columns: List[str]) -> Dict[str, List[float]]:
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"警告: 无法读取文件 {csv_file}: {e}")
        return {}
    result = {}
    for col in columns:
        if col not in df.columns:
            # 不打印列缺失警告，以保持原有输出简洁
            result[col] = []
            continue
        valid_data = df[col].dropna().astype(float).tolist()
        result[col] = valid_data
    return result


def prepare_plot_positions(all_data, file_labels, columns):
    width_between_files = 3
    width_between_columns = 1.0
    all_positions = []
    all_data_list = []
    all_colors = []
    column_colors = ['#FFB3BA', '#BAE1FF', '#D4F1BE'] 
    
    for file_idx, file_label in enumerate(file_labels):
        base_position = file_idx * width_between_files + 1
        for col_idx, column in enumerate(columns):
            if column in all_data[file_label] and len(all_data[file_label][column]) > 0:
                position = base_position + col_idx * width_between_columns
                all_positions.append(position)
                all_data_list.append(all_data[file_label][column])
                all_colors.append(column_colors[col_idx % len(column_colors)])
    return all_positions, all_data_list, all_colors


def add_annotations_and_style(ax, all_data_list, all_positions, columns, file_labels):
    width_between_files = 3
    width_between_columns = 1.0
    for i, (data, position) in enumerate(zip(all_data_list, all_positions)):
        if len(data) > 0:
            mean_val = np.mean(data)
            min_val = np.min(data)
            y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
            col_in_group_idx = i % len(columns)
            current_col_name = columns[col_in_group_idx].upper()
            ax.text(position, mean_val + y_range * 0.02, f'Mean: {mean_val:.3f}',
                    ha='center', va='bottom', fontsize=8, color='red', fontweight='bold')
            if "GAP" in current_col_name:
                ax.text(position, min_val - y_range * 0.02, f'Min: {min_val:.3f}',
                        ha='center', va='top', fontsize=8, color='darkgreen', fontweight='bold')
    file_center_positions = [i * width_between_files + 1 + (len(columns)-1)*width_between_columns/2 for i in range(len(file_labels))]
    ax.set_xticks(file_center_positions)
    ax.set_xticklabels(file_labels, rotation=45, ha='right')
    ax.set_title('Frontier Orbital Energies & Band Gap', fontsize=12, fontweight='bold', pad=15)
    ax.set_ylabel('Energy (eV)', fontsize=10, fontweight='bold')
    ax.set_xlabel('Scaffold', fontsize=10, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

def create_legend(ax, columns, column_colors):
    def format_name(n):
        u = n.upper()
        if 'HOMO' in u: return 'HOMO'
        if 'LUMO' in u: return 'LUMO'
        if 'GAP' in u: return 'Gap'
        return n
    legend_elements = []
    for idx, col in enumerate(columns):
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=column_colors[idx%len(column_colors)], 
                                             alpha=0.5, edgecolor='black', label=format_name(col)))
    legend_elements.append(mlines.Line2D([0], [0], color='red', linestyle='--', linewidth=2, label='Mean'))
    legend_elements.append(mlines.Line2D([0], [0], color='blue', linestyle='-', linewidth=2, label='Median'))
    ax.legend(handles=legend_elements, loc='upper right', framealpha=0.9)


def plot_violinplot(all_data, output_path, columns, figsize, file_labels,
                    pictos_labels: Optional[Dict[str, str]] = None, 
                    script_dir: Optional[Path] = None):
    fig, ax = plt.subplots(figsize=figsize)
    # 过滤掉非能级列
    target_cols = [c for c in columns if 'HOMO' in c.upper() or 'LUMO' in c.upper() or 'GAP' in c.upper()]
    all_positions, all_data_list, all_colors = prepare_plot_positions(all_data, file_labels, target_cols)
    parts = ax.violinplot(all_data_list, positions=all_positions, showmeans=True, showmedians=True, widths=0.6)
    for pc, color in zip(parts['bodies'], all_colors):
        pc.set_facecolor(color); pc.set_alpha(0.5); pc.set_edgecolor('black'); pc.set_linewidth(1.0)
    if 'cmeans' in parts: parts['cmeans'].set_color('red'); parts['cmeans'].set_linewidth(2.0); parts['cmeans'].set_linestyle('--')
    if 'cmedians' in parts: parts['cmedians'].set_color('blue'); parts['cmedians'].set_linewidth(2.0)
    for partname in ('cbars', 'cmins', 'cmaxes'):
        if partname in parts: parts[partname].set_color('black')
    
    # 加载并标记参考点
    if pictos_labels and script_dir:
        import re
        paths_to_try = [
            script_dir.parent / "data" / "publishinpdftv22.csv",
            script_dir / "data" / "publishinpdftv22.csv",
            Path("data/publishinpdftv22.csv")
        ]
        publish_file = None
        for p in paths_to_try:
            if p.exists():
                publish_file = p
                break
        
        if publish_file:
            try:
                ref_df = pd.read_csv(publish_file)
                width_between_files = 3
                width_between_columns = 1.0
                
                for _, row in ref_df.iterrows():
                    index_val = str(row.get('index', '')).strip()
                    if not index_val or not index_val.startswith('pic'):
                        continue
                    
                    match = re.search(r'pic(\d+)', index_val, re.IGNORECASE)
                    if not match:
                        continue
                    pic_id = match.group(1)
                    pic_key = f"pic{pic_id}".lower()
                    
                    raw_label = pictos_labels.get(pic_key, None)
                    if not raw_label:
                        continue
                    
                    if raw_label.lower().startswith('s'):
                        batch_label = f"S{raw_label[1:].upper()}"
                    else:
                        batch_label = raw_label
                    
                    # 检查 batch_label 是否在 file_labels 中
                    if batch_label not in file_labels:
                        continue
                    
                    file_idx = file_labels.index(batch_label)
                    base_position = file_idx * width_between_files + 1
                    
                    # 标记 HOMO, LUMO, Gap
                    for col_name in ['HOMO(eV)', 'LUMO(eV)', 'Gap(eV)']:
                        if col_name not in target_cols:
                            continue
                        
                        ref_val = row.get(col_name, np.nan)
                        if pd.notna(ref_val):
                            try:
                                ref_val = float(ref_val)
                                col_idx = target_cols.index(col_name)
                                position = base_position + col_idx * width_between_columns
                                # 稍微侧移（向右偏移 0.2）
                                ax.scatter(position - 0.2, ref_val, s=160, color='#e74c3c', 
                                          edgecolors='black', linewidth=0.8, marker='*', 
                                          zorder=10)
                            except (ValueError, TypeError):
                                continue
            except Exception as e:
                print(f"  ⚠️ 加载参考点数据时出错: {e}")
    
    add_annotations_and_style(ax, all_data_list, all_positions, target_cols, file_labels)
    create_legend(ax, target_cols, ['#FFB3BA', '#BAE1FF', '#D4F1BE'])
    plt.tight_layout(); plt.savefig(output_path, dpi=300, bbox_inches='tight'); plt.close()

def plot_boxplot(all_data, output_path, columns, figsize, file_labels):
    fig, ax = plt.subplots(figsize=figsize)
    target_cols = [c for c in columns if 'HOMO' in c.upper() or 'LUMO' in c.upper() or 'GAP' in c.upper()]
    all_positions, all_data_list, all_colors = prepare_plot_positions(all_data, file_labels, target_cols)
    bp = ax.boxplot(all_data_list, positions=all_positions, widths=0.6, patch_artist=True,
                    showmeans=True, meanline=True, 
                    medianprops={'color': 'blue', 'linewidth': 2.0},
                    meanprops={'color': 'red', 'linestyle': '--', 'linewidth': 2.0},
                    flierprops={'marker': 'o', 'markerfacecolor': 'gray', 'markersize': 4, 'alpha': 0.5})
    for patch, color in zip(bp['boxes'], all_colors):
        patch.set_facecolor(color); patch.set_alpha(0.5); patch.set_edgecolor('black')
    add_annotations_and_style(ax, all_data_list, all_positions, target_cols, file_labels)
    create_legend(ax, target_cols, ['#FFB3BA', '#BAE1FF', '#D4F1BE'])
    plt.tight_layout(); plt.savefig(output_path, dpi=300, bbox_inches='tight'); plt.close()


def plot_scatter_homo_lumo(all_data, output_path, file_labels):
    fig, ax = plt.subplots(figsize=(10, 10))
    cmap = plt.get_cmap('tab10')
    homo_col, lumo_col = None, None
    for label in all_data:
        for col in all_data[label]:
            if 'HOMO' in col.upper(): homo_col = col
            if 'LUMO' in col.upper(): lumo_col = col
        if homo_col and lumo_col: break
    if not homo_col or not lumo_col: return

    all_h, all_l, all_g = [], [], []
    for idx, label in enumerate(file_labels):
        if label in all_data:
            h = all_data[label].get(homo_col, [])
            l = all_data[label].get(lumo_col, [])
            # 只保留同时有 HOMO 和 LUMO 数据的点对
            if len(h) > 0 and len(l) > 0:
                # 取最小长度，确保数据点数量一致
                min_len = min(len(h), len(l))
                h_matched = h[:min_len]
                l_matched = l[:min_len]
                if min_len > 0:
                    ax.scatter(h_matched, l_matched, label=label, color=cmap(idx % 10), alpha=0.5, edgecolors='none', s=25, zorder=3)
                    all_h.extend(h_matched); all_l.extend(l_matched)
                    all_g.extend(np.array(l_matched) - np.array(h_matched))
                    if len(h) != len(l):
                        print(f"  ⚠ {label}: HOMO数据点={len(h)}, LUMO数据点={len(l)}, 已匹配{min_len}个点对")
    if not all_h: return

    ax.set_xlim(-8.0, -4.0); ax.set_ylim(-5.0, -1.0)
    x_range = np.array(ax.get_xlim()); target_gap = 1
    ax.plot(x_range, x_range + target_gap, color='gray', linestyle='--', linewidth=1.2, alpha=0.6, zorder=1)
    ax.text(x_range[1]-0.3, x_range[1] + target_gap, f' Gap={target_gap}eV', fontsize=10, color='gray', va='center', fontweight='bold', rotation=45)

    g_range = np.linspace(min(all_g)-0.2, max(all_g)+0.2, 200)
    kde = gaussian_kde(all_g); density = kde(g_range); density = density / max(density) * 0.82 
    C = (min(all_h) + max(all_l)) * 0.98
    line_x = (C - g_range) / 2; line_y = (C + g_range) / 2
    ax.plot(line_x, line_y, color='black', linewidth=1.5, alpha=0.9, zorder=4)
    ticks_val = np.arange(0.5, 4.5 , 0.5); tick_len = 0.05 
    for tv in ticks_val:
        tx0 = (C - tv) / 2; ty0 = (C + tv) / 2
        ax.plot([tx0, tx0 - tick_len], [ty0, ty0 - tick_len], color='black', linewidth=1, zorder=5)
        ax.text(tx0 - tick_len*1.5, ty0 - tick_len*1.5, f'{tv}', fontsize=7, ha='right', va='top', color='black')

    offset_x = line_x + density * 0.5; offset_y = line_y + density * 0.5
    ax.plot(offset_x, offset_y, color='#789661', linewidth=1, alpha=0.8, zorder=5)
    ax.set_title('HOMO vs LUMO Energy Distribution', fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('HOMO (eV)', fontsize=12, fontweight='bold'); ax.set_ylabel('LUMO (eV)', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.3, zorder=0); ax.legend(title='Scaffold', bbox_to_anchor=(0.98, 1), loc='upper left', framealpha=0.8)
    plt.tight_layout(); plt.savefig(output_path, dpi=300, bbox_inches='tight'); plt.close()


def plot_gap_dipole_scatter(all_data: Dict, output_path: Path, file_labels: List[str], 
                            pictos_labels: Optional[Dict[str, str]] = None, 
                            script_dir: Optional[Path] = None):
    """
    更新：Gap (LUMO-HOMO) vs. Dipole 散点图
    1. 生成汇总图（All Scaffolds）：保留图例，边际图显示总体分布（灰色）。
    2. 生成每一个骨架的单图：不带图例，边际图颜色与该骨架点的颜色一致，标题注明骨架。
    3. 在所有图中标记参考点（来自 publishinpdftv22.csv）
    """
    import matplotlib.gridspec as gridspec
    from scipy.stats import gaussian_kde
    import re

    # 1. 获取列名逻辑 (保持不变)
    gap_col, dipole_col = None, None
    for label in all_data:
        for col in all_data[label]:
            u_col = col.upper()
            if 'GAP' in u_col: gap_col = col
            if 'DIPOLE' in u_col: dipole_col = col
        if gap_col and dipole_col: break
    
    if not gap_col or not dipole_col:
        return

    # 2. 加载参考点数据
    ref_points_dict = {}  # {batch_label: {'Gap': float, 'Dipole': float}}
    if pictos_labels and script_dir:
        paths_to_try = [
            script_dir.parent / "data" / "publishinpdftv22.csv",
            script_dir / "data" / "publishinpdftv22.csv",
            Path("data/publishinpdftv22.csv")
        ]
        publish_file = None
        for p in paths_to_try:
            if p.exists():
                publish_file = p
                break
        
        if publish_file:
            try:
                ref_df = pd.read_csv(publish_file)
                for _, row in ref_df.iterrows():
                    index_val = str(row.get('index', '')).strip()
                    if not index_val or not index_val.startswith('pic'):
                        continue
                    
                    match = re.search(r'pic(\d+)', index_val, re.IGNORECASE)
                    if not match:
                        continue
                    pic_id = match.group(1)
                    pic_key = f"pic{pic_id}".lower()
                    
                    raw_label = pictos_labels.get(pic_key, None)
                    if not raw_label:
                        continue
                    
                    if raw_label.lower().startswith('s'):
                        batch_label = f"S{raw_label[1:].upper()}"
                    else:
                        batch_label = raw_label
                    
                    gap_val = row.get('Gap(eV)', np.nan)
                    dipole_val = row.get('Dipole(D)', np.nan)
                    if pd.notna(gap_val) and pd.notna(dipole_val):
                        try:
                            ref_points_dict[batch_label] = {
                                'Gap': float(gap_val),
                                'Dipole': float(dipole_val)
                            }
                            print(f"  ✅ 已加载参考点 {index_val} -> {batch_label}: Gap={float(gap_val):.3f}, Dipole={float(dipole_val):.2f}")
                        except (ValueError, TypeError):
                            continue
            except Exception as e:
                print(f"  ⚠️ 加载参考点数据时出错: {e}")

    # 定义内部绘图逻辑
    def create_single_plot(data_to_plot, scaffold_title, save_path, marginal_color, show_legend=False, ref_points_to_plot=None):
        fig = plt.figure(figsize=(10, 10))
        # 布局：边际图占 1/7 比例
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 6], width_ratios=[6, 1], 
                               hspace=0.05, wspace=0.05)
        
        ax_main = fig.add_subplot(gs[1, 0])      # 主散点图
        ax_marg_x = fig.add_subplot(gs[0, 0], sharex=ax_main)  # 顶部分布
        ax_marg_y = fig.add_subplot(gs[1, 1], sharey=ax_main)  # 右侧分布

        cmap = plt.get_cmap('tab10')
        combined_g, combined_d = [], []

        # 绘制散点
        for lbl, vals in data_to_plot.items():
            g_vals = vals.get(gap_col, [])
            d_vals = vals.get(dipole_col, [])
            if len(g_vals) > 0 and len(d_vals) > 0:
                min_len = min(len(g_vals), len(d_vals))
                g_m = np.array(g_vals[:min_len])
                d_m = np.array(d_vals[:min_len])
                
                # 确定颜色：如果是汇总图则按索引从 tab10 取色，单图则直接使用 marginal_color
                idx_in_labels = file_labels.index(lbl) if lbl in file_labels else 0
                point_color = cmap(idx_in_labels % 10) if len(data_to_plot) > 1 else marginal_color
                
                ax_main.scatter(g_m, d_m, label=lbl, color=point_color, alpha=0.6, edgecolors='none', s=35, zorder=3)
                combined_g.extend(g_m)
                combined_d.extend(d_m)

        # 绘制边际分布 (KDE)
        if combined_g and combined_d:
            gx = np.linspace(min(combined_g), max(combined_g), 200)
            gk = gaussian_kde(combined_g)(gx)
            ax_marg_x.plot(gx, gk, color=marginal_color, lw=1.5)
            ax_marg_x.fill_between(gx, 0, gk, color=marginal_color, alpha=0.3)
            
            dy = np.linspace(min(combined_d), max(combined_d), 200)
            dk = gaussian_kde(combined_d)(dy)
            ax_marg_y.plot(dk, dy, color=marginal_color, lw=1.5)
            ax_marg_y.fill_betweenx(dy, 0, dk, color=marginal_color, alpha=0.3)

        # 隐藏边际图轴线
        ax_marg_x.axis('off')
        ax_marg_y.axis('off')

        # 主图样式
        ax_main.set_xlabel('Band Gap (eV)', fontsize=12, fontweight='bold')
        ax_main.set_ylabel('Dipole Moment (D)', fontsize=12, fontweight='bold')
        ax_main.grid(True, linestyle='--', alpha=0.3, zorder=0)
        
        # 单个骨架图设置坐标轴范围
        if scaffold_title != "All Scaffolds":
            ax_main.set_xlim(left=0, right=5)
            ax_main.set_ylim(bottom=0, top=140)
        
        fig.suptitle(f'Band Gap vs. Dipole Moment - {scaffold_title}', 
                     fontsize=14, fontweight='bold', y=0.9)
        
        # 标记参考点（红色五角星）
        if ref_points_to_plot:
            ref_point_label_added = False
            for ref_label, ref_data in ref_points_to_plot.items():
                gap_val = ref_data.get('Gap', None)
                dipole_val = ref_data.get('Dipole', None)
                if gap_val is not None and dipole_val is not None:
                    label = 'Ref Point' if not ref_point_label_added else ''
                    ax_main.scatter(gap_val, dipole_val, s=180, color='#e74c3c', 
                                   edgecolors='black', linewidth=0.8, marker='*', 
                                   zorder=10, label=label)
                    ref_point_label_added = True
        
        # 只有汇总图才显示图例
        if show_legend:
            ax_main.legend(title='Scaffold', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

    # --- 1. 生成汇总图 (保留图例，边际图设为中性灰以突出总体分布) ---
    combined_dict = {lbl: all_data[lbl] for lbl in file_labels if lbl in all_data}
    # 汇总图显示所有参考点
    create_single_plot(combined_dict, "All Scaffolds", output_path, marginal_color='#4a69bd', 
                      show_legend=True, ref_points_to_plot=ref_points_dict)

    # --- 2. 遍历生成每个骨架的独立图 (移除图例，颜色同步) ---
    for idx, label in enumerate(file_labels):
        if label in all_data:
            single_scaffold_dict = {label: all_data[label]}
            individual_path = output_path.parent / f"{output_path.stem}_{label}{output_path.suffix}"
            
            # 使用对应骨架的特定颜色
            scaffold_color = plt.get_cmap('tab10')(idx % 10)
            
            # 单个骨架图只显示对应的参考点
            single_ref_point = {label: ref_points_dict[label]} if label in ref_points_dict else None
            
            create_single_plot(
                single_scaffold_dict, 
                f"Scaffold {label}", 
                individual_path, 
                marginal_color=scaffold_color, 
                show_legend=False,
                ref_points_to_plot=single_ref_point
            )


def plot_gap_s1t1_scatter(all_data: Dict, output_path: Path, file_labels: List[str], 
                          pictos_labels: Optional[Dict[str, str]] = None, 
                          script_dir: Optional[Path] = None):
    """
    更新：Gap (LUMO-HOMO) vs. S1-T1 散点图
    1. 生成汇总图（All Scaffolds）：保留图例，边际图显示总体分布（灰色）。
    2. 生成每一个骨架的单图：不带图例，边际图颜色与该骨架点的颜色一致，标题注明骨架。
    3. 在所有图中标记参考点（来自 publishinpdftv22.csv）
    """
    import matplotlib.gridspec as gridspec
    from scipy.stats import gaussian_kde
    import re

    # 1. 获取列名逻辑
    gap_col, s1t1_col = None, None
    all_columns = set()
    for label in all_data:
        for col in all_data[label]:
            all_columns.add(col)
            u_col = col.upper()
            if 'GAP' in u_col and gap_col is None: 
                gap_col = col
            # 查找 S1-T1 列：精确匹配 "S1-T1" 或包含 "S1" 和 "T1" 且中间有 "-"
            if s1t1_col is None:
                if 'S1-T1' in col or (col.upper().startswith('S1') and 'T1' in col and '-' in col):
                    s1t1_col = col
        if gap_col and s1t1_col: break
    
    # 添加调试信息
    if not gap_col:
        print(f"  ⚠️ plot_gap_s1t1_scatter: 未找到 Gap 列")
        if all_columns:
            print(f"     可用列名: {sorted(all_columns)}")
        return
    if not s1t1_col:
        print(f"  ⚠️ plot_gap_s1t1_scatter: 未找到 S1-T1 列")
        if all_columns:
            print(f"     可用列名: {sorted(all_columns)}")
            # 尝试查找类似的列
            similar_cols = [c for c in all_columns if 'S1' in c.upper() or 'T1' in c.upper()]
            if similar_cols:
                print(f"     相关列名: {similar_cols}")
        return
    
    print(f"  ✅ plot_gap_s1t1_scatter: 找到列 - Gap: {gap_col}, S1-T1: {s1t1_col}")

    # 2. 加载参考点数据
    ref_points_dict = {}  # {batch_label: {'Gap': float, 'S1T1': float}}
    if pictos_labels and script_dir:
        paths_to_try = [
            script_dir.parent / "data" / "publishinpdftv22.csv",
            script_dir / "data" / "publishinpdftv22.csv",
            Path("data/publishinpdftv22.csv")
        ]
        publish_file = None
        for p in paths_to_try:
            if p.exists():
                publish_file = p
                break
        
        if publish_file:
            try:
                ref_df = pd.read_csv(publish_file)
                for _, row in ref_df.iterrows():
                    index_val = str(row.get('index', '')).strip()
                    if not index_val or not index_val.startswith('pic'):
                        continue
                    
                    match = re.search(r'pic(\d+)', index_val, re.IGNORECASE)
                    if not match:
                        continue
                    pic_id = match.group(1)
                    pic_key = f"pic{pic_id}".lower()
                    
                    raw_label = pictos_labels.get(pic_key, None)
                    if not raw_label:
                        continue
                    
                    if raw_label.lower().startswith('s'):
                        batch_label = f"S{raw_label[1:].upper()}"
                    else:
                        batch_label = raw_label
                    
                    gap_val = row.get('Gap(eV)', np.nan)
                    s1t1_val = row.get('S1-T1', np.nan)
                    if pd.notna(gap_val) and pd.notna(s1t1_val):
                        try:
                            ref_points_dict[batch_label] = {
                                'Gap': float(gap_val),
                                'S1T1': float(s1t1_val)
                            }
                            print(f"  ✅ 已加载参考点 {index_val} -> {batch_label}: Gap={float(gap_val):.3f}, S1-T1={float(s1t1_val):.3f}")
                        except (ValueError, TypeError):
                            continue
            except Exception as e:
                print(f"  ⚠️ 加载参考点数据时出错: {e}")

    # 定义内部绘图逻辑
    def create_single_plot(data_to_plot, scaffold_title, save_path, marginal_color, show_legend=False, ref_points_to_plot=None):
        fig = plt.figure(figsize=(10, 10))
        # 布局：边际图占 1/7 比例
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 6], width_ratios=[6, 1], 
                               hspace=0.05, wspace=0.05)
        
        ax_main = fig.add_subplot(gs[1, 0])      # 主散点图
        ax_marg_x = fig.add_subplot(gs[0, 0], sharex=ax_main)  # 顶部分布
        ax_marg_y = fig.add_subplot(gs[1, 1], sharey=ax_main)  # 右侧分布

        cmap = plt.get_cmap('tab10')
        combined_g, combined_s1t1 = [], []

        # 绘制散点
        for lbl, vals in data_to_plot.items():
            g_vals = vals.get(gap_col, [])
            s1t1_vals = vals.get(s1t1_col, [])
            if len(g_vals) > 0 and len(s1t1_vals) > 0:
                min_len = min(len(g_vals), len(s1t1_vals))
                g_m = np.array(g_vals[:min_len])
                s1t1_m = np.array(s1t1_vals[:min_len])
                
                # 确定颜色：如果是汇总图则按索引从 tab10 取色，单图则直接使用 marginal_color
                idx_in_labels = file_labels.index(lbl) if lbl in file_labels else 0
                point_color = cmap(idx_in_labels % 10) if len(data_to_plot) > 1 else marginal_color
                
                ax_main.scatter(g_m, s1t1_m, label=lbl, color=point_color, alpha=0.6, edgecolors='none', s=35, zorder=3)
                combined_g.extend(g_m)
                combined_s1t1.extend(s1t1_m)

        # 绘制边际分布 (KDE)
        if combined_g and combined_s1t1:
            gx = np.linspace(min(combined_g), max(combined_g), 200)
            gk = gaussian_kde(combined_g)(gx)
            ax_marg_x.plot(gx, gk, color=marginal_color, lw=1.5)
            ax_marg_x.fill_between(gx, 0, gk, color=marginal_color, alpha=0.3)
            
            dy = np.linspace(min(combined_s1t1), max(combined_s1t1), 200)
            dk = gaussian_kde(combined_s1t1)(dy)
            ax_marg_y.plot(dk, dy, color=marginal_color, lw=1.5)
            ax_marg_y.fill_betweenx(dy, 0, dk, color=marginal_color, alpha=0.3)

        # 隐藏边际图轴线
        ax_marg_x.axis('off')
        ax_marg_y.axis('off')

        # 主图样式
        ax_main.set_xlabel('Band Gap (eV)', fontsize=12, fontweight='bold')
        ax_main.set_ylabel('S1-T1 (eV)', fontsize=12, fontweight='bold')
        ax_main.grid(True, linestyle='--', alpha=0.3, zorder=0)
        
        # 单个骨架图设置坐标轴范围
        if scaffold_title != "All Scaffolds":
            ax_main.set_xlim(left=0, right=5)
            # y轴范围需要根据实际数据调整，暂时不设置固定上限
        
        fig.suptitle(f'Band Gap vs. S1-T1 - {scaffold_title}', 
                     fontsize=14, fontweight='bold', y=0.9)
        
        # 标记参考点（红色五角星）
        if ref_points_to_plot:
            ref_point_label_added = False
            for ref_label, ref_data in ref_points_to_plot.items():
                gap_val = ref_data.get('Gap', None)
                s1t1_val = ref_data.get('S1T1', None)
                if gap_val is not None and s1t1_val is not None:
                    label = 'Ref Point' if not ref_point_label_added else ''
                    ax_main.scatter(gap_val, s1t1_val, s=180, color='#e74c3c', 
                                   edgecolors='black', linewidth=0.8, marker='*', 
                                   zorder=10, label=label)
                    ref_point_label_added = True
        
        # 只有汇总图才显示图例
        if show_legend:
            ax_main.legend(title='Scaffold', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

    # --- 1. 生成汇总图 (保留图例，边际图设为中性灰以突出总体分布) ---
    combined_dict = {lbl: all_data[lbl] for lbl in file_labels if lbl in all_data}
    # 汇总图显示所有参考点
    create_single_plot(combined_dict, "All Scaffolds", output_path, marginal_color='#4a69bd', 
                      show_legend=True, ref_points_to_plot=ref_points_dict)

    # --- 2. 遍历生成每个骨架的独立图 (移除图例，颜色同步) ---
    for idx, label in enumerate(file_labels):
        if label in all_data:
            single_scaffold_dict = {label: all_data[label]}
            individual_path = output_path.parent / f"{output_path.stem}_{label}{output_path.suffix}"
            
            # 使用对应骨架的特定颜色
            scaffold_color = plt.get_cmap('tab10')(idx % 10)
            
            # 单个骨架图只显示对应的参考点
            single_ref_point = {label: ref_points_dict[label]} if label in ref_points_dict else None
            
            create_single_plot(
                single_scaffold_dict, 
                f"Scaffold {label}", 
                individual_path, 
                marginal_color=scaffold_color, 
                show_legend=False,
                ref_points_to_plot=single_ref_point
            )

def plot_charge_dipole_scatter(all_data: Dict, output_path: Path, file_labels: List[str],
                                pictos_labels: Optional[Dict[str, str]] = None, 
                                script_dir: Optional[Path] = None):
    """
    新增：Charge vs. Dipole 散点分布图
    1. 生成汇总图（All Scaffolds）：Charge 增加抖动(Jitter)以便于观察。
    2. 生成每一个骨架的单图。
    3. 在每个 Charge (0-6) 位置增加垂直方向的 Dipole 边际分布。
    4. 在所有图中标记参考点（来自 publishinpdftv22.csv）
    """
    import matplotlib.gridspec as gridspec
    from scipy.stats import gaussian_kde
    import numpy as np
    import re
    
    # 查找列名
    charge_col, dipole_col = None, None
    for label in all_data:
        for col in all_data[label]:
            u_col = col.upper()
            if 'CHARGE' in u_col and charge_col is None: charge_col = col
            if 'DIPOLE' in u_col and dipole_col is None: dipole_col = col
        if charge_col and dipole_col: break
    
    # 添加调试信息
    if not charge_col:
        print(f"  ⚠️ plot_charge_dipole_scatter: 未找到 Charge 列")
        return
    if not dipole_col:
        print(f"  ⚠️ plot_charge_dipole_scatter: 未找到 Dipole 列")
        return
    
    print(f"  ✅ plot_charge_dipole_scatter: 找到列 - Charge: {charge_col}, Dipole: {dipole_col}")

    # 加载参考点数据
    ref_points_dict = {}  # {batch_label: {'Charge': float, 'Dipole': float}}
    if pictos_labels and script_dir:
        paths_to_try = [
            script_dir.parent / "data" / "publishinpdftv22.csv",
            script_dir / "data" / "publishinpdftv22.csv",
            Path("data/publishinpdftv22.csv")
        ]
        publish_file = None
        for p in paths_to_try:
            if p.exists():
                publish_file = p
                break
        
        if publish_file:
            try:
                ref_df = pd.read_csv(publish_file)
                for _, row in ref_df.iterrows():
                    index_val = str(row.get('index', '')).strip()
                    if not index_val or not index_val.startswith('pic'):
                        continue
                    
                    match = re.search(r'pic(\d+)', index_val, re.IGNORECASE)
                    if not match:
                        continue
                    pic_id = match.group(1)
                    pic_key = f"pic{pic_id}".lower()
                    
                    raw_label = pictos_labels.get(pic_key, None)
                    if not raw_label:
                        continue
                    
                    if raw_label.lower().startswith('s'):
                        batch_label = f"S{raw_label[1:].upper()}"
                    else:
                        batch_label = raw_label
                    
                    charge_val = row.get('charge', np.nan)
                    dipole_val = row.get('Dipole(D)', np.nan)
                    if pd.notna(charge_val) and pd.notna(dipole_val):
                        try:
                            ref_points_dict[batch_label] = {
                                'Charge': float(charge_val),
                                'Dipole': float(dipole_val)
                            }
                            print(f"  ✅ 已加载参考点 {index_val} -> {batch_label}: Charge={float(charge_val):.0f}, Dipole={float(dipole_val):.2f}")
                        except (ValueError, TypeError):
                            continue
            except Exception as e:
                print(f"  ⚠️ 加载参考点数据时出错: {e}")

    def create_single_plot(data_to_plot, scaffold_title, save_path, marginal_color, show_legend=False, ref_points_to_plot=None):
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 6], width_ratios=[6, 1], hspace=0.05, wspace=0.05)
        ax_main = fig.add_subplot(gs[1, 0]); ax_marg_x = fig.add_subplot(gs[0, 0], sharex=ax_main); ax_marg_y = fig.add_subplot(gs[1, 1], sharey=ax_main)
        cmap = plt.get_cmap('tab10'); combined_c, combined_d = [], []

        # 1. 绘制散点
        for lbl, vals in data_to_plot.items():
            c_vals = vals.get(charge_col, []); d_vals = vals.get(dipole_col, [])
            if len(c_vals) > 0 and len(d_vals) > 0:
                min_len = min(len(c_vals), len(d_vals))
                c_m = np.array(c_vals[:min_len]); d_m = np.array(d_vals[:min_len])
                # 增加抖动
                jitter = np.random.normal(0, 0.08, size=len(c_m))
                idx = file_labels.index(lbl) if lbl in file_labels else 0
                point_color = cmap(idx % 10) if len(data_to_plot) > 1 else marginal_color
                ax_main.scatter(c_m + jitter, d_m, label=lbl, color=point_color, alpha=0.5, edgecolors='none', s=35, zorder=3)
                combined_c.extend(c_m); combined_d.extend(d_m)

        # 2. 在每个电荷位置 (0-6) 增加垂直边际分布 (KDE)
        if combined_c and combined_d:
            unique_charges = np.arange(0, 7)
            for x_pos in unique_charges:
                # 提取当前电荷下的所有偶极矩数据
                current_d = []
                for lbl, vals in data_to_plot.items():
                    c_vals = vals.get(charge_col, [])
                    d_vals = vals.get(dipole_col, [])
                    if len(c_vals) > 0 and len(d_vals) > 0:
                        # 确保长度一致
                        min_len = min(len(c_vals), len(d_vals))
                        c_arr = np.array(c_vals[:min_len])
                        d_arr = np.array(d_vals[:min_len])
                        # 匹配当前电荷值
                        mask = (np.round(c_arr) == x_pos)
                        current_d.extend(d_arr[mask].tolist())
                
                # 如果当前电荷有足够的数据点（>10）且数据有方差，则绘制 KDE
                if len(current_d) > 10:
                    current_d_arr = np.array(current_d)
                    # 检查数据是否有方差（不是所有值都相同）
                    if np.std(current_d_arr) > 1e-6:  # 避免方差为0的情况
                        try:
                            d_eval = np.linspace(min(combined_d), max(combined_d), 200)
                            kde = gaussian_kde(current_d_arr, bw_method=0.3)
                            kde_val = kde(d_eval)
                            # 缩放 KDE 宽度使其不重叠（最大宽度 0.4）
                            if kde_val.max() > 0:
                                kde_val = (kde_val / kde_val.max()) * 0.4
                                
                                # 绘制在散点图背景上
                                ax_main.fill_betweenx(d_eval, x_pos, x_pos + kde_val, 
                                                      color=marginal_color, alpha=0.15, zorder=1)
                                ax_main.plot(x_pos + kde_val, d_eval, color=marginal_color, lw=1, alpha=0.3, zorder=2)
                        except Exception as e:
                            # 如果 KDE 失败（例如数据共线），跳过
                            pass

        # 3. 绘制顶部和右侧的总体边际分布
        if combined_c and combined_d:
            cx = np.linspace(-0.5, 6.5, 300); ck = gaussian_kde(combined_c, bw_method=0.2)(cx)
            ax_marg_x.plot(cx, ck, color=marginal_color, lw=1.5); ax_marg_x.fill_between(cx, 0, ck, color=marginal_color, alpha=0.3)
            dy = np.linspace(min(combined_d), max(combined_d), 200); dk = gaussian_kde(combined_d)(dy)
            ax_marg_y.plot(dk, dy, color=marginal_color, lw=1.5); ax_marg_y.fill_betweenx(dy, 0, dk, color=marginal_color, alpha=0.3)

        ax_marg_x.axis('off'); ax_marg_y.axis('off')
        ax_main.set_xlabel('Charge State', fontweight='bold'); ax_main.set_ylabel('Dipole Moment (D)', fontweight='bold')
        ax_main.set_xlim(left=-0.5, right=6.5)
        ax_main.set_xticks(np.arange(0, 7))
        ax_main.grid(True, linestyle='--', alpha=0.3)
        fig.suptitle(f'Charge vs. Dipole Moment - {scaffold_title}', fontsize=14, fontweight='bold', y=0.9)
        
        # 标记参考点（红色五角星）
        if ref_points_to_plot:
            ref_point_label_added = False
            for ref_label, ref_data in ref_points_to_plot.items():
                charge_val = ref_data.get('Charge', None)
                dipole_val = ref_data.get('Dipole', None)
                if charge_val is not None and dipole_val is not None:
                    label = 'Ref Point' if not ref_point_label_added else ''
                    ax_main.scatter(charge_val, dipole_val, s=180, color='#e74c3c', 
                                   edgecolors='black', linewidth=0.8, marker='*', 
                                   zorder=10, label=label)
                    ref_point_label_added = True
        
        if show_legend: ax_main.legend(title='Scaffold', bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.savefig(save_path, dpi=300, bbox_inches='tight'); plt.close()

    # --- 1. 生成汇总图 (保留图例，边际图设为蓝色以突出总体分布) ---
    combined_dict = {lbl: all_data[lbl] for lbl in file_labels if lbl in all_data}
    # 汇总图显示所有参考点
    create_single_plot(combined_dict, "All Scaffolds", output_path, '#4a69bd', True, ref_points_to_plot=ref_points_dict)

    # --- 2. 遍历生成每个骨架的独立图 (移除图例，颜色同步) ---
    for idx, label in enumerate(file_labels):
        if label in all_data:
            single_scaffold_dict = {label: all_data[label]}
            individual_path = output_path.parent / f"{output_path.stem}_{label}{output_path.suffix}"
            
            # 单个骨架图只显示对应的参考点
            single_ref_point = {label: ref_points_dict[label]} if label in ref_points_dict else None
            
            create_single_plot(
                single_scaffold_dict, 
                f"Scaffold {label}", 
                individual_path, 
                plt.get_cmap('tab10')(idx % 10), 
                False,
                ref_points_to_plot=single_ref_point
            )
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_files", nargs='+', type=str)
    # 增加 Dipole(D)
    parser.add_argument("--columns", nargs='+', default=["HOMO(eV)", "LUMO(eV)", "Gap(eV)", "Dipole(D)"])
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--width", type=int, default=14)
    parser.add_argument("--height", type=int, default=8)
    
    args = parser.parse_args()
    script_dir = Path(__file__).parent.absolute()
    pictos_labels = load_pictos_labels(script_dir.parent / "data" / "pictos.txt")
    
    input_files = [Path(f) for f in args.input_files if Path(f).exists()]
    if not input_files: return 1
    
    all_data = {}
    print("\n[文件读取与标签映射]")
    # 确保 S1-T1 和 charge 列也被读取（如果存在）
    columns_to_read = list(args.columns)
    if "S1-T1" not in columns_to_read:
        columns_to_read.append("S1-T1")
    if "charge" not in columns_to_read:
        columns_to_read.append("charge")
    
    for csv_file in input_files:
        label = get_label_from_filename(csv_file.name, pictos_labels)
        data = extract_data_from_csv(csv_file, columns_to_read)
        # 检查数据是否有效（使用 args.columns 计算，因为 S1-T1 可能为空）
        total_points = sum(len(data.get(col, [])) for col in args.columns)
        if total_points > 0:
            if label in all_data:
                for col in columns_to_read: all_data[label][col].extend(data[col])
            else: all_data[label] = data
            # 显示每个列的详细数据点数量（包括 S1-T1）
            col_info = ", ".join([f"{col}={len(data.get(col, []))}" for col in columns_to_read])
            print(f"  ✓ {csv_file.name} → {label} (总计: {total_points}, {col_info})")
        else:
            print(f"  ⚠ {csv_file.name} → {label} (无有效数据)")
            
    if not all_data: 
        print("错误: 没有有效数据")
        return 1
    sorted_labels = sort_labels_by_s_order(list(all_data.keys()))
    print(f"\n[最终标签列表] {sorted_labels}")
    
    # 检查每个标签的数据完整性
    print("\n[数据完整性检查]")
    for label in sorted_labels:
        if label in all_data:
            data_info = []
            for col in args.columns:
                count = len(all_data[label].get(col, []))
                data_info.append(f"{col}={count}")
            print(f"  {label}: {', '.join(data_info)}")
    
    base_out = Path(args.output) if args.output else script_dir / "frontier_orbital_energies"
    
    # 路径定义
    out_violin = base_out.parent / f"{base_out.stem}_violin.png"
    out_box = base_out.parent / f"{base_out.stem}_box.png"
    out_scatter = base_out.parent / f"{base_out.stem}_scatter.png"
    out_gap_dipole = base_out.parent / f"{base_out.stem}_gap_dipole.png"
    out_gap_s1t1 = base_out.parent / f"{base_out.stem}_gap_s1t1.png"
    out_charge_dipole = base_out.parent / f"{base_out.stem}_charge_dipole.png"

    # 绘图逻辑
    plot_violinplot(all_data, out_violin, args.columns, (args.width, args.height*0.83), sorted_labels,
                   pictos_labels=pictos_labels, script_dir=script_dir)
    plot_boxplot(all_data, out_box, args.columns, (args.width, args.height), sorted_labels)
    plot_scatter_homo_lumo(all_data, out_scatter, sorted_labels)
    plot_gap_dipole_scatter(all_data, out_gap_dipole, sorted_labels, 
                           pictos_labels=pictos_labels, script_dir=script_dir)
    plot_gap_s1t1_scatter(all_data, out_gap_s1t1, sorted_labels, 
                         pictos_labels=pictos_labels, script_dir=script_dir)
    plot_charge_dipole_scatter(all_data, out_charge_dipole, sorted_labels, 
                              pictos_labels=pictos_labels, script_dir=script_dir)
    
    # 保留原有打印
    print(f"成功生成：小提琴图、箱线图、带 1eV 投影的散点图。")
    
    # ---------------------------------------------------------
    # 新增：明确告知文件保存位置的 print
    # ---------------------------------------------------------
    print(f"\n[文件保存详情]")
    print(f"1. 小提琴图: {out_violin.resolve()}")
    print(f"2. 箱线图  : {out_box.resolve()}")
    print(f"3. 散点图  : {out_scatter.resolve()}")
    print(f"4. Gap-Dipole图: {out_gap_dipole.resolve()}")
    print(f"5. Gap-S1T1图: {out_gap_s1t1.resolve()}")
    print(f"6. Charge-Dipole图: {out_charge_dipole.resolve()}")
    
    return 0

if __name__ == "__main__":
    raise SystemExit(main())