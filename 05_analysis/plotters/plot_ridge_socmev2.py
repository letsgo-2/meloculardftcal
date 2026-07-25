#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成左右分栏的Ridge Plot（山脊图） - 统计逻辑修正版

更新说明：
1. 修正统计逻辑：Max和Median现在基于**全量有效数据**计算，不再受限于绘图范围。
2. 优化标注位置：如果真实Max超出显示范围，文字会吸附在最右侧，但显示真实数值。
3. 保持数据过滤：依然自动过滤 > 10000 的异常值。
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# 检查库依赖
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def load_pictos_labels(pictos_file: Path) -> Dict[str, str]:
    labels = {}
    if not pictos_file.exists(): return labels
    try:
        with open(pictos_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue
                if ':' in line:
                    parts = line.split(':', 1)
                    labels[parts[0].strip()] = parts[1].strip()
                else:
                    parts = line.split()
                    if len(parts) >= 2: labels[parts[0]] = parts[1]
    except Exception: pass
    return labels

def get_label_from_filename(filename: str, pictos_labels: Dict[str, str]) -> str:
    pic_name = None
    for pic in pictos_labels.keys():
        if filename.startswith(pic):
            pic_name = pic
            break
    if pic_name: return pictos_labels[pic_name]
    if filename.startswith('pic') and any(c.isdigit() for c in filename):
        match = filename.split('_')[0] if '_' in filename else filename.split('.')[0]
        if match in pictos_labels: return pictos_labels[match]
    return Path(filename).stem

def sort_labels_by_s_order(labels: List[str]) -> List[str]:
    def get_s_number(label):
        label_upper = label.upper()
        if label_upper.startswith('S'):
            try: return int(label_upper[1:])
            except: return 999
        return 999
    
    sorted_labels = sorted(labels, key=get_s_number)
    formatted_labels = []
    for label in sorted_labels:
        if label.lower().startswith('s') and label[1:].isdigit():
             formatted_labels.append(f"S{label[1:].upper()}")
        else:
             formatted_labels.append(label)
    return formatted_labels

def extract_data_from_csv(csv_file: Path, columns: List[str]) -> Dict[str, List[float]]:
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"警告: 无法读取 {csv_file}: {e}")
        return {}
    
    result = {}
    for col in columns:
        if col not in df.columns:
            result[col] = []
            continue
        valid_data = df[col].dropna().astype(float)
        # 过滤极端的非物理值 (>10000)
        valid_data = valid_data[valid_data <= 10000].tolist()
        result[col] = valid_data
    return result

def plot_split_ridge_plot(
    all_data: Dict[str, Dict[str, List[float]]],
    output_path: Path,
    columns: List[str],
    figsize: Tuple[int, int] = (12, 10),
    sort_order: Optional[List[str]] = None,
    overlap: float = 0.7,
    scale_factor: float = 0.4,
    x_labels: List[str] = None,
    xmax_custom: List[float] = None,
    ref_points: Optional[Dict[str, Dict[str, float]]] = None,
):
    if not all_data: return

    if len(columns) > 2:
        print("提示: 仅绘制前两列数据进行对比。")
        columns = columns[:2]
    
    if sort_order is not None:
        file_labels = [label for label in sort_order if label in all_data]
        for label in all_data.keys():
            if label not in file_labels: file_labels.append(label)
    else:
        file_labels = sort_labels_by_s_order(list(all_data.keys()))
    
    n_labels = len(file_labels)
    n_cols = len(columns)
    
    # 淡色系配色
    col_colors = ['#FFB3BA', '#BAE1FF'] 
    
    # -----------------------------------------------------------
    # 计算每一列的X轴范围
    # -----------------------------------------------------------
    col_limits = []
    
    for idx, col in enumerate(columns):
        all_vals = []
        for lab in file_labels:
            if col in all_data[lab]: all_vals.extend(all_data[lab][col])
        
        if not all_vals:
            col_limits.append((-0.5, 1))
            continue
            
        vmin = -0.5 
        p99 = np.percentile(all_vals, 99)
        vmax = max(vmin + 0.1, p99 * 1.2)
        
        # 应用自定义范围
        if xmax_custom and idx < len(xmax_custom):
            user_max = xmax_custom[idx]
            if user_max > 0:
                print(f"  列 '{col}' 使用自定义最大值: {user_max}")
                vmax = user_max
        else:
            print(f"  列 '{col}' 使用自动最大值: {vmax:.2f}")

        col_limits.append((vmin, vmax))

    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans']})
    fig = plt.figure(figsize=figsize)
    
    hspace = -overlap 
    # 调整 GridSpec 的 top 参数，减少上方空白
    gs = GridSpec(n_labels, n_cols, figure=fig, hspace=hspace, wspace=0.15, top=0.92, bottom=0.05, left=0.08, right=0.95)
    
    # 添加默认标题，位置在 GridSpec 下方
    fig.suptitle('SOCME Distribution Analysis', fontsize=12, fontweight='bold', y=0.86)
    
    axes_grid = []
    
    for row_idx, file_label in enumerate(file_labels):
        row_axes = []
        for col_idx, column in enumerate(columns):
            sharex = axes_grid[0][col_idx] if row_idx > 0 else None
            ax = fig.add_subplot(gs[row_idx, col_idx], sharex=sharex)
            row_axes.append(ax)
            
            if column in all_data[file_label] and len(all_data[file_label][column]) > 0:
                # 1. 获取全量数据 (已过滤 >10000)
                data_full = np.array(all_data[file_label][column])
                c_min, c_max = col_limits[col_idx]
                
                if len(data_full) > 0:
                    # ---------------------------------------------------
                    # 关键修改：统计数据基于全量数据计算
                    # ---------------------------------------------------
                    real_median_val = np.median(data_full)
                    real_max_val = np.max(data_full)

                    # ---------------------------------------------------
                    # 绘图数据：仅截取显示范围内的数据用于生成平滑曲线
                    # (避免极远处的点导致主峰被压缩成一条线)
                    # ---------------------------------------------------
                    data_for_plot = data_full[data_full <= c_max * 1.5]
                    
                    if len(data_for_plot) > 1:
                        color = col_colors[col_idx % len(col_colors)]
                        x_grid = np.linspace(c_min, c_max, 500)
                        
                        y_median_height = 0
                        y_max_density = 1.0  # 默认值，避免未定义
                        
                        if HAS_SCIPY:
                            try:
                                kde = stats.gaussian_kde(data_for_plot, bw_method=0.2)
                                y = kde(x_grid)
                                y_max_density = y.max()
                                if y_max_density > 0:
                                    y = y / y_max_density * scale_factor
                                    # 计算中位数高度 (如果中位数在绘图范围内)
                                    if real_median_val <= c_max:
                                        y_median_height = (kde(real_median_val) / y_max_density * scale_factor)[0]
                                    else:
                                        # 如果中位数都超出范围了，就画在最右边或不画高度
                                        y_median_height = 0.05
                                
                                ax.fill_between(x_grid, 0, y, alpha=0.5, color=color, zorder=10)
                                ax.plot(x_grid, y, color='black', linewidth=0.5, alpha=1.0, zorder=11)
                            except: 
                                y_max_density = 1.0  # 如果 KDE 失败，使用默认值
                        else:
                            counts, bins = np.histogram(data_for_plot, bins=50, density=True)
                            centers = 0.5*(bins[:-1]+bins[1:])
                            y_max_density = counts.max()
                            if y_max_density > 0: 
                                counts = counts / y_max_density * scale_factor
                                if real_median_val <= c_max:
                                    idx_median = (np.abs(centers - real_median_val)).argmin()
                                    y_median_height = counts[idx_median]

                            ax.fill_between(centers, 0, counts, alpha=0.5, color=color, step='mid', zorder=10)
                            ax.step(centers, counts, where='mid', color=color, linewidth=1.5, zorder=11)

                        # ---------------------------------------------------
                        # 绘制统计标注
                        # ---------------------------------------------------
                        
                        # 1. 中位数 (Median)
                        # 如果中位数在显示范围内，正常画线；如果超出，吸附在边缘
                        median_display_x = min(real_median_val, c_max)
                        
                        if real_median_val <= c_max:
                            ax.vlines(x=median_display_x, ymin=0, ymax=y_median_height, 
                                     color='blue', linestyle='-', linewidth=1.2, alpha=0.8, zorder=15)
                            # 计算右移距离（x轴范围的2%）
                            x_range = c_max - c_min
                            median_text_x = median_display_x + x_range * 0.02
                            ax.text(median_text_x, y_median_height + 0.05 * scale_factor, f'{real_median_val:.1f}', 
                                    ha='left', va='bottom', fontsize=11, color='blue', fontweight='bold', zorder=20)
                        
                        # 2. 最大值 (Max)
                        # 位置吸附逻辑：如果真实Max超出坐标轴，文字停在坐标轴最右侧，但显示真实数值
                        # 如果Max没有超出坐标系，标签右移一点以便更好显示
                        if real_max_val <= c_max:
                            # 没有超出：右移一点（移动距离为坐标轴范围的2%）
                            x_range = c_max - c_min
                            max_display_x = real_max_val + x_range * 0.1
                            # 确保不会超出太多
                            max_display_x = min(max_display_x, c_max * 1.05)
                        else:
                            # 超出坐标系：吸附在右边界
                            max_display_x = c_max
                        
                        # 稍微往上抬一点(0.02)避免和轴线重叠
                        ax.text(max_display_x, 0.02, f'Max:{real_max_val:.1f}', 
                                ha='right', va='bottom', fontsize=7, color='darkgreen', fontweight='bold', zorder=20)
                        
                        # 标记参考点（红色五角星）
                        if ref_points and file_label in ref_points:
                            ref_data = ref_points[file_label]
                            if column in ref_data:
                                ref_val = ref_data[column]
                                if ref_val is not None:
                                    # 计算参考点在 KDE 曲线上的 y 值
                                    ref_y = 0.02  # 默认值
                                    if HAS_SCIPY and len(data_for_plot) > 1:
                                        try:
                                            kde = stats.gaussian_kde(data_for_plot, bw_method=0.2)
                                            kde_val = kde(ref_val)[0]
                                            if y_max_density > 0:
                                                ref_y = (kde_val / y_max_density * scale_factor)
                                            else:
                                                ref_y = 0.02
                                        except:
                                            ref_y = 0.02
                                    else:
                                        ref_y = 0.02
                                    
                                    # 确保参考点在显示范围内
                                    if ref_val <= c_max:
                                        ax.scatter(ref_val, ref_y, s=160, color='#e74c3c', 
                                                  edgecolors='black', linewidth=0.8, marker='*', 
                                                  zorder=20)

            if col_idx == 0:
                ax.text(-0.02, 0.05, file_label, 
                        ha='right', va='center', fontsize=12, fontweight='bold', 
                        color='#4A4A4A', transform=ax.transAxes)
            
            ax.set_yticks([])
            ax.set_ylim(0, 1.0)
            ax.set_xlim(col_limits[col_idx])
            
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(True)
            ax.spines['bottom'].set_color('#4A4A4A')
            ax.spines['bottom'].set_linewidth(1.0)
            ax.patch.set_alpha(0)
            
            if row_idx < n_labels - 1:
                ax.tick_params(axis='x', labelbottom=False, bottom=False)
            else:
                ax.tick_params(axis='x', labelbottom=True, bottom=True, labelsize=10, color='#4A4A4A', labelcolor='#4A4A4A')
                label_text = x_labels[col_idx] if x_labels and col_idx < len(x_labels) else column
                ax.set_xlabel(label_text, fontsize=11, fontweight='bold', labelpad=10, color='#333333')

        axes_grid.append(row_axes)

    # -----------------------------------------------------------
    # 添加图例说明
    # -----------------------------------------------------------
    if axes_grid and axes_grid[0]:
        # 使用右上角的子图（第一行最后一列）来添加图例
        top_right_ax = axes_grid[0][-1]
        fig = top_right_ax.figure
        
        # 获取列标签
        label_socme = x_labels[0] if x_labels and len(x_labels) > 0 else "SOCME"
        label_socme_de = x_labels[1] if x_labels and len(x_labels) > 1 else "SOCME/ΔE"
        
        # 创建图例元素：蓝色线表示中位数，以及两列数据的颜色填充块
        legend_elements = [
            mlines.Line2D([0], [0], color='blue', linestyle='-', linewidth=1.5, 
                         alpha=0.8, label='Median'),
            mpatches.Patch(facecolor=col_colors[0], alpha=0.5, edgecolor='none', label=label_socme),
            mpatches.Patch(facecolor=col_colors[1], alpha=0.5, edgecolor='none', label=label_socme_de),
        ]
        
        # 创建图例
        top_right_ax.legend(handles=legend_elements, 
                           loc='upper right', 
                           bbox_to_anchor=(1.0, 0.6),
                           frameon=True, 
                           fancybox=True, 
                           shadow=False, 
                           fontsize=9,
                           framealpha=0.95, 
                           edgecolor='#CCCCCC',
                           facecolor='white')

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # 添加黑色边框
    # fig.patch.set_edgecolor('black')
    # fig.patch.set_linewidth(2.0)
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='black', pad_inches=0.1)
    print(f"分栏Ridge plot已保存至: {output_path}")
    plt.close()

def main() -> int:
    parser = argparse.ArgumentParser(description="生成双栏 Ridge Plot（带准确统计标注）")
    parser.add_argument("input_files", nargs='+', type=str, help="输入CSV文件")
    parser.add_argument("--columns", nargs='+', default=["max_SOCME_S1_Ti_value", "max_SOCME_over_deltaE_S1_Ti_value"], help="列名")
    parser.add_argument("--output", "-o", type=str, default=None, help="输出路径")
    parser.add_argument("--width", type=int, default=12, help="宽")
    parser.add_argument("--height", type=int, default=10, help="高")
    parser.add_argument("--overlap", type=float, default=0.5, help="重叠比例")
    parser.add_argument("--scale", type=float, default=0.4, help="山峰高度 (0.1-1.0)")
    parser.add_argument("--xlabels", nargs='+', default=["SOCME (cm$^{-1}$)", "SOCME/ΔE "], help="两列的X轴标签")
    parser.add_argument("--xmax", nargs='+', type=float, default=None, help="手动指定X轴最大值，例如: --xmax 50 0.8")
    
    args = parser.parse_args()
    
    script_dir = Path(__file__).parent.absolute()
    pictos_file = script_dir.parent / "data" / "pictos.txt"
    pictos_labels = load_pictos_labels(pictos_file)
    
    input_files = [Path(f) for f in args.input_files if Path(f).exists()]
    if not input_files: return 1
        
    all_data = {}
    for csv_file in input_files:
        file_label = get_label_from_filename(csv_file.name, pictos_labels)
        if file_label.lower().startswith('s'): file_label = f"S{file_label[1:].upper()}"
        
        data = extract_data_from_csv(csv_file, args.columns)
        if file_label in all_data:
            for col in args.columns:
                if col in data: all_data[file_label][col].extend(data[col])
        else:
            all_data[file_label] = data
            
    has_data = any(len(d.get(c, [])) > 0 for d in all_data.values() for c in args.columns)
    if not has_data: return 1
    
    # 打印每个骨架的统计信息
    print("\n" + "="*80)
    print("每个骨架的统计信息 (平均值和最大值):")
    print("="*80)
    
    # 获取列标签（用于显示）
    col_labels = args.xlabels if args.xlabels else args.columns
    
    # 按 S 顺序排序
    sorted_labels = sort_labels_by_s_order(list(all_data.keys()))
    
    for label in sorted_labels:
        if label not in all_data:
            continue
        
        print(f"\n{label}:")
        for idx, col in enumerate(args.columns):
            if col in all_data[label] and len(all_data[label][col]) > 0:
                data = np.array(all_data[label][col])
                mean_val = np.mean(data)
                max_val = np.max(data)
                col_display = col_labels[idx] if idx < len(col_labels) else col
                print(f"  {col_display}:")
                print(f"    平均值 (Mean): {mean_val:.4f}")
                print(f"    最大值 (Max):  {max_val:.4f}")
            else:
                col_display = col_labels[idx] if idx < len(col_labels) else col
                print(f"  {col_display}: 无数据")
    
    print("\n" + "="*80 + "\n")
    
    # 加载参考点数据
    ref_points_dict: Dict[str, Dict[str, float]] = {}  # {batch_label: {'max_SOCME_S1_Ti_value': float, 'max_SOCME_over_deltaE_S1_Ti_value': float}}
    publish_file = script_dir.parent / "data" / "publishinpdftv22.csv"
    if publish_file.exists():
        try:
            import re
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
                
                # 提取 max_SOCME_S1_Ti_value 和 max_SOCME_over_deltaE_S1_Ti_value
                socme_val = row.get('max_SOCME_S1_Ti_value', np.nan)
                socme_delta_val = row.get('max_SOCME_over_deltaE_S1_Ti_value', np.nan)
                
                ref_data = {}
                if pd.notna(socme_val):
                    try:
                        ref_data['max_SOCME_S1_Ti_value'] = float(socme_val)
                    except (ValueError, TypeError):
                        pass
                if pd.notna(socme_delta_val):
                    try:
                        ref_data['max_SOCME_over_deltaE_S1_Ti_value'] = float(socme_delta_val)
                    except (ValueError, TypeError):
                        pass
                
                if ref_data:
                    ref_points_dict[batch_label] = ref_data
                    print(f"  ✅ 已加载参考点 {index_val} -> {batch_label}: {ref_data}")
        except Exception as e:
            print(f"  ⚠️ 加载参考点数据时出错: {e}")
        
    output_path = Path(args.output) if args.output else script_dir / "socme_ridge_stats2.png"
    
    plot_split_ridge_plot(
        all_data,
        output_path,
        args.columns,
        figsize=(args.width, args.height),
        overlap=args.overlap,
        scale_factor=args.scale,
        x_labels=args.xlabels,
        xmax_custom=args.xmax,
        ref_points=ref_points_dict,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())