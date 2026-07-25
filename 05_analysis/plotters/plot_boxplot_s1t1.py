#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S1和T1能量值小提琴图绘制脚本

功能说明：
    从多个CSV文件中提取S1_energy_eV和T1_energy_eV数据，生成小提琴图（曲线箱体图）。
    支持多个文件输入，将所有文件的数据绘制在同一张图中进行比较。
    小提琴图结合了箱体图和核密度估计，能更好地展示数据分布形状。

主要特性：
    1. 支持多个CSV文件同时输入
    2. 提取S1_energy_eV和T1_energy_eV列数据
    3. 为每个文件每个列生成小提琴图（曲线箱体图），放在同一张图中进行对比
    4. 自动处理缺失值和无效数据
    5. 输出文件保存在脚本所在目录
    6. 使用pictos.txt中的标签映射，按S顺序排列
    7. 小提琴图显示数据分布曲线，更直观地展示数据分布特征

输入文件格式要求：
    - CSV文件必须包含以下列：S1_energy_eV, T1_energy_eV
    - 列名必须准确匹配（大小写敏感）

使用示例：
    1. 单个文件：
       python3 plot/plot_boxplot_s1t1.py data/pic3inpout_m.csv
    
    2. 多个文件：
       python3 plot/plot_boxplot_s1t1.py data/pic3inpout_m.csv data/pic9inpout_m.csv data/pic16inpout_m.csv
    
    3. 指定输出文件名：
       python3 plot/plot_boxplot_s1t1.py data/*.csv --output my_s1t1_boxplot.png
    
    4. 自定义图片大小：
       python3 plot/plot_boxplot_s1t1.py data/*.csv --width 16 --height 10
    
    5. 只绘制S1能量：
       python3 plot/plot_boxplot_s1t1.py data/*.csv --columns S1_energy_eV

输出说明：
    - 默认输出文件名：s1t1_energy_boxplot.png（保存在脚本所在目录）
    - 图片包含多个文件的小提琴图，便于对比分析
"""
from __future__ import annotations

from typing import List, Dict, Optional
import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines


def load_pictos_labels(pictos_file: Path) -> Dict[str, str]:
    """
    从pictos.txt文件中加载pic到标签的映射
    
    Args:
        pictos_file: pictos.txt文件路径
    
    Returns:
        字典，键为pic名称（如"pic3"），值为标签（如"S9"）
    """
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
    """
    从文件名中提取pic名称，并返回对应的标签（统一为大写S格式）
    
    Args:
        filename: 文件名（如"pic3inpout_m.csv"）
        pictos_labels: pic名称到标签的映射
    
    Returns:
        标签（如"S9"），统一为大写S格式，如果找不到则返回文件名
    """
    # 提取pic名称（如"pic3inpout_m.csv" -> "pic3"）
    pic_name = None
    for pic in pictos_labels.keys():
        if filename.startswith(pic):
            pic_name = pic
            break
    
    if pic_name and pic_name in pictos_labels:
        label = pictos_labels[pic_name]
        # 统一格式为S1, S2等（大写S）
        if label.lower().startswith('s'):
            return f"S{label[1:].upper()}"
        return label
    
    # 如果找不到，尝试从文件名提取（如"pic3inpout_m" -> "pic3"）
    if filename.startswith('pic') and any(c.isdigit() for c in filename):
        match = filename.split('_')[0] if '_' in filename else filename.split('.')[0]
        if match.startswith('pic') and match in pictos_labels:
            label = pictos_labels[match]
            if label.lower().startswith('s'):
                return f"S{label[1:].upper()}"
            return label
    
    # 如果都找不到，返回文件名（不含扩展名）
    return Path(filename).stem


def sort_labels_by_s_order(labels: List[str]) -> List[str]:
    """
    按照S的顺序排序标签（S1, S2, S3, ..., S9）
    统一转换为大写S
    
    Args:
        labels: 标签列表
    
    Returns:
        排序后的标签列表（统一为大写S）
    """
    def get_s_number(label):
        # 提取S后面的数字，统一处理大小写
        label_upper = label.upper()
        if label_upper.startswith('S'):
            try:
                return int(label_upper[1:])
            except:
                return 999  # 无法解析的放在最后
        return 999
    
    # 排序并统一转换为大写S格式
    sorted_labels = sorted(labels, key=get_s_number)
    # 统一格式为S1, S2等（大写S）
    formatted_labels = []
    for label in sorted_labels:
        if label.lower().startswith('s'):
            formatted_labels.append(f"S{label[1:].upper()}")
        else:
            formatted_labels.append(label)
    
    return formatted_labels


def extract_data_from_csv(csv_file: Path, columns: List[str]) -> Dict[str, List[float]]:
    """
    从CSV文件中提取指定列的数据
    
    Args:
        csv_file: CSV文件路径
        columns: 要提取的列名列表
    
    Returns:
        字典，键为列名，值为该列的有效数据列表
    """
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"警告: 无法读取文件 {csv_file}: {e}")
        return {}
    
    result = {}
    for col in columns:
        if col not in df.columns:
            print(f"警告: 文件 {csv_file.name} 中找不到列 '{col}'")
            result[col] = []
            continue
    
    # 提取有效数据（非NaN，非空）
        valid_data = df[col].dropna().astype(float)
        
        # 如果列名是 T1_energy_eV，剔除负数值
        if col == 'T1_energy_eV':
            valid_data = valid_data[valid_data >= 0]
            if len(valid_data) < len(df[col].dropna()):
                removed_count = len(df[col].dropna()) - len(valid_data)
                print(f"    注意: {col} 中剔除了 {removed_count} 个负数值")
        
        result[col] = valid_data.tolist()
    
    return result


def _draw_energy_violin(
    ax,
    all_data: Dict[str, Dict[str, List[float]]],
    columns: List[str],
    sort_order: Optional[List[str]] = None,
    ref_points: Optional[Dict[str, Dict[str, float]]] = None,
):
    """
    在给定的 ax 上绘制 S1/T1 能量小提琴图（不创建/保存 figure）
    """
    n_columns = len(columns)
    # 如果提供了排序顺序，使用它；否则按照S的顺序排序
    if sort_order is not None:
        file_labels = [label for label in sort_order if label in all_data]
        # 添加不在排序列表中的标签
        for label in all_data.keys():
            if label not in file_labels:
                file_labels.append(label)
    else:
        file_labels = sort_labels_by_s_order(list(all_data.keys()))
    
    # 准备数据：每个文件每个列都要绘制
    # 位置计算：每个文件有两个位置（一个S1_energy_eV，一个T1_energy_eV）
    n_files = len(file_labels)
    width_between_files = 3  # 文件之间的间距
    width_between_columns = 1.0  # 同一文件内两列之间的间距
    
    all_positions: List[float] = []
    all_data_list: List[List[float]] = []
    all_labels: List[str] = []
    all_colors: List = []
    
    # 为每个列定义颜色（使用淡色系）
    column_colors = ['#FFB3BA', '#BAE1FF']  # S1用淡粉色，T1用淡蓝色
    if len(columns) > 2:
        # 如果有更多列，使用淡色颜色映射
        column_colors = plt.cm.Pastel1(np.linspace(0, 1, len(columns)))
    
    for file_idx, file_label in enumerate(file_labels):
        base_position = file_idx * width_between_files + 1
        
        for col_idx, column in enumerate(columns):
            if column in all_data[file_label] and len(all_data[file_label][column]) > 0:
                position = base_position + col_idx * width_between_columns
                all_positions.append(position)
                all_data_list.append(all_data[file_label][column])
                # 标签格式：文件名_列名 或 列名
                all_labels.append(f"{file_label}\n{column}")
                all_colors.append(column_colors[col_idx] if isinstance(column_colors[col_idx], str) else column_colors[col_idx])
    
    if not all_data_list:
        ax.text(0.5, 0.5, 'No data available', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Energy Distribution', fontsize=12, fontweight='bold')
        return [], []
    
    # 绘制所有小提琴图
    parts = ax.violinplot(
        all_data_list,
        positions=all_positions,
        showmeans=True,
        showmedians=True,
        widths=0.6
    )
    
    # 设置小提琴图颜色（使用更淡的透明度，与plot_boxplot_s1_t1.py保持一致）
    for pc, color in zip(parts['bodies'], all_colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.5)  # 降低透明度，使颜色更淡
        pc.set_edgecolor('black')
        pc.set_linewidth(1.0)  # 稍微减小边框宽度
    
    # 设置其他部分（均值线、中位数线等）的颜色和样式（与plot_boxplot_s1_t1.py保持一致）
    # cmeans: 均值线（LineCollection），cmedians: 中位数线（LineCollection）
    if 'cmeans' in parts:
        cmeans = parts['cmeans']
        # cmeans 是 LineCollection 对象
        cmeans.set_color('red')  # 均值线用红色，更明显
        cmeans.set_linewidth(2.5)  # 增加线宽
        cmeans.set_linestyle('--')  # 虚线样式
    
    if 'cmedians' in parts:
        cmedians = parts['cmedians']
        # cmedians 是 LineCollection 对象
        cmedians.set_color('blue')  # 中位数线用蓝色
        cmedians.set_linewidth(2.5)  # 增加线宽
        cmedians.set_linestyle('-')  # 实线样式
    
    # 其他线（cbars: 中位数延伸线, cmins/cmaxes: 最小/最大值线）
    for partname in ('cbars', 'cmins', 'cmaxes'):
        if partname in parts:
            part = parts[partname]
            # 这些也可能是列表
            if isinstance(part, list):
                for line in part:
                    line.set_color('black')
                    line.set_linewidth(1.0)
            else:
                part.set_color('black')
                part.set_linewidth(1.0)
    
    # 设置x轴标签和刻度
    # 在每个文件组的中心位置显示文件名
    file_center_positions = []
    file_center_labels = []
    for file_idx, file_label in enumerate(file_labels):
        base_position = file_idx * width_between_files + 1
        center_position = base_position + (len(columns) - 1) * width_between_columns / 2
        file_center_positions.append(center_position)
        file_center_labels.append(file_label)
    
    ax.set_xticks(file_center_positions)
    ax.set_xticklabels(file_center_labels, rotation=45, ha='right')
    
    # 添加次要刻度线标记每个数据点（可选）
    ax.set_xticks(all_positions, minor=True)
    ax.tick_params(axis='x', which='minor', length=4)
    
    # 在每个小提琴图上标注均值和最小值
    for i, (data, position) in enumerate(zip(all_data_list, all_positions)):
        if len(data) > 0:
            mean_val = np.mean(data)
            min_val = np.min(data)
            
            # 获取y轴范围，用于确定标注位置
            y_max = ax.get_ylim()[1]
            
            # 在均值位置上方标注均值
            ax.text(position, mean_val + (y_max - mean_val) * 0.05, 
                   f'Mean: {mean_val:.3f}',
                   ha='center', va='bottom', fontsize=8, color='red', fontweight='bold')
            
            # 在最小值位置标注最小值
            ax.text(position, min_val - (min_val - ax.get_ylim()[0]) * 0.05,
                   f'Min: {min_val:.3f}',
                   ha='center', va='top', fontsize=8, color='darkgreen', fontweight='bold')
    
    
    # 设置标题和标签
    # 格式化列名为显示名称
    def format_column_name(col_name):
        if col_name == 'S1_energy_eV':
            return 'S1 Energy'
        elif col_name == 'T1_energy_eV':
            return 'T1 Energy'
        return col_name
    
    formatted_columns = [format_column_name(col) for col in columns]
    title = ' vs '.join(formatted_columns)
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_ylabel('Energy (eV)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Scaffold', fontsize=12, fontweight='bold')
    
    # 添加参考线，并在旁边标注数字
    reference_lines = [0.98, 2.0, 2.5, 3.26]
    reference_labels = ['0.98', '2.0', '2.5', '3.26']
    
    # 获取x轴范围，用于确定标注位置
    x_min, x_max = ax.get_xlim()
    x_range = x_max - x_min
    
    for ref_val, ref_label in zip(reference_lines, reference_labels):
        # 绘制参考线
        ax.axhline(y=ref_val, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, zorder=2)
        # 在参考线右侧添加文字标注
        ax.text(x_max - x_range * 0.02, ref_val, f'{ref_label} eV', 
                ha='right', va='center', fontsize=9, color='gray', 
                fontweight='bold', bbox=dict(boxstyle='round,pad=0.2', 
                facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5),
                zorder=3)
    
    # 添加图例：包括数据系列、均值和中位数线
    legend_elements = []
    # 数据系列图例
    for col_idx, column in enumerate(columns):
        color = column_colors[col_idx] if isinstance(column_colors[col_idx], str) else column_colors[col_idx]
        display_name = format_column_name(column)
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.5, edgecolor='black', label=display_name))
    
    # 统计线图例
    legend_elements.append(mlines.Line2D([0], [0], color='red', linestyle='--', linewidth=2, label='Mean'))
    legend_elements.append(mlines.Line2D([0], [0], color='blue', linestyle='-', linewidth=2, label='Median'))
    
    ax.legend(handles=legend_elements, loc='upper right', framealpha=0.9)
    
    # 添加网格
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # 标记参考点（红色五角星）
    if ref_points:
        for file_idx, file_label in enumerate(file_labels):
            if file_label in ref_points:
                base_position = file_idx * width_between_files + 1
                ref_data = ref_points[file_label]
                
                for col_idx, column in enumerate(columns):
                    if column in ref_data:
                        ref_val = ref_data[column]
                        if ref_val is not None:
                            position = base_position + col_idx * width_between_columns
                            # 稍微侧移（向左偏移 0.2）
                            ax.scatter(position - 0.2, ref_val, s=160, color='#e74c3c', 
                                      edgecolors='black', linewidth=0.8, marker='*', 
                                      zorder=10)
    
    # 返回每个骨架中心位置和对应标签，供振子强度共用
    return file_center_positions, file_center_labels


def plot_boxplot(
    all_data: Dict[str, Dict[str, List[float]]],
    output_path: Path,
    columns: List[str],
    figsize: tuple = (14, 8),
    sort_order: Optional[List[str]] = None,
    osc_data: Optional[Dict[str, List[float]]] = None,
    osc_label: str = "Oscillator Strength",
    ref_points: Optional[Dict[str, Dict[str, float]]] = None,
):
    """
    在一个图中绘制 S1/T1 能量小提琴图；
    如提供 osc_data，则在同一张图的右侧用第二个 y 轴绘制振子强度分布（小提琴），
    并在右侧显示振子强度的刻度和 y 轴标签。
    """
    if not all_data:
        print("错误: 没有有效数据可绘制")
        return

    fig, ax_energy = plt.subplots(1, 1, figsize=figsize)

    # 能量图，返回每个骨架的中心位置 & 标签
    file_center_positions, file_center_labels = _draw_energy_violin(
        ax_energy, all_data, columns, sort_order=sort_order, ref_points=ref_points
    )

    # 在同一个图中，用右侧 y 轴绘制振子强度（如果提供）
    enable_osc = osc_data is not None and len(osc_data) > 0
    if enable_osc and file_center_positions and file_center_labels:
        # 按能量图的骨架顺序取出振子强度数据
        osc_values_list: List[List[float]] = []
        valid_positions: List[float] = []
        valid_labels: List[str] = []
        for pos, label in zip(file_center_positions, file_center_labels):
            vals = osc_data.get(label, [])
            if vals:
                osc_values_list.append(vals)
                valid_positions.append(pos)
                valid_labels.append(label)

        if osc_values_list:
            ax_osc = ax_energy.twinx()
            parts_osc = ax_osc.violinplot(
                osc_values_list,
                positions=valid_positions,
                showmeans=True,
                showmedians=True,
                widths=0.5,
            )
            base_color = "#C2E3FF"
            for pc in parts_osc["bodies"]:
                pc.set_facecolor(base_color)
                pc.set_alpha(0.4)
                pc.set_edgecolor("black")
                pc.set_linewidth(0.8)

            if "cmeans" in parts_osc:
                parts_osc["cmeans"].set_color("darkred")
                parts_osc["cmeans"].set_linewidth(1.5)
                parts_osc["cmeans"].set_linestyle("--")
            if "cmedians" in parts_osc:
                parts_osc["cmedians"].set_color("darkblue")
                parts_osc["cmedians"].set_linewidth(1.5)
                parts_osc["cmedians"].set_linestyle("-")
            for partname in ("cbars", "cmins", "cmaxes"):
                if partname in parts_osc:
                    parts_osc[partname].set_color("gray")
                    parts_osc[partname].set_linewidth(0.8)

            # 右侧 y 轴标签和刻度（即振子强度数值）
            ax_osc.set_ylabel(osc_label, fontsize=12, fontweight="bold")
            ax_osc.tick_params(axis="y", labelsize=10)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    if enable_osc:
        print(f"能量 + 振子强度小提琴图已保存至: {output_path}")
    else:
        print(f"小提琴图已保存至: {output_path}")
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="从多个CSV文件中提取S1_energy_eV和T1_energy_eV数据，生成小提琴图（曲线箱体图）进行对比",
        epilog="""
使用示例：
  单个文件：
    %(prog)s data/pic3inpout_m.csv
  
  多个文件：
    %(prog)s data/pic3inpout_m.csv data/pic9inpout_m.csv data/pic16inpout_m.csv
  
  指定输出文件：
    %(prog)s data/*.csv --output my_s1t1_boxplot.png
  
  自定义图片大小：
    %(prog)s data/*.csv --width 16 --height 10
  
  只绘制S1能量：
    %(prog)s data/*.csv --columns S1_energy_eV

注意：
  - CSV文件必须包含列：S1_energy_eV, T1_energy_eV
  - 输出文件默认保存在脚本所在目录
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "input_files",
        nargs='+',
        type=str,
        help="输入CSV文件路径（可指定多个文件）"
    )
    parser.add_argument(
        "--columns",
        nargs='+',
        default=["S1_energy_eV", "T1_energy_eV"],
        help="要提取和绘制的列名（默认: S1_energy_eV T1_energy_eV）"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="输出图片路径（可选）。如果不指定，默认保存在脚本所在目录，文件名为s1t1_energy_boxplot.png"
    )
    parser.add_argument(
        "--width",
        type=int,
        default=14,
        help="图片宽度，单位：英寸（默认: 14）"
    )
    parser.add_argument(
        "--height",
        type=int,
        default=8,
        help="图片高度，单位：英寸（默认: 8）"
    )
    parser.add_argument(
        "--enable-osc",
        action="store_true",
        help="是否额外绘制振子强度小提琴图（默认不绘制）",
    )
    parser.add_argument(
        "--osc-column",
        type=str,
        default="S1_fosc_D2",
        help="振子强度所在列名（默认: Oscillator_strength）",
    )
    
    args = parser.parse_args()
    
    # 加载pictos标签映射
    script_dir = Path(__file__).parent.absolute()
    pictos_file = script_dir.parent / "data" / "pictos.txt"
    pictos_labels = load_pictos_labels(pictos_file)
    
    if pictos_labels:
        print(f"已加载 {len(pictos_labels)} 个pic标签映射")
    
    # 处理输入文件
    input_files = []
    for file_path in args.input_files:
        path = Path(file_path)
        if not path.exists():
            print(f"警告: 文件不存在，跳过: {path}")
            continue
        input_files.append(path)
    
    if not input_files:
        print("错误: 没有找到有效的输入文件")
        return 1
    
    print(f"找到 {len(input_files)} 个有效输入文件:")
    for f in input_files:
        print(f"  - {f}")
    
    # 提取数据
    print(f"\n正在从文件中提取数据...")
    all_data = {}
    osc_all_data: Dict[str, List[float]] = {}
    file_to_label = {}  # 记录文件到标签的映射
    
    for csv_file in input_files:
        # 使用pictos标签，如果没有则使用文件名
        file_label = get_label_from_filename(csv_file.name, pictos_labels)
        file_to_label[csv_file.name] = file_label
        print(f"  处理: {csv_file.name} -> {file_label}")
        
        data = extract_data_from_csv(csv_file, args.columns)
        
        # 如果该标签已存在（可能多个文件映射到同一标签），合并数据
        if file_label in all_data:
            print(f"    注意: 标签 {file_label} 已存在，合并数据")
            for col in args.columns:
                if col in data:
                    all_data[file_label][col].extend(data[col])
        else:
            all_data[file_label] = data
        
        # 如果需要绘制振子强度，单独提取一列
        if args.enable_osc and args.osc_column:
            osc_dict = extract_data_from_csv(csv_file, [args.osc_column])
            osc_values = osc_dict.get(args.osc_column, [])
            if osc_values:
                if file_label in osc_all_data:
                    osc_all_data[file_label].extend(osc_values)
                else:
                    osc_all_data[file_label] = list(osc_values)
        
        # 统计每个列的数据量
        for col in args.columns:
            count = len(data.get(col, []))
        if count > 0:
                print(f"    {col}: {count} 个有效数据点")
        else:
                print(f"    {col}: 无有效数据")
    
    # 检查是否有有效能量数据
    has_energy_data = False
    for file_label, data in all_data.items():
        for col in args.columns:
            if col in data and len(data[col]) > 0:
                has_energy_data = True
                break
        if has_energy_data:
            break
    
    if not has_energy_data:
        print("错误: 所有文件中都没有找到有效能量数据")
        return 1
    
    # 确定输出路径
    if args.output:
        output_path = Path(args.output)
    else:
        # 获取脚本所在目录
        script_dir = Path(__file__).parent.absolute()
        output_path = script_dir / "s1t1_energy_boxplot.png"
    
    # 按照S的顺序排序标签
    sorted_labels = sort_labels_by_s_order(list(all_data.keys()))
    print(f"\n标签排序（按S顺序）: {sorted_labels}")
    
    # 加载参考点数据
    ref_points_dict: Dict[str, Dict[str, float]] = {}  # {batch_label: {'S1_energy_eV': float, 'T1_energy_eV': float}}
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
                
                # 提取 S1_energy_eV 和 T1_energy_eV
                s1_val = row.get('S1_energy_eV', np.nan)
                t1_val = row.get('T1_energy_eV', np.nan)
                
                ref_data = {}
                if pd.notna(s1_val):
                    try:
                        ref_data['S1_energy_eV'] = float(s1_val)
                    except (ValueError, TypeError):
                        pass
                if pd.notna(t1_val):
                    try:
                        ref_data['T1_energy_eV'] = float(t1_val)
                    except (ValueError, TypeError):
                        pass
                
                if ref_data:
                    ref_points_dict[batch_label] = ref_data
                    print(f"  ✅ 已加载参考点 {index_val} -> {batch_label}: {ref_data}")
        except Exception as e:
            print(f"  ⚠️ 加载参考点数据时出错: {e}")
    
    # 绘制能量 + 振子强度小提琴图（在同一张图里）
    print(f"\n正在生成小提琴图（曲线箱体图）...")
    if args.enable_osc and osc_all_data:
        print("  已开启振子强度绘制，将在同一张图中显示能量和振子强度。")
        plot_boxplot(
            all_data,
            output_path,
            args.columns,
            figsize=(args.width, args.height),
            sort_order=sorted_labels,
            osc_data=osc_all_data,
            osc_label=args.osc_column,
            ref_points=ref_points_dict,
        )
    else:
        if args.enable_osc and not osc_all_data:
            print("  提示: 开启了振子强度绘制，但没有找到有效振子强度数据，仅绘制能量图。")
        plot_boxplot(
            all_data,
            output_path,
            args.columns,
            figsize=(args.width, args.height),
            sort_order=sorted_labels,
            osc_data=None,
            ref_points=ref_points_dict,
        )
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

