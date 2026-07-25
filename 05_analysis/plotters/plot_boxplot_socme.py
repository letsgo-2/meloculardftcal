#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SOCME 小提琴图绘制脚本

功能说明：
    从多个CSV文件中提取SOCME相关数据，生成小提琴图（曲线箱体图）。
    支持多个文件输入，将所有文件的数据绘制在同一张图中进行比较。
    小提琴图结合了箱体图和核密度估计，能更好地展示数据分布形状。

主要特性：
    1. 支持多个CSV文件同时输入
    2. 提取max_SOCME_S1_Ti_value和max_SOCME_over_deltaE_S1_Ti_value两列数据
    3. 为每个文件生成小提琴图（曲线箱体图），放在同一张图中进行对比
    4. 自动处理缺失值和无效数据
    5. 输出文件保存在脚本所在目录
    6. 小提琴图显示数据分布曲线，更直观地展示数据分布特征
    7. 使用pictos.txt中的标签映射，按S顺序排列

输入文件格式要求：
    - CSV文件必须包含以下列：max_SOCME_S1_Ti_value, max_SOCME_over_deltaE_S1_Ti_value
    - 列名必须准确匹配（大小写敏感）

使用示例：
    1. 单个文件：
       python3 plot_boxplot_socme.py data/pic3inpout_m.csv
    
    2. 多个文件：
       python3 plot_boxplot_socme.py data/pic3inpout_m.csv data/pic9inpout_m.csv data/pic16inpout_m.csv
    
    3. 指定输出文件名：
       python3 plot_boxplot_socme.py data/*.csv --output my_socme_boxplot.png
    
    4. 只绘制SOCME：
       python3 plot_boxplot_socme.py data/*.csv --columns max_SOCME_S1_Ti_value
    
    5. 自定义图片大小：
       python3 plot_boxplot_socme.py data/*.csv --width 16 --height 10

输出说明：
    - 默认输出文件名：socme_boxplot.png（保存在脚本所在目录）
    - 图片包含多个文件的箱体图，便于对比分析
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
        
        # 提取有效数据（非NaN，非空），并过滤掉大于10000的值
        valid_data = df[col].dropna().astype(float)
        # 过滤掉大于10000的值
        valid_data = valid_data[valid_data <= 10000].tolist()
        result[col] = valid_data
    
    return result


def plot_boxplot(
    all_data: Dict[str, Dict[str, List[float]]],
    output_path: Path,
    columns: List[str],
    figsize: tuple = (14, 8),
    sort_order: Optional[List[str]] = None,
    ref_points: Optional[Dict[str, Dict[str, float]]] = None,
):
    """
    绘制小提琴图（曲线箱体图）
    
    Args:
        all_data: 嵌套字典，结构为 {标签: {列名: [数据]}}
        output_path: 输出图片路径
        columns: 要绘制的列名列表
        figsize: 图片大小
        sort_order: 标签的排序顺序（如果提供，将按照此顺序排列）
    """
    n_columns = len(columns)
    n_files = len(all_data)
    
    if n_files == 0:
        print("错误: 没有有效数据可绘制")
        return
    
    # 创建单个图
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
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
    # 位置计算：每个文件有两个位置（一个列1，一个列2）
    n_files = len(file_labels)
    width_between_files = 3  # 文件之间的间距
    width_between_columns = 1.0  # 同一文件内两列之间的间距
    
    all_positions = []
    all_data_list = []
    all_labels = []
    all_colors = []
    
    # 为每个列定义颜色（使用淡色系）
    column_colors = ['#FFB3BA', '#BAE1FF']  # 第一个用淡粉色，第二个用淡蓝色
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
                all_labels.append(f"{file_label}\n{column}")
                all_colors.append(column_colors[col_idx] if isinstance(column_colors[col_idx], str) else column_colors[col_idx])
    
    if not all_data_list:
        ax.text(0.5, 0.5, 'No data available', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Data Comparison', fontsize=12, fontweight='bold')
        plt.tight_layout()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # 绘制所有小提琴图
    parts = ax.violinplot(
        all_data_list,
        positions=all_positions,
        showmeans=True,
        showmedians=True,
        widths=0.6
    )
    
    # 设置小提琴图颜色（使用更淡的透明度）
    for pc, color in zip(parts['bodies'], all_colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.5)  # 降低透明度，使颜色更淡
        pc.set_edgecolor('black')
        pc.set_linewidth(1.0)  # 稍微减小边框宽度
    
    # 设置其他部分（均值线、中位数线等）的颜色和样式
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
    
    # 设置y轴范围（在标注之前设置，以便标注位置正确）
    ax.set_ylim(bottom=0, top=500)
    
    # 在每个小提琴图上标注均值和最大值
    for i, (data, position) in enumerate(zip(all_data_list, all_positions)):
        if len(data) > 0:
            mean_val = np.mean(data)
            max_val = np.max(data)
            
            # 获取y轴范围
            y_max = ax.get_ylim()[1]  # 500
            y_min = ax.get_ylim()[0]  # 0
            y_range = y_max - y_min
            
            # 计算标注位置的偏移量（使用y轴范围的5%作为偏移）
            offset = y_range * 0.05
            
            # 标注均值：如果均值接近顶部，标注放在下方；否则放在上方
            if mean_val > y_max - offset * 3:  # 太接近顶部
                mean_text_y = mean_val - offset
                mean_va = 'top'
            else:
                mean_text_y = mean_val + offset
                mean_va = 'bottom'
            
            # 确保标注在图表范围内
            mean_text_y = max(y_min + offset * 0.5, min(y_max - offset * 0.5, mean_text_y))
            
            ax.text(position, mean_text_y, 
                   f'Mean: {mean_val:.3f}',
                   ha='center', va=mean_va, fontsize=8, color='red', fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='none'))
            
            # 标注最大值：如果最大值接近顶部，标注放在下方；否则放在上方
            if max_val > y_max - offset * 3:  # 太接近顶部
                max_text_y = max_val - offset
                max_va = 'top'
            else:
                max_text_y = max_val + offset
                max_va = 'bottom'
            
            # 确保标注在图表范围内，并且不与均值标注重叠
            max_text_y = max(y_min + offset * 0.5, min(y_max - offset * 0.5, max_text_y))
            
            # 如果最大值标注与均值标注太接近，调整位置
            if abs(max_text_y - mean_text_y) < offset * 0.8:
                if max_val > mean_val:
                    # 最大值在均值上方，将最大值标注再往上移
                    max_text_y = max(mean_text_y + offset * 1.5, min(y_max - offset * 0.5, max_val + offset))
                else:
                    # 最大值在均值下方，将最大值标注再往下移
                    max_text_y = min(mean_text_y - offset * 1.5, max(y_min + offset * 0.5, max_val - offset))
            
            ax.text(position, max_text_y,
                   f'Max: {max_val:.3f}',
                   ha='center', va=max_va, fontsize=8, color='darkgreen', fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # 设置标题和标签
    # 格式化列名为显示名称
    def format_column_name(col_name):
        if col_name == 'max_SOCME_S1_Ti_value':
            return 'SOCME'
        elif col_name == 'max_SOCME_over_deltaE_S1_Ti_value':
            return 'SOCME/ΔE'
        return col_name
    
    formatted_columns = [format_column_name(col) for col in columns]
    title = ' vs '.join(formatted_columns)
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_ylabel('Value', fontsize=12, fontweight='bold')
    ax.set_xlabel('moleculars', fontsize=12, fontweight='bold')
    
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
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"小提琴图已保存至: {output_path}")
    
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="从多个CSV文件中提取SOCME数据，生成小提琴图（曲线箱体图）进行对比",
        epilog="""
使用示例：
  单个文件：
    %(prog)s data/pic3inpout_m.csv
  
  多个文件：
    %(prog)s data/pic3inpout_m.csv data/pic9inpout_m.csv data/pic16inpout_m.csv
  
  指定输出文件：
    %(prog)s data/*.csv --output my_socme_boxplot.png
  
  只绘制SOCME：
    %(prog)s data/*.csv --columns max_SOCME_S1_Ti_value

注意：
  - CSV文件必须包含列：max_SOCME_S1_Ti_value, max_SOCME_over_deltaE_S1_Ti_value
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
        default=["max_SOCME_S1_Ti_value", "max_SOCME_over_deltaE_S1_Ti_value"],
        help="要提取和绘制的列名（默认: max_SOCME_S1_Ti_value max_SOCME_over_deltaE_S1_Ti_value）"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="输出图片路径（可选）。如果不指定，默认保存在脚本所在目录，文件名为socme_boxplot.png"
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
        
        # 统计每个列的数据量
        for col in args.columns:
            count = len(data.get(col, []))
            if count > 0:
                print(f"    {col}: {count} 个有效数据点")
            else:
                print(f"    {col}: 无有效数据")
    
    # 检查是否有有效数据
    has_data = False
    for file_label, data in all_data.items():
        for col in args.columns:
            if col in data and len(data[col]) > 0:
                has_data = True
                break
        if has_data:
            break
    
    if not has_data:
        print("错误: 所有文件中都没有找到有效数据")
        return 1
    
    # 确定输出路径
    if args.output:
        output_path = Path(args.output)
    else:
        # 获取脚本所在目录
        script_dir = Path(__file__).parent.absolute()
        output_path = script_dir / "socme_boxplot.png"
    
    # 按照S的顺序排序标签
    sorted_labels = sort_labels_by_s_order(list(all_data.keys()))
    print(f"\n标签排序（按S顺序）: {sorted_labels}")
    
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
    
    # 绘制小提琴图
    print(f"\n正在生成小提琴图（曲线箱体图）...")
    plot_boxplot(
        all_data,
        output_path,
        args.columns,
        figsize=(args.width, args.height),
        sort_order=sorted_labels,
        ref_points=ref_points_dict,
    )
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

