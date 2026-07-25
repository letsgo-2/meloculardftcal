#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制 min(S1-Ti) 和 SOCME 的散点图与帕累托前沿

功能说明：
    从CSV文件中提取 min(S1-Ti) 和 SOCME 数据，绘制散点图。
    计算并显示帕累托前沿（Pareto Front），其中：
    - min(S1-Ti) 越小越好
    - SOCME 越大越好
    
    将数据分成三段进行可视化（根据数据范围自动划分）。

使用方法：
    python plot/plot_pareto_s1ti_socme.py data/pic16inpout_m.csv
    
    指定输出文件：
    python plot/plot_pareto_s1ti_socme.py data/pic16inpout_m.csv --output my_pareto.png
"""

import argparse
import sys
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.spatial import ConvexHull


def load_pictos_labels(pictos_file: Path) -> Dict[str, str]:
    """
    从pictos.txt文件加载pic名称到标签的映射
    
    Args:
        pictos_file: pictos.txt文件路径
    
    Returns:
        字典，键为pic名称，值为标签（如"S9"）
    """
    labels = {}
    if not pictos_file.exists():
        return labels
    
    try:
        with open(pictos_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                # 支持冒号或空格分隔
                if ':' in line:
                    parts = line.split(':', 1)
                    pic_name = parts[0].strip()
                    label = parts[1].strip()
                    labels[pic_name] = label
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        pic_name = parts[0]
                        label = parts[1]
                        labels[pic_name] = label
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
    filename_lower = filename.lower()
    for pic in pictos_labels.keys():
        if filename_lower.startswith(pic.lower()):
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

def extract_data_from_csv(csv_file: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    从CSV文件中提取数据
    
    Args:
        csv_file: CSV文件路径
    
    Returns:
        (min_s1_ti, socme, ids) - 三个数组
    """
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"错误: 无法读取文件 {csv_file}: {e}")
        sys.exit(1)
    
    # 检查必要的列
    if 'min_positive_S1-Ti' not in df.columns:
        print("错误: CSV文件中没有找到 'min_positive_S1-Ti' 列")
        sys.exit(1)
    
    if 'max_SOCME_over_deltaE_S1_Ti_value' not in df.columns:
        print("错误: CSV文件中没有找到 'max_SOCME_over_deltaE_S1_Ti_value' 列")
        sys.exit(1)
    
    # 提取数据，过滤掉缺失值
    mask = df['min_positive_S1-Ti'].notna() & df['max_SOCME_over_deltaE_S1_Ti_value'].notna()
    min_s1_ti = df.loc[mask, 'min_positive_S1-Ti'].values.astype(float)
    socme = df.loc[mask, 'max_SOCME_over_deltaE_S1_Ti_value'].values.astype(float)
    
    # 获取ID（如果有）
    if 'ID' in df.columns:
        ids = df.loc[mask, 'ID'].values
    else:
        ids = np.arange(len(min_s1_ti))
    
    # 排除 SOCME > 10000 的数据点
    socme_filter = socme <= 10000
    min_s1_ti = min_s1_ti[socme_filter]
    socme = socme[socme_filter]
    ids = ids[socme_filter]
    
    excluded_count = np.sum(~socme_filter)
    if excluded_count > 0:
        print(f"  排除了 {excluded_count} 个 SOCME > 10000 的数据点")
    
    print(f"从 {csv_file.name} 中提取了 {len(min_s1_ti)} 个有效数据点")
    
    return min_s1_ti, socme, ids


def load_publish_points(
    publish_csv: Path,
    pictos_labels: Dict[str, str],
) -> Dict[str, Tuple[float, float]]:
    """
    从 publishinpdftv22.csv 中读取每个 pic 的基准点坐标，用于在单个骨架图上高亮一个点。
    
    返回:
        { "pic3": (min_positive_S1-Ti, max_SOCME_over_deltaE_S1_Ti_value), ... }
    """
    publish_points: Dict[str, Tuple[float, float]] = {}
    if not publish_csv.exists():
        print(f"提示: 未找到基准点文件 {publish_csv}，将不在单图中额外标点。")
        return publish_points
    try:
        df = pd.read_csv(publish_csv)
    except Exception as e:
        print(f"警告: 无法读取基准点文件 {publish_csv}: {e}")
        return publish_points

    required_cols = [
        "index",
        "min_positive_S1-Ti",
        "max_SOCME_over_deltaE_S1_Ti_value",
    ]
    for col in required_cols:
        if col not in df.columns:
            print(f"警告: 基准点文件中缺少列 '{col}'，将不在单图中额外标点。")
            return {}

    for _, row in df.iterrows():
        idx = str(row["index"]).strip()  # 一般是 pic3, pic4 ...
        if not idx:
            continue
        try:
            x_val = float(row["min_positive_S1-Ti"])
            y_val = float(row["max_SOCME_over_deltaE_S1_Ti_value"])
        except (TypeError, ValueError):
            continue
        if np.isnan(x_val) or np.isnan(y_val):
            continue

        # 1) 以 pic 名称作为键存一份
        publish_points[idx] = (x_val, y_val)

        # 2) 根据 pictos.txt 映射到 S 标签，再存一份，便于按 S* 匹配
        s_label_raw = pictos_labels.get(idx)
        s_label_formatted = None
        if s_label_raw:
            # 仿照 get_label_from_filename 的格式化规则，统一成 S1 / S9V1 这种
            if s_label_raw.lower().startswith("s"):
                s_label_formatted = f"S{s_label_raw[1:].upper()}"
            else:
                s_label_formatted = s_label_raw
            publish_points[s_label_formatted] = (x_val, y_val)

        if s_label_formatted:
            print(
                f"  基准点: index={idx} -> label={s_label_formatted}, "
                f"min_positive_S1-Ti={x_val:.4f}, max_SOCME_over_deltaE_S1_Ti_value={y_val:.4f}"
            )
        else:
            print(
                f"  基准点: index={idx}, min_positive_S1-Ti={x_val:.4f}, "
                f"max_SOCME_over_deltaE_S1_Ti_value={y_val:.4f} (无 pictos 标签映射)"
            )

    if publish_points:
        print(f"已从 {publish_csv.name} 读取到 {len(publish_points)} 个基准点，将在对应骨架单图上高亮显示。")
    else:
        print(f"注意: 从 {publish_csv.name} 中没有成功读取到任何有效基准点。")
    return publish_points


def find_pareto_front(min_s1_ti: np.ndarray, socme: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    找到帕累托前沿
    
    帕累托前沿定义：
    - min(S1-Ti) 越小越好
    - SOCME 越大越好
    
    Args:
        min_s1_ti: min(S1-Ti) 值数组
        socme: SOCME 值数组
    
    Returns:
        (pareto_min_s1_ti, pareto_socme) - 帕累托前沿上的点
    """
    # 合并数据
    points = np.column_stack([min_s1_ti, -socme])  # 将SOCME取负，转换为最小化问题
    
    # 找到非支配点（Pareto最优解）
    pareto_indices = []
    n = len(points)
    
    for i in range(n):
        is_pareto = True
        for j in range(n):
            if i != j:
                # 如果点j在所有目标上都优于点i，则i不是帕累托最优
                if (points[j, 0] <= points[i, 0] and points[j, 1] <= points[i, 1]) and \
                   (points[j, 0] < points[i, 0] or points[j, 1] < points[i, 1]):
                    is_pareto = False
                    break
        
        if is_pareto:
            pareto_indices.append(i)
    
    pareto_indices = np.array(pareto_indices)
    
    if len(pareto_indices) == 0:
        return np.array([]), np.array([])
    
    # 按 min_s1_ti 排序
    sorted_indices = np.argsort(min_s1_ti[pareto_indices])
    pareto_indices_sorted = pareto_indices[sorted_indices]
    
    return min_s1_ti[pareto_indices_sorted], socme[pareto_indices_sorted]


def classify_data_into_segments(min_s1_ti: np.ndarray, socme: np.ndarray) -> np.ndarray:
    """
    将数据分成三段（基于综合评分）
    
    综合评分：考虑 min(S1-Ti) 越小越好，SOCME 越大越好
    使用归一化的综合得分进行分类
    
    Args:
        min_s1_ti: min(S1-Ti) 值数组
        socme: SOCME 值数组
    
    Returns:
        segment_labels - 每个点的分段标签 (0=低, 1=中, 2=高)
    """
    # 归一化数据（0-1范围）
    min_s1_ti_norm = (min_s1_ti - min_s1_ti.min()) / (min_s1_ti.max() - min_s1_ti.min() + 1e-10)
    socme_norm = (socme - socme.min()) / (socme.max() - socme.min() + 1e-10)
    
    # 综合评分：min(S1-Ti)越小越好（取反），SOCME越大越好
    # 得分 = -min_s1_ti_norm + socme_norm
    scores = -min_s1_ti_norm + socme_norm
    
    # 分成三段
    q33 = np.percentile(scores, 33.33)
    q67 = np.percentile(scores, 66.67)
    
    segments = np.zeros(len(scores), dtype=int)
    segments[scores > q67] = 2  # 高
    segments[(scores > q33) & (scores <= q67)] = 1  # 中
    segments[scores <= q33] = 0  # 低
    
    return segments

def plot_individual_scaffold_with_marginals(
    min_s1_ti: np.ndarray, 
    socme: np.ndarray, 
    display_label: str, 
    color: str, 
    marker: str, 
    output_dir: Path,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    highlight_point: Optional[Tuple[float, float]] = None,
):
    """
    为每个骨架生成带有边际分布图的散点图
    """
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4, hspace=0.2, wspace=0.2)
    
    # 主散点图
    ax_main = fig.add_subplot(gs[1:, :3])
    # X轴分布 (顶部)
    ax_histx = fig.add_subplot(gs[0, :3], sharex=ax_main)
    # Y轴分布 (右侧)
    ax_histy = fig.add_subplot(gs[1:, 3], sharey=ax_main)
    
    # 1. 绘制主散点图并设置坐标范围（确定范围后再画直方图以统一间隔）
    ax_main.scatter(min_s1_ti, socme, c=color, marker=marker, s=50, edgecolors='black', linewidths=0.5, alpha=0.7)
    ax_main.set_xlabel('min(S1-Ti) (eV)', fontweight='bold')
    ax_main.set_ylabel('SOCME/ΔE', fontweight='bold')
    
    if xlim: ax_main.set_xlim(xlim)
    else: ax_main.set_xlim(0, 0.6)
    
    if ylim and (ylim[0] is not None or ylim[1] is not None):
        ax_main.set_ylim(ylim)
    
    # 获取当前轴的显示范围
    curr_xlim = ax_main.get_xlim()
    curr_ylim = ax_main.get_ylim()
    
    # 高亮一个额外点（来自 publishinpdftv22.csv）
    if highlight_point is not None:
        hx, hy = highlight_point
        ax_main.scatter(
            [hx],
            [hy],
            c="red",
            marker="*",
            s=140,
            edgecolors="black",
            linewidths=0.8,
            zorder=5,
            label="publish point",
        )
        # 计算相对于基准点，满足条件的点的比例
        # 条件：min(S1-Ti) <= 基准点的 min(S1-Ti) 且 SOCME/ΔE >= 基准点的 SOCME/ΔE
        better_or_equal_mask = (min_s1_ti <= hx) & (socme >= hy)
        better_or_equal_count = np.sum(better_or_equal_mask)
        total_count = len(min_s1_ti)
        if total_count > 0:
            better_or_equal_ratio = better_or_equal_count / total_count
            print(
                f"    相对于基准点 ({hx:.4f}, {hy:.4f}): "
                f"{better_or_equal_count}/{total_count} 个点满足条件 "
                f"(min(S1-Ti) <= {hx:.4f} 且 SOCME/ΔE >= {hy:.4f}), "
                f"比例 = {better_or_equal_ratio:.2%}"
            )
        # 避免重复图例，这里不专门加 legend，由颜色和形状区分
    
    # 2. 绘制边际分布 (通过 range 参数强制直方图在整个显示范围内均匀划分 Bin)
    # X分布
    ax_histx.hist(min_s1_ti, bins=30, range=curr_xlim, color=color, alpha=0.4, density=True, edgecolor='white')
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histx.set_yticks([])
    ax_histx.spines['top'].set_visible(False)
    ax_histx.spines['right'].set_visible(False)
    ax_histx.spines['left'].set_visible(False)
    
    # Y分布
    ax_histy.hist(socme, bins=30, range=curr_ylim, color=color, alpha=0.4, orientation='horizontal', density=True, edgecolor='white')
    ax_histy.tick_params(axis="y", labelleft=False)
    ax_histy.set_xticks([])
    ax_histy.spines['top'].set_visible(False)
    ax_histy.spines['right'].set_visible(False)
    ax_histy.spines['bottom'].set_visible(False)

    fig.suptitle(f'Scaffold: {display_label}', fontsize=14, fontweight='bold', y=0.95)
    
    save_path = output_dir / f"individual_dist_{display_label}.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  单个骨架分布图已保存至: {save_path}")
def plot_total_distribution_with_marginals(
    data_dict: dict,
    display_labels: dict,
    dataset_colors: list,
    dataset_markers: list,
    sorted_labels: list,
    output_dir: Path,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None
):
    """
    新增功能：将所有数据汇总在一张图中，边际直方图显示全体数据的总分布。
    """
    fig = plt.figure(figsize=(10, 10))
    gs = GridSpec(4, 4, hspace=0.2, wspace=0.2)
    
    ax_main = fig.add_subplot(gs[1:, :3])
    ax_histx = fig.add_subplot(gs[0, :3], sharex=ax_main)
    ax_histy = fig.add_subplot(gs[1:, 3], sharey=ax_main)
    
    # 设置范围
    if xlim: ax_main.set_xlim(xlim)
    else: ax_main.set_xlim(0, 0.6)
    if ylim: ax_main.set_ylim(ylim)
    
    curr_xlim = ax_main.get_xlim()
    curr_ylim = ax_main.get_ylim()

    all_x = []
    all_y = []

    # 循环绘制散点并收集所有数据
    for idx, label in enumerate(sorted_labels):
        x, y, _ = data_dict[label]
        all_x.extend(x)
        all_y.extend(y)
        
        color = dataset_colors[idx % len(dataset_colors)]
        marker = dataset_markers[idx % len(dataset_markers)]
        
        ax_main.scatter(x, y, c=color, marker=marker, s=35, alpha=0.5, 
                        edgecolors='black', linewidths=0.2, label=display_labels[label])

    # 绘制全体数据的单一总分布直方图（灰色代表总和）
    ax_histx.hist(all_x, bins=30, range=curr_xlim, color='#7F8C8D', alpha=0.5, density=True, edgecolor='white')
    ax_histy.hist(all_y, bins=30, range=curr_ylim, color='#7F8C8D', alpha=0.5, orientation='horizontal', density=True, edgecolor='white')

    ax_main.set_xlabel('min(S1-Ti) (eV)', fontweight='bold')
    ax_main.set_ylabel('SOCME/ΔE', fontweight='bold')
    ax_main.legend(loc='best', fontsize=8, ncol=2, framealpha=0.8)
    
    ax_histx.axis('off')
    ax_histy.axis('off')

    fig.suptitle('All Scaffolds', fontsize=15, fontweight='bold', y=0.96)
    
    save_path = output_dir / "total_combined_distribution.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"全体数据汇总分布图已保存至: {save_path}")

def plot_pareto_scatter(
    data_dict: dict,  # {label: (min_s1_ti, socme, ids)}
    output_path: Path,
    figsize: Tuple[int, int] = (10, 8),
    label_mapping: Optional[Dict[str, str]] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_pareto_for: Optional[List[str]] = None,  # 指定要显示帕累托曲线的标签列表
    publish_points: Optional[Dict[str, Tuple[float, float]]] = None,  # 每个 pic 的基准点
):
    """
    绘制散点图和帕累托前沿（支持多个数据集）
    
    Args:
        data_dict: 字典，键为标签，值为(min_s1_ti, socme, ids)元组
        output_path: 输出图片路径
        figsize: 图片大小
    """
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    # 定义不同数据集的颜色（使用更明显的颜色）
    dataset_colors = [
        '#E74C3C',  # 红色
        '#3498DB',  # 蓝色
        '#2ECC71',  # 绿色
        '#F39C12',  # 橙色
        '#9B59B6',  # 紫色
        '#1ABC9C',  # 青色
        '#E67E22',  # 深橙色
        '#34495E',  # 深灰色
        '#16A085',  # 深青色
        '#D35400',  # 深红色
    ]
    dataset_markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h', 'P']
    
    # 定义三段颜色（每个数据集内部使用）
    segment_colors = ['#FFB3BA', '#BAE1FF', '#BAFFC9']  # 低、中、高（淡红、淡蓝、淡绿）
    segment_labels = ['Low', 'Medium', 'High']
    
    all_pareto_points = []  # 存储所有数据集的帕累托点
    
    # 如果有标签映射，使用映射后的标签
    if label_mapping:
        display_labels = {key: label_mapping.get(key, key) for key in data_dict.keys()}
    else:
        display_labels = {key: key for key in data_dict.keys()}
    
    # 按照标签排序（如果有S标签，按S数字排序）
    def get_s_number(label):
        if label.upper().startswith('S'):
            try:
                return int(label.upper()[1:])
            except:
                return 999
        return 999
    
    sorted_labels = sorted(data_dict.keys(), key=lambda x: get_s_number(display_labels.get(x, x)))
    
    output_dir = output_path.parent

    # --- 新增图表调用：汇总全体数据分布图 ---
    plot_total_distribution_with_marginals(
        data_dict, display_labels, dataset_colors, dataset_markers, 
        sorted_labels, output_dir, xlim, ylim
    )
    
    for idx, label in enumerate(sorted_labels):
        min_s1_ti, socme, ids = data_dict[label]
        display_label = display_labels[label]
        if len(min_s1_ti) == 0:
            continue
        
        # 找到帕累托前沿
        pareto_min, pareto_socme = find_pareto_front(min_s1_ti, socme)
        
        # 将数据分成三段
        segments = classify_data_into_segments(min_s1_ti, socme)
        
        # 选择颜色和标记
        dataset_color = dataset_colors[idx % len(dataset_colors)]
        dataset_marker = dataset_markers[idx % len(dataset_markers)]
        
        # --- 原有：为每个骨架生成带边际分布的单图 ---
        # 单个骨架图：如果有 publish 基准点，则额外高亮一个点
        highlight_pt = None
        if publish_points is not None:
            key_used = None
            # 先尝试用 pic 名（label，例如 pic3）
            if label in publish_points:
                key_used = label
            # 再尝试用 pictos 映射后的 S 标签（display_label，例如 S9V1）
            elif display_label in publish_points:
                key_used = display_label

            if key_used is not None:
                highlight_pt = publish_points[key_used]
                print(
                    f"  单图高亮点: pic={label}, display_label={display_label}, "
                    f"匹配键={key_used}, 点=({highlight_pt[0]:.4f}, {highlight_pt[1]:.4f})"
                )
            else:
                print(
                    f"  单图无基准点: pic={label}, display_label={display_label} "
                    f"在 publish 基准点中未找到对应 index/标签"
                )
        plot_individual_scaffold_with_marginals(
            min_s1_ti,
            socme,
            display_label,
            dataset_color,
            dataset_marker,
            output_dir,
            xlim,
            ylim,
            highlight_point=highlight_pt,
        )
        
        # 绘制散点图（所有数据点使用相同样式）
        ax.scatter(
            min_s1_ti,
            socme,
            c=dataset_color,
            alpha=0.5,
            s=40,
            edgecolors='black',
            linewidths=0.3,
            marker=dataset_marker,
            label=f'{display_label}',
            zorder=1
        )
        
        # 绘制帕累托前沿（根据开关决定是否显示）
        should_show_pareto = True
        if show_pareto_for is not None:
            # 如果指定了要显示帕累托曲线的标签列表，检查当前标签是否在列表中
            # 检查原始标签或显示标签
            should_show_pareto = (label in show_pareto_for or display_label in show_pareto_for)
        
        if len(pareto_min) > 0 and should_show_pareto:
            # 排序以便绘制连续线
            sort_idx = np.argsort(pareto_min)
            pareto_min_sorted = pareto_min[sort_idx]
            pareto_socme_sorted = pareto_socme[sort_idx]
            
            # 绘制帕累托前沿曲线（使用实线，线宽更细）
            ax.plot(
                pareto_min_sorted,
                pareto_socme_sorted,
                color=dataset_color,
                linewidth=1.5,
                alpha=0.9,
                linestyle='-',
                label=f'{display_label} Pareto Front',
                zorder=3
            )
            
            all_pareto_points.append((display_label, len(pareto_min)))
            print(f"{display_label} ({label}): 找到 {len(pareto_min)} 个帕累托最优解，已绘制曲线")
        elif len(pareto_min) > 0:
            print(f"{display_label} ({label}): 找到 {len(pareto_min)} 个帕累托最优解，但根据设置不显示曲线")
    
    # 设置标签和标题
    ax.set_xlabel('min(S1-Ti) (eV)', fontsize=10, fontweight='bold')
    ax.set_ylabel('SOCME (cm⁻¹)', fontsize=10, fontweight='bold')
    ax.set_title('Pareto Front:min(S1-Ti) ↓, SOCME ↑', 
                 fontsize=11, fontweight='bold', pad=15)
    
    # 设置x轴范围
    if xlim is not None:
        left = xlim[0] if xlim[0] is not None else None
        right = xlim[1] if xlim[1] is not None else None
        if left is not None or right is not None:
            ax.set_xlim(left=left, right=right)
    else:
        # 默认限制最大值为0.6
        ax.set_xlim(left=0, right=0.6)
    
    # 设置y轴范围（如果指定）
    if ylim is not None:
        bottom = ylim[0] if ylim[0] is not None else None
        top = ylim[1] if ylim[1] is not None else None
        if bottom is not None or top is not None:
            ax.set_ylim(bottom=bottom, top=top)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # 添加图例（只显示数据集和帕累托前沿，不显示帕累托点）
    # 创建自定义图例，只包含数据集和帕累托前沿
    from matplotlib.lines import Line2D
    legend_elements = []
    
    # 按照相同的顺序生成图例
    for idx, label in enumerate(sorted_labels):
        min_s1_ti, socme, ids = data_dict[label]
        display_label = display_labels[label]
        
        if len(min_s1_ti) == 0:
            continue
        
        dataset_color = dataset_colors[idx % len(dataset_colors)]
        dataset_marker = dataset_markers[idx % len(dataset_markers)]
        
        # 数据集散点图
        legend_elements.append(
            plt.Line2D([0], [0], marker=dataset_marker, color='w', 
                      markerfacecolor=dataset_color, markersize=8, 
                      markeredgecolor='black', markeredgewidth=0.5,
                      label=f'{display_label}', linestyle='None')
        )
        
        # 帕累托前沿线（如果有且允许显示）
        should_show_pareto = True
        if show_pareto_for is not None:
            should_show_pareto = (label in show_pareto_for or display_label in show_pareto_for)
        
        pareto_min, pareto_socme = find_pareto_front(min_s1_ti, socme)
        if len(pareto_min) > 0 and should_show_pareto:
            legend_elements.append(
                Line2D([0], [0], color=dataset_color, linewidth=1.5,
                      label=f'{display_label} Pareto Front')
            )
    
    ax.legend(handles=legend_elements, loc='best', framealpha=0.9, fontsize=10, ncol=2)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"总散点图和帕累托前沿已保存至: {output_path}")
    
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="绘制 min(S1-Ti) 和 SOCME 的散点图与帕累托前沿",
        epilog="""
使用示例：
  单个文件：
    %(prog)s data/pic16inpout_m.csv
  
  多个文件：
    %(prog)s data/pic3inpout_m.csv data/pic16inpout_m.csv
  
  指定输出文件：
    %(prog)s data/pic16inpout_m.csv --output my_pareto.png
  
  自定义图片大小：
    %(prog)s data/pic16inpout_m.csv --width 18 --height 10
  
  设置横轴范围：
    %(prog)s data/pic16inpout_m.csv --xmin 0 --xmax 0.5
  
  设置纵轴范围：
    %(prog)s data/pic16inpout_m.csv --ymin 0 --ymax 1000
  
  只显示特定文件的帕累托曲线（支持 pic* 或 S* 格式）：
    %(prog)s data/pic3inpout_m.csv data/pic16inpout_m.csv --show-pareto pic3 pic16
    %(prog)s data/pic3inpout_m.csv data/pic16inpout_m.csv --show-pareto S9 S16

注意：
  - CSV文件必须包含列：min_positive_S1-Ti, max_SOCME_over_deltaE_S1_Ti_value
  - 输出文件默认保存在脚本所在目录
  - 支持多个文件输入，会在同一张图上显示
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
        "--output", "-o",
        type=str,
        default=None,
        help="输出图片路径（可选）。如果不指定，默认保存在脚本所在目录，文件名为pareto_s1ti_socmegap.png"
    )
    parser.add_argument(
        "--width",
        type=int,
        default=10,
        help="图片宽度，单位：英寸（默认: 10）"
    )
    parser.add_argument(
        "--height",
        type=int,
        default=8,
        help="图片高度，单位：英寸（默认: 8）"
    )
    parser.add_argument(
        "--ymin",
        type=float,
        default=None,
        help="纵轴（SOCME）最小值（默认: 自动）"
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="纵轴（SOCME）最大值（默认: 自动）"
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=None,
        help="横轴（min(S1-Ti)）最小值（默认: 0）"
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="横轴（min(S1-Ti)）最大值（默认: 0.6）"
    )
    parser.add_argument(
        "--show-pareto",
        nargs='+',
        type=str,
        default=None,
        help="指定要显示帕累托曲线的文件标签（例如: pic3 pic16 或 S9 S16）。支持 pic* 和 S* 格式。如果不指定，默认显示所有文件的帕累托曲线。散点图始终显示。"
    )
    
    args = parser.parse_args()
    
    # 加载pictos标签映射
    script_dir = Path(__file__).parent.absolute()
    pictos_file = script_dir.parent / "data" / "pictos.txt"
    pictos_labels = load_pictos_labels(pictos_file)
    
    if pictos_labels:
        print(f"已加载 {len(pictos_labels)} 个pic标签映射")

    # 读取 publish 基准点（用于在单个骨架图上高亮一个点）
    publish_csv = script_dir.parent / "data" / "publishinpdftv22.csv"
    publish_points = load_publish_points(publish_csv, pictos_labels)
    
    # 处理输入文件
    data_dict = {}
    label_mapping = {}  # pic名称到显示标签的映射
    
    for input_file in args.input_files:
        input_path = Path(input_file)
        if not input_path.exists():
            print(f"警告: 文件不存在，跳过: {input_path}")
            continue
        
        # 从文件名提取pic名称（如 pic16inpout_m.csv -> pic16）
        pic_name = input_path.stem.replace('inpout_m', '').replace('inpout', '')
        if not pic_name:
            pic_name = input_path.stem
        
        # 使用pictos标签，如果没有则使用pic名称
        display_label = get_label_from_filename(pic_name, pictos_labels)
        label_mapping[pic_name] = display_label
        
        # 提取数据
        print(f"正在从 {input_path.name} 中提取数据...")
        min_s1_ti, socme, ids = extract_data_from_csv(input_path)
        
        if len(min_s1_ti) == 0:
            print(f"警告: {input_path.name} 中没有找到有效数据，跳过")
            continue
        
        data_dict[pic_name] = (min_s1_ti, socme, ids)
    
    if len(data_dict) == 0:
        print("错误: 没有找到有效数据")
        return 1
    
    # 确定输出路径
    if args.output:
        output_path = Path(args.output)
    else:
        # 获取脚本所在目录
        script_dir = Path(__file__).parent.absolute()
        output_path = script_dir / "pareto_s1ti_socmegap.png"
    
    # 准备x轴和y轴范围参数
    xlim = None
    if args.xmin is not None or args.xmax is not None:
        xlim = (args.xmin if args.xmin is not None else None, 
                args.xmax if args.xmax is not None else None)
    
    ylim = None
    if args.ymin is not None or args.ymax is not None:
        ylim = (args.ymin if args.ymin is not None else None, 
                args.ymax if args.ymax is not None else None)
    
    # 处理帕累托曲线显示开关
    show_pareto_labels = []
    if args.show_pareto:
        # 将用户输入的标签转换为对应的 pic_name 和 display_label
        show_pareto_labels = []
        for user_label in args.show_pareto:
            # 标准化用户输入的标签（统一大小写）
            user_label_normalized = user_label.strip()
            # 如果输入的是 S* 格式，统一为大写
            if user_label_normalized.lower().startswith('s') and len(user_label_normalized) > 1:
                user_label_normalized = f"S{user_label_normalized[1:].upper()}"
            
            # 尝试匹配 pic_name 或 display_label（大小写不敏感）
            found = False
            for pic_name, display_label in label_mapping.items():
                # 匹配 pic_name（如 pic3）
                if user_label_normalized.lower() == pic_name.lower():
                    show_pareto_labels.append(pic_name)  # 添加原始标签
                    show_pareto_labels.append(display_label)  # 添加显示标签
                    found = True
                    break
                # 匹配 display_label（如 S9）
                elif user_label_normalized.upper() == display_label.upper():
                    show_pareto_labels.append(pic_name)  # 添加原始标签
                    show_pareto_labels.append(display_label)  # 添加显示标签
                    found = True
                    break
            
            if not found:
                # 如果没找到，直接使用用户输入的标签（可能是文件名）
                show_pareto_labels.append(user_label_normalized)
                print(f"警告: 未找到标签 '{user_label}' 的匹配项，将尝试直接使用")
        
        print(f"将显示以下标签的帕累托曲线: {args.show_pareto}")
        print(f"  匹配到的标签列表: {show_pareto_labels}")
    
    # 绘制散点图和帕累托前沿
    print(f"\n正在生成散点图和帕累托前沿...")
    plot_pareto_scatter(
        data_dict,
        output_path,
        figsize=(args.width, args.height),
        label_mapping=label_mapping,
        xlim=xlim,
        ylim=ylim,
        show_pareto_for=show_pareto_labels,
        publish_points=publish_points,
    )
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())