#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
装饰片段统计与可视化脚本

功能说明：
    从CSV或xlsx文件中读取scaffold字段，统计装饰片段出现频率，并生成柱状图。
    柱状图显示前N个最常见的装饰片段，每个柱子上方显示对应的分子结构图片。

主要特性：
    1. 自动规范化SMILES字符串，确保相同分子被正确统计
    2. 在柱状图上显示分子结构图（透明背景，无边框）
    3. 支持CSV和xlsx格式输入
    4. 输出文件自动按照输入文件名命名

输入文件格式要求：
    - 输入文件必须包含"Scaffold"列（或通过--scaffold-col指定列名）
    - Scaffold列的格式：骨架|装饰1|装饰2|装饰3|...
      第一个片段是骨架，用"|"分隔的后续片段都是装饰片段

使用示例：
    1. 基本用法（使用默认设置，输出文件自动命名为输入文件名.png）：
       python3 plot_decorations_statistics.py data/pic16inpout_m.csv
    
    2. 显示前5个最常见装饰片段：
       python3 plot_decorations_statistics.py data/pic16inpout_m.csv --top-n 5
    
    3. 指定输出路径：
       python3 plot_decorations_statistics.py data/pic16inpout_m.csv --output plot/my_result.png
    
    4. 自定义图片大小：
       python3 plot_decorations_statistics.py data/pic16inpout_m.csv --width 16 --height 10
    
    5. 指定scaffold列名（如果列名不是"Scaffold"）：
       python3 plot_decorations_statistics.py data/myfile.csv --scaffold-col "MyScaffold"
    
    6. 组合使用多个参数：
       python3 plot_decorations_statistics.py data/pic16inpout_m.csv --top-n 15 --width 18 --height 10 --output plot/top15.png

输出说明：
    - 如果不指定--output参数，输出文件将保存在脚本所在目录（plot目录）
    - 输出文件名 = 输入文件名（不含扩展名）+ ".png"
    - 例如：输入 data/pic16inpout_m.csv → 输出 plot/pic16inpout_m.png

注意事项：
    - 需要安装：pandas, matplotlib, rdkit, pillow
    - 相同分子的不同SMILES表示会被规范化后合并统计
    - 分子图片背景透明，无边框
"""
from __future__ import annotations

from typing import Dict, List, Tuple, Optional
import argparse
import io
from pathlib import Path
from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image


def _strip_numbers_for_drawing(m: Chem.Mol) -> None:
    """防止 RDKit 把 map number / isotope 当成标签画出来"""
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
        a.SetIsotope(0)


def draw_molecule_to_image(smiles: str, size: Tuple[int, int] = (200, 200)) -> Optional[Image.Image]:
    """
    将SMILES字符串绘制为PIL Image对象
    返回装饰片段的图片（保留*标记）
    """
    if not smiles or not smiles.strip():
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None
        
        _strip_numbers_for_drawing(mol)
        
        # 获取dummy原子并标记为"*"
        dummy_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        
        try:
            rdDepictor.SetPreferCoordGen(True)
        except Exception:
            pass
        rdDepictor.Compute2DCoords(mol)
        
        drawer = rdMolDraw2D.MolDraw2DCairo(int(size[0]), int(size[1]))
        opts = drawer.drawOptions()
        
        # 不显示编号
        opts.addAtomIndices = False
        opts.addBondIndices = False
        
        # 关闭可能显示数字的选项
        for attr in ("includeAtomTags", "includeAtomMapNumbers", "addAtomMapNumbers", "explicitMethyl"):
            if hasattr(opts, attr):
                try:
                    setattr(opts, attr, False)
                except Exception:
                    pass
        
        # 标记dummy原子为"*"
        if dummy_atoms:
            for idx in dummy_atoms:
                if 0 <= idx < mol.GetNumAtoms():
                    opts.atomLabels[idx] = "*"
        
        drawer.DrawMolecule(mol, highlightAtoms=dummy_atoms)
        drawer.FinishDrawing()
        
        img_bytes = drawer.GetDrawingText()
        img = Image.open(io.BytesIO(img_bytes))
        
        # 将白色背景转换为透明
        if img.mode == 'RGBA':
            # 如果已经是RGBA模式，直接返回
            return img
        elif img.mode == 'RGB':
            # 转换为RGBA，并将白色像素设为透明
            img = img.convert('RGBA')
            data = img.getdata()
            new_data = []
            for item in data:
                # 将白色（或接近白色）的像素设为透明
                if item[0] > 240 and item[1] > 240 and item[2] > 240:
                    new_data.append((255, 255, 255, 0))  # 完全透明
                else:
                    new_data.append(item)
            img.putdata(new_data)
        return img
        
    except Exception as e:
        print(f"警告: 无法绘制 {smiles}: {e}")
        return None


def extract_decorations_from_scaffold(scaffold: str) -> List[str]:
    """
    从scaffold字段中提取装饰片段
    scaffold格式: 骨架|装饰1|装饰2|...
    返回装饰片段列表（不包括骨架）
    """
    if not scaffold or pd.isna(scaffold):
        return []
    
    parts = str(scaffold).strip().split('|')
    if len(parts) <= 1:
        return []
    
    # 第一个是骨架，后面的都是装饰片段
    decorations = [part.strip() for part in parts[1:] if part.strip()]
    return decorations


def canonicalize_smiles(smiles: str) -> Optional[str]:
    """
    将SMILES字符串规范化为标准形式
    如果无法解析，返回None
    """
    if not smiles or not smiles.strip():
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None
        # 规范化SMILES（相同的分子会有相同的规范化SMILES）
        canon_smiles = Chem.MolToSmiles(mol)
        return canon_smiles
    except Exception:
        return None


def count_decorations(df: pd.DataFrame, scaffold_col: str = "Scaffold") -> Tuple[Counter, Dict[str, str]]:
    """
    统计所有装饰片段的数量（使用规范化的SMILES）
    返回:
        Counter对象：键为规范化SMILES，值为出现次数
        Dict对象：规范化SMILES到原始SMILES的映射（用于显示）
    """
    all_decorations = []
    original_smiles_map = {}  # 规范化SMILES -> 原始SMILES（保存第一个遇到的原始形式）
    
    for scaffold in df[scaffold_col]:
        decorations = extract_decorations_from_scaffold(scaffold)
        for dec in decorations:
            # 规范化SMILES
            canon_smiles = canonicalize_smiles(dec)
            if canon_smiles:
                all_decorations.append(canon_smiles)
                # 保存第一个遇到的原始SMILES作为代表
                if canon_smiles not in original_smiles_map:
                    original_smiles_map[canon_smiles] = dec
    
    return Counter(all_decorations), original_smiles_map


def plot_top_decorations(
    decoration_counts: List[Tuple[str, int]],
    original_smiles_map: Dict[str, str],
    output_path: Path,
    top_n: int = 10,
    figsize: Tuple[int, int] = (14, 8),
    bar_width: float = 0.6,
):
    """
    绘制前N个最常见装饰片段的柱状图，每个柱子上方显示装饰片段图片
    
    Args:
        decoration_counts: (规范化SMILES, count) 元组列表，已按count降序排序
        original_smiles_map: 规范化SMILES到原始SMILES的映射
        output_path: 输出图片路径
        top_n: 显示前N个
        figsize: 图片大小
        bar_width: 柱子宽度
    """
    # 只取前N个
    top_decorations = decoration_counts[:top_n]
    
    if not top_decorations:
        print("警告: 没有装饰片段可显示")
        return
    
    # 准备数据（使用规范化SMILES来绘制，因为它们代表实际相同的分子）
    canon_smiles_list = [item[0] for item in top_decorations]
    counts_list = [item[1] for item in top_decorations]
    # 获取原始SMILES用于显示（如果规范化版本无法解析，则尝试原始版本）
    smiles_list = []
    for canon_smiles in canon_smiles_list:
        # 优先使用规范化SMILES，如果原始SMILES在映射中，也可以尝试使用它
        smiles_to_use = canon_smiles
        if canon_smiles in original_smiles_map:
            # 两个都尝试，优先用能画的
            smiles_list.append((canon_smiles, original_smiles_map[canon_smiles]))
        else:
            smiles_list.append((canon_smiles, canon_smiles))
    
    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)
    
    # 绘制柱状图
    x_pos = np.arange(len(smiles_list))
    bars = ax.bar(x_pos, counts_list, width=bar_width, color='steelblue', edgecolor='black', linewidth=1.2)
    
    # 设置x轴标签（使用简短标签）
    labels = [f"Dec {i+1}" for i in range(len(smiles_list))]
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=0, ha='center')
    
    # 设置y轴和标题
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_xlabel('Decoration Fragments', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {len(top_decorations)} Most Common Decoration Fragments', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # 在柱子上方显示数值
    for i, (bar, count) in enumerate(zip(bars, counts_list)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{count}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # 在柱子上方添加分子图片
    mol_img_size = (100, 100)  # 每个分子图片的大小
    
    for i, (canon_smiles, orig_smiles) in enumerate(smiles_list):
        # 尝试用规范化SMILES画图，如果失败则用原始SMILES
        img = draw_molecule_to_image(canon_smiles, size=mol_img_size)
        if img is None and orig_smiles != canon_smiles:
            img = draw_molecule_to_image(orig_smiles, size=mol_img_size)
        
        if img is not None:
            # 计算图片位置（柱子顶部上方）
            bar = bars[i]
            x_coord = bar.get_x() + bar.get_width() / 2
            y_coord = bar.get_height()
            
            # 创建OffsetImage
            imagebox = OffsetImage(img, zoom=0.6)
            
            # 创建AnnotationBbox，位置在柱子顶部上方（无边框，透明背景）
            ab = AnnotationBbox(
                imagebox,
                (x_coord, y_coord),
                frameon=False,  # 无边框
                box_alignment=(0.5, 0.5),
                xybox=(0, 25),  # 图片位置相对于柱顶向上偏移25像素
                xycoords='data',
                boxcoords='offset points',
            )
            ax.add_artist(ab)
        else:
            print(f"警告: 无法为 {canon_smiles[:50]}... 生成图片")
    
    # 调整布局，为顶部图片留出空间
    plt.subplots_adjust(top=0.85, bottom=0.15)
    
    # 添加网格
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # 保存图片
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"柱状图已保存至: {output_path}")
    
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="统计装饰片段出现频率并生成可视化柱状图",
        epilog="""
使用示例：
  基本用法：
    %(prog)s data/pic16inpout_m.csv
  
  显示前5个最常见装饰片段：
    %(prog)s data/pic16inpout_m.csv --top-n 5
  
  指定输出路径：
    %(prog)s data/pic16inpout_m.csv --output plot/my_result.png
  
  自定义图片大小：
    %(prog)s data/pic16inpout_m.csv --width 16 --height 10
  
  指定scaffold列名：
    %(prog)s data/myfile.csv --scaffold-col "MyScaffold"

输入文件格式：
  - 必须包含Scaffold列（格式：骨架|装饰1|装饰2|...）
  - 第一个片段是骨架，后续用"|"分隔的都是装饰片段
  
输出说明：
  - 默认输出：与输入文件同名，扩展名为.png，保存在脚本所在目录
  - 例如：输入 data/file.csv → 输出 plot/file.png
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "input_file", 
        type=str, 
        help="输入文件路径（支持.csv或.xlsx格式）。文件必须包含Scaffold列，格式为：骨架|装饰1|装饰2|..."
    )
    parser.add_argument(
        "--scaffold-col", 
        default="Scaffold", 
        help="Scaffold列的名称（默认: 'Scaffold'）。如果输入文件的列名不是'Scaffold'，使用此参数指定"
    )
    parser.add_argument(
        "--output", "-o", 
        type=str, 
        default=None, 
        help="输出图片路径（可选）。如果不指定，输出文件名将与输入文件相同，扩展名改为.png，保存在脚本所在目录（plot目录）"
    )
    parser.add_argument(
        "--top-n", 
        type=int, 
        default=10, 
        help="显示前N个最常见装饰片段（默认: 10）。按出现频率从高到低排序"
    )
    parser.add_argument(
        "--width", 
        type=int, 
        default=14, 
        help="输出图片的宽度，单位：英寸（默认: 14）。可根据需要显示的装饰片段数量调整"
    )
    parser.add_argument(
        "--height", 
        type=int, 
        default=8, 
        help="输出图片的高度，单位：英寸（默认: 8）。可根据需要调整"
    )
    
    args = parser.parse_args()
    
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"错误: 文件不存在: {input_path}")
        return 1
    
    # 读取文件
    try:
        if input_path.suffix.lower() == '.xlsx':
            df = pd.read_excel(input_path)
        elif input_path.suffix.lower() == '.csv':
            df = pd.read_csv(input_path)
        else:
            print(f"错误: 不支持的文件格式 {input_path.suffix}，请使用.xlsx或.csv")
            return 1
    except Exception as e:
        print(f"错误: 无法读取文件: {e}")
        return 1
    
    # 检查scaffold列是否存在
    if args.scaffold_col not in df.columns:
        print(f"错误: 找不到列 '{args.scaffold_col}'")
        print(f"可用列: {', '.join(df.columns)}")
        return 1
    
    # 统计装饰片段（使用规范化SMILES）
    print(f"正在统计装饰片段（规范化SMILES）...")
    decoration_counter, original_smiles_map = count_decorations(df, scaffold_col=args.scaffold_col)
    
    if not decoration_counter:
        print("警告: 没有找到任何装饰片段")
        return 0
    
    print(f"共找到 {len(decoration_counter)} 种不同的装饰片段（规范化后）")
    print(f"装饰片段总出现次数: {sum(decoration_counter.values())}")
    
    # 获取最常见的装饰片段
    top_decorations = decoration_counter.most_common(args.top_n)
    
    print(f"\n前 {len(top_decorations)} 个最常见装饰片段（规范化后）:")
    for i, (canon_smiles, count) in enumerate(top_decorations, 1):
        orig_smiles = original_smiles_map.get(canon_smiles, canon_smiles)
        print(f"  {i}. Count={count}: {canon_smiles[:80]}{'...' if len(canon_smiles) > 80 else ''}")
        if orig_smiles != canon_smiles:
            print(f"      原始: {orig_smiles[:80]}{'...' if len(orig_smiles) > 80 else ''}")
    
    # 确定输出路径（默认保存在脚本所在目录，使用输入文件名，扩展名改为.png）
    if args.output:
        output_path = Path(args.output)
    else:
        # 获取脚本所在目录
        script_dir = Path(__file__).parent.absolute()
        # 使用输入文件名（不含扩展名）加上.png作为输出文件名，保存在脚本目录
        output_path = script_dir / f"{input_path.stem}.png"
    
    # 绘制柱状图
    print(f"\n正在生成柱状图...")
    plot_top_decorations(
        top_decorations,
        original_smiles_map,
        output_path,
        top_n=args.top_n,
        figsize=(args.width, args.height),
    )
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

