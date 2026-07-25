#!/usr/bin/env python3
"""
从CSV文件读取SMILES，计算分子量，并生成分子量分布图

用法:
    # 处理单个CSV文件（会在原文件中添加Molecular_Weight列）
    python calculate_molecular_weight.py data/pic3.csv
    python calculate_molecular_weight.py data/pic3.csv --output-dir custom_dir
    
    # 处理多个CSV文件（不会修改原文件，只生成总分布图）
    python calculate_molecular_weight.py data/pic3.csv data/pic4.csv data/pic8.csv
    python calculate_molecular_weight.py data/*.csv --output-dir wmimage
    
说明:
    - 单个文件: 会在CSV文件中添加Molecular_Weight列，并生成该文件的分布图
    - 多个文件: 不会修改原CSV文件，只计算所有文件的分子量并生成一个总分布图
    - 输出目录: 默认使用 wmimage/ 目录存放分布图
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
from collections import Counter

# 配置matplotlib字体（使用英文标签，避免中文字体问题）
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# 使用英文标签
LABEL_MW = 'Molecular Weight (g/mol)'
LABEL_FREQ = 'Frequency'
LABEL_TITLE = 'Molecular Weight Distribution'
LABEL_TOTAL = 'Total'
LABEL_MEAN = 'Mean'
LABEL_MEDIAN = 'Median'
LABEL_MIN = 'Min'
LABEL_MAX = 'Max'
LABEL_STD = 'Std'

# 尝试导入RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: RDKit未安装，将尝试使用Open Babel")
    print("建议安装RDKit以获得更好的性能: pip install rdkit")

# 尝试导入Open Babel作为备选
try:
    import openbabel
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False
    if not RDKIT_AVAILABLE:
        print("错误: RDKit和Open Babel都不可用，请至少安装其中一个")
        print("安装RDKit: pip install rdkit")
        print("安装Open Babel: pip install openbabel-wheel 或 conda install -c conda-forge openbabel")
        sys.exit(1)

# 原子量字典（用于备用计算）
ATOMIC_WEIGHTS = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 27.982, 'Si': 28.086, 'P': 30.974,
    'S': 32.066, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Fe': 55.845, 'Br': 79.904, 'I': 126.904
}

def calculate_molecular_weight_rdkit(smiles):
    """
    使用RDKit从SMILES计算分子量
    
    Args:
        smiles: SMILES字符串
        
    Returns:
        分子量（g/mol），如果失败返回None
    """
    if not RDKIT_AVAILABLE:
        return None
    
    try:
        # 从SMILES创建分子对象
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # 使用RDKit的分子量描述符
        mw = Descriptors.MolWt(mol)
        return mw
    except Exception as e:
        return None

def calculate_molecular_weight_openbabel(smiles):
    """
    使用Open Babel从SMILES计算分子量
    
    Args:
        smiles: SMILES字符串
        
    Returns:
        分子量（g/mol），如果失败返回None
    """
    if not OPENBABEL_AVAILABLE:
        return None
    
    try:
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("smi")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, smiles)
        
        if mol.NumAtoms() == 0:
            return None
        
        # 计算分子量
        mw = 0.0
        for atom in openbabel.OBMolAtomIter(mol):
            atomic_num = atom.GetAtomicNum()
            atomic_weight = openbabel.GetAtomicMass(atomic_num)
            mw += atomic_weight
        
        return mw
    except Exception as e:
        return None

def calculate_molecular_weight(smiles):
    """
    计算分子量（优先使用RDKit，失败则尝试Open Babel）
    
    Args:
        smiles: SMILES字符串
        
    Returns:
        分子量（g/mol），如果失败返回None
    """
    # 优先使用RDKit
    if RDKIT_AVAILABLE:
        mw = calculate_molecular_weight_rdkit(smiles)
        if mw is not None:
            return mw
    
    # 如果RDKit失败或不可用，尝试Open Babel
    if OPENBABEL_AVAILABLE:
        mw = calculate_molecular_weight_openbabel(smiles)
        if mw is not None:
            return mw
    
    return None

def calculate_molecular_weight_from_formula(formula_dict):
    """
    从分子式字典计算分子量（备用方法）
    
    Args:
        formula_dict: 元素符号到数量的字典，如 {'C': 20, 'H': 25, 'N': 3, 'O': 2}
        
    Returns:
        分子量（g/mol）
    """
    total_weight = 0.0
    for element, count in formula_dict.items():
        if element in ATOMIC_WEIGHTS:
            total_weight += ATOMIC_WEIGHTS[element] * count
        else:
            print(f"警告: 未知元素 {element}，使用默认原子量")
            total_weight += 1.0 * count  # 默认值
    return total_weight

def get_molecular_formula_rdkit(smiles):
    """
    使用RDKit从SMILES获取分子式
    
    Args:
        smiles: SMILES字符串
        
    Returns:
        元素计数字典，如 {'C': 20, 'H': 25, 'N': 3, 'O': 2}
    """
    if not RDKIT_AVAILABLE:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # 统计原子
        atom_counts = Counter()
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] += 1
        
        return dict(atom_counts)
    except Exception as e:
        print(f"警告: 获取分子式失败 {smiles[:50]}: {e}")
        return None

def format_molecular_formula(formula_dict):
    """
    格式化分子式为字符串（C、H在前，其他按字母顺序）
    
    Args:
        formula_dict: 元素计数字典
        
    Returns:
        格式化的分子式字符串，如 "C20H25N3O2"
    """
    if not formula_dict:
        return ""
    
    # 定义元素顺序：C、H在前，其他按字母顺序
    priority_elements = ['C', 'H']
    other_elements = sorted([e for e in formula_dict.keys() if e not in priority_elements])
    ordered_elements = priority_elements + other_elements
    
    formula_parts = []
    for element in ordered_elements:
        if element in formula_dict and formula_dict[element] > 0:
            count = formula_dict[element]
            if count == 1:
                formula_parts.append(element)
            else:
                formula_parts.append(f"{element}{count}")
    
    return "".join(formula_parts)

def process_csv_file(csv_file_path, output_dir=None, save_to_csv=True):
    """
    处理CSV文件：计算分子量并生成分布图
    
    Args:
        csv_file_path: 输入的CSV文件路径
        output_dir: 输出目录（用于存放图片），如果为None则使用默认目录
        save_to_csv: 是否将分子量保存到CSV文件（单个文件时True，多个文件时False）
    
    Returns:
        分子量列表（成功计算的分子量）
    """
    csv_path = Path(csv_file_path)
    if not csv_path.exists():
        print(f"错误: 文件不存在: {csv_file_path}")
        return
    
    print(f"正在读取文件: {csv_path}")
    
    # 读取CSV文件
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"错误: 无法读取CSV文件: {e}")
        return
    
    # 查找SMILES列
    smiles_col = None
    for col in df.columns:
        if col.upper() == 'SMILES':
            smiles_col = col
            break
    
    if smiles_col is None:
        print("错误: 未找到SMILES列")
        print(f"可用列: {list(df.columns)}")
        return
    
    print(f"找到SMILES列: {smiles_col}")
    print(f"总行数: {len(df)}")
    
    # 计算分子量
    print("\n正在计算分子量...")
    molecular_weights = []
    molecular_formulas = []
    failed_count = 0
    
    for idx, smiles in enumerate(df[smiles_col], 1):
        if pd.isna(smiles) or str(smiles).strip() == '':
            molecular_weights.append(None)
            molecular_formulas.append('')
            failed_count += 1
            continue
        
        smiles_str = str(smiles).strip()
        
        # 计算分子量（自动选择可用的方法）
        mw = calculate_molecular_weight(smiles_str)
        molecular_weights.append(mw)
        
        # 获取分子式（可选）
        formula_dict = get_molecular_formula_rdkit(smiles_str)
        if formula_dict:
            formula_str = format_molecular_formula(formula_dict)
        else:
            formula_str = ''
        molecular_formulas.append(formula_str)
        
        if mw is None:
            failed_count += 1
        
        if idx % 100 == 0:
            print(f"  已处理: {idx}/{len(df)}")
    
    # 如果指定保存，则添加分子量列到DataFrame并保存
    if save_to_csv:
        df['Molecular_Weight'] = molecular_weights
        print(f"\n正在保存更新后的CSV文件...")
        df.to_csv(csv_path, index=False)
        print(f"已保存: {csv_path}")
    
    # 统计信息
    valid_weights = [mw for mw in molecular_weights if mw is not None]
    if valid_weights:
        print(f"\n分子量统计 ({csv_path.name}):")
        print(f"  成功计算: {len(valid_weights)}")
        print(f"  失败: {failed_count}")
        print(f"  最小值: {min(valid_weights):.2f} g/mol")
        print(f"  最大值: {max(valid_weights):.2f} g/mol")
        print(f"  平均值: {np.mean(valid_weights):.2f} g/mol")
        print(f"  中位数: {np.median(valid_weights):.2f} g/mol")
    else:
        print(f"\n警告: {csv_path.name} 没有成功计算任何分子量")
        return []
    
    # 如果save_to_csv为False（多个文件模式），不生成单个文件的分布图
    if not save_to_csv:
        return valid_weights
    
    # 生成分子量分布图
    print("\n正在生成分子量分布图...")
    
    # 确定输出目录
    if output_dir is None:
        output_dir = Path("wmimage")
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 生成图片文件名（基于CSV文件名）
    csv_stem = csv_path.stem  # 不含扩展名的文件名
    image_path = output_dir / f"{csv_stem}_molecular_weight_distribution.png"
    
    # 创建分布图
    plt.figure(figsize=(10, 6))
    
    # 直方图
    plt.hist(valid_weights, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    plt.xlabel(LABEL_MW, fontsize=12)
    plt.ylabel(LABEL_FREQ, fontsize=12)
    
    plt.title(f'{LABEL_TITLE} - {csv_stem}\n'
              f'{LABEL_TOTAL}: {len(valid_weights)}, '
              f'{LABEL_MEAN}: {np.mean(valid_weights):.1f} g/mol, '
              f'{LABEL_MEDIAN}: {np.median(valid_weights):.1f} g/mol', 
              fontsize=14)
    stats_text = (f'{LABEL_MIN}: {min(valid_weights):.1f} g/mol\n'
                  f'{LABEL_MAX}: {max(valid_weights):.1f} g/mol\n'
                  f'{LABEL_MEAN}: {np.mean(valid_weights):.1f} g/mol\n'
                  f'{LABEL_MEDIAN}: {np.median(valid_weights):.1f} g/mol\n'
                  f'{LABEL_STD}: {np.std(valid_weights):.1f} g/mol')
    
    plt.grid(True, alpha=0.3)
    plt.text(0.98, 0.98, stats_text, 
             transform=plt.gca().transAxes,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             fontsize=10)
    
    plt.tight_layout()
    plt.savefig(image_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"已保存分布图: {image_path}")
    print("\n完成！")
    
    return valid_weights

def process_multiple_csv_files(csv_file_paths, output_dir=None):
    """
    处理多个CSV文件：计算所有文件的分子量并生成总分布图
    
    Args:
        csv_file_paths: CSV文件路径列表
        output_dir: 输出目录（用于存放图片），如果为None则使用默认目录
    """
    print(f"检测到 {len(csv_file_paths)} 个CSV文件")
    print("=" * 60)
    
    all_weights = []
    file_stats = []
    
    # 处理每个CSV文件
    for idx, csv_file in enumerate(csv_file_paths, 1):
        print(f"\n[{idx}/{len(csv_file_paths)}] 处理文件: {csv_file}")
        print("-" * 60)
        
        # 不保存到CSV，只收集分子量
        weights = process_csv_file(csv_file, output_dir=None, save_to_csv=False)
        all_weights.extend(weights)
        
        if weights:
            file_stats.append({
                'file': Path(csv_file).name,
                'count': len(weights),
                'min': min(weights),
                'max': max(weights),
                'mean': np.mean(weights),
                'median': np.median(weights)
            })
    
    # 统计所有文件的总体信息
    if not all_weights:
        print("\n错误: 所有文件都没有成功计算任何分子量")
        return
    
    print("\n" + "=" * 60)
    print("总体统计信息:")
    print("=" * 60)
    print(f"  总文件数: {len(csv_file_paths)}")
    print(f"  总分子数: {len(all_weights)}")
    print(f"  最小值: {min(all_weights):.2f} g/mol")
    print(f"  最大值: {max(all_weights):.2f} g/mol")
    print(f"  平均值: {np.mean(all_weights):.2f} g/mol")
    print(f"  中位数: {np.median(all_weights):.2f} g/mol")
    print(f"  标准差: {np.std(all_weights):.2f} g/mol")
    
    # 生成总分布图
    print("\n正在生成总分子量分布图...")
    
    # 确定输出目录
    if output_dir is None:
        output_dir = Path("wmimage")
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 生成图片文件名
    image_path = output_dir / "combined_molecular_weight_distribution.png"
    
    # 创建分布图
    plt.figure(figsize=(12, 7))
    
    # 直方图
    plt.hist(all_weights, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    plt.xlabel(LABEL_MW, fontsize=12)
    plt.ylabel(LABEL_FREQ, fontsize=12)
    
    plt.title(f'{LABEL_TITLE} - Combined ({len(csv_file_paths)} files)\n'
              f'{LABEL_TOTAL}: {len(all_weights)}, '
              f'{LABEL_MEAN}: {np.mean(all_weights):.1f} g/mol, '
              f'{LABEL_MEDIAN}: {np.median(all_weights):.1f} g/mol', 
              fontsize=14)
    
    stats_text = (f'{LABEL_MIN}: {min(all_weights):.1f} g/mol\n'
                  f'{LABEL_MAX}: {max(all_weights):.1f} g/mol\n'
                  f'{LABEL_MEAN}: {np.mean(all_weights):.1f} g/mol\n'
                  f'{LABEL_MEDIAN}: {np.median(all_weights):.1f} g/mol\n'
                  f'{LABEL_STD}: {np.std(all_weights):.1f} g/mol\n'
                  f'Files: {len(csv_file_paths)}')
    
    plt.grid(True, alpha=0.3)
    plt.text(0.98, 0.98, stats_text, 
             transform=plt.gca().transAxes,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             fontsize=10)
    
    plt.tight_layout()
    plt.savefig(image_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"已保存总分布图: {image_path}")
    
    # 打印各文件统计信息
    if file_stats:
        print("\n各文件统计信息:")
        print("-" * 60)
        for stat in file_stats:
            print(f"  {stat['file']}: {stat['count']} molecules, "
                  f"mean={stat['mean']:.1f} g/mol")
    
    print("\n完成！")

def main():
    parser = argparse.ArgumentParser(
        description='从CSV文件读取SMILES，计算分子量，并生成分布图',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 处理单个CSV文件（会在原文件中添加Molecular_Weight列）
  %(prog)s data/pic3.csv
  
  # 处理多个CSV文件（不会修改原文件，只生成总分布图）
  %(prog)s data/pic3.csv data/pic4.csv data/pic8.csv
  
  # 使用通配符处理多个文件
  %(prog)s data/*.csv
  
  # 指定输出目录
  %(prog)s data/pic3.csv --output-dir custom_dir
        """
    )
    parser.add_argument('csv_files', 
                       nargs='+',
                       help='输入的CSV文件路径（可以指定多个文件）')
    parser.add_argument('--output-dir', '-o',
                       default=None,
                       help='输出图片目录（默认: wmimage）')
    
    args = parser.parse_args()
    
    # 检查是单个文件还是多个文件
    if len(args.csv_files) == 1:
        # 单个文件：添加分子量列并生成分布图
        process_csv_file(args.csv_files[0], args.output_dir, save_to_csv=True)
    else:
        # 多个文件：只计算分子量并生成总分布图
        process_multiple_csv_files(args.csv_files, args.output_dir)

if __name__ == '__main__':
    main()

