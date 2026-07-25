#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
合并 TD-DFT summary CSV 和基态优化 CSV 文件

这个脚本用于将两个CSV文件合并成一个文件：
1. {input_dir}out/{input_dir}out_summary.csv (TD-DFT数据)
2. {input_dir}out/{input_dir}out_opt.csv (基态优化数据)

合并后的文件保存为：{input_dir}out/{input_dir}dft.csv

使用方法：
    python merge_dft_csv.py pic3inp
    python merge_dft_csv.py pic3inp -o custom_output_dir

命令行参数：
    input               输入目录名（如 pic3inp）
    -o, --output        输出目录（可选，默认：{input_dir}out）
"""

import argparse
import csv
import os
import sys
import io

def ensure_unicode(s):
    """确保字符串是 unicode（Python 2）或 str（Python 3）"""
    if sys.version_info[0] >= 3:
        return str(s) if s is not None else ""
    else:
        if isinstance(s, unicode):
            return s
        elif isinstance(s, str):
            return s.decode('utf-8')
        else:
            return unicode(s) if s is not None else u""


def read_csv_to_dict(csv_path):
    """
    读取CSV文件，返回以filename为key的字典
    
    Args:
        csv_path: CSV文件路径
    
    Returns:
        dict: {filename: {column: value, ...}, ...} 或 None 如果失败
    """
    try:
        if not os.path.exists(csv_path):
            sys.stderr.write("Warning: CSV file not found: {}\n".format(csv_path))
            return None
        
        data = {}
        with io.open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                filename = row.get('filename', '').strip()
                if filename:
                    data[filename] = row
        
        return data
    except Exception as e:
        sys.stderr.write("Error reading CSV file {}: {}\n".format(csv_path, e))
        return None


def merge_csv_files(summary_data, opt_data, output_path):
    """
    合并两个CSV数据字典
    
    Args:
        summary_data: TD-DFT summary数据字典
        opt_data: 基态优化数据字典
        output_path: 输出文件路径
    
    Returns:
        True 如果成功，False 如果失败
    """
    try:
        # 获取所有唯一的filename
        all_filenames = set()
        if summary_data:
            all_filenames.update(summary_data.keys())
        if opt_data:
            all_filenames.update(opt_data.keys())
        
        if not all_filenames:
            sys.stderr.write("Error: No data to merge\n")
            return False
        
        # 确定所有列名
        summary_columns = []
        opt_columns = []
        
        if summary_data:
            # 从第一个条目获取列名（排除filename）
            first_summary = next(iter(summary_data.values()))
            summary_columns = [col for col in first_summary.keys() if col != 'filename']
        
        if opt_data:
            # 从第一个条目获取列名（排除filename）
            first_opt = next(iter(opt_data.values()))
            opt_columns = [col for col in first_opt.keys() if col != 'filename']
        
        # 合并列名：filename + summary列 + opt列
        all_columns = ['filename'] + summary_columns + opt_columns
        
        # 写入合并后的CSV
        if sys.version_info[0] >= 3:
            newline_arg = ""
        else:
            newline_arg = None
        
        with io.open(output_path, "w", newline=newline_arg, encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=all_columns)
            writer.writeheader()
            
            # 按filename排序写入
            for filename in sorted(all_filenames):
                merged_row = {'filename': ensure_unicode(filename)}
                
                # 添加summary数据
                if summary_data and filename in summary_data:
                    for col in summary_columns:
                        merged_row[col] = ensure_unicode(summary_data[filename].get(col, ''))
                else:
                    # 如果summary中没有这个filename，填充空值
                    for col in summary_columns:
                        merged_row[col] = ensure_unicode("")
                
                # 添加opt数据
                if opt_data and filename in opt_data:
                    for col in opt_columns:
                        merged_row[col] = ensure_unicode(opt_data[filename].get(col, ''))
                else:
                    # 如果opt中没有这个filename，填充空值
                    for col in opt_columns:
                        merged_row[col] = ensure_unicode("")
                
                writer.writerow(merged_row)
        
        return True
    except Exception as e:
        sys.stderr.write("Error merging CSV files: {}\n".format(e))
        return False


def main():
    parser = argparse.ArgumentParser(description="Merge TD-DFT summary and ground-state optimization CSV files")
    parser.add_argument("input", help="Input directory name (e.g., pic3inp)")
    parser.add_argument("-o", "--output", help="Output directory (optional, default: {input_dir}out)")
    args = parser.parse_args()
    
    input_dir_name = args.input
    
    # 确定输入目录的绝对路径
    input_abs = os.path.abspath(input_dir_name)
    input_basename = os.path.basename(input_abs)
    input_parent = os.path.dirname(input_abs)
    
    # 确定输出目录
    if args.output:
        output_dir = args.output
    else:
        output_dir_name = input_basename + "out"
        output_dir = os.path.join(input_parent, output_dir_name)
    
    if not os.path.exists(output_dir):
        sys.stderr.write("Error: Output directory does not exist: {}\n".format(output_dir))
        sys.exit(1)
    
    # 确定两个输入CSV文件路径
    summary_csv_name = input_basename + "out_tdv2.csv"
    opt_csv_name = input_basename + "out_opt.csv"
    
    summary_csv_path = os.path.join(output_dir, summary_csv_name)
    opt_csv_path = os.path.join(output_dir, opt_csv_name)
    
    # 确定输出文件路径
    output_csv_name = input_basename + "dftv20.csv"
    output_csv_path = os.path.join(output_dir, output_csv_name)
    
    print("Input directory: {}".format(input_dir_name))
    print("Output directory: {}".format(output_dir))
    print("Summary CSV: {}".format(summary_csv_path))
    print("Opt CSV: {}".format(opt_csv_path))
    print("Output CSV: {}".format(output_csv_path))
    
    # 读取两个CSV文件
    print("\nReading CSV files...")
    summary_data = read_csv_to_dict(summary_csv_path)
    opt_data = read_csv_to_dict(opt_csv_path)
    
    if summary_data is None and opt_data is None:
        sys.stderr.write("Error: Both CSV files are missing or cannot be read\n")
        sys.exit(1)
    
    if summary_data:
        print("  Found {} entries in summary CSV".format(len(summary_data)))
    else:
        print("  Warning: Summary CSV not found or empty")
    
    if opt_data:
        print("  Found {} entries in opt CSV".format(len(opt_data)))
    else:
        print("  Warning: Opt CSV not found or empty")
    
    # 合并文件
    print("\nMerging CSV files...")
    success = merge_csv_files(summary_data, opt_data, output_csv_path)
    
    if success:
        print("Successfully merged CSV files!")
        print("Output file: {}".format(output_csv_path))
        sys.exit(0)
    else:
        sys.stderr.write("Error: Failed to merge CSV files\n")
        sys.exit(1)


if __name__ == "__main__":
    main()

