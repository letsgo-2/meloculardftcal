#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Parse ORCA TD-DFT/SOC logs and create per-molecule and summary CSV files.

Supported modes:

1. Single log:
   ``python3 03_result_parsing/parse_orca_td_socv2.py <log> --smiles <SMILES>``
2. Batch directory:
   ``python3 03_result_parsing/parse_orca_td_socv2.py work/<batch>inp
   --csv data/raw/<batch>.csv --output work/<batch>inpout``
3. Rebuild a summary from existing per-molecule ``*td.csv`` files with
   ``--summary-only``.

Batch logs must follow ``*_<id>td.log`` or ``*_<id>_td.log``. The integer ID
maps to source CSV data row ``id - 1``. Batch mode writes
``<output_dir_name>_tdv2.csv``, matching the merge-stage filename contract.

The parser extracts singlet/triplet energies, oscillator strengths,
wavelengths, S1-Ti gaps, SOCME magnitudes, and SOCME/abs(deltaE). It expects
ORCA sections headed by ``ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE
MOMENTS`` and ``CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS``.
"""

import argparse
import csv
import os
import re
import math
import io
import sys
import glob

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


def get_charge_from_inp(inp_path):
    """
    从 ORCA 输入文件读取体系总电荷。
    支持两种格式：
    - * xyz <charge> <multiplicity>
    - * xyzfile <charge> <multiplicity> <filename>
    """
    if not os.path.exists(inp_path):
        return None

    try:
        with io.open(inp_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        for line in lines:
            s = line.strip()
            # 支持 * xyz 和 * xyzfile 两种格式
            if s.lower().startswith("* xyz") or s.startswith("*xyz"):
                # 去掉开头的 *，然后分割
                # 例如: "* xyz 4 1" -> ["xyz", "4", "1"]
                # 例如: "* xyzfile 4 1 file.xyz" -> ["xyzfile", "4", "1", "file.xyz"]
                parts = s.lstrip("*").strip().split()
                if len(parts) >= 3:  # 至少要有 xyz/xyzfile, charge, multiplicity
                    try:
                        return int(parts[1])  # parts[0] 是 "xyz" 或 "xyzfile", parts[1] 是电荷
                    except (ValueError, IndexError):
                        continue
        return None
    except Exception as e:
        sys.stderr.write("Warning: Failed to parse charge from {}: {}\n".format(inp_path, e))
        return None

def get_smiles_from_csv(log_filename, csv_path):
    """
    从 CSV 文件中根据 log 文件名获取对应的 SMILES。
    
    log 文件名格式: <batch>_XXXXtd.log，其中 XXXX 是数字编号（如 0001）
    CSV 文件行号从1开始（不包括标题行），对应文件名中的数字。
    
    例如: pic16_0001td.log -> CSV 第2行（第1行是标题）
          pic3_0012td.log  -> CSV 第13行（第1行是标题）
    """
    # 从文件名提取数字 ID
    # 匹配任意批次的 <batch>_XXXXtd.log 或 <batch>_XXXX_td.log，
    # 提取最后一个下划线后的数字编号。
    match = re.search(r'.*_(\d+)_?td\.log$', log_filename, re.IGNORECASE)
    if not match:
        return None

    if not csv_path:
        return None
    
    row_num = int(match.group(1))  # 转换为整数，例如 0001 -> 1
    
    try:
        if not os.path.exists(csv_path):
            sys.stderr.write("Warning: CSV file not found: {}\n".format(csv_path))
            return None
            
        with io.open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            
            # row_num 从1开始，对应 rows[0]
            if 1 <= row_num <= len(rows):
                smiles = rows[row_num - 1].get('SMILES', '')
                return smiles if smiles else None
            else:
                sys.stderr.write("Warning: Row number {} out of range (CSV has {} rows)\n".format(row_num, len(rows)))
                return None
    except (IOError, ValueError, KeyError) as e:
        sys.stderr.write("Warning: Failed to read SMILES from CSV: {}\n".format(e))
        return None

def parse_electric_spectrum(text):
    """
    解析 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' 部分。
    返回:
        s_states: { index: (E_eV, E_cm, lambda_nm, fosc) }
        t_states: { index: (E_eV, E_cm, lambda_nm, fosc) }
    """
    header = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    start = text.find(header)
    if start < 0:
        raise ValueError("Electric dipole absorption spectrum section not found")

    next_header = "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS"
    end = text.find(next_header, start)
    if end < 0:
        end = len(text)

    block = text[start:end].splitlines()

    s_states = {}
    t_states = {}

    for line in block:
        line = line.rstrip()
        if "0-1A" not in line or "->" not in line:
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        if parts[0] != "0-1A" or parts[1] != "->":
            continue
        final_state = parts[2]
        try:
            e_ev = float(parts[3])
            e_cm = float(parts[4])
            wl_nm = float(parts[5])
            fosc = float(parts[6])
        except ValueError:
            continue

        m = re.match(r"(\d+)-([13])A", final_state)
        if not m:
            continue

        index = int(m.group(1))
        mult = int(m.group(2))
        if mult == 1:
            s_states[index] = (e_ev, e_cm, wl_nm, fosc)
        elif mult == 3:
            t_states[index] = (e_ev, e_cm, wl_nm, fosc)

    if not s_states or not t_states:
        raise ValueError("Failed to parse singlet or triplet states from electric spectrum")

    return s_states, t_states


def parse_socme(text):
    """
    解析 'CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS' 部分。
    返回:
        soc: { T_index: { S_index: |SOCME| (cm^-1) } }
    """
    header = "CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS"
    start = text.find(header)
    if start < 0:
        raise ValueError("SOCME section not found")

    tail = text[start:].splitlines()
    soc = {}
    data_started = False

    pattern = re.compile(
        r"\s*(\d+)\s+(\d+)\s+"
        r"\(\s*([+-]?\d+\.\d+)\s*,\s*([+-]?\d+\.\d+)\s*\)\s+"
        r"\(\s*([+-]?\d+\.\d+)\s*,\s*([+-]?\d+\.\d+)\s*\)\s+"
        r"\(\s*([+-]?\d+\.\d+)\s*,\s*([+-]?\d+\.\d+)\s*\)"
    )

    for line in tail:
        m = pattern.match(line)
        if m:
            data_started = True
            t_idx = int(m.group(1))
            s_idx = int(m.group(2))
            z_re, z_im = float(m.group(3)), float(m.group(4))
            x_re, x_im = float(m.group(5)), float(m.group(6))
            y_re, y_im = float(m.group(7)), float(m.group(8))

            # 每个方向的复数模长
            mz = math.hypot(z_re, z_im)
            mx = math.hypot(x_re, x_im)
            my = math.hypot(y_re, y_im)
            # 三个方向合成总模长
            mag = math.sqrt(mz**2 + mx**2 + my**2)

            if t_idx not in soc:
                soc[t_idx] = {}
            soc[t_idx][s_idx] = mag
        else:
            if data_started:
                break

    if not soc:
        raise ValueError("No SOCME data parsed")

    return soc


def build_csv_rows(basename, smiles, s_states, t_states, soc, charge=None):
    rows = []

    # 第一行：文件名（不带扩展名）
    rows.append([ensure_unicode(basename)])
    # 第二行：SMILES
    rows.append([ensure_unicode(smiles if smiles is not None else "")])
    # 第三行：体系总电荷（可选）
    rows.append([ensure_unicode(charge if charge is not None else "")])

    # ===== S 态块 =====
    rows.append([])  # 空行分隔
    rows.append([ensure_unicode("S_state"), ensure_unicode("energy_eV"), ensure_unicode("fosc_D2"), ensure_unicode("wavelength_nm")])
    # S0：基态，能量设 0，fosc 设 0
    rows.append([ensure_unicode("S0"), ensure_unicode("0.000000"), ensure_unicode("0.000000000"), ensure_unicode("0.0")])
    for idx in sorted(s_states.keys()):
        e_ev, e_cm, wl_nm, fosc = s_states[idx]
        rows.append([
            ensure_unicode("S{}".format(idx)),
            ensure_unicode("{:.6f}".format(e_ev)),
            ensure_unicode("{:.9f}".format(fosc)),
            ensure_unicode("{:.1f}".format(wl_nm)),
        ])

    # ===== T 态块 =====
    rows.append([])
    rows.append([ensure_unicode("T_state"), ensure_unicode("energy_eV"), ensure_unicode("fosc_D2"), ensure_unicode("wavelength_nm")])
    for idx in sorted(t_states.keys()):
        e_ev, e_cm, wl_nm, fosc = t_states[idx]
        rows.append([
            ensure_unicode("T{}".format(idx)),
            ensure_unicode("{:.6f}".format(e_ev)),
            ensure_unicode("{:.9f}".format(fosc)),
            ensure_unicode("{:.1f}".format(wl_nm)),
        ])

    # ===== S1–Tn 能量差块（单位 eV）=====
    rows.append([])
    rows.append([ensure_unicode("S1_Tn_deltaE_eV")])
    rows.append([ensure_unicode("S_label"), ensure_unicode("T_label"), ensure_unicode("deltaE_eV")])  # deltaE = E(S1) - E(Tn)
    if 1 in s_states:
        e_s1 = s_states[1][0]
        for t_idx in sorted(t_states.keys()):
            e_t = t_states[t_idx][0]
            delta = e_s1 - e_t
            rows.append([
                ensure_unicode("S1"),
                ensure_unicode("T{}".format(t_idx)),
                ensure_unicode("{:.6f}".format(delta)),
            ])

    # ===== SOCME 模长矩阵 =====
    rows.append([])
    rows.append([ensure_unicode("SOCME_abs_cm^-1")])

    all_s_idx = set()
    all_t_idx = set()
    for t_idx, s_map in soc.items():
        all_t_idx.add(t_idx)
        all_s_idx.update(s_map.keys())

    if not all_s_idx or not all_t_idx:
        return rows

    s_indices = sorted(all_s_idx)
    t_indices = sorted(all_t_idx)

    header = [ensure_unicode("T\\S")] + [ensure_unicode("S{}".format(s)) for s in s_indices]
    rows.append(header)

    for t_idx in t_indices:
        row = [ensure_unicode("T{}".format(t_idx))]
        s_map = soc.get(t_idx, {})
        for s_idx in s_indices:
            val = s_map.get(s_idx, None)
            if val is None:
                row.append(ensure_unicode(""))
            else:
                row.append(ensure_unicode("{:.6f}".format(val)))
        rows.append(row)

    # ===== SOCME / |ΔE_ST| （ΔE 用 eV 的绝对值）=====
    rows.append([])
    rows.append([ensure_unicode("SOCME_over_deltaE_(cm^-1_per_eV)")])
    rows.append(header)

    # 各个 S 的能量（S0 设 0）
    s_energy_eV = {0: 0.0}
    for s_idx, (e_ev, e_cm, wl_nm, fosc) in s_states.items():
        s_energy_eV[s_idx] = e_ev

    for t_idx in t_indices:
        e_t = t_states.get(t_idx, (None, None, None, None))[0]
        row = [ensure_unicode("T{}".format(t_idx))]
        for s_idx in s_indices:
            soc_val = soc.get(t_idx, {}).get(s_idx, None)
            e_s = s_energy_eV.get(s_idx)
            if soc_val is None or e_t is None or e_s is None:
                row.append(ensure_unicode(""))
                continue
            delta_e = e_s - e_t
            if delta_e == 0.0:
                row.append(ensure_unicode(""))
            else:
                ratio = soc_val / abs(delta_e)
                row.append(ensure_unicode("{:.6f}".format(ratio)))
        rows.append(row)

    return rows


def parse_generated_csv(csv_path):
    """
    解析已生成的 CSV 文件，提取 SMILES、S1-T1 和 S1-Ti 数据。
    
    Returns:
        (filename, s1_energy, t1_energy, smiles, s1_t1, min_positive_delta, min_positive_ti, charge, s1_fosc, s1_wavelength, 
         max_socme_s1_ti_value, max_socme_s1_ti_index, max_socme_over_delta_s1_ti_value, max_socme_over_delta_s1_ti_index) 或 None 如果失败
    """
    try:
        with io.open(csv_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            return None
        
        # 第1行是文件名（去掉 td 后缀）
        basename_with_td = lines[0].strip()
        # 去掉 td 后缀（可能是 td 或 _td）
        filename = basename_with_td
        if filename.endswith("_td"):
            filename = filename[:-3]
        elif filename.endswith("td"):
            filename = filename[:-2]
        
        # 第2行是 SMILES
        smiles = lines[1].strip()
        # 第3行（可选）是体系总电荷
        charge = None
        if len(lines) >= 3:
            line_e = lines[2].strip()
            if line_e:
                try:
                    charge = int(line_e)
                except ValueError:
                    pass
        
        # 查找 S1 的能量、fosc、波长（在 S_state 部分）
        s1_energy = None
        s1_fosc = None
        s1_wavelength = None
        in_s_state_section = False
        for line in lines:
            line = line.strip()
            if line == "S_state,energy_eV,fosc_D2,wavelength_nm":
                in_s_state_section = True
                continue
            if in_s_state_section:
                if line.startswith("S1,"):
                    parts = line.split(",")
                    if len(parts) >= 4:
                        try:
                            s1_energy = float(parts[1].strip())
                            s1_fosc = float(parts[2].strip())
                            s1_wavelength = float(parts[3].strip())
                        except ValueError:
                            pass
                    break
                if line == "" or line.startswith("T_state"):
                    break
        
        # 查找 T1 的能量（在 T_state 部分）
        t1_energy = None
        in_t_state_section = False
        for line in lines:
            line = line.strip()
            if line == "T_state,energy_eV,fosc_D2,wavelength_nm":
                in_t_state_section = True
                continue
            if in_t_state_section:
                if line.startswith("T1,"):
                    parts = line.split(",")
                    if len(parts) >= 2:
                        try:
                            t1_energy = float(parts[1].strip())
                        except ValueError:
                            pass
                    break
                if line == "" or line.startswith("S1_Tn_deltaE_eV"):
                    break
        
        # 查找 S1_Tn_deltaE_eV 部分
        s1_t1 = None
        min_positive_delta = None
        min_positive_ti = None
        
        in_delta_section = False
        header_passed = False
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            if line == "S1_Tn_deltaE_eV":
                in_delta_section = True
                header_passed = False
                continue
            
            if in_delta_section:
                # 跳过标题行
                if line.startswith("S_label"):
                    header_passed = True
                    continue
                
                # 如果还没过标题行，继续
                if not header_passed:
                    continue
                
                # 处理数据行
                if line.startswith("S1,"):
                    parts = line.split(",")
                    if len(parts) >= 3:
                        t_label = parts[1].strip()
                        delta_str = parts[2].strip()
                        try:
                            delta = float(delta_str)
                            
                            # 提取 T 态的索引
                            if t_label.startswith("T"):
                                t_idx = int(t_label[1:])
                                
                                # S1-T1
                                if t_idx == 1:
                                    s1_t1 = delta
                                
                                # 找正数的最小值
                                if delta > 0:
                                    if min_positive_delta is None or delta < min_positive_delta:
                                        min_positive_delta = delta
                                        min_positive_ti = t_idx
                        except (ValueError, IndexError):
                            continue
                
                # 如果遇到空行且已经处理过数据，停止
                if line == "" and header_passed and s1_t1 is not None:
                    break
                
                # 如果遇到下一个部分（SOCME），停止
                if line.startswith("SOCME"):
                    break
        
        if s1_t1 is None:
            return None
        
        # 解析 SOCME_abs_cm^-1 部分，找到 S1 列的最大值和对应的 Ti
        max_socme_s1_ti_value = None
        max_socme_s1_ti_index = None
        max_socme_over_delta_s1_ti_value = None
        max_socme_over_delta_s1_ti_index = None
        
        in_socme_abs_section = False
        in_socme_over_delta_section = False
        socme_abs_header_passed = False
        socme_over_delta_header_passed = False
        s1_col_index = None
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # 首先检查是否是新的部分标题（优先级最高）
            # 查找 SOCME_abs_cm^-1 部分
            if line == "SOCME_abs_cm^-1":
                in_socme_abs_section = True
                in_socme_over_delta_section = False
                socme_abs_header_passed = False
                socme_over_delta_header_passed = False
                s1_col_index = None
                continue
            
            # 查找 SOCME_over_deltaE_(cm^-1_per_eV) 部分
            if line == "SOCME_over_deltaE_(cm^-1_per_eV)":
                in_socme_abs_section = False
                in_socme_over_delta_section = True
                socme_abs_header_passed = False
                socme_over_delta_header_passed = False
                s1_col_index = None
                continue
            
            # 处理 SOCME_abs_cm^-1 部分
            if in_socme_abs_section:
                # 解析标题行，找到 S1 列的位置
                if line.startswith("T\\S") or line.startswith("T/S"):
                    parts = line.split(",")
                    for idx, col in enumerate(parts):
                        col = col.strip()
                        if col == "S1":
                            s1_col_index = idx
                            break
                    socme_abs_header_passed = True
                    continue
                
                # 处理数据行
                if socme_abs_header_passed and s1_col_index is not None:
                    if line.startswith("T"):
                        parts = line.split(",")
                        if len(parts) > s1_col_index:
                            try:
                                t_label = parts[0].strip()
                                if t_label.startswith("T"):
                                    t_idx = int(t_label[1:])

                                    # 【核心修改点】：增加能级范围过滤
                                    # 只有当 T_i 能量低于 S1 时（t_idx <= min_positive_ti），才参与最大值统计
                                    if min_positive_ti is not None and t_idx <= min_positive_ti:
                                        socme_value = float(parts[s1_col_index].strip())
                                    
                                    
                                        if max_socme_s1_ti_value is None or socme_value > max_socme_s1_ti_value:
                                            max_socme_s1_ti_value = socme_value
                                            max_socme_s1_ti_index = t_idx
                            except (ValueError, IndexError):
                                continue
            
            # 处理 SOCME_over_deltaE_(cm^-1_per_eV) 部分
            if in_socme_over_delta_section:
                # 解析标题行，找到 S1 列的位置
                if line.startswith("T\\S") or line.startswith("T/S"):
                    parts = line.split(",")
                    for idx, col in enumerate(parts):
                        col = col.strip()
                        if col == "S1":
                            s1_col_index = idx
                            break
                    socme_over_delta_header_passed = True
                    continue
                
                # 处理数据行
                if socme_over_delta_header_passed and s1_col_index is not None:
                    if line.startswith("T"):
                        parts = line.split(",")
                        if len(parts) > s1_col_index:
                            try:
                                t_label = parts[0].strip()
                                if t_label.startswith("T"):
                                    t_idx = int(t_label[1:])
                                    # 【核心修改点】：同样的过滤逻辑
                                    if min_positive_ti is not None and t_idx <= min_positive_ti:
                                        socme_over_delta_value = float(parts[s1_col_index].strip())
                                    
                                        if max_socme_over_delta_s1_ti_value is None or socme_over_delta_value > max_socme_over_delta_s1_ti_value:
                                            max_socme_over_delta_s1_ti_value = socme_over_delta_value
                                            max_socme_over_delta_s1_ti_index = t_idx
                            except (ValueError, IndexError):
                                continue
        
        return (filename, s1_energy, t1_energy, smiles, s1_t1, min_positive_delta, min_positive_ti, charge, s1_fosc, s1_wavelength,
                max_socme_s1_ti_value, max_socme_s1_ti_index, max_socme_over_delta_s1_ti_value, max_socme_over_delta_s1_ti_index)
    
    except Exception as e:
        sys.stderr.write("Warning: Failed to parse CSV {}: {}\n".format(csv_path, e))
        return None


def generate_summary_csv(input_dir, output_dir=None, summary_filename=None):
    """
    从输入目录中的所有 CSV 文件生成汇总 CSV。
    
    Args:
        input_dir: 包含 CSV 文件的目录
        output_dir: 输出目录（如果为 None，则使用 input_dir）
        summary_filename: 汇总文件名（如果为 None，则使用"文件夹名_summary.csv"）
    
    Returns:
        True 如果成功，False 如果失败
    """
    try:
        # 查找所有 CSV 文件
        pattern = os.path.join(input_dir, "*td.csv")
        csv_files = sorted(glob.glob(pattern))
        
        if not csv_files:
            sys.stderr.write("No CSV files found in directory: {}\n".format(input_dir))
            return False
        
        summary_rows = []
        summary_rows.append([ensure_unicode("filename"), 
                            ensure_unicode("S1_energy_eV"),
                            ensure_unicode("T1_energy_eV"),
                            ensure_unicode("S1_fosc_D2"),
                            ensure_unicode("S1_wavelength_nm"),
                            ensure_unicode("charge"),
                            ensure_unicode("SMILES"), 
                            ensure_unicode("S1-T1"), 
                            ensure_unicode("min_positive_S1-Ti"), 
                            ensure_unicode("Ti_for_min_positive"),
                            ensure_unicode("max_SOCME_S1_Ti_value"),
                            ensure_unicode("max_SOCME_S1_Ti_index"),
                            ensure_unicode("max_SOCME_over_deltaE_S1_Ti_value"),
                            ensure_unicode("max_SOCME_over_deltaE_S1_Ti_index")])
        
        for csv_file in csv_files:
            result = parse_generated_csv(csv_file)
            if result:
                (filename, s1_energy, t1_energy, smiles, s1_t1, 
                 min_positive_delta, min_positive_ti, charge,
                 s1_fosc, s1_wavelength,
                 max_socme_s1_ti_value, max_socme_s1_ti_index,
                 max_socme_over_delta_s1_ti_value, max_socme_over_delta_s1_ti_index) = result
                
                # 格式化数据
                s1_energy_str = "{:.6f}".format(s1_energy) if s1_energy is not None else ""
                t1_energy_str = "{:.6f}".format(t1_energy) if t1_energy is not None else ""
                s1_fosc_str = "{:.9f}".format(s1_fosc) if s1_fosc is not None else ""
                s1_wl_str = "{:.1f}".format(s1_wavelength) if s1_wavelength is not None else ""
                charge_str = "{}".format(charge) if charge is not None else ""
                s1_t1_str = "{:.6f}".format(s1_t1) if s1_t1 is not None else ""
                min_delta_str = "{:.6f}".format(min_positive_delta) if min_positive_delta is not None else ""
                min_ti_str = "T{}".format(min_positive_ti) if min_positive_ti is not None else ""
                max_socme_value_str = "{:.6f}".format(max_socme_s1_ti_value) if max_socme_s1_ti_value is not None else ""
                max_socme_ti_str = "T{}".format(max_socme_s1_ti_index) if max_socme_s1_ti_index is not None else ""
                max_socme_over_delta_value_str = "{:.6f}".format(max_socme_over_delta_s1_ti_value) if max_socme_over_delta_s1_ti_value is not None else ""
                max_socme_over_delta_ti_str = "T{}".format(max_socme_over_delta_s1_ti_index) if max_socme_over_delta_s1_ti_index is not None else ""
                
                summary_rows.append([
                    ensure_unicode(filename),
                    ensure_unicode(s1_energy_str),
                    ensure_unicode(t1_energy_str),
                    ensure_unicode(s1_fosc_str),
                    ensure_unicode(s1_wl_str),
                    ensure_unicode(charge_str),
                    ensure_unicode(smiles),
                    ensure_unicode(s1_t1_str),
                    ensure_unicode(min_delta_str),
                    ensure_unicode(min_ti_str),
                    ensure_unicode(max_socme_value_str),
                    ensure_unicode(max_socme_ti_str),
                    ensure_unicode(max_socme_over_delta_value_str),
                    ensure_unicode(max_socme_over_delta_ti_str)
                ])
        
        # 确定输出目录
        if output_dir is None:
            output_dir = input_dir
        else:
            # 确保输出目录存在
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        
        # 确定summary文件名
        if summary_filename is None:
            # 使用输入目录的文件夹名作为summary文件名前缀
            dir_basename = os.path.basename(os.path.abspath(input_dir))
            summary_filename = "{}_summary.csv".format(dir_basename)
        
        # 写入汇总文件
        summary_path = os.path.join(output_dir, summary_filename)
        
        if sys.version_info[0] >= 3:
            newline_arg = ""
        else:
            newline_arg = None
        
        with io.open(summary_path, "w", newline=newline_arg, encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerows(summary_rows)
        
        print("Summary CSV generated: {} (from {} CSV files)".format(summary_path, len(csv_files)))
        return True
    
    except Exception as e:
        sys.stderr.write("Error generating summary CSV: {}\n".format(e))
        return False


def process_single_log_file(log_path, output_path, smiles_csv_path, smiles_override=None):
    """
    处理单个 log 文件，提取信息并保存到 CSV。
    
    Args:
        log_path: log 文件路径
        output_path: 输出 CSV 文件路径
        smiles_csv_path: SMILES CSV 文件路径
        smiles_override: 手动指定的 SMILES（可选）
    
    Returns:
        True 如果成功，False 如果失败
    """
    try:
        # Python 2/3 兼容的读文件方式
        with io.open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()

        s_states, t_states = parse_electric_spectrum(text)
        soc = parse_socme(text)

        basename = os.path.splitext(os.path.basename(log_path))[0]
        # 读取体系总电荷：从对应的 td.inp 文件中读取（如 xhx01td.log -> xhx01td.inp）
        charge = None
        log_dir = os.path.dirname(log_path)
        # 首先尝试同目录下的 td.inp 文件
        inp_candidate = os.path.join(log_dir, basename + ".inp")
        if not os.path.exists(inp_candidate):
            # 如果找不到，尝试父目录（如输出在 pic16inpout，输入在 pic16inp）
            inp_candidate = os.path.join(os.path.dirname(os.path.abspath(log_path)), "..", basename + ".inp")
            inp_candidate = os.path.abspath(inp_candidate)
        
        if os.path.exists(inp_candidate):
            charge = get_charge_from_inp(inp_candidate)
        
        # 获取 SMILES：优先使用手动指定，否则从 CSV 文件读取
        smiles = smiles_override
        if smiles is None:
            log_filename = os.path.basename(log_path)
            csv_path = smiles_csv_path
            # 如果 CSV 路径是相对路径，尝试多个位置
            if not os.path.isabs(csv_path):
                # 首先尝试相对于输入文件所在目录
                input_dir = os.path.dirname(os.path.abspath(log_path))
                test_path = os.path.join(input_dir, csv_path)
                if os.path.exists(test_path):
                    csv_path = test_path
                else:
                    # 尝试相对于脚本所在目录
                    script_dir = os.path.dirname(os.path.abspath(__file__))
                    test_path = os.path.join(script_dir, csv_path)
                    if os.path.exists(test_path):
                        csv_path = test_path
                    else:
                        # 最后尝试相对于当前工作目录
                        csv_path = os.path.abspath(csv_path)
            smiles = get_smiles_from_csv(log_filename, csv_path)
            if smiles is None:
                smiles = ""  # 如果找不到，使用空字符串
        
        rows = build_csv_rows(basename, smiles, s_states, t_states, soc, charge)

        # 确保输出目录存在
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Python 2/3 兼容的写 CSV
        if sys.version_info[0] >= 3:
            newline_arg = ""
        else:
            newline_arg = None

        with io.open(output_path, "w", newline=newline_arg, encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerows(rows)
        
        return True
    except Exception as e:
        sys.stderr.write("Error processing {}: {}\n".format(log_path, e))
        return False


def main():
    parser = argparse.ArgumentParser(description="Parse ORCA TD-DFT log (with SOC) into CSV.")
    parser.add_argument("input", help="ORCA TD-DFT log file or directory containing *td.log files (or *td.csv files if --summary-only)")
    parser.add_argument("-o", "--output", help="Output CSV path (for single file) or output directory (for directory input)")
    parser.add_argument("--smiles", help="SMILES string to write on the second line (overrides CSV lookup)", default=None)
    parser.add_argument(
        "--csv",
        help="Source CSV used for ID-to-SMILES row lookup; required unless --smiles is supplied",
        default=None,
    )
    parser.add_argument("--summary-only", action="store_true", help="Only generate summary CSV from existing CSV files (skip log processing)")
    args = parser.parse_args()

    in_path = args.input

    # 如果只是生成 summary，直接处理
    if args.summary_only:
        if not os.path.isdir(in_path):
            sys.stderr.write("Error: --summary-only requires a directory containing CSV files: {}\n".format(in_path))
            sys.exit(1)
        
        # --- 优化后的输出处理逻辑 ---
        summary_out_dir = None
        summary_out_name = None

        if args.output:
            # 判断用户输入的是否是具体的 .csv 文件
            if args.output.lower().endswith('.csv'):
                summary_out_dir = os.path.dirname(args.output)
                summary_out_name = os.path.basename(args.output)
                # 如果用户只写了 -o test.csv，dirname 会是空，设为当前目录
                if not summary_out_dir:
                    summary_out_dir = "."
            else:
                # 如果输入的是目录
                summary_out_dir = args.output
                summary_out_name = None
        else:
            # 默认：输出到输入目录，文件名自动生成
            summary_out_dir = in_path
            summary_out_name = None
        
        print("Generating summary CSV...")
        print("Input directory: {}".format(in_path))
        
        # 调用函数，传入具体的文件名（如果有的话）
        success = generate_summary_csv(in_path, summary_out_dir, summary_out_name)
        sys.exit(0 if success else 1)

    # 检查输入是文件还是目录
    if os.path.isfile(in_path):
        # 单文件模式
        out_path = args.output
        if not out_path:
            out_path = os.path.splitext(in_path)[0] + ".csv"
        
        success = process_single_log_file(in_path, out_path, args.csv, args.smiles)
        sys.exit(0 if success else 1)
    
    elif os.path.isdir(in_path):
        # 目录模式：批量处理
        # 查找所有 *td.log 文件
        pattern = os.path.join(in_path, "*td.log")
        log_files = glob.glob(pattern)
        
        if not log_files:
            sys.stderr.write("No *td.log files found in directory: {}\n".format(in_path))
            sys.exit(1)
        
        # 确定输出目录
        if args.output:
            output_dir = args.output
        else:
            # 自动创建：输入文件夹名 + "out"
            input_basename = os.path.basename(os.path.abspath(in_path))
            input_dir = os.path.dirname(os.path.abspath(in_path))
            output_dir = os.path.join(input_dir, input_basename + "out")
        
        # 确保输出目录存在
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        print("Processing {} log files...".format(len(log_files)))
        print("Output directory: {}".format(output_dir))
        
        # 处理每个文件
        success_count = 0
        fail_count = 0
        
        for log_file in sorted(log_files):
            log_basename = os.path.basename(log_file)
            csv_basename = os.path.splitext(log_basename)[0] + ".csv"
            output_path = os.path.join(output_dir, csv_basename)
            
            print("Processing: {} -> {}".format(log_basename, csv_basename))
            
            if process_single_log_file(log_file, output_path, args.csv, args.smiles):
                success_count += 1
            else:
                fail_count += 1
        
        print("\nCompleted: {} succeeded, {} failed".format(success_count, fail_count))
        
        # 生成汇总 CSV
        if success_count > 0:
            print("\nGenerating summary CSV...")
            # 与 merge_dft_csv.py 的输入契约保持一致：
            # work/pic16inpout -> pic16inpout_tdv2.csv
            summary_name = "{}_tdv2.csv".format(
                os.path.basename(os.path.abspath(output_dir))
            )
            generate_summary_csv(output_dir, output_dir, summary_name)
        
        sys.exit(0 if fail_count == 0 else 1)
    
    else:
        sys.stderr.write("Error: Input path does not exist: {}\n".format(in_path))
        sys.exit(1)


if __name__ == "__main__":
    main()
