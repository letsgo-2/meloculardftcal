#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import re
import io
import sys
import glob

def parse_ground_state_properties(text):
    """从 ORCA Log 中解析基态核心物理量，确保只取最后一次出现的值"""
    results = {
        "total_energy": None, "homo": None, "lumo": None, "gap": None,
        "dipole": None, "gibbs": None, "enthalpy": None, "entropy": None,
        "imaginary": "No",
        "nuclear_repulsion": None,
        "electronic_energy": None,
        "one_electron_energy": None,
        "two_electron_energy": None,
        "cpcm_dielectric": None
    }

    # 1. 总电子能量 - 使用 findall 并取 [-1]，确保是最后一次优化后的能量
    energy_matches = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-]?\d+\.\d+)', text)
    if energy_matches:
        results["total_energy"] = float(energy_matches[-1])

    # 2. 详细能量组成 - 使用 rfind 从后往前定位最后一个 TOTAL SCF ENERGY 块
    scf_block_start = text.rfind("TOTAL SCF ENERGY")
    if scf_block_start != -1:
        # 截取该块后的一段文本进行解析
        scf_text = text[scf_block_start:scf_block_start+2500]
        
        nuc_m = re.search(r'Nuclear Repulsion\s+:\s+([-]?\d+\.\d+)', scf_text)
        ele_m = re.search(r'Electronic Energy\s+:\s+([-]?\d+\.\d+)', scf_text)
        one_m = re.search(r'One Electron Energy:\s+([-]?\d+\.\d+)', scf_text)
        two_m = re.search(r'Two Electron Energy:\s+([-]?\d+\.\d+)', scf_text)
        cpc_m = re.search(r'CPCM Dielectric\s+:\s+([-]?\d+\.\d+)', scf_text)
        
        if nuc_m: results["nuclear_repulsion"] = float(nuc_m.group(1))
        if ele_m: results["electronic_energy"] = float(ele_m.group(1))
        if one_m: results["one_electron_energy"] = float(one_m.group(1))
        if two_m: results["two_electron_energy"] = float(two_m.group(1))
        if cpc_m: results["cpcm_dielectric"] = float(cpc_m.group(1))

    # 3. 轨道能量 (HOMO / LUMO) - 使用 rfind 锁定最后一次出现的轨道块
    orb_marker = "ORBITAL ENERGIES"
    orb_start = text.rfind(orb_marker)
    if orb_start != -1:
        orb_lines = text[orb_start:].splitlines()
        row_pat = re.compile(r'^\s*(\d+)\s+([0-2]\.\d+)\s+[-]?\d+\.\d+\s+([-]?\d+\.\d+)')
        orbitals = []
        for line in orb_lines:
            m = row_pat.match(line)
            if m:
                orbitals.append({'occ': float(m.group(2)), 'ev': float(m.group(3))})
        
        for i in range(len(orbitals)):
            if orbitals[i]['occ'] > 0.1:
                results["homo"] = orbitals[i]['ev']
                if i + 1 < len(orbitals) and orbitals[i+1]['occ'] < 0.1:
                    results["lumo"] = orbitals[i+1]['ev']
                    results["gap"] = results["lumo"] - results["homo"]
                    break

    # 4. 偶极矩 - 取最后一次匹配
    dip_matches = re.findall(r'Magnitude \(Debye\)\s+:\s+([\d\.]+)', text)
    if dip_matches:
        results["dipole"] = float(dip_matches[-1])

    # 5. 热力学量 (通常在 Freq 任务最后只出现一次)
    g_match = re.search(r'Final Gibbs free energy\s+\.\.\.\s+([-]?[\d\.]+)\s+Eh', text)
    if g_match: results["gibbs"] = float(g_match.group(1))
    
    h_match = re.search(r'Total enthalpy\s+\.\.\.\s+([-]?[\d\.]+)\s+Eh', text)
    if h_match: results["enthalpy"] = float(h_match.group(1))
    
    s_match = re.search(r'Final entropy term\s+\.\.\.\s+([\d\.]+)\s+Eh', text)
    if s_match: results["entropy"] = float(s_match.group(1))

    # 6. 虚频检测
    freq_start = text.rfind("VIBRATIONAL FREQUENCIES")
    if freq_start != -1:
        freq_block = text[freq_start:freq_start+1500]
        if re.search(r'\d+:\s+-\d+\.\d+\s+cm\*\*-1', freq_block):
            results["imaginary"] = "Yes"

    return results

def process_file(log_path, out_dir):
    basename = os.path.basename(log_path)
    target_name = os.path.splitext(basename)[0] + ".csv"
    out_path = os.path.join(out_dir, target_name)

    try:
        with io.open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
        
        d = parse_ground_state_properties(content)
        if d["homo"] is None: return None

        # 写入详细 CSV
        rows = [
            ["Property", "Value", "Unit"],
            ["SourceFile", basename, "-"],
            ["Total_Energy", "{:.12f}".format(d["total_energy"]) if d["total_energy"] else "N/A", "Eh"],
            ["Nuclear_Repulsion", "{:.12f}".format(d["nuclear_repulsion"]) if d["nuclear_repulsion"] else "N/A", "Eh"],
            ["Electronic_Energy", "{:.12f}".format(d["electronic_energy"]) if d["electronic_energy"] else "N/A", "Eh"],
            ["One_Electron_Energy", "{:.12f}".format(d["one_electron_energy"]) if d["one_electron_energy"] else "N/A", "Eh"],
            ["Two_Electron_Energy", "{:.12f}".format(d["two_electron_energy"]) if d["two_electron_energy"] else "N/A", "Eh"],
            ["CPCM_Solvation_Energy", "{:.12f}".format(d["cpcm_dielectric"]) if d["cpcm_dielectric"] else "N/A", "Eh"],
            ["HOMO", "{:.6f}".format(d["homo"]) if d["homo"] else "N/A", "eV"],
            ["LUMO", "{:.6f}".format(d["lumo"]) if d["lumo"] else "N/A", "eV"],
            ["Energy_Gap", "{:.6f}".format(d["gap"]) if d["gap"] else "N/A", "eV"],
            ["Dipole_Moment", "{:.6f}".format(d["dipole"]) if d["dipole"] else "N/A", "Debye"],
            ["Gibbs_Free_Energy", "{:.12f}".format(d["gibbs"]) if d["gibbs"] else "N/A", "Eh"],
            ["Enthalpy", "{:.12f}".format(d["enthalpy"]) if d["enthalpy"] else "N/A", "Eh"],
            ["Entropy_TS", "{:.12f}".format(d["entropy"]) if d["entropy"] else "N/A", "Eh"],
            ["Has_Imaginary_Freq", d["imaginary"], "Bool"]
        ]
        with io.open(out_path, "w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerows(rows)
        
        return {
            "filename": os.path.splitext(basename)[0],
            "Total_Energy(Eh)": d["total_energy"],
            "Nuclear_Repulsion(Eh)": d["nuclear_repulsion"],
            "Solvation(Eh)": d["cpcm_dielectric"],
            "HOMO(eV)": d["homo"],
            "LUMO(eV)": d["lumo"],
            "Gap(eV)": d["gap"],
            "Dipole(D)": d["dipole"],
            "Gibbs(Eh)": d["gibbs"],
            "Imaginary": d["imaginary"]
        }
    except: return None

def main():
    parser = argparse.ArgumentParser(description="ORCA 性质提取 V12 - 确保提取最后一次收敛结果")
    parser.add_argument("input", help="输入目录 (如 pic3inp)")
    parser.add_argument("-o", "--output", help="输出目录 (可选)")
    args = parser.parse_args()

    in_path = args.input
    if not os.path.isdir(in_path):
        print("Error: Input directory not found.")
        sys.exit(1)

    input_abs = os.path.abspath(in_path)
    input_basename = os.path.basename(input_abs)
    input_parent = os.path.dirname(input_abs)
    
    # pic3inp -> pic3inpout
    output_dir = args.output if args.output else os.path.join(input_parent, input_basename + "out")
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    # 汇总表名：pic3inpout_opt.csv，放在 pic3inpout 文件夹内
    summary_csv_name = os.path.basename(output_dir) + "_opt.csv"
    summary_path = os.path.join(output_dir, summary_csv_name)

    log_files = [f for f in glob.glob(os.path.join(in_path, "*.log")) if not f.lower().endswith("td.log")]

    summary_rows = []
    print("Processing {} ground-state files...".format(len(log_files)))
    for f in sorted(log_files):
        res = process_file(f, output_dir)
        if res:
            summary_rows.append(res)
            print("  [Converged] {}".format(os.path.basename(f)))

    if summary_rows:
        header = ["filename", "Total_Energy(Eh)", "Nuclear_Repulsion(Eh)", "Solvation(Eh)", "HOMO(eV)", "LUMO(eV)", "Gap(eV)", "Dipole(D)", "Gibbs(Eh)", "Imaginary"]
        with io.open(summary_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(summary_rows)
        print("\nAll data successfully saved to: {}".format(output_dir))

if __name__ == "__main__":
    main()