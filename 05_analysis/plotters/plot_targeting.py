#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3.6 (Targeting-relevant analysis): charge-driven mitochondrial targeting proxy vs physicochemical balance

This script follows (and reuses) your existing "read-many-CSV + pictos label mapping" pattern,
but generates the recommended Fig 3.6 a–d and an additional multi-panel figure.

Input CSV must contain (at least):
  - charge, LogP, TPSA, Dipole(D)  (optional: Molecular_Weight, SA_Score, Imaginary)

Outputs (PNG by default; optional PDF):
  - fig3_6a_charge_heatmap.png
  - fig3_6b_positive_fraction.png
  - fig3_6c_logp_tpsa_window.png
  - fig3_6d_charge_vs_tpsa_strip.png
  - fig3_6d2_facets_charge_vs_tpsa.png  (optional, larger)
  - fig3_6_abcd_multipanel.png

Usage:
  python3 plot_fig3_6_targeting.py data/pic3inpout_m.csv data/pic9inpout_m.csv ...

Optional:
  --filter-imaginary            keep only Imaginary == 0 (recommended for publication)
  --logp-min 1 --logp-max 5     set feasible window for LogP
  --tpsa-min 20 --tpsa-max 120  set feasible window for TPSA
  --save-pdf                    also save PDF copies
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap

# 设置标准字体，避免中文字体警告
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False


# -----------------------------
# Label mapping (same style as your script)
# -----------------------------
def load_pictos_labels(pictos_file: Path) -> Dict[str, str]:
    labels: Dict[str, str] = {}
    if not pictos_file.exists():
        print(f"Warning: pictos file not found: {pictos_file}")
        return labels
    try:
        with open(pictos_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or ":" not in line:
                    continue
                pic_name, label = line.split(":", 1)
                labels[pic_name.strip()] = label.strip()
    except Exception as e:
        print(f"Warning: failed to read pictos file: {e}")
    return labels


def get_label_from_filename(filename: str, pictos_labels: Dict[str, str]) -> str:
    pic_name = None
    for pic in pictos_labels.keys():
        if filename.startswith(pic):
            pic_name = pic
            break

    if pic_name and pic_name in pictos_labels:
        label = pictos_labels[pic_name]
        if label.lower().startswith("s"):
            return f"S{label[1:].upper()}"
        return label

    # fallback
    return Path(filename).stem


def sort_labels_by_s_order(labels: List[str]) -> List[str]:
    # 【优化】：实现 S1-S9 顺序排列，并将 pic0/REF 放在最后
    def get_sort_key(label: str):
        u = label.upper()
        # 处理 pic0 或包含 REF 字符的标签
        if "PIC0" in u or "REF" in u:
            return (2, u)  # 优先级最低，排在最后
        if u.startswith("S"):
            try:
                # 提取数字，如 "S1" -> 1
                return (0, int(u[1:]))
            except:
                return (1, u)
        return (1, u)

    return sorted(labels, key=get_sort_key)


def load_reference_points_from_publish(
    script_dir: Path,
    pictos_labels: Dict[str, str]
) -> Dict[str, Dict[str, float]]:
    """从 publishinpdftv22.csv 加载参考点数据
    
    Returns:
        Dict[str, Dict[str, float]]: {batch_label: {'charge': float, 'TPSA': float}}
    """
    import re
    
    # 尝试多个可能的路径
    paths_to_try = [
        script_dir / "data" / "publishinpdftv22.csv",
        script_dir.parent / "data" / "publishinpdftv22.csv",
        Path("data/publishinpdftv22.csv")
    ]
    
    ref_points = {}
    publish_file = None
    for p in paths_to_try:
        if p.exists():
            publish_file = p
            break
    
    if not publish_file:
        print("  ⚠️  未找到 publishinpdftv22.csv 文件")
        return ref_points
    
    try:
        df = pd.read_csv(publish_file)
        for _, row in df.iterrows():
            index_val = str(row.get('index', '')).strip()
            if not index_val or not index_val.startswith('pic'):
                continue
            
            # 提取 pic 序号
            match = re.search(r'pic(\d+)', index_val, re.IGNORECASE)
            if not match:
                continue
            pic_id = match.group(1)
            pic_key = f"pic{pic_id}".lower()
            
            # 通过 pictos 映射获取 batch 标签
            raw_label = pictos_labels.get(pic_key, None)
            if not raw_label:
                continue
            
            # 格式化标签（与 get_label_from_filename 保持一致）
            if raw_label.lower().startswith("s"):
                batch_label = f"S{raw_label[1:].upper()}"
            else:
                batch_label = raw_label
            
            # 提取 charge、TPSA 和 LogP
            charge = row.get('charge', np.nan)
            tpsa = row.get('TPSA', np.nan)
            logp = row.get('LogP', np.nan)
            
            # 构建参考点数据（charge 和 TPSA 必须存在，LogP 可选）
            point_data = {}
            if pd.notna(charge):
                try:
                    point_data['charge'] = float(charge)
                except (ValueError, TypeError):
                    pass
            if pd.notna(tpsa):
                try:
                    point_data['TPSA'] = float(tpsa)
                except (ValueError, TypeError):
                    pass
            if pd.notna(logp):
                try:
                    point_data['LogP'] = float(logp)
                except (ValueError, TypeError):
                    pass
            
            # 至少需要有 charge 或 TPSA 或 LogP 中的一个
            if point_data:
                ref_points[batch_label] = point_data
                charge_str = f"charge={point_data.get('charge', 'N/A')}" if 'charge' in point_data else ""
                tpsa_str = f"TPSA={point_data.get('TPSA', 'N/A'):.1f}" if 'TPSA' in point_data else ""
                logp_str = f"LogP={point_data.get('LogP', 'N/A'):.2f}" if 'LogP' in point_data else ""
                info_list = [s for s in [charge_str, tpsa_str, logp_str] if s]
                print(f"  ✅ 已加载参考点 {index_val} -> {batch_label}: {', '.join(info_list)}")
    except Exception as e:
        print(f"  ⚠️  加载参考点数据时出错: {e}")
    
    return ref_points


# -----------------------------
# Data loading
# -----------------------------
def _to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def read_and_concat_csvs(
    csv_files: List[Path],
    pictos_labels: Dict[str, str],
    filter_imaginary: bool = False,
) -> Tuple[pd.DataFrame, List[str]]:
    frames: List[pd.DataFrame] = []
    batch_labels: List[str] = []

    for f in csv_files:
        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f"Warning: cannot read {f}: {e}")
            continue

        batch = get_label_from_filename(f.name, pictos_labels)
        batch_labels.append(batch)
        df = df.copy()
        df["Batch"] = batch  # unify group key
        # 标记是否为背景参照组
        df["Is_Ref"] = True if ("PIC0" in f.name.upper() or "REF" in batch.upper()) else False

        # numeric conversions
        for col in ["charge", "LogP", "TPSA", "Dipole(D)", "Molecular_Weight", "SA_Score"]:
            if col in df.columns:
                df[col] = _to_num(df[col])

        # enforce charge as integer-like where possible
        if "charge" in df.columns:
            df["charge"] = np.floor(df["charge"].astype(float) + 1e-9)  # tolerate 1.0 etc.
            df.loc[(df["charge"] < 0) | (df["charge"] > 50), "charge"] = np.nan  # sanity

        if filter_imaginary and "Imaginary" in df.columns:
            # ORCA parser output is Yes/No; accept legacy numeric/bool forms too.
            # Missing or unrecognized values are excluded from this quality gate.
            imag_raw = df["Imaginary"]
            imag_text = imag_raw.astype(str).str.strip().str.lower()
            imag_num = pd.to_numeric(imag_raw, errors="coerce")
            no_imaginary = (
                imag_text.isin({"no", "false", "0", "0.0"})
                | imag_num.eq(0.0)
            )
            df = df.loc[no_imaginary].copy()

        frames.append(df)

    if not frames:
        return pd.DataFrame(), []

    all_df = pd.concat(frames, ignore_index=True)

    # Determine plotting order from the batches we actually have
    uniq_batches = list(pd.unique(all_df["Batch"]))
    ordered = sort_labels_by_s_order(uniq_batches)
    return all_df, ordered


# -----------------------------
# Plot helpers
# -----------------------------
def _save(fig: plt.Figure, out: Path, save_pdf: bool = False) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    if save_pdf:
        fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}" + (" (+PDF)" if save_pdf else ""))


def _charge_cmap_norm(charge_max: int = 6):
    # 0..charge_max discrete colormap
    base = plt.get_cmap("viridis", charge_max + 1)
    colors = [base(i) for i in range(charge_max + 1)]
    cmap = ListedColormap(colors)
    boundaries = np.arange(-0.5, charge_max + 1.5, 1.0)
    norm = BoundaryNorm(boundaries, cmap.N)
    return cmap, norm


# -----------------------------
# Fig 3.6a: Heatmap of charge composition per Batch
# -----------------------------
def plot_36a_charge_heatmap(
    df: pd.DataFrame,
    batch_order: List[str],
    charge_max: int,
    out: Path,
    save_pdf: bool = False,
    pictos_labels: Optional[Dict[str, str]] = None,
    script_dir: Optional[Path] = None,
) -> None:
    if df.empty or "charge" not in df.columns:
        print("Skip 3.6a: empty data or missing 'charge'.")
        return

    charges = list(range(0, charge_max + 1))

    # compute proportions
    mat = np.zeros((len(charges), len(batch_order)), dtype=float)
    cnt = np.zeros((len(charges), len(batch_order)), dtype=int)

    for j, b in enumerate(batch_order):
        sub = df.loc[df["Batch"] == b, "charge"].dropna()
        n = len(sub)
        if n == 0:
            continue
        for i, c in enumerate(charges):
            k = int((sub == c).sum())
            cnt[i, j] = k
            mat[i, j] = k / n

    # 加载参考点数据（从 publishinpdftv22.csv）
    ref_points = {}
    if pictos_labels and script_dir:
        ref_points = load_reference_points_from_publish(script_dir, pictos_labels)
        if ref_points:
            print(f"  [3.6a] 已加载 {len(ref_points)} 个参考点用于标记")

    fig, ax = plt.subplots(figsize=(1.0 * max(8, len(batch_order)), 8))
    
    # 【优化】：美化热图，0电荷对应白色
    colors = ["#ffffff", "#4da6ff", "#004d99"] # 白 -> 蓝渐变
    custom_cmap = LinearSegmentedColormap.from_list("charge_heat", colors)
    im = ax.imshow(mat, aspect="auto", origin="lower", cmap=custom_cmap)

    ax.set_xticks(np.arange(len(batch_order)))
    ax.set_xticklabels(batch_order, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(charges)))
    ax.set_yticklabels([str(c) for c in charges])
    # ax.set_xlabel("Scaffold / Batch", fontweight="bold")
    ax.set_ylabel("Charge Count", fontweight="bold")
    ax.set_title("Charge composition within each scaffold", fontweight="bold", pad=12)

    # 获取主 axes 的位置，创建 colorbar axes 使高度与 y 轴一致
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.02, pos.y0, 0.02, pos.height])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Proportion within scaffold", fontweight="bold")

    # annotate with %
    for i in range(len(charges)):
        for j in range(len(batch_order)):
            if cnt[i, j] == 0 or mat[i, j] == 0:
                continue
            # 如果比例小于1%但大于0，至少显示一个数字（显示为"<1%"）
            if mat[i, j] < 0.01:
                text = "<1%"
            else:
                text = f"{mat[i, j]*100:.0f}%"
            ax.text(j, i, text,
                    ha="center", va="center",
                    fontsize=11,
                    color="white" if mat[i, j] > 0.45 else "black")
    
    # 标记参考点（红色五角星）- 在对应的 batch 列和 charge 行位置
    ref_point_color = "#e74c3c"
    ref_point_label_added = False
    for j, batch_name in enumerate(batch_order):
        if batch_name in ref_points:
            ref_pt = ref_points[batch_name]
            ref_charge = ref_pt.get('charge', None)
            if ref_charge is not None and 0 <= int(ref_charge) <= charge_max:
                charge_idx = int(ref_charge)
                # 找到该 charge 在 charges 列表中的索引
                if charge_idx in charges:
                    i = charges.index(charge_idx)
                    # 在热图的对应位置标记红色五角星（稍微向右偏移，避免挡住百分比）
                    label = 'Reference Point' if not ref_point_label_added else ''
                    ax.scatter(j + 0.35, i, s=80, color=ref_point_color, 
                             edgecolors='black', linewidth=0.8, marker='*', 
                             zorder=10, label=label)
                    ref_point_label_added = True
    
    # 如果有参考点，添加图例
    if ref_point_label_added:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles=handles, labels=labels, loc='upper right', fontsize=9, framealpha=0.9)

    ax.grid(False)
    _save(fig, out, save_pdf=save_pdf)


# -----------------------------
# Fig 3.6b: Positive fraction (charge > 0) per Batch
# -----------------------------
def plot_36b_positive_fraction(
    df: pd.DataFrame,
    batch_order: List[str],
    out: Path,
    save_pdf: bool = False,
    positive_threshold: int = 0,
    sort_desc: bool = False,
    pictos_labels: Optional[Dict[str, str]] = None,
    script_dir: Optional[Path] = None,
) -> None:
    """
    Fig 3.6b | 阳离子比例累积堆叠图 (支持 Charge 0-6)
    优化：0电荷为灰色，阳离子总量为0时不标注百分比，美化图表宽度。
    """
    if df.empty or "charge" not in df.columns:
        print("Skip 3.6b: 缺少电荷数据。")
        return

    # 设置最大电荷为 6
    charge_max = 6 
    data_rows = []
    
    # 1. 提取各批次电荷分布数据
    for b in batch_order:
        sub = df.loc[df["Batch"] == b, "charge"].dropna()
        n = len(sub)
        if n == 0: continue
        
        row = {"Batch": b}
        # 计算 0 到 6 的分布比例
        for c in range(0, charge_max + 1):
            row[f"C{c}"] = (sub == c).sum() / n
        # 计算总阳离子比例 (Charge > 0)
        row["pos_total"] = (sub > 0).sum() / n
        data_rows.append(row)

    stat = pd.DataFrame(data_rows)

    # 2. 加载参考点数据（从 publishinpdftv22.csv）
    ref_points = {}
    if pictos_labels and script_dir:
        ref_points = load_reference_points_from_publish(script_dir, pictos_labels)
        if ref_points:
            print(f"  [3.6b] 已加载 {len(ref_points)} 个参考点用于标记")

    # 3. 动态调整画布宽度，防止过宽 (每个柱子约 0.6 英寸)
    fig_width = max(5, len(stat) * 0.6 + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, 5), dpi=300)
    
    # 4. 配色方案：Charge 0 为灰色，1-6 使用 Pastel1 离散浅色（非常浅和柔和）
    cmap = plt.get_cmap("Pastel1")
    # colors[0] 为灰色，其余为 Pastel1 映射色（非常浅的离散颜色）
    colors = ["#F2F2F2"] + [cmap(i % cmap.N) for i in range(1, charge_max + 1)]
    ref_point_color = "#e74c3c"  # 红色参考点
    
    # 5. 绘制堆叠柱状图，同时记录每个 charge 的底部位置（用于标记参考点）
    bottom = np.zeros(len(stat))
    bottom_positions = {}  # {batch_name: {charge: bottom_y}}
    batch_total_heights = {}  # {batch_name: total_height}
    
    for c in range(0, charge_max + 1):
        label = f"Charge +{c}" if c > 0 else "Neutral (0)"
        # 增加极细白边分割色块
        ax.bar(stat["Batch"], stat[f"C{c}"], bottom=bottom, color=colors[c], 
               edgecolor="white", linewidth=0.5, label=label)
        
        # 记录每个 batch 的每个 charge 的底部位置（用于标记参考点）
        for i, batch_name in enumerate(stat["Batch"]):
            if batch_name not in bottom_positions:
                bottom_positions[batch_name] = {}
            bottom_positions[batch_name][c] = bottom[i]
        
        bottom += stat[f"C{c}"]
    
    # 计算每个 batch 的柱子总高度（所有 charge 段的总和）
    for i, batch_name in enumerate(stat["Batch"]):
        batch_total_heights[batch_name] = bottom[i]

    # 6. 标记参考点（红色点）- 在对应 charge 堆叠段的中心位置
    ref_point_label_added = False
    for i, batch_name in enumerate(stat["Batch"]):
        if batch_name in ref_points:
            ref_pt = ref_points[batch_name]
            ref_charge = ref_pt.get('charge', None)
            if ref_charge is not None and 0 <= int(ref_charge) <= charge_max:
                charge_idx = int(ref_charge)
                # 获取该 charge 堆叠段的底部和高度
                if batch_name in bottom_positions and charge_idx in bottom_positions[batch_name]:
                    segment_bottom = bottom_positions[batch_name][charge_idx]
                    segment_height = stat.loc[stat["Batch"] == batch_name, f"C{charge_idx}"].values[0]
                    
                    # 确保堆叠段高度大于 0（即该 batch 有该 charge 的数据）
                    if segment_height > 0:
                        # 在堆叠段的中心位置标记红色五角星
                        y_pos = segment_bottom + segment_height / 2
                        # 只为第一个参考点添加 label，用于图例
                        label = 'Reference Point' if not ref_point_label_added else ''
                        ax.scatter(i, y_pos, s=80, color=ref_point_color, 
                                 edgecolors='black', linewidth=0.8, marker='*', 
                                 zorder=10, label=label)
                        ref_point_label_added = True
                    else:
                        # 如果该 batch 没有该 charge 的数据，在柱子顶部标记
                        total_height = batch_total_heights.get(batch_name, 1.0)
                        y_pos = total_height + 0.05  # 稍微高于柱子顶部
                        label = 'Reference Point' if not ref_point_label_added else ''
                        ax.scatter(i, y_pos, s=80, color=ref_point_color, 
                                 edgecolors='black', linewidth=0.8, marker='*', 
                                 zorder=10, label=label)
                        ref_point_label_added = True

    # 7. 在柱子顶部标注阳离子总比例 (%)
    # 逻辑：如果 val == 0 (即只有灰色部分)，则不进行标注
    for i, val in enumerate(stat["pos_total"]):
        if val > 0:  # 只有当存在阳离子时才标注百分比
            ax.text(i, 1.02, f"{val*100:.1f}%", ha='center', va='bottom', 
                    fontsize=8.5, fontweight='bold', color="#2c3e50")

    # 8. 学术化视觉美化
    ax.set_ylim(0, 1.15) # 预留顶部百分比空间
    ax.set_ylabel("Population Fraction", fontweight="bold", fontsize=10)
    ax.set_title("Fig 3.6b | Targeted Charge Distribution Composition", 
                 fontweight="bold", fontsize=12, pad=20)
    
    # 图例设置：放在右侧，无边框
    ax.legend(title="Charge State", bbox_to_anchor=(1.02, 1), loc='upper left', 
              fontsize=8.5, frameon=False)
    
    # 隐藏上方和右侧边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#cccccc')
    ax.spines['bottom'].set_color('#cccccc')

    # 辅助网格线
    ax.grid(axis="y", alpha=0.3, linestyle="--", color='#bdc3c7', zorder=0)
    ax.set_axisbelow(True)
    
    # 刻度美化
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(np.arange(0, 1.1, 0.2))

    plt.tight_layout()
    _save(fig, out, save_pdf=save_pdf)


# -----------------------------
# Fig 3.6c: LogP vs TPSA scatter with feasible window; color by charge
# -----------------------------
def plot_36c_logp_tpsa_window(
    df: pd.DataFrame,
    batch_order: List[str],
    out: Path,
    save_pdf: bool = False,
    charge_max: int = 6,
    logp_window: Tuple[float, float] = (1.0, 5.0),
    tpsa_window: Tuple[float, float] = (20.0, 120.0),
    script_dir: Optional[Path] = None,
    pictos_labels: Optional[Dict[str, str]] = None,
) -> None:
    need_cols = ["charge", "LogP", "TPSA", "Batch"]
    for c in need_cols:
        if c not in df.columns:
            print(f"Skip 3.6c: missing column '{c}'.")
            return

    # =================================================================
    # 第一部分：完全保留您原来的逻辑，生成总图 (Total Plot)
    # =================================================================
    df_ref = df[df["Is_Ref"] == True].dropna(subset=["LogP", "TPSA"])
    df_gen = df[df["Is_Ref"] == False].dropna(subset=["LogP", "TPSA", "charge"])

    # 为总图使用 tab10 colormap
    base_tab10 = plt.get_cmap("tab10")
    tab10_colors = [base_tab10(i % base_tab10.N) for i in range(charge_max + 1)]
    cmap_tab10 = ListedColormap(tab10_colors)
    boundaries = np.arange(-0.5, charge_max + 1.5, 1.0)
    norm_tab10 = BoundaryNorm(boundaries, cmap_tab10.N)
    
    # 【备用方案】最初的 viridis colormap（取消注释下面的代码，注释掉上面的 tab10 代码即可切换）
    # cmap_total, norm_total = _charge_cmap_norm(charge_max=charge_max)
    
    # 其他图使用原来的 colormap
    cmap, norm = _charge_cmap_norm(charge_max=charge_max)

    fig, ax = plt.subplots(figsize=(7.5, 6))
    if not df_ref.empty:
        ax.scatter(df_ref["LogP"], df_ref["TPSA"], color="lightgrey", s=12, alpha=0.3, label="RefDataset", zorder=1)

    if not df_gen.empty:
        # 使用 tab10 colormap
        sc = ax.scatter(df_gen["LogP"], df_gen["TPSA"], c=df_gen["charge"], cmap=cmap_tab10, norm=norm_tab10,
                        s=25, alpha=0.5, edgecolors="none", zorder=2)
        # 【备用方案】使用 viridis colormap（取消注释下面两行，注释掉上面两行即可切换）
        # sc = ax.scatter(df_gen["LogP"], df_gen["TPSA"], c=df_gen["charge"], cmap=cmap_total, norm=norm_total,
        #                 s=25, alpha=0.5, edgecolors="none", zorder=2)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.01, ticks=np.arange(0, charge_max + 1))
        cbar.set_label("charge count", fontweight="bold")

    ref_points_dict = {}
    ref_points_list = []
    if script_dir and pictos_labels:
        paths_to_try = [
            script_dir / "data" / "publishinpdftv22.csv",
            script_dir.parent / "data" / "publishinpdftv22.csv",
            Path("data/publishinpdftv22.csv")
        ]
        publish_file = None
        for p in paths_to_try:
            if p.exists():
                publish_file = p
                break
        if publish_file:
            try:
                import re
                ref_df = pd.read_csv(publish_file)
                for _, row in ref_df.iterrows():
                    index_val = str(row.get('index', '')).strip()
                    if not index_val or not index_val.startswith('pic'):
                        continue
                    match = re.search(r'pic(\d+)', index_val, re.IGNORECASE)
                    if not match: continue
                    pic_id = match.group(1); pic_key = f"pic{pic_id}".lower()
                    raw_label = pictos_labels.get(pic_key, None)
                    if not raw_label: continue
                    if raw_label.lower().startswith("s"): batch_label = f"S{raw_label[1:].upper()}"
                    else: batch_label = raw_label
                    logp_val = row.get('LogP', np.nan); tpsa_val = row.get('TPSA', np.nan)
                    if pd.notna(logp_val) and pd.notna(tpsa_val):
                        point_data = {'LogP': float(logp_val), 'TPSA': float(tpsa_val)}
                        ref_points_dict[batch_label] = point_data
                        ref_points_list.append(point_data)
                if ref_points_list:
                    ax.scatter([pt['LogP'] for pt in ref_points_list], [pt['TPSA'] for pt in ref_points_list], 
                               s=80, color='#e74c3c', edgecolors='black', linewidth=0.8, marker='*', zorder=10, label='Reference Points')
            except Exception as e: print(f"  ⚠️ [3.6c] 加载参考点出错: {e}")

    ax.set_xlabel("LogP", fontweight="bold"); ax.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold")
    ax.set_title("LogP vs TPSA", fontweight="bold", pad=12); ax.grid(alpha=0.25, linestyle="--"); ax.set_axisbelow(True)
    handles, labels = ax.get_legend_handles_labels()
    if handles: ax.legend(handles=handles, labels=labels, loc='upper right', fontsize=9, framealpha=0.9)
    _save(fig, out, save_pdf=save_pdf)

    # =================================================================
    # 第二部分：批量生成额外图片 (保留逻辑)
    # =================================================================
    print(f"  [3.6c] 开始生成额外图片...")
    # A. RefDataset 图
    if not df_ref.empty:
        fig_ref, ax_ref = plt.subplots(figsize=(7.5, 6))
        sc_ref = ax_ref.scatter(df_ref["LogP"], df_ref["TPSA"], c=df_ref["charge"], cmap=cmap, norm=norm, s=20, alpha=0.8, linewidth=0.5, zorder=5)
        fig_ref.colorbar(sc_ref, ax=ax_ref, fraction=0.04, pad=0.02, ticks=np.arange(0, charge_max + 1)).set_label("charge", fontweight="bold")
        ax_ref.set_xlabel("LogP", fontweight="bold"); ax_ref.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold")
        ax_ref.set_title("RefDataset (Baseline)", fontweight="bold"); ax_ref.grid(alpha=0.25, linestyle="--")
        _save(fig_ref, out.parent / f"{out.stem}_RefOnly{out.suffix}", save_pdf=save_pdf)

    # B. Generation (All) 图
    if not df_gen.empty:
        fig_gen, ax_gen = plt.subplots(figsize=(7.5, 6))
        sc_gen = ax_gen.scatter(df_gen["LogP"], df_gen["TPSA"], c=df_gen["charge"], cmap=cmap, norm=norm, s=25, alpha=0.7, edgecolors="none", zorder=2)
        fig_gen.colorbar(sc_gen, ax=ax_gen, fraction=0.04, pad=0.02, ticks=np.arange(0, charge_max + 1)).set_label("charge (integer)", fontweight="bold")
        if ref_points_list:
            ax_gen.scatter([pt['LogP'] for pt in ref_points_list], [pt['TPSA'] for pt in ref_points_list], s=80, color='#e74c3c', edgecolors='black', marker='*', zorder=10, label='Reference Points')
        ax_gen.set_xlabel("LogP", fontweight="bold"); ax_gen.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold")
        ax_gen.set_title("Generation (All)", fontweight="bold"); ax_gen.grid(alpha=0.25, linestyle="--"); ax_gen.set_axisbelow(True)
        _save(fig_gen, out.parent / f"{out.stem}_GenOnly{out.suffix}", save_pdf=save_pdf)

    # C. 骨架循环图
    for batch in batch_order:
        if batch == "RefDaset" or "Ref" in batch: continue
        df_batch = df[df["Batch"] == batch].dropna(subset=["LogP", "TPSA"])
        if df_batch.empty: continue
        fig_b, ax_b = plt.subplots(figsize=(7.5, 6))
        if not df_ref.empty: ax_b.scatter(df_ref["LogP"], df_ref["TPSA"], color="lightgrey", s=12, alpha=0.3, label="RefDataset", zorder=1)
        sc_b = ax_b.scatter(df_batch["LogP"], df_batch["TPSA"], c=df_batch["charge"], cmap=cmap, norm=norm, s=40, alpha=0.8, edgecolors="white", linewidth=0.5, zorder=5)
        fig_b.colorbar(sc_b, ax=ax_b, fraction=0.04, pad=0.02, ticks=np.arange(0, charge_max + 1)).set_label("charge", fontweight="bold")
        if batch in ref_points_dict:
            ref_pt = ref_points_dict[batch]
            ax_b.scatter(ref_pt['LogP'], ref_pt['TPSA'], s=80, color='#e74c3c', edgecolors='black', marker='*', zorder=10, label='Reference Point')
        ax_b.set_xlabel("LogP", fontweight="bold"); ax_b.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold")
        ax_b.set_title(f"Scaffold {batch} balance", fontweight="bold"); ax_b.grid(alpha=0.25, linestyle="--"); ax_b.legend(loc='upper right', fontsize=8)
        _save(fig_b, out.parent / f"{out.stem}_{batch}{out.suffix}", save_pdf=save_pdf)

    # =================================================================
    # E. 新增：三维理化性质展示 (X=LogP, Y=TPSA, Z=Charge)
    # =================================================================
    if not df_gen.empty:
        from mpl_toolkits.mplot3d import Axes3D
        print(f"  [3.6c] 正在生成三维理化空间展示图...")
        fig_3d = plt.figure(figsize=(10, 8))
        ax_3d = fig_3d.add_subplot(111, projection='3d')
        
        # 绘制三维散点
        # 使用相同的 cmap 和 norm 保证颜色含义一致
        sc_3d = ax_3d.scatter(
            df_gen["LogP"], 
            df_gen["TPSA"], 
            df_gen["charge"], 
            c=df_gen["charge"], 
            cmap=cmap, 
            norm=norm,
            s=35, 
            alpha=0.7, 
            edgecolors='white', 
            linewidth=0.3
        )
        
        # 设置轴标签
        ax_3d.set_xlabel("LogP", fontweight="bold", labelpad=10)
        ax_3d.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold", labelpad=10)
        ax_3d.set_zlabel("Charge Count", fontweight="bold", labelpad=10)
        ax_3d.set_title("3D Physicochemical Distribution", fontweight="bold", pad=20)
        
        # 调整刻度
        ax_3d.set_zticks(np.arange(0, charge_max + 1))
        
        # 角度设置：仰角调低（向下一点），方位角稍微旋转
        # elev=15 使视角更平视，能更好地看清电荷层级
        ax_3d.view_init(elev=15, azim=-60)
        
        # 添加颜色条
        cax_3d = fig_3d.add_axes([0.88, 0.25, 0.02, 0.5]) # 手动调整颜色条位置避免重叠
        fig_3d.colorbar(sc_3d, cax=cax_3d, ticks=np.arange(0, charge_max + 1)).set_label("Charge", fontweight="bold")
        
        _save(fig_3d, out.parent / f"{out.stem}_3D_View{out.suffix}", save_pdf=save_pdf)


# -----------------------------
# Fig 3.6d: charge vs TPSA (strip with jitter); color by Batch
# -----------------------------
def plot_36d_charge_vs_tpsa_strip(
    df: pd.DataFrame,
    batch_order: List[str],
    out: Path,
    save_pdf: bool = False,
    charge_max: int = 6,
    jitter: float = 0.18,
    seed: int = 0,
) -> None:
    need_cols = ["charge", "TPSA", "Batch"]
    for c in need_cols:
        if c not in df.columns:
            print(f"Skip 3.6d: missing column '{c}'.")
            return

    sub = df[need_cols].dropna()
    if sub.empty:
        print("Skip 3.6d: no valid rows after dropna.")
        return

    # keep charge range
    sub = sub.loc[(sub["charge"] >= 0) & (sub["charge"] <= charge_max)].copy()
    if sub.empty:
        print("Skip 3.6d: no rows within charge range.")
        return

    rng = np.random.default_rng(seed)

    # color per batch (9 batches): use tab10 for readability
    tab = plt.get_cmap("tab10")
    batch_to_color = {b: tab(i % 10) for i, b in enumerate(batch_order)}

    fig, ax = plt.subplots(figsize=(7.5, 6))
    for b in batch_order:
        sb = sub.loc[sub["Batch"] == b]
        if sb.empty:
            continue
        x = sb["charge"].values.astype(float)
        xj = x + rng.uniform(-jitter, jitter, size=len(x))
        ax.scatter(
            xj,
            sb["TPSA"].values,
            s=14,
            alpha=0.55,
            color=batch_to_color[b],
            edgecolors="none",
            label=b,
        )

    ax.set_xlim(-0.6, charge_max + 0.6)
    ax.set_xticks(np.arange(0, charge_max + 1))
    ax.set_xlabel("charge", fontweight="bold")
    ax.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold")
    ax.set_title("Charge vs TPSA (jittered strip), colored by scaffold batch", fontweight="bold", pad=12)
    ax.grid(alpha=0.25, linestyle="--")
    ax.set_axisbelow(True)

    # compact legend
    handles = [mlines.Line2D([], [], color=batch_to_color[b], marker="o", linestyle="None", markersize=6, label=b)
               for b in batch_order if b in batch_to_color]
    ax.legend(handles=handles, loc="upper right", framealpha=0.9, fontsize=9)

    _save(fig, out, save_pdf=save_pdf)


# -----------------------------
# Optional: Fig 3.6d2 facets (3x3): charge vs TPSA per Batch (jitter + box)
# -----------------------------
def plot_36d2_facets_charge_vs_tpsa(
    df: pd.DataFrame,
    batch_order: List[str],
    out: Path,
    save_pdf: bool = False,
    charge_max: int = 6,
    jitter: float = 0.15,
    seed: int = 0,
    pictos_labels: Optional[Dict[str, str]] = None,
    script_dir: Optional[Path] = None,
) -> None:
    need_cols = ["charge", "TPSA", "Batch"]
    for c in need_cols:
        if c not in df.columns:
            print(f"Skip 3.6d2: 缺少核心列 '{c}'")
            return

    # 1. 核心标记逻辑
    if "Is_Ref" not in df.columns:
        df["Is_Ref"] = df["Batch"].str.contains("Ref|pic0", case=False, na=False)

    # 2. 数据清洗：保留有效数值
    sub = df.dropna(subset=["charge", "TPSA"]).copy()
    sub = sub.loc[(sub["charge"] >= 0) & (sub["charge"] <= charge_max)].copy()
    
    if sub.empty:
        print("Skip 3.6d2: 清洗后无有效数据。请检查 pic0 是否已添加 charge 列。")
        return

    # 3. 构建绘图序列 (保持 S1-S9 -> RefDaset -> Generation Total)
    all_avail = sub["Batch"].unique()
    gen_individual = [b for b in batch_order if b in all_avail and 
                      not sub.loc[sub["Batch"] == b, "Is_Ref"].any()]
    ref_group = [b for b in all_avail if sub.loc[sub["Batch"] == b, "Is_Ref"].any()]
    gen_total_label = "Generation (All)"
    plot_sequence = gen_individual + ref_group + [gen_total_label]

    # 4. 画布配置
    n = len(plot_sequence)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4.2 * nrows), sharex=False, sharey=True)
    axes = np.array(axes).reshape(nrows, ncols)

    rng = np.random.default_rng(seed)
    x_ticks = np.arange(0, charge_max + 1)
    x_lim = (-0.6, charge_max + 0.6)
    
    # 5. 加载参考点数据（从 publishinpdftv22.csv）
    ref_points = {}
    if pictos_labels and script_dir:
        print("\n--- 正在加载参考点数据 (publishinpdftv22.csv) ---")
        ref_points = load_reference_points_from_publish(script_dir, pictos_labels)
        if ref_points:
            print(f"  已加载 {len(ref_points)} 个参考点")

    # 6. 统一配色方案 (所有子图共用一套颜色)
    # 使用统一的蓝色系，或者你可以根据需要更换为 plt.cm.viridis
    standard_box_color = "#74b9ff"  # 浅蓝色填充
    standard_pt_color = "#0984ff"   # 深蓝色点
    median_line_color = "#d63031"   # 红色中位线
    ref_point_color = "#e74c3c"     # 红色参考点

    for idx, b in enumerate(plot_sequence):
        r, c = divmod(idx, ncols)
        ax = axes[r, c]
        
        # 数据选择
        if b == gen_total_label:
            sb = sub.loc[sub["Is_Ref"] == False]
        else:
            sb = sub.loc[sub["Batch"] == b]

        if sb.empty:
            ax.axis("off")
            continue

        # 6. 绘制抖动散点 (Strip)
        x_raw = sb["charge"].values.astype(float)
        xj = x_raw + rng.uniform(-jitter, jitter, size=len(x_raw))
        ax.scatter(xj, sb["TPSA"].values, s=15, alpha=0.4, color=standard_pt_color, 
                   edgecolors="none", zorder=2)

        # 7. 绘制箱线图 (Boxplot)
        data_list, pos_list = [], []
        for ch in x_ticks:
            vals = sb.loc[sb["charge"] == ch, "TPSA"].values
            if len(vals) > 0:
                data_list.append(vals)
                pos_list.append(ch)
        
        if data_list:
            bp = ax.boxplot(data_list, positions=pos_list, widths=0.45, patch_artist=True, 
                            showfliers=False, zorder=3, 
                            medianprops=dict(color=median_line_color, linewidth=1.5))
            for box in bp["boxes"]:
                box.set_facecolor(standard_box_color)
                box.set_alpha(0.4)
                box.set_edgecolor("#2d3436")
                box.set_linewidth(0.8)

        # 7.5. 绘制参考点（红色点）- 如果该 batch 有参考点
        if b in ref_points and b != gen_total_label:
            ref_pt = ref_points[b]
            ref_charge = ref_pt.get('charge', None)
            ref_tpsa = ref_pt.get('TPSA', None)
            if ref_charge is not None and ref_tpsa is not None:
                # 确保参考点的 charge 在有效范围内（0 到 charge_max）
                if 0 <= ref_charge <= charge_max and pd.notna(ref_tpsa):
                    ax.scatter(ref_charge, ref_tpsa, s=120, color=ref_point_color, 
                             edgecolors='black', linewidth=1.5, marker='o', 
                             zorder=10, label='Reference Point' if idx == 0 else '')

        # 8. 统一坐标轴与标题
        ax.set_title(b, fontweight="bold", fontsize=13)
        ax.grid(axis='y', alpha=0.2, linestyle="--", zorder=1)
        ax.set_xlim(x_lim)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([str(int(t)) for t in x_ticks])
        ax.set_xlabel("Charge", fontweight="bold")
        if c == 0:
            ax.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold")

    # 关闭空子图
    for k in range(n, nrows * ncols):
        axes.flatten()[k].axis("off")

    fig.suptitle("Charge vs TPSA", 
                 fontweight="bold", fontsize=15, y=0.99)
    
    plt.subplots_adjust(left=0.08, right=0.95, top=0.94, bottom=0.08, wspace=0.25, hspace=0.45)
    _save(fig, out, save_pdf=save_pdf)
# -----------------------------
# Multi-panel: a–d
# -----------------------------
def plot_36_abcd_multipanel(
    df: pd.DataFrame,
    batch_order: List[str],
    out: Path,
    save_pdf: bool = False,
    charge_max: int = 6,
    positive_threshold: int = 0,
    logp_window: Tuple[float, float] = (1.0, 5.0),
    tpsa_window: Tuple[float, float] = (20.0, 120.0),
) -> None:
    if df.empty:
        print("Skip multipanel: empty data.")
        return

    # Prepare essentials
    if "charge" not in df.columns or "Batch" not in df.columns:
        print("Skip multipanel: missing 'charge' or 'Batch'.")
        return

    charges = list(range(0, charge_max + 1))

    # a) heatmap matrix
    mat = np.zeros((len(charges), len(batch_order)), dtype=float)
    for j, b in enumerate(batch_order):
        subc = df.loc[df["Batch"] == b, "charge"].dropna()
        n = len(subc)
        if n == 0:
            continue
        for i, c in enumerate(charges):
            mat[i, j] = float((subc == c).sum()) / n

    # b) 采用堆叠逻辑
    stat_list = []
    for b in batch_order:
        subc = df.loc[df["Batch"] == b, "charge"].dropna()
        n = len(subc)
        val = np.nan if n == 0 else float((subc > positive_threshold).sum()) / n
        stat_list.append({"Batch": b, "pos_frac": val})
    stat = pd.DataFrame(stat_list) # 保持 batch_order 顺序

    cmap, norm = _charge_cmap_norm(charge_max=charge_max)

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, wspace=0.18, hspace=0.22)

    # (a)
    axa = fig.add_subplot(gs[0, 0])
    # 统一美化配色
    im = axa.imshow(mat, aspect="auto", origin="lower", cmap=LinearSegmentedColormap.from_list("w_blue", ["#ffffff", "#004d99"]))
    axa.set_xticks(np.arange(len(batch_order)))
    axa.set_xticklabels(batch_order, rotation=45, ha="right")
    axa.set_yticks(np.arange(len(charges)))
    axa.set_yticklabels([str(c) for c in charges])
    axa.set_title("(a) Charge composition (proportion)", fontweight="bold")
    axa.set_xlabel("Scaffold / Batch", fontweight="bold")
    axa.set_ylabel("charge", fontweight="bold")
    cbar = fig.colorbar(im, ax=axa, fraction=0.03, pad=0.02)
    cbar.set_label("Proportion", fontweight="bold")

    # (b)
    axb = fig.add_subplot(gs[0, 1])
    x = np.arange(len(stat))
    # 为多面板图保持简单性，只画阳离子总量，但保持排序
    axb.bar(x, stat["pos_frac"], edgecolor="black", linewidth=0.8, color="#4da6ff")
    axb.set_xticks(x)
    axb.set_xticklabels(stat["Batch"].tolist(), rotation=45, ha="right")
    axb.set_ylim(0, 1.0)
    axb.set_title("(b) Cationic fraction P(charge>0)", fontweight="bold")
    axb.set_ylabel("P(charge > 0)", fontweight="bold")
    axb.set_xlabel("Scaffold / Batch", fontweight="bold")
    axb.grid(axis="y", alpha=0.3, linestyle="--")
    axb.set_axisbelow(True)

    # (c)
    axc = fig.add_subplot(gs[1, 0])
    df_ref = df[df["Is_Ref"] == True]
    df_gen = df[df["Is_Ref"] == False]
    
    if not df_ref.empty:
        axc.scatter(df_ref["LogP"], df_ref["TPSA"], color="lightgrey", s=10, alpha=0.2, zorder=1)
    
    if not df_gen.empty:
        sc = axc.scatter(df_gen["LogP"], df_gen["TPSA"], c=df_gen["charge"], 
                        cmap=cmap, norm=norm, s=16, alpha=0.55, edgecolors="none", zorder=2)
        cb = fig.colorbar(sc, ax=axc, fraction=0.03, pad=0.02, ticks=np.arange(0, charge_max + 1))
        cb.set_label("charge", fontweight="bold")
        
    x0, x1 = logp_window
    y0, y1 = tpsa_window
    rect = mpatches.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, linewidth=2.0, linestyle="--", edgecolor="black")
    axc.add_patch(rect)
    axc.set_title("(c) LogP vs TPSA (Colored by charge)", fontweight="bold")
    axc.set_xlabel("LogP", fontweight="bold")
    axc.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold")
    axc.grid(alpha=0.25, linestyle="--")

    # (d)
    axd = fig.add_subplot(gs[1, 1])
    if "TPSA" in df.columns:
        subd = df[["charge", "TPSA", "Batch"]].dropna()
        subd = subd.loc[(subd["charge"] >= 0) & (subd["charge"] <= charge_max)]
        rng = np.random.default_rng(0)
        tab = plt.get_cmap("tab10")
        b_to_c = {b: tab(i % 10) for i, b in enumerate(batch_order)}
        for b in batch_order:
            sb = subd.loc[subd["Batch"] == b]
            if not sb.empty:
                xj = sb["charge"].values.astype(float) + rng.uniform(-0.18, 0.18, size=len(sb))
                axd.scatter(xj, sb["TPSA"].values, s=12, alpha=0.5, color=b_to_c[b], edgecolors="none")
        axd.set_xlim(-0.6, charge_max + 0.6)
        axd.set_xticks(np.arange(0, charge_max + 1))
        axd.set_title("(d) Charge vs TPSA (jittered strip)", fontweight="bold")
        axd.set_xlabel("charge", fontweight="bold")
        axd.set_ylabel(r"TPSA ($\AA^2$)", fontweight="bold")
        axd.grid(alpha=0.25, linestyle="--")

    fig.suptitle("Fig 3.6 | Targeting-relevant analysis (charge-driven proxy vs physicochemical balance)",
                  fontweight="bold", y=0.995)
    plt.subplots_adjust(left=0.08, right=0.95, top=0.96, bottom=0.08, wspace=0.18, hspace=0.22)

    _save(fig, out, save_pdf=save_pdf)


# -----------------------------
# Main
# -----------------------------
def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot Fig 3.6 (a–d) for mitochondrial-targeting proxy analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("input_files", nargs="+", type=str, help="Input CSV files (multiple)")
    parser.add_argument("--output-dir", type=str, default=None, help="Output directory (default: script dir)")
    parser.add_argument("--pictos", type=str, default=None, help="pictos.txt path (default: ../data/pictos.txt)")

    parser.add_argument("--filter-imaginary", action="store_true",
                        help="Keep only no-imaginary rows: No/False/0 (recommended)")
    parser.add_argument("--charge-max", type=int, default=6, help="Max charge to display (default: 6)")
    parser.add_argument("--positive-threshold", type=int, default=0, help="Positive if charge > threshold (default: 0)")

    parser.add_argument("--logp-min", type=float, default=1.0, help="Feasible window LogP min (default: 1.0)")
    parser.add_argument("--logp-max", type=float, default=5.0, help="Feasible window LogP max (default: 5.0)")
    parser.add_argument("--tpsa-min", type=float, default=20.0, help="Feasible window TPSA min (default: 20.0)")
    parser.add_argument("--tpsa-max", type=float, default=120.0, help="Feasible window TPSA max (default: 120.0)")

    parser.add_argument("--save-pdf", action="store_true", help="Also save PDF copies")

    parser.add_argument("--no-facets", action="store_true",
                        help="Skip the large facet figure (3x3) to save time/space")
    args = parser.parse_args()

    # Resolve paths
    script_dir = Path(__file__).parent.absolute()
    out_dir = Path(args.output_dir) if args.output_dir else script_dir

    pictos_path = Path(args.pictos) if args.pictos else (script_dir.parent / "data" / "pictos.txt")
    pictos_labels = load_pictos_labels(pictos_path)
    if pictos_labels:
        print(f"Loaded pictos mapping: {len(pictos_labels)} entries from {pictos_path}")
    else:
        print("No pictos mapping loaded (will label batches by filename stem).")

    # Input files
    csv_files: List[Path] = []
    for p in args.input_files:
        fp = Path(p)
        if not fp.exists():
            print(f"Warning: file not found: {fp}")
            continue
        csv_files.append(fp)
    if not csv_files:
        print("Error: no valid input files.")
        return 1

    # Load & concat
    df, batch_order = read_and_concat_csvs(
        csv_files,
        pictos_labels=pictos_labels,
        filter_imaginary=args.filter_imaginary,
    )
    if df.empty:
        print("Error: empty dataframe after loading.")
        return 1

    print(f"Loaded rows: {len(df)}")
    print(f"Batches: {batch_order}")

    # Output paths
    fig_a = out_dir / "fig3_6a_charge_heatmap.png"
    fig_b = out_dir / "fig3_6b_positive_fraction.png"
    fig_c = out_dir / "fig3_6c_logp_tpsa_window.png"
    fig_d = out_dir / "fig3_6d_charge_vs_tpsa_strip.png"
    fig_d2 = out_dir / "fig3_6d2_facets_charge_vs_tpsa.png"
    fig_abcd = out_dir / "fig3_6_abcd_multipanel.png"

    # Plot (a–d)
    plot_36a_charge_heatmap(
        df, batch_order, charge_max=args.charge_max, out=fig_a, save_pdf=args.save_pdf,
        pictos_labels=pictos_labels,
        script_dir=script_dir
    )
    plot_36b_positive_fraction(
        df, batch_order, out=fig_b, save_pdf=args.save_pdf,
        positive_threshold=args.positive_threshold, sort_desc=False, # 强制保持 Batch 排序
        pictos_labels=pictos_labels,
        script_dir=script_dir
    )
    plot_36c_logp_tpsa_window(
        df, batch_order, out=fig_c, save_pdf=args.save_pdf,
        charge_max=args.charge_max,
        logp_window=(args.logp_min, args.logp_max),
        tpsa_window=(args.tpsa_min, args.tpsa_max),
        script_dir=script_dir,
        pictos_labels=pictos_labels
    )
    plot_36d_charge_vs_tpsa_strip(
        df, batch_order, out=fig_d, save_pdf=args.save_pdf,
        charge_max=args.charge_max
    )

    if not args.no_facets:
        plot_36d2_facets_charge_vs_tpsa(
            df, batch_order, out=fig_d2, save_pdf=args.save_pdf,
            charge_max=args.charge_max,
            pictos_labels=pictos_labels,
            script_dir=script_dir
        )

    # Multipanel
    plot_36_abcd_multipanel(
        df, batch_order, out=fig_abcd, save_pdf=args.save_pdf,
        charge_max=args.charge_max,
        positive_threshold=args.positive_threshold,
        logp_window=(args.logp_min, args.logp_max),
        tpsa_window=(args.tpsa_min, args.tpsa_max),
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
