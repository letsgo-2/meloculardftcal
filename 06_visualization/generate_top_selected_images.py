#!/usr/bin/env python3
"""
从候选CSV中筛选出满足条件的数据并绘图：
1) min(S1-Ti) 比参考CSV某个 index（如 pic3）更小
2) SOCME/δE 比该 index 更大

然后按排序规则取前N（默认10）生成带标注分子图。
"""

import argparse
import os
import re
from pathlib import Path

import pandas as pd

import generate_molecule_images_with_props as gm


def normalize_key(text):
    if text is None:
        return ""
    s = str(text).replace("\ufeff", "").strip().lower()
    for ch in [" ", "-", "_", "(", ")", "[", "]", "{", "}", "/", "\\"]:
        s = s.replace(ch, "")
    return s


def find_column(df, aliases):
    cols = list(df.columns)
    exact = {c: c for c in cols}
    norm_map = {normalize_key(c): c for c in cols}
    for a in aliases:
        if a in exact:
            return exact[a]
    for a in aliases:
        key = normalize_key(a)
        if key in norm_map:
            return norm_map[key]
    return None


def infer_target_index(candidate_path, cand):
    stem = candidate_path.stem
    m = re.match(r"^([A-Za-z]+\d+)", stem)
    if m:
        return m.group(1)
    if "ID" in cand.columns:
        for v in cand["ID"]:
            if pd.notna(v):
                return str(v).split("_")[0].strip()
    if "index" in cand.columns:
        for v in cand["index"]:
            if pd.notna(v):
                return str(v).strip()
    raise ValueError("无法自动推断 target index（请检查候选文件名或ID列）")


def main():
    parser = argparse.ArgumentParser(description="筛选并绘制前N个满足条件的分子图")
    parser.add_argument("candidate_csv", help="候选CSV路径，例如: data/pic3inpdftv22.csv")
    parser.add_argument(
        "--reference-csv",
        required=True,
        help="参考CSV路径（必填，必须包含性质列），例如: data/publishinpdftv22.csv",
    )
    parser.add_argument("--top-n", type=int, default=10, help="取前N个（默认10）")
    parser.add_argument(
        "--tool",
        choices=["rdkit", "openbabel"],
        default="rdkit",
        help="绘图工具（默认rdkit）",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="输出目录（默认: data/image/<候选文件名>_selected）",
    )
    args = parser.parse_args()

    candidate_path = Path(args.candidate_csv)
    reference_path = Path(args.reference_csv)
    if not candidate_path.exists():
        raise FileNotFoundError(f"候选文件不存在: {candidate_path}")
    if not reference_path.exists():
        raise FileNotFoundError(f"参考文件不存在: {reference_path}")

    cand = pd.read_csv(candidate_path)
    ref = pd.read_csv(reference_path)

    smiles_col = find_column(cand, ["SMILES", "smiles"])
    cand_min_col = find_column(cand, ["min(S1-Ti)", "min_positive_S1-Ti", "min_positive_S1_Ti"])
    cand_soc_col = find_column(cand, ["SOCME/δE", "max_SOCME_over_deltaE_S1_Ti_value", "SOCME_over_deltaE"])

    ref_index_col = find_column(ref, ["index", "Index"])
    ref_min_col = find_column(ref, ["min(S1-Ti)", "min_positive_S1-Ti", "min_positive_S1_Ti"])
    ref_soc_col = find_column(ref, ["SOCME/δE", "max_SOCME_over_deltaE_S1_Ti_value", "SOCME_over_deltaE"])

    required = {
        "smiles_col": smiles_col,
        "cand_min_col": cand_min_col,
        "cand_soc_col": cand_soc_col,
        "ref_index_col": ref_index_col,
        "ref_min_col": ref_min_col,
        "ref_soc_col": ref_soc_col,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise ValueError(
            "参考CSV必须包含 index、min(S1-Ti)、SOCME/δE 等性质列；"
            f"当前缺少: {missing}。请使用 data/publishinpdftv22.csv 作为参考文件。"
        )

    target_index = infer_target_index(candidate_path, cand)
    ref_one = ref[ref[ref_index_col].astype(str).str.strip() == target_index]
    if ref_one.empty:
        raise ValueError(f"参考CSV中未找到 index={target_index} 的条目")
    ref_row = ref_one.iloc[0]
    r_min = float(ref_row[ref_min_col])
    r_soc = float(ref_row[ref_soc_col])

    selected = []
    for _, row in cand.iterrows():
        try:
            c_min = float(row[cand_min_col])
            c_soc = float(row[cand_soc_col])
        except Exception:
            continue

        if (c_min < r_min) and (c_soc > r_soc):
            delta_soc = c_soc - r_soc
            delta_min = r_min - c_min
            selected.append(
                {
                    "row": row,
                    "c_min": c_min,
                    "c_soc": c_soc,
                    "delta_soc": delta_soc,
                    "delta_min": delta_min,
                }
            )

    total = len(cand)
    hit_total = len(selected)
    ratio = (hit_total / total) if total else 0.0
    print("=== 筛选统计 ===")
    print(f"参考文件: {reference_path}")
    print(f"目标index: {target_index}")
    print(f"参考阈值: min(S1-Ti)={r_min:.6g}, SOCME/δE={r_soc:.6g}")
    print(f"候选总数: {total}")
    print(f"满足条件数量: {hit_total}")
    print(f"满足条件比例: {ratio:.4%}")

    if not selected:
        print("没有找到满足条件的数据。")
        return

    # 优先 SOCME/δE 提升更大，再看 min(S1-Ti) 降低更大
    selected.sort(key=lambda x: (x["delta_soc"], x["delta_min"]), reverse=True)
    top_n = max(1, args.top_n)
    selected = selected[:top_n]

    if args.output_dir:
        out_dir = Path(args.output_dir)
    else:
        out_dir = candidate_path.parent / "image" / f"{candidate_path.stem}_selected"
    os.makedirs(out_dir, exist_ok=True)

    success = 0
    for i, item in enumerate(selected, start=1):
        row = item["row"]
        smiles = row[smiles_col]
        if pd.isna(smiles) or not str(smiles).strip():
            continue

        rid = row.get("ID", f"{target_index}_{i:04d}")
        filename = f"{rid}.png"
        out_path = out_dir / filename

        ok = gm.generate_molecule_image(str(smiles), str(out_path), tool=args.tool)
        if not ok:
            continue

        props = gm.collect_properties(row)
        props = [("No.", str(i))] + props + [("ΔSOCME/δE", gm.fmt_value(item["delta_soc"]))]
        gm.add_properties_text_to_image(out_path, props)
        success += 1
        print(f"[{success:02d}] {filename}")

    print("\n完成!")
    print(f"筛选命中总数: {hit_total}")
    print(f"输出 top-n: {len(selected)}")
    print(f"成功绘图: {success}")
    print(f"输出目录: {out_dir}")


if __name__ == "__main__":
    main()

