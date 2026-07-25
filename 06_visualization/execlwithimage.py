#!/usr/bin/env python3
"""
从 CSV 读取每一行的 SMILES，保持所有原始数据不变，
只在最后一列插入用 obabel 生成的分子结构图片，输出 Excel (.xlsx)。

功能说明：
    - 读取 CSV 文件，保持所有原始列不变
    - 从 SMILES 列读取分子结构字符串
    - 使用 obabel 为每个分子生成 2D 结构图片（PNG）
    - 将图片插入到 Excel 文件的最后一列
    - 自动调整行高和列宽以适应图片

示例用法：
    # 基本用法（使用默认设置）
    python execlwithimage.py \
        --csv data/pic16.csv \
        --output data/pic16_obabel_img.xlsx

    # 指定 SMILES 列名（如果列名不是 SMILES）
    python execlwithimage.py \
        --csv data/pic16.csv \
        --output data/pic16_obabel_img.xlsx \
        --smiles-col SMILES

    # 自定义图片尺寸
    python execlwithimage.py \
        --csv data/pic16.csv \
        --output data/pic16_obabel_img.xlsx \
        --img-width 300 \
        --img-height 300

    # 指定图片输出目录（默认会在 CSV 同级创建 <stem>_img_ob 目录）
    python execlwithimage.py \
        --csv data/pic16.csv \
        --output data/pic16_obabel_img.xlsx \
        --img-dir data/image/pic16_ob

    # 自定义图片列名
    python execlwithimage.py \
        --csv data/pic16.csv \
        --output data/pic16_obabel_img.xlsx \
        --img-col-name "分子结构"

    # 完整示例（所有参数）
    python execlwithimage.py \
        --csv data/pic16.csv \
        --output data/pic16_obabel_img.xlsx \
        --smiles-col SMILES \
        --img-col-name Structure \
        --img-dir data/image/pic16_ob \
        --img-width 250 \
        --img-height 250

依赖要求：
    - pandas: pip install pandas
    - openpyxl: pip install openpyxl
    - obabel: apt-get install openbabel (系统命令)
"""

import argparse
import subprocess
import os
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.utils import get_column_letter
from openpyxl.drawing.image import Image as XLImage


def check_obabel():
    """检查 obabel 命令是否可用"""
    try:
        result = subprocess.run(
            ["obabel", "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        # obabel --version 会输出版本信息到stdout，但返回非零退出码
        # 所以检查输出中是否包含"Open Babel"来判断命令是否可用
        if 'Open Babel' in result.stdout or 'Open Babel' in result.stderr:
            return True
        return False
    except Exception:
        return False


def smiles_to_png_with_obabel(smiles: str, png_path: Path, size=(250, 250)) -> bool:
    """
    使用 obabel 从 SMILES 生成 PNG 图片（2D 构型）

    Args:
        smiles: SMILES 字符串
        png_path: 输出 PNG 路径
        size: (width, height) 像素，注意：obabel 只使用宽度，高度会按比例自动调整
    """
    png_path = Path(png_path)
    # 构造命令
    # 注意：obabel 的 -xp 参数只设置图片宽度，高度会根据分子结构自动按比例调整
    # cmd = [
    #     "obabel",
    #     f"-:{smiles}",          # 从命令行参数读 SMILES
    #     "-O", str(png_path),    # 输出 PNG 文件
    #     "--gen2d",              # 生成 2D 坐标
    #     "-xp", str(size[0]),    # 图片宽度（像素），高度会自动按比例调整
    # ]
    cmd = [
        "obabel",
        f"-:{smiles}",          # 从命令行参数读 SMILES
        "-O", str(png_path),    # 输出 PNG 文件
        "--gen2d",              # 生成 2D 坐标
        "-xp", str(1000),    # 图片宽度（像素），高度会自动按比例调整
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"[obabel 错误] SMILES: {smiles}")
        print(result.stderr.strip())
        return False

    if not png_path.exists() or os.path.getsize(png_path) == 0:
        print(f"[警告] obabel 生成的图片为空: {png_path}")
        return False

    return True


def csv_to_excel_with_obabel_images(
    csv_path,
    xlsx_path,
    smiles_col="SMILES",
    img_col_name="Structure",
    img_dir=None,
    img_size=(200, 200),
):
    csv_path = Path(csv_path)
    xlsx_path = Path(xlsx_path)

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV 不存在: {csv_path}")

    if not check_obabel():
        raise RuntimeError("找不到 obabel 命令，请先安装并加入 PATH")

    # 读取 CSV
    df = pd.read_csv(csv_path)

    if smiles_col not in df.columns:
        raise ValueError(f"CSV 中找不到 SMILES 列: {smiles_col}")

    # 图片输出目录
    if img_dir is None:
        # 默认：和 CSV 同级建一个 子目录，例如 data/pic16_img_ob
        img_dir = csv_path.parent / f"{csv_path.stem}_img_ob"
    img_dir = Path(img_dir)
    img_dir.mkdir(parents=True, exist_ok=True)
    print(f"图片输出目录: {img_dir}")

    # 创建工作簿
    wb = Workbook()
    ws = wb.active
    ws.title = csv_path.stem

    # 表头：原有列 + 最后一列图片列
    headers = list(df.columns) + [img_col_name]
    ws.append(headers)

    img_col_idx = len(headers)             # 图片在最后一列
    img_col_letter = get_column_letter(img_col_idx)

    # 逐行写入
    for i, (_, row) in enumerate(df.iterrows(), start=2):  # Excel 行号从 2 开始
        # 原始数据一字不动写进去
        row_values = list(row.values) + [""]
        ws.append(row_values)

        smiles_val = row[smiles_col]
        if pd.isna(smiles_val):
            continue
        smiles = str(smiles_val).strip()
        if not smiles:
            continue

        # 每行生成一个 PNG：mol_0001.png 这种
        img_name = f"mol_{i-1:04d}.png"    # i-1 对应原始数据行号
        img_path = img_dir / img_name

        ok = smiles_to_png_with_obabel(smiles, img_path, size=img_size)
        if not ok:
            print(f"[跳过] 第 {i-1} 行图片生成失败")
            continue

        try:
            xl_img = XLImage(str(img_path))
            # 设置图片显示大小（Excel 里）
            xl_img.width, xl_img.height = img_size

            # 调整行高，避免图片被压扁
            ws.row_dimensions[i].height = img_size[1] * 0.75

            # 插入到最后一列
            cell_ref = f"{img_col_letter}{i}"
            ws.add_image(xl_img, cell_ref)
        except Exception as e:
            print(f"[插入失败] 行 {i-1}, 图片 {img_path}: {e}")

    # 简单调整列宽
    for col_idx in range(1, img_col_idx + 1):
        col_letter = get_column_letter(col_idx)
        ws.column_dimensions[col_letter].width = 18

    wb.save(xlsx_path)
    print(f"完成！已生成 Excel: {xlsx_path}")


def main():
    parser = argparse.ArgumentParser(
        description="从 CSV 读取 SMILES，用 obabel 生成结构图，插入到 Excel 最后一列"
    )
    parser.add_argument("--csv", required=True, help="输入 CSV 路径，如 data/pic16.csv")
    parser.add_argument(
        "--output",
        required=True,
        help="输出 Excel 路径，如 data/pic16_obabel_img.xlsx",
    )
    parser.add_argument(
        "--smiles-col",
        default="SMILES",
        help="CSV 中 SMILES 列名（默认: SMILES）",
    )
    parser.add_argument(
        "--img-col-name",
        default="Structure",
        help="图片列列名（默认: Structure）",
    )
    parser.add_argument(
        "--img-dir",
        default=None,
        help="图片输出目录（默认: 和 CSV 同级的 <stem>_img_ob）",
    )
    parser.add_argument(
        "--img-width",
        type=int,
        default=200,
        help="图片宽度（像素，默认 200）。注意：obabel 只使用宽度，高度会根据分子结构自动按比例调整",
    )
    parser.add_argument(
        "--img-height",
        type=int,
        default=200,
        help="图片高度（像素，默认 200）。注意：此参数仅用于 Excel 中设置行高，实际图片高度由 obabel 自动按比例调整",
    )

    args = parser.parse_args()
    csv_to_excel_with_obabel_images(
        csv_path=args.csv,
        xlsx_path=args.output,
        smiles_col=args.smiles_col,
        img_col_name=args.img_col_name,
        img_dir=args.img_dir,
        img_size=(args.img_width, args.img_height),
    )


if __name__ == "__main__":
    main()
