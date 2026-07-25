#!/usr/bin/env python3
"""
从指定的CSV文件读取SMILES，生成分子结构图片并叠加关键性质文本。

叠加字段（按展示标签）:
- HOMO-LUMO(eV)
- S1(eV)
- T1(eV)
- S1-T1(eV)
- min(S1-Ti)
- SOCME(cm -1)
- SOCME/δE

用法:
    # 默认使用RDKit
    python generate_molecule_images_with_props.py data/publishinpdftv22.csv

    # 指定工具
    python generate_molecule_images_with_props.py data/publishinpdftv22.csv --tool openbabel

    # 指定输出基础目录
    python generate_molecule_images_with_props.py data/publishinpdftv22.csv --output data/image

    # 测试模式：只生成指定序号
    python generate_molecule_images_with_props.py data/publishinpdftv22.csv --test 1
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: RDKit未安装，将无法使用rdkit工具")


def check_obabel_command():
    """检查obabel命令行工具是否可用。"""
    try:
        result = subprocess.run(
            ["obabel", "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        return "Open Babel" in result.stdout or "Open Babel" in result.stderr
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
    except Exception:
        return False


def generate_molecule_image_rdkit(smiles, output_path, size=(400, 400)):
    """使用RDKit从SMILES生成分子结构图片。"""
    if not RDKIT_AVAILABLE:
        return False

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"警告: 无法解析SMILES: {smiles}")
            return False
        img = Draw.MolToImage(mol, size=size)
        img.save(output_path)
        return True
    except Exception as e:
        print(f"错误: RDKit生成图片失败 {output_path}: {e}")
        return False


def generate_molecule_image_openbabel(smiles, output_path, size=(400, 400)):
    """使用Open Babel命令行工具从SMILES生成分子结构图片。"""
    try:
        cmd = [
            "obabel",
            "-:" + smiles,
            "-O",
            str(output_path),
            "--gen2d",
            "-xp",
            str(size[0]),
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0 and os.path.exists(output_path):
            return os.path.getsize(output_path) > 1000
        return False
    except FileNotFoundError:
        print("错误: 找不到obabel命令，请确保Open Babel已安装并在PATH中")
        return False
    except Exception as e:
        print(f"错误: Open Babel生成图片失败 {output_path}: {e}")
        return False


def generate_molecule_image(smiles, output_path, tool="rdkit", size=(400, 400)):
    """根据工具选择生成分子结构图片。"""
    if tool == "rdkit":
        return generate_molecule_image_rdkit(smiles, output_path, size)
    if tool == "openbabel":
        if check_obabel_command():
            return generate_molecule_image_openbabel(smiles, output_path, size)
        print("错误: Open Babel命令行工具不可用，请安装: apt-get install openbabel")
        return False
    print(f"错误: 未知的工具: {tool}")
    return False


def get_value_by_aliases(row, aliases):
    """按候选列名顺序读取值。"""
    for col in aliases:
        if col in row.index and pd.notna(row[col]):
            return row[col]
    return None


def normalize_key(text):
    """统一列名格式，减少因空格/大小写/符号导致的匹配失败。"""
    if text is None:
        return ""
    s = str(text).replace("\ufeff", "").strip().lower()
    for ch in [" ", "-", "_", "(", ")", "[", "]", "{", "}", "/", "\\"]:
        s = s.replace(ch, "")
    return s


def get_value_by_aliases_fuzzy(row, aliases):
    """先精确匹配，再按规范化列名模糊匹配。"""
    exact = get_value_by_aliases(row, aliases)
    if exact is not None:
        return exact

    norm_to_real = {}
    for col in row.index:
        norm_to_real[normalize_key(col)] = col

    for alias in aliases:
        key = normalize_key(alias)
        if key in norm_to_real:
            real_col = norm_to_real[key]
            val = row[real_col]
            if pd.notna(val):
                return val
    return None


def fmt_value(v, digits=6):
    """将值格式化为紧凑字符串。"""
    if v is None:
        return "NA"
    if isinstance(v, str):
        return v
    try:
        return f"{float(v):.{digits}g}"
    except Exception:
        return str(v)


def collect_properties(row):
    """提取并组织用于图片展示的性质字段。"""
    field_map = [
        ("HOMO-LUMO(eV)", ["HOMO-LUMO(eV)", "Gap(eV)", "gap_eV"]),
        ("S1(eV)", ["S1(eV)", "S1_energy_eV", "S1_energy"]),
        ("T1(eV)", ["T1(eV)", "T1_energy_eV", "T1_energy"]),
        ("S1-T1(eV)", ["S1-T1(eV)", "S1-T1", "S1_minus_T1"]),
        ("min(S1-Ti)", ["min(S1-Ti)", "min_positive_S1-Ti", "min_positive_S1_Ti"]),
        ("SOCME(cm -1)", ["SOCME(cm -1)", "max_SOCME_S1_Ti_value", "max_SOCME"]),
        ("SOCME/δE", ["SOCME/δE", "max_SOCME_over_deltaE_S1_Ti_value", "SOCME_over_deltaE"]),
    ]

    props = []
    for label, aliases in field_map:
        props.append((label, fmt_value(get_value_by_aliases_fuzzy(row, aliases))))

    # 兜底: 若未直接读到HOMO-LUMO，尝试由HOMO和LUMO计算
    homo_lumo_idx = 0
    if props[homo_lumo_idx][1] == "NA":
        homo = get_value_by_aliases_fuzzy(row, ["HOMO(eV)", "HOMO", "homo"])
        lumo = get_value_by_aliases_fuzzy(row, ["LUMO(eV)", "LUMO", "lumo"])
        try:
            if homo is not None and lumo is not None:
                gap = float(lumo) - float(homo)
                props[homo_lumo_idx] = ("HOMO-LUMO(eV)", fmt_value(gap))
        except Exception:
            pass

    return props


def get_readable_font(size=18):
    """优先加载常见TrueType字体，失败时回退默认字体。"""
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    ]
    for p in candidates:
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, size=size)
            except Exception:
                continue
    return ImageFont.load_default()


def add_properties_text_to_image(image_path, properties, margin=14):
    """在分子图片底部追加文本区域并绘制性质信息。"""
    img = Image.open(image_path).convert("RGB")
    font = get_readable_font(size=18)
    draw = ImageDraw.Draw(img)

    lines = [f"{k}: {v}" for k, v in properties]
    line_heights = []
    max_width = 0
    for line in lines:
        bbox = draw.textbbox((0, 0), line, font=font)
        w = bbox[2] - bbox[0]
        h = bbox[3] - bbox[1]
        max_width = max(max_width, w)
        line_heights.append(h)

    line_gap = 8
    title = "Properties"
    tb = draw.textbbox((0, 0), title, font=font)
    title_h = tb[3] - tb[1]
    text_height = title_h + 8 + sum(line_heights) + line_gap * (len(lines) - 1)
    new_width = max(img.width, max_width + margin * 2)
    new_height = img.height + text_height + margin * 2

    canvas = Image.new("RGB", (new_width, new_height), "white")
    canvas.paste(img, ((new_width - img.width) // 2, 0))
    cdraw = ImageDraw.Draw(canvas)

    # 用浅灰背景区突出文字信息，避免“看不见”
    panel_top = img.height
    cdraw.rectangle(
        [(0, panel_top), (new_width, new_height)],
        fill=(242, 242, 242),
        outline=(180, 180, 180),
        width=1,
    )
    cdraw.line([(0, panel_top), (new_width, panel_top)], fill=(120, 120, 120), width=2)

    y = img.height + margin
    cdraw.text((margin, y), title, fill=(20, 20, 20), font=font)
    y += title_h + 8
    for line in lines:
        cdraw.text((margin, y), line, fill=(20, 20, 20), font=font)
        bbox = cdraw.textbbox((margin, y), line, font=font)
        y += (bbox[3] - bbox[1]) + line_gap

    canvas.save(image_path)


def resolve_output_dir(csv_path, output_base_dir, tool):
    """确定输出目录。"""
    file_base = csv_path.stem
    if output_base_dir is None:
        base_output = csv_path.parent / "image" / file_base
        if tool == "openbabel":
            return Path(str(base_output) + "ob")
        return base_output
    return Path(output_base_dir) / file_base


def process_csv_file(csv_file_path, output_base_dir=None, tool="rdkit", test_index=None, verbose=False):
    """处理CSV并生成带性质标注的图片。"""
    csv_path = Path(csv_file_path)
    if not csv_path.exists():
        print(f"错误: 文件不存在: {csv_file_path}")
        return

    if tool == "rdkit" and not RDKIT_AVAILABLE:
        print("错误: RDKit不可用，请安装: pip install rdkit")
        return
    if tool == "openbabel" and not check_obabel_command():
        print("错误: Open Babel命令行工具不可用，请安装: apt-get install openbabel")
        return

    output_dir = resolve_output_dir(csv_path, output_base_dir, tool)
    os.makedirs(output_dir, exist_ok=True)

    print(f"处理文件: {csv_path.name}")
    print(f"使用工具: {tool}")
    print(f"输出目录: {output_dir}")

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"错误: 读取CSV失败: {e}")
        return

    if "SMILES" not in df.columns:
        print(f"错误: {csv_path.name} 中没有SMILES列")
        return

    if test_index is not None:
        if test_index < 1 or test_index > len(df):
            print(f"错误: 序号 {test_index} 超出范围（1-{len(df)}）")
            return
        target_indices = [test_index - 1]
    else:
        target_indices = list(df.index)

    total_images = 0
    failed_images = 0
    file_base = csv_path.stem

    for idx in target_indices:
        row = df.iloc[idx]
        smiles = row.get("SMILES")
        if pd.isna(smiles) or not smiles:
            failed_images += 1
            continue

        image_filename = f"{file_base}_{idx + 1:04d}.png"
        output_path = output_dir / image_filename

        ok = generate_molecule_image(smiles, str(output_path), tool=tool)
        if not ok:
            failed_images += 1
            continue

        props = collect_properties(row)
        try:
            add_properties_text_to_image(output_path, props)
            if verbose:
                print(f"[{idx + 1}] {image_filename}")
                print("  " + " | ".join([f"{k}={v}" for k, v in props]))
            total_images += 1
        except Exception as e:
            print(f"警告: 图片已生成，但叠加文字失败 {output_path}: {e}")
            failed_images += 1

    print("\n完成!")
    print(f"成功生成 {total_images} 张带性质标注的图片")
    if failed_images > 0:
        print(f"失败 {failed_images} 张图片")
    print(f"图片保存在: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="从CSV读取SMILES并生成带性质标注的分子结构图")
    parser.add_argument("csv_file", help="CSV文件路径，例如: data/publishinpdftv22.csv")
    parser.add_argument("--output", "-o", default=None, help="输出基础目录 (默认: data/image)")
    parser.add_argument(
        "--tool",
        "-t",
        choices=["rdkit", "openbabel"],
        default="rdkit",
        help="使用的工具: rdkit 或 openbabel (默认: rdkit)",
    )
    parser.add_argument(
        "--test",
        type=int,
        metavar="INDEX",
        help="测试模式: 只生成指定序号（从1开始），例如: --test 1",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="打印每张图片实际叠加的属性值，便于核对是否正确读取CSV",
    )
    args = parser.parse_args()

    process_csv_file(
        csv_file_path=args.csv_file,
        output_base_dir=args.output,
        tool=args.tool,
        test_index=args.test,
        verbose=args.verbose,
    )
