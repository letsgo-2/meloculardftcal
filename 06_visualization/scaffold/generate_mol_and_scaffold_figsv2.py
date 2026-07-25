#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""主要用来处理pic4，因为不符合规则"""
from __future__ import annotations

from typing import Tuple, List, Dict, Optional, Set
import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AdjustQueryParameters, AdjustQueryProperties


# =========================
# Utilities
# =========================

def _strip_numbers_for_drawing(m: Chem.Mol) -> None:
    """防止 RDKit 把 map number / isotope / tag 画出来（例如 C:12 或 13C）"""
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
        a.SetIsotope(0)


# =========================
# Scaffold parsing & match
# =========================

def split_scaffold_spec(scaffold_spec: str) -> str:
    """scaffold_spec: core|frag1|frag2|... 这里只取 core（第一个 | 前）"""
    return scaffold_spec.split("|", 1)[0].strip()


def core_query_and_attachment_points(
    core_with_stars: str
) -> Tuple[Chem.Mol, Chem.Mol, List[int], List[int]]:
    """
    返回：
      core_full:     含 [*] 的 core mol（用于画 scaffold 图）
      core_backbone: 去掉 [*] 的 core mol（用于在整分子里做匹配）
      dummy_idxs:    core_full 中 dummy([*]) 原子 idx 列表（保序）
      attach_in_backbone: backbone 中作为连接位点的原子 idx 列表（与 dummy 相邻；保序，与 dummy 对应）
    """
    core_full = Chem.MolFromSmiles(core_with_stars)
    if core_full is None:
        raise ValueError(f"无法解析 scaffold core：{core_with_stars}")

    dummy_idxs = [a.GetIdx() for a in core_full.GetAtoms() if a.GetAtomicNum() == 0]
    if not dummy_idxs:
        core_backbone = Chem.Mol(core_full)
        return core_full, core_backbone, [], []

    # dummy 的邻居（在 core_full 里的 idx）——顺序与 dummy 一致
    attach_old: List[int] = []
    for d in dummy_idxs:
        nbrs = [n.GetIdx() for n in core_full.GetAtomWithIdx(d).GetNeighbors()]
        if len(nbrs) != 1:
            raise ValueError("发现 [*] 不是单连接位点（可能是桥连写法），需要更复杂处理。")
        attach_old.append(nbrs[0])

    # 删除 dummy 得 backbone
    rw = Chem.RWMol(core_full)
    for d in sorted(dummy_idxs, reverse=True):
        rw.RemoveAtom(d)
    core_backbone = rw.GetMol()
    Chem.SanitizeMol(core_backbone)

    # old->new 映射：core_full idx -> backbone idx
    old_to_new: Dict[int, int] = {}
    new_i = 0
    dummy_set = set(dummy_idxs)
    for old_i in range(core_full.GetNumAtoms()):
        if old_i in dummy_set:
            continue
        old_to_new[old_i] = new_i
        new_i += 1

    attach_in_backbone: List[int] = [old_to_new[old_i] for old_i in attach_old]
    return core_full, core_backbone, dummy_idxs, attach_in_backbone


def make_generic_element_query(core_backbone: Chem.Mol) -> Chem.Mol:
    """
    把 core_backbone 转成“鲁棒 query”：
      - 原子：只约束元素 [#Z]（不约束电荷/芳香/氢数）
      - 键：generic bond（不约束单/双/芳香）
    用于解决共振/芳香表示差异导致的匹配失败。
    """
    params = AdjustQueryParameters()
    params.makeBondsGeneric = True
    # bond 会变成 query bond；atom 仍是普通 Atom，需要我们换成 QueryAtom
    q = AdjustQueryProperties(core_backbone, params)

    rw = Chem.RWMol(q)
    for i, a in enumerate(q.GetAtoms()):
        z = a.GetAtomicNum()
        if z == 0:
            qatom = Chem.AtomFromSmarts("[*]")
        else:
            qatom = Chem.AtomFromSmarts(f"[#{z}]")
        rw.ReplaceAtom(i, qatom, preserveProps=False)

    q2 = rw.GetMol()
    return q2


def get_core_match(mol: Chem.Mol, core_backbone: Chem.Mol) -> Tuple[int, ...]:
    """
    优先用“鲁棒 query”匹配（忽略 bond order/charge/aromatic），失败再尝试精确匹配。
    """
    query = make_generic_element_query(core_backbone)
    matches = mol.GetSubstructMatches(query, useChirality=False, maxMatches=20)
    if matches:
        return matches[0]

    # fallback：精确匹配（便于你发现是不是 scaffold 写错了）
    matches2 = mol.GetSubstructMatches(core_backbone, useChirality=False, maxMatches=20)
    if matches2:
        return matches2[0]

    raise ValueError("目标分子中找不到该 scaffold core（去掉[*]后）。")


# =========================
# Drawing
# =========================

def draw_png(
    mol: Chem.Mol,
    out_png: Path,
    size: Tuple[int, int] = (900, 700),
    atom_labels: Optional[Dict[int, str]] = None,
    highlight_atoms: Optional[List[int]] = None,
    highlight_bonds: Optional[List[int]] = None,
):
    """
    约束：
    - 不显示编号（atom index / bond index）
    - 不写 C（RDKit 默认骨架线式不写 C；我们也不对非 dummy 写 label）
    - 只对指定原子写 label（用于连接点 "*")
    """
    mol = Chem.Mol(mol)
    _strip_numbers_for_drawing(mol)

    try:
        rdDepictor.SetPreferCoordGen(True)
    except Exception:
        pass
    rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(int(size[0]), int(size[1]))
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    opts.addBondIndices = False

    # 关掉可能导致 tag/map 画出来的选项（不同 RDKit 版本字段名略不同）
    for attr in ("includeAtomTags", "includeAtomMapNumbers", "addAtomMapNumbers", "explicitMethyl"):
        if hasattr(opts, attr):
            try:
                setattr(opts, attr, False)
            except Exception:
                pass

    if atom_labels:
        for idx, lab in atom_labels.items():
            if 0 <= idx < mol.GetNumAtoms():
                opts.atomLabels[idx] = str(lab)

    if highlight_atoms is None:
        highlight_atoms = []
    if highlight_bonds is None:
        highlight_bonds = []

    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightBonds=highlight_bonds)
    drawer.FinishDrawing()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_png.write_bytes(drawer.GetDrawingText())


# =========================
# Decorations extraction
# =========================

def extract_and_draw_decorations(
    mol: Chem.Mol,
    core_atoms_in_mol: List[int],
    attach_atoms_in_mol_inorder: List[int],
    out_prefix: Path,
    img_size: Tuple[int, int],
):
    """
    切断“骨架(core)–装饰(R)”连接键，提取所有装饰片段并画出来。
    - 每个装饰片段单独成图：{prefix}_dec_1.png, {prefix}_dec_2.png ...
    - 装饰里保留 "*" 并高亮（切断键后 dummy 原子）
    - 装饰里不显示编号、不写 C（只给 dummy 写 "*"，且清掉 map/isotope）
    """
    core_set: Set[int] = set(core_atoms_in_mol)

    mol2 = Chem.Mol(mol)
    _strip_numbers_for_drawing(mol2)

    bond_indices: List[int] = []
    dummy_labels: List[Tuple[int, int]] = []

    for aidx in attach_atoms_in_mol_inorder:
        if not (0 <= aidx < mol2.GetNumAtoms()):
            continue
        atom = mol2.GetAtomWithIdx(aidx)
        outside_nbrs = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() not in core_set]
        if not outside_nbrs:
            continue

        for nb in outside_nbrs:
            b = mol2.GetBondBetweenAtoms(aidx, nb)
            if b is None:
                continue
            bond_indices.append(b.GetIdx())
            dummy_labels.append((0, 0))  # 不给 isotope 编号，避免任何数字被画出来

    if not bond_indices:
        print("No decoration bonds found (core may cover whole molecule or no substituents).")
        return

    fragmol = Chem.FragmentOnBonds(
        mol2,
        bondIndices=bond_indices,
        addDummies=True,
        dummyLabels=dummy_labels,
    )
    _strip_numbers_for_drawing(fragmol)

    # atom id 分组（用于判断是否包含 core）
    frag_atom_tuples = Chem.GetMolFrags(fragmol, asMols=False, sanitizeFrags=False)
    frag_mols = Chem.GetMolFrags(fragmol, asMols=True, sanitizeFrags=False)

    deco_frags: List[Chem.Mol] = []
    for atom_ids, fm in zip(frag_atom_tuples, frag_mols):
        if set(atom_ids) & core_set:
            continue
        try:
            Chem.SanitizeMol(fm)
        except Exception:
            pass
        _strip_numbers_for_drawing(fm)
        deco_frags.append(fm)

    if not deco_frags:
        print("No decoration fragments extracted.")
        return

    saved: List[Path] = []
    combined: Optional[Chem.Mol] = None

    for i, fm in enumerate(deco_frags, start=1):
        dummy_atoms = [a.GetIdx() for a in fm.GetAtoms() if a.GetAtomicNum() == 0]
        labels = {idx: "*" for idx in dummy_atoms}

        out = out_prefix.parent.joinpath(f"{out_prefix.name}_dec_{i}.png")
        draw_png(
            fm,
            out_png=out,
            size=img_size,
            atom_labels=labels,
            highlight_atoms=dummy_atoms,
        )
        saved.append(out)
        combined = fm if combined is None else Chem.CombineMols(combined, fm)

    if combined is not None:
        _strip_numbers_for_drawing(combined)
        dummy_atoms = [a.GetIdx() for a in combined.GetAtoms() if a.GetAtomicNum() == 0]
        labels = {idx: "*" for idx in dummy_atoms}

        out_all = out_prefix.parent.joinpath(f"{out_prefix.name}_decorations_all.png")
        draw_png(
            combined,
            out_png=out_all,
            size=img_size,
            atom_labels=labels,
            highlight_atoms=dummy_atoms,
        )
        saved.append(out_all)

    print("Saved decoration figures:")
    for p in saved:
        print(f"  {p}")


# =========================
# Main pipeline
# =========================

def generate_molecule_and_scaffold_figs(
    smiles: str,
    scaffold_spec: str,
    out_prefix: Path,
    img_size: Tuple[int, int] = (900, 700),
):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Bad SMILES for molecule.")
    _strip_numbers_for_drawing(mol)

    core_smiles = split_scaffold_spec(scaffold_spec)
    core_full, core_backbone, dummy_idxs, attach_in_backbone = core_query_and_attachment_points(core_smiles)

    match = get_core_match(mol, core_backbone)
    core_atoms_in_mol_sorted = sorted(set(match))

    # 连接位点（整分子里对应的母核原子）——保序（对应 core 的 [*] 顺序）
    attach_atoms_in_mol_inorder = [match[i] for i in attach_in_backbone]
    attach_atoms_unique = sorted(set(attach_atoms_in_mol_inorder))

    # 1) 整分子图：只标母核连接点 "*" 并高亮（不标红母核）
    mol_labels = {aidx: "*" for aidx in attach_atoms_unique}
    draw_png(
        mol,
        out_png=out_prefix.parent.joinpath(f"{out_prefix.name}_molecule.png"),
        size=img_size,
        atom_labels=mol_labels,
        highlight_atoms=attach_atoms_unique,
    )

    # 2) scaffold(core) 图：只标 dummy [*] 为 "*" 并高亮
    _strip_numbers_for_drawing(core_full)
    core_labels = {didx: "*" for didx in dummy_idxs}
    draw_png(
        core_full,
        out_png=out_prefix.parent.joinpath(f"{out_prefix.name}_scaffold.png"),
        size=img_size,
        atom_labels=core_labels,
        highlight_atoms=dummy_idxs,
    )

    # 3) 装饰片段
    extract_and_draw_decorations(
        mol=mol,
        core_atoms_in_mol=core_atoms_in_mol_sorted,
        attach_atoms_in_mol_inorder=attach_atoms_in_mol_inorder,
        out_prefix=out_prefix,
        img_size=img_size,
    )

    print("Saved base figures:")
    print(f"  {out_prefix.parent.joinpath(out_prefix.name + '_molecule.png')}")
    print(f"  {out_prefix.parent.joinpath(out_prefix.name + '_scaffold.png')}")


def process_csv(
    csv_file: Path,
    smiles_col: str = "SMILES",
    scaffold_col: str = "Scaffold",
    output_dir: Path = Path("."),
    row_1based: int = 1,
    out_name: Optional[str] = None,
    img_size: Tuple[int, int] = (900, 700),
):
    df = pd.read_csv(csv_file)
    if df.empty:
        raise ValueError("CSV 是空的。")

    if row_1based <= 0:
        row_1based = 1
    row0 = row_1based - 1
    if row0 >= len(df):
        raise ValueError(f"CSV 只有 {len(df)} 行，无法选择第 {row_1based} 行。")

    if smiles_col not in df.columns:
        raise ValueError(f"CSV 缺少列：{smiles_col}")
    if scaffold_col not in df.columns:
        raise ValueError(f"CSV 缺少列：{scaffold_col}")

    smiles = str(df.iloc[row0][smiles_col]).strip()
    scaffold = str(df.iloc[row0][scaffold_col]).strip()

    base = out_name.strip() if (out_name and out_name.strip()) else csv_file.stem
    out_prefix = output_dir.joinpath(base)

    print(f"Processing CSV: {csv_file}")
    print(f"Row: {row_1based}/{len(df)}")
    generate_molecule_and_scaffold_figs(smiles, scaffold, out_prefix, img_size=img_size)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="整分子图（仅标母核连接点*）+ scaffold图（仅标[*]点*）+ 装饰片段图（每片段保留*并高亮；不写C不编号；core 匹配更鲁棒）"
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles", "-s", help="输入的 SMILES 字符串")
    group.add_argument("--csv", help="输入 CSV 文件路径")

    parser.add_argument("--scaffold", "-c", help="scaffold 规范字符串（仅在 --smiles 模式下需要）")
    parser.add_argument("--smiles-col", default="SMILES", help="CSV 中 SMILES 列名（默认: SMILES）")
    parser.add_argument("--scaffold-col", default="Scaffold", help="CSV 中 Scaffold 列名（默认: Scaffold）")
    parser.add_argument("--row", "-r", type=int, default=1, help="选择 CSV 的第几行（1 表示第一行；默认: 1）")

    parser.add_argument("--output-dir", default=".", help="输出目录（默认: 当前目录）")
    parser.add_argument("--output-name", default=None, help="输出前缀名字（不含扩展名）。CSV 模式默认用 CSV 文件名。")

    parser.add_argument("--w", type=int, default=900, help="图片宽度（默认: 900）")
    parser.add_argument("--h", type=int, default=700, help="图片高度（默认: 700）")

    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    img_size = (int(args.w), int(args.h))

    try:
        if args.csv:
            process_csv(
                csv_file=Path(args.csv),
                smiles_col=args.smiles_col,
                scaffold_col=args.scaffold_col,
                output_dir=output_dir,
                row_1based=int(args.row),
                out_name=args.output_name,
                img_size=img_size,
            )
        else:
            if not args.scaffold:
                raise ValueError("使用 --smiles 时必须提供 --scaffold 参数。")
            out_name = args.output_name.strip() if args.output_name else "out"
            out_prefix = output_dir.joinpath(out_name)
            generate_molecule_and_scaffold_figs(
                smiles=args.smiles.strip(),
                scaffold_spec=args.scaffold.strip(),
                out_prefix=out_prefix,
                img_size=img_size,
            )
    except Exception as e:
        print(f"错误: {e}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
