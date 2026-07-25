#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, List, Tuple, Set, Dict

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


# ===================== utils =====================

def unique_preserve(seq: List[int]) -> List[int]:
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


# ===================== Scaffold parsing =====================

def split_scaffold_spec(scaffold_spec: str) -> str:
    """scaffold_spec 可能是 core|frag1|frag2|...，这里只取 core（第一个 | 前）。"""
    return scaffold_spec.split("|", 1)[0].strip()


def core_query_and_attachment_points(core_with_stars: str) -> Tuple[Chem.Mol, Chem.Mol, List[int], List[int]]:
    """
    输入：含 [*] 的 core SMILES
    输出：
      core_full:         含 [*] 的 core mol（用于画 scaffold 图）
      core_backbone:     去掉 [*] 的 core mol（用于在整分子里做匹配）
      dummy_idxs:        core_full 中 dummy([*]) 原子 idx 列表（顺序用于 *1,*2...）
      attach_in_backbone: backbone 中作为连接位点的原子 idx 列表（顺序与 dummy 对应 -> AP1/AP2...）
    """
    core_full = Chem.MolFromSmiles(core_with_stars)
    if core_full is None:
        raise ValueError("无法解析 scaffold core（第一个 | 前的部分）。请确认 SMILES/芳香性写法。")

    dummy_idxs = [a.GetIdx() for a in core_full.GetAtoms() if a.GetAtomicNum() == 0]  # [*] -> atomic num 0
    if not dummy_idxs:
        core_backbone = Chem.Mol(core_full)
        return core_full, core_backbone, [], []

    # dummy 的邻居（这才是母核的连接位点原子，在 core_full 的 idx）
    attach_old = []
    for d in dummy_idxs:
        nbrs = [n.GetIdx() for n in core_full.GetAtomWithIdx(d).GetNeighbors()]
        if len(nbrs) != 1:
            raise ValueError("发现 [*] 不是单连接位点（可能是桥连写法），需要更复杂处理。")
        attach_old.append(nbrs[0])

    # 删除 dummy 得到 backbone
    rw = Chem.RWMol(core_full)
    for d in sorted(dummy_idxs, reverse=True):
        rw.RemoveAtom(d)
    core_backbone = rw.GetMol()
    Chem.SanitizeMol(core_backbone)

    # old(core_full idx) -> new(backbone idx) 映射
    old_to_new = {}
    new_i = 0
    dummy_set = set(dummy_idxs)
    for old_i in range(core_full.GetNumAtoms()):
        if old_i in dummy_set:
            continue
        old_to_new[old_i] = new_i
        new_i += 1

    # 关键：按 dummy 顺序得到 attachment（与 *1,*2... 一致），并去重（正常不应重复）
    attach_in_backbone = unique_preserve([old_to_new[i] for i in attach_old])
    return core_full, core_backbone, dummy_idxs, attach_in_backbone


def get_core_match(mol: Chem.Mol, core_backbone: Chem.Mol) -> Tuple[int, ...]:
    """
    在整分子中找到 core_backbone（去掉[*]后的母核）的一个匹配。
    如有多个匹配，这里默认取第一个。
    """
    matches = mol.GetSubstructMatches(core_backbone)
    if not matches:
        raise ValueError("目标分子中找不到该 scaffold core（去掉[*]后）。")
    return matches[0]


# ===================== Decoration signature =====================

def _collect_branch_atoms(mol: Chem.Mol, start: int, forbidden: Set[int]) -> Set[int]:
    """从 start 出发 BFS，禁止进入 forbidden（母核集合）。"""
    q = [start]
    seen = {start}
    while q:
        x = q.pop()
        for nb in mol.GetAtomWithIdx(x).GetNeighbors():
            j = nb.GetIdx()
            if j in forbidden or j in seen:
                continue
            seen.add(j)
            q.append(j)
    return seen


def decoration_signature(mol: Chem.Mol, core_set: Set[int], attach_atom: int) -> Tuple[str, ...]:
    """
    attach_atom 是母核原子。
    signature = 该母核原子向外（不进入 core）的每个支链片段 canonical SMILES，排序后 tuple。
    """
    atom = mol.GetAtomWithIdx(attach_atom)
    outside = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() not in core_set]
    if not outside:
        return tuple()

    frags = []
    for s in outside:
        branch = sorted(_collect_branch_atoms(mol, s, forbidden=core_set))
        smi = Chem.MolFragmentToSmiles(mol, atomsToUse=branch, canonical=True)
        frags.append(smi)
    return tuple(sorted(frags))


# ===================== Symmetry (AP-based) =====================

def is_involution(perm: Tuple[int, ...]) -> bool:
    """二阶置换：p(p(i)) = i"""
    return all(perm[perm[i]] == i for i in range(len(perm)))


@dataclass
class SymmetryResult:
    axis_symmetric: bool
    center_symmetric: bool
    core_atoms: List[int]
    attachment_atoms: List[int]  # 顺序与 *1,*2... 一致 -> AP1/AP2...
    attachment_signature: Dict[int, Tuple[str, ...]]
    example_axis_perm: Optional[List[int]]
    example_center_perm: Optional[List[int]]
    num_automorphisms_seen: int


def _prepare_core_and_ap(
    smiles: str,
    scaffold_spec: str,
) -> Tuple[Chem.Mol, Set[int], List[int], Dict[int, Tuple[str, ...]], Tuple[int, ...], Chem.Mol, List[int], List[int]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Bad SMILES for molecule.")

    core_smiles = split_scaffold_spec(scaffold_spec)
    core_full, core_backbone, dummy_idxs, attach_in_backbone = core_query_and_attachment_points(core_smiles)
    match = get_core_match(mol, core_backbone)

    core_set = set(match)

    # AP：按 scaffold dummy 顺序映射到整分子 idx（即 AP1/2/3/4 对应 *1/2/3/4）
    attachment_atoms = [match[i] for i in attach_in_backbone]
    attachment_atoms = unique_preserve(attachment_atoms)

    sig = {a: decoration_signature(mol, core_set, a) for a in attachment_atoms}

    return mol, core_set, attachment_atoms, sig, match, core_full, dummy_idxs, attach_in_backbone


def _valid_involution_iter(
    mol: Chem.Mol,
    core_set: Set[int],
    attachment_atoms: List[int],
    sig: Dict[int, Tuple[str, ...]],
    use_chirality: bool,
    max_matches: int,
    min_core_moved_abs: int,
    min_core_moved_frac: float,
):
    autos = mol.GetSubstructMatches(
        mol,
        uniquify=False,
        useChirality=use_chirality,
        maxMatches=max_matches,
    )

    n_core = len(core_set)
    moved_threshold = max(min_core_moved_abs, int(min_core_moved_frac * n_core))

    for perm in autos:
        # 跳过恒等
        if all(i == perm[i] for i in range(len(perm))):
            continue
        # 必须二阶
        if not is_involution(perm):
            continue
        # 必须 core -> core
        if {perm[i] for i in core_set} != core_set:
            continue
        # 过滤局部对称：core 移动太少的不算
        moved_core = sum(1 for i in core_set if perm[i] != i)
        if moved_core < moved_threshold:
            continue
        # 装饰签名必须匹配：sig(AP) == sig(p(AP))
        ok = True
        for a in attachment_atoms:
            if sig[a] != sig[perm[a]]:
                ok = False
                break
        if not ok:
            continue

        yield perm, len(autos)


def check_symmetry_axis_or_center(
    smiles: str,
    scaffold_spec: str,
    use_chirality: bool = False,
    max_matches: int = 100000,
    min_core_moved_abs: int = 6,
    min_core_moved_frac: float = 0.25,
) -> SymmetryResult:
    """
    AP-based 判别：
      - axis:   存在有效 involution，且 AP 上有 fixed 点
      - center: 存在有效 involution，且 AP 上无 fixed 点
    """
    mol, core_set, attachment_atoms, sig, *_ = _prepare_core_and_ap(smiles, scaffold_spec)

    if len(attachment_atoms) == 0:
        return SymmetryResult(
            axis_symmetric=False,
            center_symmetric=False,
            core_atoms=sorted(core_set),
            attachment_atoms=[],
            attachment_signature={},
            example_axis_perm=None,
            example_center_perm=None,
            num_automorphisms_seen=0,
        )

    axis = False
    center = False
    best_axis = None
    best_center = None
    autos_seen = 0

    for perm, autos_seen in _valid_involution_iter(
        mol, core_set, attachment_atoms, sig,
        use_chirality=use_chirality,
        max_matches=max_matches,
        min_core_moved_abs=min_core_moved_abs,
        min_core_moved_frac=min_core_moved_frac,
    ):
        fixed_ap = sum(1 for a in attachment_atoms if perm[a] == a)
        if fixed_ap == 0:
            center = True
            if best_center is None:
                best_center = list(perm)
        else:
            axis = True
            if best_axis is None:
                best_axis = list(perm)
            else:
                old_fixed = sum(1 for a in attachment_atoms if best_axis[a] == a)
                if fixed_ap > old_fixed:
                    best_axis = list(perm)

        if axis and center:
            break

    return SymmetryResult(
        axis_symmetric=axis,
        center_symmetric=center,
        core_atoms=sorted(core_set),
        attachment_atoms=attachment_atoms,
        attachment_signature=sig,
        example_axis_perm=best_axis,
        example_center_perm=best_center,
        num_automorphisms_seen=autos_seen,
    )


# ===================== AP grouping enumeration =====================

def ap_grouping_from_involution(perm: Tuple[int, ...], ap_atoms: List[int]) -> List[Tuple[int, ...]]:
    """
    perm 为 involution。返回 AP 的分组：[(a,), (b,c), ...]
    其中 (a,) 表示 fixed，(b,c) 表示互换。
    """
    ap_set = set(ap_atoms)
    seen = set()
    groups: List[Tuple[int, ...]] = []
    for a in ap_atoms:
        if a in seen:
            continue
        b = perm[a]
        if b == a:
            groups.append((a,))
            seen.add(a)
        else:
            # 理论上 b 也应是 AP（否则不是 AP->AP 的对称，直接按实际输出也能看出来）
            if b in ap_set:
                pair = tuple(sorted((a, b)))
                groups.append(pair)
                seen.add(a)
                seen.add(b)
            else:
                # 极端情况：AP 被映射到非 AP（一般不会，因为你 core->core 且 sig 过滤）
                groups.append((a, b))
                seen.add(a)
                seen.add(b)
    # 规范化排序
    groups = [tuple(g) for g in groups]
    groups.sort(key=lambda x: (len(x), x))
    return groups


def enumerate_ap_groupings(
    smiles: str,
    scaffold_spec: str,
    use_chirality: bool = False,
    max_matches: int = 100000,
    min_core_moved_abs: int = 6,
    min_core_moved_frac: float = 0.25,
) -> Tuple[List[List[Tuple[int, ...]]], List[List[Tuple[int, ...]]]]:
    """
    返回：
      axis_groupings:   所有有效 axis involution 的 AP 分组（去重）
      center_groupings: 所有有效 center involution 的 AP 分组（去重）
    """
    mol, core_set, ap_atoms, sig, *_ = _prepare_core_and_ap(smiles, scaffold_spec)
    if len(ap_atoms) == 0:
        return [], []

    axis_set = set()
    center_set = set()
    axis_groupings = []
    center_groupings = []

    for perm, _ in _valid_involution_iter(
        mol, core_set, ap_atoms, sig,
        use_chirality=use_chirality,
        max_matches=max_matches,
        min_core_moved_abs=min_core_moved_abs,
        min_core_moved_frac=min_core_moved_frac,
    ):
        fixed_ap = sum(1 for a in ap_atoms if perm[a] == a)
        groups = ap_grouping_from_involution(perm, ap_atoms)
        key = tuple(groups)

        if fixed_ap == 0:
            if key not in center_set:
                center_set.add(key)
                center_groupings.append(groups)
        else:
            if key not in axis_set:
                axis_set.add(key)
                axis_groupings.append(groups)

    return axis_groupings, center_groupings


# ===================== Drawing =====================

def draw_png(
    mol: Chem.Mol,
    out_png: str,
    size: Tuple[int, int] = (900, 700),
    atom_labels: Optional[Dict[int, str]] = None,
    highlight_atoms: Optional[List[int]] = None,
):
    mol = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.addAtomIndices = False

    if atom_labels:
        for idx, lab in atom_labels.items():
            if 0 <= idx < mol.GetNumAtoms():
                opts.atomLabels[idx] = lab

    drawer.DrawMolecule(mol, highlightAtoms=(highlight_atoms or []))
    drawer.FinishDrawing()
    with open(out_png, "wb") as f:
        f.write(drawer.GetDrawingText())


def generate_molecule_and_scaffold_figs(
    smiles: str,
    scaffold_spec: str,
    out_prefix: str,
) -> None:
    """
    生成两张图：
      - {prefix}_molecule.png：整分子，core 高亮 + AP1/AP2... 标注（与 *1/*2... 对应）
      - {prefix}_scaffold.png：scaffold(core 含[*])，dummy 标 *1,*2...，邻居标 P1,P2...
    """
    mol, core_set, ap_atoms, _, match, core_full, dummy_idxs, attach_in_backbone = _prepare_core_and_ap(smiles, scaffold_spec)

    core_atoms_in_mol = sorted(core_set)

    # 分子图：APk 标注（k 与 scaffold 的 *k 对应）
    mol_labels = {aidx: f"AP{k}" for k, aidx in enumerate(ap_atoms, start=1)}
    draw_png(
        mol,
        out_png=f"{out_prefix}_molecule.png",
        atom_labels=mol_labels,
        highlight_atoms=core_atoms_in_mol,
    )

    # scaffold 图：标 *1/*2...，并给邻居标 P1/P2...
    core_labels: Dict[int, str] = {}
    for k, didx in enumerate(dummy_idxs, start=1):
        core_labels[didx] = f"*{k}"
        nbr = core_full.GetAtomWithIdx(didx).GetNeighbors()[0].GetIdx()
        core_labels[nbr] = (core_labels.get(nbr, "") + ("|" if core_labels.get(nbr) else "") + f"P{k}")

    draw_png(
        core_full,
        out_png=f"{out_prefix}_scaffold.png",
        atom_labels=core_labels,
        highlight_atoms=dummy_idxs,
    )


# ===================== Formatting outputs =====================

def format_grouping(groups: List[Tuple[int, ...]], ap_label: Dict[int, str]) -> str:
    parts = []
    for g in groups:
        if len(g) == 1:
            a = g[0]
            parts.append(f"{ap_label.get(a, str(a))} (fixed)")
        elif len(g) == 2:
            a, b = g
            parts.append(f"{ap_label.get(a, str(a))} <-> {ap_label.get(b, str(b))}")
        else:
            parts.append(str(tuple(ap_label.get(x, str(x)) for x in g)))
    return " | ".join(parts)


def show_attach_mapping(perm: Optional[List[int]], ap_atoms: List[int], ap_label: Dict[int, str], title: str):
    if not perm:
        return
    fixed_ap = sum(1 for a in ap_atoms if perm[a] == a)
    print(f"\n--- {title}: AP mapping ---")
    print(f"fixed_AP = {fixed_ap} / {len(ap_atoms)}")
    for a in ap_atoms:
        print(f"{ap_label.get(a, str(a))} ({a}) -> {ap_label.get(perm[a], str(perm[a]))} ({perm[a]})")


# ===================== Main =====================

def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="AP-based axis/center symmetry (scaffold+decorations), enumerate AP groupings, and draw PNGs."
    )
    ap.add_argument("--smiles", required=True, help="Molecule SMILES")
    ap.add_argument("--scaffold", required=True, help="Scaffold spec: core|frag1|... (we use only core)")
    ap.add_argument("--prefix", default="out", help="Output prefix for PNGs")
    ap.add_argument("--use-chirality", action="store_true", help="Use chirality when enumerating automorphisms")
    ap.add_argument("--max-matches", type=int, default=100000, help="Max automorphisms to enumerate")
    ap.add_argument("--min-core-moved-abs", type=int, default=6, help="Filter: min number of moved core atoms")
    ap.add_argument("--min-core-moved-frac", type=float, default=0.25, help="Filter: min fraction of moved core atoms")
    ap.add_argument("--max-groupings", type=int, default=50, help="Max groupings to print per kind")
    args = ap.parse_args()

    # 1) symmetry
    res = check_symmetry_axis_or_center(
        args.smiles,
        args.scaffold,
        use_chirality=args.use_chirality,
        max_matches=args.max_matches,
        min_core_moved_abs=args.min_core_moved_abs,
        min_core_moved_frac=args.min_core_moved_frac,
    )

    ap_atoms = res.attachment_atoms
    ap_label = {a: f"AP{k}" for k, a in enumerate(ap_atoms, start=1)}

    print("=== Symmetry (AP-based, topological) ===")
    print(f"axis_symmetric   : {res.axis_symmetric}")
    print(f"center_symmetric : {res.center_symmetric}")
    print(f"core_atoms       : {len(res.core_atoms)} atoms")
    print(f"attachment_atoms : {len(res.attachment_atoms)} (AP1..APn aligns with *1..*n)")
    print("AP index mapping :", ", ".join([f"{ap_label[a]}={a}" for a in ap_atoms]))
    print(f"automorphisms_seen (capped) : {res.num_automorphisms_seen}")

    show_attach_mapping(res.example_axis_perm, ap_atoms, ap_label, "Example AXIS involution")
    show_attach_mapping(res.example_center_perm, ap_atoms, ap_label, "Example CENTER involution")

    # 2) enumerate all AP groupings (dedup)
    axis_groups, center_groups = enumerate_ap_groupings(
        args.smiles,
        args.scaffold,
        use_chirality=args.use_chirality,
        max_matches=args.max_matches,
        min_core_moved_abs=args.min_core_moved_abs,
        min_core_moved_frac=args.min_core_moved_frac,
    )

    print("\n=== AP symmetry groupings (all valid involutions, deduped) ===")

    if center_groups:
        print(f"\nCENTER groupings (fixed_AP=0): {len(center_groups)}")
        for i, g in enumerate(center_groups[: args.max_groupings], start=1):
            print(f"  [{i}] {format_grouping(g, ap_label)}")
        if len(center_groups) > args.max_groupings:
            print(f"  ... ({len(center_groups) - args.max_groupings} more, increase --max-groupings to show)")
    else:
        print("\nCENTER groupings: none found")

    if axis_groups:
        print(f"\nAXIS groupings (fixed_AP>0): {len(axis_groups)}")
        for i, g in enumerate(axis_groups[: args.max_groupings], start=1):
            print(f"  [{i}] {format_grouping(g, ap_label)}")
        if len(axis_groups) > args.max_groupings:
            print(f"  ... ({len(axis_groups) - args.max_groupings} more, increase --max-groupings to show)")
    else:
        print("\nAXIS groupings: none found")

    # 3) draw figures
    generate_molecule_and_scaffold_figs(args.smiles, args.scaffold, out_prefix=args.prefix)
    print("\n=== Saved PNGs ===")
    print(f"{args.prefix}_molecule.png")
    print(f"{args.prefix}_scaffold.png")


if __name__ == "__main__":
    main()


"""
python symmetry_and_drawv2.py \
  --smiles "CCN(CC)c1ccc(-c2ccc(-c3c4nc(c(-c5ccc(N(CC)CC)cc5)c5ccc([nH]5)c(-c5ccc(-c6ccc(N(CC)CC)cc6)cc5)c5nc(c(-c6ccc(N(CC)CC)cc6)c6ccc3[nH]6)C=C5)C=C4)cc2)cc1" \
  --scaffold "[*]c1c2nc(c([*])c3ccc([nH]3)c([*])c3nc(c([*])c4ccc1[nH]4)C=C3)C=C2|c1([*])[n+](C)cc2ccccc2c1|c1(N(C)C)ccc([*])cc1|c1[n+](C)c([*])cc2ccccc12|c1cc(N(C)C)ccc1[*]" \
  --prefix demo

python symmetry_and_drawv2.py \
  --smiles "Cn1c(-c2c3nc(c(-c4n(C)c5ccccc5[n+]4C)c4ccc([nH]4)c(-c4n(C)c5ccccc5[n+]4C)c4nc(c(-c5n(C)c6ccccc6[n+]5C)c5ccc2[nH]5)C=C4)C=C3)[n+](C)c2ccccc21" \
  --scaffold "[*]c1c2nc(c([*])c3ccc([nH]3)c([*])c3nc(c([*])c4ccc1[nH]4)C=C3)C=C2|c1([*])[n+](C)c2ccccc2n1C|c1cc2[n+](C)c([*])n(C)c2cc1|c12n(C)c([*])[n+](C)c1cccc2|C[n+]1c([*])n(C)c2ccccc12" \
  --prefix demo

python symmetry_and_drawv2.py \
  --smiles "CC1(C)c2ccccc2N(c2ccc(-c3c4nc(c(-c5ccc(N6c7ccccc7C(C)(C)c7ccccc76)cc5)c5ccc([nH]5)c(-c5ccc(N6c7ccccc7C(C)(C)c7ccccc76)cc5)c5nc(c(-c6ccc(N7c8ccccc8C(C)(C)c8ccccc87)cc6)c6ccc3[nH]6)C=C5)C=C4)cc2)c2ccccc21 " \
  --scaffold "[*]c1c2nc(c([*])c3ccc([nH]3)c([*])c3nc(c([*])c4ccc1[nH]4)C=C3)C=C2|[*]c1ccc(N2c3ccccc3C(C)(C)c3ccccc32)cc1|c1c2c(ccc1)C(C)(C)c1c(cccc1)N2c1ccc([*])cc1|c1([*])ccc(N2c3c(cccc3)C(C)(C)c3ccccc32)cc1|c1cc([*])ccc1N1c2c(cccc2)C(C)(C)c2ccccc21" \
  --prefix demo 

"""