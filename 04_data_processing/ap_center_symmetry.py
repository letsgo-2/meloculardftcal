# ap_center_symmetry.py
# -*- coding: utf-8 -*-
"""
AP-based center symmetry (your definition) for nAP=4.

Definition (your rule):
- For 4 attachment points (AP1..AP4), center_symmetric is True iff
  all three perfect pairings can make "within-pair decorations equal".
  This is equivalent to: signature(AP1) == signature(AP2) == signature(AP3) == signature(AP4).

AP order is aligned with the order of [*] in the scaffold core (the first part before '|').
"""

from __future__ import annotations
from typing import List, Tuple, Set, Optional
from rdkit import Chem
from rdkit import RDLogger

# 抑制 RDKit 的警告信息（如化合价超限等）
RDLogger.DisableLog('rdApp.*')


def _split_scaffold_core(scaffold_spec: str) -> str:
    """scaffold_spec may be 'core|frag1|frag2|...'; we only use core."""
    return scaffold_spec.split("|", 1)[0].strip()


def _unique_preserve(seq: List[int]) -> List[int]:
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _core_backbone_and_attach_order(core_with_stars: str) -> Tuple[Chem.Mol, List[int]]:
    """
    Parse core with [*], remove [*] to get backbone.
    Return:
      - backbone mol
      - attach indices in backbone order aligned to [*] order (AP1..)
    """
    core_full = Chem.MolFromSmiles(core_with_stars)
    if core_full is None:
        raise ValueError("Cannot parse scaffold core SMILES.")

    dummy_idxs = [a.GetIdx() for a in core_full.GetAtoms() if a.GetAtomicNum() == 0]  # [*]
    if not dummy_idxs:
        # no [*]
        return Chem.Mol(core_full), []

    # For each dummy in order, record its neighbor (attachment atom in core_full)
    attach_old = []
    for d in dummy_idxs:
        nbrs = [n.GetIdx() for n in core_full.GetAtomWithIdx(d).GetNeighbors()]
        if len(nbrs) != 1:
            raise ValueError("[*] is not singly-connected; bridged patterns not supported.")
        attach_old.append(nbrs[0])

    # Remove dummy atoms -> backbone
    rw = Chem.RWMol(core_full)
    for d in sorted(dummy_idxs, reverse=True):
        rw.RemoveAtom(d)
    backbone = rw.GetMol()
    Chem.SanitizeMol(backbone)

    # Map old index (core_full) -> new index (backbone)
    dummy_set = set(dummy_idxs)
    old_to_new = {}
    new_i = 0
    for old_i in range(core_full.GetNumAtoms()):
        if old_i in dummy_set:
            continue
        old_to_new[old_i] = new_i
        new_i += 1

    attach_in_backbone = _unique_preserve([old_to_new[i] for i in attach_old])
    return backbone, attach_in_backbone


def _get_core_match(mol: Chem.Mol, backbone: Chem.Mol) -> Tuple[int, ...]:
    matches = mol.GetSubstructMatches(backbone)
    if not matches:
        raise ValueError("Backbone (scaffold core) not found in molecule.")
    return matches[0]


def _collect_branch_atoms(mol: Chem.Mol, start: int, forbidden: Set[int]) -> Set[int]:
    """BFS from start without entering forbidden(core atoms)."""
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


def _decoration_signature(mol: Chem.Mol, core_set: Set[int], ap_atom: int) -> Tuple[str, ...]:
    """
    signature(AP) = canonical SMILES of each branch leaving core at AP (excluding core), sorted as tuple.
    """
    atom = mol.GetAtomWithIdx(ap_atom)
    outside = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() not in core_set]
    if not outside:
        return tuple()

    frags = []
    for s in outside:
        branch = sorted(_collect_branch_atoms(mol, s, forbidden=core_set))
        smi = Chem.MolFragmentToSmiles(mol, atomsToUse=branch, canonical=True)
        frags.append(smi)
    return tuple(sorted(frags))


def is_center_symmetric_ap4(
    smiles: str,
    scaffold_spec: str,
    *,
    strict: bool = False,
) -> bool:
    """
    Return True/False only.

    True iff:
      - exactly 4 APs are defined by scaffold core ([*] count == 4)
      - all four AP decoration signatures are identical

    If strict=False (default), any parsing/matching failure returns False.
    If strict=True, raises ValueError on parsing/matching problems.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Bad molecule SMILES.")

        core = _split_scaffold_core(scaffold_spec)
        backbone, attach_in_backbone = _core_backbone_and_attach_order(core)

        # must be 4 APs by your definition
        if len(attach_in_backbone) != 4:
            return False

        match = _get_core_match(mol, backbone)
        core_set = set(match)

        # AP atoms in molecule, ordered to align with *1..*4
        ap_atoms = _unique_preserve([match[i] for i in attach_in_backbone])
        if len(ap_atoms) != 4:
            return False

        sigs = [_decoration_signature(mol, core_set, a) for a in ap_atoms]
        return len(set(sigs)) == 1  # all equal => 3 pairings all valid => center symmetric

    except Exception:
        if strict:
            raise
        return False


__all__ = ["is_center_symmetric_ap4"]


if __name__ == "__main__":
    # Optional quick CLI test:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--smiles", required=True)
    p.add_argument("--scaffold", required=True)
    p.add_argument("--strict", action="store_true")
    args = p.parse_args()
    print(is_center_symmetric_ap4(args.smiles, args.scaffold, strict=args.strict))
