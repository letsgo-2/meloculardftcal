#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
addcenter.py

Read a CSV containing columns: Scaffold, SMILES
Append/overwrite the last column "SameS" with "true"/"false" computed by:
    from ap_center_symmetry import is_center_symmetric_ap4

Usage:
    # Update input file in-place (no -o)
    python addcenter.py input.csv

    # Write to a new output file
    python addcenter.py input.csv -o output.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

from ap_center_symmetry import is_center_symmetric_ap4


def detect_dialect(sample: str) -> csv.Dialect:
    try:
        return csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", "|"])
    except Exception:
        return csv.excel  # fallback (comma)


def process_csv(in_path: Path, out_path: Path) -> None:
    # Read a small sample to infer delimiter
    with in_path.open("r", encoding="utf-8-sig", newline="") as f:
        sample = f.read(4096)
        dialect = detect_dialect(sample)

    with in_path.open("r", encoding="utf-8-sig", newline="") as fin, out_path.open(
        "w", encoding="utf-8", newline=""
    ) as fout:
        reader = csv.DictReader(fin, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError("Empty CSV or missing header row.")

        fieldnames = list(reader.fieldnames)
        if "Scaffold" not in fieldnames or "SMILES" not in fieldnames:
            raise ValueError(f"CSV must contain columns 'Scaffold' and 'SMILES'. Got: {fieldnames}")

        # Ensure SameS is the last column (overwrite if already exists)
        out_fieldnames = [c for c in fieldnames if c != "SameS"] + ["SameS"]

        writer = csv.DictWriter(fout, fieldnames=out_fieldnames, dialect=dialect)
        writer.writeheader()

        for row in reader:
            scaffold = (row.get("Scaffold") or "").strip()
            smiles = (row.get("SMILES") or "").strip()

            if not scaffold or not smiles:
                same = False
            else:
                same = is_center_symmetric_ap4(smiles, scaffold, strict=False)

            row["SameS"] = "true" if same else "false"
            writer.writerow({k: row.get(k, "") for k in out_fieldnames})


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("csv_in", type=Path, help="Input CSV file (must contain Scaffold and SMILES columns)")
    p.add_argument("-o", "--out", type=Path, default=None, help="Output CSV file path (optional)")
    args = p.parse_args()

    in_path: Path = args.csv_in
    if not in_path.exists():
        raise FileNotFoundError(f"Input not found: {in_path}")

    # If -o is provided -> write to out file.
    # If not provided -> update input file in-place via temp file + atomic replace.
    if args.out is not None:
        out_path = args.out
        process_csv(in_path, out_path)
        print(f"Saved: {out_path}")
        return 0

    # In-place update
    tmp_fd = None
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, dir=str(in_path.parent), suffix=".tmp", encoding="utf-8", newline=""
        ) as tf:
            tmp_path = Path(tf.name)

        process_csv(in_path, tmp_path)

        # atomic replace
        os.replace(str(tmp_path), str(in_path))
        print(f"Updated in-place: {in_path}")
        return 0
    finally:
        # If something failed before os.replace, cleanup temp file
        if tmp_path is not None and tmp_path.exists():
            # If we already replaced successfully, tmp_path won't exist
            try:
                tmp_path.unlink()
            except Exception:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
