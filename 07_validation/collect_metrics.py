#!/usr/bin/env python3
"""Create machine-readable quality and scientific-performance metrics."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Union

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = "data/core/pic16_s1/input/pic16.csv"
DEFAULT_DATASET = "data/core/pic16_s1/final/pic16inpdftv22.csv"
DEFAULT_TD = "data/core/pic16_s1/intermediate/pic16inpout_tdv2.csv"
DEFAULT_V20 = "data/core/pic16_s1/intermediate/pic16inpdftv20.csv"
DEFAULT_V21 = "data/core/pic16_s1/intermediate/pic16inpdftv21.csv"
DEFAULT_REFERENCE = "05_analysis/data/publishinpdftv22.csv"
DEFAULT_OUTPUT = "reports/pic16_s1/metrics.json"

REQUIRED_COLUMNS = [
    "ID",
    "SMILES",
    "S1_energy_eV",
    "T1_energy_eV",
    "S1_fosc_D2",
    "S1_wavelength_nm",
    "min_positive_S1-Ti",
    "max_SOCME_S1_Ti_value",
    "max_SOCME_over_deltaE_S1_Ti_value",
    "HOMO(eV)",
    "LUMO(eV)",
    "Gap(eV)",
    "Dipole(D)",
    "Imaginary",
    "LogP",
    "TPSA",
    "SA_Score",
]

SUMMARY_COLUMNS = [
    "S1_energy_eV",
    "T1_energy_eV",
    "S1_fosc_D2",
    "S1_wavelength_nm",
    "S1-T1",
    "min_positive_S1-Ti",
    "max_SOCME_S1_Ti_value",
    "max_SOCME_over_deltaE_S1_Ti_value",
    "HOMO(eV)",
    "LUMO(eV)",
    "Gap(eV)",
    "Dipole(D)",
    "LogP",
    "TPSA",
    "SA_Score",
    "Molecular_Weight",
]


def repository_path(value: str, parser: argparse.ArgumentParser) -> tuple[Path, str]:
    path = Path(value)
    if path.is_absolute():
        parser.error("all input/output paths must be repository-relative")
    return ROOT / path, path.as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def count_rows(path: Optional[Path]) -> Optional[int]:
    if path is None or not path.exists():
        return None
    return int(len(pd.read_csv(path)))


def no_imaginary_mask(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip().str.lower()
    numeric = pd.to_numeric(series, errors="coerce")
    return text.isin({"no", "false", "0", "0.0"}) | numeric.eq(0.0)


def pareto_count(frame: pd.DataFrame, directions: list[str]) -> int:
    values = frame.to_numpy(dtype=float, copy=True)
    for index, direction in enumerate(directions):
        if direction == "maximize":
            values[:, index] *= -1
    nondominated = np.ones(len(values), dtype=bool)
    for index, point in enumerate(values):
        dominates_point = np.all(values <= point, axis=1) & np.any(
            values < point, axis=1
        )
        nondominated[index] = not bool(np.any(dominates_point))
    return int(nondominated.sum())


def describe(series: pd.Series) -> Dict[str, Union[float, int, None]]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return {
            "count": 0,
            "min": None,
            "p10": None,
            "p25": None,
            "median": None,
            "mean": None,
            "p75": None,
            "p90": None,
            "max": None,
        }
    quantiles = values.quantile([0.1, 0.25, 0.5, 0.75, 0.9])
    return {
        "count": int(values.count()),
        "min": float(values.min()),
        "p10": float(quantiles.loc[0.1]),
        "p25": float(quantiles.loc[0.25]),
        "median": float(quantiles.loc[0.5]),
        "mean": float(values.mean()),
        "p75": float(quantiles.loc[0.75]),
        "p90": float(quantiles.loc[0.9]),
        "max": float(values.max()),
    }


def safe_rate(
    numerator: Optional[int], denominator: Optional[int]
) -> Optional[float]:
    if numerator is None or denominator in (None, 0):
        return None
    return float(numerator / denominator)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize pipeline coverage, data quality, Pareto size, and property distributions."
    )
    parser.add_argument("--source", default=DEFAULT_SOURCE)
    parser.add_argument("--dataset", default=DEFAULT_DATASET)
    parser.add_argument("--td-summary", default=DEFAULT_TD)
    parser.add_argument("--v20", default=DEFAULT_V20)
    parser.add_argument("--v21", default=DEFAULT_V21)
    parser.add_argument("--reference", default=DEFAULT_REFERENCE)
    parser.add_argument("--reference-index", default="pic16")
    parser.add_argument(
        "--no-reference",
        action="store_true",
        help="Skip comparison with a reference molecule.",
    )
    parser.add_argument("--dataset-id", default="pic16_s1")
    parser.add_argument("--display-label", default="S1")
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--logp-min", type=float, default=1.0)
    parser.add_argument("--logp-max", type=float, default=5.0)
    parser.add_argument("--tpsa-min", type=float, default=20.0)
    parser.add_argument("--tpsa-max", type=float, default=120.0)
    args = parser.parse_args()

    source_path, source_rel = repository_path(args.source, parser)
    dataset_path, dataset_rel = repository_path(args.dataset, parser)
    td_path, td_rel = repository_path(args.td_summary, parser)
    v20_path, v20_rel = repository_path(args.v20, parser)
    v21_path, v21_rel = repository_path(args.v21, parser)
    if args.no_reference:
        reference_path: Optional[Path] = None
        reference_rel: Optional[str] = None
    else:
        reference_path, reference_rel = repository_path(args.reference, parser)
    output_path, output_rel = repository_path(args.output, parser)

    for required in (source_path, dataset_path):
        if not required.exists():
            parser.error("required input not found: {}".format(required.relative_to(ROOT)))

    source = pd.read_csv(source_path)
    data = pd.read_csv(dataset_path)
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in data]
    complete_columns = [column for column in REQUIRED_COLUMNS if column in data]
    complete_mask = (
        data[complete_columns].notna().all(axis=1)
        if complete_columns
        else pd.Series(False, index=data.index)
    )
    no_imaginary = (
        no_imaginary_mask(data["Imaginary"])
        if "Imaginary" in data
        else pd.Series(False, index=data.index)
    )

    pareto_columns_2d = [
        "min_positive_S1-Ti",
        "max_SOCME_over_deltaE_S1_Ti_value",
    ]
    pareto_2d_data = data[pareto_columns_2d].apply(
        pd.to_numeric, errors="coerce"
    ).dropna()
    pareto_2d = pareto_count(pareto_2d_data, ["minimize", "maximize"])

    pareto_columns_3d = [
        "min_positive_S1-Ti",
        "max_SOCME_S1_Ti_value",
        "max_SOCME_over_deltaE_S1_Ti_value",
    ]
    pareto_3d_data = data[pareto_columns_3d].apply(
        pd.to_numeric, errors="coerce"
    ).dropna()
    pareto_3d = pareto_count(
        pareto_3d_data, ["minimize", "maximize", "maximize"]
    )

    logp = pd.to_numeric(data.get("LogP"), errors="coerce")
    tpsa = pd.to_numeric(data.get("TPSA"), errors="coerce")
    in_window = (
        logp.between(args.logp_min, args.logp_max, inclusive="both")
        & tpsa.between(args.tpsa_min, args.tpsa_max, inclusive="both")
    )

    reference_benchmark: dict[str, Any]
    if reference_path is not None and reference_path.exists():
        reference_table = pd.read_csv(reference_path)
        if "index" not in reference_table:
            parser.error("reference table must contain an 'index' column")
        matches = reference_table[
            reference_table["index"].astype(str) == str(args.reference_index)
        ]
        if len(matches) != 1:
            parser.error(
                "reference index {!r} must match exactly one row; found {}".format(
                    args.reference_index, len(matches)
                )
            )
        reference = matches.iloc[0]
        ref_gap = float(reference["min_positive_S1-Ti"])
        ref_soc = float(reference["max_SOCME_S1_Ti_value"])
        ref_ratio = float(reference["max_SOCME_over_deltaE_S1_Ti_value"])
        candidate_gap = pd.to_numeric(
            data["min_positive_S1-Ti"], errors="coerce"
        )
        candidate_soc = pd.to_numeric(
            data["max_SOCME_S1_Ti_value"], errors="coerce"
        )
        candidate_ratio = pd.to_numeric(
            data["max_SOCME_over_deltaE_S1_Ti_value"], errors="coerce"
        )
        beats_gap_ratio = (candidate_gap <= ref_gap) & (
            candidate_ratio >= ref_ratio
        )
        beats_all_three = beats_gap_ratio & (candidate_soc >= ref_soc)
        reference_imaginary = str(reference.get("Imaginary", ""))
        reference_has_imaginary = reference_imaginary.strip().lower() in {
            "yes",
            "true",
            "1",
            "1.0",
        }
        reference_benchmark = {
            "reference_path": reference_rel,
            "reference_index": str(args.reference_index),
            "reference_values": {
                "min_positive_S1-Ti": ref_gap,
                "max_SOCME_S1_Ti_value": ref_soc,
                "max_SOCME_over_deltaE_S1_Ti_value": ref_ratio,
                "Imaginary": reference_imaginary,
            },
            "gap_and_ratio_better_or_equal": {
                "all_rows": int(beats_gap_ratio.sum()),
                "all_fraction": safe_rate(int(beats_gap_ratio.sum()), len(data)),
                "no_imaginary_rows": int((beats_gap_ratio & no_imaginary).sum()),
                "no_imaginary_fraction": safe_rate(
                    int((beats_gap_ratio & no_imaginary).sum()),
                    int(no_imaginary.sum()),
                ),
            },
            "gap_soc_and_ratio_better_or_equal": {
                "all_rows": int(beats_all_three.sum()),
                "all_fraction": safe_rate(int(beats_all_three.sum()), len(data)),
                "no_imaginary_rows": int((beats_all_three & no_imaginary).sum()),
                "no_imaginary_fraction": safe_rate(
                    int((beats_all_three & no_imaginary).sum()),
                    int(no_imaginary.sum()),
                ),
            },
            "comparison_contract": "Compare only data produced with a scientifically compatible quantum protocol.",
            "warning": (
                "The selected reference has an imaginary frequency and is an electronic-property benchmark, not a validated minimum."
                if reference_has_imaginary
                else None
            ),
        }
    else:
        reference_benchmark = {
            "reference_path": reference_rel,
            "reference_index": str(args.reference_index),
            "available": False,
            "reason": (
                "disabled_by_cli"
                if args.no_reference
                else "reference_file_not_found"
            ),
        }

    n_input = int(len(source))
    n_td = count_rows(td_path)
    n_v20 = count_rows(v20_path)
    n_v21 = count_rows(v21_path)
    n_v22 = int(len(data))
    n_complete = int(complete_mask.sum()) if not missing_columns else 0
    n_no_imaginary = int(no_imaginary.sum())
    input_hashes = {
        source_rel: sha256(source_path),
        dataset_rel: sha256(dataset_path),
    }
    if reference_path is not None and reference_path.exists():
        input_hashes[reference_rel] = sha256(reference_path)

    report: dict[str, Any] = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_id": args.dataset_id,
        "display_label": args.display_label,
        "inputs": {
            "source": source_rel,
            "dataset": dataset_rel,
            "td_summary": td_rel,
            "v20": v20_rel,
            "v21": v21_rel,
            "reference": reference_rel,
        },
        "sha256": input_hashes,
        "pipeline_counts": {
            "N_input": n_input,
            "N_td_summary": n_td,
            "N_v20": n_v20,
            "N_v21": n_v21,
            "N_v22": n_v22,
            "N_complete_v22": n_complete,
            "N_no_imaginary": n_no_imaginary,
            "N_opt_normal": None,
            "N_td_normal": None,
            "normal_termination_note": "Raw ORCA logs are intentionally not bundled; normal-termination counts require work/<batch>inp/*.log.",
        },
        "pipeline_rates": {
            "td_summary_coverage_of_input": safe_rate(n_td, n_input),
            "v20_coverage_of_input": safe_rate(n_v20, n_input),
            "v21_coverage_of_input": safe_rate(n_v21, n_input),
            "v22_coverage_of_input": safe_rate(n_v22, n_input),
            "complete_v22_rate": safe_rate(n_complete, n_input),
            "no_imaginary_rate": safe_rate(n_no_imaginary, n_v22),
            "opt_normal_termination_rate": None,
            "td_normal_termination_rate": None,
        },
        "quality": {
            "required_columns": REQUIRED_COLUMNS,
            "missing_required_columns": missing_columns,
            "all_ids_unique": bool(data["ID"].is_unique) if "ID" in data else False,
            "all_smiles_present": bool(data["SMILES"].notna().all())
            if "SMILES" in data
            else False,
        },
        "scientific_performance": {
            "pareto_gap_vs_soc_ratio": {
                "valid_rows": int(len(pareto_2d_data)),
                "nondominated_rows": pareto_2d,
                "fraction": safe_rate(pareto_2d, int(len(pareto_2d_data))),
                "directions": {
                    "min_positive_S1-Ti": "minimize",
                    "max_SOCME_over_deltaE_S1_Ti_value": "maximize",
                },
            },
            "pareto_gap_soc_ratio_3d": {
                "valid_rows": int(len(pareto_3d_data)),
                "nondominated_rows": pareto_3d,
                "fraction": safe_rate(pareto_3d, int(len(pareto_3d_data))),
                "directions": {
                    "min_positive_S1-Ti": "minimize",
                    "max_SOCME_S1_Ti_value": "maximize",
                    "max_SOCME_over_deltaE_S1_Ti_value": "maximize",
                },
            },
            "targeting_window": {
                "logp": [args.logp_min, args.logp_max],
                "tpsa_angstrom2": [args.tpsa_min, args.tpsa_max],
                "matching_rows": int(in_window.sum()),
                "fraction": safe_rate(int(in_window.sum()), n_v22),
                "interpretation": (
                    "Visualization default only; zero hits indicates that this window is not suitable for this chemical space."
                    if int(in_window.sum()) == 0
                    else "Visualization thresholds only; validate them against the application before using them as scientific constraints."
                ),
            },
            "reference_benchmark": reference_benchmark,
        },
        "property_summary": {
            column: describe(data[column])
            for column in SUMMARY_COLUMNS
            if column in data
        },
        "output": output_rel,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(output_rel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
