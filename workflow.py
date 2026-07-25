#!/usr/bin/env python3
"""Read-only command-line access to the repository's agent interface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "config" / "agent_interface.json"


def load_config() -> dict[str, Any]:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def emit(value: Any) -> None:
    print(json.dumps(value, ensure_ascii=False, indent=2))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Inspect workflow stages, paths, and adjustable parameters as JSON."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("show", help="Print the complete machine-readable interface.")
    subparsers.add_parser("stages", help="List all stable stage IDs.")
    subparsers.add_parser(
        "software", help="Show required software and reference hardware."
    )
    subparsers.add_parser("parameters", help="Show resource, quantum, and edit interfaces.")

    stage_parser = subparsers.add_parser("stage", help="Show one stage definition.")
    stage_parser.add_argument("stage_id", help="Stable ID such as S41_PARSE_TD")

    paths_parser = subparsers.add_parser("paths", help="Show repository-relative paths.")
    paths_parser.add_argument(
        "--check",
        action="store_true",
        help="Include whether every concrete path currently exists.",
    )

    args = parser.parse_args()
    config = load_config()

    if args.command == "show":
        emit(config)
        return 0

    if args.command == "stages":
        emit(
            [
                {
                    "stage_id": stage_id,
                    "entrypoint": stage["entrypoint"],
                    "external_mutation": stage["external_mutation"],
                }
                for stage_id, stage in config["stages"].items()
            ]
        )
        return 0

    if args.command == "software":
        emit(
            {
                "software": config["software"],
                "reference_hardware": config["reference_hardware"],
            }
        )
        return 0

    if args.command == "stage":
        stage = config["stages"].get(args.stage_id)
        if stage is None:
            parser.error(
                "unknown stage {!r}; choose from {}".format(
                    args.stage_id, ", ".join(config["stages"])
                )
            )
        emit({"stage_id": args.stage_id, **stage})
        return 0

    if args.command == "paths":
        if not args.check:
            emit(config["paths"])
            return 0
        rows = {}
        for name, relative_path in config["paths"].items():
            path = Path(relative_path)
            rows[name] = {
                "path": relative_path,
                "repository_relative": not path.is_absolute(),
                "exists": (ROOT / path).exists() if not path.is_absolute() else False,
            }
        emit(rows)
        return 0

    if args.command == "parameters":
        emit(
            {
                "scheduler_defaults": config["scheduler_defaults"],
                "scheduler_tools": config["scheduler_tools"],
                "quantum_profiles": config["quantum_profiles"],
                "parameter_edit_map": config["parameter_edit_map"],
                "analysis_entrypoints": config["analysis_entrypoints"],
                "analysis_parameter_interfaces": config[
                    "analysis_parameter_interfaces"
                ],
                "quality_contract": config["quality_contract"],
                "default_analysis_windows": config["default_analysis_windows"],
            }
        )
        return 0

    parser.error("unsupported command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
