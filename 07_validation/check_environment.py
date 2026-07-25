#!/usr/bin/env python3
"""Collect a non-sensitive software and hardware capability report."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


ROOT = Path(__file__).resolve().parents[1]


def package_version(distribution: str) -> Optional[str]:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def command_version(command: str, args: list[str]) -> dict[str, Any]:
    executable = shutil.which(command)
    if executable is None:
        return {"available": False, "version": None}
    try:
        result = subprocess.run(
            [executable, *args],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        combined = (result.stdout or result.stderr).strip().splitlines()
        return {
            "available": result.returncode == 0,
            "version": combined[0] if combined else None,
            "returncode": result.returncode,
        }
    except Exception as exc:
        return {"available": False, "version": None, "error": str(exc)}


def cpu_model() -> Optional[str]:
    cpuinfo = Path("/proc/cpuinfo")
    if not cpuinfo.exists():
        return platform.processor() or None
    for line in cpuinfo.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.lower().startswith("model name"):
            return line.split(":", 1)[1].strip()
    return platform.processor() or None


def memory_gib() -> Optional[float]:
    meminfo = Path("/proc/meminfo")
    if not meminfo.exists():
        return None
    for line in meminfo.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("MemTotal:"):
            kib = float(line.split()[1])
            return round(kib / 1024 / 1024, 2)
    return None


def orca_version() -> dict[str, Any]:
    configured = os.environ.get("ORCA_BIN")
    # Do not fall back to a generic "orca" command: desktop systems often
    # provide the unrelated GNOME screen reader under that name.
    if not configured:
        return {
            "available": False,
            "configured_by": "ORCA_BIN",
            "version": None,
            "reason": "ORCA_BIN is not set",
        }
    executable = configured
    try:
        result = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        combined = (result.stdout or result.stderr).strip().splitlines()
        return {
            "available": result.returncode == 0,
            "configured_by": "ORCA_BIN",
            "version": combined[0] if combined else None,
            "returncode": result.returncode,
        }
    except Exception as exc:
        return {
            "available": False,
            "configured_by": "ORCA_BIN",
            "version": None,
            "error": str(exc),
        }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write a non-sensitive JSON report of the current environment."
    )
    parser.add_argument(
        "--output",
        default="reports/environment.json",
        help="Repository-relative output path.",
    )
    args = parser.parse_args()

    output = Path(args.output)
    if output.is_absolute():
        parser.error("--output must be relative to the repository root")
    output_path = ROOT / output

    report = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_runtime_capabilities",
        "reference_environment": "config/agent_interface.json",
        "python": {
            "version": platform.python_version(),
            "executable_kind": Path(sys.executable).name,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "cpu_model": cpu_model(),
            "logical_cpus": os.cpu_count(),
            "memory_gib": memory_gib(),
        },
        "python_packages": {
            name: package_version(distribution)
            for name, distribution in {
                "numpy": "numpy",
                "pandas": "pandas",
                "scipy": "scipy",
                "matplotlib": "matplotlib",
                "seaborn": "seaborn",
                "rdkit": "rdkit",
                "openpyxl": "openpyxl",
                "scikit_learn": "scikit-learn",
                "umap_learn": "umap-learn",
            }.items()
        },
        "commands": {
            "obabel": command_version("obabel", ["-V"]),
            "mpirun": command_version("mpirun", ["--version"]),
            "sbatch": command_version("sbatch", ["--version"]),
            "squeue": command_version("squeue", ["--version"]),
            "orca": orca_version(),
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(output.as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
