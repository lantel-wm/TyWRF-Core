#!/usr/bin/env python
"""Inspect WRF NetCDF reference files and report dimensions plus core fields."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import netCDF4


CORE_VARIABLES = (
    "Times",
    "XLAT",
    "XLONG",
    "HGT",
    "U",
    "V",
    "W",
    "PH",
    "PHB",
    "T",
    "MU",
    "MUB",
    "P",
    "PB",
    "QVAPOR",
    "QCLOUD",
    "QRAIN",
    "QICE",
    "QSNOW",
    "QGRAUP",
    "QNICE",
    "QNRAIN",
    "PSFC",
    "U10",
    "V10",
    "T2",
    "Q2",
    "RAINC",
    "RAINNC",
)


@dataclass(frozen=True)
class VariableSummary:
    dimensions: tuple[str, ...]
    shape: tuple[int, ...]
    dtype: str


@dataclass(frozen=True)
class FileSummary:
    path: str
    status: str
    dimensions: dict[str, int]
    present_core_variables: list[str]
    missing_core_variables: list[str]
    variables: dict[str, VariableSummary]


def summarize(path: Path, core_variables: Iterable[str] = CORE_VARIABLES) -> FileSummary:
    requested = tuple(core_variables)
    with netCDF4.Dataset(path) as dataset:
        dimensions = {name: len(dim) for name, dim in dataset.dimensions.items()}
        present = [name for name in requested if name in dataset.variables]
        missing = [name for name in requested if name not in dataset.variables]
        variables = {
            name: VariableSummary(
                dimensions=tuple(dataset.variables[name].dimensions),
                shape=tuple(int(v) for v in dataset.variables[name].shape),
                dtype=str(dataset.variables[name].dtype),
            )
            for name in present
        }
    return FileSummary(
        path=str(path),
        status="ok" if not missing else "missing_core_variables",
        dimensions=dimensions,
        present_core_variables=present,
        missing_core_variables=missing,
        variables=variables,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path, help="NetCDF files to inspect")
    parser.add_argument(
        "--core-variables",
        nargs="+",
        default=list(CORE_VARIABLES),
        help="Core variables expected in each reference file",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summaries = [summarize(path, args.core_variables) for path in args.paths]
    report = [asdict(summary) for summary in summaries]
    report_json = json.dumps(report, indent=2 if args.pretty else None)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0 if all(summary.status == "ok" for summary in summaries) else 1


if __name__ == "__main__":
    raise SystemExit(main())
