#!/usr/bin/env python
"""Write a reduced WRF-compatible core wrfout file from a source NetCDF file."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import netCDF4


CORE_WRFOUT_VARIABLES = (
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
class CoreWriteSummary:
    source: str
    destination: str
    copied_variables: list[str]
    missing_variables: list[str]
    copied_dimensions: list[str]


def _all_dimensions(source: netCDF4.Dataset) -> list[str]:
    return list(source.dimensions)


def _copy_attributes(source, destination) -> None:
    attrs = {name: source.getncattr(name) for name in source.ncattrs() if name != "_FillValue"}
    destination.setncatts(attrs)


def copy_core_wrfout(
    source_path: Path,
    destination_path: Path,
    variables: Iterable[str] = CORE_WRFOUT_VARIABLES,
    *,
    allow_missing: bool = False,
) -> CoreWriteSummary:
    requested = list(variables)
    destination_path.parent.mkdir(parents=True, exist_ok=True)

    with netCDF4.Dataset(source_path) as source:
        copied = [name for name in requested if name in source.variables]
        missing = [name for name in requested if name not in source.variables]
        if missing and not allow_missing:
            raise KeyError(f"source file is missing required variables: {', '.join(missing)}")

        dimensions = _all_dimensions(source)

        with netCDF4.Dataset(destination_path, "w", format=source.data_model) as destination:
            _copy_attributes(source, destination)
            destination.setncattr("TYWRF_CORE_WRFOUT", "true")
            destination.setncattr("TYWRF_CORE_SOURCE", str(source_path))
            destination.setncattr("TYWRF_CORE_VARIABLE_COUNT", len(copied))

            for name in dimensions:
                dimension = source.dimensions[name]
                size = None if dimension.isunlimited() else len(dimension)
                destination.createDimension(name, size)

            for name in copied:
                source_variable = source.variables[name]
                fill_value = getattr(source_variable, "_FillValue", None)
                destination_variable = destination.createVariable(
                    name,
                    source_variable.datatype,
                    source_variable.dimensions,
                    fill_value=fill_value,
                )
                _copy_attributes(source_variable, destination_variable)
                destination_variable[:] = source_variable[:]

    return CoreWriteSummary(
        source=str(source_path),
        destination=str(destination_path),
        copied_variables=copied,
        missing_variables=missing,
        copied_dimensions=dimensions,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Source WRF NetCDF file")
    parser.add_argument("destination", type=Path, help="Destination reduced wrfout file")
    parser.add_argument("--variables", nargs="+", default=list(CORE_WRFOUT_VARIABLES))
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = copy_core_wrfout(
        args.source,
        args.destination,
        variables=args.variables,
        allow_missing=args.allow_missing,
    )
    print(json.dumps(asdict(summary), indent=2 if args.pretty else None))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
