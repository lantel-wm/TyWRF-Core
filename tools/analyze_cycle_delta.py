#!/usr/bin/env python
"""Analyze 10 minute WRF reference field deltas from cycle start to cycle end."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime
import json
import math
from pathlib import Path
from typing import Iterable

import netCDF4
import numpy as np

try:
    from tools.compare_wrfout import STRICT_CORE_VARIABLES
except ModuleNotFoundError:
    from compare_wrfout import STRICT_CORE_VARIABLES


WRF_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"
DEFAULT_DOMAIN = "d02"
DEFAULT_THRESHOLD = 0.05
DEFAULT_VARIABLES = STRICT_CORE_VARIABLES
DEFAULT_THRESHOLDS = {name: DEFAULT_THRESHOLD for name in DEFAULT_VARIABLES}


@dataclass(frozen=True)
class CycleDeltaVariable:
    variable: str
    status: str
    shape: tuple[int, ...] | None = None
    start_shape: tuple[int, ...] | None = None
    end_shape: tuple[int, ...] | None = None
    rmse_persistence: float | None = None
    normalized_rmse: float | None = None
    mean_abs_delta: float | None = None
    max_abs_delta: float | None = None
    finite_pair_count: int | None = None
    total_pair_count: int | None = None
    threshold: float | None = None
    exceeds_threshold: bool | None = None
    message: str | None = None


@dataclass(frozen=True)
class CycleDeltaReport:
    status: str
    start_file: str
    end_file: str
    domain: str | None
    start_time: str | None
    end_time: str | None
    variables: list[CycleDeltaVariable]
    summary: dict[str, object]


def parse_wrf_time(value: str) -> datetime:
    for fmt in (WRF_TIME_FORMAT, "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            pass
    raise ValueError(f"time must match YYYY-MM-DD_HH:MM:SS: {value}")


def format_wrf_time(value: datetime) -> str:
    return value.strftime(WRF_TIME_FORMAT)


def wrfout_filename(domain: str, valid_time: datetime) -> str:
    return f"wrfout_{domain}_{format_wrf_time(valid_time)}"


def resolve_cycle_files(
    *,
    start_file: Path | None = None,
    end_file: Path | None = None,
    reference_dir: Path | None = None,
    domain: str = DEFAULT_DOMAIN,
    start: str | datetime | None = None,
    end: str | datetime | None = None,
) -> tuple[Path, Path, str | None, str | None, str | None]:
    explicit = start_file is not None or end_file is not None
    named = reference_dir is not None or start is not None or end is not None
    if explicit and named:
        raise ValueError(
            "use either --start-file/--end-file or --reference-dir/--domain/--start/--end"
        )
    if explicit:
        if start_file is None or end_file is None:
            raise ValueError("--start-file and --end-file must be provided together")
        return start_file, end_file, None, None, None

    if reference_dir is None or start is None or end is None:
        raise ValueError(
            "provide --start-file/--end-file or --reference-dir, --start, and --end"
        )
    if domain not in {"d01", "d02"}:
        raise ValueError(f"unsupported domain: {domain}")

    start_time = parse_wrf_time(start) if isinstance(start, str) else start
    end_time = parse_wrf_time(end) if isinstance(end, str) else end
    if end_time <= start_time:
        raise ValueError("end time must be after start time")

    return (
        reference_dir / wrfout_filename(domain, start_time),
        reference_dir / wrfout_filename(domain, end_time),
        domain,
        format_wrf_time(start_time),
        format_wrf_time(end_time),
    )


def _valid_numeric_pairs(
    start: np.ndarray,
    end: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    start_masked = np.ma.masked_invalid(np.ma.asarray(start, dtype=np.float64))
    end_masked = np.ma.masked_invalid(np.ma.asarray(end, dtype=np.float64))

    start_data = np.asarray(start_masked.filled(np.nan), dtype=np.float64)
    end_data = np.asarray(end_masked.filled(np.nan), dtype=np.float64)
    combined_mask = np.ma.getmaskarray(start_masked) | np.ma.getmaskarray(end_masked)
    finite_mask = np.isfinite(start_data) & np.isfinite(end_data) & ~combined_mask

    return start_data[finite_mask], end_data[finite_mask], int(end_data.size)


def analyze_variable_delta(
    start_dataset: netCDF4.Dataset,
    end_dataset: netCDF4.Dataset,
    variable: str,
    *,
    thresholds: dict[str, float] | None = DEFAULT_THRESHOLDS,
) -> CycleDeltaVariable:
    threshold = thresholds.get(variable) if thresholds else None
    if variable not in start_dataset.variables:
        return CycleDeltaVariable(
            variable=variable,
            status="missing_start",
            threshold=threshold,
            message="variable missing from cycle start file",
        )
    if variable not in end_dataset.variables:
        return CycleDeltaVariable(
            variable=variable,
            status="missing_end",
            start_shape=tuple(int(v) for v in start_dataset.variables[variable].shape),
            threshold=threshold,
            message="variable missing from cycle end file",
        )

    start_data = start_dataset.variables[variable][:]
    end_data = end_dataset.variables[variable][:]
    start_shape = tuple(int(v) for v in start_data.shape)
    end_shape = tuple(int(v) for v in end_data.shape)
    if start_shape != end_shape:
        return CycleDeltaVariable(
            variable=variable,
            status="shape_mismatch",
            start_shape=start_shape,
            end_shape=end_shape,
            threshold=threshold,
            message=f"cycle start shape {start_shape} != cycle end shape {end_shape}",
        )

    try:
        start_values, end_values, total_pair_count = _valid_numeric_pairs(start_data, end_data)
    except (TypeError, ValueError) as exc:
        return CycleDeltaVariable(
            variable=variable,
            status="non_numeric",
            shape=start_shape,
            start_shape=start_shape,
            end_shape=end_shape,
            threshold=threshold,
            message=f"variable cannot be compared numerically: {exc}",
        )

    finite_pair_count = int(start_values.size)
    if finite_pair_count == 0:
        return CycleDeltaVariable(
            variable=variable,
            status="no_valid_pairs",
            shape=start_shape,
            start_shape=start_shape,
            end_shape=end_shape,
            finite_pair_count=0,
            total_pair_count=total_pair_count,
            threshold=threshold,
            message="no finite unmasked sample pairs available for delta analysis",
        )

    delta = start_values - end_values
    rmse = float(math.sqrt(np.mean(delta * delta)))
    scale = float(math.sqrt(np.mean(end_values * end_values)))
    if scale == 0.0:
        normalized = 0.0 if rmse == 0.0 else math.inf
    else:
        normalized = rmse / scale
    mean_abs = float(np.mean(np.abs(delta)))
    max_abs = float(np.max(np.abs(delta)))
    exceeds_threshold = (
        None
        if threshold is None
        else math.isnan(normalized) or normalized > threshold
    )
    message = (
        f"normalized persistence RMSE {normalized:g} exceeds strict threshold {threshold:g}"
        if exceeds_threshold
        else None
    )

    return CycleDeltaVariable(
        variable=variable,
        status="computed",
        shape=start_shape,
        start_shape=start_shape,
        end_shape=end_shape,
        rmse_persistence=rmse,
        normalized_rmse=normalized,
        mean_abs_delta=mean_abs,
        max_abs_delta=max_abs,
        finite_pair_count=finite_pair_count,
        total_pair_count=total_pair_count,
        threshold=threshold,
        exceeds_threshold=exceeds_threshold,
        message=message,
    )


def analyze_cycle_delta(
    start_path: Path,
    end_path: Path,
    *,
    variables: Iterable[str] = DEFAULT_VARIABLES,
    thresholds: dict[str, float] | None = DEFAULT_THRESHOLDS,
    domain: str | None = None,
    start_time: str | None = None,
    end_time: str | None = None,
) -> CycleDeltaReport:
    with netCDF4.Dataset(start_path) as start_dataset, netCDF4.Dataset(end_path) as end_dataset:
        analyzed_variables = [
            analyze_variable_delta(
                start_dataset,
                end_dataset,
                variable,
                thresholds=thresholds,
            )
            for variable in variables
        ]

    summary = _summary(analyzed_variables)
    unavailable_count = sum(1 for variable in analyzed_variables if variable.status != "computed")
    return CycleDeltaReport(
        status="computed" if unavailable_count == 0 else "incomplete",
        start_file=str(start_path),
        end_file=str(end_path),
        domain=domain,
        start_time=start_time,
        end_time=end_time,
        variables=analyzed_variables,
        summary=summary,
    )


def analyze_reference_cycle(
    reference_dir: Path,
    *,
    domain: str = DEFAULT_DOMAIN,
    start: str | datetime,
    end: str | datetime,
    variables: Iterable[str] = DEFAULT_VARIABLES,
    thresholds: dict[str, float] | None = DEFAULT_THRESHOLDS,
) -> CycleDeltaReport:
    start_path, end_path, resolved_domain, start_time, end_time = resolve_cycle_files(
        reference_dir=reference_dir,
        domain=domain,
        start=start,
        end=end,
    )
    return analyze_cycle_delta(
        start_path,
        end_path,
        variables=variables,
        thresholds=thresholds,
        domain=resolved_domain,
        start_time=start_time,
        end_time=end_time,
    )


def _summary(variables: list[CycleDeltaVariable]) -> dict[str, object]:
    computed = [variable for variable in variables if variable.status == "computed"]
    threshold_exceeded = [
        variable for variable in computed if variable.exceeds_threshold is True
    ]
    largest = sorted(
        computed,
        key=lambda variable: (
            variable.normalized_rmse is not None,
            float("-inf")
            if variable.normalized_rmse is None
            else variable.normalized_rmse,
        ),
        reverse=True,
    )
    first_unavailable = next(
        (variable.variable for variable in variables if variable.status != "computed"),
        None,
    )
    return {
        "total": len(variables),
        "computed": len(computed),
        "unavailable": len(variables) - len(computed),
        "strict_threshold_exceeded": len(threshold_exceeded),
        "first_failing_variable": (
            threshold_exceeded[0].variable if threshold_exceeded else None
        ),
        "first_unavailable_variable": first_unavailable,
        "largest_normalized_deltas": [
            {
                "variable": variable.variable,
                "normalized_rmse": variable.normalized_rmse,
                "rmse_persistence": variable.rmse_persistence,
            }
            for variable in largest
        ],
    }


def _strict_json_value(value):
    if is_dataclass(value) and not isinstance(value, type):
        return _strict_json_value(asdict(value))
    if isinstance(value, dict):
        return {key: _strict_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_strict_json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def report_to_dict(report: CycleDeltaReport) -> dict:
    return _strict_json_value(report)


def report_to_json(report: CycleDeltaReport, *, pretty: bool = False) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def parse_threshold_overrides(overrides: Iterable[str]) -> dict[str, float]:
    thresholds = dict(DEFAULT_THRESHOLDS)
    for override in overrides:
        if "=" not in override:
            raise ValueError(f"threshold override must use VARIABLE=VALUE: {override}")
        variable, raw_value = override.split("=", 1)
        variable = variable.strip()
        if not variable:
            raise ValueError(f"threshold override is missing a variable name: {override}")
        thresholds[variable] = float(raw_value)
    return thresholds


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_argument_group("input selection")
    input_group.add_argument("--start-file", type=Path, help="Cycle start WRF NetCDF file")
    input_group.add_argument("--end-file", type=Path, help="Cycle end WRF NetCDF file")
    input_group.add_argument("--reference-dir", type=Path, help="Directory with WRF reference files")
    input_group.add_argument("--domain", choices=("d01", "d02"), default=DEFAULT_DOMAIN)
    input_group.add_argument("--start", help="Cycle start time, for example 2025-07-26_00:00:00")
    input_group.add_argument("--end", help="Cycle end time, for example 2025-07-26_00:10:00")
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_VARIABLES),
        help="Variables to analyze; defaults to strict field gate variables",
    )
    parser.add_argument(
        "--threshold",
        action="append",
        default=[],
        metavar="VARIABLE=VALUE",
        help="Override or add a normalized RMSE threshold; repeat as needed",
    )
    parser.add_argument(
        "--no-thresholds",
        action="store_true",
        help="Report deltas without strict-threshold exceedance flags",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        start_file, end_file, domain, start_time, end_time = resolve_cycle_files(
            start_file=args.start_file,
            end_file=args.end_file,
            reference_dir=args.reference_dir,
            domain=args.domain,
            start=args.start,
            end=args.end,
        )
        thresholds = None if args.no_thresholds else parse_threshold_overrides(args.threshold)
    except ValueError as exc:
        parser.error(str(exc))

    report = analyze_cycle_delta(
        start_file,
        end_file,
        variables=args.variables,
        thresholds=thresholds,
        domain=domain,
        start_time=start_time,
        end_time=end_time,
    )
    report_json = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
