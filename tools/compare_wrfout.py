#!/usr/bin/env python
"""Compare selected variables between two WRF-compatible NetCDF files."""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from dataclasses import asdict, dataclass, is_dataclass
from pathlib import Path
from typing import Iterable

import netCDF4
import numpy as np

try:
    from tools.tc_diagnostics import compare_tc_diagnostics
except ModuleNotFoundError:
    from tc_diagnostics import compare_tc_diagnostics


STRICT_CORE_VARIABLES = (
    "U",
    "V",
    "T",
    "PH",
    "MU",
    "P",
    "QVAPOR",
)

W_AND_HYDROMETEOR_VARIABLES = (
    "W",
    "QCLOUD",
    "QRAIN",
    "QICE",
    "QSNOW",
    "QGRAUP",
    "QNICE",
    "QNRAIN",
)

NEAR_SURFACE_VARIABLES = (
    "PSFC",
    "U10",
    "V10",
    "T2",
    "Q2",
)

PRECIPITATION_VARIABLES = (
    "RAINC",
    "RAINNC",
)

DEFAULT_VARIABLES = (
    *STRICT_CORE_VARIABLES,
    *W_AND_HYDROMETEOR_VARIABLES,
    *NEAR_SURFACE_VARIABLES,
    *PRECIPITATION_VARIABLES,
)

DEFAULT_THRESHOLDS = {
    **{name: 0.05 for name in STRICT_CORE_VARIABLES},
    **{name: 0.10 for name in W_AND_HYDROMETEOR_VARIABLES},
    **{name: 0.10 for name in NEAR_SURFACE_VARIABLES},
}


@dataclass(frozen=True)
class VariableComparison:
    variable: str
    status: str
    reference_shape: tuple[int, ...] | None = None
    candidate_shape: tuple[int, ...] | None = None
    rmse: float | None = None
    normalized_rmse: float | None = None
    max_abs_error: float | None = None
    threshold: float | None = None
    valid_count: int | None = None
    total_count: int | None = None
    message: str | None = None


@dataclass(frozen=True)
class ComparisonReport:
    reference: str
    candidate: str
    status: str
    summary: dict[str, int]
    variables: list[VariableComparison]
    diagnostics: dict[str, object]


@dataclass(frozen=True)
class ErrorMetrics:
    rmse: float
    normalized_rmse: float
    max_abs_error: float
    valid_count: int
    total_count: int


def _valid_numeric_pairs(reference: np.ndarray, candidate: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    ref = np.ma.masked_invalid(np.ma.asarray(reference, dtype=np.float64))
    cand = np.ma.masked_invalid(np.ma.asarray(candidate, dtype=np.float64))

    ref_data = np.asarray(ref.filled(np.nan), dtype=np.float64)
    cand_data = np.asarray(cand.filled(np.nan), dtype=np.float64)
    combined_mask = np.ma.getmaskarray(ref) | np.ma.getmaskarray(cand)
    finite_mask = np.isfinite(ref_data) & np.isfinite(cand_data) & ~combined_mask

    return ref_data[finite_mask], cand_data[finite_mask], int(ref_data.size)


def error_metrics(reference: np.ndarray, candidate: np.ndarray) -> ErrorMetrics:
    """Return RMSE, normalized RMSE, max absolute error, and sample counts."""

    ref_values, cand_values, total_count = _valid_numeric_pairs(reference, candidate)
    valid_count = int(ref_values.size)
    if valid_count == 0:
        return ErrorMetrics(
            rmse=math.nan,
            normalized_rmse=math.nan,
            max_abs_error=math.nan,
            valid_count=0,
            total_count=total_count,
        )

    diff = cand_values - ref_values
    rmse = float(math.sqrt(np.mean(diff * diff)))
    scale = float(math.sqrt(np.mean(ref_values * ref_values)))
    if scale == 0.0:
        norm = 0.0 if rmse == 0.0 else math.inf
    else:
        norm = rmse / scale
    max_abs = float(np.max(np.abs(diff)))
    return ErrorMetrics(
        rmse=rmse,
        normalized_rmse=norm,
        max_abs_error=max_abs,
        valid_count=valid_count,
        total_count=total_count,
    )


def normalized_rmse(reference: np.ndarray, candidate: np.ndarray) -> tuple[float, float, float]:
    """Return RMSE, normalized RMSE, and max absolute error."""

    metrics = error_metrics(reference, candidate)
    return metrics.rmse, metrics.normalized_rmse, metrics.max_abs_error


def compare_variable(
    reference: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    variable: str,
    thresholds: dict[str, float] | None = None,
) -> VariableComparison:
    threshold = thresholds.get(variable) if thresholds else None
    if variable not in reference.variables:
        return VariableComparison(
            variable=variable,
            status="missing_reference",
            threshold=threshold,
            message="variable missing from reference",
        )
    if variable not in candidate.variables:
        return VariableComparison(
            variable=variable,
            status="missing_candidate",
            reference_shape=tuple(int(v) for v in reference.variables[variable].shape),
            threshold=threshold,
            message="variable missing from candidate",
        )

    ref_data = reference.variables[variable][:]
    cand_data = candidate.variables[variable][:]
    if ref_data.shape != cand_data.shape:
        return VariableComparison(
            variable,
            "shape_mismatch",
            reference_shape=tuple(int(v) for v in ref_data.shape),
            candidate_shape=tuple(int(v) for v in cand_data.shape),
            threshold=threshold,
            message=f"reference shape {ref_data.shape} != candidate shape {cand_data.shape}",
        )

    try:
        metrics = error_metrics(ref_data, cand_data)
    except (TypeError, ValueError) as exc:
        shape = tuple(int(v) for v in ref_data.shape)
        return VariableComparison(
            variable=variable,
            status="non_numeric",
            reference_shape=shape,
            candidate_shape=shape,
            threshold=threshold,
            message=f"variable cannot be compared numerically: {exc}",
        )

    shape = tuple(int(v) for v in ref_data.shape)
    if metrics.valid_count == 0:
        return VariableComparison(
            variable=variable,
            status="no_valid_samples",
            reference_shape=shape,
            candidate_shape=shape,
            threshold=threshold,
            valid_count=metrics.valid_count,
            total_count=metrics.total_count,
            message="no finite unmasked sample pairs available for comparison",
        )

    status = "ok"
    message = None
    if threshold is not None and (
        math.isnan(metrics.normalized_rmse) or metrics.normalized_rmse > threshold
    ):
        status = "threshold_exceeded"
        message = f"normalized RMSE {metrics.normalized_rmse:g} exceeds threshold {threshold:g}"

    return VariableComparison(
        variable=variable,
        status=status,
        reference_shape=shape,
        candidate_shape=shape,
        rmse=metrics.rmse,
        normalized_rmse=metrics.normalized_rmse,
        max_abs_error=metrics.max_abs_error,
        threshold=threshold,
        valid_count=metrics.valid_count,
        total_count=metrics.total_count,
        message=message,
    )


def compare_files(
    reference_path: Path,
    candidate_path: Path,
    variables: Iterable[str] = DEFAULT_VARIABLES,
    thresholds: dict[str, float] | None = DEFAULT_THRESHOLDS,
    *,
    include_tc_diagnostics: bool = False,
    diagnostic_time_index: int = -1,
) -> ComparisonReport:
    with netCDF4.Dataset(reference_path) as reference, netCDF4.Dataset(candidate_path) as candidate:
        comparisons = [
            compare_variable(reference, candidate, variable, thresholds=thresholds)
            for variable in variables
        ]

    counts = Counter(comparison.status for comparison in comparisons)
    summary = dict(sorted(counts.items()))
    summary["total"] = len(comparisons)
    summary["failed"] = sum(1 for comparison in comparisons if comparison.status != "ok")
    status = "ok" if summary["failed"] == 0 else "failed"
    tc_diagnostics = (
        compare_tc_diagnostics(
            reference_path,
            candidate_path,
            time_index=diagnostic_time_index,
        )
        if include_tc_diagnostics
        else {
            "status": "pending",
            "message": "Run with --tc-diagnostics to report TC center, PSFC-min proxy MSLP, Vmax, and rainfall diagnostics.",
        }
    )
    return ComparisonReport(
        reference=str(reference_path),
        candidate=str(candidate_path),
        status=status,
        summary=summary,
        variables=comparisons,
        diagnostics={"tc": tc_diagnostics},
    )


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


def report_to_dict(report: ComparisonReport) -> dict:
    return _strict_json_value(report)


def report_to_json(report: ComparisonReport, *, pretty: bool = False) -> str:
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
    parser.add_argument("reference", type=Path, help="WRF reference NetCDF file")
    parser.add_argument("candidate", type=Path, help="TyWRF-Core candidate NetCDF file")
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_VARIABLES),
        help="Variables to compare",
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
        help="Report metrics without failing variables that exceed default thresholds",
    )
    parser.add_argument(
        "--tc-diagnostics",
        action="store_true",
        help="Include TC center, PSFC-min proxy MSLP, Vmax, and rainfall diagnostics",
    )
    parser.add_argument(
        "--diagnostic-time-index",
        type=int,
        default=-1,
        help="Time index used for optional TC diagnostics; defaults to the last record",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        thresholds = None if args.no_thresholds else parse_threshold_overrides(args.threshold)
    except ValueError as exc:
        parser.error(str(exc))

    report = compare_files(
        args.reference,
        args.candidate,
        args.variables,
        thresholds=thresholds,
        include_tc_diagnostics=args.tc_diagnostics,
        diagnostic_time_index=args.diagnostic_time_index,
    )
    report_json = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0 if report.status == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
