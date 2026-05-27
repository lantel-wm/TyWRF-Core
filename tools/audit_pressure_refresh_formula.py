#!/usr/bin/env python
"""Audit pressure-refresh perturbation/full-pressure formula semantics."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any

import netCDF4
import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.audit_pressure_refresh_candidate import (
    _candidate_metadata,
    _infer_target_region,
    _metrics_to_dict,
    _read_variable_data,
    _region_mask_for_data,
    error_metrics,
)


START_CLOSE_NORMALIZED_TOLERANCE = 1.0e-6
START_CLOSE_RATIO = 0.10
FULL_PRESSURE_BETTER_RATIO = 0.90


@dataclass(frozen=True)
class FormulaExpression:
    label: str
    status: str
    present: bool
    values: np.ndarray | None = None
    shape: tuple[int, ...] | None = None
    missing_inputs: tuple[str, ...] = ()
    message: str | None = None


@dataclass(frozen=True)
class FormulaComparison:
    name: str
    status: str
    reference_expression: str
    candidate_expression: str
    reference_present: bool
    candidate_present: bool
    reference_shape: tuple[int, ...] | None = None
    candidate_shape: tuple[int, ...] | None = None
    missing_inputs: list[str] | None = None
    rmse: float | None = None
    normalized_rmse: float | None = None
    max_abs_error: float | None = None
    valid_count: int | None = None
    total_count: int | None = None
    reference_finite_count: int | None = None
    candidate_finite_count: int | None = None
    sum_squared_error: float | None = None
    bias_mean: float | None = None
    region_split: dict[str, Any] | None = None
    message: str | None = None


@dataclass(frozen=True)
class PressureRefreshFormulaAudit:
    reference_end: str
    candidate: str
    candidate_start: str | None
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    summary: dict[str, Any]
    metadata: dict[str, Any]
    region: dict[str, Any]
    comparisons: list[FormulaComparison]
    diagnosis: dict[str, Any]


def _finite_bias_mean(reference: np.ndarray, candidate: np.ndarray) -> float | None:
    ref_values = np.ma.asarray(reference, dtype=np.float64)
    cand_values = np.ma.asarray(candidate, dtype=np.float64)
    ref_mask = np.ma.getmaskarray(np.ma.masked_invalid(ref_values))
    cand_mask = np.ma.getmaskarray(np.ma.masked_invalid(cand_values))
    values = np.asarray(cand_values.filled(np.nan) - ref_values.filled(np.nan))
    finite = np.isfinite(values) & ~ref_mask & ~cand_mask
    if not np.any(finite):
        return None
    return float(np.mean(values[finite]))


def _load_variable(
    dataset: netCDF4.Dataset | None,
    owner: str,
    name: str,
    *,
    time_index: int,
) -> FormulaExpression:
    label = f"{owner}.{name}"
    if dataset is None:
        return FormulaExpression(
            label=label,
            status="not_available",
            present=False,
            missing_inputs=(label,),
            message=f"{owner} dataset was not supplied",
        )
    if name not in dataset.variables:
        return FormulaExpression(
            label=label,
            status="not_available",
            present=False,
            missing_inputs=(label,),
            message=f"{label} is missing",
        )
    try:
        values = np.asarray(
            _read_variable_data(dataset, name, time_index=time_index),
            dtype=np.float64,
        )
    except (IndexError, TypeError, ValueError) as exc:
        return FormulaExpression(
            label=label,
            status="not_available",
            present=True,
            missing_inputs=(label,),
            message=f"{label} cannot be read: {exc}",
        )
    return FormulaExpression(
        label=label,
        status="available",
        present=True,
        values=values,
        shape=tuple(int(value) for value in values.shape),
    )


def _pressure_expression(
    dataset: netCDF4.Dataset | None,
    owner: str,
    *,
    full_pressure: bool,
    time_index: int,
) -> FormulaExpression:
    p = _load_variable(dataset, owner, "P", time_index=time_index)
    if not full_pressure:
        return p

    pb = _load_variable(dataset, owner, "PB", time_index=time_index)
    label = f"{owner}.P+PB"
    missing = [*p.missing_inputs, *pb.missing_inputs]
    if p.status != "available" or pb.status != "available":
        messages = [item for item in (p.message, pb.message) if item]
        return FormulaExpression(
            label=label,
            status="not_available",
            present=p.present and pb.present,
            missing_inputs=tuple(missing),
            message=(
                "; ".join(messages) if messages else "P/PB inputs are not available"
            ),
        )
    if p.values is None or pb.values is None:
        return FormulaExpression(
            label=label,
            status="not_available",
            present=p.present and pb.present,
            missing_inputs=tuple(missing) or (label,),
            message="P/PB values are not available",
        )
    if p.values.shape != pb.values.shape:
        return FormulaExpression(
            label=label,
            status="shape_mismatch",
            present=True,
            shape=tuple(int(value) for value in p.values.shape),
            message=(
                f"{p.label} shape {p.values.shape} != "
                f"{pb.label} shape {pb.values.shape}"
            ),
        )
    values = p.values + pb.values
    return FormulaExpression(
        label=label,
        status="available",
        present=True,
        values=values,
        shape=tuple(int(value) for value in values.shape),
    )


def _target_region_bias(
    reference: np.ndarray,
    candidate: np.ndarray,
    mask_2d: np.ndarray,
) -> float | None:
    target_mask = _region_mask_for_data(reference, mask_2d)
    return _finite_bias_mean(reference[target_mask], candidate[target_mask])


def _comparison_region_split(
    reference: np.ndarray,
    candidate: np.ndarray,
    region: dict[str, Any],
    global_sum_squared_error: float,
) -> dict[str, Any]:
    if reference.ndim < 3:
        return {
            "status": "not_applicable",
            "message": "target-region split requires a 3D pressure field",
        }
    if region.get("status") not in {"available", "available_count_mismatch"}:
        return {
            "status": "not_available",
            "message": region.get("message", "target region is not available"),
        }

    mask_2d = region.get("mask")
    if not isinstance(mask_2d, np.ndarray) or mask_2d.shape != reference.shape[-2:]:
        return {
            "status": "not_available",
            "message": "target region mask does not match expression horizontal shape",
        }

    target_mask = _region_mask_for_data(reference, mask_2d)
    non_target_mask = ~target_mask
    target_metrics = error_metrics(reference, candidate, selection_mask=target_mask)
    non_target_metrics = error_metrics(
        reference,
        candidate,
        selection_mask=non_target_mask,
    )
    target_sse = target_metrics.sum_squared_error
    if (
        math.isfinite(target_sse)
        and math.isfinite(global_sum_squared_error)
        and global_sum_squared_error > 0.0
    ):
        target_error_fraction = target_sse / global_sum_squared_error
    else:
        target_error_fraction = None

    return {
        "status": region["status"],
        "target_region": _metrics_to_dict(target_metrics),
        "non_target_region": _metrics_to_dict(non_target_metrics),
        "target_error_fraction": target_error_fraction,
        "target_region_dominates_global_error": (
            None if target_error_fraction is None else target_error_fraction > 0.5
        ),
        "target_region_bias_mean": _target_region_bias(
            reference,
            candidate,
            mask_2d,
        ),
        "non_target_region_bias_mean": _target_region_bias(
            reference,
            candidate,
            ~mask_2d,
        ),
        "target_column_count_matches_metadata": region.get(
            "target_column_count_matches_metadata"
        ),
        "message": region.get("message"),
    }


def _compare_expressions(
    name: str,
    reference: FormulaExpression,
    candidate: FormulaExpression,
    region: dict[str, Any],
) -> FormulaComparison:
    missing = [*reference.missing_inputs, *candidate.missing_inputs]
    if reference.status != "available" or candidate.status != "available":
        status = (
            "shape_mismatch"
            if "shape_mismatch" in {reference.status, candidate.status}
            else "not_available"
        )
        messages = [item for item in (reference.message, candidate.message) if item]
        return FormulaComparison(
            name=name,
            status=status,
            reference_expression=reference.label,
            candidate_expression=candidate.label,
            reference_present=reference.present,
            candidate_present=candidate.present,
            reference_shape=reference.shape,
            candidate_shape=candidate.shape,
            missing_inputs=missing,
            region_split={
                "status": "not_available",
                "message": "comparison inputs are not available",
            },
            message=(
                "; ".join(messages)
                if messages
                else "comparison inputs are not available"
            ),
        )
    if reference.values is None or candidate.values is None:
        return FormulaComparison(
            name=name,
            status="not_available",
            reference_expression=reference.label,
            candidate_expression=candidate.label,
            reference_present=reference.present,
            candidate_present=candidate.present,
            reference_shape=reference.shape,
            candidate_shape=candidate.shape,
            missing_inputs=missing,
            region_split={
                "status": "not_available",
                "message": "comparison values are not available",
            },
            message="comparison values are not available",
        )
    if reference.values.shape != candidate.values.shape:
        return FormulaComparison(
            name=name,
            status="shape_mismatch",
            reference_expression=reference.label,
            candidate_expression=candidate.label,
            reference_present=reference.present,
            candidate_present=candidate.present,
            reference_shape=tuple(int(value) for value in reference.values.shape),
            candidate_shape=tuple(int(value) for value in candidate.values.shape),
            missing_inputs=missing,
            region_split={
                "status": "not_available",
                "message": "comparison shapes do not match",
            },
            message=(
                f"{reference.label} shape {reference.values.shape} != "
                f"{candidate.label} shape {candidate.values.shape}"
            ),
        )

    metrics = error_metrics(reference.values, candidate.values)
    split = _comparison_region_split(
        reference.values,
        candidate.values,
        region,
        metrics.sum_squared_error,
    )
    return FormulaComparison(
        name=name,
        status="computed" if metrics.valid_count else "no_valid_samples",
        reference_expression=reference.label,
        candidate_expression=candidate.label,
        reference_present=reference.present,
        candidate_present=candidate.present,
        reference_shape=tuple(int(value) for value in reference.values.shape),
        candidate_shape=tuple(int(value) for value in candidate.values.shape),
        missing_inputs=missing,
        rmse=metrics.rmse,
        normalized_rmse=metrics.normalized_rmse,
        max_abs_error=metrics.max_abs_error,
        valid_count=metrics.valid_count,
        total_count=metrics.total_count,
        reference_finite_count=metrics.reference_finite_count,
        candidate_finite_count=metrics.candidate_finite_count,
        sum_squared_error=metrics.sum_squared_error,
        bias_mean=_finite_bias_mean(reference.values, candidate.values),
        region_split=split,
    )


def _horizontal_shape_for_p(
    reference_end: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    *,
    time_index: int,
) -> tuple[int, int] | None:
    if "P" not in reference_end.variables or "P" not in candidate.variables:
        return None
    try:
        ref_p = _read_variable_data(reference_end, "P", time_index=time_index)
        cand_p = _read_variable_data(candidate, "P", time_index=time_index)
    except (IndexError, TypeError, ValueError):
        return None
    if ref_p.shape != cand_p.shape or ref_p.ndim < 3:
        return None
    return int(ref_p.shape[-2]), int(ref_p.shape[-1])


def _comparison_by_name(
    comparisons: list[FormulaComparison],
    name: str,
) -> FormulaComparison | None:
    return next((item for item in comparisons if item.name == name), None)


def _metric_for_diagnosis(comparison: FormulaComparison | None) -> float | None:
    if comparison is None or comparison.status != "computed":
        return None
    split = comparison.region_split or {}
    if split.get("status") in {"available", "available_count_mismatch"}:
        target = split.get("target_region")
        if isinstance(target, dict):
            value = target.get("normalized_rmse")
            if isinstance(value, (int, float)) and math.isfinite(float(value)):
                return float(value)
    value = comparison.normalized_rmse
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return float(value)
    return None


def _rmse(comparison: FormulaComparison | None) -> float | None:
    if comparison is None or comparison.status != "computed":
        return None
    if isinstance(comparison.rmse, (int, float)) and math.isfinite(
        float(comparison.rmse)
    ):
        return float(comparison.rmse)
    return None


def _rmse_ratio(candidate: float | None, baseline: float | None) -> float | None:
    if candidate is None or baseline is None:
        return None
    if baseline == 0.0:
        return 0.0 if candidate == 0.0 else math.inf
    return candidate / baseline


def _build_diagnosis(
    comparisons: list[FormulaComparison],
    *,
    start_close_normalized_tolerance: float,
    start_close_ratio: float,
) -> dict[str, Any]:
    p_vs_ref = _comparison_by_name(comparisons, "candidate_p_vs_reference_p")
    full_vs_ref = _comparison_by_name(
        comparisons,
        "candidate_full_pressure_vs_reference_full_pressure",
    )
    p_vs_start = _comparison_by_name(comparisons, "candidate_p_vs_start_p")
    p_vs_ref_full = _comparison_by_name(
        comparisons,
        "candidate_p_vs_reference_full_pressure",
    )
    full_vs_ref_p = _comparison_by_name(
        comparisons,
        "candidate_full_pressure_vs_reference_p",
    )

    p_region_split = p_vs_ref.region_split if p_vs_ref else None
    if (
        p_region_split
        and p_region_split.get("status") in {"available", "available_count_mismatch"}
    ):
        target_region_p_error_dominates = p_region_split.get(
            "target_region_dominates_global_error"
        )
        target_region_error_fraction = p_region_split.get("target_error_fraction")
        target_region_bias_mean = p_region_split.get("target_region_bias_mean")
    else:
        target_region_p_error_dominates = None
        target_region_error_fraction = None
        target_region_bias_mean = None

    p_rmse = _rmse(p_vs_ref)
    full_rmse = _rmse(full_vs_ref)
    full_pressure_to_perturbation_rmse_ratio = _rmse_ratio(full_rmse, p_rmse)
    full_pressure_better = (
        None
        if full_pressure_to_perturbation_rmse_ratio is None
        else full_pressure_to_perturbation_rmse_ratio <= FULL_PRESSURE_BETTER_RATIO
    )

    start_metric = _metric_for_diagnosis(p_vs_start)
    reference_metric = _metric_for_diagnosis(p_vs_ref)
    if start_metric is None:
        candidate_p_close_to_start = None
    else:
        candidate_p_close_to_start = start_metric <= start_close_normalized_tolerance
        if not candidate_p_close_to_start and reference_metric is not None:
            candidate_p_close_to_start = (
                start_metric <= reference_metric * start_close_ratio
            )

    ref_full_rmse = _rmse(p_vs_ref_full)
    candidate_p_to_reference_full_rmse_ratio = _rmse_ratio(ref_full_rmse, p_rmse)
    candidate_p_better_matches_reference_full_pressure = (
        None
        if candidate_p_to_reference_full_rmse_ratio is None
        else candidate_p_to_reference_full_rmse_ratio <= FULL_PRESSURE_BETTER_RATIO
    )

    full_vs_ref_p_rmse = _rmse(full_vs_ref_p)
    candidate_full_to_reference_p_rmse_ratio = _rmse_ratio(full_vs_ref_p_rmse, p_rmse)
    candidate_full_better_matches_reference_p = (
        None
        if candidate_full_to_reference_p_rmse_ratio is None
        else candidate_full_to_reference_p_rmse_ratio <= FULL_PRESSURE_BETTER_RATIO
    )

    return {
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "target_region_p_error_dominates": target_region_p_error_dominates,
        "p_target_region_error_fraction": target_region_error_fraction,
        "full_pressure_better_than_perturbation": full_pressure_better,
        "full_pressure_to_perturbation_rmse_ratio": (
            full_pressure_to_perturbation_rmse_ratio
        ),
        "candidate_p_close_to_start_p": candidate_p_close_to_start,
        "candidate_p_bias_mean": None if p_vs_ref is None else p_vs_ref.bias_mean,
        "target_region_bias_mean": target_region_bias_mean,
        "candidate_p_target_or_global_start_normalized_rmse": start_metric,
        "candidate_p_target_or_global_reference_normalized_rmse": reference_metric,
        "candidate_p_better_matches_reference_full_pressure": (
            candidate_p_better_matches_reference_full_pressure
        ),
        "candidate_p_to_reference_full_rmse_ratio": (
            candidate_p_to_reference_full_rmse_ratio
        ),
        "candidate_full_pressure_better_matches_reference_p": (
            candidate_full_better_matches_reference_p
        ),
        "candidate_full_to_reference_p_rmse_ratio": (
            candidate_full_to_reference_p_rmse_ratio
        ),
        "message": (
            "diagnostic-only formula/source audit; reference comparisons explain "
            "pressure semantics and never generate a validation pass"
        ),
    }


def _build_comparisons(
    reference_end: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    candidate_start: netCDF4.Dataset | None,
    region: dict[str, Any],
    *,
    time_index: int,
    start_time_index: int,
) -> list[FormulaComparison]:
    reference_p = _pressure_expression(
        reference_end,
        "reference_end",
        full_pressure=False,
        time_index=time_index,
    )
    candidate_p = _pressure_expression(
        candidate,
        "candidate",
        full_pressure=False,
        time_index=time_index,
    )
    reference_full = _pressure_expression(
        reference_end,
        "reference_end",
        full_pressure=True,
        time_index=time_index,
    )
    candidate_full = _pressure_expression(
        candidate,
        "candidate",
        full_pressure=True,
        time_index=time_index,
    )
    start_p = _pressure_expression(
        candidate_start,
        "candidate_start",
        full_pressure=False,
        time_index=start_time_index,
    )
    start_full = _pressure_expression(
        candidate_start,
        "candidate_start",
        full_pressure=True,
        time_index=start_time_index,
    )

    return [
        _compare_expressions(
            "candidate_p_vs_reference_p",
            reference_p,
            candidate_p,
            region,
        ),
        _compare_expressions(
            "candidate_full_pressure_vs_reference_full_pressure",
            reference_full,
            candidate_full,
            region,
        ),
        _compare_expressions(
            "candidate_p_vs_start_p",
            start_p,
            candidate_p,
            region,
        ),
        _compare_expressions(
            "candidate_full_pressure_vs_start_full_pressure",
            start_full,
            candidate_full,
            region,
        ),
        _compare_expressions(
            "candidate_full_pressure_vs_reference_p",
            reference_p,
            candidate_full,
            region,
        ),
        _compare_expressions(
            "candidate_p_vs_reference_full_pressure",
            reference_full,
            candidate_p,
            region,
        ),
    ]


def audit_pressure_refresh_formula(
    reference_end_path: Path,
    candidate_path: Path,
    candidate_start_path: Path | None = None,
    *,
    time_index: int = -1,
    start_time_index: int = -1,
    start_close_normalized_tolerance: float = START_CLOSE_NORMALIZED_TOLERANCE,
    start_close_ratio: float = START_CLOSE_RATIO,
) -> PressureRefreshFormulaAudit:
    with netCDF4.Dataset(reference_end_path) as reference_end, netCDF4.Dataset(
        candidate_path
    ) as candidate:
        if candidate_start_path is None:
            candidate_start = None
            metadata = _candidate_metadata(candidate)
            horizontal_shape = _horizontal_shape_for_p(
                reference_end,
                candidate,
                time_index=time_index,
            )
            region = _infer_target_region(candidate, horizontal_shape)
            comparisons = _build_comparisons(
                reference_end,
                candidate,
                candidate_start,
                region,
                time_index=time_index,
                start_time_index=start_time_index,
            )
        else:
            with netCDF4.Dataset(candidate_start_path) as candidate_start_dataset:
                metadata = _candidate_metadata(candidate)
                horizontal_shape = _horizontal_shape_for_p(
                    reference_end,
                    candidate,
                    time_index=time_index,
                )
                region = _infer_target_region(candidate, horizontal_shape)
                comparisons = _build_comparisons(
                    reference_end,
                    candidate,
                    candidate_start_dataset,
                    region,
                    time_index=time_index,
                    start_time_index=start_time_index,
                )

    public_region = {key: value for key, value in region.items() if key != "mask"}
    counts: dict[str, int] = {}
    for comparison in comparisons:
        counts[comparison.status] = counts.get(comparison.status, 0) + 1

    diagnosis = _build_diagnosis(
        comparisons,
        start_close_normalized_tolerance=start_close_normalized_tolerance,
        start_close_ratio=start_close_ratio,
    )
    summary = {
        **dict(sorted(counts.items())),
        "total": len(comparisons),
        "computed": counts.get("computed", 0),
        "not_available": counts.get("not_available", 0),
        "shape_mismatch": counts.get("shape_mismatch", 0),
        "no_valid_samples": counts.get("no_valid_samples", 0),
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "target_region_p_error_dominates": diagnosis[
            "target_region_p_error_dominates"
        ],
        "full_pressure_better_than_perturbation": diagnosis[
            "full_pressure_better_than_perturbation"
        ],
        "candidate_p_close_to_start_p": diagnosis["candidate_p_close_to_start_p"],
        "candidate_p_bias_mean": diagnosis["candidate_p_bias_mean"],
        "target_region_bias_mean": diagnosis["target_region_bias_mean"],
    }
    status = "computed" if counts.get("computed", 0) else "computed_with_flags"
    if counts.get("shape_mismatch", 0) or counts.get("no_valid_samples", 0):
        status = "computed_with_flags"

    return PressureRefreshFormulaAudit(
        reference_end=str(reference_end_path),
        candidate=str(candidate_path),
        candidate_start=(
            None if candidate_start_path is None else str(candidate_start_path)
        ),
        status=status,
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        summary=summary,
        metadata=metadata,
        region=public_region,
        comparisons=comparisons,
        diagnosis=diagnosis,
    )


def _strict_json_value(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _strict_json_value(asdict(value))
    if isinstance(value, dict):
        return {key: _strict_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_strict_json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return _strict_json_value(value.tolist())
    if isinstance(value, np.generic):
        return _strict_json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def report_to_dict(report: PressureRefreshFormulaAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: PressureRefreshFormulaAudit,
    *,
    pretty: bool = False,
) -> str:
    return json.dumps(
        report_to_dict(report),
        indent=2 if pretty else None,
        allow_nan=False,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-end",
        required=True,
        type=Path,
        help="WRF reference wrfout at the candidate end time",
    )
    parser.add_argument(
        "--candidate",
        required=True,
        type=Path,
        help="TyWRF candidate wrfout at the same end time",
    )
    parser.add_argument(
        "--candidate-start",
        type=Path,
        help="Optional TyWRF/WRF candidate start-state wrfout for P-stasis checks",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for end-time Time-leading variables",
    )
    parser.add_argument(
        "--start-time-index",
        type=int,
        default=-1,
        help="Time index used for candidate-start Time-leading variables",
    )
    parser.add_argument(
        "--start-close-normalized-tolerance",
        type=float,
        default=START_CLOSE_NORMALIZED_TOLERANCE,
        help="Absolute normalized-RMSE tolerance for candidate P vs start P closeness",
    )
    parser.add_argument(
        "--start-close-ratio",
        type=float,
        default=START_CLOSE_RATIO,
        help="Relative tolerance versus candidate P end-reference error for start closeness",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Write the JSON report to this path",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_pressure_refresh_formula(
            args.reference_end,
            args.candidate,
            args.candidate_start,
            time_index=args.time_index,
            start_time_index=args.start_time_index,
            start_close_normalized_tolerance=args.start_close_normalized_tolerance,
            start_close_ratio=args.start_close_ratio,
        )
    except (OSError, ValueError, IndexError) as exc:
        parser.error(str(exc))

    text = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
