#!/usr/bin/env python
"""Audit target-region vertical pressure bias for pressure-refresh diagnostics."""

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
    _read_variable_data,
    _region_mask_for_data,
    error_metrics,
)


DEFAULT_WORST_LEVEL_COUNT = 5
COMPANION_FIELDS = ("PB", "P+PB", "MU", "MUB", "PH", "PHB", "HGT", "PSFC")


@dataclass(frozen=True)
class PerLevelPressureMetrics:
    level: int
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    mean_candidate_p: float | None = None
    mean_reference_p: float | None = None
    mean_diff: float | None = None
    rmse: float | None = None
    normalized_rmse: float | None = None
    max_abs_error: float | None = None
    valid_count: int | None = None
    total_count: int | None = None
    reference_finite_count: int | None = None
    candidate_finite_count: int | None = None
    message: str | None = None


@dataclass(frozen=True)
class CompanionFieldSummary:
    field: str
    level: int
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    mean_candidate: float | None = None
    mean_reference: float | None = None
    mean_diff: float | None = None
    valid_count: int | None = None
    reference_shape: tuple[int, ...] | None = None
    candidate_shape: tuple[int, ...] | None = None
    alignment: str | None = None
    message: str | None = None


@dataclass(frozen=True)
class SourceEvolutionSummary:
    level: int
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    mean_source_p: float | None = None
    mean_candidate_p: float | None = None
    mean_reference_p: float | None = None
    mean_candidate_minus_source_p: float | None = None
    mean_reference_minus_source_p: float | None = None
    mean_candidate_evolution_minus_reference_evolution_p: float | None = None
    valid_count: int | None = None
    source_shape: tuple[int, ...] | None = None
    message: str | None = None


@dataclass(frozen=True)
class WorstLevelAttribution:
    level: int
    rank_by_rmse: int | None
    rank_by_abs_bias: int | None
    diagnostic_only: bool
    candidate_model_pass: str
    pressure_metrics: PerLevelPressureMetrics
    companion_fields: list[CompanionFieldSummary]
    source_evolution: SourceEvolutionSummary


@dataclass(frozen=True)
class PressureRefreshVerticalBiasAudit:
    reference_end: str
    candidate_end: str
    source_start: str | None
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    summary: dict[str, Any]
    metadata: dict[str, Any]
    region: dict[str, Any]
    per_level_p: list[PerLevelPressureMetrics]
    worst_levels_by_rmse: list[int]
    worst_levels_by_abs_bias: list[int]
    worst_level_attribution: list[WorstLevelAttribution]
    diagnosis: dict[str, Any]


def _finite_mean(values: np.ndarray, mask: np.ndarray) -> float | None:
    masked = np.ma.masked_invalid(np.ma.asarray(values, dtype=np.float64))
    array = np.asarray(masked.filled(np.nan), dtype=np.float64)
    finite = np.isfinite(array) & ~np.ma.getmaskarray(masked) & mask
    if not np.any(finite):
        return None
    return float(np.mean(array[finite]))


def _finite_diff_mean(
    reference: np.ndarray,
    candidate: np.ndarray,
    mask: np.ndarray,
) -> tuple[float | None, int]:
    ref_masked = np.ma.masked_invalid(np.ma.asarray(reference, dtype=np.float64))
    cand_masked = np.ma.masked_invalid(np.ma.asarray(candidate, dtype=np.float64))
    ref_values = np.asarray(ref_masked.filled(np.nan), dtype=np.float64)
    cand_values = np.asarray(cand_masked.filled(np.nan), dtype=np.float64)
    finite = (
        np.isfinite(ref_values)
        & np.isfinite(cand_values)
        & ~np.ma.getmaskarray(ref_masked)
        & ~np.ma.getmaskarray(cand_masked)
        & mask
    )
    if not np.any(finite):
        return None, 0
    return float(np.mean(cand_values[finite] - ref_values[finite])), int(
        np.count_nonzero(finite)
    )


def _load_optional_variable(
    dataset: netCDF4.Dataset,
    name: str,
    *,
    time_index: int,
) -> np.ndarray | None:
    if name not in dataset.variables:
        return None
    return np.asarray(
        _read_variable_data(dataset, name, time_index=time_index),
        dtype=np.float64,
    )


def _read_p_pair(
    reference_end: netCDF4.Dataset,
    candidate_end: netCDF4.Dataset,
    *,
    time_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    if "P" not in reference_end.variables:
        raise ValueError("reference-end file is missing P")
    if "P" not in candidate_end.variables:
        raise ValueError("candidate-end file is missing P")
    reference_p = np.asarray(
        _read_variable_data(reference_end, "P", time_index=time_index),
        dtype=np.float64,
    )
    candidate_p = np.asarray(
        _read_variable_data(candidate_end, "P", time_index=time_index),
        dtype=np.float64,
    )
    if reference_p.shape != candidate_p.shape:
        raise ValueError(
            f"P shape mismatch: reference {reference_p.shape} != "
            f"candidate {candidate_p.shape}"
        )
    if reference_p.ndim < 3:
        raise ValueError(
            f"P must be at least 3D after Time selection, got {reference_p.shape}"
        )
    return reference_p, candidate_p


def _target_mask_for_p_level(
    p: np.ndarray,
    mask_2d: np.ndarray,
    level: int,
) -> np.ndarray:
    return _region_mask_for_data(p, mask_2d)[level]


def _per_level_metrics(
    reference_p: np.ndarray,
    candidate_p: np.ndarray,
    mask_2d: np.ndarray,
) -> list[PerLevelPressureMetrics]:
    levels: list[PerLevelPressureMetrics] = []
    for level in range(int(reference_p.shape[-3])):
        ref_slice = reference_p[level]
        cand_slice = candidate_p[level]
        mask = _target_mask_for_p_level(reference_p, mask_2d, level)
        metrics = error_metrics(ref_slice, cand_slice, selection_mask=mask)
        mean_diff, _ = _finite_diff_mean(ref_slice, cand_slice, mask)
        status = "computed" if metrics.valid_count else "no_valid_samples"
        levels.append(
            PerLevelPressureMetrics(
                level=level,
                status=status,
                diagnostic_only=True,
                candidate_model_pass="not_applicable",
                mean_candidate_p=_finite_mean(cand_slice, mask),
                mean_reference_p=_finite_mean(ref_slice, mask),
                mean_diff=mean_diff,
                rmse=metrics.rmse,
                normalized_rmse=metrics.normalized_rmse,
                max_abs_error=metrics.max_abs_error,
                valid_count=metrics.valid_count,
                total_count=metrics.total_count,
                reference_finite_count=metrics.reference_finite_count,
                candidate_finite_count=metrics.candidate_finite_count,
                message=(
                    None
                    if metrics.valid_count
                    else "no finite target-region P sample pairs are available"
                ),
            )
        )
    return levels


def _finite_metric_value(value: float | None) -> float | None:
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return float(value)
    return None


def _rank_levels_by_rmse(
    levels: list[PerLevelPressureMetrics],
    count: int,
) -> list[int]:
    ranked = [
        item
        for item in levels
        if item.valid_count and _finite_metric_value(item.rmse) is not None
    ]
    ranked.sort(key=lambda item: (-float(item.rmse), item.level))
    return [item.level for item in ranked[:count]]


def _rank_levels_by_abs_bias(
    levels: list[PerLevelPressureMetrics],
    count: int,
) -> list[int]:
    ranked = [
        item
        for item in levels
        if item.valid_count and _finite_metric_value(item.mean_diff) is not None
    ]
    ranked.sort(key=lambda item: (-abs(float(item.mean_diff)), item.level))
    return [item.level for item in ranked[:count]]


def _level_rank_map(levels: list[int]) -> dict[int, int]:
    return {level: index + 1 for index, level in enumerate(levels)}


def _align_field_to_p_level(
    field: np.ndarray,
    *,
    level: int,
    p_shape: tuple[int, ...],
) -> tuple[np.ndarray | None, str | None, str | None]:
    if field.shape[-2:] != p_shape[-2:]:
        return None, None, (
            f"horizontal shape {field.shape[-2:]} does not match P {p_shape[-2:]}"
        )
    p_nz = int(p_shape[-3])
    if field.ndim >= 3:
        field_nz = int(field.shape[-3])
        if field_nz == p_nz:
            return field[level], "mass_level_same_index", None
        if field_nz == p_nz + 1 and level + 1 < field_nz:
            return (
                0.5 * (field[level] + field[level + 1]),
                "staggered_vertical_adjacent_mean",
                None,
            )
        return None, None, (
            f"vertical shape {field_nz} cannot align to P level count {p_nz}"
        )
    if field.ndim == 2:
        return field, "horizontal_2d_reused_for_level", None
    return None, None, f"unsupported field shape {field.shape}"


def _companion_expression(
    dataset: netCDF4.Dataset,
    field: str,
    *,
    time_index: int,
) -> tuple[np.ndarray | None, str | None]:
    if field == "P+PB":
        p = _load_optional_variable(dataset, "P", time_index=time_index)
        pb = _load_optional_variable(dataset, "PB", time_index=time_index)
        if p is None or pb is None:
            missing = [name for name, value in (("P", p), ("PB", pb)) if value is None]
            return None, f"missing inputs: {', '.join(missing)}"
        if p.shape != pb.shape:
            return None, f"P shape {p.shape} != PB shape {pb.shape}"
        return p + pb, None

    values = _load_optional_variable(dataset, field, time_index=time_index)
    if values is None:
        return None, f"{field} is missing"
    return values, None


def _companion_summary(
    reference_end: netCDF4.Dataset,
    candidate_end: netCDF4.Dataset,
    field: str,
    *,
    level: int,
    p_shape: tuple[int, ...],
    mask_2d: np.ndarray,
    time_index: int,
) -> CompanionFieldSummary:
    ref_values, ref_error = _companion_expression(
        reference_end,
        field,
        time_index=time_index,
    )
    cand_values, cand_error = _companion_expression(
        candidate_end,
        field,
        time_index=time_index,
    )
    messages = [message for message in (ref_error, cand_error) if message]
    if ref_values is None or cand_values is None:
        return CompanionFieldSummary(
            field=field,
            level=level,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_shape=None if ref_values is None else tuple(ref_values.shape),
            candidate_shape=None if cand_values is None else tuple(cand_values.shape),
            message="; ".join(messages) if messages else "companion field unavailable",
        )
    if ref_values.shape != cand_values.shape:
        return CompanionFieldSummary(
            field=field,
            level=level,
            status="shape_mismatch",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_shape=tuple(ref_values.shape),
            candidate_shape=tuple(cand_values.shape),
            message=(
                f"reference shape {ref_values.shape} != "
                f"candidate shape {cand_values.shape}"
            ),
        )

    ref_aligned, ref_alignment, ref_align_error = _align_field_to_p_level(
        ref_values,
        level=level,
        p_shape=p_shape,
    )
    cand_aligned, cand_alignment, cand_align_error = _align_field_to_p_level(
        cand_values,
        level=level,
        p_shape=p_shape,
    )
    messages = [message for message in (ref_align_error, cand_align_error) if message]
    if ref_aligned is None or cand_aligned is None:
        return CompanionFieldSummary(
            field=field,
            level=level,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_shape=tuple(ref_values.shape),
            candidate_shape=tuple(cand_values.shape),
            message=(
                "; ".join(messages) if messages else "field cannot align to P level"
            ),
        )

    mean_diff, valid_count = _finite_diff_mean(ref_aligned, cand_aligned, mask_2d)
    return CompanionFieldSummary(
        field=field,
        level=level,
        status="computed" if valid_count else "no_valid_samples",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        mean_candidate=_finite_mean(cand_aligned, mask_2d),
        mean_reference=_finite_mean(ref_aligned, mask_2d),
        mean_diff=mean_diff,
        valid_count=valid_count,
        reference_shape=tuple(ref_values.shape),
        candidate_shape=tuple(cand_values.shape),
        alignment=(
            ref_alignment if ref_alignment == cand_alignment else "mixed_alignment"
        ),
        message=(
            None
            if valid_count
            else "no finite target-region sample pairs are available"
        ),
    )


def _source_evolution_summary(
    reference_p: np.ndarray,
    candidate_p: np.ndarray,
    source_start: netCDF4.Dataset | None,
    *,
    level: int,
    mask_2d: np.ndarray,
    source_time_index: int,
) -> SourceEvolutionSummary:
    if source_start is None:
        return SourceEvolutionSummary(
            level=level,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            message="source/start file was not supplied",
        )
    if "P" not in source_start.variables:
        return SourceEvolutionSummary(
            level=level,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            message="source/start file is missing P",
        )
    source_p = np.asarray(
        _read_variable_data(source_start, "P", time_index=source_time_index),
        dtype=np.float64,
    )
    if source_p.shape != reference_p.shape:
        return SourceEvolutionSummary(
            level=level,
            status="shape_mismatch",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            source_shape=tuple(source_p.shape),
            message=(
                f"source/start P shape {source_p.shape} != "
                f"end P shape {reference_p.shape}"
            ),
        )

    mask = _target_mask_for_p_level(reference_p, mask_2d, level)
    source_slice = source_p[level]
    candidate_slice = candidate_p[level]
    reference_slice = reference_p[level]
    candidate_evolution, candidate_count = _finite_diff_mean(
        source_slice,
        candidate_slice,
        mask,
    )
    reference_evolution, reference_count = _finite_diff_mean(
        source_slice,
        reference_slice,
        mask,
    )
    if candidate_evolution is None or reference_evolution is None:
        evolution_delta = None
    else:
        evolution_delta = candidate_evolution - reference_evolution

    valid_count = min(candidate_count, reference_count)
    return SourceEvolutionSummary(
        level=level,
        status="computed" if valid_count else "no_valid_samples",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        mean_source_p=_finite_mean(source_slice, mask),
        mean_candidate_p=_finite_mean(candidate_slice, mask),
        mean_reference_p=_finite_mean(reference_slice, mask),
        mean_candidate_minus_source_p=candidate_evolution,
        mean_reference_minus_source_p=reference_evolution,
        mean_candidate_evolution_minus_reference_evolution_p=evolution_delta,
        valid_count=valid_count,
        source_shape=tuple(source_p.shape),
        message=(
            None
            if valid_count
            else "no finite target-region source/end sample pairs are available"
        ),
    )


def _build_worst_level_attribution(
    reference_end: netCDF4.Dataset,
    candidate_end: netCDF4.Dataset,
    source_start: netCDF4.Dataset | None,
    reference_p: np.ndarray,
    candidate_p: np.ndarray,
    per_level: list[PerLevelPressureMetrics],
    worst_by_rmse: list[int],
    worst_by_bias: list[int],
    mask_2d: np.ndarray,
    *,
    time_index: int,
    source_time_index: int,
) -> list[WorstLevelAttribution]:
    rmse_ranks = _level_rank_map(worst_by_rmse)
    bias_ranks = _level_rank_map(worst_by_bias)
    selected_levels = sorted(set(worst_by_rmse) | set(worst_by_bias))
    by_level = {item.level: item for item in per_level}
    attribution: list[WorstLevelAttribution] = []
    for level in selected_levels:
        attribution.append(
            WorstLevelAttribution(
                level=level,
                rank_by_rmse=rmse_ranks.get(level),
                rank_by_abs_bias=bias_ranks.get(level),
                diagnostic_only=True,
                candidate_model_pass="not_applicable",
                pressure_metrics=by_level[level],
                companion_fields=[
                    _companion_summary(
                        reference_end,
                        candidate_end,
                        field,
                        level=level,
                        p_shape=reference_p.shape,
                        mask_2d=mask_2d,
                        time_index=time_index,
                    )
                    for field in COMPANION_FIELDS
                ],
                source_evolution=_source_evolution_summary(
                    reference_p,
                    candidate_p,
                    source_start,
                    level=level,
                    mask_2d=mask_2d,
                    source_time_index=source_time_index,
                ),
            )
        )
    return attribution


def audit_pressure_refresh_vertical_bias(
    reference_end_path: Path,
    candidate_end_path: Path,
    source_start_path: Path | None = None,
    *,
    time_index: int = -1,
    source_time_index: int = -1,
    worst_level_count: int = DEFAULT_WORST_LEVEL_COUNT,
) -> PressureRefreshVerticalBiasAudit:
    if worst_level_count <= 0:
        raise ValueError("worst-level-count must be positive")

    with netCDF4.Dataset(reference_end_path) as reference_end, netCDF4.Dataset(
        candidate_end_path
    ) as candidate_end:
        reference_p, candidate_p = _read_p_pair(
            reference_end,
            candidate_end,
            time_index=time_index,
        )
        metadata = _candidate_metadata(candidate_end)
        region = _infer_target_region(
            candidate_end,
            (int(reference_p.shape[-2]), int(reference_p.shape[-1])),
        )

        source_dataset: netCDF4.Dataset | None = None
        if source_start_path is not None:
            source_dataset = netCDF4.Dataset(source_start_path)
        try:
            region_available = region.get("status") in {
                "available",
                "available_count_mismatch",
            }
            if region_available and isinstance(region.get("mask"), np.ndarray):
                mask_2d = region["mask"]
                per_level = _per_level_metrics(reference_p, candidate_p, mask_2d)
                worst_by_rmse = _rank_levels_by_rmse(per_level, worst_level_count)
                worst_by_bias = _rank_levels_by_abs_bias(per_level, worst_level_count)
                attribution = _build_worst_level_attribution(
                    reference_end,
                    candidate_end,
                    source_dataset,
                    reference_p,
                    candidate_p,
                    per_level,
                    worst_by_rmse,
                    worst_by_bias,
                    mask_2d,
                    time_index=time_index,
                    source_time_index=source_time_index,
                )
                region_status = region["status"]
            else:
                per_level = []
                worst_by_rmse = []
                worst_by_bias = []
                attribution = []
                region_status = "not_available"
        finally:
            if source_dataset is not None:
                source_dataset.close()

    public_region = {key: value for key, value in region.items() if key != "mask"}
    computed_levels = sum(1 for item in per_level if item.status == "computed")
    max_rmse = max(
        (
            float(item.rmse)
            for item in per_level
            if _finite_metric_value(item.rmse) is not None
        ),
        default=None,
    )
    max_abs_bias = max(
        (
            abs(float(item.mean_diff))
            for item in per_level
            if _finite_metric_value(item.mean_diff) is not None
        ),
        default=None,
    )
    summary = {
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "p_variable": "P",
        "total_levels": int(reference_p.shape[-3]),
        "computed_levels": computed_levels,
        "region_status": region_status,
        "worst_level_count": worst_level_count,
        "worst_level_by_rmse": worst_by_rmse[0] if worst_by_rmse else None,
        "worst_level_by_abs_bias": worst_by_bias[0] if worst_by_bias else None,
        "max_target_region_rmse": max_rmse,
        "max_target_region_abs_bias": max_abs_bias,
        "source_start_available": source_start_path is not None,
    }
    status = "computed" if computed_levels else "not_available"
    if any(item.status != "computed" for item in per_level):
        status = "computed_with_flags" if computed_levels else "not_available"

    return PressureRefreshVerticalBiasAudit(
        reference_end=str(reference_end_path),
        candidate_end=str(candidate_end_path),
        source_start=None if source_start_path is None else str(source_start_path),
        status=status,
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        summary=summary,
        metadata=metadata,
        region=public_region,
        per_level_p=per_level,
        worst_levels_by_rmse=worst_by_rmse,
        worst_levels_by_abs_bias=worst_by_bias,
        worst_level_attribution=attribution,
        diagnosis={
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": (
                "diagnostic-only vertical pressure-bias attribution; reference-end "
                "comparisons explain pressure-refresh error structure and never "
                "create a validation gate pass"
            ),
        },
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


def report_to_dict(report: PressureRefreshVerticalBiasAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: PressureRefreshVerticalBiasAudit,
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
        "--candidate-end",
        "--candidate",
        dest="candidate_end",
        required=True,
        type=Path,
        help="TyWRF candidate wrfout at the same end time",
    )
    parser.add_argument(
        "--source-start",
        "--start",
        dest="source_start",
        type=Path,
        help="Optional source/start wrfout used for P evolution attribution",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for end-time Time-leading variables",
    )
    parser.add_argument(
        "--source-time-index",
        type=int,
        default=-1,
        help="Time index used for source/start Time-leading variables",
    )
    parser.add_argument(
        "--worst-level-count",
        type=int,
        default=DEFAULT_WORST_LEVEL_COUNT,
        help="Number of worst levels to rank by RMSE and absolute mean bias",
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
        report = audit_pressure_refresh_vertical_bias(
            args.reference_end,
            args.candidate_end,
            args.source_start,
            time_index=args.time_index,
            source_time_index=args.source_time_index,
            worst_level_count=args.worst_level_count,
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
