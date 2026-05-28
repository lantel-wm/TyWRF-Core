#!/usr/bin/env python
"""Audit selected field evolution errors against WRF reference evolution."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

import netCDF4
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if __package__ in {None, ""} and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.audit_selected_field_state import (
    DEFAULT_STATE_VARIABLES,
    ErrorMetrics,
    _candidate_metadata,
    _field_region_mask,
    _infer_target_region,
    _mass_horizontal_shape,
    _metrics_to_dict,
    _read_variable_data,
    _region_mask_for_data,
    error_metrics,
)


@dataclass(frozen=True)
class EvolutionVariableAudit:
    variable: str
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    reference_end_present: bool
    candidate_end_present: bool
    reference_start_present: bool
    candidate_start_present: bool
    reference_end_shape: tuple[int, ...] | None = None
    candidate_end_shape: tuple[int, ...] | None = None
    reference_start_shape: tuple[int, ...] | None = None
    candidate_start_shape: tuple[int, ...] | None = None
    rmse: float | None = None
    normalized_rmse: float | None = None
    max_abs_error: float | None = None
    valid_count: int | None = None
    total_count: int | None = None
    reference_finite_count: int | None = None
    candidate_finite_count: int | None = None
    evolution_amplitude: dict[str, Any] | None = None
    region_split: dict[str, Any] | None = None
    level_summary: dict[str, Any] | None = None
    message: str | None = None


@dataclass(frozen=True)
class SelectedFieldEvolutionAudit:
    reference_end: str
    candidate_end: str
    reference_start: str
    candidate_start: str
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    summary: dict[str, Any]
    metadata: dict[str, Any]
    region: dict[str, Any]
    variables: list[EvolutionVariableAudit]
    diagnosis: dict[str, Any]


def _finite_pair(
    reference_evolution: np.ndarray,
    candidate_evolution: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reference = np.ma.masked_invalid(np.ma.asarray(reference_evolution, dtype=np.float64))
    candidate = np.ma.masked_invalid(np.ma.asarray(candidate_evolution, dtype=np.float64))
    ref_values = np.asarray(reference.filled(np.nan), dtype=np.float64)
    cand_values = np.asarray(candidate.filled(np.nan), dtype=np.float64)
    valid = (
        np.isfinite(ref_values)
        & np.isfinite(cand_values)
        & ~np.ma.getmaskarray(reference)
        & ~np.ma.getmaskarray(candidate)
    )
    return ref_values, cand_values, valid


def evolution_amplitude(
    reference_evolution: np.ndarray,
    candidate_evolution: np.ndarray,
    *,
    selection_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    ref_values, cand_values, valid = _finite_pair(reference_evolution, candidate_evolution)
    if ref_values.shape != cand_values.shape:
        raise ValueError(f"shape mismatch: {ref_values.shape} != {cand_values.shape}")
    if selection_mask is not None:
        selected = np.asarray(selection_mask, dtype=bool)
        if selected.shape != ref_values.shape:
            selected = np.broadcast_to(selected, ref_values.shape)
        valid = valid & selected

    valid_count = int(np.count_nonzero(valid))
    if valid_count == 0:
        return {
            "status": "no_valid_samples",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "reference_evolution_rms": None,
            "candidate_evolution_rms": None,
            "amplitude_ratio": None,
            "capture_fraction": None,
            "valid_count": 0,
        }

    ref_sample = ref_values[valid]
    cand_sample = cand_values[valid]
    ref_sse = float(np.sum(ref_sample * ref_sample))
    cand_sse = float(np.sum(cand_sample * cand_sample))
    ref_rms = float(math.sqrt(ref_sse / valid_count))
    cand_rms = float(math.sqrt(cand_sse / valid_count))
    if ref_rms == 0.0:
        amplitude_ratio = None
    else:
        amplitude_ratio = cand_rms / ref_rms
    if ref_sse == 0.0:
        capture_fraction = None
    else:
        capture_fraction = float(np.sum(cand_sample * ref_sample) / ref_sse)

    return {
        "status": "available",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "reference_evolution_rms": ref_rms,
        "candidate_evolution_rms": cand_rms,
        "amplitude_ratio": amplitude_ratio,
        "capture_fraction": capture_fraction,
        "valid_count": valid_count,
    }


def _evolution_region_split(
    variable: str,
    reference_evolution: np.ndarray,
    candidate_evolution: np.ndarray,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    global_metrics: ErrorMetrics,
) -> dict[str, Any]:
    if reference_evolution.ndim < 2:
        return {
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "variable has no horizontal dimensions",
        }

    mask_2d, mask_info = _field_region_mask(
        variable,
        (int(reference_evolution.shape[-2]), int(reference_evolution.shape[-1])),
        mass_shape,
        region,
    )
    if mask_2d is None:
        return {
            **mask_info,
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
        }

    target_mask = _region_mask_for_data(reference_evolution, mask_2d)
    overlap_mask = ~target_mask
    target_metrics = error_metrics(
        reference_evolution,
        candidate_evolution,
        selection_mask=target_mask,
    )
    overlap_metrics = error_metrics(
        reference_evolution,
        candidate_evolution,
        selection_mask=overlap_mask,
    )
    target_sse = target_metrics.sum_squared_error
    global_sse = global_metrics.sum_squared_error
    if math.isfinite(target_sse) and math.isfinite(global_sse) and global_sse > 0.0:
        target_error_fraction = target_sse / global_sse
    else:
        target_error_fraction = None

    return {
        **mask_info,
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "target_region": {
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            **_metrics_to_dict(target_metrics),
            "evolution_amplitude": evolution_amplitude(
                reference_evolution,
                candidate_evolution,
                selection_mask=target_mask,
            ),
        },
        "overlap_region": {
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            **_metrics_to_dict(overlap_metrics),
            "evolution_amplitude": evolution_amplitude(
                reference_evolution,
                candidate_evolution,
                selection_mask=overlap_mask,
            ),
        },
        "target_error_fraction": target_error_fraction,
        "target_region_dominates_global_error": (
            None if target_error_fraction is None else target_error_fraction > 0.5
        ),
    }


def _level_summaries(
    variable: str,
    reference_evolution: np.ndarray,
    candidate_evolution: np.ndarray,
    *,
    enabled: bool,
    limit: int,
) -> dict[str, Any]:
    if not enabled:
        return {
            "status": "not_requested",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
        }
    if variable not in {"U", "V"}:
        return {
            "status": "not_applicable",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "per-level summaries are limited to U/V",
        }
    if reference_evolution.ndim != 3:
        return {
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "expected a 3D field after selecting the Time index",
        }

    levels = []
    for level in range(reference_evolution.shape[0]):
        metrics = error_metrics(
            reference_evolution[level, ...],
            candidate_evolution[level, ...],
        )
        levels.append(
            {
                "level": level,
                **_metrics_to_dict(metrics),
                "evolution_amplitude": evolution_amplitude(
                    reference_evolution[level, ...],
                    candidate_evolution[level, ...],
                ),
            }
        )

    def sort_key(item: dict[str, Any]) -> float:
        value = item["normalized_rmse"]
        if value is None or math.isnan(value):
            return -math.inf
        return float(value)

    worst = sorted(levels, key=sort_key, reverse=True)[: max(0, limit)]
    return {
        "status": "available",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "vertical_axis": 0,
        "level_count": len(levels),
        "worst_levels": worst,
    }


def _missing_variable_audit(
    variable: str,
    *,
    reference_end_present: bool,
    candidate_end_present: bool,
    reference_start_present: bool,
    candidate_start_present: bool,
) -> EvolutionVariableAudit:
    missing_owner = []
    if not reference_end_present:
        missing_owner.append("reference-end")
    if not candidate_end_present:
        missing_owner.append("candidate-end")
    if not reference_start_present:
        missing_owner.append("reference-start")
    if not candidate_start_present:
        missing_owner.append("candidate-start")
    return EvolutionVariableAudit(
        variable=variable,
        status="not_available",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        reference_end_present=reference_end_present,
        candidate_end_present=candidate_end_present,
        reference_start_present=reference_start_present,
        candidate_start_present=candidate_start_present,
        region_split={
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "variable missing",
        },
        level_summary={
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "variable missing",
        },
        message=f"variable missing from {' and '.join(missing_owner)}",
    )


def compare_evolution_variable(
    reference_end: netCDF4.Dataset,
    candidate_end: netCDF4.Dataset,
    reference_start: netCDF4.Dataset,
    candidate_start: netCDF4.Dataset,
    variable: str,
    *,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    time_index: int,
    start_time_index: int,
    include_level_summary: bool,
    level_summary_count: int,
) -> EvolutionVariableAudit:
    reference_end_present = variable in reference_end.variables
    candidate_end_present = variable in candidate_end.variables
    reference_start_present = variable in reference_start.variables
    candidate_start_present = variable in candidate_start.variables
    if not all(
        (
            reference_end_present,
            candidate_end_present,
            reference_start_present,
            candidate_start_present,
        )
    ):
        return _missing_variable_audit(
            variable,
            reference_end_present=reference_end_present,
            candidate_end_present=candidate_end_present,
            reference_start_present=reference_start_present,
            candidate_start_present=candidate_start_present,
        )

    try:
        ref_end = _read_variable_data(reference_end, variable, time_index=time_index)
        cand_end = _read_variable_data(candidate_end, variable, time_index=time_index)
        ref_start = _read_variable_data(
            reference_start,
            variable,
            time_index=start_time_index,
        )
        cand_start = _read_variable_data(
            candidate_start,
            variable,
            time_index=start_time_index,
        )
    except (IndexError, TypeError, ValueError) as exc:
        return EvolutionVariableAudit(
            variable=variable,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_end_present=True,
            candidate_end_present=True,
            reference_start_present=True,
            candidate_start_present=True,
            message=f"variable cannot be read: {exc}",
        )

    shapes = {
        "reference_end_shape": tuple(int(value) for value in ref_end.shape),
        "candidate_end_shape": tuple(int(value) for value in cand_end.shape),
        "reference_start_shape": tuple(int(value) for value in ref_start.shape),
        "candidate_start_shape": tuple(int(value) for value in cand_start.shape),
    }
    unique_shapes = set(shapes.values())
    if len(unique_shapes) != 1:
        return EvolutionVariableAudit(
            variable=variable,
            status="shape_mismatch",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_end_present=True,
            candidate_end_present=True,
            reference_start_present=True,
            candidate_start_present=True,
            **shapes,
            message="end/start shapes differ for evolution comparison",
        )

    try:
        reference_evolution = np.asarray(ref_end, dtype=np.float64) - np.asarray(
            ref_start,
            dtype=np.float64,
        )
        candidate_evolution = np.asarray(cand_end, dtype=np.float64) - np.asarray(
            cand_start,
            dtype=np.float64,
        )
        metrics = error_metrics(reference_evolution, candidate_evolution)
        amplitude = evolution_amplitude(reference_evolution, candidate_evolution)
    except (TypeError, ValueError) as exc:
        return EvolutionVariableAudit(
            variable=variable,
            status="non_numeric",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_end_present=True,
            candidate_end_present=True,
            reference_start_present=True,
            candidate_start_present=True,
            **shapes,
            message=f"variable cannot be compared numerically: {exc}",
        )

    if metrics.valid_count == 0:
        status = "no_valid_samples"
        message = "no finite evolution sample pairs are available"
    else:
        status = "computed"
        message = None

    split = _evolution_region_split(
        variable,
        reference_evolution,
        candidate_evolution,
        region,
        mass_shape,
        metrics,
    )
    level_summary = _level_summaries(
        variable,
        reference_evolution,
        candidate_evolution,
        enabled=include_level_summary,
        limit=level_summary_count,
    )

    return EvolutionVariableAudit(
        variable=variable,
        status=status,
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        reference_end_present=True,
        candidate_end_present=True,
        reference_start_present=True,
        candidate_start_present=True,
        **shapes,
        rmse=metrics.rmse,
        normalized_rmse=metrics.normalized_rmse,
        max_abs_error=metrics.max_abs_error,
        valid_count=metrics.valid_count,
        total_count=metrics.total_count,
        reference_finite_count=metrics.reference_finite_count,
        candidate_finite_count=metrics.candidate_finite_count,
        evolution_amplitude=amplitude,
        region_split=split,
        level_summary=level_summary,
        message=message,
    )


def audit_selected_field_evolution(
    reference_end_path: Path,
    candidate_end_path: Path,
    start_path: Path | None = None,
    *,
    reference_start_path: Path | None = None,
    candidate_start_path: Path | None = None,
    variables: Iterable[str] = DEFAULT_STATE_VARIABLES,
    time_index: int = -1,
    start_time_index: int = -1,
    include_level_summary: bool = False,
    level_summary_count: int = 5,
) -> SelectedFieldEvolutionAudit:
    if start_path is not None:
        if reference_start_path is None:
            reference_start_path = start_path
        if candidate_start_path is None:
            candidate_start_path = start_path
    if reference_start_path is None or candidate_start_path is None:
        raise ValueError(
            "provide a shared start file or both --reference-start and --candidate-start"
        )

    variables = tuple(variables)
    with netCDF4.Dataset(reference_end_path) as reference_end, netCDF4.Dataset(
        candidate_end_path
    ) as candidate_end, netCDF4.Dataset(
        reference_start_path
    ) as reference_start, netCDF4.Dataset(candidate_start_path) as candidate_start:
        metadata = _candidate_metadata(candidate_end)
        region = {
            **_infer_target_region(candidate_end),
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
        }
        mass_shape = _mass_horizontal_shape(
            reference_end,
            candidate_end,
            time_index=time_index,
        )
        audits = [
            compare_evolution_variable(
                reference_end,
                candidate_end,
                reference_start,
                candidate_start,
                variable,
                region=region,
                mass_shape=mass_shape,
                time_index=time_index,
                start_time_index=start_time_index,
                include_level_summary=include_level_summary,
                level_summary_count=level_summary_count,
            )
            for variable in variables
        ]

    counts = Counter(audit.status for audit in audits)
    target_fractions = {}
    target_dominates = {}
    amplitude_ratios = {}
    capture_fractions = {}
    for audit in audits:
        split = audit.region_split or {}
        target_fractions[audit.variable] = (
            split.get("target_error_fraction")
            if split.get("status") == "available"
            else None
        )
        target_dominates[audit.variable] = (
            split.get("target_region_dominates_global_error")
            if split.get("status") == "available"
            else None
        )
        amplitude = audit.evolution_amplitude or {}
        amplitude_ratios[audit.variable] = amplitude.get("amplitude_ratio")
        capture_fractions[audit.variable] = amplitude.get("capture_fraction")

    numeric_flags = sum(
        1
        for audit in audits
        if audit.status in {"no_valid_samples", "non_numeric", "shape_mismatch"}
    )
    status = "computed_with_flags" if numeric_flags or counts.get("not_available") else "computed"
    summary = {
        **dict(sorted(counts.items())),
        "total": len(audits),
        "available": sum(1 for audit in audits if audit.status != "not_available"),
        "not_available": counts.get("not_available", 0),
        "numeric_failures": numeric_flags,
        "shape_mismatches": counts.get("shape_mismatch", 0),
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "shared_start_file": reference_start_path == candidate_start_path,
        "target_region_error_fraction": target_fractions,
        "target_region_dominates_error": target_dominates,
        "amplitude_ratio": amplitude_ratios,
        "capture_fraction": capture_fractions,
    }
    return SelectedFieldEvolutionAudit(
        reference_end=str(reference_end_path),
        candidate_end=str(candidate_end_path),
        reference_start=str(reference_start_path),
        candidate_start=str(candidate_start_path),
        status=status,
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        summary=summary,
        metadata=metadata,
        region=region,
        variables=audits,
        diagnosis={
            "target_region_error_fraction": target_fractions,
            "target_region_dominates_error": target_dominates,
            "amplitude_ratio": amplitude_ratios,
            "capture_fraction": capture_fractions,
            "candidate_model_pass": "not_applicable",
            "message": (
                "diagnostic-only selected-field evolution audit; reference-end "
                "comparison is error attribution, not validation pass generation"
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


def report_to_dict(report: SelectedFieldEvolutionAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: SelectedFieldEvolutionAudit,
    *,
    pretty: bool = False,
) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_end", type=Path, help="WRF reference-end wrfout file")
    parser.add_argument("candidate_end", type=Path, help="TyWRF candidate-end wrfout file")
    parser.add_argument(
        "start",
        type=Path,
        nargs="?",
        help="Shared candidate/reference start-state wrfout file",
    )
    parser.add_argument(
        "--reference-start",
        type=Path,
        help="Reference start-state wrfout file when starts are separate",
    )
    parser.add_argument(
        "--candidate-start",
        type=Path,
        help="Candidate start-state wrfout file when starts are separate",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_STATE_VARIABLES),
        help="State variables to audit",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for reference-end and candidate-end Time-leading variables",
    )
    parser.add_argument(
        "--start-time-index",
        type=int,
        default=-1,
        help="Time index used for start-state Time-leading variables",
    )
    parser.add_argument(
        "--level-summary",
        action="store_true",
        help="Include worst-level summaries for 3D U/V fields",
    )
    parser.add_argument(
        "--level-count",
        type=int,
        default=5,
        help="Number of worst vertical levels to report when --level-summary is set",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_selected_field_evolution(
            args.reference_end,
            args.candidate_end,
            args.start,
            reference_start_path=args.reference_start,
            candidate_start_path=args.candidate_start,
            variables=args.variables,
            time_index=args.time_index,
            start_time_index=args.start_time_index,
            include_level_summary=args.level_summary,
            level_summary_count=args.level_count,
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
