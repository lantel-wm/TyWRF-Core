#!/usr/bin/env python
"""Diagnostic-only audit for selected-field horizontal spatial alignment."""

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
class SpatialAlignmentVariableAudit:
    variable: str
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    reference_end_present: bool
    candidate_end_present: bool
    start_present: bool
    reference_end_shape: tuple[int, ...] | None = None
    candidate_end_shape: tuple[int, ...] | None = None
    start_shape: tuple[int, ...] | None = None
    max_shift: int | None = None
    end_state: dict[str, Any] | None = None
    evolution: dict[str, Any] | None = None
    message: str | None = None


@dataclass(frozen=True)
class SelectedFieldSpatialAlignmentAudit:
    reference_end: str
    candidate_end: str
    start: str | None
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    summary: dict[str, Any]
    metadata: dict[str, Any]
    region: dict[str, Any]
    variables: list[SpatialAlignmentVariableAudit]
    diagnosis: dict[str, Any]


def _shifted_candidate(candidate: np.ndarray, *, di: int, dj: int) -> np.ndarray:
    values = np.asarray(candidate, dtype=np.float64)
    if values.ndim < 2:
        if di == 0 and dj == 0:
            return values.copy()
        shifted = np.full(values.shape, np.nan, dtype=np.float64)
        return shifted

    ny, nx = int(values.shape[-2]), int(values.shape[-1])
    shifted = np.full(values.shape, np.nan, dtype=np.float64)
    src_i0 = max(0, -di)
    src_i1 = nx - max(0, di)
    dst_i0 = max(0, di)
    dst_i1 = nx - max(0, -di)
    src_j0 = max(0, -dj)
    src_j1 = ny - max(0, dj)
    dst_j0 = max(0, dj)
    dst_j1 = ny - max(0, -dj)
    if src_i0 >= src_i1 or src_j0 >= src_j1:
        return shifted

    shifted[..., dst_j0:dst_j1, dst_i0:dst_i1] = values[
        ..., src_j0:src_j1, src_i0:src_i1
    ]
    return shifted


def _best_shift_entry(entries: list[dict[str, Any]]) -> dict[str, Any] | None:
    available = [entry for entry in entries if entry["metrics"]["valid_count"] > 0]
    if not available:
        return None

    def sort_key(entry: dict[str, Any]) -> tuple[float, float, int, int, int]:
        metrics = entry["metrics"]
        normalized = metrics["normalized_rmse"]
        rmse = metrics["rmse"]
        if normalized is None or not math.isfinite(float(normalized)):
            normalized_key = math.inf
        else:
            normalized_key = float(normalized)
        if rmse is None or not math.isfinite(float(rmse)):
            rmse_key = math.inf
        else:
            rmse_key = float(rmse)
        shift = entry["shift"]
        distance = abs(int(shift["di"])) + abs(int(shift["dj"]))
        return normalized_key, rmse_key, distance, int(shift["dj"]), int(shift["di"])

    return min(available, key=sort_key)


def _improvement_summary(
    baseline: dict[str, Any],
    best: dict[str, Any] | None,
) -> dict[str, Any]:
    if best is None:
        return {
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "no shift had valid sample pairs",
        }

    base_metrics = baseline["metrics"]
    best_metrics = best["metrics"]

    def ratio(base: Any, shifted: Any) -> float | None:
        if base is None or shifted is None:
            return None
        base_value = float(base)
        shifted_value = float(shifted)
        if not math.isfinite(base_value) or not math.isfinite(shifted_value):
            return None
        if shifted_value == 0.0:
            return math.inf if base_value > 0.0 else 1.0
        return base_value / shifted_value

    def reduction(base: Any, shifted: Any) -> float | None:
        if base is None or shifted is None:
            return None
        base_value = float(base)
        shifted_value = float(shifted)
        if not math.isfinite(base_value) or not math.isfinite(shifted_value):
            return None
        if base_value == 0.0:
            return 0.0 if shifted_value == 0.0 else -math.inf
        return (base_value - shifted_value) / base_value

    norm_reduction = reduction(
        base_metrics["normalized_rmse"],
        best_metrics["normalized_rmse"],
    )
    return {
        "status": "available",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "best_shift": best["shift"],
        "baseline": base_metrics,
        "best": best_metrics,
        "rmse_improvement_ratio": ratio(base_metrics["rmse"], best_metrics["rmse"]),
        "normalized_rmse_improvement_ratio": ratio(
            base_metrics["normalized_rmse"],
            best_metrics["normalized_rmse"],
        ),
        "rmse_reduction_fraction": reduction(base_metrics["rmse"], best_metrics["rmse"]),
        "normalized_rmse_reduction_fraction": norm_reduction,
        "improved": None if norm_reduction is None else norm_reduction > 0.0,
        "message": (
            "integer shifts are diagnostic-only error attribution and must not be "
            "applied to candidate fields"
        ),
    }


def _scan_shift_alignment(
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    max_shift: int,
    selection_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    for dj in range(-max_shift, max_shift + 1):
        for di in range(-max_shift, max_shift + 1):
            shifted = _shifted_candidate(candidate, di=di, dj=dj)
            metrics = error_metrics(reference, shifted, selection_mask=selection_mask)
            entries.append(
                {
                    "shift": {"di": di, "dj": dj},
                    "diagnostic_only": True,
                    "candidate_model_pass": "not_applicable",
                    "metrics": _metrics_to_dict(metrics),
                }
            )

    baseline = next(
        entry
        for entry in entries
        if entry["shift"]["di"] == 0 and entry["shift"]["dj"] == 0
    )
    best = _best_shift_entry(entries)
    return {
        **_improvement_summary(baseline, best),
        "shift_count": len(entries),
        "tested_shifts": entries,
    }


def _scan_shifted_end_evolution_alignment(
    reference_evolution: np.ndarray,
    candidate_end: np.ndarray,
    candidate_start: np.ndarray,
    *,
    max_shift: int,
    selection_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    for dj in range(-max_shift, max_shift + 1):
        for di in range(-max_shift, max_shift + 1):
            shifted_end = _shifted_candidate(candidate_end, di=di, dj=dj)
            shifted_evolution = shifted_end - candidate_start
            metrics = error_metrics(
                reference_evolution,
                shifted_evolution,
                selection_mask=selection_mask,
            )
            entries.append(
                {
                    "shift": {"di": di, "dj": dj},
                    "diagnostic_only": True,
                    "candidate_model_pass": "not_applicable",
                    "metrics": _metrics_to_dict(metrics),
                }
            )

    baseline = next(
        entry
        for entry in entries
        if entry["shift"]["di"] == 0 and entry["shift"]["dj"] == 0
    )
    best = _best_shift_entry(entries)
    return {
        **_improvement_summary(baseline, best),
        "shift_count": len(entries),
        "tested_shifts": entries,
    }


def _region_alignment(
    variable: str,
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    max_shift: int,
) -> dict[str, Any]:
    if reference.ndim < 2:
        return {
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "variable has no horizontal dimensions",
        }

    mask_2d, mask_info = _field_region_mask(
        variable,
        (int(reference.shape[-2]), int(reference.shape[-1])),
        mass_shape,
        region,
    )
    if mask_2d is None:
        return {
            **mask_info,
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
        }

    target_mask = _region_mask_for_data(reference, mask_2d)
    overlap_mask = ~target_mask
    return {
        **mask_info,
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "target_region": _scan_shift_alignment(
            reference,
            candidate,
            max_shift=max_shift,
            selection_mask=target_mask,
        ),
        "overlap_region": _scan_shift_alignment(
            reference,
            candidate,
            max_shift=max_shift,
            selection_mask=overlap_mask,
        ),
    }


def _alignment_block(
    variable: str,
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    max_shift: int,
) -> dict[str, Any]:
    global_alignment = _scan_shift_alignment(
        reference,
        candidate,
        max_shift=max_shift,
    )
    return {
        "status": "computed",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "global": global_alignment,
        "region_split": _region_alignment(
            variable,
            reference,
            candidate,
            region=region,
            mass_shape=mass_shape,
            max_shift=max_shift,
        ),
    }


def _evolution_region_alignment(
    variable: str,
    reference_evolution: np.ndarray,
    candidate_end: np.ndarray,
    candidate_start: np.ndarray,
    *,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    max_shift: int,
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
    return {
        **mask_info,
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "target_region": _scan_shifted_end_evolution_alignment(
            reference_evolution,
            candidate_end,
            candidate_start,
            max_shift=max_shift,
            selection_mask=target_mask,
        ),
        "overlap_region": _scan_shifted_end_evolution_alignment(
            reference_evolution,
            candidate_end,
            candidate_start,
            max_shift=max_shift,
            selection_mask=overlap_mask,
        ),
    }


def _evolution_alignment_block(
    variable: str,
    reference_evolution: np.ndarray,
    candidate_end: np.ndarray,
    candidate_start: np.ndarray,
    *,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    max_shift: int,
) -> dict[str, Any]:
    return {
        "status": "computed",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "global": _scan_shifted_end_evolution_alignment(
            reference_evolution,
            candidate_end,
            candidate_start,
            max_shift=max_shift,
        ),
        "region_split": _evolution_region_alignment(
            variable,
            reference_evolution,
            candidate_end,
            candidate_start,
            region=region,
            mass_shape=mass_shape,
            max_shift=max_shift,
        ),
    }


def _missing_variable_audit(
    variable: str,
    *,
    reference_end_present: bool,
    candidate_end_present: bool,
    start_present: bool,
    max_shift: int,
) -> SpatialAlignmentVariableAudit:
    missing_owner = []
    if not reference_end_present:
        missing_owner.append("reference-end")
    if not candidate_end_present:
        missing_owner.append("candidate-end")
    if not start_present:
        missing_owner.append("start")
    return SpatialAlignmentVariableAudit(
        variable=variable,
        status="not_available",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        reference_end_present=reference_end_present,
        candidate_end_present=candidate_end_present,
        start_present=start_present,
        max_shift=max_shift,
        end_state={
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "variable missing",
        },
        evolution={
            "status": "not_available" if start_present else "not_requested",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "variable missing" if start_present else "start file not provided",
        },
        message=f"variable missing from {' and '.join(missing_owner)}",
    )


def compare_spatial_alignment_variable(
    reference_end: netCDF4.Dataset,
    candidate_end: netCDF4.Dataset,
    start: netCDF4.Dataset | None,
    variable: str,
    *,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    max_shift: int,
    time_index: int,
    start_time_index: int,
) -> SpatialAlignmentVariableAudit:
    reference_present = variable in reference_end.variables
    candidate_present = variable in candidate_end.variables
    start_present = start is not None and variable in start.variables
    if not reference_present or not candidate_present:
        return _missing_variable_audit(
            variable,
            reference_end_present=reference_present,
            candidate_end_present=candidate_present,
            start_present=start_present,
            max_shift=max_shift,
        )

    try:
        ref_end = _read_variable_data(reference_end, variable, time_index=time_index)
        cand_end = _read_variable_data(candidate_end, variable, time_index=time_index)
    except (IndexError, TypeError, ValueError) as exc:
        return SpatialAlignmentVariableAudit(
            variable=variable,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_end_present=True,
            candidate_end_present=True,
            start_present=start_present,
            max_shift=max_shift,
            message=f"end-state variable cannot be read: {exc}",
        )

    ref_shape = tuple(int(value) for value in ref_end.shape)
    cand_shape = tuple(int(value) for value in cand_end.shape)
    if ref_shape != cand_shape:
        return SpatialAlignmentVariableAudit(
            variable=variable,
            status="shape_mismatch",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_end_present=True,
            candidate_end_present=True,
            start_present=start_present,
            reference_end_shape=ref_shape,
            candidate_end_shape=cand_shape,
            max_shift=max_shift,
            message="reference-end and candidate-end shapes differ",
        )

    try:
        end_state = _alignment_block(
            variable,
            np.asarray(ref_end, dtype=np.float64),
            np.asarray(cand_end, dtype=np.float64),
            region=region,
            mass_shape=mass_shape,
            max_shift=max_shift,
        )
    except (TypeError, ValueError) as exc:
        return SpatialAlignmentVariableAudit(
            variable=variable,
            status="non_numeric",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            reference_end_present=True,
            candidate_end_present=True,
            start_present=start_present,
            reference_end_shape=ref_shape,
            candidate_end_shape=cand_shape,
            max_shift=max_shift,
            message=f"end-state variable cannot be compared numerically: {exc}",
        )

    evolution: dict[str, Any]
    start_shape: tuple[int, ...] | None = None
    if start is None:
        evolution = {
            "status": "not_requested",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "start file not provided",
        }
    elif not start_present:
        evolution = {
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": "variable missing from start file",
        }
    else:
        try:
            start_data = _read_variable_data(
                start,
                variable,
                time_index=start_time_index,
            )
        except (IndexError, TypeError, ValueError) as exc:
            evolution = {
                "status": "not_available",
                "diagnostic_only": True,
                "candidate_model_pass": "not_applicable",
                "message": f"start variable cannot be read: {exc}",
            }
        else:
            start_shape = tuple(int(value) for value in start_data.shape)
            if start_shape != ref_shape:
                evolution = {
                    "status": "shape_mismatch",
                    "diagnostic_only": True,
                    "candidate_model_pass": "not_applicable",
                    "reference_end_shape": ref_shape,
                    "start_shape": start_shape,
                    "message": "start shape differs from end-state shape",
                }
            else:
                reference_evolution = np.asarray(ref_end, dtype=np.float64) - np.asarray(
                    start_data,
                    dtype=np.float64,
                )
                candidate_start = np.asarray(start_data, dtype=np.float64)
                evolution = _evolution_alignment_block(
                    variable,
                    reference_evolution,
                    np.asarray(cand_end, dtype=np.float64),
                    candidate_start,
                    region=region,
                    mass_shape=mass_shape,
                    max_shift=max_shift,
                )
                evolution["message"] = (
                    "computed as shift(candidate_end) - candidate_start versus "
                    "reference_end - reference_start; shared start is diagnostic-only"
                )

    status = "computed"
    if evolution.get("status") in {"not_available", "shape_mismatch"}:
        status = "computed_with_flags"

    return SpatialAlignmentVariableAudit(
        variable=variable,
        status=status,
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        reference_end_present=True,
        candidate_end_present=True,
        start_present=start_present,
        reference_end_shape=ref_shape,
        candidate_end_shape=cand_shape,
        start_shape=start_shape,
        max_shift=max_shift,
        end_state=end_state,
        evolution=evolution,
    )


def audit_selected_field_spatial_alignment(
    reference_end_path: Path,
    candidate_end_path: Path,
    start_path: Path | None = None,
    *,
    variables: Iterable[str] = DEFAULT_STATE_VARIABLES,
    max_shift: int = 2,
    time_index: int = -1,
    start_time_index: int = -1,
) -> SelectedFieldSpatialAlignmentAudit:
    if max_shift < 0:
        raise ValueError("--max-shift must be non-negative")

    variables = tuple(variables)
    with netCDF4.Dataset(reference_end_path) as reference_end, netCDF4.Dataset(
        candidate_end_path
    ) as candidate_end:
        start_context = netCDF4.Dataset(start_path) if start_path else None
        try:
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
                compare_spatial_alignment_variable(
                    reference_end,
                    candidate_end,
                    start_context,
                    variable,
                    region=region,
                    mass_shape=mass_shape,
                    max_shift=max_shift,
                    time_index=time_index,
                    start_time_index=start_time_index,
                )
                for variable in variables
            ]
        finally:
            if start_context is not None:
                start_context.close()

    counts = Counter(audit.status for audit in audits)
    end_best_shifts = {}
    evolution_best_shifts = {}
    end_improved = {}
    evolution_improved = {}
    for audit in audits:
        end_global = ((audit.end_state or {}).get("global") or {})
        evolution_global = ((audit.evolution or {}).get("global") or {})
        end_best_shifts[audit.variable] = end_global.get("best_shift")
        evolution_best_shifts[audit.variable] = evolution_global.get("best_shift")
        end_improved[audit.variable] = end_global.get("improved")
        evolution_improved[audit.variable] = evolution_global.get("improved")

    status = (
        "computed_with_flags"
        if (
            counts.get("not_available")
            or counts.get("shape_mismatch")
            or counts.get("computed_with_flags")
        )
        else "computed"
    )
    summary = {
        **dict(sorted(counts.items())),
        "total": len(audits),
        "available": sum(1 for audit in audits if audit.status != "not_available"),
        "not_available": counts.get("not_available", 0),
        "shape_mismatches": counts.get("shape_mismatch", 0),
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "max_shift": max_shift,
        "shift_window": {
            "di": [-max_shift, max_shift],
            "dj": [-max_shift, max_shift],
        },
        "start_provided": start_path is not None,
        "end_state_best_shift": end_best_shifts,
        "evolution_best_shift": evolution_best_shifts,
        "end_state_improved_by_shift": end_improved,
        "evolution_improved_by_shift": evolution_improved,
    }
    return SelectedFieldSpatialAlignmentAudit(
        reference_end=str(reference_end_path),
        candidate_end=str(candidate_end_path),
        start=None if start_path is None else str(start_path),
        status=status,
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        summary=summary,
        metadata=metadata,
        region=region,
        variables=audits,
        diagnosis={
            "end_state_best_shift": end_best_shifts,
            "evolution_best_shift": evolution_best_shifts,
            "end_state_improved_by_shift": end_improved,
            "evolution_improved_by_shift": evolution_improved,
            "candidate_model_pass": "not_applicable",
            "diagnostic_only": True,
            "message": (
                "diagnostic-only spatial alignment audit; reference-end output is "
                "used only for error attribution. Best shifts must not be applied "
                "to candidate fields and do not generate validation gate passes."
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


def report_to_dict(report: SelectedFieldSpatialAlignmentAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: SelectedFieldSpatialAlignmentAudit,
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
        help="Optional shared reference/candidate start-state wrfout file",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_STATE_VARIABLES),
        help="State variables to audit",
    )
    parser.add_argument(
        "--max-shift",
        type=int,
        default=2,
        help="Maximum absolute integer i/j shift to test",
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
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_selected_field_spatial_alignment(
            args.reference_end,
            args.candidate_end,
            args.start,
            variables=args.variables,
            max_shift=args.max_shift,
            time_index=args.time_index,
            start_time_index=args.start_time_index,
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
