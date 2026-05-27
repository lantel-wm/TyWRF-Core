#!/usr/bin/env python
"""Audit pressure-refresh candidate pressure errors against a WRF reference file."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable

import netCDF4
import numpy as np


PRESSURE_REFRESH_VARIABLES = ("P", "PB", "MUB", "PHB", "SLP")
P_LIKE_3D_VARIABLES = {"P", "PB", "PHB"}
DEFAULT_THRESHOLDS = {
    "P": 0.05,
    "PB": 0.05,
    "MUB": 0.05,
    "PHB": 0.05,
    "SLP": 0.10,
}
BOOL_TRUE_VALUES = {"1", "true", "t", "yes", "y", "on"}
BOOL_FALSE_VALUES = {"0", "false", "f", "no", "n", "off"}
BLOCKING_CANDIDATE_KIND_TOKENS = ("diagnostic", "closure", "remap", "oracle")
METADATA_ATTRS = (
    "TYWRF_DIAGNOSTIC_ONLY",
    "TYWRF_GATE_CANDIDATE",
    "TYWRF_INTEGRATOR_OUTPUT",
    "TYWRF_VALIDATION_GATE_ONLY",
    "TYWRF_CANDIDATE_KIND",
    "TYWRF_FROM_PARENT_START",
    "TYWRF_TO_PARENT_START",
    "TYWRF_PARENT_GRID_RATIO",
    "TYWRF_PRESSURE_REFRESH_APPLIED",
    "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS",
    "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY",
    "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY",
    "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT",
    "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT",
    "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT",
    "TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT",
    "TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT",
    "TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS",
    "TYWRF_PRESSURE_REFRESH_TOUCHED_HALO_CELLS",
    "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT",
    "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO",
    "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED",
)


@dataclass(frozen=True)
class ErrorMetrics:
    rmse: float
    normalized_rmse: float
    max_abs_error: float
    valid_count: int
    total_count: int
    reference_finite_count: int
    candidate_finite_count: int
    sum_squared_error: float


@dataclass(frozen=True)
class VariableAudit:
    variable: str
    status: str
    reference_present: bool
    candidate_present: bool
    reference_shape: tuple[int, ...] | None = None
    candidate_shape: tuple[int, ...] | None = None
    threshold: float | None = None
    rmse: float | None = None
    normalized_rmse: float | None = None
    max_abs_error: float | None = None
    valid_count: int | None = None
    total_count: int | None = None
    reference_finite_count: int | None = None
    candidate_finite_count: int | None = None
    region_split: dict[str, Any] | None = None
    message: str | None = None


@dataclass(frozen=True)
class PressureRefreshCandidateAudit:
    reference: str
    candidate: str
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    summary: dict[str, Any]
    metadata: dict[str, Any]
    region: dict[str, Any]
    variables: list[VariableAudit]
    diagnosis: dict[str, Any]


def _netcdf_attr_to_json(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, (list, tuple)):
        return [_netcdf_attr_to_json(item) for item in value]
    if hasattr(value, "tolist"):
        return _netcdf_attr_to_json(value.tolist())
    if hasattr(value, "item"):
        return _netcdf_attr_to_json(value.item())
    return str(value).strip()


def _read_attrs(dataset: netCDF4.Dataset) -> dict[str, Any]:
    attrs = {}
    for name in dataset.ncattrs():
        if name in METADATA_ATTRS or name.startswith("TYWRF_PRESSURE_REFRESH_"):
            attrs[name] = _netcdf_attr_to_json(dataset.getncattr(name))
    return attrs


def _coerce_bool(value: Any) -> bool | None:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in BOOL_TRUE_VALUES:
        return True
    if text in BOOL_FALSE_VALUES:
        return False
    return None


def _coerce_int(value: Any) -> int | None:
    if isinstance(value, (bool, np.bool_)) or value is None:
        return None
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if numeric.is_integer():
            return int(numeric)
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        try:
            numeric = float(text)
        except ValueError:
            return None
        if numeric.is_integer():
            return int(numeric)
    return None


def _coerce_pair(value: Any) -> tuple[int, int] | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)) and len(value) == 2:
        first = _coerce_int(value[0])
        second = _coerce_int(value[1])
        if first is not None and second is not None:
            return first, second
        return None

    parts = [part for part in re.split(r"[\s,]+", str(value).strip()) if part]
    if len(parts) != 2:
        return None
    first = _coerce_int(parts[0])
    second = _coerce_int(parts[1])
    if first is None or second is None:
        return None
    return first, second


def _candidate_metadata(candidate: netCDF4.Dataset) -> dict[str, Any]:
    attrs = _read_attrs(candidate)
    gate_candidate = _coerce_bool(attrs.get("TYWRF_GATE_CANDIDATE"))
    candidate_diagnostic_only = _coerce_bool(attrs.get("TYWRF_DIAGNOSTIC_ONLY"))
    integrator_output = _coerce_bool(attrs.get("TYWRF_INTEGRATOR_OUTPUT"))
    validation_gate_only = _coerce_bool(attrs.get("TYWRF_VALIDATION_GATE_ONLY"))
    candidate_kind = attrs.get("TYWRF_CANDIDATE_KIND")
    pressure_refresh_applied = _coerce_bool(
        attrs.get("TYWRF_PRESSURE_REFRESH_APPLIED")
    )
    integration_status = attrs.get("TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS")
    pressure_refresh_experimental_apply = _coerce_bool(
        attrs.get("TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY")
    )
    experimental_pressure_refresh_apply = _coerce_bool(
        attrs.get("TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY")
    )

    blockers: list[str] = []
    if candidate_diagnostic_only is True:
        blockers.append("TYWRF_DIAGNOSTIC_ONLY=true")
    if gate_candidate is False:
        blockers.append("TYWRF_GATE_CANDIDATE=false")
    elif gate_candidate is not True:
        blockers.append("TYWRF_GATE_CANDIDATE is not true")
    if integrator_output is False:
        blockers.append("TYWRF_INTEGRATOR_OUTPUT=false")
    elif integrator_output is not True:
        blockers.append("TYWRF_INTEGRATOR_OUTPUT is not true")
    if validation_gate_only is True:
        blockers.append("TYWRF_VALIDATION_GATE_ONLY=true")

    kind = str(candidate_kind).lower() if candidate_kind else ""
    if kind and any(token in kind for token in BLOCKING_CANDIDATE_KIND_TOKENS):
        blockers.append(f"TYWRF_CANDIDATE_KIND={candidate_kind}")
    if pressure_refresh_experimental_apply is True:
        blockers.append("TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY=true")
    if experimental_pressure_refresh_apply is True:
        blockers.append("TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY=true")
    if str(integration_status).strip().lower() == "experimental_apply_test_only":
        blockers.append(
            "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS=experimental_apply_test_only"
        )

    hidden_seam = (
        pressure_refresh_experimental_apply is True
        or experimental_pressure_refresh_apply is True
        or str(integration_status).strip().lower() == "experimental_apply_test_only"
    )
    gate_eligible = not blockers
    normal_candidate = (
        pressure_refresh_applied is True
        and integration_status == "applied_to_candidate"
        and gate_eligible
        and candidate_diagnostic_only is False
        and gate_candidate is True
        and integrator_output is True
        and validation_gate_only is not True
        and not hidden_seam
    )

    if pressure_refresh_applied is True and hidden_seam:
        disposition_status = "hidden_seam_diagnostic_evidence_only"
    elif pressure_refresh_applied is True and normal_candidate:
        disposition_status = "normal_candidate_metadata_only"
    elif pressure_refresh_applied is True and not gate_eligible:
        disposition_status = "diagnostic_helper_evidence_only"
    elif pressure_refresh_applied is True:
        disposition_status = "candidate_metadata_only"
    elif pressure_refresh_applied is False:
        disposition_status = "not_applied"
    else:
        disposition_status = "not_reported"

    disposition = {
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "gate_candidate": gate_candidate,
        "candidate_diagnostic_only": candidate_diagnostic_only,
        "integrator_output": integrator_output,
        "validation_gate_only": validation_gate_only,
        "candidate_kind": candidate_kind,
        "candidate_gate_eligible": gate_eligible,
        "candidate_gate_blockers": blockers,
        "pressure_refresh_applied": pressure_refresh_applied,
        "pressure_refresh_integration_status": integration_status,
        "pressure_refresh_experimental_apply": pressure_refresh_experimental_apply,
        "experimental_pressure_refresh_apply": experimental_pressure_refresh_apply,
        "normal_candidate_pressure_refresh": normal_candidate,
        "hidden_seam_pressure_refresh": hidden_seam,
        "pressure_refresh_disposition": {
            "status": disposition_status,
            "candidate_model_pass": "not_applicable",
            "message": (
                "pressure-refresh metadata is diagnostic evidence only; "
                "this audit never creates a model validation pass"
            ),
        },
        "from_parent_start": _pair_to_dict(
            _coerce_pair(attrs.get("TYWRF_FROM_PARENT_START"))
        ),
        "to_parent_start": _pair_to_dict(
            _coerce_pair(attrs.get("TYWRF_TO_PARENT_START"))
        ),
        "parent_grid_ratio": _coerce_int(attrs.get("TYWRF_PARENT_GRID_RATIO")),
        "pressure_refresh_target_column_count": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT")
        ),
        "pressure_refresh_refreshed_column_count": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT")
        ),
        "pressure_refresh_refreshed_point_count": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT")
        ),
        "pressure_refresh_skipped_point_count": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT")
        ),
        "pressure_refresh_invalid_point_count": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT")
        ),
        "pressure_refresh_touched_overlap_cells": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS")
        ),
        "pressure_refresh_touched_halo_cells": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_TOUCHED_HALO_CELLS")
        ),
        "pressure_refresh_refreshed_p_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS")
        ),
        "pressure_refresh_changed_p_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS")
        ),
        "pressure_refresh_changed_pb_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS")
        ),
        "pressure_refresh_changed_mub_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS")
        ),
        "pressure_refresh_changed_phb_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS")
        ),
        "pressure_refresh_changed_p_matches_refreshed_point_count": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT")
        ),
        "pressure_refresh_invalid_and_skipped_points_zero": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO")
        ),
        "pressure_refresh_overlap_halo_untouched": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED")
        ),
    }
    return {
        "status": "available",
        "attrs": attrs,
        "present_attrs": sorted(attrs),
        "missing_attrs": [name for name in METADATA_ATTRS if name not in attrs],
        "disposition": disposition,
    }


def _pair_to_dict(value: tuple[int, int] | None) -> dict[str, int] | None:
    if value is None:
        return None
    return {"i": value[0], "j": value[1]}


def _read_variable_data(
    dataset: netCDF4.Dataset,
    name: str,
    *,
    time_index: int,
) -> Any:
    variable = dataset.variables[name]
    if variable.ndim > 0 and variable.dimensions[0] == "Time":
        length = int(variable.shape[0])
        index = time_index if time_index >= 0 else length + time_index
        if index < 0 or index >= length:
            raise IndexError(
                f"{name} has Time length {length}, cannot read time_index={time_index}"
            )
        return variable[index, ...]
    return variable[:]


def _finite_array(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    masked = np.ma.masked_invalid(np.ma.asarray(data, dtype=np.float64))
    values = np.asarray(masked.filled(np.nan), dtype=np.float64)
    finite = np.isfinite(values) & ~np.ma.getmaskarray(masked)
    return values, finite


def error_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    selection_mask: np.ndarray | None = None,
) -> ErrorMetrics:
    ref_values, ref_finite = _finite_array(reference)
    cand_values, cand_finite = _finite_array(candidate)
    if ref_values.shape != cand_values.shape:
        raise ValueError(f"shape mismatch: {ref_values.shape} != {cand_values.shape}")

    if selection_mask is None:
        selected = np.ones(ref_values.shape, dtype=bool)
    else:
        selected = np.asarray(selection_mask, dtype=bool)
        if selected.shape != ref_values.shape:
            selected = np.broadcast_to(selected, ref_values.shape)

    pair_mask = selected & ref_finite & cand_finite
    ref_count = int(np.count_nonzero(selected & ref_finite))
    cand_count = int(np.count_nonzero(selected & cand_finite))
    total_count = int(np.count_nonzero(selected))
    valid_count = int(np.count_nonzero(pair_mask))
    if valid_count == 0:
        return ErrorMetrics(
            rmse=math.nan,
            normalized_rmse=math.nan,
            max_abs_error=math.nan,
            valid_count=0,
            total_count=total_count,
            reference_finite_count=ref_count,
            candidate_finite_count=cand_count,
            sum_squared_error=math.nan,
        )

    ref_sample = ref_values[pair_mask]
    diff = cand_values[pair_mask] - ref_sample
    sum_squared_error = float(np.sum(diff * diff))
    rmse = float(math.sqrt(sum_squared_error / valid_count))
    scale = float(math.sqrt(np.mean(ref_sample * ref_sample)))
    if scale == 0.0:
        normalized = 0.0 if rmse == 0.0 else math.inf
    else:
        normalized = rmse / scale
    return ErrorMetrics(
        rmse=rmse,
        normalized_rmse=normalized,
        max_abs_error=float(np.max(np.abs(diff))),
        valid_count=valid_count,
        total_count=total_count,
        reference_finite_count=ref_count,
        candidate_finite_count=cand_count,
        sum_squared_error=sum_squared_error,
    )


def _metrics_to_dict(metrics: ErrorMetrics) -> dict[str, Any]:
    return {
        "rmse": metrics.rmse,
        "normalized_rmse": metrics.normalized_rmse,
        "max_abs_error": metrics.max_abs_error,
        "valid_count": metrics.valid_count,
        "total_count": metrics.total_count,
        "reference_finite_count": metrics.reference_finite_count,
        "candidate_finite_count": metrics.candidate_finite_count,
        "sum_squared_error": metrics.sum_squared_error,
    }


def _infer_target_region(
    candidate: netCDF4.Dataset,
    horizontal_shape: tuple[int, int] | None,
) -> dict[str, Any]:
    if horizontal_shape is None:
        return {
            "status": "not_available",
            "message": "no P-like horizontal shape is available",
        }
    attrs = _read_attrs(candidate)
    from_start = _coerce_pair(attrs.get("TYWRF_FROM_PARENT_START"))
    to_start = _coerce_pair(attrs.get("TYWRF_TO_PARENT_START"))
    ratio = _coerce_int(attrs.get("TYWRF_PARENT_GRID_RATIO"))
    target_count = _coerce_int(attrs.get("TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT"))
    missing = []
    if from_start is None:
        missing.append("TYWRF_FROM_PARENT_START")
    if to_start is None:
        missing.append("TYWRF_TO_PARENT_START")
    if ratio is None:
        missing.append("TYWRF_PARENT_GRID_RATIO")
    if target_count is None:
        missing.append("TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT")
    if missing:
        return {
            "status": "not_available",
            "missing_metadata": missing,
            "message": "cannot infer refreshed target region from candidate metadata",
        }
    if ratio <= 0:
        return {
            "status": "not_available",
            "message": "TYWRF_PARENT_GRID_RATIO must be positive",
        }

    ny, nx = horizontal_shape
    child_delta_i = (to_start[0] - from_start[0]) * ratio
    child_delta_j = (to_start[1] - from_start[1]) * ratio
    if abs(child_delta_i) >= nx or abs(child_delta_j) >= ny:
        return {
            "status": "not_available",
            "from_parent_start": _pair_to_dict(from_start),
            "to_parent_start": _pair_to_dict(to_start),
            "parent_grid_ratio": ratio,
            "child_delta": {"i": child_delta_i, "j": child_delta_j},
            "horizontal_shape": {"ny": ny, "nx": nx},
            "metadata_target_column_count": target_count,
            "message": "inferred child shift exceeds the candidate horizontal shape",
        }

    mask = np.zeros((ny, nx), dtype=bool)
    if child_delta_i > 0:
        mask[:, nx - child_delta_i :] = True
    elif child_delta_i < 0:
        mask[:, : -child_delta_i] = True
    if child_delta_j > 0:
        mask[ny - child_delta_j :, :] = True
    elif child_delta_j < 0:
        mask[: -child_delta_j, :] = True

    inferred_count = int(np.count_nonzero(mask))
    count_matches = inferred_count == target_count
    return {
        "status": "available" if count_matches else "available_count_mismatch",
        "from_parent_start": _pair_to_dict(from_start),
        "to_parent_start": _pair_to_dict(to_start),
        "parent_grid_ratio": ratio,
        "child_delta": {"i": child_delta_i, "j": child_delta_j},
        "horizontal_shape": {"ny": ny, "nx": nx},
        "metadata_target_column_count": target_count,
        "inferred_target_column_count": inferred_count,
        "target_column_count_matches_metadata": count_matches,
        "mask": mask,
        "message": (
            None
            if count_matches
            else "inferred exposed column count does not match metadata target count"
        ),
    }


def _region_mask_for_data(data: np.ndarray, mask_2d: np.ndarray) -> np.ndarray:
    leading = (1,) * (data.ndim - 2)
    return np.broadcast_to(mask_2d.reshape((*leading, *mask_2d.shape)), data.shape)


def _region_split(
    variable: str,
    reference: np.ndarray,
    candidate: np.ndarray,
    region: dict[str, Any],
    global_metrics: ErrorMetrics,
) -> dict[str, Any]:
    if variable not in P_LIKE_3D_VARIABLES or reference.ndim < 3:
        return {
            "status": "not_applicable",
            "message": "region split is only computed for P-like 3D pressure fields",
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
            "message": "target region mask does not match variable horizontal shape",
        }

    target_mask = _region_mask_for_data(reference, mask_2d)
    non_target_mask = ~target_mask
    target_metrics = error_metrics(
        reference,
        candidate,
        selection_mask=target_mask,
    )
    non_target_metrics = error_metrics(
        reference,
        candidate,
        selection_mask=non_target_mask,
    )
    target_sse = target_metrics.sum_squared_error
    global_sse = global_metrics.sum_squared_error
    if math.isfinite(target_sse) and math.isfinite(global_sse) and global_sse > 0.0:
        target_error_fraction = target_sse / global_sse
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
        "target_column_count_matches_metadata": region.get(
            "target_column_count_matches_metadata"
        ),
        "message": region.get("message"),
    }


def _horizontal_shape_for_region(
    reference: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    variables: Iterable[str],
    *,
    time_index: int,
) -> tuple[int, int] | None:
    for variable in variables:
        if variable not in P_LIKE_3D_VARIABLES:
            continue
        if variable not in reference.variables or variable not in candidate.variables:
            continue
        try:
            ref_data = _read_variable_data(reference, variable, time_index=time_index)
            cand_data = _read_variable_data(candidate, variable, time_index=time_index)
        except (IndexError, TypeError, ValueError):
            continue
        if ref_data.shape == cand_data.shape and ref_data.ndim >= 3:
            return int(ref_data.shape[-2]), int(ref_data.shape[-1])
    return None


def compare_pressure_variable(
    reference: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    variable: str,
    *,
    thresholds: dict[str, float] | None,
    region: dict[str, Any],
    time_index: int,
) -> VariableAudit:
    threshold = thresholds.get(variable) if thresholds else None
    reference_present = variable in reference.variables
    candidate_present = variable in candidate.variables
    if not reference_present or not candidate_present:
        missing_owner = []
        if not reference_present:
            missing_owner.append("reference")
        if not candidate_present:
            missing_owner.append("candidate")
        return VariableAudit(
            variable=variable,
            status="not_available",
            reference_present=reference_present,
            candidate_present=candidate_present,
            threshold=threshold,
            region_split={"status": "not_available", "message": "variable missing"},
            message=f"variable missing from {' and '.join(missing_owner)}",
        )

    try:
        ref_data = _read_variable_data(reference, variable, time_index=time_index)
        cand_data = _read_variable_data(candidate, variable, time_index=time_index)
    except (IndexError, TypeError, ValueError) as exc:
        return VariableAudit(
            variable=variable,
            status="not_available",
            reference_present=True,
            candidate_present=True,
            threshold=threshold,
            message=f"variable cannot be read: {exc}",
        )

    ref_shape = tuple(int(value) for value in ref_data.shape)
    cand_shape = tuple(int(value) for value in cand_data.shape)
    if ref_shape != cand_shape:
        return VariableAudit(
            variable=variable,
            status="shape_mismatch",
            reference_present=True,
            candidate_present=True,
            reference_shape=ref_shape,
            candidate_shape=cand_shape,
            threshold=threshold,
            message=f"reference shape {ref_shape} != candidate shape {cand_shape}",
        )

    try:
        metrics = error_metrics(ref_data, cand_data)
    except (TypeError, ValueError) as exc:
        return VariableAudit(
            variable=variable,
            status="non_numeric",
            reference_present=True,
            candidate_present=True,
            reference_shape=ref_shape,
            candidate_shape=cand_shape,
            threshold=threshold,
            message=f"variable cannot be compared numerically: {exc}",
        )

    if metrics.valid_count == 0:
        status = "no_valid_samples"
        message = "no finite sample pairs are available"
    elif threshold is not None and (
        math.isnan(metrics.normalized_rmse) or metrics.normalized_rmse > threshold
    ):
        status = "threshold_exceeded"
        message = (
            f"normalized RMSE {metrics.normalized_rmse:g} exceeds "
            f"threshold {threshold:g}"
        )
    else:
        status = "ok"
        message = None

    split = _region_split(variable, ref_data, cand_data, region, metrics)
    return VariableAudit(
        variable=variable,
        status=status,
        reference_present=True,
        candidate_present=True,
        reference_shape=ref_shape,
        candidate_shape=cand_shape,
        threshold=threshold,
        rmse=metrics.rmse,
        normalized_rmse=metrics.normalized_rmse,
        max_abs_error=metrics.max_abs_error,
        valid_count=metrics.valid_count,
        total_count=metrics.total_count,
        reference_finite_count=metrics.reference_finite_count,
        candidate_finite_count=metrics.candidate_finite_count,
        region_split=split,
        message=message,
    )


def _first_failing_pressure_variable(
    variables: list[VariableAudit],
) -> dict[str, Any] | None:
    for audit in variables:
        if audit.status in {"threshold_exceeded", "no_valid_samples", "non_numeric"}:
            return {
                "variable": audit.variable,
                "status": audit.status,
                "normalized_rmse": audit.normalized_rmse,
                "threshold": audit.threshold,
                "message": audit.message,
            }
        if audit.status == "shape_mismatch":
            return {
                "variable": audit.variable,
                "status": audit.status,
                "normalized_rmse": None,
                "threshold": audit.threshold,
                "message": audit.message,
            }
    return None


def _p_region_diagnosis(variables: list[VariableAudit]) -> dict[str, Any]:
    p_audit = next((audit for audit in variables if audit.variable == "P"), None)
    split = p_audit.region_split if p_audit else None
    if not split or split.get("status") not in {"available", "available_count_mismatch"}:
        return {
            "status": "not_available",
            "refreshed_region_p_dominates_global_error": None,
            "p_refreshed_region_error_fraction": None,
            "message": "P region split is not available",
        }
    return {
        "status": split["status"],
        "refreshed_region_p_dominates_global_error": split.get(
            "target_region_dominates_global_error"
        ),
        "p_refreshed_region_error_fraction": split.get("target_error_fraction"),
        "message": split.get("message"),
    }


def audit_pressure_refresh_candidate(
    reference_path: Path,
    candidate_path: Path,
    variables: Iterable[str] = PRESSURE_REFRESH_VARIABLES,
    thresholds: dict[str, float] | None = DEFAULT_THRESHOLDS,
    *,
    time_index: int = -1,
) -> PressureRefreshCandidateAudit:
    variables = tuple(variables)
    with netCDF4.Dataset(reference_path) as reference, netCDF4.Dataset(
        candidate_path
    ) as candidate:
        metadata = _candidate_metadata(candidate)
        horizontal_shape = _horizontal_shape_for_region(
            reference,
            candidate,
            variables,
            time_index=time_index,
        )
        region = _infer_target_region(candidate, horizontal_shape)
        audits = [
            compare_pressure_variable(
                reference,
                candidate,
                variable,
                thresholds=thresholds,
                region=region,
                time_index=time_index,
            )
            for variable in variables
        ]

    public_region = {key: value for key, value in region.items() if key != "mask"}
    counts = Counter(audit.status for audit in audits)
    failing = _first_failing_pressure_variable(audits)
    p_diagnosis = _p_region_diagnosis(audits)
    status = "computed_with_flags" if failing else "computed"
    summary = {
        **dict(sorted(counts.items())),
        "total": len(audits),
        "available": sum(1 for audit in audits if audit.status != "not_available"),
        "not_available": counts.get("not_available", 0),
        "numeric_failures": sum(
            1
            for audit in audits
            if audit.status in {"threshold_exceeded", "no_valid_samples", "non_numeric"}
        ),
        "shape_mismatches": counts.get("shape_mismatch", 0),
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "first_failing_pressure_variable": (
            None if failing is None else failing["variable"]
        ),
        "refreshed_region_p_dominates_global_error": p_diagnosis[
            "refreshed_region_p_dominates_global_error"
        ],
        "p_refreshed_region_error_fraction": p_diagnosis[
            "p_refreshed_region_error_fraction"
        ],
    }
    return PressureRefreshCandidateAudit(
        reference=str(reference_path),
        candidate=str(candidate_path),
        status=status,
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        summary=summary,
        metadata=metadata,
        region=public_region,
        variables=audits,
        diagnosis={
            "first_failing_pressure_variable": failing,
            "p_region_dominance": p_diagnosis,
            "candidate_model_pass": "not_applicable",
            "message": (
                "diagnostic-only pressure-refresh candidate audit; compare-to-WRF "
                "results are error analysis, not validation pass generation"
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


def report_to_dict(report: PressureRefreshCandidateAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: PressureRefreshCandidateAudit,
    *,
    pretty: bool = False,
) -> str:
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
    parser.add_argument("reference", type=Path, help="WRF reference wrfout file")
    parser.add_argument("candidate", type=Path, help="TyWRF candidate wrfout file")
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(PRESSURE_REFRESH_VARIABLES),
        help="Pressure-refresh-related variables to audit",
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
        help="Report metrics without threshold-exceeded status",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for Time-leading variables; defaults to the last record",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        thresholds = None if args.no_thresholds else parse_threshold_overrides(args.threshold)
        report = audit_pressure_refresh_candidate(
            args.reference,
            args.candidate,
            variables=args.variables,
            thresholds=thresholds,
            time_index=args.time_index,
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
