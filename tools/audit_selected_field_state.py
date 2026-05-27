#!/usr/bin/env python
"""Audit selected state-field errors for U/V/MU against a WRF reference file."""

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


DEFAULT_STATE_VARIABLES = ("U", "V", "MU")
DEFAULT_THRESHOLDS = {name: 0.05 for name in DEFAULT_STATE_VARIABLES}
METADATA_ATTRS = (
    "TYWRF_DIAGNOSTIC_ONLY",
    "TYWRF_GATE_CANDIDATE",
    "TYWRF_INTEGRATOR_OUTPUT",
    "TYWRF_VALIDATION_GATE_ONLY",
    "TYWRF_CANDIDATE_KIND",
    "TYWRF_FROM_PARENT_START",
    "TYWRF_TO_PARENT_START",
    "TYWRF_PARENT_GRID_RATIO",
)
BOOL_TRUE_VALUES = {"1", "true", "t", "yes", "y", "on"}
BOOL_FALSE_VALUES = {"0", "false", "f", "no", "n", "off"}


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
class StateVariableAudit:
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
    persistence: dict[str, Any] | None = None
    message: str | None = None


@dataclass(frozen=True)
class SelectedFieldStateAudit:
    reference_end: str
    candidate: str
    candidate_start: str | None
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    summary: dict[str, Any]
    metadata: dict[str, Any]
    region: dict[str, Any]
    variables: list[StateVariableAudit]
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
        if name in METADATA_ATTRS:
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


def _pair_to_dict(value: tuple[int, int] | None) -> dict[str, int] | None:
    if value is None:
        return None
    return {"i": value[0], "j": value[1]}


def _candidate_metadata(candidate: netCDF4.Dataset) -> dict[str, Any]:
    attrs = _read_attrs(candidate)
    return {
        "status": "available",
        "attrs": attrs,
        "present_attrs": sorted(attrs),
        "missing_attrs": [name for name in METADATA_ATTRS if name not in attrs],
        "disposition": {
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "gate_candidate": _coerce_bool(attrs.get("TYWRF_GATE_CANDIDATE")),
            "candidate_diagnostic_only": _coerce_bool(
                attrs.get("TYWRF_DIAGNOSTIC_ONLY")
            ),
            "integrator_output": _coerce_bool(attrs.get("TYWRF_INTEGRATOR_OUTPUT")),
            "validation_gate_only": _coerce_bool(
                attrs.get("TYWRF_VALIDATION_GATE_ONLY")
            ),
            "candidate_kind": attrs.get("TYWRF_CANDIDATE_KIND"),
            "from_parent_start": _pair_to_dict(
                _coerce_pair(attrs.get("TYWRF_FROM_PARENT_START"))
            ),
            "to_parent_start": _pair_to_dict(
                _coerce_pair(attrs.get("TYWRF_TO_PARENT_START"))
            ),
            "parent_grid_ratio": _coerce_int(attrs.get("TYWRF_PARENT_GRID_RATIO")),
        },
    }


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


def _infer_target_region(candidate: netCDF4.Dataset) -> dict[str, Any]:
    attrs = _read_attrs(candidate)
    from_start = _coerce_pair(attrs.get("TYWRF_FROM_PARENT_START"))
    to_start = _coerce_pair(attrs.get("TYWRF_TO_PARENT_START"))
    ratio = _coerce_int(attrs.get("TYWRF_PARENT_GRID_RATIO"))
    missing = []
    if from_start is None:
        missing.append("TYWRF_FROM_PARENT_START")
    if to_start is None:
        missing.append("TYWRF_TO_PARENT_START")
    if ratio is None:
        missing.append("TYWRF_PARENT_GRID_RATIO")
    if missing:
        return {
            "status": "not_available",
            "missing_metadata": missing,
            "message": "cannot infer moving-nest exposed region from candidate metadata",
        }
    if ratio <= 0:
        return {
            "status": "not_available",
            "message": "TYWRF_PARENT_GRID_RATIO must be positive",
        }

    child_delta_i = (to_start[0] - from_start[0]) * ratio
    child_delta_j = (to_start[1] - from_start[1]) * ratio
    return {
        "status": "available",
        "from_parent_start": _pair_to_dict(from_start),
        "to_parent_start": _pair_to_dict(to_start),
        "parent_grid_ratio": ratio,
        "child_delta": {"i": child_delta_i, "j": child_delta_j},
        "message": None,
    }


def _mass_horizontal_shape(
    reference: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    *,
    time_index: int,
) -> tuple[int, int] | None:
    if "MU" not in reference.variables or "MU" not in candidate.variables:
        return None
    try:
        ref_mu = _read_variable_data(reference, "MU", time_index=time_index)
        cand_mu = _read_variable_data(candidate, "MU", time_index=time_index)
    except (IndexError, TypeError, ValueError):
        return None
    if ref_mu.shape != cand_mu.shape or ref_mu.ndim < 2:
        return None
    return int(ref_mu.shape[-2]), int(ref_mu.shape[-1])


def _edge_mask(shape: tuple[int, int], delta_i: int, delta_j: int) -> np.ndarray:
    ny, nx = shape
    if abs(delta_i) > nx or abs(delta_j) > ny:
        raise ValueError("moving-nest child delta exceeds horizontal shape")
    mask = np.zeros((ny, nx), dtype=bool)
    if delta_i > 0:
        mask[:, nx - delta_i :] = True
    elif delta_i < 0:
        mask[:, : -delta_i] = True
    if delta_j > 0:
        mask[ny - delta_j :, :] = True
    elif delta_j < 0:
        mask[: -delta_j, :] = True
    return mask


def _field_region_mask(
    variable: str,
    horizontal_shape: tuple[int, int],
    mass_shape: tuple[int, int] | None,
    region: dict[str, Any],
) -> tuple[np.ndarray | None, dict[str, Any]]:
    if region.get("status") != "available":
        return None, {
            "status": "not_available",
            "message": region.get("message", "moving-nest exposed region unavailable"),
        }

    child_delta = region.get("child_delta") or {}
    delta_i = _coerce_int(child_delta.get("i"))
    delta_j = _coerce_int(child_delta.get("j"))
    if delta_i is None or delta_j is None:
        return None, {
            "status": "not_available",
            "message": "moving-nest child delta is not available",
        }

    ny, nx = horizontal_shape
    shape_kind = "mass_grid"
    mask_delta_i = abs(delta_i)
    mask_delta_j = abs(delta_j)
    if variable == "U":
        if mass_shape is None:
            return None, {
                "status": "not_available",
                "message": "mass-grid shape is required to split staggered U errors",
            }
        mass_ny, mass_nx = mass_shape
        if (ny, nx) == (mass_ny, mass_nx):
            shape_kind = "mass_grid"
        elif (ny, nx) == (mass_ny, mass_nx + 1):
            shape_kind = "x_staggered"
            if delta_i != 0:
                mask_delta_i += 1
        else:
            return None, {
                "status": "not_available",
                "horizontal_shape": {"ny": ny, "nx": nx},
                "mass_horizontal_shape": {"ny": mass_ny, "nx": mass_nx},
                "message": "U horizontal shape is not mass-grid or x-staggered",
            }
    elif variable == "V":
        if mass_shape is None:
            return None, {
                "status": "not_available",
                "message": "mass-grid shape is required to split staggered V errors",
            }
        mass_ny, mass_nx = mass_shape
        if (ny, nx) == (mass_ny, mass_nx):
            shape_kind = "mass_grid"
        elif (ny, nx) == (mass_ny + 1, mass_nx):
            shape_kind = "y_staggered"
            if delta_j != 0:
                mask_delta_j += 1
        else:
            return None, {
                "status": "not_available",
                "horizontal_shape": {"ny": ny, "nx": nx},
                "mass_horizontal_shape": {"ny": mass_ny, "nx": mass_nx},
                "message": "V horizontal shape is not mass-grid or y-staggered",
            }

    signed_i = mask_delta_i if delta_i > 0 else -mask_delta_i
    signed_j = mask_delta_j if delta_j > 0 else -mask_delta_j
    try:
        mask = _edge_mask((ny, nx), signed_i, signed_j)
    except ValueError as exc:
        return None, {
            "status": "not_available",
            "horizontal_shape": {"ny": ny, "nx": nx},
            "message": str(exc),
        }

    return mask, {
        "status": "available",
        "shape_kind": shape_kind,
        "horizontal_shape": {"ny": ny, "nx": nx},
        "target_point_count_2d": int(np.count_nonzero(mask)),
        "overlap_point_count_2d": int(mask.size - np.count_nonzero(mask)),
        "message": None,
    }


def _region_mask_for_data(data: np.ndarray, mask_2d: np.ndarray) -> np.ndarray:
    leading = (1,) * (data.ndim - 2)
    return np.broadcast_to(mask_2d.reshape((*leading, *mask_2d.shape)), data.shape)


def _region_split(
    variable: str,
    reference: np.ndarray,
    candidate: np.ndarray,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    global_metrics: ErrorMetrics,
) -> dict[str, Any]:
    if reference.ndim < 2:
        return {
            "status": "not_available",
            "message": "variable has no horizontal dimensions",
        }

    mask_2d, mask_info = _field_region_mask(
        variable,
        (int(reference.shape[-2]), int(reference.shape[-1])),
        mass_shape,
        region,
    )
    if mask_2d is None:
        return mask_info

    target_mask = _region_mask_for_data(reference, mask_2d)
    overlap_mask = ~target_mask
    target_metrics = error_metrics(
        reference,
        candidate,
        selection_mask=target_mask,
    )
    overlap_metrics = error_metrics(
        reference,
        candidate,
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
        "target_region": _metrics_to_dict(target_metrics),
        "overlap_region": _metrics_to_dict(overlap_metrics),
        "target_error_fraction": target_error_fraction,
        "target_region_dominates_global_error": (
            None if target_error_fraction is None else target_error_fraction > 0.5
        ),
    }


def _persistence_comparison(
    variable: str,
    reference_end: np.ndarray,
    candidate: np.ndarray,
    candidate_start_dataset: netCDF4.Dataset | None,
    *,
    start_time_index: int,
    close_to_start_threshold: float,
    start_closeness_ratio: float,
) -> dict[str, Any]:
    if candidate_start_dataset is None:
        return {
            "status": "not_requested",
            "candidate_close_to_start": None,
            "message": "candidate-start file was not provided",
        }
    if variable not in candidate_start_dataset.variables:
        return {
            "status": "not_available",
            "candidate_start_present": False,
            "candidate_close_to_start": None,
            "message": "variable missing from candidate-start file",
        }
    try:
        start = _read_variable_data(
            candidate_start_dataset,
            variable,
            time_index=start_time_index,
        )
    except (IndexError, TypeError, ValueError) as exc:
        return {
            "status": "not_available",
            "candidate_start_present": True,
            "candidate_close_to_start": None,
            "message": f"candidate-start variable cannot be read: {exc}",
        }

    start_shape = tuple(int(value) for value in start.shape)
    if start_shape != tuple(int(value) for value in reference_end.shape):
        return {
            "status": "shape_mismatch",
            "candidate_start_present": True,
            "candidate_start_shape": start_shape,
            "reference_shape": tuple(int(value) for value in reference_end.shape),
            "candidate_close_to_start": None,
            "message": "candidate-start shape differs from reference-end shape",
        }

    try:
        candidate_vs_start = error_metrics(start, candidate)
        reference_end_vs_start = error_metrics(start, reference_end)
    except (TypeError, ValueError) as exc:
        return {
            "status": "non_numeric",
            "candidate_start_present": True,
            "candidate_start_shape": start_shape,
            "candidate_close_to_start": None,
            "message": f"candidate-start comparison failed: {exc}",
        }

    cand_norm = candidate_vs_start.normalized_rmse
    ref_norm = reference_end_vs_start.normalized_rmse
    if not math.isfinite(cand_norm):
        candidate_close = False
    elif cand_norm <= close_to_start_threshold:
        candidate_close = True
    elif math.isfinite(ref_norm) and ref_norm > close_to_start_threshold:
        candidate_close = cand_norm <= ref_norm * start_closeness_ratio
    else:
        candidate_close = False

    if math.isfinite(cand_norm) and math.isfinite(ref_norm) and ref_norm > 0.0:
        candidate_start_distance_fraction = cand_norm / ref_norm
    else:
        candidate_start_distance_fraction = None

    return {
        "status": "available",
        "candidate_start_present": True,
        "candidate_start_shape": start_shape,
        "candidate_vs_start": _metrics_to_dict(candidate_vs_start),
        "reference_end_vs_start": _metrics_to_dict(reference_end_vs_start),
        "candidate_start_distance_fraction_of_reference_evolution": (
            candidate_start_distance_fraction
        ),
        "candidate_close_to_start": candidate_close,
        "close_to_start_threshold": close_to_start_threshold,
        "start_closeness_ratio": start_closeness_ratio,
    }


def compare_state_variable(
    reference_end: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    candidate_start: netCDF4.Dataset | None,
    variable: str,
    *,
    thresholds: dict[str, float] | None,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    time_index: int,
    start_time_index: int,
    close_to_start_threshold: float,
    start_closeness_ratio: float,
) -> StateVariableAudit:
    threshold = thresholds.get(variable) if thresholds else None
    reference_present = variable in reference_end.variables
    candidate_present = variable in candidate.variables
    if not reference_present or not candidate_present:
        missing_owner = []
        if not reference_present:
            missing_owner.append("reference-end")
        if not candidate_present:
            missing_owner.append("candidate")
        return StateVariableAudit(
            variable=variable,
            status="not_available",
            reference_present=reference_present,
            candidate_present=candidate_present,
            threshold=threshold,
            region_split={"status": "not_available", "message": "variable missing"},
            persistence={
                "status": "not_available",
                "candidate_close_to_start": None,
                "message": "variable missing",
            },
            message=f"variable missing from {' and '.join(missing_owner)}",
        )

    try:
        ref_data = _read_variable_data(reference_end, variable, time_index=time_index)
        cand_data = _read_variable_data(candidate, variable, time_index=time_index)
    except (IndexError, TypeError, ValueError) as exc:
        return StateVariableAudit(
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
        return StateVariableAudit(
            variable=variable,
            status="shape_mismatch",
            reference_present=True,
            candidate_present=True,
            reference_shape=ref_shape,
            candidate_shape=cand_shape,
            threshold=threshold,
            message=f"reference-end shape {ref_shape} != candidate shape {cand_shape}",
        )

    try:
        metrics = error_metrics(ref_data, cand_data)
    except (TypeError, ValueError) as exc:
        return StateVariableAudit(
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

    split = _region_split(variable, ref_data, cand_data, region, mass_shape, metrics)
    persistence = _persistence_comparison(
        variable,
        ref_data,
        cand_data,
        candidate_start,
        start_time_index=start_time_index,
        close_to_start_threshold=close_to_start_threshold,
        start_closeness_ratio=start_closeness_ratio,
    )
    return StateVariableAudit(
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
        persistence=persistence,
        message=message,
    )


def _first_failing_state_variable(
    variables: list[StateVariableAudit],
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


def _field_diagnosis(variables: list[StateVariableAudit]) -> dict[str, Any]:
    target_fractions: dict[str, float | None] = {}
    target_dominates: dict[str, bool | None] = {}
    close_to_start: dict[str, bool | None] = {}
    dominated_fields: list[str] = []

    for audit in variables:
        split = audit.region_split or {}
        if split.get("status") == "available":
            fraction = split.get("target_error_fraction")
            dominates = split.get("target_region_dominates_global_error")
        else:
            fraction = None
            dominates = None
        target_fractions[audit.variable] = fraction
        target_dominates[audit.variable] = dominates
        if dominates is True:
            dominated_fields.append(audit.variable)

        persistence = audit.persistence or {}
        close_to_start[audit.variable] = persistence.get("candidate_close_to_start")

    available_dominance = [
        value for value in target_dominates.values() if value is not None
    ]
    if available_dominance:
        any_dominates: bool | None = any(available_dominance)
    else:
        any_dominates = None

    return {
        "target_region_error_fraction": target_fractions,
        "target_region_dominates_error": target_dominates,
        "candidate_close_to_start": close_to_start,
        "fields_with_target_region_dominated_error": dominated_fields,
        "exposed_or_interpolated_cells_dominate_state_error": any_dominates,
    }


def audit_selected_field_state(
    reference_end_path: Path,
    candidate_path: Path,
    *,
    candidate_start_path: Path | None = None,
    variables: Iterable[str] = DEFAULT_STATE_VARIABLES,
    thresholds: dict[str, float] | None = DEFAULT_THRESHOLDS,
    time_index: int = -1,
    start_time_index: int = -1,
    close_to_start_threshold: float = 1.0e-6,
    start_closeness_ratio: float = 0.25,
) -> SelectedFieldStateAudit:
    variables = tuple(variables)
    with netCDF4.Dataset(reference_end_path) as reference_end, netCDF4.Dataset(
        candidate_path
    ) as candidate:
        candidate_start_context = (
            netCDF4.Dataset(candidate_start_path) if candidate_start_path else None
        )
        try:
            metadata = _candidate_metadata(candidate)
            region = _infer_target_region(candidate)
            mass_shape = _mass_horizontal_shape(
                reference_end,
                candidate,
                time_index=time_index,
            )
            audits = [
                compare_state_variable(
                    reference_end,
                    candidate,
                    candidate_start_context,
                    variable,
                    thresholds=thresholds,
                    region=region,
                    mass_shape=mass_shape,
                    time_index=time_index,
                    start_time_index=start_time_index,
                    close_to_start_threshold=close_to_start_threshold,
                    start_closeness_ratio=start_closeness_ratio,
                )
                for variable in variables
            ]
        finally:
            if candidate_start_context is not None:
                candidate_start_context.close()

    counts = Counter(audit.status for audit in audits)
    failing = _first_failing_state_variable(audits)
    field_diagnosis = _field_diagnosis(audits)
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
        "first_failing_state_variable": (
            None if failing is None else failing["variable"]
        ),
        "candidate_start_provided": candidate_start_path is not None,
        "target_region_error_fraction": field_diagnosis[
            "target_region_error_fraction"
        ],
        "candidate_close_to_start": field_diagnosis["candidate_close_to_start"],
        "target_region_dominates_error": field_diagnosis[
            "target_region_dominates_error"
        ],
    }
    return SelectedFieldStateAudit(
        reference_end=str(reference_end_path),
        candidate=str(candidate_path),
        candidate_start=None if candidate_start_path is None else str(candidate_start_path),
        status=status,
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        summary=summary,
        metadata=metadata,
        region=region,
        variables=audits,
        diagnosis={
            "first_failing_state_variable": failing,
            **field_diagnosis,
            "candidate_model_pass": "not_applicable",
            "message": (
                "diagnostic-only selected-field state audit; compare-to-WRF "
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


def report_to_dict(report: SelectedFieldStateAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: SelectedFieldStateAudit,
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
    parser.add_argument("reference_end", type=Path, help="WRF reference-end wrfout file")
    parser.add_argument("candidate", type=Path, help="TyWRF candidate wrfout file")
    parser.add_argument(
        "--candidate-start",
        type=Path,
        help="Optional candidate/start-state wrfout file for persistence diagnostics",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_STATE_VARIABLES),
        help="State variables to audit",
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
        help="Time index used for reference-end and candidate Time-leading variables",
    )
    parser.add_argument(
        "--start-time-index",
        type=int,
        default=-1,
        help="Time index used for candidate-start Time-leading variables",
    )
    parser.add_argument(
        "--close-to-start-threshold",
        type=float,
        default=1.0e-6,
        help="Absolute normalized RMSE threshold for exact/near persistence",
    )
    parser.add_argument(
        "--start-closeness-ratio",
        type=float,
        default=0.25,
        help=(
            "Also mark close-to-start when candidate-start distance is at most this "
            "fraction of the reference-end evolution distance"
        ),
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        thresholds = None if args.no_thresholds else parse_threshold_overrides(args.threshold)
        report = audit_selected_field_state(
            args.reference_end,
            args.candidate,
            candidate_start_path=args.candidate_start,
            variables=args.variables,
            thresholds=thresholds,
            time_index=args.time_index,
            start_time_index=args.start_time_index,
            close_to_start_threshold=args.close_to_start_threshold,
            start_closeness_ratio=args.start_closeness_ratio,
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
