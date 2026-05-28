#!/usr/bin/env python
"""Audit pressure-refresh formula and staging inputs.

This tool is observation-only. Reference-end fields are used only as
diagnostic evidence and never as candidate-generation guidance.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable

import netCDF4
import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.audit_pressure_refresh_candidate import (  # noqa: E402
    _candidate_metadata,
    _coerce_bool,
    _coerce_int,
    _infer_target_region,
    _metrics_to_dict,
    _netcdf_attr_to_json,
    _read_variable_data,
    _region_mask_for_data,
    error_metrics,
)


FORMULA_FIELDS = (
    "P",
    "PB",
    "P+PB",
    "MU",
    "MUB",
    "MU+MUB",
    "PH",
    "PHB",
    "PH+PHB",
    "T",
    "QVAPOR",
    "HGT",
)
REQUIRED_CANDIDATE_FIELDS = ("P", "PB", "MU", "MUB", "PH", "PHB", "T")
LOW_LEVELS = tuple(range(5))
WORST_LEVEL_COUNT = 5
LOW_LEVEL_P_BIAS_PA = 500.0
WEAK_COMPANION_DIFF = 100.0
INPUT_NORMALIZED_RMSE_WARNING = 0.05
SOURCE_RELATIVE_PATHS = (
    "src/dynamics/pressure_refresh.cpp",
    "src/dynamics/pressure_refresh_hook.cpp",
    "src/dynamics/pressure_refresh_staging.cpp",
    "include/tywrf/dynamics/pressure_refresh.hpp",
    "include/tywrf/dynamics/pressure_refresh_hook.hpp",
    "include/tywrf/dynamics/pressure_refresh_staging.hpp",
)
FORMULA_TERMS = (
    "MU+MUB",
    "PH+PHB",
    "PB subtraction",
    "ALB",
    "C3F/C4F",
    "C3H/C4H",
    "P_TOP",
    "theta_m",
)
PRESSURE_COUNT_ATTRS = (
    "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT",
    "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT",
    "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT",
    "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS",
    "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS",
    "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS",
    "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS",
    "TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT",
    "TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT",
    "TYWRF_WOULD_REFRESH_P_POINT_COUNT",
)
PRESSURE_FLAG_ATTRS = (
    "TYWRF_PRESSURE_REFRESH_APPLIED",
    "TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED",
    "TYWRF_PRESSURE_REFRESH_PROVIDER_OK",
    "TYWRF_PRESSURE_REFRESH_STAGING_OK",
    "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT",
    "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO",
    "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED",
    "TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT",
)
CORE_ATTRS = (
    "TYWRF_DIAGNOSTIC_ONLY",
    "TYWRF_GATE_CANDIDATE",
    "TYWRF_INTEGRATOR_OUTPUT",
    "TYWRF_VALIDATION_GATE_ONLY",
    "TYWRF_CANDIDATE_KIND",
    "TYWRF_FROM_PARENT_START",
    "TYWRF_TO_PARENT_START",
    "TYWRF_PARENT_GRID_RATIO",
    "TYWRF_PRESSURE_REFRESH_HELPER_NAME",
    "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE",
    "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE",
    "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE",
    "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_SOURCE",
    "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_PROVENANCE",
)


@dataclass(frozen=True)
class FieldExpression:
    name: str
    owner: str
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    present: bool
    values: np.ndarray | None = None
    shape: tuple[int, ...] | None = None
    missing_inputs: tuple[str, ...] = ()
    message: str | None = None


@dataclass(frozen=True)
class RiskFlag:
    code: str
    severity: str
    message: str
    evidence: dict[str, Any]
    diagnostic_only: bool = True
    candidate_model_pass: str = "not_applicable"


@dataclass(frozen=True)
class PressureRefreshFormulaInputsAudit:
    candidate_end: str
    reference_end: str | None
    source_start: str | None
    producer_json: str | None
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    metadata_contract: dict[str, Any]
    formula_input_metrics: dict[str, Any]
    source_formula_inventory: dict[str, Any]
    risk_flags: list[RiskFlag]
    summary: dict[str, Any]


def _diag(**values: Any) -> dict[str, Any]:
    return {
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        **values,
    }


def _finite_float(value: Any) -> float | None:
    if isinstance(value, (int, float, np.integer, np.floating)):
        numeric = float(value)
        if math.isfinite(numeric):
            return numeric
    return None


def _read_relevant_attrs(dataset: netCDF4.Dataset) -> dict[str, Any]:
    attrs: dict[str, Any] = {}
    for name in dataset.ncattrs():
        if (
            name in CORE_ATTRS
            or name in PRESSURE_COUNT_ATTRS
            or name in PRESSURE_FLAG_ATTRS
            or name.startswith("TYWRF_PRESSURE_REFRESH_")
            or name.startswith("TYWRF_PRESSURE_COMPUTE_")
            or name.startswith("TYWRF_PROVIDER_")
            or name.startswith("TYWRF_WOULD_")
            or name.startswith("TYWRF_SYNC_")
            or "TERRAIN" in name
            or "HGT" in name
        ):
            attrs[name] = _netcdf_attr_to_json(dataset.getncattr(name))
    return attrs


def _producer_json_summary(
    producer_json_path: Path | None,
    payload: dict[str, Any] | None,
) -> dict[str, Any]:
    if producer_json_path is None:
        return _diag(
            status="not_available",
            path=None,
            message="D58 producer JSON was not supplied",
        )
    if payload is None:
        return _diag(
            status="not_available",
            path=str(producer_json_path),
            message="D58 producer JSON could not be loaded",
        )

    flags = payload.get("formula_risk_flags", [])
    if not isinstance(flags, list):
        flags = []
    flag_summaries = [
        {
            "code": item.get("code"),
            "severity": item.get("severity"),
            "message": item.get("message"),
        }
        for item in flags
        if isinstance(item, dict)
    ]
    return _diag(
        status="available",
        path=str(producer_json_path),
        producer_status=payload.get("status"),
        risk_flag_count=len(flag_summaries),
        risk_flags=flag_summaries,
        metadata_status=(payload.get("metadata_contract") or {}).get("status")
        if isinstance(payload.get("metadata_contract"), dict)
        else None,
        inventory_status=(payload.get("code_producer_inventory") or {}).get("status")
        if isinstance(payload.get("code_producer_inventory"), dict)
        else None,
    )


def _metadata_contract(
    candidate_end: netCDF4.Dataset,
    candidate_path: Path,
    producer_json_path: Path | None,
    producer_payload: dict[str, Any] | None,
) -> dict[str, Any]:
    attrs = _read_relevant_attrs(candidate_end)
    counts = {name: _coerce_int(attrs.get(name)) for name in PRESSURE_COUNT_ATTRS}
    flags = {name: _coerce_bool(attrs.get(name)) for name in PRESSURE_FLAG_ATTRS}
    terrain_attrs = {
        name: value
        for name, value in attrs.items()
        if "TERRAIN" in name or "HGT" in name or "ALB" in name
    }
    return _diag(
        status="available",
        candidate_end=str(candidate_path),
        attrs=attrs,
        present_attrs=sorted(attrs),
        counts=counts,
        flags=flags,
        terrain_metadata=terrain_attrs,
        candidate_metadata=_candidate_metadata(candidate_end),
        producer_json=_producer_json_summary(producer_json_path, producer_payload),
        reference_usage_policy=(
            "reference_end inputs are diagnostic-only and are never used to "
            "generate, tune, patch, or select candidate NetCDF fields"
        ),
    )


def _load_json(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("producer JSON must contain an object")
    return payload


def _load_variable(
    dataset: netCDF4.Dataset | None,
    owner: str,
    name: str,
    *,
    time_index: int,
) -> FieldExpression:
    if dataset is None:
        return FieldExpression(
            name=name,
            owner=owner,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=False,
            missing_inputs=(f"{owner}.{name}",),
            message=f"{owner} dataset was not supplied",
        )
    if name not in dataset.variables:
        return FieldExpression(
            name=name,
            owner=owner,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=False,
            missing_inputs=(f"{owner}.{name}",),
            message=f"{owner}.{name} is missing",
        )
    try:
        values = np.asarray(
            _read_variable_data(dataset, name, time_index=time_index),
            dtype=np.float64,
        )
    except (IndexError, TypeError, ValueError) as exc:
        return FieldExpression(
            name=name,
            owner=owner,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=True,
            missing_inputs=(f"{owner}.{name}",),
            message=f"{owner}.{name} cannot be read: {exc}",
        )
    return FieldExpression(
        name=name,
        owner=owner,
        status="available",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        present=True,
        values=values,
        shape=tuple(int(value) for value in values.shape),
    )


def _expression(
    dataset: netCDF4.Dataset | None,
    owner: str,
    name: str,
    *,
    time_index: int,
) -> FieldExpression:
    if "+" not in name:
        return _load_variable(dataset, owner, name, time_index=time_index)

    parts = tuple(part.strip() for part in name.split("+"))
    loaded = [
        _load_variable(dataset, owner, part, time_index=time_index)
        for part in parts
    ]
    missing = tuple(
        missing for item in loaded for missing in item.missing_inputs
    )
    label = f"{owner}.{name}"
    if any(item.status != "available" for item in loaded):
        messages = [item.message for item in loaded if item.message]
        return FieldExpression(
            name=name,
            owner=owner,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=all(item.present for item in loaded),
            missing_inputs=missing,
            message=(
                "; ".join(messages)
                if messages
                else f"{label} inputs are not available"
            ),
        )

    values = [item.values for item in loaded]
    if any(value is None for value in values):
        return FieldExpression(
            name=name,
            owner=owner,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=True,
            missing_inputs=missing or (label,),
            message=f"{label} values are not available",
        )

    shapes = [tuple(value.shape) for value in values if value is not None]
    if len(set(shapes)) != 1:
        return FieldExpression(
            name=name,
            owner=owner,
            status="shape_mismatch",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=True,
            shape=shapes[0],
            message=f"{label} input shapes differ: {shapes}",
        )

    summed = np.zeros(shapes[0], dtype=np.float64)
    for value in values:
        summed = summed + np.asarray(value, dtype=np.float64)
    return FieldExpression(
        name=name,
        owner=owner,
        status="available",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        present=True,
        values=summed,
        shape=tuple(int(value) for value in summed.shape),
    )


def _finite_values(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    masked = np.ma.masked_invalid(np.ma.asarray(values, dtype=np.float64))
    array = np.asarray(masked.filled(np.nan), dtype=np.float64)
    finite = np.isfinite(array) & ~np.ma.getmaskarray(masked)
    return array, finite


def _mean(values: np.ndarray, selection_mask: np.ndarray | None = None) -> float | None:
    array, finite = _finite_values(values)
    if selection_mask is not None:
        selected = np.asarray(selection_mask, dtype=bool)
        if selected.shape != array.shape:
            selected = np.broadcast_to(selected, array.shape)
        finite = finite & selected
    if not np.any(finite):
        return None
    return float(np.mean(array[finite]))


def _mean_diff(
    reference: np.ndarray,
    candidate: np.ndarray,
    selection_mask: np.ndarray | None = None,
) -> float | None:
    ref_values, ref_finite = _finite_values(reference)
    cand_values, cand_finite = _finite_values(candidate)
    finite = ref_finite & cand_finite
    if selection_mask is not None:
        selected = np.asarray(selection_mask, dtype=bool)
        if selected.shape != ref_values.shape:
            selected = np.broadcast_to(selected, ref_values.shape)
        finite = finite & selected
    if not np.any(finite):
        return None
    return float(np.mean(cand_values[finite] - ref_values[finite]))


def _metric_report(
    reference: np.ndarray,
    candidate: np.ndarray,
    selection_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    metrics = error_metrics(reference, candidate, selection_mask=selection_mask)
    return _diag(
        status="computed" if metrics.valid_count else "no_valid_samples",
        **_metrics_to_dict(metrics),
        mean_reference=_mean(reference, selection_mask),
        mean_candidate=_mean(candidate, selection_mask),
        mean_diff=_mean_diff(reference, candidate, selection_mask),
    )


def _target_mask_report(
    values: np.ndarray,
    mask_2d: np.ndarray | None,
) -> tuple[np.ndarray | None, dict[str, Any]]:
    if mask_2d is None:
        return None, _diag(status="not_available", message="target region is missing")
    if values.ndim < 2 or tuple(values.shape[-2:]) != tuple(mask_2d.shape):
        return None, _diag(
            status="not_available",
            message=(
                f"target mask shape {mask_2d.shape} does not match field "
                f"horizontal shape {values.shape[-2:] if values.ndim >= 2 else None}"
            ),
        )
    return _region_mask_for_data(values, mask_2d), _diag(status="available")


def _expression_summary(
    expression: FieldExpression,
    mask_2d: np.ndarray | None,
) -> dict[str, Any]:
    if expression.status != "available" or expression.values is None:
        return _diag(
            status=expression.status,
            present=expression.present,
            shape=expression.shape,
            missing_inputs=list(expression.missing_inputs),
            message=expression.message,
        )

    values = expression.values
    array, finite = _finite_values(values)
    target_mask, mask_report = _target_mask_report(values, mask_2d)
    if target_mask is not None:
        target_finite = finite & target_mask
        target_finite_count = int(np.count_nonzero(target_finite))
        target_mean = float(np.mean(array[target_finite])) if target_finite_count else None
    else:
        target_finite_count = None
        target_mean = None

    finite_count = int(np.count_nonzero(finite))
    return _diag(
        status="available",
        present=True,
        shape=tuple(int(value) for value in values.shape),
        total_count=int(values.size),
        finite_count=finite_count,
        min=None if finite_count == 0 else float(np.min(array[finite])),
        max=None if finite_count == 0 else float(np.max(array[finite])),
        mean=None if finite_count == 0 else float(np.mean(array[finite])),
        target_region=_diag(
            status=mask_report["status"],
            finite_count=target_finite_count,
            mean=target_mean,
            message=mask_report.get("message"),
        ),
    )


def _level_metric(
    reference: np.ndarray,
    candidate: np.ndarray,
    level: int,
    mask_2d: np.ndarray,
) -> dict[str, Any]:
    if reference.ndim < 3:
        return _diag(
            status="not_applicable",
            level=level,
            message="field has no vertical level axis",
        )
    nz = int(reference.shape[-3])
    if level < 0 or level >= nz:
        return _diag(
            status="out_of_range",
            level=level,
            vertical_level_count=nz,
        )
    ref_slice = np.take(reference, level, axis=-3)
    cand_slice = np.take(candidate, level, axis=-3)
    return _diag(level=level, **_metric_report(ref_slice, cand_slice, mask_2d))


def _low_level_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    mask_2d: np.ndarray | None,
    levels: Iterable[int],
) -> list[dict[str, Any]]:
    if mask_2d is None:
        return [
            _diag(status="not_available", level=level, message="target region is missing")
            for level in levels
        ]
    return [
        _level_metric(reference, candidate, int(level), mask_2d)
        for level in levels
    ]


def _worst_level_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    mask_2d: np.ndarray | None,
    count: int,
) -> list[dict[str, Any]]:
    if reference.ndim < 3:
        return []
    if mask_2d is None:
        return []
    levels: list[dict[str, Any]] = []
    for level in range(int(reference.shape[-3])):
        item = _level_metric(reference, candidate, level, mask_2d)
        rmse = _finite_float(item.get("rmse"))
        if item.get("status") == "computed" and rmse is not None:
            levels.append(item)
    levels.sort(key=lambda item: (-float(item["rmse"]), int(item["level"])))
    return levels[:count]


def _comparison_report(
    baseline: FieldExpression,
    candidate: FieldExpression,
    mask_2d: np.ndarray | None,
) -> dict[str, Any]:
    missing = [*baseline.missing_inputs, *candidate.missing_inputs]
    if baseline.status != "available" or candidate.status != "available":
        status = (
            "shape_mismatch"
            if "shape_mismatch" in {baseline.status, candidate.status}
            else "not_available"
        )
        messages = [item for item in (baseline.message, candidate.message) if item]
        return _diag(
            status=status,
            baseline=baseline.owner,
            candidate=candidate.owner,
            baseline_shape=baseline.shape,
            candidate_shape=candidate.shape,
            missing_inputs=missing,
            low_levels=[],
            worst_levels_by_target_rmse=[],
            message=(
                "; ".join(messages)
                if messages
                else "comparison inputs are not available"
            ),
        )
    if baseline.values is None or candidate.values is None:
        return _diag(
            status="not_available",
            baseline=baseline.owner,
            candidate=candidate.owner,
            missing_inputs=missing,
            low_levels=[],
            worst_levels_by_target_rmse=[],
            message="comparison values are not available",
        )
    if baseline.values.shape != candidate.values.shape:
        return _diag(
            status="shape_mismatch",
            baseline=baseline.owner,
            candidate=candidate.owner,
            baseline_shape=tuple(int(value) for value in baseline.values.shape),
            candidate_shape=tuple(int(value) for value in candidate.values.shape),
            missing_inputs=missing,
            low_levels=[],
            worst_levels_by_target_rmse=[],
            message=(
                f"{baseline.owner}.{baseline.name} shape {baseline.values.shape} "
                f"!= {candidate.owner}.{candidate.name} shape {candidate.values.shape}"
            ),
        )

    target_mask, mask_report = _target_mask_report(candidate.values, mask_2d)
    return _diag(
        **_metric_report(baseline.values, candidate.values),
        baseline=baseline.owner,
        candidate=candidate.owner,
        baseline_shape=tuple(int(value) for value in baseline.values.shape),
        candidate_shape=tuple(int(value) for value in candidate.values.shape),
        missing_inputs=missing,
        target_region=(
            _metric_report(baseline.values, candidate.values, target_mask)
            if target_mask is not None
            else mask_report
        ),
        low_levels=_low_level_metrics(
            baseline.values,
            candidate.values,
            mask_2d,
            LOW_LEVELS,
        ),
        worst_levels_by_target_rmse=_worst_level_metrics(
            baseline.values,
            candidate.values,
            mask_2d,
            WORST_LEVEL_COUNT,
        ),
    )


def _horizontal_shape(
    candidate: netCDF4.Dataset,
    *,
    time_index: int,
) -> tuple[int, int] | None:
    for name in ("P", "MU", "HGT", "PB", "PH"):
        if name not in candidate.variables:
            continue
        try:
            values = np.asarray(
                _read_variable_data(candidate, name, time_index=time_index)
            )
        except (IndexError, TypeError, ValueError):
            continue
        if values.ndim >= 2:
            return int(values.shape[-2]), int(values.shape[-1])
    return None


def _public_region(region: dict[str, Any]) -> dict[str, Any]:
    public = {key: value for key, value in region.items() if key != "mask"}
    return _diag(**public)


def _field_metric_entry(
    name: str,
    candidate_end: netCDF4.Dataset,
    reference_end: netCDF4.Dataset | None,
    source_start: netCDF4.Dataset | None,
    mask_2d: np.ndarray | None,
    *,
    time_index: int,
    source_time_index: int,
) -> dict[str, Any]:
    candidate_expr = _expression(
        candidate_end,
        "candidate_end",
        name,
        time_index=time_index,
    )
    reference_expr = _expression(
        reference_end,
        "reference_end",
        name,
        time_index=time_index,
    )
    source_expr = _expression(
        source_start,
        "source_start",
        name,
        time_index=source_time_index,
    )
    return _diag(
        field=name,
        candidate_end=_expression_summary(candidate_expr, mask_2d),
        reference_end=_expression_summary(reference_expr, mask_2d),
        source_start=_expression_summary(source_expr, mask_2d),
        candidate_vs_reference_end=_comparison_report(
            reference_expr,
            candidate_expr,
            mask_2d,
        ),
        candidate_vs_source_start=_comparison_report(
            source_expr,
            candidate_expr,
            mask_2d,
        ),
        reference_end_vs_source_start=_comparison_report(
            source_expr,
            reference_expr,
            mask_2d,
        ),
    )


def _formula_input_metrics(
    candidate_end: netCDF4.Dataset,
    reference_end: netCDF4.Dataset | None,
    source_start: netCDF4.Dataset | None,
    metadata_contract: dict[str, Any],
    *,
    time_index: int,
    source_time_index: int,
) -> dict[str, Any]:
    horizontal_shape = _horizontal_shape(candidate_end, time_index=time_index)
    region = _infer_target_region(candidate_end, horizontal_shape)
    mask_2d = region.get("mask") if isinstance(region.get("mask"), np.ndarray) else None
    fields = {
        name: _field_metric_entry(
            name,
            candidate_end,
            reference_end,
            source_start,
            mask_2d,
            time_index=time_index,
            source_time_index=source_time_index,
        )
        for name in FORMULA_FIELDS
    }
    return _diag(
        status="computed",
        time_index=time_index,
        source_time_index=source_time_index,
        target_region=_public_region(region),
        pressure_refresh_metadata_counts=metadata_contract.get("counts", {}),
        fields=fields,
    )


_FUNCTION_RE = re.compile(
    r"(?P<name>[A-Za-z_][\w:]*)\s*\([^;{}]*\)\s*(?:noexcept\s*)?(?:const\s*)?\{",
    re.MULTILINE,
)


def _line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def _brace_end(text: str, open_brace_index: int) -> int:
    depth = 0
    for index in range(open_brace_index, len(text)):
        char = text[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return index + 1
    return len(text)


def _term_present(term: str, body: str) -> bool:
    text = body.lower()
    if term == "MU+MUB":
        return "inputs.mu" in text and "inputs.mub" in text
    if term == "PH+PHB":
        return "inputs.ph" in text and "inputs.phb" in text
    if term == "PB subtraction":
        return "total_pressure" in text and "- pb" in text
    if term == "ALB":
        return "inputs.alb" in text or "external_alb" in text or " alb" in text
    if term == "C3F/C4F":
        return "c3f" in text and "c4f" in text
    if term == "C3H/C4H":
        return "c3h" in text and "c4h" in text
    if term == "P_TOP":
        return "p_top" in text
    if term == "theta_m":
        return "theta" in text and (
            "base_potential_temperature" in text or "use_theta_m" in text
        )
    return False


def _source_function_records(path: Path, relative_path: str, text: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for match in _FUNCTION_RE.finditer(text):
        name = match.group("name")
        open_brace = text.find("{", match.end() - 1)
        if open_brace < 0:
            continue
        body = text[open_brace : _brace_end(text, open_brace)]
        formula_terms = {
            term: _term_present(term, body)
            for term in FORMULA_TERMS
        }
        signals = {
            "is_compute_krosa_pressure": name.endswith("compute_krosa_pressure"),
            "is_refresh_krosa_moving_nest_pressure": name.endswith(
                "refresh_krosa_moving_nest_pressure"
            ),
            "is_staging_input_builder": name.endswith(
                "make_krosa_pressure_refresh_inputs"
            ),
            "writes_perturbation_p": "inputs.p(" in body and "=" in body,
            "syncs_base_state_inputs": all(
                token in body for token in ("sync_exposed_3d", "pb", "phb")
            )
            and "mub" in body,
            "uses_base_state_provider_alb": "provider_views.alb" in body
            or "external_alb" in body,
        }
        if any(signals.values()) or "pressure" in name.lower():
            records.append(
                _diag(
                    name=name,
                    relative_path=relative_path,
                    path=str(path),
                    line=_line_number(text, match.start()),
                    formula_terms=formula_terms,
                    signals=signals,
                )
            )
    return records


def _source_formula_inventory(source_root: Path) -> dict[str, Any]:
    files: list[dict[str, Any]] = []
    functions: list[dict[str, Any]] = []
    for relative_path in SOURCE_RELATIVE_PATHS:
        path = source_root / relative_path
        if not path.exists():
            files.append(
                _diag(
                    status="missing",
                    relative_path=relative_path,
                    path=str(path),
                    functions=[],
                )
            )
            continue
        text = path.read_text(encoding="utf-8")
        file_functions = _source_function_records(path, relative_path, text)
        functions.extend(file_functions)
        files.append(
            _diag(
                status="scanned",
                relative_path=relative_path,
                path=str(path),
                functions=file_functions,
            )
        )

    compute_functions = [
        item
        for item in functions
        if (item.get("signals") or {}).get("is_compute_krosa_pressure")
    ]
    terms = {
        term: _diag(
            status="present" if any(
                (function.get("formula_terms") or {}).get(term)
                for function in compute_functions
            ) else "missing",
            present_in=[
                {
                    "name": function["name"],
                    "relative_path": function["relative_path"],
                    "line": function["line"],
                }
                for function in compute_functions
                if (function.get("formula_terms") or {}).get(term)
            ],
        )
        for term in FORMULA_TERMS
    }
    missing_terms = [
        term for term, item in terms.items() if item["status"] != "present"
    ]
    return _diag(
        status="scanned" if any(item["status"] == "scanned" for item in files) else "not_available",
        source_root=str(source_root),
        files=files,
        compute_krosa_pressure=_diag(
            status="available" if compute_functions else "not_available",
            definitions=compute_functions,
            formula_terms=terms,
            missing_formula_terms=missing_terms,
        ),
        staging_hooks=_diag(
            refresh_functions=[
                item for item in functions
                if (item.get("signals") or {}).get("is_refresh_krosa_moving_nest_pressure")
            ],
            staging_input_builders=[
                item for item in functions
                if (item.get("signals") or {}).get("is_staging_input_builder")
            ],
            base_state_sync_functions=[
                item for item in functions
                if (item.get("signals") or {}).get("syncs_base_state_inputs")
            ],
        ),
        summary=_diag(
            source_files_scanned=sum(1 for item in files if item["status"] == "scanned"),
            source_files_missing=sum(1 for item in files if item["status"] == "missing"),
            compute_function_count=len(compute_functions),
            missing_formula_terms=missing_terms,
        ),
    )


def _add_flag(
    flags: list[RiskFlag],
    code: str,
    severity: str,
    message: str,
    evidence: dict[str, Any],
) -> None:
    flags.append(
        RiskFlag(
            code=code,
            severity=severity,
            message=message,
            evidence=evidence,
        )
    )


def _field(metrics: dict[str, Any], name: str) -> dict[str, Any]:
    fields = metrics.get("fields")
    if not isinstance(fields, dict):
        return {}
    item = fields.get(name)
    return item if isinstance(item, dict) else {}


def _first_low_level_flag_item(
    comparison: dict[str, Any],
    threshold: float,
) -> dict[str, Any] | None:
    levels = comparison.get("low_levels")
    if not isinstance(levels, list):
        return None
    for item in levels:
        if not isinstance(item, dict):
            continue
        mean_diff = _finite_float(item.get("mean_diff"))
        rmse = _finite_float(item.get("rmse"))
        if (
            mean_diff is not None
            and abs(mean_diff) >= threshold
        ) or (rmse is not None and rmse >= threshold):
            return item
    return None


def _comparison_low_level(
    metrics: dict[str, Any],
    field_name: str,
    comparison_name: str,
    level: int,
) -> dict[str, Any] | None:
    comparison = _field(metrics, field_name).get(comparison_name)
    if not isinstance(comparison, dict):
        return None
    for item in comparison.get("low_levels", []):
        if isinstance(item, dict) and item.get("level") == level:
            return item
    return None


def _close(a: float, b: float) -> bool:
    return abs(a - b) <= max(50.0, 0.10 * max(abs(a), abs(b), 1.0))


def _metadata_risks(
    metadata: dict[str, Any],
    flags: list[RiskFlag],
) -> None:
    terrain_values = [
        value
        for value in (metadata.get("terrain_metadata") or {}).values()
    ]
    typed_flags = metadata.get("flags") or {}
    if typed_flags.get("TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT") is True or any(
        "moved_candidate_HGT" in str(value) for value in terrain_values
    ):
        _add_flag(
            flags,
            "terrain_provenance_uses_moved_candidate_hgt",
            "warning",
            "pressure-refresh metadata indicates terrain came from moved candidate HGT",
            {
                "terrain_metadata": metadata.get("terrain_metadata"),
                "provider_terrain_uses_moved_candidate_hgt": typed_flags.get(
                    "TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT"
                ),
            },
        )

    counts = metadata.get("counts") or {}
    changed_p = _coerce_int(counts.get("TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS"))
    companion_changed = [
        _coerce_int(counts.get("TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS")),
        _coerce_int(counts.get("TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS")),
        _coerce_int(counts.get("TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS")),
    ]
    if changed_p and all((value or 0) < changed_p for value in companion_changed):
        _add_flag(
            flags,
            "pressure_refresh_changes_p_more_than_base_state_companions",
            "info",
            "metadata shows more changed P points than changed PB/MUB/PHB companion points",
            {
                "changed_p_points": changed_p,
                "changed_pb_points": companion_changed[0],
                "changed_mub_points": companion_changed[1],
                "changed_phb_points": companion_changed[2],
            },
        )

    producer = metadata.get("producer_json") or {}
    producer_flags = producer.get("risk_flags")
    if isinstance(producer_flags, list):
        warning_codes = [
            item.get("code")
            for item in producer_flags
            if isinstance(item, dict)
            and item.get("severity") in {"warning", "error", "critical"}
        ]
        if warning_codes:
            _add_flag(
                flags,
                "producer_json_pressure_formula_risks_present",
                "warning",
                "D58 producer JSON carries pressure-refresh formula/staging risks",
                {"producer_risk_codes": warning_codes},
            )


def _metric_risks(metrics: dict[str, Any], flags: list[RiskFlag]) -> None:
    for field_name in REQUIRED_CANDIDATE_FIELDS:
        field = _field(metrics, field_name)
        candidate = field.get("candidate_end") if isinstance(field, dict) else {}
        if isinstance(candidate, dict) and candidate.get("status") != "available":
            _add_flag(
                flags,
                "missing_candidate_formula_input",
                "warning",
                f"candidate is missing required formula input {field_name}",
                {
                    "field": field_name,
                    "status": candidate.get("status"),
                    "message": candidate.get("message"),
                },
            )

    p_comparison = _field(metrics, "P").get("candidate_vs_reference_end")
    p_level = (
        _first_low_level_flag_item(p_comparison, LOW_LEVEL_P_BIAS_PA)
        if isinstance(p_comparison, dict)
        else None
    )
    if p_level is not None:
        _add_flag(
            flags,
            "low_level_target_p_bias",
            "warning",
            "low-level target-region perturbation P differs strongly from reference",
            {
                "level": p_level.get("level"),
                "mean_diff": p_level.get("mean_diff"),
                "rmse": p_level.get("rmse"),
                "valid_count": p_level.get("valid_count"),
            },
        )

    if p_level is not None:
        level = int(p_level["level"])
        p_diff = _finite_float(p_level.get("mean_diff"))
        full = _comparison_low_level(
            metrics,
            "P+PB",
            "candidate_vs_reference_end",
            level,
        )
        pb = _comparison_low_level(
            metrics,
            "PB",
            "candidate_vs_reference_end",
            level,
        )
        full_diff = _finite_float((full or {}).get("mean_diff"))
        pb_diff = _finite_float((pb or {}).get("mean_diff"))
        if (
            p_diff is not None
            and full_diff is not None
            and abs(p_diff) >= LOW_LEVEL_P_BIAS_PA
            and _close(full_diff, p_diff)
            and (pb_diff is None or abs(pb_diff) <= WEAK_COMPANION_DIFF)
        ):
            _add_flag(
                flags,
                "p_plus_pb_inherits_perturbation_p_bias",
                "warning",
                "P+PB low-level target-region bias closely tracks perturbation P",
                {
                    "level": level,
                    "p_mean_diff": p_diff,
                    "p_plus_pb_mean_diff": full_diff,
                    "pb_mean_diff": pb_diff,
                },
            )

        weak_companions: dict[str, float | None] = {}
        for name in ("PB", "MUB", "PHB"):
            item = _comparison_low_level(
                metrics,
                name,
                "candidate_vs_reference_end",
                level,
            )
            weak_companions[name] = _finite_float((item or {}).get("mean_diff"))
        if (
            p_diff is not None
            and abs(p_diff) >= LOW_LEVEL_P_BIAS_PA
            and all(
                value is not None and abs(value) <= WEAK_COMPANION_DIFF
                for value in weak_companions.values()
            )
        ):
            _add_flag(
                flags,
                "p_bias_with_weak_base_state_companion_diffs",
                "warning",
                "large P bias appears while PB/MUB/PHB input differences are weak",
                {"level": level, "p_mean_diff": p_diff, **weak_companions},
            )

    for field_name in ("MU+MUB", "PH+PHB", "T", "HGT"):
        comparison = _field(metrics, field_name).get("candidate_vs_reference_end")
        target = comparison.get("target_region") if isinstance(comparison, dict) else {}
        normalized = _finite_float(target.get("normalized_rmse") if isinstance(target, dict) else None)
        if normalized is not None and normalized > INPUT_NORMALIZED_RMSE_WARNING:
            _add_flag(
                flags,
                "formula_input_target_region_rmse_high",
                "warning",
                f"{field_name} target-region normalized RMSE exceeds diagnostic threshold",
                {"field": field_name, "target_normalized_rmse": normalized},
            )


def _source_risks(inventory: dict[str, Any], flags: list[RiskFlag]) -> None:
    compute = inventory.get("compute_krosa_pressure") or {}
    missing = compute.get("missing_formula_terms") if isinstance(compute, dict) else []
    if missing:
        _add_flag(
            flags,
            "source_formula_terms_missing",
            "warning",
            "static inventory did not find all expected compute_krosa_pressure formula terms",
            {"missing_formula_terms": missing},
        )
    if isinstance(compute, dict) and compute.get("status") == "not_available":
        _add_flag(
            flags,
            "compute_krosa_pressure_not_found",
            "warning",
            "static inventory did not find compute_krosa_pressure",
            {"source_root": inventory.get("source_root")},
        )


def _risk_flags(
    metadata: dict[str, Any],
    metrics: dict[str, Any],
    inventory: dict[str, Any],
) -> list[RiskFlag]:
    flags: list[RiskFlag] = []
    _metadata_risks(metadata, flags)
    _metric_risks(metrics, flags)
    _source_risks(inventory, flags)
    return flags


def audit_pressure_refresh_formula_inputs(
    candidate_end_path: Path,
    *,
    reference_end_path: Path | None = None,
    source_start_path: Path | None = None,
    producer_json_path: Path | None = None,
    source_root: Path = Path("."),
    time_index: int = -1,
    source_time_index: int = -1,
) -> PressureRefreshFormulaInputsAudit:
    producer_payload = _load_json(producer_json_path)
    with ExitStack() as stack:
        candidate_end = stack.enter_context(netCDF4.Dataset(candidate_end_path))
        reference_end = (
            stack.enter_context(netCDF4.Dataset(reference_end_path))
            if reference_end_path is not None
            else None
        )
        source_start = (
            stack.enter_context(netCDF4.Dataset(source_start_path))
            if source_start_path is not None
            else None
        )
        metadata = _metadata_contract(
            candidate_end,
            candidate_end_path,
            producer_json_path,
            producer_payload,
        )
        metrics = _formula_input_metrics(
            candidate_end,
            reference_end,
            source_start,
            metadata,
            time_index=time_index,
            source_time_index=source_time_index,
        )

    inventory = _source_formula_inventory(source_root)
    flags = _risk_flags(metadata, metrics, inventory)
    computed_fields = sum(
        1
        for field in metrics.get("fields", {}).values()
        if isinstance(field, dict)
        and (field.get("candidate_end") or {}).get("status") == "available"
    )
    summary = _diag(
        status="computed_with_flags" if flags else "computed",
        formula_field_count=len(FORMULA_FIELDS),
        available_candidate_field_count=computed_fields,
        risk_flag_count=len(flags),
        risk_flag_codes=[flag.code for flag in flags],
        target_region_status=(metrics.get("target_region") or {}).get("status"),
        strict_gate_status="not_evaluated",
        message=(
            "diagnostic-only pressure-refresh formula/input audit; no candidate "
            "NetCDF values are created, changed, tuned, patched, or selected"
        ),
    )
    return PressureRefreshFormulaInputsAudit(
        candidate_end=str(candidate_end_path),
        reference_end=None if reference_end_path is None else str(reference_end_path),
        source_start=None if source_start_path is None else str(source_start_path),
        producer_json=None if producer_json_path is None else str(producer_json_path),
        status=summary["status"],
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        metadata_contract=metadata,
        formula_input_metrics=metrics,
        source_formula_inventory=inventory,
        risk_flags=flags,
        summary=summary,
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


def report_to_dict(report: PressureRefreshFormulaInputsAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: PressureRefreshFormulaInputsAudit,
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
        "--candidate-end",
        required=True,
        type=Path,
        help="TyWRF candidate/end NetCDF file to inspect",
    )
    parser.add_argument(
        "--reference-end",
        type=Path,
        help="Optional WRF reference-end NetCDF file for diagnostic comparison only",
    )
    parser.add_argument(
        "--source-start",
        type=Path,
        help="Optional source/start NetCDF file for staging/evolution diagnostics",
    )
    parser.add_argument(
        "--producer-json",
        type=Path,
        help="Optional D58 audit_pressure_refresh_producer.py JSON report",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=Path("."),
        help="Repository/source root containing pressure-refresh C++ sources",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for candidate/reference Time-leading variables",
    )
    parser.add_argument(
        "--source-time-index",
        type=int,
        default=-1,
        help="Time index used for source/start Time-leading variables",
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
        report = audit_pressure_refresh_formula_inputs(
            args.candidate_end,
            reference_end_path=args.reference_end,
            source_start_path=args.source_start,
            producer_json_path=args.producer_json,
            source_root=args.source_root,
            time_index=args.time_index,
            source_time_index=args.source_time_index,
        )
    except (OSError, ValueError, IndexError, json.JSONDecodeError) as exc:
        parser.error(str(exc))

    text = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
