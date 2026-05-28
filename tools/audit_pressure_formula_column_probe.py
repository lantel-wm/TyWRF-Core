#!/usr/bin/env python
"""Probe per-column pressure formula terms in moved-nest target cells.

This tool is observation-only. It never creates, modifies, patches, tunes,
selects, or writes candidate NetCDF fields. Reference-end inputs are used only
as diagnostic evidence.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

import netCDF4
import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.audit_pressure_refresh_candidate import (  # noqa: E402
    _candidate_metadata,
    _infer_target_region,
    _netcdf_attr_to_json,
    _read_variable_data,
)
from tools.audit_pressure_refresh_formula_inputs import (  # noqa: E402
    _source_formula_inventory as _d59_source_formula_inventory,
)


DEFAULT_LEVELS = tuple(range(5))
DEFAULT_MAX_COLUMNS = 5
PRESSURE_ERROR_WARNING_PA = 500.0
WEAK_COMPANION_WARNING_ABS = 100.0
FIELD_NAMES = (
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
REQUIRED_FORMULA_INPUTS = ("P", "PB", "MU", "MUB", "PH", "PHB", "T")
SOURCE_FORMULA_TERMS = (
    "MU+MUB",
    "PH+PHB",
    "PB subtraction",
    "ALB",
    "C3F/C4F",
    "C3H/C4H",
    "P_TOP",
    "theta_m",
)
COEFFICIENT_TERMS = ("C3F/C4F", "C3H/C4H", "P_TOP", "ALB", "theta_m")
THETA_METADATA_ATTRS = (
    "USE_THETA_M",
    "TYWRF_USE_THETA_M",
    "TYWRF_PRESSURE_REFRESH_USE_THETA_M",
    "TYWRF_PRESSURE_REFRESH_THERMODYNAMIC_MODE",
    "TYWRF_PRESSURE_REFRESH_FORMULA",
)


@dataclass(frozen=True)
class FieldData:
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
class PressureFormulaColumnProbeAudit:
    candidate_end: str
    reference_end: str | None
    source_start: str | None
    formula_input_json: str | None
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    column_selection: dict[str, Any]
    column_formula_terms: dict[str, Any]
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
    if isinstance(value, (np.ma.core.MaskedConstant,)):
        return None
    if isinstance(value, np.ndarray) and value.shape == ():
        value = value.item()
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        numeric = float(value)
        if math.isfinite(numeric):
            return numeric
    return None


def _load_json(path: Path | None, label: str) -> dict[str, Any] | None:
    if path is None:
        return None
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"{label} JSON must contain an object")
    return payload


def _formula_input_json_summary(
    formula_input_json_path: Path | None,
    payload: dict[str, Any] | None,
) -> dict[str, Any]:
    if formula_input_json_path is None:
        return _diag(
            status="not_available",
            path=None,
            message="D59 formula/input JSON was not supplied",
        )
    if payload is None:
        return _diag(
            status="not_available",
            path=str(formula_input_json_path),
            message="D59 formula/input JSON could not be loaded",
        )

    flags = payload.get("risk_flags", [])
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
    metrics = payload.get("formula_input_metrics")
    if not isinstance(metrics, dict):
        metrics = {}
    return _diag(
        status="available",
        path=str(formula_input_json_path),
        audit_status=payload.get("status"),
        risk_flag_count=len(flag_summaries),
        risk_flags=flag_summaries,
        target_region=metrics.get("target_region"),
        pressure_refresh_metadata_counts=metrics.get(
            "pressure_refresh_metadata_counts"
        ),
    )


def _read_variable(
    dataset: netCDF4.Dataset | None,
    owner: str,
    name: str,
    *,
    time_index: int,
) -> FieldData:
    if dataset is None:
        return FieldData(
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
        return FieldData(
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
            np.ma.asarray(_read_variable_data(dataset, name, time_index=time_index)),
            dtype=np.float64,
        )
    except (IndexError, TypeError, ValueError) as exc:
        return FieldData(
            name=name,
            owner=owner,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=True,
            missing_inputs=(f"{owner}.{name}",),
            message=f"{owner}.{name} cannot be read: {exc}",
        )
    return FieldData(
        name=name,
        owner=owner,
        status="available",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        present=True,
        values=values,
        shape=tuple(int(value) for value in values.shape),
    )


def _field_expression(
    dataset: netCDF4.Dataset | None,
    owner: str,
    name: str,
    *,
    time_index: int,
) -> FieldData:
    if "+" not in name:
        return _read_variable(dataset, owner, name, time_index=time_index)

    parts = tuple(part.strip() for part in name.split("+"))
    loaded = [
        _read_variable(dataset, owner, part, time_index=time_index)
        for part in parts
    ]
    missing = tuple(item for field in loaded for item in field.missing_inputs)
    if any(field.status != "available" for field in loaded):
        messages = [field.message for field in loaded if field.message]
        return FieldData(
            name=name,
            owner=owner,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=all(field.present for field in loaded),
            missing_inputs=missing,
            message=(
                "; ".join(messages)
                if messages
                else f"{owner}.{name} inputs are not available"
            ),
        )

    arrays = [field.values for field in loaded]
    if any(array is None for array in arrays):
        return FieldData(
            name=name,
            owner=owner,
            status="not_available",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=True,
            missing_inputs=missing or (f"{owner}.{name}",),
            message=f"{owner}.{name} values are unavailable",
        )
    shapes = [tuple(array.shape) for array in arrays if array is not None]
    if len(set(shapes)) != 1:
        return FieldData(
            name=name,
            owner=owner,
            status="shape_mismatch",
            diagnostic_only=True,
            candidate_model_pass="not_applicable",
            present=True,
            shape=shapes[0],
            missing_inputs=missing,
            message=f"{owner}.{name} input shapes differ: {shapes}",
        )

    summed = np.zeros(shapes[0], dtype=np.float64)
    for array in arrays:
        summed = summed + np.asarray(array, dtype=np.float64)
    return FieldData(
        name=name,
        owner=owner,
        status="available",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        present=True,
        values=summed,
        shape=tuple(int(value) for value in summed.shape),
    )


def _horizontal_shape(
    candidate_end: netCDF4.Dataset,
    *,
    time_index: int,
) -> tuple[int, int] | None:
    for name in ("P", "PB", "MU", "HGT", "PH"):
        if name not in candidate_end.variables:
            continue
        try:
            values = np.asarray(_read_variable_data(candidate_end, name, time_index=time_index))
        except (IndexError, TypeError, ValueError):
            continue
        if values.ndim >= 2:
            return int(values.shape[-2]), int(values.shape[-1])
    return None


def _public_region(region: dict[str, Any]) -> dict[str, Any]:
    return _diag(**{key: value for key, value in region.items() if key != "mask"})


def _target_mask(candidate_end: netCDF4.Dataset, *, time_index: int) -> tuple[np.ndarray | None, dict[str, Any]]:
    shape = _horizontal_shape(candidate_end, time_index=time_index)
    region = _infer_target_region(candidate_end, shape)
    mask = region.get("mask")
    if isinstance(mask, np.ndarray):
        return np.asarray(mask, dtype=bool), _public_region(region)
    return None, _public_region(region)


def _read_p(
    dataset: netCDF4.Dataset | None,
    owner: str,
    *,
    time_index: int,
) -> tuple[np.ndarray | None, dict[str, Any]]:
    field = _read_variable(dataset, owner, "P", time_index=time_index)
    if field.status != "available" or field.values is None:
        return None, _diag(
            status=field.status,
            missing_inputs=list(field.missing_inputs),
            message=field.message,
        )
    if field.values.ndim < 3:
        return None, _diag(
            status="not_available",
            shape=field.shape,
            message=f"{owner}.P must have a vertical level axis",
        )
    return field.values, _diag(status="available", shape=field.shape)


def _value_at(values: np.ndarray, level: int, j: int, i: int) -> tuple[str, float | None, str | None]:
    if values.ndim == 2:
        if j < 0 or j >= values.shape[-2] or i < 0 or i >= values.shape[-1]:
            return "out_of_range", None, "horizontal index is outside the field"
        return _scalar_status(values[j, i])
    if values.ndim == 3:
        if level < 0 or level >= values.shape[-3]:
            return "out_of_range", None, "vertical level is outside the field"
        if j < 0 or j >= values.shape[-2] or i < 0 or i >= values.shape[-1]:
            return "out_of_range", None, "horizontal index is outside the field"
        return _scalar_status(values[level, j, i])
    return "not_available", None, f"unsupported field rank {values.ndim}"


def _scalar_status(value: Any) -> tuple[str, float | None, str | None]:
    numeric = _finite_float(value)
    if numeric is None:
        return "not_finite", None, "value is missing, masked, or non-finite"
    return "available", numeric, None


def _sample_field(field: FieldData, level: int, j: int, i: int) -> dict[str, Any]:
    if field.status != "available" or field.values is None:
        return _diag(
            status=field.status,
            value=None,
            present=field.present,
            shape=field.shape,
            missing_inputs=list(field.missing_inputs),
            message=field.message,
        )
    status, value, message = _value_at(field.values, level, j, i)
    return _diag(
        status=status,
        value=value,
        shape=tuple(int(item) for item in field.values.shape),
        level_axis="none_2d" if field.values.ndim == 2 else "vertical",
        message=message,
    )


def _diff_record(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
    left_value = _finite_float(left.get("value"))
    right_value = _finite_float(right.get("value"))
    if left_value is None or right_value is None:
        return _diag(
            status="not_available",
            value=None,
            message="one or both values are unavailable",
        )
    return _diag(status="available", value=left_value - right_value)


def _field_level_record(
    name: str,
    candidate_field: FieldData,
    reference_field: FieldData,
    source_field: FieldData,
    *,
    level: int,
    j: int,
    i: int,
) -> dict[str, Any]:
    candidate = _sample_field(candidate_field, level, j, i)
    reference = _sample_field(reference_field, level, j, i)
    source = _sample_field(source_field, level, j, i)
    return _diag(
        field=name,
        level=level,
        i=i,
        j=j,
        candidate=candidate,
        reference=reference,
        source=source,
        candidate_minus_reference=_diff_record(candidate, reference),
        candidate_minus_source=_diff_record(candidate, source),
        reference_minus_source=_diff_record(reference, source),
    )


def _rank_target_columns(
    candidate_end: netCDF4.Dataset,
    reference_end: netCDF4.Dataset | None,
    mask_2d: np.ndarray | None,
    region_report: dict[str, Any],
    *,
    levels: tuple[int, ...],
    max_columns: int,
    time_index: int,
) -> dict[str, Any]:
    candidate_metadata = _diag(**_candidate_metadata(candidate_end))
    if mask_2d is None:
        return _diag(
            status="not_available",
            ranking_method="none",
            levels=list(levels),
            max_columns=max_columns,
            target_region=region_report,
            candidate_metadata=candidate_metadata,
            selected_columns=[],
            message="target-region metadata is not available",
        )

    target_indices = [(int(j), int(i)) for j, i in np.argwhere(mask_2d)]
    if not target_indices:
        return _diag(
            status="not_available",
            ranking_method="none",
            levels=list(levels),
            max_columns=max_columns,
            target_region=region_report,
            candidate_metadata=candidate_metadata,
            selected_columns=[],
            message="target-region mask has no selected columns",
        )

    candidate_p, candidate_p_status = _read_p(
        candidate_end,
        "candidate_end",
        time_index=time_index,
    )
    if reference_end is None:
        selected = [
            _diag(
                rank=rank + 1,
                i=i,
                j=j,
                ranking_score=None,
                max_abs_low_level_p_error=None,
                mean_abs_low_level_p_error=None,
                low_level_p_errors=[],
                selection_basis="target_region_order_no_reference",
            )
            for rank, (j, i) in enumerate(target_indices[:max_columns])
        ]
        return _diag(
            status="selected_without_reference",
            ranking_method="target_region_row_major",
            levels=list(levels),
            max_columns=max_columns,
            target_region=region_report,
            candidate_metadata=candidate_metadata,
            candidate_p_status=candidate_p_status,
            reference_p_status=_diag(
                status="not_available",
                message="reference-end file was not supplied",
            ),
            selected_columns=selected,
        )

    reference_p, reference_p_status = _read_p(
        reference_end,
        "reference_end",
        time_index=time_index,
    )
    if candidate_p is None or reference_p is None:
        return _diag(
            status="not_available",
            ranking_method="low_level_p_abs_error",
            levels=list(levels),
            max_columns=max_columns,
            target_region=region_report,
            candidate_metadata=candidate_metadata,
            candidate_p_status=candidate_p_status,
            reference_p_status=reference_p_status,
            selected_columns=[],
            message="P is unavailable for candidate or reference",
        )
    if candidate_p.shape != reference_p.shape:
        return _diag(
            status="not_available",
            ranking_method="low_level_p_abs_error",
            levels=list(levels),
            max_columns=max_columns,
            target_region=region_report,
            candidate_metadata=candidate_metadata,
            candidate_p_status=candidate_p_status,
            reference_p_status=reference_p_status,
            selected_columns=[],
            message=(
                f"P shape mismatch: candidate {candidate_p.shape} != "
                f"reference {reference_p.shape}"
            ),
        )

    ranked: list[dict[str, Any]] = []
    for j, i in target_indices:
        per_level: list[dict[str, Any]] = []
        abs_errors: list[float] = []
        for level in levels:
            if level < 0 or level >= candidate_p.shape[-3]:
                per_level.append(
                    _diag(
                        status="out_of_range",
                        level=level,
                        candidate_p=None,
                        reference_p=None,
                        candidate_minus_reference=None,
                    )
                )
                continue
            cand_status, cand_value, cand_message = _value_at(candidate_p, level, j, i)
            ref_status, ref_value, ref_message = _value_at(reference_p, level, j, i)
            diff = (
                cand_value - ref_value
                if cand_value is not None and ref_value is not None
                else None
            )
            abs_error = None if diff is None else abs(diff)
            if abs_error is not None:
                abs_errors.append(abs_error)
            per_level.append(
                _diag(
                    status="available" if diff is not None else "not_available",
                    level=level,
                    candidate_p=cand_value,
                    reference_p=ref_value,
                    candidate_minus_reference=diff,
                    abs_error=abs_error,
                    message=cand_message or ref_message,
                    candidate_status=cand_status,
                    reference_status=ref_status,
                )
            )
        if not abs_errors:
            score = None
            mean_abs = None
            max_abs = None
        else:
            max_abs = float(max(abs_errors))
            mean_abs = float(np.mean(abs_errors))
            score = max_abs
        ranked.append(
            _diag(
                rank=None,
                i=i,
                j=j,
                ranking_score=score,
                max_abs_low_level_p_error=max_abs,
                mean_abs_low_level_p_error=mean_abs,
                low_level_p_errors=per_level,
                selection_basis="low_level_p_abs_error",
            )
        )

    ranked.sort(
        key=lambda item: (
            -float(item["ranking_score"] or -math.inf),
            -float(item["mean_abs_low_level_p_error"] or -math.inf),
            int(item["j"]),
            int(item["i"]),
        )
    )
    selected: list[dict[str, Any]] = []
    for rank, item in enumerate(ranked[:max_columns], start=1):
        selected.append({**item, "rank": rank})

    return _diag(
        status="ranked_by_reference",
        ranking_method="low_level_p_abs_error",
        levels=list(levels),
        max_columns=max_columns,
        target_region=region_report,
        candidate_metadata=candidate_metadata,
        candidate_p_status=candidate_p_status,
        reference_p_status=reference_p_status,
        selected_columns=selected,
    )


def _fields_for_owner(
    dataset: netCDF4.Dataset | None,
    owner: str,
    *,
    time_index: int,
) -> dict[str, FieldData]:
    return {
        name: _field_expression(dataset, owner, name, time_index=time_index)
        for name in FIELD_NAMES
    }


def _source_term_status(
    inventory: dict[str, Any],
    term: str,
) -> dict[str, Any]:
    compute = inventory.get("compute_krosa_pressure")
    if not isinstance(compute, dict):
        return _diag(status="not_available", term=term)
    terms = compute.get("formula_terms")
    if not isinstance(terms, dict):
        return _diag(status="not_available", term=term)
    item = terms.get(term)
    if isinstance(item, dict):
        return item
    return _diag(status="not_available", term=term)


def _read_coeff_array(
    dataset: netCDF4.Dataset | None,
    owner: str,
    name: str,
    *,
    time_index: int,
) -> tuple[np.ndarray | None, dict[str, Any]]:
    field = _read_variable(dataset, owner, name, time_index=time_index)
    if field.status != "available" or field.values is None:
        return None, _diag(
            status=field.status,
            missing_inputs=list(field.missing_inputs),
            message=field.message,
        )
    if field.values.ndim != 1:
        return None, _diag(
            status="not_available",
            shape=field.shape,
            message=f"{owner}.{name} must be a 1D coefficient vector after Time selection",
        )
    return field.values, _diag(status="available", shape=field.shape)


def _sample_full_coeff_pair(
    dataset: netCDF4.Dataset | None,
    owner: str,
    level: int,
    *,
    time_index: int,
) -> dict[str, Any]:
    c3f, c3f_status = _read_coeff_array(dataset, owner, "C3F", time_index=time_index)
    c4f, c4f_status = _read_coeff_array(dataset, owner, "C4F", time_index=time_index)
    missing = []
    if c3f is None:
        missing.append("C3F")
    if c4f is None:
        missing.append("C4F")
    if missing:
        return _diag(
            status="not_available",
            missing_names=missing,
            C3F=c3f_status,
            C4F=c4f_status,
        )
    if level < 0 or level + 1 >= c3f.shape[0] or level + 1 >= c4f.shape[0]:
        return _diag(
            status="out_of_range",
            missing_names=[],
            level=level,
            C3F_shape=tuple(int(value) for value in c3f.shape),
            C4F_shape=tuple(int(value) for value in c4f.shape),
        )
    return _diag(
        status="available",
        values={
            "C3F_k": float(c3f[level]),
            "C3F_k_plus_1": float(c3f[level + 1]),
            "C4F_k": float(c4f[level]),
            "C4F_k_plus_1": float(c4f[level + 1]),
        },
    )


def _sample_mass_coeff_pair(
    dataset: netCDF4.Dataset | None,
    owner: str,
    level: int,
    *,
    time_index: int,
) -> dict[str, Any]:
    c3h, c3h_status = _read_coeff_array(dataset, owner, "C3H", time_index=time_index)
    c4h, c4h_status = _read_coeff_array(dataset, owner, "C4H", time_index=time_index)
    missing = []
    if c3h is None:
        missing.append("C3H")
    if c4h is None:
        missing.append("C4H")
    if missing:
        return _diag(
            status="not_available",
            missing_names=missing,
            C3H=c3h_status,
            C4H=c4h_status,
        )
    if level < 0 or level >= c3h.shape[0] or level >= c4h.shape[0]:
        return _diag(
            status="out_of_range",
            missing_names=[],
            level=level,
            C3H_shape=tuple(int(value) for value in c3h.shape),
            C4H_shape=tuple(int(value) for value in c4h.shape),
        )
    return _diag(
        status="available",
        values={
            "C3H_k": float(c3h[level]),
            "C4H_k": float(c4h[level]),
        },
    )


def _read_p_top(
    dataset: netCDF4.Dataset | None,
    owner: str,
    *,
    time_index: int,
) -> dict[str, Any]:
    if dataset is None:
        return _diag(
            status="not_available",
            missing_names=["P_TOP"],
            message=f"{owner} dataset was not supplied",
        )
    if "P_TOP" in dataset.ncattrs():
        value = _finite_float(_netcdf_attr_to_json(dataset.getncattr("P_TOP")))
        return _diag(
            status="available" if value is not None else "not_finite",
            value=value,
            source="global_attribute",
        )
    if "P_TOP" not in dataset.variables:
        return _diag(
            status="not_available",
            missing_names=["P_TOP"],
            message=f"{owner}.P_TOP is not present as a global attribute or variable",
        )
    try:
        raw = _read_variable_data(dataset, "P_TOP", time_index=time_index)
    except (IndexError, TypeError, ValueError) as exc:
        return _diag(
            status="not_available",
            missing_names=["P_TOP"],
            message=f"{owner}.P_TOP cannot be read: {exc}",
        )
    array = np.asarray(raw, dtype=np.float64)
    if array.shape not in {(), (1,)}:
        return _diag(
            status="not_available",
            missing_names=[],
            shape=tuple(int(value) for value in array.shape),
            message=f"{owner}.P_TOP is not scalar after Time selection",
        )
    value = _finite_float(array.reshape(-1)[0])
    return _diag(
        status="available" if value is not None else "not_finite",
        value=value,
        source="variable",
    )


def _sample_alb(
    dataset: netCDF4.Dataset | None,
    owner: str,
    level: int,
    j: int,
    i: int,
    *,
    time_index: int,
) -> dict[str, Any]:
    field = _read_variable(dataset, owner, "ALB", time_index=time_index)
    return _sample_field(field, level, j, i)


def _read_theta_metadata(dataset: netCDF4.Dataset | None, owner: str) -> dict[str, Any]:
    if dataset is None:
        return _diag(
            status="not_available",
            missing_names=["theta_m"],
            message=f"{owner} dataset was not supplied",
        )
    attrs = {
        name: _netcdf_attr_to_json(dataset.getncattr(name))
        for name in THETA_METADATA_ATTRS
        if name in dataset.ncattrs()
    }
    if attrs:
        return _diag(status="available_metadata", attrs=attrs)
    variable_names = [name for name in ("theta_m", "THETA_M") if name in dataset.variables]
    if variable_names:
        return _diag(
            status="available_variable",
            variables=variable_names,
            message="theta_m is represented as a NetCDF variable, but this probe does not evaluate it as a pressure coefficient",
        )
    return _diag(
        status="not_available",
        missing_names=["theta_m"],
        message="theta_m is not available as explicit NetCDF metadata",
    )


def _coefficient_records(
    candidate_end: netCDF4.Dataset,
    reference_end: netCDF4.Dataset | None,
    source_start: netCDF4.Dataset | None,
    source_inventory: dict[str, Any],
    *,
    level: int,
    j: int,
    i: int,
    time_index: int,
    source_time_index: int,
) -> dict[str, Any]:
    return {
        "C3F/C4F": _diag(
            term="C3F/C4F",
            source_inventory=_source_term_status(source_inventory, "C3F/C4F"),
            candidate=_sample_full_coeff_pair(
                candidate_end,
                "candidate_end",
                level,
                time_index=time_index,
            ),
            reference=_sample_full_coeff_pair(
                reference_end,
                "reference_end",
                level,
                time_index=time_index,
            ),
            source=_sample_full_coeff_pair(
                source_start,
                "source_start",
                level,
                time_index=source_time_index,
            ),
        ),
        "C3H/C4H": _diag(
            term="C3H/C4H",
            source_inventory=_source_term_status(source_inventory, "C3H/C4H"),
            candidate=_sample_mass_coeff_pair(
                candidate_end,
                "candidate_end",
                level,
                time_index=time_index,
            ),
            reference=_sample_mass_coeff_pair(
                reference_end,
                "reference_end",
                level,
                time_index=time_index,
            ),
            source=_sample_mass_coeff_pair(
                source_start,
                "source_start",
                level,
                time_index=source_time_index,
            ),
        ),
        "P_TOP": _diag(
            term="P_TOP",
            source_inventory=_source_term_status(source_inventory, "P_TOP"),
            candidate=_read_p_top(candidate_end, "candidate_end", time_index=time_index),
            reference=_read_p_top(reference_end, "reference_end", time_index=time_index),
            source=_read_p_top(source_start, "source_start", time_index=source_time_index),
        ),
        "ALB": _diag(
            term="ALB",
            source_inventory=_source_term_status(source_inventory, "ALB"),
            candidate=_sample_alb(
                candidate_end,
                "candidate_end",
                level,
                j,
                i,
                time_index=time_index,
            ),
            reference=_sample_alb(
                reference_end,
                "reference_end",
                level,
                j,
                i,
                time_index=time_index,
            ),
            source=_sample_alb(
                source_start,
                "source_start",
                level,
                j,
                i,
                time_index=source_time_index,
            ),
        ),
        "theta_m": _diag(
            term="theta_m",
            source_inventory=_source_term_status(source_inventory, "theta_m"),
            candidate=_read_theta_metadata(candidate_end, "candidate_end"),
            reference=_read_theta_metadata(reference_end, "reference_end"),
            source=_read_theta_metadata(source_start, "source_start"),
        ),
    }


def _column_formula_terms(
    candidate_end: netCDF4.Dataset,
    reference_end: netCDF4.Dataset | None,
    source_start: netCDF4.Dataset | None,
    source_inventory: dict[str, Any],
    selected_columns: list[dict[str, Any]],
    *,
    levels: tuple[int, ...],
    time_index: int,
    source_time_index: int,
) -> dict[str, Any]:
    candidate_fields = _fields_for_owner(
        candidate_end,
        "candidate_end",
        time_index=time_index,
    )
    reference_fields = _fields_for_owner(
        reference_end,
        "reference_end",
        time_index=time_index,
    )
    source_fields = _fields_for_owner(
        source_start,
        "source_start",
        time_index=source_time_index,
    )

    columns: list[dict[str, Any]] = []
    for column in selected_columns:
        i = int(column["i"])
        j = int(column["j"])
        level_records = []
        for level in levels:
            fields = {
                name: _field_level_record(
                    name,
                    candidate_fields[name],
                    reference_fields[name],
                    source_fields[name],
                    level=level,
                    j=j,
                    i=i,
                )
                for name in FIELD_NAMES
            }
            level_records.append(
                _diag(
                    level=level,
                    i=i,
                    j=j,
                    fields=fields,
                    coefficient_terms=_coefficient_records(
                        candidate_end,
                        reference_end,
                        source_start,
                        source_inventory,
                        level=level,
                        j=j,
                        i=i,
                        time_index=time_index,
                        source_time_index=source_time_index,
                    ),
                )
            )
        columns.append(
            _diag(
                rank=column.get("rank"),
                i=i,
                j=j,
                ranking_score=column.get("ranking_score"),
                max_abs_low_level_p_error=column.get(
                    "max_abs_low_level_p_error"
                ),
                mean_abs_low_level_p_error=column.get(
                    "mean_abs_low_level_p_error"
                ),
                low_level_p_errors=column.get("low_level_p_errors", []),
                levels=level_records,
            )
        )

    required_status = {
        name: {
            "candidate_end": candidate_fields[name].status
            if name in candidate_fields
            else "not_available",
            "reference_end": reference_fields[name].status
            if name in reference_fields
            else "not_available",
            "source_start": source_fields[name].status
            if name in source_fields
            else "not_available",
        }
        for name in REQUIRED_FORMULA_INPUTS
    }
    return _diag(
        status="computed",
        levels=list(levels),
        field_names=list(FIELD_NAMES),
        coefficient_terms=list(COEFFICIENT_TERMS),
        required_formula_input_status=required_status,
        selected_column_count=len(columns),
        columns=columns,
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


def _diff_value(field_record: dict[str, Any], diff_name: str) -> float | None:
    diff = field_record.get(diff_name)
    if isinstance(diff, dict):
        return _finite_float(diff.get("value"))
    return None


def _field_value(field_record: dict[str, Any], owner: str) -> float | None:
    owner_record = field_record.get(owner)
    if isinstance(owner_record, dict):
        return _finite_float(owner_record.get("value"))
    return None


def _candidate_coeff_missing(column_terms: dict[str, Any]) -> list[str]:
    missing: set[str] = set()
    columns = column_terms.get("columns")
    if not isinstance(columns, list):
        return []
    for column in columns:
        if not isinstance(column, dict):
            continue
        for level in column.get("levels", []):
            if not isinstance(level, dict):
                continue
            coeffs = level.get("coefficient_terms")
            if not isinstance(coeffs, dict):
                continue
            for term in COEFFICIENT_TERMS:
                record = coeffs.get(term)
                candidate = record.get("candidate") if isinstance(record, dict) else None
                if isinstance(candidate, dict) and candidate.get("status") not in {
                    "available",
                    "available_metadata",
                    "available_variable",
                }:
                    missing.add(term)
    return sorted(missing)


def _risk_flags(
    selection: dict[str, Any],
    terms: dict[str, Any],
    source_inventory: dict[str, Any],
    formula_input_summary: dict[str, Any],
    *,
    reference_supplied: bool,
    source_supplied: bool,
) -> list[RiskFlag]:
    flags: list[RiskFlag] = []
    if not reference_supplied:
        _add_flag(
            flags,
            "reference_end_missing_ranked_by_target_order",
            "warning",
            "reference-end file was not supplied, so target columns are not ranked by P error",
            {"selection_status": selection.get("status")},
        )
    if not source_supplied:
        _add_flag(
            flags,
            "source_start_missing",
            "info",
            "source-start file was not supplied; candidate-source and reference-source diffs are unavailable",
            {"selected_column_count": len(selection.get("selected_columns", []))},
        )
    if selection.get("status") == "not_available":
        _add_flag(
            flags,
            "target_region_or_ranking_not_available",
            "warning",
            "target-region metadata or P ranking inputs are unavailable",
            {"message": selection.get("message"), "target_region": selection.get("target_region")},
        )

    for column in terms.get("columns", []):
        if not isinstance(column, dict):
            continue
        for level in column.get("levels", []):
            if not isinstance(level, dict):
                continue
            fields = level.get("fields")
            if not isinstance(fields, dict):
                continue
            p = fields.get("P")
            if not isinstance(p, dict):
                continue
            p_diff = _diff_value(p, "candidate_minus_reference")
            p_candidate_source = _diff_value(p, "candidate_minus_source")
            if p_diff is not None and abs(p_diff) >= PRESSURE_ERROR_WARNING_PA:
                _add_flag(
                    flags,
                    "low_level_column_p_error_large",
                    "warning",
                    "selected target column has large low-level perturbation-P error",
                    {
                        "rank": column.get("rank"),
                        "i": column.get("i"),
                        "j": column.get("j"),
                        "level": level.get("level"),
                        "candidate_minus_reference_p": p_diff,
                    },
                )

            companion_diffs = {
                name: _diff_value(fields.get(name, {}), "candidate_minus_reference")
                for name in ("PB", "MU+MUB", "PH+PHB")
            }
            t_diff = _diff_value(fields.get("T", {}), "candidate_minus_reference")
            q_diff = _diff_value(
                fields.get("QVAPOR", {}),
                "candidate_minus_reference",
            )
            if (
                p_diff is not None
                and abs(p_diff) >= PRESSURE_ERROR_WARNING_PA
                and all(
                    value is not None and abs(value) <= WEAK_COMPANION_WARNING_ABS
                    for value in companion_diffs.values()
                )
            ):
                _add_flag(
                    flags,
                    "p_bias_with_weak_column_companion_terms",
                    "warning",
                    "per-column P error is large while base-state companion formula terms differ weakly",
                    {
                        "rank": column.get("rank"),
                        "i": column.get("i"),
                        "j": column.get("j"),
                        "level": level.get("level"),
                        "candidate_minus_reference_p": p_diff,
                        "candidate_minus_reference_t": t_diff,
                        "candidate_minus_reference_qvapor": q_diff,
                        **{
                            f"candidate_minus_reference_{name}": value
                            for name, value in companion_diffs.items()
                        },
                    },
                )
            if p_diff is not None and p_candidate_source is not None and abs(p_diff - p_candidate_source) <= max(50.0, 0.1 * abs(p_diff)):
                _add_flag(
                    flags,
                    "column_candidate_source_delta_matches_p_error",
                    "info",
                    "candidate-source P delta closely matches candidate-reference P error in a selected column",
                    {
                        "rank": column.get("rank"),
                        "i": column.get("i"),
                        "j": column.get("j"),
                        "level": level.get("level"),
                        "candidate_minus_reference_p": p_diff,
                        "candidate_minus_source_p": p_candidate_source,
                    },
                )

    missing_coefficients = _candidate_coeff_missing(terms)
    if missing_coefficients:
        _add_flag(
            flags,
            "coefficient_terms_missing_from_candidate_netcdf",
            "info",
            "one or more pressure-formula coefficient terms are not explicitly available in candidate NetCDF",
            {"missing_candidate_coefficient_terms": missing_coefficients},
        )

    compute = source_inventory.get("compute_krosa_pressure")
    missing_terms = []
    if isinstance(compute, dict):
        raw_missing = compute.get("missing_formula_terms")
        if isinstance(raw_missing, list):
            missing_terms = raw_missing
    if missing_terms:
        _add_flag(
            flags,
            "source_formula_terms_missing",
            "warning",
            "static source inventory did not find all expected compute_krosa_pressure formula terms",
            {"missing_formula_terms": missing_terms},
        )
    if isinstance(compute, dict) and compute.get("status") == "not_available":
        _add_flag(
            flags,
            "compute_krosa_pressure_not_found",
            "warning",
            "static source inventory did not find compute_krosa_pressure",
            {"source_root": source_inventory.get("source_root")},
        )

    prior_flags = formula_input_summary.get("risk_flags")
    if isinstance(prior_flags, list):
        warning_codes = [
            item.get("code")
            for item in prior_flags
            if isinstance(item, dict)
            and item.get("severity") in {"warning", "error", "critical"}
        ]
        if warning_codes:
            _add_flag(
                flags,
                "formula_input_json_prior_risks_present",
                "warning",
                "D59 formula/input JSON already reports pressure formula risks",
                {"prior_risk_codes": warning_codes},
            )
    return flags


def _parse_levels(value: str) -> tuple[int, ...]:
    levels: list[int] = []
    for token in value.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        try:
            level = int(stripped)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"invalid level {stripped!r}; expected comma-separated integers"
            ) from exc
        if level < 0:
            raise argparse.ArgumentTypeError("levels must be non-negative")
        levels.append(level)
    if not levels:
        raise argparse.ArgumentTypeError("at least one level is required")
    return tuple(levels)


def audit_pressure_formula_column_probe(
    candidate_end_path: Path,
    *,
    reference_end_path: Path | None = None,
    source_start_path: Path | None = None,
    formula_input_json_path: Path | None = None,
    source_root: Path = Path("."),
    levels: Iterable[int] = DEFAULT_LEVELS,
    max_columns: int = DEFAULT_MAX_COLUMNS,
    time_index: int = -1,
    source_time_index: int = -1,
) -> PressureFormulaColumnProbeAudit:
    if max_columns <= 0:
        raise ValueError("--max-columns must be positive")
    level_tuple = tuple(int(level) for level in levels)
    if any(level < 0 for level in level_tuple):
        raise ValueError("levels must be non-negative")
    if not level_tuple:
        raise ValueError("at least one level is required")

    formula_input_payload = _load_json(
        formula_input_json_path,
        "formula-input",
    )
    formula_input_summary = _formula_input_json_summary(
        formula_input_json_path,
        formula_input_payload,
    )
    source_inventory = _d59_source_formula_inventory(source_root)

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
        mask_2d, region_report = _target_mask(candidate_end, time_index=time_index)
        selection = _rank_target_columns(
            candidate_end,
            reference_end,
            mask_2d,
            region_report,
            levels=level_tuple,
            max_columns=max_columns,
            time_index=time_index,
        )
        terms = _column_formula_terms(
            candidate_end,
            reference_end,
            source_start,
            source_inventory,
            selection.get("selected_columns", []),
            levels=level_tuple,
            time_index=time_index,
            source_time_index=source_time_index,
        )

    selection = _diag(
        **selection,
        prior_formula_input_json=formula_input_summary,
    )
    flags = _risk_flags(
        selection,
        terms,
        source_inventory,
        formula_input_summary,
        reference_supplied=reference_end_path is not None,
        source_supplied=source_start_path is not None,
    )
    summary = _diag(
        status="computed_with_flags" if flags else "computed",
        selected_column_count=len(selection.get("selected_columns", [])),
        level_count=len(level_tuple),
        max_columns=max_columns,
        risk_flag_count=len(flags),
        risk_flag_codes=[flag.code for flag in flags],
        strict_gate_status="not_evaluated",
        message=(
            "diagnostic-only per-column pressure formula-term probe; no "
            "candidate NetCDF values are created, changed, tuned, patched, or selected"
        ),
    )
    return PressureFormulaColumnProbeAudit(
        candidate_end=str(candidate_end_path),
        reference_end=None if reference_end_path is None else str(reference_end_path),
        source_start=None if source_start_path is None else str(source_start_path),
        formula_input_json=(
            None if formula_input_json_path is None else str(formula_input_json_path)
        ),
        status=summary["status"],
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        column_selection=selection,
        column_formula_terms=terms,
        source_formula_inventory=source_inventory,
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


def report_to_dict(report: PressureFormulaColumnProbeAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: PressureFormulaColumnProbeAudit,
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
        help="Optional WRF reference-end NetCDF file for diagnostic ranking/comparison",
    )
    parser.add_argument(
        "--source-start",
        type=Path,
        help="Optional source/start NetCDF file for evolution diagnostics",
    )
    parser.add_argument(
        "--formula-input-json",
        type=Path,
        help="Optional D59 audit_pressure_refresh_formula_inputs.py JSON report",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=Path("."),
        help="Repository/source root containing pressure-refresh C++ sources",
    )
    parser.add_argument(
        "--levels",
        type=_parse_levels,
        default=DEFAULT_LEVELS,
        help="Comma-separated low levels to probe, default 0,1,2,3,4",
    )
    parser.add_argument(
        "--max-columns",
        type=int,
        default=DEFAULT_MAX_COLUMNS,
        help="Maximum target columns to report",
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
        report = audit_pressure_formula_column_probe(
            args.candidate_end,
            reference_end_path=args.reference_end,
            source_start_path=args.source_start,
            formula_input_json_path=args.formula_input_json,
            source_root=args.source_root,
            levels=args.levels,
            max_columns=args.max_columns,
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
