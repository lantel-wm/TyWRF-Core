#!/usr/bin/env python
"""Audit D61 pressure-column runtime probe NetCDF attributes.

This tool is diagnostic-only. It reads candidate NetCDF metadata, structures
runtime pressure observations, and never creates, modifies, selects, patches,
or tunes candidate fields.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

import netCDF4
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if __package__ in {None, ""} and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


ATTR_PREFIX = "TYWRF_PRESSURE_COLUMN_PROBE_"
FORMULA_ATTR_PREFIX = "TYWRF_PRESSURE_FORMULA_OBSERVATION_"
REQUIRED_ATTRS = (
    "VERSION",
    "ENABLED",
    "EVIDENCE_ONLY",
    "INDEX_BASE",
    "COLUMN_COUNT",
    "LEVEL_COUNT",
    "PHASE_COUNT",
    "RECORD_COUNT",
    "COLUMNS",
    "LEVELS",
    "PHASES",
    "FIELDS",
    "NOT_AVAILABLE",
    "VALUES",
)
OPTIONAL_PROBE_ATTRS = ("FORMULA_OBSERVATION_ENABLED",)
FORMULA_REQUIRED_ATTRS = (
    "VERSION",
    "ENABLED",
    "EVIDENCE_ONLY",
    "INDEX_BASE",
    "REQUEST_COUNT",
    "RECORD_COUNT",
    "VALID_COUNT",
    "INVALID_COUNT",
    "OUT_OF_BOUNDS_COUNT",
    "OUTSIDE_TARGET_REGION_COUNT",
    "FIELDS",
    "VALUES",
)
EXPECTED_PHASES = ("post_static_refresh", "post_pressure_refresh")
FORMULA_TERMS = ("ALB", "C3F", "C4F", "C3H", "C4H", "P_TOP", "theta_m")
FORMULA_ID_FIELDS = {"status", "valid", "i", "j", "k"}
FORMULA_RECORDED_STATUS = "recorded"
FORMULA_OUT_OF_BOUNDS_STATUS = "request_out_of_bounds"
FORMULA_OUTSIDE_TARGET_REGION_STATUS = "request_outside_target_region"
FORMULA_SKIPPED_STATUSES = {
    FORMULA_OUT_OF_BOUNDS_STATUS,
    FORMULA_OUTSIDE_TARGET_REGION_STATUS,
}
FORMULA_VALUE_FIELDS = (
    "mu_total",
    "pfu",
    "pfd",
    "phm",
    "log_ratio",
    "phi_lower",
    "phi_upper",
    "delta_phi",
    "ALB",
    "PB",
    "theta",
    "alpha_total",
    "alpha_perturbation",
    "alpha_from_wrf_branch",
    "pressure_base",
    "total_pressure",
    "perturbation_pressure_pa",
)
VALUE_KEY_TO_FIELD = {
    "P_PLUS_PB": "P+PB",
    "MU_PLUS_MUB": "MU+MUB",
    "PH_PLUS_PHB": "PH+PHB",
}
LARGE_P_DROP_PA = 500.0
PRESSURE_MATCH_TOLERANCE_PA = 1.0e-3
LARGE_FRACTIONAL_TOTAL_PRESSURE_INCREASE = 0.01
MOST_RECORDS_FRACTION = 0.8
SENSITIVITY_DRIVER_FIELDS = (
    "theta",
    "delta_phi",
    "mu_total",
    "phm",
    "pfu",
    "pfd",
    "alpha_total",
    "pressure_base",
    "PB",
)


def _diag(**values: Any) -> dict[str, Any]:
    return {
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        **values,
    }


def _risk(
    code: str,
    severity: str,
    message: str,
    evidence: dict[str, Any],
) -> dict[str, Any]:
    return _diag(code=code, severity=severity, message=message, evidence=evidence)


def _attr_to_json(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray):
        return _strict_json_value(value.tolist())
    if isinstance(value, np.generic):
        return _strict_json_value(value.item())
    return value


def _strict_json_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _strict_json_value(item) for key, item in value.items()}
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


def _coerce_int(value: Any) -> int | None:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            numeric = float(text)
        except ValueError:
            return None
        if numeric.is_integer():
            return int(numeric)
    return None


def _coerce_float(value: Any) -> float | None:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bool):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _coerce_bool(value: Any) -> bool | None:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in {0, 1}:
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "y"}:
            return True
        if normalized in {"false", "0", "no", "n"}:
            return False
    return None


def _safe_divide(numerator: float, denominator: float) -> float | None:
    if denominator == 0.0:
        return None
    result = numerator / denominator
    return result if math.isfinite(result) else None


def _csv_items(value: Any) -> list[str]:
    if value is None:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _parse_columns(value: Any) -> tuple[list[dict[str, int]], list[str]]:
    errors: list[str] = []
    columns: list[dict[str, int]] = []
    if value is None or str(value).strip() == "":
        return columns, errors
    for index, raw_column in enumerate(str(value).split(";"), start=1):
        raw_column = raw_column.strip()
        if not raw_column:
            continue
        parts = [part.strip() for part in raw_column.split(",")]
        if len(parts) != 2:
            errors.append(f"column {index} must be i,j: {raw_column}")
            continue
        try:
            i_value = int(parts[0])
            j_value = int(parts[1])
        except ValueError:
            errors.append(f"column {index} has non-integer coordinates: {raw_column}")
            continue
        columns.append({"i": i_value, "j": j_value})
    return columns, errors


def _parse_levels(value: Any) -> tuple[list[int], list[str]]:
    errors: list[str] = []
    levels: list[int] = []
    for index, raw_level in enumerate(_csv_items(value), start=1):
        try:
            levels.append(int(raw_level))
        except ValueError:
            errors.append(f"level {index} is not an integer: {raw_level}")
    return levels, errors


def _read_probe_attrs(dataset: netCDF4.Dataset) -> tuple[dict[str, Any], list[str]]:
    attrs: dict[str, Any] = {}
    missing: list[str] = []
    dataset_attrs = set(dataset.ncattrs())
    for suffix in REQUIRED_ATTRS:
        name = f"{ATTR_PREFIX}{suffix}"
        if name not in dataset_attrs:
            missing.append(name)
            attrs[suffix.lower()] = None
            continue
        attrs[suffix.lower()] = _attr_to_json(dataset.getncattr(name))
    for suffix in OPTIONAL_PROBE_ATTRS:
        name = f"{ATTR_PREFIX}{suffix}"
        if name in dataset_attrs:
            attrs[suffix.lower()] = _attr_to_json(dataset.getncattr(name))
    return attrs, missing


def _read_formula_attrs(
    dataset: netCDF4.Dataset,
) -> tuple[dict[str, Any], list[str], list[str]]:
    attrs: dict[str, Any] = {}
    missing: list[str] = []
    dataset_attrs = set(dataset.ncattrs())
    present = sorted(
        name for name in dataset_attrs if name.startswith(FORMULA_ATTR_PREFIX)
    )
    for suffix in FORMULA_REQUIRED_ATTRS:
        name = f"{FORMULA_ATTR_PREFIX}{suffix}"
        key = suffix.lower()
        if name not in dataset_attrs:
            attrs[key] = None
            if present:
                missing.append(name)
            continue
        attrs[key] = _attr_to_json(dataset.getncattr(name))
    return attrs, missing, present


def _parse_record_fields(raw_record: str) -> tuple[dict[str, str], list[str]]:
    fields: dict[str, str] = {}
    errors: list[str] = []
    for raw_part in raw_record.split(";"):
        part = raw_part.strip()
        if not part:
            continue
        if "=" not in part:
            errors.append(f"record part lacks '=': {part}")
            continue
        key, value = part.split("=", 1)
        key = key.strip()
        if not key:
            errors.append(f"record part has empty key: {part}")
            continue
        fields[key] = value.strip()
    return fields, errors


def _parse_probe_values(
    raw_values: Any,
    declared_fields: list[str],
) -> tuple[list[dict[str, Any]], list[str], int]:
    if raw_values is None:
        return [], [], 0

    records: list[dict[str, Any]] = []
    errors: list[str] = []
    raw_text = str(raw_values).strip()
    if raw_text == "":
        return records, errors, 0

    raw_records = [item.strip() for item in raw_text.split("|")]
    for record_index, raw_record in enumerate(raw_records, start=1):
        if not raw_record:
            errors.append(f"record {record_index} is empty")
            continue

        fields, field_errors = _parse_record_fields(raw_record)
        errors.extend(f"record {record_index}: {error}" for error in field_errors)

        phase = fields.get("phase")
        i_value = _coerce_int(fields.get("i"))
        j_value = _coerce_int(fields.get("j"))
        k_value = _coerce_int(fields.get("k"))
        if not phase:
            errors.append(f"record {record_index}: missing phase")
        if i_value is None:
            errors.append(f"record {record_index}: missing or invalid i")
        if j_value is None:
            errors.append(f"record {record_index}: missing or invalid j")
        if k_value is None:
            errors.append(f"record {record_index}: missing or invalid k")

        values: dict[str, float] = {}
        for key, value in fields.items():
            if key in {"phase", "i", "j", "k"}:
                continue
            field_name = VALUE_KEY_TO_FIELD.get(key, key)
            numeric = _coerce_float(value)
            if numeric is None:
                errors.append(
                    f"record {record_index}: field {field_name} is not a finite float"
                )
                continue
            values[field_name] = numeric

        if declared_fields:
            missing_fields = [field for field in declared_fields if field not in values]
            if missing_fields:
                errors.append(
                    f"record {record_index}: missing declared fields "
                    f"{','.join(missing_fields)}"
                )

        if phase is None or i_value is None or j_value is None or k_value is None:
            continue

        records.append(
            {
                "record_index": record_index,
                "phase": phase,
                "i": i_value,
                "j": j_value,
                "k": k_value,
                "values": values,
            }
        )

    return records, errors, len(raw_records)


def _parse_formula_values(
    raw_values: Any,
    declared_fields: list[str],
) -> tuple[list[dict[str, Any]], list[str], list[dict[str, Any]], int]:
    if raw_values is None:
        return [], [], [], 0

    records: list[dict[str, Any]] = []
    errors: list[str] = []
    term_errors: list[dict[str, Any]] = []
    raw_text = str(raw_values).strip()
    if raw_text == "":
        return records, errors, term_errors, 0

    expected_fields = declared_fields or [*FORMULA_ID_FIELDS, *FORMULA_VALUE_FIELDS]
    expected_numeric_fields = [
        field for field in expected_fields if field not in FORMULA_ID_FIELDS
    ]
    raw_records = [item.strip() for item in raw_text.split("|")]
    for record_index, raw_record in enumerate(raw_records, start=1):
        if not raw_record:
            errors.append(f"record {record_index} is empty")
            continue

        fields, field_errors = _parse_record_fields(raw_record)
        errors.extend(f"record {record_index}: {error}" for error in field_errors)

        status = fields.get("status")
        valid = _coerce_bool(fields.get("valid"))
        i_value = _coerce_int(fields.get("i"))
        j_value = _coerce_int(fields.get("j"))
        k_value = _coerce_int(fields.get("k"))
        if status is None or status == "":
            errors.append(f"record {record_index}: missing status")
        if valid is None:
            errors.append(f"record {record_index}: missing or invalid valid")
        if i_value is None:
            errors.append(f"record {record_index}: missing or invalid i")
        if j_value is None:
            errors.append(f"record {record_index}: missing or invalid j")
        if k_value is None:
            errors.append(f"record {record_index}: missing or invalid k")

        status_kind = _formula_status_kind(status)
        requires_formula_terms = status_kind in {"recorded", "invalid", "unknown"}
        values: dict[str, float] = {}
        for field in expected_numeric_fields:
            if field not in fields:
                if requires_formula_terms:
                    term_errors.append(
                        {
                            "record_index": record_index,
                            "field": field,
                            "reason": "missing",
                        }
                    )
                continue
            numeric = _coerce_float(fields[field])
            if numeric is None:
                if requires_formula_terms:
                    term_errors.append(
                        {
                            "record_index": record_index,
                            "field": field,
                            "raw_value": fields[field],
                            "reason": "not_finite_float",
                        }
                    )
                continue
            values[field] = numeric

        extra_numeric_fields = sorted(
            set(fields) - set(expected_numeric_fields) - FORMULA_ID_FIELDS
        )
        for field in extra_numeric_fields:
            numeric = _coerce_float(fields[field])
            if numeric is not None:
                values[field] = numeric

        if i_value is None or j_value is None or k_value is None:
            continue

        records.append(
            {
                "record_index": record_index,
                "status": status,
                "valid": valid,
                "i": i_value,
                "j": j_value,
                "k": k_value,
                "values": values,
            }
        )

    return records, errors, term_errors, len(raw_records)


def _column_key(column: dict[str, int]) -> tuple[int, int]:
    return int(column["i"]), int(column["j"])


def _observed_columns(records: list[dict[str, Any]]) -> list[dict[str, int]]:
    seen: set[tuple[int, int]] = set()
    columns: list[dict[str, int]] = []
    for record in records:
        key = (record["i"], record["j"])
        if key in seen:
            continue
        seen.add(key)
        columns.append({"i": record["i"], "j": record["j"]})
    return columns


def _observed_levels(records: list[dict[str, Any]]) -> list[int]:
    levels: list[int] = []
    seen: set[int] = set()
    for record in records:
        level = int(record["k"])
        if level in seen:
            continue
        seen.add(level)
        levels.append(level)
    return levels


def _observed_phases(records: list[dict[str, Any]]) -> list[str]:
    phases: list[str] = []
    seen: set[str] = set()
    for record in records:
        phase = str(record["phase"])
        if phase in seen:
            continue
        seen.add(phase)
        phases.append(phase)
    return phases


def _ordered_records(
    declared_columns: list[dict[str, int]],
    declared_levels: list[int],
    declared_phases: list[str],
    records: list[dict[str, Any]],
) -> tuple[list[dict[str, int]], list[int], list[str]]:
    columns = declared_columns or _observed_columns(records)
    levels = declared_levels or _observed_levels(records)
    phases = declared_phases or _observed_phases(records)
    return columns, levels, phases


def _build_observations(
    columns: list[dict[str, int]],
    levels: list[int],
    phases: list[str],
    records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[str]]:
    errors: list[str] = []
    grouped: dict[tuple[int, int, int, str], dict[str, Any]] = {}
    for record in records:
        key = (record["i"], record["j"], record["k"], record["phase"])
        if key in grouped:
            errors.append(
                "duplicate record for "
                f"i={record['i']},j={record['j']},k={record['k']},phase={record['phase']}"
            )
            continue
        grouped[key] = record

    output_columns: list[dict[str, Any]] = []
    for column in columns:
        i_value, j_value = _column_key(column)
        level_entries: list[dict[str, Any]] = []
        for level in levels:
            phase_values: dict[str, dict[str, float]] = {}
            for phase in phases:
                record = grouped.get((i_value, j_value, level, phase))
                if record is not None:
                    phase_values[phase] = record["values"]

            deltas: dict[str, dict[str, float]] = {}
            before = phase_values.get("post_static_refresh")
            after = phase_values.get("post_pressure_refresh")
            if before is not None and after is not None:
                common_fields = sorted(set(before) & set(after))
                deltas["post_pressure_refresh_minus_post_static_refresh"] = {
                    field: after[field] - before[field] for field in common_fields
                }

            level_entries.append(
                {
                    "k": level,
                    "phases": phase_values,
                    "deltas": deltas,
                    "diagnostic_only": True,
                    "candidate_model_pass": "not_applicable",
                }
            )
        output_columns.append(
            {
                "i": i_value,
                "j": j_value,
                "levels": level_entries,
                "diagnostic_only": True,
                "candidate_model_pass": "not_applicable",
            }
        )
    return output_columns, errors


def _record_count_summary(
    declared_count: int | None,
    raw_record_count: int,
    parsed_count: int,
    columns: list[dict[str, int]],
    levels: list[int],
    phases: list[str],
) -> dict[str, Any]:
    expected = None
    if columns and levels and phases:
        expected = len(columns) * len(levels) * len(phases)
    return _diag(
        declared=declared_count,
        raw_record_segments=raw_record_count,
        parsed=parsed_count,
        expected_from_declared_shape=expected,
        declared_matches_parsed=(
            None if declared_count is None else declared_count == parsed_count
        ),
        declared_matches_expected=(
            None if declared_count is None or expected is None else declared_count == expected
        ),
        expected_matches_parsed=(None if expected is None else expected == parsed_count),
    )


def _phase_missing(
    declared_phases: list[str],
    observed_phases: list[str],
) -> list[str]:
    available = set(declared_phases) | set(observed_phases)
    return [phase for phase in EXPECTED_PHASES if phase not in available]


def _negative_post_pressure(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    negatives: list[dict[str, Any]] = []
    for record in records:
        if record["phase"] != "post_pressure_refresh":
            continue
        p_value = record["values"].get("P")
        if p_value is not None and p_value < 0.0:
            negatives.append(
                {
                    "i": record["i"],
                    "j": record["j"],
                    "k": record["k"],
                    "P": p_value,
                }
            )
    return negatives


def _large_p_drops(observations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    drops: list[dict[str, Any]] = []
    for column in observations:
        for level in column["levels"]:
            delta = level["deltas"].get(
                "post_pressure_refresh_minus_post_static_refresh", {}
            ).get("P")
            if delta is not None and delta <= -LARGE_P_DROP_PA:
                drops.append(
                    {
                        "i": column["i"],
                        "j": column["j"],
                        "k": level["k"],
                        "delta_P": delta,
                    }
                )
    drops.sort(key=lambda item: item["delta_P"])
    return drops


def _formula_status_kind(status: Any) -> str:
    if status == FORMULA_RECORDED_STATUS:
        return "recorded"
    if status == FORMULA_OUT_OF_BOUNDS_STATUS:
        return "out_of_bounds"
    if status == FORMULA_OUTSIDE_TARGET_REGION_STATUS:
        return "outside_target_region"
    if isinstance(status, str) and status.startswith("invalid_"):
        return "invalid"
    return "unknown"


def _formula_record_status_problem(record: dict[str, Any]) -> str | None:
    status_kind = _formula_status_kind(record.get("status"))
    valid = record.get("valid")
    if status_kind == "recorded":
        return None if valid is True else "recorded_without_valid_flag"
    if status_kind in {"out_of_bounds", "outside_target_region"}:
        return None if valid is False else f"{status_kind}_with_valid_flag"
    if status_kind == "invalid":
        return "invalid_formula_status"
    return "unknown_formula_status"


def _formula_record_count_summary(
    attrs: dict[str, Any],
    raw_record_count: int,
    parsed_records: list[dict[str, Any]],
) -> dict[str, Any]:
    declared_count = _coerce_int(attrs.get("record_count"))
    declared_valid_count = _coerce_int(attrs.get("valid_count"))
    declared_invalid_count = _coerce_int(attrs.get("invalid_count"))
    declared_out_of_bounds_count = _coerce_int(attrs.get("out_of_bounds_count"))
    declared_outside_target_region_count = _coerce_int(
        attrs.get("outside_target_region_count")
    )
    parsed_count = len(parsed_records)
    valid_count = sum(
        1
        for record in parsed_records
        if record.get("valid") is True
        and _formula_status_kind(record.get("status")) == "recorded"
    )
    out_of_bounds_count = sum(
        1
        for record in parsed_records
        if _formula_status_kind(record.get("status")) == "out_of_bounds"
    )
    outside_target_region_count = sum(
        1
        for record in parsed_records
        if _formula_status_kind(record.get("status")) == "outside_target_region"
    )
    invalid_count = sum(
        1
        for record in parsed_records
        if _formula_status_kind(record.get("status")) in {"invalid", "unknown"}
    )
    return _diag(
        request_count=_coerce_int(attrs.get("request_count")),
        declared=declared_count,
        raw_record_segments=raw_record_count,
        parsed=parsed_count,
        declared_matches_parsed=(
            None if declared_count is None else declared_count == parsed_count
        ),
        valid_declared=declared_valid_count,
        valid_observed=valid_count,
        valid_matches_declared=(
            None
            if declared_valid_count is None
            else declared_valid_count == valid_count
        ),
        invalid_declared=declared_invalid_count,
        invalid_observed=invalid_count,
        invalid_matches_declared=(
            None
            if declared_invalid_count is None
            else declared_invalid_count == invalid_count
        ),
        out_of_bounds_declared=declared_out_of_bounds_count,
        out_of_bounds_observed=out_of_bounds_count,
        out_of_bounds_matches_declared=(
            None
            if declared_out_of_bounds_count is None
            else declared_out_of_bounds_count == out_of_bounds_count
        ),
        outside_target_region_declared=declared_outside_target_region_count,
        outside_target_region_observed=outside_target_region_count,
        outside_target_region_matches_declared=(
            None
            if declared_outside_target_region_count is None
            else declared_outside_target_region_count == outside_target_region_count
        ),
    )


def _build_formula_records_by_column_level(
    records: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        key = f"{record['i']},{record['j']},{record['k']}"
        grouped.setdefault(key, []).append(record)
    return grouped


def _formula_supplied_fields(records: list[dict[str, Any]]) -> list[str]:
    supplied: set[str] = set()
    for record in records:
        supplied.update(field for field in record["values"] if field in FORMULA_VALUE_FIELDS)
    return sorted(supplied)


def _probe_pressure_maps(
    records: list[dict[str, Any]],
) -> tuple[dict[tuple[int, int, int], float], dict[tuple[int, int, int], float]]:
    static: dict[tuple[int, int, int], float] = {}
    post: dict[tuple[int, int, int], float] = {}
    for record in records:
        p_value = record["values"].get("P")
        if p_value is None:
            continue
        key = (record["i"], record["j"], record["k"])
        if record["phase"] == "post_static_refresh":
            static[key] = p_value
        elif record["phase"] == "post_pressure_refresh":
            post[key] = p_value
    return static, post


def _range_summary(values: list[float]) -> dict[str, Any]:
    if not values:
        return _diag(count=0, min=None, max=None, mean=None)
    return _diag(
        count=len(values),
        min=min(values),
        max=max(values),
        mean=sum(values) / len(values),
    )


def _build_formula_sensitivity_record(record: dict[str, Any]) -> dict[str, Any] | None:
    values = record["values"]
    total_pressure = values.get("total_pressure")
    pb_value = values.get("PB")
    if total_pressure is None or pb_value is None:
        return None

    gap_to_pb = pb_value - total_pressure
    pressure_change_needed = max(gap_to_pb, 0.0)
    fractional_increase_needed = (
        _safe_divide(pressure_change_needed, total_pressure)
        if total_pressure > 0.0
        else None
    )
    pressure_base = values.get("pressure_base")
    pressure_base_fractional_proxy = None
    pressure_base_delta_proxy = None
    if (
        pressure_base is not None
        and total_pressure > 0.0
        and pressure_base != 0.0
    ):
        pressure_base_fractional_proxy = fractional_increase_needed
        if fractional_increase_needed is not None:
            pressure_base_delta_proxy = pressure_base * fractional_increase_needed

    drivers = {
        field: values[field]
        for field in SENSITIVITY_DRIVER_FIELDS
        if field in values
    }
    pfu = values.get("pfu")
    pfd = values.get("pfd")
    phm = values.get("phm")

    return _diag(
        i=record["i"],
        j=record["j"],
        k=record["k"],
        total_pressure_pa=total_pressure,
        PB=pb_value,
        total_pressure_gap_to_pb_pa=gap_to_pb,
        relative_total_pressure_gap_to_pb=_safe_divide(gap_to_pb, abs(pb_value)),
        pressure_change_needed_to_make_perturbation_nonnegative_pa=(
            pressure_change_needed
        ),
        approximate_fractional_total_pressure_increase_needed=(
            fractional_increase_needed
        ),
        available_drivers=drivers,
        exact_budget_directional_sensitivities=_diag(
            perturbation_pressure_pa={
                "with_respect_to_total_pressure": 1.0,
                "with_respect_to_PB": -1.0,
            },
            total_pressure_gap_to_pb_pa={
                "with_respect_to_total_pressure": -1.0,
                "with_respect_to_PB": 1.0,
            },
        ),
        diagnostic_driver_hints=_diag(
            total_pressure_direct=_diag(
                required_delta_pa=pressure_change_needed,
                fractional_change_needed=fractional_increase_needed,
            ),
            PB_threshold_direct=_diag(
                required_delta_pa=-pressure_change_needed,
                fractional_change_needed=(
                    _safe_divide(-pressure_change_needed, pb_value)
                    if pb_value != 0.0
                    else None
                ),
            ),
            pressure_base_multiplicative_proxy=_diag(
                available=pressure_base is not None,
                assumption=(
                    "Proxy only: if recorded pressure_base scales total_pressure "
                    "linearly while other formula terms stay fixed."
                ),
                required_delta_pa=pressure_base_delta_proxy,
                fractional_change_needed=pressure_base_fractional_proxy,
            ),
            sign_only_hints=_diag(
                theta=(
                    "larger theta can raise total pressure only if this record's "
                    "formula branch uses theta as a positive thermodynamic factor"
                ),
                alpha_total=(
                    "smaller alpha_total can raise total pressure only if this "
                    "record's formula branch uses alpha_total as a denominator"
                ),
                delta_phi=(
                    "delta_phi is formula-coupled; this report does not infer a "
                    "standalone pressure direction from the recorded terms"
                ),
                mu_total=(
                    "mu_total is formula-coupled; this report does not infer a "
                    "standalone pressure direction from the recorded terms"
                ),
                phm_pfu_pfd=(
                    "phm, pfu, and pfd are interpolation drivers; inspect staging "
                    "if the pressure gap is large"
                ),
            ),
            bounded_recorded_ratios=_diag(
                pfu_plus_pfd=None if pfu is None or pfd is None else pfu + pfd,
                delta_phi_over_abs_phm=(
                    None
                    if values.get("delta_phi") is None or phm in {None, 0.0}
                    else _safe_divide(values["delta_phi"], abs(phm))
                ),
            ),
        ),
    )


def _build_formula_sensitivity_summary(
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    gaps = [record["total_pressure_gap_to_pb_pa"] for record in records]
    needed = [
        record["pressure_change_needed_to_make_perturbation_nonnegative_pa"]
        for record in records
    ]
    fractional_needed = [
        record["approximate_fractional_total_pressure_increase_needed"]
        for record in records
        if record["approximate_fractional_total_pressure_increase_needed"] is not None
    ]
    records_requiring_increase = [
        record
        for record in records
        if record["pressure_change_needed_to_make_perturbation_nonnegative_pa"] > 0.0
    ]
    records_requiring_large_increase = [
        record
        for record in records
        if (
            record["approximate_fractional_total_pressure_increase_needed"]
            is not None
            and record["approximate_fractional_total_pressure_increase_needed"]
            >= LARGE_FRACTIONAL_TOTAL_PRESSURE_INCREASE
        )
    ]
    large_fraction = (
        len(records_requiring_large_increase) / len(records) if records else 0.0
    )
    return _diag(
        record_count=len(records),
        total_pressure_gap_to_pb_pa=_range_summary(gaps),
        pressure_change_needed_to_make_perturbation_nonnegative_pa=(
            _range_summary(needed)
        ),
        approximate_fractional_total_pressure_increase_needed=(
            _range_summary(fractional_needed)
        ),
        records_requiring_total_pressure_increase_count=len(
            records_requiring_increase
        ),
        records_requiring_large_total_pressure_increase_count=len(
            records_requiring_large_increase
        ),
        large_fractional_total_pressure_increase_threshold=(
            LARGE_FRACTIONAL_TOTAL_PRESSURE_INCREASE
        ),
        records_requiring_large_total_pressure_increase_fraction=large_fraction,
        most_records_require_large_total_pressure_increase=(
            bool(records)
            and large_fraction >= MOST_RECORDS_FRACTION
        ),
        all_records_require_total_pressure_increase=(
            bool(records) and len(records_requiring_increase) == len(records)
        ),
    )


def _correlate_formula_observation(
    formula_records: list[dict[str, Any]],
    probe_records: list[dict[str, Any]],
) -> dict[str, Any]:
    static_p, post_p = _probe_pressure_maps(probe_records)
    matched: list[dict[str, Any]] = []
    mismatches: list[dict[str, Any]] = []
    missing_probe_post: list[dict[str, Any]] = []
    pressure_budget_records: list[dict[str, Any]] = []
    formula_sensitivity_records: list[dict[str, Any]] = []
    total_pressure_below_pb_records: list[dict[str, Any]] = []
    pressure_drop_explained_records: list[dict[str, Any]] = []
    for record in formula_records:
        formula_p = record["values"].get("perturbation_pressure_pa")
        if formula_p is None:
            continue
        key = (record["i"], record["j"], record["k"])
        if key not in post_p:
            missing_probe_post.append(
                {
                    "i": record["i"],
                    "j": record["j"],
                    "k": record["k"],
                    "formula_perturbation_pressure_pa": formula_p,
                }
            )
            continue
        probe_p = post_p[key]
        diff = formula_p - probe_p
        probe_delta = None
        if key in static_p:
            probe_delta = probe_p - static_p[key]
        item = {
            "i": record["i"],
            "j": record["j"],
            "k": record["k"],
            "formula_perturbation_pressure_pa": formula_p,
            "probe_post_pressure_refresh_p": probe_p,
            "difference_pa": diff,
            "probe_delta_p": probe_delta,
            "matches_within_tolerance": abs(diff) <= PRESSURE_MATCH_TOLERANCE_PA,
        }
        matched.append(item)
        if not item["matches_within_tolerance"]:
            mismatches.append(item)

        total_pressure = record["values"].get("total_pressure")
        pb_value = record["values"].get("PB")
        if total_pressure is None or pb_value is None:
            continue
        total_minus_pb = total_pressure - pb_value
        perturbation_minus_total_minus_pb = formula_p - total_minus_pb
        post_matches_total_minus_pb = (
            abs(probe_p - total_minus_pb) <= PRESSURE_MATCH_TOLERANCE_PA
        )
        large_drop = probe_delta is not None and probe_delta <= -LARGE_P_DROP_PA
        drop_matches_total_minus_pb = (
            probe_delta is not None
            and key in static_p
            and abs(probe_delta - (total_minus_pb - static_p[key]))
            <= PRESSURE_MATCH_TOLERANCE_PA
        )
        drop_explained = bool(
            large_drop
            and total_pressure < pb_value
            and post_matches_total_minus_pb
            and drop_matches_total_minus_pb
        )
        budget = {
            "i": record["i"],
            "j": record["j"],
            "k": record["k"],
            "total_pressure_pa": total_pressure,
            "PB": pb_value,
            "perturbation_pressure_pa": formula_p,
            "probe_post_pressure_refresh_p": probe_p,
            "total_pressure_minus_pb_pa": total_minus_pb,
            "perturbation_minus_total_minus_pb_pa": (
                perturbation_minus_total_minus_pb
            ),
            "formula_total_pressure_below_pb": total_pressure < pb_value,
            "post_refresh_p_matches_total_minus_pb": post_matches_total_minus_pb,
            "probe_delta_p": probe_delta,
            "large_probe_drop": large_drop,
            "pressure_drop_expected_from_total_minus_pb_pa": (
                None if key not in static_p else total_minus_pb - static_p[key]
            ),
            "probe_delta_matches_total_minus_pb_minus_static_p": (
                drop_matches_total_minus_pb
            ),
            "large_drop_explained_by_formula_base_subtraction": drop_explained,
        }
        pressure_budget_records.append(budget)
        sensitivity_record = _build_formula_sensitivity_record(record)
        if sensitivity_record is not None:
            formula_sensitivity_records.append(sensitivity_record)
        if budget["formula_total_pressure_below_pb"]:
            total_pressure_below_pb_records.append(budget)
        if drop_explained:
            pressure_drop_explained_records.append(budget)
    correlation = _diag(
        match_tolerance_pa=PRESSURE_MATCH_TOLERANCE_PA,
        matched_records=matched,
        pressure_mismatches=mismatches,
        missing_post_pressure_probe_records=missing_probe_post,
    )
    if pressure_budget_records:
        correlation["pressure_budget"] = _diag(
            records=pressure_budget_records,
            total_pressure_below_pb_records=total_pressure_below_pb_records,
            pressure_drop_explained_by_base_subtraction_records=(
                pressure_drop_explained_records
            ),
        )
        if formula_sensitivity_records:
            correlation["pressure_budget"]["formula_sensitivity"] = _diag(
                records=formula_sensitivity_records,
                summary=_build_formula_sensitivity_summary(
                    formula_sensitivity_records
                ),
            )
    return correlation


def _formula_claims_enabled(
    probe_attrs: dict[str, Any],
    formula_attrs: dict[str, Any],
) -> bool:
    return (
        _coerce_bool(formula_attrs.get("enabled")) is True
        or _coerce_bool(probe_attrs.get("formula_observation_enabled")) is True
    )


def _build_risk_flags(
    missing_attrs: list[str],
    parse_errors: list[str],
    record_count: dict[str, Any],
    declared_phases: list[str],
    observed_phases: list[str],
    not_available: list[str],
    records: list[dict[str, Any]],
    observations: list[dict[str, Any]],
    formula_terms_supplied: bool,
) -> list[dict[str, Any]]:
    flags: list[dict[str, Any]] = []
    if missing_attrs:
        flags.append(
            _risk(
                "missing_probe_attrs",
                "error",
                "One or more pressure-column runtime probe attributes are missing.",
                {"missing_attrs": missing_attrs},
            )
        )

    if parse_errors:
        flags.append(
            _risk(
                "malformed_values",
                "error",
                "Pressure-column runtime probe attributes or records are malformed.",
                {"errors": parse_errors[:50], "error_count": len(parse_errors)},
            )
        )

    mismatch_keys = (
        "declared_matches_parsed",
        "declared_matches_expected",
        "expected_matches_parsed",
    )
    if any(record_count.get(key) is False for key in mismatch_keys):
        flags.append(
            _risk(
                "record_count_mismatch",
                "error",
                "Pressure-column runtime probe record counts are inconsistent.",
                {"record_count": record_count},
            )
        )

    missing_phases = _phase_missing(declared_phases, observed_phases)
    if missing_phases:
        flags.append(
            _risk(
                "missing_expected_phases",
                "error",
                "Expected before/after pressure-column probe phases are missing.",
                {
                    "expected_phases": list(EXPECTED_PHASES),
                    "missing_phases": missing_phases,
                    "declared_phases": declared_phases,
                    "observed_phases": observed_phases,
                },
            )
        )

    negative_records = _negative_post_pressure(records)
    if negative_records:
        flags.append(
            _risk(
                "post_pressure_refresh_p_negative",
                "error",
                "post_pressure_refresh contains negative perturbation pressure.",
                {
                    "count": len(negative_records),
                    "min_P": min(item["P"] for item in negative_records),
                    "examples": negative_records[:10],
                },
            )
        )

    drops = _large_p_drops(observations)
    if drops:
        flags.append(
            _risk(
                "large_p_drop_magnitude",
                "warning",
                "P drops by at least the configured threshold after pressure refresh.",
                {
                    "threshold_pa": LARGE_P_DROP_PA,
                    "count": len(drops),
                    "largest_drop": drops[0],
                    "examples": drops[:10],
                },
            )
        )

    unavailable_terms = [term for term in FORMULA_TERMS if term in set(not_available)]
    if unavailable_terms and formula_terms_supplied:
        flags.append(
            _risk(
                "probe_field_terms_unavailable",
                "warning",
                "NetCDF pressure-column probe fields still mark formula terms unavailable, but formula observation supplied internal terms separately.",
                {
                    "required_formula_terms": list(FORMULA_TERMS),
                    "unavailable_probe_fields": unavailable_terms,
                    "formula_observation_supplies_internal_terms": True,
                },
            )
        )
    elif unavailable_terms:
        flags.append(
            _risk(
                "formula_terms_unavailable",
                "warning",
                "Pressure refresh formula terms remain unavailable in runtime probe metadata.",
                {
                    "required_formula_terms": list(FORMULA_TERMS),
                    "unavailable_terms": unavailable_terms,
                },
            )
        )

    return flags


def _build_formula_risk_flags(
    *,
    formula_present: bool,
    formula_enabled_claims_true: bool,
    formula_missing_attrs: list[str],
    formula_parse_errors: list[str],
    formula_term_errors: list[dict[str, Any]],
    formula_record_count: dict[str, Any],
    formula_records: list[dict[str, Any]],
    formula_correlation: dict[str, Any],
    not_available: list[str],
) -> list[dict[str, Any]]:
    flags: list[dict[str, Any]] = []
    unavailable_terms = [term for term in FORMULA_TERMS if term in set(not_available)]
    if formula_enabled_claims_true and unavailable_terms and (
        not formula_present or formula_missing_attrs
    ):
        flags.append(
            _risk(
                "missing_formula_observation",
                "error",
                "Formula observation was claimed enabled while probe metadata still marks formula terms unavailable, but formula observation attrs are incomplete.",
                {
                    "formula_present": formula_present,
                    "missing_attrs": formula_missing_attrs,
                    "unavailable_probe_terms": unavailable_terms,
                },
            )
        )

    if formula_parse_errors:
        flags.append(
            _risk(
                "malformed_formula_observation_values",
                "error",
                "Pressure formula observation records are malformed.",
                {
                    "errors": formula_parse_errors[:50],
                    "error_count": len(formula_parse_errors),
                },
            )
        )

    mismatch_keys = (
        "declared_matches_parsed",
        "valid_matches_declared",
        "invalid_matches_declared",
        "out_of_bounds_matches_declared",
        "outside_target_region_matches_declared",
    )
    if any(formula_record_count.get(key) is False for key in mismatch_keys):
        flags.append(
            _risk(
                "formula_observation_count_mismatch",
                "error",
                "Pressure formula observation record counts are inconsistent.",
                {"record_count": formula_record_count},
            )
        )

    invalid_records = [
        {
            "i": record["i"],
            "j": record["j"],
            "k": record["k"],
            "status": record.get("status"),
            "valid": record.get("valid"),
        }
        for record in formula_records
        if _formula_record_status_problem(record) is not None
    ]
    if invalid_records:
        flags.append(
            _risk(
                "formula_observation_invalid_status",
                "warning",
                "Formula observation contains invalid or internally inconsistent statuses.",
                {
                    "count": len(invalid_records),
                    "examples": invalid_records[:20],
                },
            )
        )

    if formula_term_errors:
        flags.append(
            _risk(
                "formula_observation_nonfinite_or_invalid_terms",
                "error",
                "Formula observation contains missing or nonfinite formula terms.",
                {
                    "errors": formula_term_errors[:50],
                    "error_count": len(formula_term_errors),
                },
            )
        )

    pressure_mismatches = formula_correlation.get("pressure_mismatches", [])
    if pressure_mismatches:
        flags.append(
            _risk(
                "formula_pressure_probe_mismatch",
                "error",
                "Formula observation perturbation pressure differs from post_pressure_refresh probe P for the same i,j,k.",
                {
                    "match_tolerance_pa": PRESSURE_MATCH_TOLERANCE_PA,
                    "count": len(pressure_mismatches),
                    "examples": pressure_mismatches[:20],
                },
            )
        )

    pressure_budget = formula_correlation.get("pressure_budget", {})
    below_pb_records = pressure_budget.get("total_pressure_below_pb_records", [])
    if below_pb_records:
        flags.append(
            _risk(
                "formula_total_pressure_below_base_pressure",
                "warning",
                "Formula observation total pressure is below the base pressure, so perturbation pressure is negative by construction.",
                {
                    "match_tolerance_pa": PRESSURE_MATCH_TOLERANCE_PA,
                    "count": len(below_pb_records),
                    "examples": below_pb_records[:20],
                },
            )
        )

    explained_records = pressure_budget.get(
        "pressure_drop_explained_by_base_subtraction_records", []
    )
    if explained_records:
        flags.append(
            _risk(
                "formula_pressure_drop_explained_by_base_subtraction",
                "warning",
                "Large post-refresh probe P drops are consistent with formula total-pressure minus base-pressure subtraction.",
                {
                    "match_tolerance_pa": PRESSURE_MATCH_TOLERANCE_PA,
                    "large_drop_threshold_pa": LARGE_P_DROP_PA,
                    "count": len(explained_records),
                    "examples": explained_records[:20],
                },
            )
        )

    sensitivity_summary = (
        pressure_budget.get("formula_sensitivity", {})
        .get("summary", {})
    )
    if sensitivity_summary.get("most_records_require_large_total_pressure_increase"):
        flags.append(
            _risk(
                "formula_total_pressure_gap_requires_staging_diagnosis",
                "warning",
                "Most formula observations require a large total-pressure increase to make perturbation pressure nonnegative; inspect pressure-formula staging before changing the PB subtraction.",
                {
                    "record_count": sensitivity_summary.get("record_count"),
                    "large_fractional_total_pressure_increase_threshold": (
                        sensitivity_summary.get(
                            "large_fractional_total_pressure_increase_threshold"
                        )
                    ),
                    "records_requiring_large_total_pressure_increase_count": (
                        sensitivity_summary.get(
                            "records_requiring_large_total_pressure_increase_count"
                        )
                    ),
                    "records_requiring_large_total_pressure_increase_fraction": (
                        sensitivity_summary.get(
                            "records_requiring_large_total_pressure_increase_fraction"
                        )
                    ),
                    "fractional_increase_range": sensitivity_summary.get(
                        "approximate_fractional_total_pressure_increase_needed"
                    ),
                    "gap_to_pb_range_pa": sensitivity_summary.get(
                        "total_pressure_gap_to_pb_pa"
                    ),
                },
            )
        )

    return flags


def audit_pressure_column_runtime_probe(candidate: Path) -> dict[str, Any]:
    candidate = Path(candidate)
    with netCDF4.Dataset(candidate) as dataset:
        attrs, missing_attrs = _read_probe_attrs(dataset)
        formula_attrs, formula_missing_attrs, formula_present_attrs = (
            _read_formula_attrs(dataset)
        )

    declared_columns, column_errors = _parse_columns(attrs.get("columns"))
    declared_levels, level_errors = _parse_levels(attrs.get("levels"))
    declared_phases = _csv_items(attrs.get("phases"))
    declared_fields = _csv_items(attrs.get("fields"))
    not_available = _csv_items(attrs.get("not_available"))
    declared_record_count = _coerce_int(attrs.get("record_count"))

    records, value_errors, raw_record_count = _parse_probe_values(
        attrs.get("values"), declared_fields
    )
    parse_errors = [*column_errors, *level_errors, *value_errors]

    observed_columns = _observed_columns(records)
    observed_levels = _observed_levels(records)
    observed_phases = _observed_phases(records)
    columns, levels, phases = _ordered_records(
        declared_columns, declared_levels, declared_phases, records
    )
    observations, observation_errors = _build_observations(
        columns, levels, phases, records
    )
    parse_errors.extend(observation_errors)

    record_count = _record_count_summary(
        declared_record_count,
        raw_record_count,
        len(records),
        declared_columns,
        declared_levels,
        declared_phases,
    )

    formula_declared_fields = _csv_items(formula_attrs.get("fields"))
    (
        formula_records,
        formula_value_errors,
        formula_term_errors,
        formula_raw_record_count,
    ) = _parse_formula_values(formula_attrs.get("values"), formula_declared_fields)
    formula_record_count = _formula_record_count_summary(
        formula_attrs, formula_raw_record_count, formula_records
    )
    formula_records_by_column_level = _build_formula_records_by_column_level(
        formula_records
    )
    formula_supplied_fields = _formula_supplied_fields(formula_records)
    formula_correlation = _correlate_formula_observation(formula_records, records)
    formula_present = bool(formula_present_attrs)
    formula_enabled_claims_true = _formula_claims_enabled(attrs, formula_attrs)

    risk_flags = _build_risk_flags(
        missing_attrs,
        parse_errors,
        record_count,
        declared_phases,
        observed_phases,
        not_available,
        records,
        observations,
        bool(formula_supplied_fields),
    )
    risk_flags.extend(
        _build_formula_risk_flags(
            formula_present=formula_present,
            formula_enabled_claims_true=formula_enabled_claims_true,
            formula_missing_attrs=formula_missing_attrs,
            formula_parse_errors=formula_value_errors,
            formula_term_errors=formula_term_errors,
            formula_record_count=formula_record_count,
            formula_records=formula_records,
            formula_correlation=formula_correlation,
            not_available=not_available,
        )
    )

    return {
        "candidate": str(candidate),
        "status": "computed_with_flags" if risk_flags else "computed",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "usage_boundary": _diag(
            statement=(
                "This report is diagnostic evidence only; it is not a validation "
                "gate pass and must not be used as candidate generation input."
            ),
            can_serve_as_gate_pass=False,
            can_generate_candidate=False,
            reads_candidate_only=True,
        ),
        "attrs": _diag(
            version=attrs.get("version"),
            enabled=attrs.get("enabled"),
            evidence_only=attrs.get("evidence_only"),
            index_base=attrs.get("index_base"),
            missing_attrs=missing_attrs,
        ),
        "columns": _diag(declared=declared_columns, observed=observed_columns),
        "levels": _diag(declared=declared_levels, observed=observed_levels),
        "phases": _diag(declared=declared_phases, observed=observed_phases),
        "fields": _diag(declared=declared_fields, not_available=not_available),
        "record_count": record_count,
        "records": _diag(parsed=records),
        "observations": _diag(columns=observations),
        "formula_observation": _diag(
            present=formula_present,
            enabled_claims_true=formula_enabled_claims_true,
            attrs=_diag(
                version=formula_attrs.get("version"),
                enabled=formula_attrs.get("enabled"),
                evidence_only=formula_attrs.get("evidence_only"),
                index_base=formula_attrs.get("index_base"),
                present_attrs=formula_present_attrs,
                missing_attrs=formula_missing_attrs,
            ),
            fields=_diag(
                declared=formula_declared_fields,
                supplied=formula_supplied_fields,
            ),
            record_count=formula_record_count,
            records=_diag(parsed=formula_records),
            records_by_column_level=_diag(items=formula_records_by_column_level),
            correlation=formula_correlation,
        ),
        "summary": _diag(
            column_count=len(columns),
            level_count=len(levels),
            phase_count=len(phases),
            parsed_record_count=len(records),
            formula_observation_present=formula_present,
            formula_record_count=len(formula_records),
            risk_flag_count=len(risk_flags),
            risk_codes=[flag["code"] for flag in risk_flags],
        ),
        "risk_flags": risk_flags,
    }


def report_to_json(report: dict[str, Any], *, pretty: bool = False) -> str:
    return json.dumps(
        _strict_json_value(report),
        indent=2 if pretty else None,
        allow_nan=False,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path, help="Candidate NetCDF file")
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_pressure_column_runtime_probe(args.candidate)
    except OSError as exc:
        parser.error(str(exc))

    text = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
