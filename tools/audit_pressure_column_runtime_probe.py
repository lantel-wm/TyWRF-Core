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
EXPECTED_PHASES = ("post_static_refresh", "post_pressure_refresh")
FORMULA_TERMS = ("ALB", "C3F", "C4F", "C3H", "C4H", "P_TOP", "theta_m")
VALUE_KEY_TO_FIELD = {
    "P_PLUS_PB": "P+PB",
    "MU_PLUS_MUB": "MU+MUB",
    "PH_PLUS_PHB": "PH+PHB",
}
LARGE_P_DROP_PA = 500.0


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


def _coerce_float(value: str) -> float | None:
    try:
        numeric = float(value)
    except ValueError:
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


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
    for suffix in REQUIRED_ATTRS:
        name = f"{ATTR_PREFIX}{suffix}"
        if name not in dataset.ncattrs():
            missing.append(name)
            attrs[suffix.lower()] = None
            continue
        attrs[suffix.lower()] = _attr_to_json(dataset.getncattr(name))
    return attrs, missing


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


def _build_risk_flags(
    missing_attrs: list[str],
    parse_errors: list[str],
    record_count: dict[str, Any],
    declared_phases: list[str],
    observed_phases: list[str],
    not_available: list[str],
    records: list[dict[str, Any]],
    observations: list[dict[str, Any]],
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
    if unavailable_terms:
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


def audit_pressure_column_runtime_probe(candidate: Path) -> dict[str, Any]:
    candidate = Path(candidate)
    with netCDF4.Dataset(candidate) as dataset:
        attrs, missing_attrs = _read_probe_attrs(dataset)

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

    risk_flags = _build_risk_flags(
        missing_attrs,
        parse_errors,
        record_count,
        declared_phases,
        observed_phases,
        not_available,
        records,
        observations,
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
        "summary": _diag(
            column_count=len(columns),
            level_count=len(levels),
            phase_count=len(phases),
            parsed_record_count=len(records),
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
