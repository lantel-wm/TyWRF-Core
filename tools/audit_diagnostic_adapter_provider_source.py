#!/usr/bin/env python
"""Audit D77 diagnostic-adapter provider-source metadata in a NetCDF output."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
from typing import Any

import netCDF4
import numpy as np


PROVIDER_PREFIX = "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_"
SOURCE_STAGING_PREFIX = "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_"
SOURCE_CHILD_DELTA_PREFIX = "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_"
DELTA_FIELDS = ("PHB", "MUB", "HT", "PB", "T_INIT", "ALB")
DELTA_FIELD_SUFFIXES = (
    "COMPARED_VALUE_COUNT",
    "DIFFERING_VALUE_COUNT",
    "MAX_ABS_DIFF",
)

REQUIRED_PROVIDER_KEYS = (
    "VERSION",
    "SOURCE_ORIGIN",
    "TERRAIN_SOURCE",
    "TERRAIN_PROVENANCE",
    "PROVIDER_OK",
    "DIAGNOSTIC_ONLY",
    "GATE_CANDIDATE",
    "INTEGRATOR_OUTPUT",
    "WRITES_CANDIDATE",
    "USES_REFERENCE_END_TRUTH",
    "USES_DIRECT_P_SHORTCUT",
    "READS_DIRECT_P",
    "WROTE_PB",
    "WROTE_T_INIT",
    "WROTE_MUB",
    "WROTE_ALB",
    "WROTE_PHB",
    "PROVIDER_RECONSTRUCTED_PHB_NOT_WRF_REBALANCE_VALIDATED",
)
REQUIRED_STAGING_KEYS = (
    "OK",
    "DIAGNOSTIC_ONLY",
    "GATE_CANDIDATE",
    "INTEGRATOR_OUTPUT",
    "WRITES_CANDIDATE",
    "USES_REFERENCE_END_TRUTH",
    "USES_DIRECT_P_SHORTCUT",
)
REQUIRED_DELTA_KEYS = (
    "VERSION",
    "DIAGNOSTIC_ONLY",
    "GATE_CANDIDATE",
    "INTEGRATOR_OUTPUT",
    "WRITES_CANDIDATE",
    "WRITES_NETCDF",
    "VALUES_IDENTICAL",
    "COMPARED_VALUE_COUNT",
    "DIFFERING_VALUE_COUNT",
    "MAX_ABS_DIFF",
)


@dataclass(frozen=True)
class AuditResult:
    status: str
    path: str
    gate_evidence: bool
    provider_source: dict[str, Any]
    source_staging: dict[str, Any]
    source_child_delta: dict[str, Any]
    largest_delta_field: dict[str, Any] | None
    reasons: list[str]


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


def _coerce_bool(value: Any) -> bool | None:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and value in (0, 1):
        return bool(value)
    if isinstance(value, (float, np.floating)) and value in (0.0, 1.0):
        return bool(value)
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "f", "no", "n", "off"}:
        return False
    return None


def _coerce_number(value: Any) -> float | None:
    if isinstance(value, (bool, np.bool_)) or value is None:
        return None
    if isinstance(value, (int, float, np.integer, np.floating)):
        numeric = float(value)
        return numeric if math.isfinite(numeric) else None
    text = str(value).strip()
    if not text:
        return None
    try:
        numeric = float(text)
    except ValueError:
        return None
    return numeric if math.isfinite(numeric) else None


def _read_prefixed_attrs(
    dataset: netCDF4.Dataset,
    prefix: str,
) -> dict[str, Any]:
    return {
        name.removeprefix(prefix): _netcdf_attr_to_json(dataset.getncattr(name))
        for name in dataset.ncattrs()
        if name.startswith(prefix)
    }


def _missing(attrs: dict[str, Any], keys: tuple[str, ...]) -> list[str]:
    return [key for key in keys if key not in attrs]


def _bool_summary(
    attrs: dict[str, Any], keys: tuple[str, ...]
) -> dict[str, bool | None]:
    return {key.lower(): _coerce_bool(attrs.get(key)) for key in keys}


def _staging_counts(attrs: dict[str, Any]) -> dict[str, Any]:
    count_keys = tuple(key for key in attrs if key.endswith("_COUNT"))
    counts = {key.lower(): _coerce_number(attrs.get(key)) for key in sorted(count_keys)}
    return {
        "all_counts": counts,
        "input_count": counts.get("input_count"),
        "staged_count": counts.get("staged_count"),
        "missing_count": counts.get("missing_count"),
        "invalid_count": counts.get("invalid_count"),
        "staged_value_count": counts.get("staged_value_count"),
        "invalid_exposed_value_count": counts.get("invalid_exposed_value_count"),
        "exposed_region_count": counts.get("exposed_region_count"),
        "exposed_counts": {
            key: value
            for key, value in counts.items()
            if key.startswith("exposed_") or key == "exposed_region_count"
        },
        "staged_counts": {
            key: value for key, value in counts.items() if key.startswith("staged_")
        },
        "invalid_counts": {
            key: value for key, value in counts.items() if key.startswith("invalid_")
        },
    }


def _provider_summary(attrs: dict[str, Any]) -> dict[str, Any]:
    bools = _bool_summary(
        attrs,
        (
            "PROVIDER_OK",
            "DIAGNOSTIC_ONLY",
            "GATE_CANDIDATE",
            "INTEGRATOR_OUTPUT",
            "WRITES_CANDIDATE",
            "WRITES_NETCDF",
            "NO_CANDIDATE_WRITE",
            "TERRAIN_OVERRIDE_USED",
            "USES_REFERENCE_END_TRUTH",
            "NO_REFERENCE_END_TRUTH",
            "USES_DIRECT_P_SHORTCUT",
            "NO_DIRECT_P_SHORTCUT",
            "READS_DIRECT_P",
            "WROTE_PB",
            "WROTE_T_INIT",
            "WROTE_MUB",
            "WROTE_ALB",
            "WROTE_PHB",
            "PROVIDER_RECONSTRUCTED_PHB_NOT_WRF_REBALANCE_VALIDATED",
        ),
    )
    return {
        "present": bool(attrs),
        "version": attrs.get("VERSION"),
        "origin": attrs.get("ORIGIN"),
        "source_origin": attrs.get("SOURCE_ORIGIN"),
        "provider_source": attrs.get("PROVIDER_SOURCE"),
        "terrain_source": attrs.get("TERRAIN_SOURCE"),
        "terrain_provenance": attrs.get("TERRAIN_PROVENANCE"),
        "terrain_override_used": bools["terrain_override_used"],
        "ht_source": attrs.get("HT_SOURCE"),
        "ht_provenance": attrs.get("HT_PROVENANCE"),
        "result_message": attrs.get("RESULT_MESSAGE"),
        "provider_ok": bools["provider_ok"],
        "diagnostic_only": bools["diagnostic_only"],
        "gate_candidate": bools["gate_candidate"],
        "integrator_output": bools["integrator_output"],
        "writes_candidate": bools["writes_candidate"],
        "writes_netcdf": bools["writes_netcdf"],
        "no_candidate_write": bools["no_candidate_write"],
        "uses_reference_end_truth": bools["uses_reference_end_truth"],
        "no_reference_end_truth": bools["no_reference_end_truth"],
        "uses_direct_p_shortcut": bools["uses_direct_p_shortcut"],
        "no_direct_p_shortcut": bools["no_direct_p_shortcut"],
        "reads_direct_p": bools["reads_direct_p"],
        "wrote_pb": bools["wrote_pb"],
        "wrote_t_init": bools["wrote_t_init"],
        "wrote_mub": bools["wrote_mub"],
        "wrote_alb": bools["wrote_alb"],
        "wrote_phb": bools["wrote_phb"],
        "phb_diagnostic_only_marker": bools[
            "provider_reconstructed_phb_not_wrf_rebalance_validated"
        ],
        "missing_required_attrs": _missing(attrs, REQUIRED_PROVIDER_KEYS),
        "raw_attrs": attrs,
    }


def _source_staging_summary(attrs: dict[str, Any]) -> dict[str, Any]:
    bools = _bool_summary(
        attrs,
        (
            "OK",
            "DIAGNOSTIC_ONLY",
            "GATE_CANDIDATE",
            "INTEGRATOR_OUTPUT",
            "WRITES_CANDIDATE",
            "WRITES_NETCDF",
            "CANDIDATE_BUFFERS_PRESERVED",
            "OWNS_STAGING_BUFFERS",
            "ALLOCATED_BUFFERS",
            "USES_REFERENCE_END_TRUTH",
            "USES_DIRECT_P_SHORTCUT",
            "READS_DIRECT_P",
            "ALIASES_CHILD",
        ),
    )
    counts = _staging_counts(attrs)
    return {
        "present": bool(attrs),
        "version": attrs.get("VERSION"),
        "provider_kind": attrs.get("PROVIDER_KIND"),
        "source": attrs.get("SOURCE"),
        "disposition": attrs.get("DISPOSITION"),
        "source_shape": attrs.get("SOURCE_SHAPE"),
        "result_message": attrs.get("RESULT_MESSAGE"),
        "ok": bools["ok"],
        "diagnostic_only": bools["diagnostic_only"],
        "gate_candidate": bools["gate_candidate"],
        "integrator_output": bools["integrator_output"],
        "writes_candidate": bools["writes_candidate"],
        "writes_netcdf": bools["writes_netcdf"],
        "candidate_buffers_preserved": bools["candidate_buffers_preserved"],
        "owns_staging_buffers": bools["owns_staging_buffers"],
        "allocated_buffers": bools["allocated_buffers"],
        "uses_reference_end_truth": bools["uses_reference_end_truth"],
        "uses_direct_p_shortcut": bools["uses_direct_p_shortcut"],
        "reads_direct_p": bools["reads_direct_p"],
        "aliases_child": bools["aliases_child"],
        "counts": counts,
        "missing_required_attrs": _missing(attrs, REQUIRED_STAGING_KEYS),
        "raw_attrs": attrs,
    }


def _delta_field_summary(attrs: dict[str, Any], field: str) -> dict[str, Any] | None:
    prefix = f"{field}_"
    present_keys = [
        prefix + suffix for suffix in DELTA_FIELD_SUFFIXES if prefix + suffix in attrs
    ]
    if not present_keys:
        return None
    return {
        "compared_value_count": _coerce_number(
            attrs.get(prefix + "COMPARED_VALUE_COUNT")
        ),
        "differing_value_count": _coerce_number(
            attrs.get(prefix + "DIFFERING_VALUE_COUNT")
        ),
        "max_abs_diff": _coerce_number(attrs.get(prefix + "MAX_ABS_DIFF")),
        "missing_attrs": [
            prefix + suffix for suffix in DELTA_FIELD_SUFFIXES if prefix + suffix not in attrs
        ],
    }


def _source_child_delta_summary(attrs: dict[str, Any]) -> dict[str, Any]:
    bools = _bool_summary(
        attrs,
        (
            "DIAGNOSTIC_ONLY",
            "GATE_CANDIDATE",
            "INTEGRATOR_OUTPUT",
            "WRITES_CANDIDATE",
            "WRITES_NETCDF",
            "VALUES_IDENTICAL",
        ),
    )
    fields = {
        field: summary
        for field in DELTA_FIELDS
        if (summary := _delta_field_summary(attrs, field)) is not None
    }
    return {
        "present": bool(attrs),
        "version": attrs.get("VERSION"),
        "source": attrs.get("SOURCE"),
        "scope": attrs.get("SCOPE"),
        "fields_attr": attrs.get("FIELDS"),
        "diagnostic_only": bools["diagnostic_only"],
        "gate_candidate": bools["gate_candidate"],
        "integrator_output": bools["integrator_output"],
        "writes_candidate": bools["writes_candidate"],
        "writes_netcdf": bools["writes_netcdf"],
        "values_identical": bools["values_identical"],
        "compared_value_count": _coerce_number(attrs.get("COMPARED_VALUE_COUNT")),
        "differing_value_count": _coerce_number(attrs.get("DIFFERING_VALUE_COUNT")),
        "max_abs_diff": _coerce_number(attrs.get("MAX_ABS_DIFF")),
        "fields": fields,
        "missing_required_attrs": _missing(attrs, REQUIRED_DELTA_KEYS),
        "raw_attrs": attrs,
    }


def _largest_delta_field(
    delta: dict[str, Any],
) -> dict[str, Any] | None:
    largest: tuple[str, float] | None = None
    for field, summary in delta["fields"].items():
        value = summary.get("max_abs_diff")
        if value is None:
            continue
        if largest is None or value > largest[1]:
            largest = (field, value)
    if largest is None:
        return None
    return {"field": largest[0], "max_abs_diff": largest[1]}


def _require_bool(
    reasons: list[str],
    summary: dict[str, Any],
    section: str,
    key: str,
    expected: bool,
) -> None:
    value = summary.get(key)
    if value is None:
        reasons.append(f"{section}.{key} is missing or not a boolean")
    elif value is not expected:
        reasons.append(f"{section}.{key} is {value}, expected {expected}")


def _validate_delta(reasons: list[str], delta: dict[str, Any]) -> None:
    compared = delta["compared_value_count"]
    differing = delta["differing_value_count"]
    max_abs_diff = delta["max_abs_diff"]
    for key, value in (
        ("compared_value_count", compared),
        ("differing_value_count", differing),
        ("max_abs_diff", max_abs_diff),
    ):
        if value is None or value < 0.0:
            reasons.append(f"source_child_delta.{key} is missing or invalid")
    if compared is not None and compared <= 0.0:
        reasons.append("source_child_delta.compared_value_count must be positive")
    if compared is not None and differing is not None and differing > compared:
        reasons.append(
            "source_child_delta.differing_value_count exceeds compared_value_count"
        )

    for field in DELTA_FIELDS:
        field_summary = delta["fields"].get(field)
        if field_summary is None:
            reasons.append(f"source_child_delta.{field} field metadata is missing")
            continue
        if field_summary["missing_attrs"]:
            reasons.append(
                f"source_child_delta.{field} missing attrs: "
                + ", ".join(field_summary["missing_attrs"])
            )
        field_compared = field_summary["compared_value_count"]
        field_differing = field_summary["differing_value_count"]
        field_max = field_summary["max_abs_diff"]
        for key, value in (
            ("compared_value_count", field_compared),
            ("differing_value_count", field_differing),
            ("max_abs_diff", field_max),
        ):
            if value is None or value < 0.0:
                reasons.append(f"source_child_delta.{field}.{key} is missing or invalid")
        if (
            field_compared is not None
            and field_differing is not None
            and field_differing > field_compared
        ):
            reasons.append(
                f"source_child_delta.{field}.differing_value_count exceeds "
                "compared_value_count"
            )

    if len(delta["fields"]) == len(DELTA_FIELDS) and all(
        field["compared_value_count"] is not None
        and field["differing_value_count"] is not None
        and field["max_abs_diff"] is not None
        for field in delta["fields"].values()
    ):
        field_compared_sum = sum(
            field["compared_value_count"] for field in delta["fields"].values()
        )
        field_differing_sum = sum(
            field["differing_value_count"] for field in delta["fields"].values()
        )
        field_max = max(field["max_abs_diff"] for field in delta["fields"].values())
        if compared is not None and not math.isclose(compared, field_compared_sum):
            reasons.append(
                "source_child_delta.compared_value_count does not match per-field sum"
            )
        if differing is not None and not math.isclose(differing, field_differing_sum):
            reasons.append(
                "source_child_delta.differing_value_count does not match per-field sum"
            )
        if max_abs_diff is not None and not math.isclose(max_abs_diff, field_max):
            reasons.append(
                "source_child_delta.max_abs_diff does not match largest per-field max"
            )


def _validate_report(
    provider: dict[str, Any],
    staging: dict[str, Any],
    delta: dict[str, Any],
) -> list[str]:
    reasons: list[str] = []
    if provider["missing_required_attrs"]:
        reasons.append(
            "provider_source missing required attrs: "
            + ", ".join(provider["missing_required_attrs"])
        )
    if staging["missing_required_attrs"]:
        reasons.append(
            "source_staging missing required attrs: "
            + ", ".join(staging["missing_required_attrs"])
        )
    if delta["missing_required_attrs"]:
        reasons.append(
            "source_child_delta missing required attrs: "
            + ", ".join(delta["missing_required_attrs"])
        )

    for section, summary in (
        ("provider_source", provider),
        ("source_staging", staging),
        ("source_child_delta", delta),
    ):
        _require_bool(reasons, summary, section, "diagnostic_only", True)
        _require_bool(reasons, summary, section, "gate_candidate", False)
        _require_bool(reasons, summary, section, "integrator_output", False)
        _require_bool(reasons, summary, section, "writes_candidate", False)
        if "writes_netcdf" in summary and summary["writes_netcdf"] is not None:
            _require_bool(reasons, summary, section, "writes_netcdf", False)

    for section, summary in (("provider_source", provider), ("source_staging", staging)):
        _require_bool(reasons, summary, section, "uses_reference_end_truth", False)
        _require_bool(reasons, summary, section, "uses_direct_p_shortcut", False)
        _require_bool(reasons, summary, section, "reads_direct_p", False)

    _require_bool(reasons, provider, "provider_source", "provider_ok", True)
    _require_bool(reasons, provider, "provider_source", "wrote_pb", True)
    _require_bool(reasons, provider, "provider_source", "wrote_t_init", True)
    _require_bool(reasons, provider, "provider_source", "wrote_mub", True)
    _require_bool(reasons, provider, "provider_source", "wrote_alb", True)
    _require_bool(reasons, provider, "provider_source", "wrote_phb", True)
    if provider["no_candidate_write"] is False:
        reasons.append("provider_source.no_candidate_write is False, expected True")
    if provider["no_reference_end_truth"] is False:
        reasons.append("provider_source.no_reference_end_truth is False, expected True")
    if provider["no_direct_p_shortcut"] is False:
        reasons.append("provider_source.no_direct_p_shortcut is False, expected True")
    _require_bool(
        reasons,
        provider,
        "provider_source",
        "phb_diagnostic_only_marker",
        True,
    )
    _require_bool(reasons, staging, "source_staging", "ok", True)
    if staging["aliases_child"] is True:
        reasons.append("source_staging.aliases_child is true, expected isolated staging")

    _validate_delta(reasons, delta)
    return reasons


def audit_diagnostic_adapter_provider_source(path: Path) -> AuditResult:
    with netCDF4.Dataset(path) as dataset:
        provider_attrs = _read_prefixed_attrs(dataset, PROVIDER_PREFIX)
        staging_attrs = _read_prefixed_attrs(dataset, SOURCE_STAGING_PREFIX)
        delta_attrs = _read_prefixed_attrs(dataset, SOURCE_CHILD_DELTA_PREFIX)

    provider = _provider_summary(provider_attrs)
    staging = _source_staging_summary(staging_attrs)
    delta = _source_child_delta_summary(delta_attrs)
    largest_delta = _largest_delta_field(delta)
    reasons = _validate_report(provider, staging, delta)
    status = "failed" if reasons else "passed"
    reasons.append(
        "D77 provider-source/source-staging/source-child-delta metadata is "
        "diagnostic-only and is not strict-gate evidence."
    )
    return AuditResult(
        status=status,
        path=str(path),
        gate_evidence=False,
        provider_source=provider,
        source_staging=staging,
        source_child_delta=delta,
        largest_delta_field=largest_delta,
        reasons=reasons,
    )


def report_to_json(report: AuditResult, *, pretty: bool = False) -> str:
    return json.dumps(
        _strict_json_value(report),
        indent=2 if pretty else None,
        allow_nan=False,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "path",
        type=Path,
        help="Selected-field hidden diagnostic adapter NetCDF output",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_diagnostic_adapter_provider_source(args.path)
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
