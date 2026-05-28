#!/usr/bin/env python
"""Diagnostic-only audit for selected-field U/V wind errors."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable

import netCDF4
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if __package__ in {None, ""} and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.audit_selected_field_state import (
    _coerce_bool,
    _coerce_int,
    _coerce_pair,
    _field_region_mask,
    _infer_target_region,
    _mass_horizontal_shape,
    _netcdf_attr_to_json,
    _pair_to_dict,
    _read_variable_data,
    _region_mask_for_data,
)


DEFAULT_WIND_VARIABLES = ("U", "V")
KEY_METADATA_ATTRS = (
    "TYWRF_DIAGNOSTIC_ONLY",
    "TYWRF_GATE_CANDIDATE",
    "TYWRF_INTEGRATOR_OUTPUT",
    "TYWRF_VALIDATION_GATE_ONLY",
    "TYWRF_CANDIDATE_KIND",
    "TYWRF_CANDIDATE_DOMAIN",
    "TYWRF_CYCLE_START",
    "TYWRF_CYCLE_END",
    "TYWRF_D02_RESOLUTION_CHECK",
    "DX",
    "DY",
    "TYWRF_FROM_PARENT_START",
    "TYWRF_TO_PARENT_START",
    "TYWRF_REMAP_FROM_PARENT_START",
    "TYWRF_REMAP_TO_PARENT_START",
    "TYWRF_PARENT_GRID_RATIO",
    "I_PARENT_START",
    "J_PARENT_START",
    "CEN_LAT",
    "CEN_LON",
    "TYWRF_SELECTED_FIELD_TIMELINE_VERSION",
    "TYWRF_SELECTED_FIELD_TIMELINE_EVIDENCE_ONLY",
    "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_COUNT",
    "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_NAMES",
    "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS",
    "TYWRF_SELECTED_FIELD_CHANGED_POINTS",
)


@dataclass(frozen=True)
class WindErrorMetrics:
    rmse: float
    normalized_rmse: float
    max_abs: float
    valid_count: int
    total_count: int
    reference_finite_count: int
    candidate_finite_count: int
    max_abs_location_index: tuple[int, ...] | None
    reference_value_at_max_abs: float | None
    candidate_value_at_max_abs: float | None
    sum_squared_error: float


@dataclass(frozen=True)
class WindVariableAudit:
    variable: str
    status: str
    diagnostic_only: bool
    gate_evidence: bool
    advances_00_20: bool
    reference_present: bool
    candidate_present: bool
    reference_shape: tuple[int, ...] | None = None
    candidate_shape: tuple[int, ...] | None = None
    whole_domain: WindErrorMetrics | None = None
    region_breakdown_status: str | None = None
    region_breakdown: dict[str, Any] | None = None
    message: str | None = None


@dataclass(frozen=True)
class SelectedFieldWindErrorAudit:
    reference: str
    candidate: str
    domain: str
    status: str
    diagnostic_only: bool
    gate_evidence: bool
    advances_00_20: bool
    summary: dict[str, Any]
    metadata: dict[str, Any]
    movement_region: dict[str, Any]
    variables: list[WindVariableAudit]
    diagnosis: dict[str, Any]


def _finite_array(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    masked = np.ma.masked_invalid(np.ma.asarray(data, dtype=np.float64))
    values = np.asarray(masked.filled(np.nan), dtype=np.float64)
    finite = np.isfinite(values) & ~np.ma.getmaskarray(masked)
    return values, finite


def wind_error_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    selection_mask: np.ndarray | None = None,
) -> WindErrorMetrics:
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
    reference_finite_count = int(np.count_nonzero(selected & ref_finite))
    candidate_finite_count = int(np.count_nonzero(selected & cand_finite))
    total_count = int(np.count_nonzero(selected))
    valid_count = int(np.count_nonzero(pair_mask))
    if valid_count == 0:
        return WindErrorMetrics(
            rmse=math.nan,
            normalized_rmse=math.nan,
            max_abs=math.nan,
            valid_count=0,
            total_count=total_count,
            reference_finite_count=reference_finite_count,
            candidate_finite_count=candidate_finite_count,
            max_abs_location_index=None,
            reference_value_at_max_abs=None,
            candidate_value_at_max_abs=None,
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

    abs_error = np.abs(cand_values - ref_values)
    masked_abs_error = np.where(pair_mask, abs_error, -math.inf)
    flat_index = int(np.argmax(masked_abs_error))
    location = tuple(int(value) for value in np.unravel_index(flat_index, ref_values.shape))
    return WindErrorMetrics(
        rmse=rmse,
        normalized_rmse=float(normalized),
        max_abs=float(masked_abs_error[location]),
        valid_count=valid_count,
        total_count=total_count,
        reference_finite_count=reference_finite_count,
        candidate_finite_count=candidate_finite_count,
        max_abs_location_index=location,
        reference_value_at_max_abs=float(ref_values[location]),
        candidate_value_at_max_abs=float(cand_values[location]),
        sum_squared_error=sum_squared_error,
    )


def _read_candidate_attrs(dataset: netCDF4.Dataset) -> dict[str, Any]:
    attrs: dict[str, Any] = {}
    for name in dataset.ncattrs():
        if (
            name in KEY_METADATA_ATTRS
            or name.startswith("TYWRF_SELECTED_FIELD_")
            or name.startswith("TYWRF_REMAP_")
        ):
            attrs[name] = _netcdf_attr_to_json(dataset.getncattr(name))
    return attrs


def _coerce_float(value: Any) -> float | None:
    if value is None or isinstance(value, (bool, np.bool_)):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _split_csv_attr(value: Any) -> list[str]:
    if value is None:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _parse_timeline_events(value: Any) -> list[dict[str, Any]]:
    if value is None:
        return []
    events: list[dict[str, Any]] = []
    for raw_event in str(value).split("|"):
        text = raw_event.strip()
        if not text:
            continue
        match = re.fullmatch(r"(\d+):([^()]+)\((.*)\)", text)
        if not match:
            events.append({"raw": text, "parsed": False})
            continue
        fields: dict[str, str] = {}
        for raw_field in match.group(3).split(","):
            if not raw_field:
                continue
            if "=" not in raw_field:
                fields[raw_field] = ""
                continue
            key, raw_value = raw_field.split("=", 1)
            fields[key] = raw_value
        events.append(
            {
                "index": int(match.group(1)),
                "name": match.group(2),
                "fields": fields,
                "parsed": True,
            }
        )
    return events


def _candidate_metadata(candidate: netCDF4.Dataset) -> dict[str, Any]:
    attrs = _read_candidate_attrs(candidate)
    timeline_events_raw = attrs.get("TYWRF_SELECTED_FIELD_TIMELINE_EVENTS")
    return {
        "status": "available",
        "diagnostic_only": True,
        "gate_evidence": False,
        "advances_00_20": False,
        "attrs": attrs,
        "present_attrs": sorted(attrs),
        "missing_key_attrs": [name for name in KEY_METADATA_ATTRS if name not in attrs],
        "parsed": {
            "diagnostic_only": _coerce_bool(attrs.get("TYWRF_DIAGNOSTIC_ONLY")),
            "gate_candidate": _coerce_bool(attrs.get("TYWRF_GATE_CANDIDATE")),
            "integrator_output": _coerce_bool(attrs.get("TYWRF_INTEGRATOR_OUTPUT")),
            "validation_gate_only": _coerce_bool(
                attrs.get("TYWRF_VALIDATION_GATE_ONLY")
            ),
            "candidate_kind": attrs.get("TYWRF_CANDIDATE_KIND"),
            "candidate_domain": attrs.get("TYWRF_CANDIDATE_DOMAIN"),
            "cycle_start": attrs.get("TYWRF_CYCLE_START"),
            "cycle_end": attrs.get("TYWRF_CYCLE_END"),
            "dx": _coerce_float(attrs.get("DX")),
            "dy": _coerce_float(attrs.get("DY")),
            "d02_resolution_check": attrs.get("TYWRF_D02_RESOLUTION_CHECK"),
            "from_parent_start": _pair_to_dict(
                _coerce_pair(attrs.get("TYWRF_FROM_PARENT_START"))
            ),
            "to_parent_start": _pair_to_dict(
                _coerce_pair(attrs.get("TYWRF_TO_PARENT_START"))
            ),
            "remap_from_parent_start": _pair_to_dict(
                _coerce_pair(attrs.get("TYWRF_REMAP_FROM_PARENT_START"))
            ),
            "remap_to_parent_start": _pair_to_dict(
                _coerce_pair(attrs.get("TYWRF_REMAP_TO_PARENT_START"))
            ),
            "parent_grid_ratio": _coerce_int(attrs.get("TYWRF_PARENT_GRID_RATIO")),
            "i_parent_start": _coerce_int(attrs.get("I_PARENT_START")),
            "j_parent_start": _coerce_int(attrs.get("J_PARENT_START")),
            "timeline_version": attrs.get("TYWRF_SELECTED_FIELD_TIMELINE_VERSION"),
            "timeline_evidence_only": _coerce_bool(
                attrs.get("TYWRF_SELECTED_FIELD_TIMELINE_EVIDENCE_ONLY")
            ),
            "timeline_event_count": _coerce_int(
                attrs.get("TYWRF_SELECTED_FIELD_TIMELINE_EVENT_COUNT")
            ),
            "timeline_event_names": _split_csv_attr(
                attrs.get("TYWRF_SELECTED_FIELD_TIMELINE_EVENT_NAMES")
            ),
            "timeline_events": _parse_timeline_events(timeline_events_raw),
            "timeline_events_raw": timeline_events_raw,
        },
    }


def _metrics_to_dict(metrics: WindErrorMetrics) -> dict[str, Any]:
    return _strict_json_value(asdict(metrics))


def _region_breakdown(
    variable: str,
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    whole_domain: WindErrorMetrics,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
) -> dict[str, Any]:
    if reference.ndim < 2:
        return {
            "status": "insufficient_metadata",
            "region_breakdown_status": "insufficient_metadata",
            "whole_domain": _metrics_to_dict(whole_domain),
            "diagnostic_only": True,
            "gate_evidence": False,
            "advances_00_20": False,
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
            "status": "insufficient_metadata",
            "region_breakdown_status": "insufficient_metadata",
            "whole_domain": _metrics_to_dict(whole_domain),
            "diagnostic_only": True,
            "gate_evidence": False,
            "advances_00_20": False,
            "message": mask_info.get(
                "message",
                "insufficient metadata for reliable overlap/exposed partition",
            ),
        }

    exposed_mask = _region_mask_for_data(reference, mask_2d)
    overlap_mask = ~exposed_mask
    exposed = wind_error_metrics(reference, candidate, selection_mask=exposed_mask)
    overlap = wind_error_metrics(reference, candidate, selection_mask=overlap_mask)
    if (
        math.isfinite(exposed.sum_squared_error)
        and math.isfinite(whole_domain.sum_squared_error)
        and whole_domain.sum_squared_error > 0.0
    ):
        exposed_error_fraction = exposed.sum_squared_error / whole_domain.sum_squared_error
    else:
        exposed_error_fraction = None

    return {
        **mask_info,
        "status": "available",
        "region_breakdown_status": "available",
        "diagnostic_only": True,
        "gate_evidence": False,
        "advances_00_20": False,
        "whole_domain": _metrics_to_dict(whole_domain),
        "exposed_region": _metrics_to_dict(exposed),
        "overlap_region": _metrics_to_dict(overlap),
        "exposed_error_fraction": exposed_error_fraction,
        "exposed_region_dominates_global_error": (
            None if exposed_error_fraction is None else exposed_error_fraction > 0.5
        ),
    }


def compare_wind_variable(
    reference: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    variable: str,
    *,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    time_index: int,
) -> WindVariableAudit:
    reference_present = variable in reference.variables
    candidate_present = variable in candidate.variables
    if not reference_present or not candidate_present:
        missing_owner = []
        if not reference_present:
            missing_owner.append("reference")
        if not candidate_present:
            missing_owner.append("candidate")
        return WindVariableAudit(
            variable=variable,
            status="not_available",
            diagnostic_only=True,
            gate_evidence=False,
            advances_00_20=False,
            reference_present=reference_present,
            candidate_present=candidate_present,
            region_breakdown_status="insufficient_metadata",
            region_breakdown={
                "status": "insufficient_metadata",
                "region_breakdown_status": "insufficient_metadata",
                "message": "variable missing",
            },
            message=f"variable missing from {' and '.join(missing_owner)}",
        )

    try:
        ref_data = _read_variable_data(reference, variable, time_index=time_index)
        cand_data = _read_variable_data(candidate, variable, time_index=time_index)
    except (IndexError, TypeError, ValueError) as exc:
        return WindVariableAudit(
            variable=variable,
            status="not_available",
            diagnostic_only=True,
            gate_evidence=False,
            advances_00_20=False,
            reference_present=True,
            candidate_present=True,
            message=f"variable cannot be read: {exc}",
        )

    ref_shape = tuple(int(value) for value in ref_data.shape)
    cand_shape = tuple(int(value) for value in cand_data.shape)
    if ref_shape != cand_shape:
        return WindVariableAudit(
            variable=variable,
            status="shape_mismatch",
            diagnostic_only=True,
            gate_evidence=False,
            advances_00_20=False,
            reference_present=True,
            candidate_present=True,
            reference_shape=ref_shape,
            candidate_shape=cand_shape,
            message=f"reference shape {ref_shape} != candidate shape {cand_shape}",
        )

    try:
        ref_values = np.asarray(ref_data, dtype=np.float64)
        cand_values = np.asarray(cand_data, dtype=np.float64)
        whole_domain = wind_error_metrics(ref_values, cand_values)
    except (TypeError, ValueError) as exc:
        return WindVariableAudit(
            variable=variable,
            status="non_numeric",
            diagnostic_only=True,
            gate_evidence=False,
            advances_00_20=False,
            reference_present=True,
            candidate_present=True,
            reference_shape=ref_shape,
            candidate_shape=cand_shape,
            message=f"variable cannot be compared numerically: {exc}",
        )

    status = "computed" if whole_domain.valid_count > 0 else "no_valid_samples"
    breakdown = _region_breakdown(
        variable,
        ref_values,
        cand_values,
        whole_domain=whole_domain,
        region=region,
        mass_shape=mass_shape,
    )
    return WindVariableAudit(
        variable=variable,
        status=status,
        diagnostic_only=True,
        gate_evidence=False,
        advances_00_20=False,
        reference_present=True,
        candidate_present=True,
        reference_shape=ref_shape,
        candidate_shape=cand_shape,
        whole_domain=whole_domain,
        region_breakdown_status=breakdown["region_breakdown_status"],
        region_breakdown=breakdown,
        message=None if status == "computed" else "no finite sample pairs are available",
    )


def audit_selected_field_wind_error(
    reference_path: Path,
    candidate_path: Path,
    *,
    domain: str = "d02",
    variables: Iterable[str] = DEFAULT_WIND_VARIABLES,
    time_index: int = -1,
) -> SelectedFieldWindErrorAudit:
    variables = tuple(variables)
    with netCDF4.Dataset(reference_path) as reference, netCDF4.Dataset(
        candidate_path
    ) as candidate:
        metadata = _candidate_metadata(candidate)
        region = {
            **_infer_target_region(candidate),
            "diagnostic_only": True,
            "gate_evidence": False,
            "advances_00_20": False,
        }
        mass_shape = _mass_horizontal_shape(
            reference,
            candidate,
            time_index=time_index,
        )
        audits = [
            compare_wind_variable(
                reference,
                candidate,
                variable,
                region=region,
                mass_shape=mass_shape,
                time_index=time_index,
            )
            for variable in variables
        ]

    counts = Counter(audit.status for audit in audits)
    region_status = {
        audit.variable: audit.region_breakdown_status for audit in audits
    }
    summary = {
        **dict(sorted(counts.items())),
        "total": len(audits),
        "available": sum(1 for audit in audits if audit.status != "not_available"),
        "not_available": counts.get("not_available", 0),
        "region_breakdown_status": region_status,
        "diagnostic_only": True,
        "gate_evidence": False,
        "advances_00_20": False,
    }
    status = "computed_with_flags" if any(
        audit.status != "computed" or audit.region_breakdown_status != "available"
        for audit in audits
    ) else "computed"
    return SelectedFieldWindErrorAudit(
        reference=str(reference_path),
        candidate=str(candidate_path),
        domain=domain,
        status=status,
        diagnostic_only=True,
        gate_evidence=False,
        advances_00_20=False,
        summary=summary,
        metadata=metadata,
        movement_region=region,
        variables=audits,
        diagnosis={
            "diagnostic_only": True,
            "gate_evidence": False,
            "advances_00_20": False,
            "message": (
                "U/V selected-field wind error audit only; results are not gate "
                "evidence and must not be used to generate or alter candidates"
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


def report_to_dict(report: SelectedFieldWindErrorAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: SelectedFieldWindErrorAudit,
    *,
    pretty: bool = False,
) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="WRF reference wrfout file")
    parser.add_argument("candidate", type=Path, help="TyWRF candidate wrfout file")
    parser.add_argument(
        "--domain",
        default="d02",
        help="Domain label to stamp in the report; default: d02",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_WIND_VARIABLES),
        help="Wind variables to audit; default: U V",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for Time-leading variables",
    )
    parser.add_argument("--output", type=Path, help="Optional JSON output path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    report = audit_selected_field_wind_error(
        args.reference,
        args.candidate,
        domain=args.domain,
        variables=args.variables,
        time_index=args.time_index,
    )
    payload = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
