#!/usr/bin/env python
"""Diagnostic-only attribution for selected-field U/V wind source errors."""

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
    _field_region_mask,
    _infer_target_region,
    _mass_horizontal_shape,
    _read_variable_data,
    _region_mask_for_data,
)
from tools.audit_selected_field_wind_error import (
    DEFAULT_WIND_VARIABLES,
    WindErrorMetrics,
    _candidate_metadata,
    wind_error_metrics,
)


@dataclass(frozen=True)
class WindSourceVariableAttribution:
    variable: str
    status: str
    diagnostic_only: bool
    gate_evidence: bool
    advances_00_20: bool
    d02_start_present: bool
    candidate_end_present: bool
    reference_end_present: bool
    d02_start_shape: tuple[int, ...] | None = None
    candidate_end_shape: tuple[int, ...] | None = None
    reference_end_shape: tuple[int, ...] | None = None
    whole_domain: dict[str, Any] | None = None
    region_breakdown_status: str | None = None
    region_breakdown: dict[str, Any] | None = None
    message: str | None = None


@dataclass(frozen=True)
class WindSourceAttributionAudit:
    d02_start: str
    candidate_end: str
    reference_end: str
    d01_parent_start: str | None
    d01_parent_end: str | None
    domain: str
    status: str
    diagnostic_only: bool
    gate_evidence: bool
    advances_00_20: bool
    summary: dict[str, Any]
    metadata: dict[str, Any]
    movement_region: dict[str, Any]
    source_pose: dict[str, Any]
    parent_context: dict[str, Any]
    variables: list[WindSourceVariableAttribution]
    diagnosis: dict[str, Any]


def _metrics_to_dict(metrics: WindErrorMetrics) -> dict[str, Any]:
    return _strict_json_value(asdict(metrics))


def _finite_values(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    masked = np.ma.masked_invalid(np.ma.asarray(data, dtype=np.float64))
    values = np.asarray(masked.filled(np.nan), dtype=np.float64)
    finite = np.isfinite(values) & ~np.ma.getmaskarray(masked)
    return values, finite


def _safe_ratio(numerator: float, denominator: float) -> float | None:
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return None
    if denominator <= 0.0:
        return None
    return float(numerator / denominator)


def _json_pair_to_tuple(value: Any) -> tuple[int, int] | None:
    if not isinstance(value, dict):
        return None
    i_value = value.get("i")
    j_value = value.get("j")
    if not isinstance(i_value, int) or not isinstance(j_value, int):
        return None
    return i_value, j_value


def _pose_check(
    name: str,
    from_value: Any,
    to_value: Any,
) -> dict[str, Any]:
    from_pair = _json_pair_to_tuple(from_value)
    to_pair = _json_pair_to_tuple(to_value)
    if from_pair is None or to_pair is None:
        missing = []
        if from_pair is None:
            missing.append(f"{name}_from_parent_start")
        if to_pair is None:
            missing.append(f"{name}_to_parent_start")
        return {
            "name": name,
            "status": "insufficient_metadata",
            "missing_metadata": missing,
            "from_parent_start": from_value,
            "to_parent_start": to_value,
            "moved": None,
            "child_delta": None,
        }
    delta = {
        "i": int(to_pair[0] - from_pair[0]),
        "j": int(to_pair[1] - from_pair[1]),
    }
    return {
        "name": name,
        "status": "available",
        "from_parent_start": {"i": from_pair[0], "j": from_pair[1]},
        "to_parent_start": {"i": to_pair[0], "j": to_pair[1]},
        "moved": delta != {"i": 0, "j": 0},
        "child_delta": delta,
    }


def _source_pose_from_metadata(metadata: dict[str, Any]) -> dict[str, Any]:
    parsed = metadata.get("parsed", {})
    checks = [
        _pose_check(
            "parent",
            parsed.get("from_parent_start"),
            parsed.get("to_parent_start"),
        ),
        _pose_check(
            "remap",
            parsed.get("remap_from_parent_start"),
            parsed.get("remap_to_parent_start"),
        ),
    ]
    available = [entry for entry in checks if entry["status"] == "available"]
    moved = [entry for entry in available if entry["moved"] is True]
    if moved:
        status = "raw_start_not_candidate_pose"
        aligned = False
        interpretation = "raw_pose_only"
        evidence = False
        message = (
            "candidate metadata indicates moving-nest pose changed; "
            "candidate_vs_d02_start and explanatory ratios compare the raw d02 "
            "start pose to candidate end pose and are not shifted/remapped start "
            "persistence evidence"
        )
    elif available:
        status = "same_pose"
        aligned = True
        interpretation = "same_pose_start_persistence_diagnostic"
        evidence = True
        message = (
            "candidate metadata indicates from/to parent starts are unchanged; "
            "raw d02 start comparison is same-pose diagnostic evidence only"
        )
    else:
        status = "insufficient_metadata"
        aligned = False
        interpretation = "metadata_insufficient_raw_pose_diagnostic"
        evidence = False
        message = (
            "candidate pose metadata is insufficient; raw d02 start comparison "
            "must not be interpreted as pose-aligned persistence evidence"
        )

    return {
        "status": status,
        "source_pose_status": status,
        "pose_aligned_with_candidate": aligned,
        "pose_aligned_start_input_provided": False,
        "candidate_vs_d02_start_interpretation": interpretation,
        "pose_aligned_start_persistence_evidence": evidence,
        "not_shifted_start_persistence_evidence": not evidence,
        "diagnostic_only": True,
        "gate_evidence": False,
        "advances_00_20": False,
        "checks": checks,
        "message": message,
    }


def _pose_annotation(source_pose: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_pose_status": source_pose["source_pose_status"],
        "pose_aligned_with_candidate": source_pose["pose_aligned_with_candidate"],
        "candidate_vs_d02_start_interpretation": source_pose[
            "candidate_vs_d02_start_interpretation"
        ],
        "pose_aligned_start_persistence_evidence": source_pose[
            "pose_aligned_start_persistence_evidence"
        ],
        "not_shifted_start_persistence_evidence": source_pose[
            "not_shifted_start_persistence_evidence"
        ],
        "source_pose_message": source_pose["message"],
    }


def attribution_ratios(
    reference_end: np.ndarray,
    d02_start: np.ndarray,
    candidate_end: np.ndarray,
    *,
    source_pose: dict[str, Any],
    selection_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    pose_fields = _pose_annotation(source_pose)
    ref_values, ref_finite = _finite_values(reference_end)
    start_values, start_finite = _finite_values(d02_start)
    candidate_values, candidate_finite = _finite_values(candidate_end)
    if ref_values.shape != start_values.shape or ref_values.shape != candidate_values.shape:
        raise ValueError(
            "shape mismatch: "
            f"reference={ref_values.shape}, start={start_values.shape}, "
            f"candidate={candidate_values.shape}"
        )

    if selection_mask is None:
        selected = np.ones(ref_values.shape, dtype=bool)
    else:
        selected = np.asarray(selection_mask, dtype=bool)
        if selected.shape != ref_values.shape:
            selected = np.broadcast_to(selected, ref_values.shape)

    valid = selected & ref_finite & start_finite & candidate_finite
    valid_count = int(np.count_nonzero(valid))
    selected_count = int(np.count_nonzero(selected))
    if valid_count == 0:
        return {
            "status": "no_valid_samples",
            **pose_fields,
            "valid_count": 0,
            "selected_count": selected_count,
            "start_state_persistence_fraction": None,
            "candidate_error_fraction_of_start_persistence": None,
            "candidate_change_fraction_of_start_persistence": None,
            "candidate_delta_fraction_of_candidate_error": None,
            "candidate_error_projection_on_start_persistence_fraction": None,
            "candidate_error_cosine_with_start_persistence": None,
            "message": "no finite reference/start/candidate triplets are available",
        }

    start_error = start_values[valid] - ref_values[valid]
    candidate_error = candidate_values[valid] - ref_values[valid]
    candidate_delta = candidate_values[valid] - start_values[valid]
    start_sse = float(np.sum(start_error * start_error))
    candidate_sse = float(np.sum(candidate_error * candidate_error))
    delta_sse = float(np.sum(candidate_delta * candidate_delta))
    dot = float(np.sum(candidate_error * start_error))
    norm_product = math.sqrt(candidate_sse * start_sse)

    return {
        "status": (
            "available"
            if source_pose["source_pose_status"] == "same_pose"
            else "raw_pose_only"
        ),
        **pose_fields,
        "valid_count": valid_count,
        "selected_count": selected_count,
        "start_persistence_sse": start_sse,
        "candidate_error_sse": candidate_sse,
        "candidate_delta_sse": delta_sse,
        "start_state_persistence_fraction": _safe_ratio(start_sse, candidate_sse),
        "candidate_error_fraction_of_start_persistence": _safe_ratio(
            candidate_sse, start_sse
        ),
        "candidate_change_fraction_of_start_persistence": _safe_ratio(
            delta_sse, start_sse
        ),
        "candidate_delta_fraction_of_candidate_error": _safe_ratio(
            delta_sse, candidate_sse
        ),
        "candidate_error_projection_on_start_persistence_fraction": _safe_ratio(
            dot, candidate_sse
        ),
        "candidate_error_cosine_with_start_persistence": _safe_ratio(dot, norm_product),
        "message": (
            "Ratios use finite reference/start/candidate triplets only. "
            + source_pose["message"]
        ),
    }


def _comparison_bundle(
    reference_end: np.ndarray,
    d02_start: np.ndarray,
    candidate_end: np.ndarray,
    *,
    source_pose: dict[str, Any],
    selection_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    candidate_vs_reference = wind_error_metrics(
        reference_end,
        candidate_end,
        selection_mask=selection_mask,
    )
    d02_start_vs_reference = wind_error_metrics(
        reference_end,
        d02_start,
        selection_mask=selection_mask,
    )
    candidate_vs_d02_start = wind_error_metrics(
        d02_start,
        candidate_end,
        selection_mask=selection_mask,
    )
    candidate_vs_d02_start_payload = {
        **_metrics_to_dict(candidate_vs_d02_start),
        **_pose_annotation(source_pose),
    }
    return {
        **_pose_annotation(source_pose),
        "candidate_vs_reference": _metrics_to_dict(candidate_vs_reference),
        "d02_start_vs_reference": _metrics_to_dict(d02_start_vs_reference),
        "candidate_vs_d02_start": candidate_vs_d02_start_payload,
        "explanatory_ratios": attribution_ratios(
            reference_end,
            d02_start,
            candidate_end,
            source_pose=source_pose,
            selection_mask=selection_mask,
        ),
    }


def _region_attribution(
    variable: str,
    reference_end: np.ndarray,
    d02_start: np.ndarray,
    candidate_end: np.ndarray,
    *,
    whole_domain: dict[str, Any],
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    source_pose: dict[str, Any],
) -> dict[str, Any]:
    if reference_end.ndim < 2:
        return {
            "status": "insufficient_metadata",
            "region_breakdown_status": "insufficient_metadata",
            "whole_domain": whole_domain,
            "diagnostic_only": True,
            "gate_evidence": False,
            "advances_00_20": False,
            "message": "variable has no horizontal dimensions",
        }

    mask_2d, mask_info = _field_region_mask(
        variable,
        (int(reference_end.shape[-2]), int(reference_end.shape[-1])),
        mass_shape,
        region,
    )
    if mask_2d is None:
        return {
            **mask_info,
            "status": "insufficient_metadata",
            "region_breakdown_status": "insufficient_metadata",
            "whole_domain": whole_domain,
            "diagnostic_only": True,
            "gate_evidence": False,
            "advances_00_20": False,
            "message": mask_info.get(
                "message",
                "insufficient metadata for reliable overlap/exposed partition",
            ),
        }

    exposed_mask = _region_mask_for_data(reference_end, mask_2d)
    overlap_mask = ~exposed_mask
    exposed = _comparison_bundle(
        reference_end,
        d02_start,
        candidate_end,
        source_pose=source_pose,
        selection_mask=exposed_mask,
    )
    overlap = _comparison_bundle(
        reference_end,
        d02_start,
        candidate_end,
        source_pose=source_pose,
        selection_mask=overlap_mask,
    )
    whole_candidate_sse = whole_domain["candidate_vs_reference"]["sum_squared_error"]
    exposed_candidate_sse = exposed["candidate_vs_reference"]["sum_squared_error"]
    overlap_candidate_sse = overlap["candidate_vs_reference"]["sum_squared_error"]
    exposed_error_fraction = (
        None
        if whole_candidate_sse in {None, 0.0}
        else _safe_ratio(float(exposed_candidate_sse), float(whole_candidate_sse))
    )
    overlap_error_fraction = (
        None
        if whole_candidate_sse in {None, 0.0}
        else _safe_ratio(float(overlap_candidate_sse), float(whole_candidate_sse))
    )

    return {
        **mask_info,
        "status": "available",
        "region_breakdown_status": "available",
        "diagnostic_only": True,
        "gate_evidence": False,
        "advances_00_20": False,
        "whole_domain": whole_domain,
        "exposed_region": exposed,
        "overlap_region": overlap,
        "exposed_candidate_error_fraction": exposed_error_fraction,
        "overlap_candidate_error_fraction": overlap_error_fraction,
        "exposed_region_dominates_candidate_error": (
            None if exposed_error_fraction is None else exposed_error_fraction > 0.5
        ),
        "overlap_region_dominates_candidate_error": (
            None if overlap_error_fraction is None else overlap_error_fraction > 0.5
        ),
    }


def _missing_variable_attribution(
    variable: str,
    *,
    d02_start_present: bool,
    candidate_end_present: bool,
    reference_end_present: bool,
) -> WindSourceVariableAttribution:
    missing = []
    if not d02_start_present:
        missing.append("d02-start")
    if not candidate_end_present:
        missing.append("candidate-end")
    if not reference_end_present:
        missing.append("reference-end")
    return WindSourceVariableAttribution(
        variable=variable,
        status="not_available",
        diagnostic_only=True,
        gate_evidence=False,
        advances_00_20=False,
        d02_start_present=d02_start_present,
        candidate_end_present=candidate_end_present,
        reference_end_present=reference_end_present,
        region_breakdown_status="insufficient_metadata",
        region_breakdown={
            "status": "insufficient_metadata",
            "region_breakdown_status": "insufficient_metadata",
            "message": "variable missing",
        },
        message=f"variable missing from {' and '.join(missing)}",
    )


def compare_wind_source_variable(
    d02_start: netCDF4.Dataset,
    candidate_end: netCDF4.Dataset,
    reference_end: netCDF4.Dataset,
    variable: str,
    *,
    region: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    source_pose: dict[str, Any],
    time_index: int,
    start_time_index: int,
) -> WindSourceVariableAttribution:
    d02_start_present = variable in d02_start.variables
    candidate_end_present = variable in candidate_end.variables
    reference_end_present = variable in reference_end.variables
    if not d02_start_present or not candidate_end_present or not reference_end_present:
        return _missing_variable_attribution(
            variable,
            d02_start_present=d02_start_present,
            candidate_end_present=candidate_end_present,
            reference_end_present=reference_end_present,
        )

    try:
        start_data = _read_variable_data(
            d02_start,
            variable,
            time_index=start_time_index,
        )
        candidate_data = _read_variable_data(
            candidate_end,
            variable,
            time_index=time_index,
        )
        reference_data = _read_variable_data(
            reference_end,
            variable,
            time_index=time_index,
        )
    except (IndexError, TypeError, ValueError) as exc:
        return WindSourceVariableAttribution(
            variable=variable,
            status="not_available",
            diagnostic_only=True,
            gate_evidence=False,
            advances_00_20=False,
            d02_start_present=True,
            candidate_end_present=True,
            reference_end_present=True,
            message=f"variable cannot be read: {exc}",
        )

    start_shape = tuple(int(value) for value in start_data.shape)
    candidate_shape = tuple(int(value) for value in candidate_data.shape)
    reference_shape = tuple(int(value) for value in reference_data.shape)
    if start_shape != reference_shape or candidate_shape != reference_shape:
        return WindSourceVariableAttribution(
            variable=variable,
            status="shape_mismatch",
            diagnostic_only=True,
            gate_evidence=False,
            advances_00_20=False,
            d02_start_present=True,
            candidate_end_present=True,
            reference_end_present=True,
            d02_start_shape=start_shape,
            candidate_end_shape=candidate_shape,
            reference_end_shape=reference_shape,
            message=(
                f"reference-end shape {reference_shape}, d02-start shape "
                f"{start_shape}, candidate-end shape {candidate_shape}"
            ),
        )

    try:
        reference_values = np.asarray(reference_data, dtype=np.float64)
        start_values = np.asarray(start_data, dtype=np.float64)
        candidate_values = np.asarray(candidate_data, dtype=np.float64)
        whole_domain = _comparison_bundle(
            reference_values,
            start_values,
            candidate_values,
            source_pose=source_pose,
        )
    except (TypeError, ValueError) as exc:
        return WindSourceVariableAttribution(
            variable=variable,
            status="non_numeric",
            diagnostic_only=True,
            gate_evidence=False,
            advances_00_20=False,
            d02_start_present=True,
            candidate_end_present=True,
            reference_end_present=True,
            d02_start_shape=start_shape,
            candidate_end_shape=candidate_shape,
            reference_end_shape=reference_shape,
            message=f"variable cannot be compared numerically: {exc}",
        )

    status = (
        "computed"
        if whole_domain["candidate_vs_reference"]["valid_count"] > 0
        else "no_valid_samples"
    )
    breakdown = _region_attribution(
        variable,
        reference_values,
        start_values,
        candidate_values,
        whole_domain=whole_domain,
        region=region,
        mass_shape=mass_shape,
        source_pose=source_pose,
    )
    return WindSourceVariableAttribution(
        variable=variable,
        status=status,
        diagnostic_only=True,
        gate_evidence=False,
        advances_00_20=False,
        d02_start_present=True,
        candidate_end_present=True,
        reference_end_present=True,
        d02_start_shape=start_shape,
        candidate_end_shape=candidate_shape,
        reference_end_shape=reference_shape,
        whole_domain=whole_domain,
        region_breakdown_status=breakdown["region_breakdown_status"],
        region_breakdown=breakdown,
        message=None if status == "computed" else "no finite sample pairs are available",
    )


def _parent_variable_summary(
    dataset: netCDF4.Dataset,
    variable: str,
    *,
    time_index: int,
) -> dict[str, Any]:
    if variable not in dataset.variables:
        return {"present": False}
    try:
        data = _read_variable_data(dataset, variable, time_index=time_index)
    except (IndexError, TypeError, ValueError) as exc:
        return {"present": True, "status": "not_available", "message": str(exc)}
    return {
        "present": True,
        "status": "available",
        "shape": tuple(int(value) for value in data.shape),
    }


def _parent_context(
    d01_parent_start_path: Path | None,
    d01_parent_end_path: Path | None,
    *,
    variables: Iterable[str],
    time_index: int,
    start_time_index: int,
) -> dict[str, Any]:
    context: dict[str, Any] = {
        "status": "not_requested",
        "diagnostic_only": True,
        "gate_evidence": False,
        "advances_00_20": False,
        "d01_parent_start": None if d01_parent_start_path is None else str(d01_parent_start_path),
        "d01_parent_end": None if d01_parent_end_path is None else str(d01_parent_end_path),
        "variables": {},
    }
    if d01_parent_start_path is None and d01_parent_end_path is None:
        return context
    if d01_parent_start_path is None or d01_parent_end_path is None:
        context["status"] = "partial"
        context["message"] = "both parent start and parent end are required for parent evolution metrics"
        return context

    with netCDF4.Dataset(d01_parent_start_path) as parent_start, netCDF4.Dataset(
        d01_parent_end_path
    ) as parent_end:
        context["status"] = "available"
        for variable in variables:
            start_summary = _parent_variable_summary(
                parent_start,
                variable,
                time_index=start_time_index,
            )
            end_summary = _parent_variable_summary(
                parent_end,
                variable,
                time_index=time_index,
            )
            entry = {
                "parent_start": start_summary,
                "parent_end": end_summary,
                "parent_end_vs_parent_start": None,
            }
            if start_summary.get("status") == "available" and end_summary.get(
                "status"
            ) == "available":
                try:
                    start_data = _read_variable_data(
                        parent_start,
                        variable,
                        time_index=start_time_index,
                    )
                    end_data = _read_variable_data(
                        parent_end,
                        variable,
                        time_index=time_index,
                    )
                    if start_data.shape == end_data.shape:
                        entry["parent_end_vs_parent_start"] = _metrics_to_dict(
                            wind_error_metrics(start_data, end_data)
                        )
                    else:
                        entry["parent_end_vs_parent_start"] = {
                            "status": "shape_mismatch",
                            "parent_start_shape": tuple(
                                int(value) for value in start_data.shape
                            ),
                            "parent_end_shape": tuple(
                                int(value) for value in end_data.shape
                            ),
                        }
                except (TypeError, ValueError, IndexError) as exc:
                    entry["parent_end_vs_parent_start"] = {
                        "status": "not_available",
                        "message": str(exc),
                    }
            context["variables"][variable] = entry
    return context


def _diagnosis(variables: list[WindSourceVariableAttribution]) -> dict[str, Any]:
    overlap_persistence: dict[str, float | None] = {}
    overlap_delta_fraction: dict[str, float | None] = {}
    exposed_candidate_fraction: dict[str, float | None] = {}
    region_status: dict[str, str | None] = {}

    for audit in variables:
        region = audit.region_breakdown or {}
        region_status[audit.variable] = audit.region_breakdown_status
        if region.get("region_breakdown_status") == "available":
            overlap_ratios = region["overlap_region"]["explanatory_ratios"]
            overlap_persistence[audit.variable] = overlap_ratios.get(
                "start_state_persistence_fraction"
            )
            overlap_delta_fraction[audit.variable] = overlap_ratios.get(
                "candidate_delta_fraction_of_candidate_error"
            )
            exposed_candidate_fraction[audit.variable] = region.get(
                "exposed_candidate_error_fraction"
            )
        else:
            overlap_persistence[audit.variable] = None
            overlap_delta_fraction[audit.variable] = None
            exposed_candidate_fraction[audit.variable] = None

    return {
        "diagnostic_only": True,
        "gate_evidence": False,
        "advances_00_20": False,
        "region_breakdown_status": region_status,
        "overlap_start_state_persistence_fraction": overlap_persistence,
        "overlap_candidate_delta_fraction_of_candidate_error": overlap_delta_fraction,
        "exposed_candidate_error_fraction": exposed_candidate_fraction,
        "message": (
            "U/V source/time-level attribution only; results are not gate evidence, "
            "must not alter candidate generation, and must not advance 00:20."
        ),
    }


def audit_selected_field_wind_source_attribution(
    d02_start_path: Path,
    candidate_end_path: Path,
    reference_end_path: Path,
    *,
    d01_parent_start_path: Path | None = None,
    d01_parent_end_path: Path | None = None,
    domain: str = "d02",
    variables: Iterable[str] = DEFAULT_WIND_VARIABLES,
    time_index: int = -1,
    start_time_index: int = -1,
) -> WindSourceAttributionAudit:
    variables = tuple(variables)
    with netCDF4.Dataset(d02_start_path) as d02_start, netCDF4.Dataset(
        candidate_end_path
    ) as candidate_end, netCDF4.Dataset(reference_end_path) as reference_end:
        metadata = _candidate_metadata(candidate_end)
        source_pose = _source_pose_from_metadata(metadata)
        region = {
            **_infer_target_region(candidate_end),
            "diagnostic_only": True,
            "gate_evidence": False,
            "advances_00_20": False,
        }
        mass_shape = _mass_horizontal_shape(
            reference_end,
            candidate_end,
            time_index=time_index,
        )
        audits = [
            compare_wind_source_variable(
                d02_start,
                candidate_end,
                reference_end,
                variable,
                region=region,
                mass_shape=mass_shape,
                source_pose=source_pose,
                time_index=time_index,
                start_time_index=start_time_index,
            )
            for variable in variables
        ]

    parent_context = _parent_context(
        d01_parent_start_path,
        d01_parent_end_path,
        variables=variables,
        time_index=time_index,
        start_time_index=start_time_index,
    )
    counts = Counter(audit.status for audit in audits)
    diagnosis = _diagnosis(audits)
    summary = {
        **dict(sorted(counts.items())),
        "total": len(audits),
        "available": sum(1 for audit in audits if audit.status != "not_available"),
        "not_available": counts.get("not_available", 0),
        "region_breakdown_status": diagnosis["region_breakdown_status"],
        "source_pose_status": source_pose["source_pose_status"],
        "pose_aligned_with_candidate": source_pose["pose_aligned_with_candidate"],
        "candidate_vs_d02_start_interpretation": source_pose[
            "candidate_vs_d02_start_interpretation"
        ],
        "source_time_levels": {
            "d02_start": {
                "path": str(d02_start_path),
                "time_index": start_time_index,
            },
            "candidate_end": {
                "path": str(candidate_end_path),
                "time_index": time_index,
            },
            "reference_end": {
                "path": str(reference_end_path),
                "time_index": time_index,
            },
            "d01_parent_start": (
                None
                if d01_parent_start_path is None
                else {
                    "path": str(d01_parent_start_path),
                    "time_index": start_time_index,
                }
            ),
            "d01_parent_end": (
                None
                if d01_parent_end_path is None
                else {
                    "path": str(d01_parent_end_path),
                    "time_index": time_index,
                }
            ),
        },
        "candidate_metadata_summary": {
            "dx": metadata["parsed"].get("dx"),
            "dy": metadata["parsed"].get("dy"),
            "d02_resolution_check": metadata["parsed"].get("d02_resolution_check"),
            "from_parent_start": metadata["parsed"].get("from_parent_start"),
            "to_parent_start": metadata["parsed"].get("to_parent_start"),
            "remap_from_parent_start": metadata["parsed"].get(
                "remap_from_parent_start"
            ),
            "remap_to_parent_start": metadata["parsed"].get("remap_to_parent_start"),
            "cycle_start": metadata["parsed"].get("cycle_start"),
            "cycle_end": metadata["parsed"].get("cycle_end"),
            "timeline_event_count": metadata["parsed"].get("timeline_event_count"),
            "timeline_event_names": metadata["parsed"].get("timeline_event_names"),
        },
        "parent_context_status": parent_context["status"],
        "diagnostic_only": True,
        "gate_evidence": False,
        "advances_00_20": False,
    }
    status = "computed_with_flags" if any(
        audit.status != "computed" or audit.region_breakdown_status != "available"
        for audit in audits
    ) else "computed"
    return WindSourceAttributionAudit(
        d02_start=str(d02_start_path),
        candidate_end=str(candidate_end_path),
        reference_end=str(reference_end_path),
        d01_parent_start=(
            None if d01_parent_start_path is None else str(d01_parent_start_path)
        ),
        d01_parent_end=None if d01_parent_end_path is None else str(d01_parent_end_path),
        domain=domain,
        status=status,
        diagnostic_only=True,
        gate_evidence=False,
        advances_00_20=False,
        summary=summary,
        metadata=metadata,
        movement_region=region,
        source_pose=source_pose,
        parent_context=parent_context,
        variables=audits,
        diagnosis=diagnosis,
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


def report_to_dict(report: WindSourceAttributionAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: WindSourceAttributionAudit,
    *,
    pretty: bool = False,
) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("d02_start", type=Path, help="d02 start-state wrfout file")
    parser.add_argument("candidate_end", type=Path, help="TyWRF candidate-end wrfout file")
    parser.add_argument("reference_end", type=Path, help="WRF reference-end wrfout file")
    parser.add_argument(
        "--d01-parent-start",
        type=Path,
        help="Optional d01 parent start wrfout for parent evolution context",
    )
    parser.add_argument(
        "--d01-parent-end",
        type=Path,
        help="Optional d01 parent end wrfout for parent evolution context",
    )
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
        help="Time index used for candidate-end/reference-end Time-leading variables",
    )
    parser.add_argument(
        "--start-time-index",
        type=int,
        default=-1,
        help="Time index used for d02-start and parent-start Time-leading variables",
    )
    parser.add_argument("--output", type=Path, help="Optional JSON output path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_selected_field_wind_source_attribution(
            args.d02_start,
            args.candidate_end,
            args.reference_end,
            d01_parent_start_path=args.d01_parent_start,
            d01_parent_end_path=args.d01_parent_end,
            domain=args.domain,
            variables=args.variables,
            time_index=args.time_index,
            start_time_index=args.start_time_index,
        )
    except (OSError, ValueError, IndexError) as exc:
        parser.error(str(exc))

    payload = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
