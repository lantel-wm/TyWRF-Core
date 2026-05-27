#!/usr/bin/env python
"""Combine KROSA 10 minute moving-nest and field-delta diagnostics."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
from typing import Any, Iterable

import netCDF4

try:
    from tools.analyze_cycle_delta import (
        DEFAULT_THRESHOLDS as DEFAULT_DELTA_THRESHOLDS,
        DEFAULT_VARIABLES as DEFAULT_DELTA_VARIABLES,
        analyze_cycle_delta,
        parse_threshold_overrides,
        report_to_dict as delta_report_to_dict,
        resolve_cycle_files,
    )
    from tools.audit_moving_nest import (
        DEFAULT_DOMAIN,
        DEFAULT_PARENT_GRID_RATIO,
        audit_moving_nest,
        report_to_dict as moving_nest_report_to_dict,
        resolve_moving_nest_files,
    )
except ModuleNotFoundError as exc:
    if exc.name != "tools":
        raise
    from analyze_cycle_delta import (
        DEFAULT_THRESHOLDS as DEFAULT_DELTA_THRESHOLDS,
        DEFAULT_VARIABLES as DEFAULT_DELTA_VARIABLES,
        analyze_cycle_delta,
        parse_threshold_overrides,
        report_to_dict as delta_report_to_dict,
        resolve_cycle_files,
    )
    from audit_moving_nest import (
        DEFAULT_DOMAIN,
        DEFAULT_PARENT_GRID_RATIO,
        audit_moving_nest,
        report_to_dict as moving_nest_report_to_dict,
        resolve_moving_nest_files,
    )


DIAGNOSTIC_ONLY_MESSAGE = (
    "diagnostic-only report; no candidate/model pass is created or evaluated"
)
PARENT_FILL_METADATA_ATTRS = (
    "TYWRF_DIAGNOSTIC_REMAP_PARENT_FILL",
    "TYWRF_MINIMUM_STATIC_REFRESH_FIELDS",
    "TYWRF_STAGGERED_STATIC_COORDS_STATUS",
    "TYWRF_P_DERIVED_REFRESH_STATUS",
    "TYWRF_PRESSURE_REFRESH_REQUIRED",
    "TYWRF_PRESSURE_REFRESH_OPT_IN",
    "TYWRF_PRESSURE_REFRESH_APPLIED",
    "TYWRF_PRESSURE_REFRESH_REQUIREMENT_STATUS",
    "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS",
    "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY",
    "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY",
    "TYWRF_PRESSURE_REFRESH_FORMULA_STATUS",
    "TYWRF_PRESSURE_REFRESH_FORMULA_STAGING_NAME",
    "TYWRF_PRESSURE_REFRESH_REGION_STAGING_NAME",
    "TYWRF_PRESSURE_REFRESH_THERMODYNAMIC_MODE",
    "TYWRF_PRESSURE_REFRESH_REQUIRED_INPUTS",
    "TYWRF_PRESSURE_REFRESH_OUTPUT_FIELD",
    "TYWRF_PRESSURE_REFRESH_HELPER_NAME",
    "TYWRF_PRESSURE_REFRESH_ALB_SOURCE",
    "TYWRF_PRESSURE_REFRESH_PROVIDER_OK",
    "TYWRF_PRESSURE_REFRESH_STAGING_OK",
    "TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED",
    "TYWRF_PRESSURE_REFRESH_TERRAIN_OVERRIDE_USED",
    "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE",
    "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE",
    "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS",
    "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS",
    "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS",
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
    "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE",
    "TYWRF_PRESSURE_REFRESH_METADATA_TIME_INDEX",
    "TYWRF_DIRECT_WRF_END_STATE_ORACLE_STATUS",
    "TYWRF_DIAGNOSTIC_ONLY",
    "TYWRF_GATE_CANDIDATE",
    "TYWRF_INTEGRATOR_OUTPUT",
    "TYWRF_VALIDATION_GATE_ONLY",
    "TYWRF_CANDIDATE_KIND",
)
BOOL_TRUE_VALUES = {"1", "true", "t", "yes", "y", "on"}
BOOL_FALSE_VALUES = {"0", "false", "f", "no", "n", "off"}
BLOCKING_CANDIDATE_KIND_TOKENS = ("diagnostic", "remap", "oracle", "closure")


@dataclass(frozen=True)
class ReportFiles:
    start_file: Path
    end_file: Path
    domain: str
    start_time: str
    end_time: str


@dataclass(frozen=True)
class TenMinuteDiagnosticsReport:
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    disposition: dict[str, Any]
    reference_dir: str
    domain: str
    start_time: str
    end_time: str
    start_file: str
    end_file: str
    candidate_file: str | None
    summary: dict[str, Any]
    candidate_metadata: dict[str, Any]
    parent_fill_metadata: dict[str, Any]
    movement_audit: dict[str, Any]
    field_delta: dict[str, Any]


def resolve_report_files(
    reference_dir: Path,
    *,
    domain: str = DEFAULT_DOMAIN,
    start: str,
    end: str,
) -> ReportFiles:
    """Resolve files through both underlying diagnostic APIs and require agreement."""
    movement_start, movement_end, movement_domain, movement_start_time, movement_end_time = (
        resolve_moving_nest_files(
            reference_dir=reference_dir,
            domain=domain,
            start=start,
            end=end,
        )
    )
    delta_start, delta_end, delta_domain, delta_start_time, delta_end_time = resolve_cycle_files(
        reference_dir=reference_dir,
        domain=domain,
        start=start,
        end=end,
    )

    if (movement_start, movement_end) != (delta_start, delta_end):
        raise ValueError(
            "moving-nest and cycle-delta diagnostics resolved different WRF files"
        )
    if (movement_domain, movement_start_time, movement_end_time) != (
        delta_domain,
        delta_start_time,
        delta_end_time,
    ):
        raise ValueError(
            "moving-nest and cycle-delta diagnostics resolved different segment metadata"
        )
    if movement_start_time is None or movement_end_time is None:
        raise ValueError("start and end times are required for the 10 minute report")

    return ReportFiles(
        start_file=movement_start,
        end_file=movement_end,
        domain=movement_domain,
        start_time=movement_start_time,
        end_time=movement_end_time,
    )


def report_10min_diagnostics(
    reference_dir: Path,
    *,
    domain: str = DEFAULT_DOMAIN,
    start: str,
    end: str,
    variables: Iterable[str] = DEFAULT_DELTA_VARIABLES,
    thresholds: dict[str, float] | None = DEFAULT_DELTA_THRESHOLDS,
    parent_grid_ratio: int = DEFAULT_PARENT_GRID_RATIO,
    log_file: Path | None = None,
    candidate_file: Path | str | None = None,
) -> TenMinuteDiagnosticsReport:
    candidate_path = Path(candidate_file) if candidate_file is not None else None
    files = resolve_report_files(
        reference_dir,
        domain=domain,
        start=start,
        end=end,
    )
    movement = audit_moving_nest(
        files.start_file,
        files.end_file,
        domain=files.domain,
        start_time=files.start_time,
        end_time=files.end_time,
        parent_grid_ratio=parent_grid_ratio,
        log_file=log_file,
    )
    delta = analyze_cycle_delta(
        files.start_file,
        files.end_file,
        variables=variables,
        thresholds=thresholds,
        domain=files.domain,
        start_time=files.start_time,
        end_time=files.end_time,
    )

    movement_payload = moving_nest_report_to_dict(movement)
    delta_payload = delta_report_to_dict(delta)
    candidate_metadata = _read_candidate_metadata(candidate_path)
    parent_fill_metadata = _parent_fill_metadata(candidate_metadata)
    summary = _combined_summary(
        movement_payload,
        delta_payload,
        candidate_metadata,
        parent_fill_metadata,
    )
    return TenMinuteDiagnosticsReport(
        status=_combined_status(movement.status, delta.status),
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        disposition={
            "mode": "diagnostic-only",
            "creates_model_candidate": False,
            "candidate_model_pass": "not_applicable",
            "candidate_gate_eligible": candidate_metadata.get(
                "candidate_gate_eligible"
            ),
            "candidate_gate_blockers": candidate_metadata.get(
                "candidate_gate_blockers"
            ),
            "pressure_refresh_disposition": parent_fill_metadata.get(
                "pressure_refresh_disposition"
            ),
            "message": DIAGNOSTIC_ONLY_MESSAGE,
        },
        reference_dir=str(reference_dir),
        domain=files.domain,
        start_time=files.start_time,
        end_time=files.end_time,
        start_file=str(files.start_file),
        end_file=str(files.end_file),
        candidate_file=str(candidate_path) if candidate_path is not None else None,
        summary=summary,
        candidate_metadata=candidate_metadata,
        parent_fill_metadata=parent_fill_metadata,
        movement_audit=movement_payload,
        field_delta=delta_payload,
    )


def _combined_status(movement_status: str, delta_status: str) -> str:
    if movement_status == "computed" and delta_status == "computed":
        return "computed"
    if "incomplete" in {movement_status, delta_status}:
        return "incomplete"
    return "computed_with_flags"


def _combined_summary(
    movement_payload: dict[str, Any],
    delta_payload: dict[str, Any],
    candidate_metadata: dict[str, Any],
    parent_fill_metadata: dict[str, Any],
) -> dict[str, Any]:
    movement_summary = movement_payload.get("summary", {})
    delta_summary = delta_payload.get("summary", {})
    resolution = movement_summary.get("resolution", {})
    return {
        "moved": movement_summary.get("moved"),
        "parent_delta": movement_summary.get("parent_delta"),
        "child_cell_delta": movement_summary.get("child_cell_delta"),
        "first_failing_variable": delta_summary.get("first_failing_variable"),
        "strict_threshold_exceeded_count": delta_summary.get(
            "strict_threshold_exceeded"
        ),
        "largest_normalized_deltas": delta_summary.get("largest_normalized_deltas"),
        "d02_2km": resolution.get("d02_2km") if isinstance(resolution, dict) else None,
        "d02_2km_status": (
            resolution.get("status") if isinstance(resolution, dict) else None
        ),
        "d02_resolution": resolution,
        "candidate_file": candidate_metadata.get("path"),
        "candidate_metadata_status": candidate_metadata.get("status"),
        "candidate_gate_candidate": candidate_metadata.get("gate_candidate"),
        "candidate_gate_eligible": candidate_metadata.get("candidate_gate_eligible"),
        "candidate_gate_blockers": candidate_metadata.get("candidate_gate_blockers"),
        "candidate_diagnostic_only": candidate_metadata.get(
            "candidate_diagnostic_only"
        ),
        "candidate_integrator_output": candidate_metadata.get("integrator_output"),
        "candidate_validation_gate_only": candidate_metadata.get(
            "validation_gate_only"
        ),
        "candidate_kind": candidate_metadata.get("candidate_kind"),
        "parent_fill_metadata_status": parent_fill_metadata.get("status"),
        "pressure_refresh_metadata_status": parent_fill_metadata.get(
            "pressure_refresh_metadata_status"
        ),
        "diagnostic_remap_parent_fill": parent_fill_metadata.get(
            "diagnostic_remap_parent_fill"
        ),
        "minimum_static_refresh_fields": parent_fill_metadata.get(
            "minimum_static_refresh_fields"
        ),
        "staggered_static_coords_status": parent_fill_metadata.get(
            "staggered_static_coords_status"
        ),
        "p_derived_refresh_status": parent_fill_metadata.get(
            "p_derived_refresh_status"
        ),
        "pressure_refresh_required": parent_fill_metadata.get(
            "pressure_refresh_required"
        ),
        "pressure_refresh_opt_in": parent_fill_metadata.get(
            "pressure_refresh_opt_in"
        ),
        "pressure_refresh_applied": parent_fill_metadata.get(
            "pressure_refresh_applied"
        ),
        "pressure_refresh_requirement_status": parent_fill_metadata.get(
            "pressure_refresh_requirement_status"
        ),
        "pressure_refresh_integration_status": parent_fill_metadata.get(
            "pressure_refresh_integration_status"
        ),
        "pressure_refresh_experimental_apply": parent_fill_metadata.get(
            "pressure_refresh_experimental_apply"
        ),
        "experimental_pressure_refresh_apply": parent_fill_metadata.get(
            "experimental_pressure_refresh_apply"
        ),
        "normal_candidate_pressure_refresh": parent_fill_metadata.get(
            "normal_candidate_pressure_refresh"
        ),
        "hidden_seam_pressure_refresh": parent_fill_metadata.get(
            "hidden_seam_pressure_refresh"
        ),
        "pressure_refresh_formula_status": parent_fill_metadata.get(
            "pressure_refresh_formula_status"
        ),
        "pressure_refresh_formula_staging_name": parent_fill_metadata.get(
            "pressure_refresh_formula_staging_name"
        ),
        "pressure_refresh_region_staging_name": parent_fill_metadata.get(
            "pressure_refresh_region_staging_name"
        ),
        "pressure_refresh_thermodynamic_mode": parent_fill_metadata.get(
            "pressure_refresh_thermodynamic_mode"
        ),
        "pressure_refresh_required_inputs": parent_fill_metadata.get(
            "pressure_refresh_required_inputs"
        ),
        "pressure_refresh_output_field": parent_fill_metadata.get(
            "pressure_refresh_output_field"
        ),
        "pressure_refresh_helper_name": parent_fill_metadata.get(
            "pressure_refresh_helper_name"
        ),
        "pressure_refresh_alb_source": parent_fill_metadata.get(
            "pressure_refresh_alb_source"
        ),
        "pressure_refresh_provider_ok": parent_fill_metadata.get(
            "pressure_refresh_provider_ok"
        ),
        "pressure_refresh_staging_ok": parent_fill_metadata.get(
            "pressure_refresh_staging_ok"
        ),
        "pressure_refresh_compute_called": parent_fill_metadata.get(
            "pressure_refresh_compute_called"
        ),
        "pressure_refresh_terrain_override_used": parent_fill_metadata.get(
            "pressure_refresh_terrain_override_used"
        ),
        "pressure_refresh_terrain_source": parent_fill_metadata.get(
            "pressure_refresh_terrain_source"
        ),
        "pressure_refresh_terrain_provenance": parent_fill_metadata.get(
            "pressure_refresh_terrain_provenance"
        ),
        "pressure_refresh_synced_pb_points": parent_fill_metadata.get(
            "pressure_refresh_synced_pb_points"
        ),
        "pressure_refresh_synced_mub_points": parent_fill_metadata.get(
            "pressure_refresh_synced_mub_points"
        ),
        "pressure_refresh_synced_phb_points": parent_fill_metadata.get(
            "pressure_refresh_synced_phb_points"
        ),
        "pressure_refresh_target_column_count": parent_fill_metadata.get(
            "pressure_refresh_target_column_count"
        ),
        "pressure_refresh_refreshed_column_count": parent_fill_metadata.get(
            "pressure_refresh_refreshed_column_count"
        ),
        "pressure_refresh_refreshed_point_count": parent_fill_metadata.get(
            "pressure_refresh_refreshed_point_count"
        ),
        "pressure_refresh_skipped_point_count": parent_fill_metadata.get(
            "pressure_refresh_skipped_point_count"
        ),
        "pressure_refresh_invalid_point_count": parent_fill_metadata.get(
            "pressure_refresh_invalid_point_count"
        ),
        "pressure_refresh_touched_overlap_cells": parent_fill_metadata.get(
            "pressure_refresh_touched_overlap_cells"
        ),
        "pressure_refresh_touched_halo_cells": parent_fill_metadata.get(
            "pressure_refresh_touched_halo_cells"
        ),
        "pressure_refresh_refreshed_p_points": parent_fill_metadata.get(
            "pressure_refresh_refreshed_p_points"
        ),
        "pressure_refresh_changed_p_points": parent_fill_metadata.get(
            "pressure_refresh_changed_p_points"
        ),
        "pressure_refresh_changed_pb_points": parent_fill_metadata.get(
            "pressure_refresh_changed_pb_points"
        ),
        "pressure_refresh_changed_mub_points": parent_fill_metadata.get(
            "pressure_refresh_changed_mub_points"
        ),
        "pressure_refresh_changed_phb_points": parent_fill_metadata.get(
            "pressure_refresh_changed_phb_points"
        ),
        "pressure_refresh_changed_p_matches_refreshed_point_count": (
            parent_fill_metadata.get(
                "pressure_refresh_changed_p_matches_refreshed_point_count"
            )
        ),
        "pressure_refresh_invalid_and_skipped_points_zero": parent_fill_metadata.get(
            "pressure_refresh_invalid_and_skipped_points_zero"
        ),
        "pressure_refresh_overlap_halo_untouched": parent_fill_metadata.get(
            "pressure_refresh_overlap_halo_untouched"
        ),
        "pressure_refresh_metadata_source": parent_fill_metadata.get(
            "pressure_refresh_metadata_source"
        ),
        "pressure_refresh_metadata_time_index": parent_fill_metadata.get(
            "pressure_refresh_metadata_time_index"
        ),
        "pressure_refresh_disposition": parent_fill_metadata.get(
            "pressure_refresh_disposition"
        ),
        "direct_wrf_end_state_oracle_status": parent_fill_metadata.get(
            "direct_wrf_end_state_oracle_status"
        ),
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "message": DIAGNOSTIC_ONLY_MESSAGE,
    }


def _read_candidate_metadata(candidate_file: Path | None) -> dict[str, Any]:
    base = {
        "path": str(candidate_file) if candidate_file is not None else None,
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
    }
    if candidate_file is None:
        return {
            **base,
            "status": "not_provided",
            "attrs": {},
            "present_attrs": [],
            "missing_attrs": [],
            "gate_candidate": None,
            "candidate_diagnostic_only": None,
            "integrator_output": None,
            "validation_gate_only": None,
            "candidate_kind": None,
            "candidate_gate_eligible": False,
            "candidate_gate_blockers": ["candidate_file_not_provided"],
            "message": "no candidate file supplied",
        }
    if not candidate_file.exists():
        return {
            **base,
            "status": "not_available",
            "attrs": {},
            "present_attrs": [],
            "missing_attrs": list(PARENT_FILL_METADATA_ATTRS),
            "gate_candidate": None,
            "candidate_diagnostic_only": None,
            "integrator_output": None,
            "validation_gate_only": None,
            "candidate_kind": None,
            "candidate_gate_eligible": False,
            "candidate_gate_blockers": ["candidate_file_not_available"],
            "message": f"candidate file does not exist: {candidate_file}",
        }

    try:
        with netCDF4.Dataset(candidate_file) as dataset:
            ncattrs = set(dataset.ncattrs())
            attrs = {
                name: _netcdf_attr_to_json(dataset.getncattr(name))
                for name in PARENT_FILL_METADATA_ATTRS
                if name in ncattrs
            }
    except OSError as exc:
        return {
            **base,
            "status": "not_available",
            "attrs": {},
            "present_attrs": [],
            "missing_attrs": list(PARENT_FILL_METADATA_ATTRS),
            "gate_candidate": None,
            "candidate_diagnostic_only": None,
            "integrator_output": None,
            "validation_gate_only": None,
            "candidate_kind": None,
            "candidate_gate_eligible": False,
            "candidate_gate_blockers": ["candidate_metadata_unreadable"],
            "message": f"candidate metadata check failed: {exc}",
        }

    result = {
        **base,
        "status": "available",
        "attrs": attrs,
        "present_attrs": sorted(attrs),
        "missing_attrs": [
            name for name in PARENT_FILL_METADATA_ATTRS if name not in attrs
        ],
        "gate_candidate": _coerce_bool(attrs.get("TYWRF_GATE_CANDIDATE")),
        "candidate_diagnostic_only": _coerce_bool(
            attrs.get("TYWRF_DIAGNOSTIC_ONLY")
        ),
        "integrator_output": _coerce_bool(attrs.get("TYWRF_INTEGRATOR_OUTPUT")),
        "validation_gate_only": _coerce_bool(
            attrs.get("TYWRF_VALIDATION_GATE_ONLY")
        ),
        "candidate_kind": attrs.get("TYWRF_CANDIDATE_KIND"),
        "message": (
            "candidate metadata surfaced for diagnostic context only; "
            "no model pass is created"
        ),
    }
    return {**result, **_candidate_gate_assessment(result)}


def _candidate_gate_assessment(candidate_metadata: dict[str, Any]) -> dict[str, Any]:
    status = candidate_metadata.get("status")
    attrs = candidate_metadata.get("attrs", {})
    if status == "not_provided":
        blockers = ["candidate_file_not_provided"]
    elif status == "not_available":
        blockers = ["candidate_file_not_available"]
    elif status != "available":
        blockers = [f"candidate_metadata_status={status}"]
    else:
        blockers = []

    if candidate_metadata.get("candidate_diagnostic_only") is True:
        blockers.append("TYWRF_DIAGNOSTIC_ONLY=true")
    if candidate_metadata.get("gate_candidate") is False:
        blockers.append("TYWRF_GATE_CANDIDATE=false")
    elif candidate_metadata.get("gate_candidate") is not True:
        blockers.append("TYWRF_GATE_CANDIDATE is not true")
    if candidate_metadata.get("integrator_output") is False:
        blockers.append("TYWRF_INTEGRATOR_OUTPUT=false")
    elif candidate_metadata.get("integrator_output") is not True:
        blockers.append("TYWRF_INTEGRATOR_OUTPUT is not true")
    if candidate_metadata.get("validation_gate_only") is True:
        blockers.append("TYWRF_VALIDATION_GATE_ONLY=true")

    candidate_kind = candidate_metadata.get("candidate_kind")
    kind = str(candidate_kind).lower() if candidate_kind else ""
    if kind and any(token in kind for token in BLOCKING_CANDIDATE_KIND_TOKENS):
        blockers.append(f"TYWRF_CANDIDATE_KIND={candidate_kind}")

    if _coerce_bool(attrs.get("TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY")) is True:
        blockers.append("TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY=true")
    if _coerce_bool(attrs.get("TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY")) is True:
        blockers.append("TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY=true")
    integration_status = attrs.get("TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS")
    if str(integration_status).strip().lower() == "experimental_apply_test_only":
        blockers.append(
            "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS=experimental_apply_test_only"
        )

    return {
        "candidate_gate_eligible": status == "available" and not blockers,
        "candidate_gate_blockers": blockers,
    }


def _parent_fill_metadata(candidate_metadata: dict[str, Any]) -> dict[str, Any]:
    attrs = candidate_metadata.get("attrs", {})
    parent_fill = _coerce_bool(attrs.get("TYWRF_DIAGNOSTIC_REMAP_PARENT_FILL"))
    pressure_refresh_metadata_available = any(
        name.startswith("TYWRF_PRESSURE_REFRESH_")
        or name == "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY"
        for name in attrs
    )
    if candidate_metadata.get("status") != "available":
        status = candidate_metadata.get("status")
    elif parent_fill is True:
        status = "available"
    elif parent_fill is False:
        status = "not_parent_fill"
    else:
        status = "missing_parent_fill_flag"

    result = {
        "status": status,
        "path": candidate_metadata.get("path"),
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "candidate_gate_eligible": candidate_metadata.get("candidate_gate_eligible"),
        "candidate_gate_blockers": candidate_metadata.get("candidate_gate_blockers"),
        "pressure_refresh_metadata_status": (
            "available" if pressure_refresh_metadata_available else "not_reported"
        ),
        "diagnostic_remap_parent_fill": parent_fill,
        "minimum_static_refresh_fields": _split_csv_attr(
            attrs.get("TYWRF_MINIMUM_STATIC_REFRESH_FIELDS")
        ),
        "staggered_static_coords_status": attrs.get(
            "TYWRF_STAGGERED_STATIC_COORDS_STATUS"
        ),
        "p_derived_refresh_status": attrs.get("TYWRF_P_DERIVED_REFRESH_STATUS"),
        "pressure_refresh_required": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_REQUIRED")
        ),
        "pressure_refresh_opt_in": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_OPT_IN")
        ),
        "pressure_refresh_applied": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_APPLIED")
        ),
        "pressure_refresh_requirement_status": attrs.get(
            "TYWRF_PRESSURE_REFRESH_REQUIREMENT_STATUS"
        ),
        "pressure_refresh_integration_status": attrs.get(
            "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS"
        ),
        "pressure_refresh_experimental_apply": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY")
        ),
        "experimental_pressure_refresh_apply": _coerce_bool(
            attrs.get("TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY")
        ),
        "pressure_refresh_formula_status": attrs.get(
            "TYWRF_PRESSURE_REFRESH_FORMULA_STATUS"
        ),
        "pressure_refresh_formula_staging_name": attrs.get(
            "TYWRF_PRESSURE_REFRESH_FORMULA_STAGING_NAME"
        ),
        "pressure_refresh_region_staging_name": attrs.get(
            "TYWRF_PRESSURE_REFRESH_REGION_STAGING_NAME"
        ),
        "pressure_refresh_thermodynamic_mode": attrs.get(
            "TYWRF_PRESSURE_REFRESH_THERMODYNAMIC_MODE"
        ),
        "pressure_refresh_required_inputs": _split_csv_attr(
            attrs.get("TYWRF_PRESSURE_REFRESH_REQUIRED_INPUTS")
        ),
        "pressure_refresh_output_field": attrs.get(
            "TYWRF_PRESSURE_REFRESH_OUTPUT_FIELD"
        ),
        "pressure_refresh_helper_name": attrs.get(
            "TYWRF_PRESSURE_REFRESH_HELPER_NAME"
        ),
        "pressure_refresh_alb_source": attrs.get(
            "TYWRF_PRESSURE_REFRESH_ALB_SOURCE"
        ),
        "pressure_refresh_provider_ok": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_PROVIDER_OK")
        ),
        "pressure_refresh_staging_ok": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_STAGING_OK")
        ),
        "pressure_refresh_compute_called": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED")
        ),
        "pressure_refresh_terrain_override_used": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_TERRAIN_OVERRIDE_USED")
        ),
        "pressure_refresh_terrain_source": attrs.get(
            "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE"
        ),
        "pressure_refresh_terrain_provenance": attrs.get(
            "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE"
        ),
        "pressure_refresh_synced_pb_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS")
        ),
        "pressure_refresh_synced_mub_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS")
        ),
        "pressure_refresh_synced_phb_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS")
        ),
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
        "pressure_refresh_metadata_source": attrs.get(
            "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE"
        ),
        "pressure_refresh_metadata_time_index": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_METADATA_TIME_INDEX")
        ),
        "direct_wrf_end_state_oracle_status": attrs.get(
            "TYWRF_DIRECT_WRF_END_STATE_ORACLE_STATUS"
        ),
        "gate_candidate": candidate_metadata.get("gate_candidate"),
    }
    result["hidden_seam_pressure_refresh"] = _hidden_seam_pressure_refresh(result)
    result["normal_candidate_pressure_refresh"] = _normal_candidate_pressure_refresh(
        result,
        candidate_metadata,
    )
    result["pressure_refresh_disposition"] = _pressure_refresh_disposition(
        result,
        candidate_metadata,
    )
    return result


def _hidden_seam_pressure_refresh(parent_fill_metadata: dict[str, Any]) -> bool:
    integration_status = parent_fill_metadata.get("pressure_refresh_integration_status")
    return (
        parent_fill_metadata.get("pressure_refresh_experimental_apply") is True
        or parent_fill_metadata.get("experimental_pressure_refresh_apply") is True
        or str(integration_status).strip().lower() == "experimental_apply_test_only"
    )


def _normal_candidate_pressure_refresh(
    parent_fill_metadata: dict[str, Any],
    candidate_metadata: dict[str, Any],
) -> bool:
    return (
        parent_fill_metadata.get("pressure_refresh_applied") is True
        and parent_fill_metadata.get("pressure_refresh_integration_status")
        == "applied_to_candidate"
        and candidate_metadata.get("candidate_gate_eligible") is True
        and candidate_metadata.get("candidate_diagnostic_only") is False
        and candidate_metadata.get("gate_candidate") is True
        and candidate_metadata.get("integrator_output") is True
        and candidate_metadata.get("validation_gate_only") is not True
        and not _hidden_seam_pressure_refresh(parent_fill_metadata)
    )


def _pressure_refresh_disposition(
    parent_fill_metadata: dict[str, Any],
    candidate_metadata: dict[str, Any],
) -> dict[str, Any]:
    applied = parent_fill_metadata.get("pressure_refresh_applied")
    gate_eligible = bool(candidate_metadata.get("candidate_gate_eligible"))
    blockers = candidate_metadata.get("candidate_gate_blockers") or []
    hidden_seam = parent_fill_metadata.get("hidden_seam_pressure_refresh") is True
    normal_candidate = (
        parent_fill_metadata.get("normal_candidate_pressure_refresh") is True
    )

    if applied is True and hidden_seam:
        status = "hidden_seam_diagnostic_evidence_only"
        message = (
            "pressure refresh came from an experimental hidden seam; "
            "it is not a model validation gate pass"
        )
    elif applied is True and normal_candidate:
        status = "normal_candidate_metadata_only"
        message = (
            "pressure_refresh_applied=true is attached to a gate-eligible "
            "integrator candidate, but this report remains diagnostic-only"
        )
    elif applied is True and not gate_eligible:
        status = "diagnostic_helper_evidence_only"
        message = (
            "pressure_refresh_applied=true is diagnostic helper evidence only; "
            "it is not a model validation gate pass"
        )
    elif applied is True:
        status = "candidate_metadata_only"
        message = (
            "pressure_refresh_applied=true is candidate metadata only; "
            "candidate_model_pass remains not_applicable in this diagnostic report"
        )
    elif applied is False:
        status = "not_applied"
        message = "pressure_refresh_applied=false"
    else:
        status = "not_reported"
        message = "pressure_refresh_applied metadata is not present"

    return {
        "status": status,
        "pressure_refresh_applied": applied,
        "pressure_refresh_integration_status": parent_fill_metadata.get(
            "pressure_refresh_integration_status"
        ),
        "pressure_refresh_experimental_apply": parent_fill_metadata.get(
            "pressure_refresh_experimental_apply"
        ),
        "experimental_pressure_refresh_apply": parent_fill_metadata.get(
            "experimental_pressure_refresh_apply"
        ),
        "normal_candidate_pressure_refresh": normal_candidate,
        "hidden_seam_pressure_refresh": hidden_seam,
        "candidate_gate_eligible": gate_eligible,
        "candidate_gate_blockers": blockers,
        "candidate_model_pass": "not_applicable",
        "message": message,
    }


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


def _coerce_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
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
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if value.is_integer():
            return int(value)
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _split_csv_attr(value: Any) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _strict_json_value(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _strict_json_value(asdict(value))
    if isinstance(value, dict):
        return {key: _strict_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_strict_json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def report_to_dict(report: TenMinuteDiagnosticsReport) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(report: TenMinuteDiagnosticsReport, *, pretty: bool = False) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-dir",
        type=Path,
        required=True,
        help="Directory with WRF reference files",
    )
    parser.add_argument("--domain", choices=(DEFAULT_DOMAIN,), default=DEFAULT_DOMAIN)
    parser.add_argument(
        "--candidate-file",
        type=Path,
        help="Optional TyWRF diagnostic output file to mine for candidate metadata",
    )
    parser.add_argument(
        "--start",
        required=True,
        help="Cycle start time, for example 2025-07-26_00:00:00",
    )
    parser.add_argument(
        "--end",
        required=True,
        help="Cycle end time, for example 2025-07-26_00:10:00",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_DELTA_VARIABLES),
        help="Variables to analyze; defaults to strict field gate variables",
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
        help="Report deltas without strict-threshold exceedance flags",
    )
    parser.add_argument(
        "--parent-grid-ratio",
        type=int,
        default=DEFAULT_PARENT_GRID_RATIO,
        help="Parent-to-child horizontal grid ratio; KROSA d02 uses 5",
    )
    parser.add_argument("--log-file", type=Path, help="Optional WRF rsl log to mine for move lines")
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        thresholds = None if args.no_thresholds else parse_threshold_overrides(args.threshold)
        report = report_10min_diagnostics(
            args.reference_dir,
            domain=args.domain,
            start=args.start,
            end=args.end,
            variables=args.variables,
            thresholds=thresholds,
            parent_grid_ratio=args.parent_grid_ratio,
            log_file=args.log_file,
            candidate_file=args.candidate_file,
        )
    except ValueError as exc:
        parser.error(str(exc))

    report_json = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
