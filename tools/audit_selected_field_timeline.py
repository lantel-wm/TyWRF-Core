#!/usr/bin/env python
"""Diagnostic-only audit of selected-field moving-nest timeline/remap evidence."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
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

from tools.audit_selected_field_pipeline import (
    summarize_evolution_json,
    summarize_spatial_alignment_json,
)
from tools.audit_selected_field_state import (
    _coerce_bool,
    _coerce_int,
    _coerce_pair,
    _field_region_mask,
    _netcdf_attr_to_json,
    _pair_to_dict,
)


DEFAULT_SOURCE_FILES = (
    "src/tools/selected_field_cycle.cpp",
    "src/dynamics/cycle_schedule.cpp",
    "include/tywrf/dynamics/cycle_schedule.hpp",
    "src/nest/state_remap.cpp",
    "src/nest/parent_child_interpolation.cpp",
    "src/nest/state_exchange.cpp",
    "src/nest/static_fields.cpp",
)
DEFAULT_TIMELINE_VARIABLES = ("U", "V", "T", "PH", "MU", "P", "QVAPOR")
STRICT_SELECTED_STATE_VARIABLES = frozenset(DEFAULT_TIMELINE_VARIABLES)
LARGE_CHILD_DELTA_I = 50
LARGE_CHILD_DELTA_J = 30
SMALL_SHIFT_LIMIT = 5

SOURCE_HOOKS = {
    "selected_field_tool": (
        "build_candidate_state",
        "remap_child_state_overlap_only",
        "build_exposed_child_state_exchange_plan",
        "interpolate_parent_to_exposed_child",
        "refresh_static_fields",
        "probe_pressure_refresh_provider_readiness",
        "probe_pressure_refresh_dry_run_contract",
        "apply_pressure_refresh",
        "stamp_gate_metadata",
    ),
    "schedule": (
        "CycleScheduleCallKind::moving_nest_move_check",
        "CycleScheduleCallKind::moving_nest_post_move_parent_fill",
        "CycleScheduleCallKind::parent_child_interpolation",
        "CycleScheduleCallKind::two_way_feedback",
        "snap_to_next_parent_step",
        "make_krosa_10min_cycle_schedule_config",
    ),
    "remap": (
        "build_remap_plan",
        "remap_child_state_overlap_only",
        "remap_child_state_overlap_with_parent_fill",
        "needs_parent_fill",
        "needs_derived_pressure_refresh",
    ),
    "state_exchange": (
        "build_exposed_child_state_exchange_plan",
        "StateExchangeField::u",
        "StateExchangeField::v",
        "StateExchangeField::mu",
        "StateExchangeField::qvapor",
        "StateExchangeField::t",
        "StateExchangeField::ph",
    ),
    "parent_child_interpolation": (
        "interpolate_parent_to_exposed_child",
        "interpolate_regions_2d",
        "interpolate_regions_3d",
        "ParentChildInterpolationMethod::bilinear",
    ),
    "static_refresh": (
        "refresh_moving_nest_static_fields",
        "local_linear_extrapolated_2d",
        "bilinear_clamped_2d",
        "uses_reference_end",
    ),
}

TIMELINE_ATTRS = (
    "TYWRF_DIAGNOSTIC_ONLY",
    "TYWRF_GATE_CANDIDATE",
    "TYWRF_INTEGRATOR_OUTPUT",
    "TYWRF_VALIDATION_GATE_ONLY",
    "TYWRF_CANDIDATE_KIND",
    "TYWRF_CANDIDATE_DOMAIN",
    "TYWRF_D02_RESOLUTION_CHECK",
    "TYWRF_CYCLE_START",
    "TYWRF_CYCLE_END",
    "TYWRF_FROM_PARENT_START",
    "TYWRF_TO_PARENT_START",
    "TYWRF_PARENT_GRID_RATIO",
    "TYWRF_STATE_VARIABLES",
    "TYWRF_PARENT_INTERPOLATED_STATE_VARIABLES",
    "TYWRF_SELECTED_FIELD_CHANGED_POINTS",
    "TYWRF_EXPOSED_EXCHANGE_POINTS",
    "TYWRF_INTERPOLATED_POINTS",
    "TYWRF_STATIC_REFRESH_APPLIED",
    "TYWRF_STATIC_REFRESH_METHOD",
    "TYWRF_STATIC_REFRESH_USES_REFERENCE_END",
    "TYWRF_STATIC_REFRESH_OVERLAP_CELLS",
    "TYWRF_STATIC_REFRESH_EXPOSED_CELLS",
    "TYWRF_STATIC_REFRESH_COORD_EXTRAPOLATED_CELLS",
    "TYWRF_STATIC_REFRESH_HGT_PARENT_INTERPOLATED_CELLS",
    "TYWRF_STATIC_REFRESH_CHANGED_TEMPLATE_POINTS",
    "TYWRF_PRESSURE_REFRESH_OPT_IN",
    "TYWRF_PRESSURE_REFRESH_APPLIED",
    "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS",
    "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS",
    "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT",
    "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS",
    "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS",
    "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS",
    "I_PARENT_START",
    "J_PARENT_START",
    "DX",
    "DY",
)


def _report_base(**extra: Any) -> dict[str, Any]:
    return {
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        **extra,
    }


def _numeric_attr(value: Any) -> float | None:
    if value is None or isinstance(value, (bool, np.bool_)):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _integerish_attr(value: Any) -> int | None:
    parsed = _coerce_int(value)
    if parsed is not None:
        return parsed
    numeric = _numeric_attr(value)
    if numeric is None or not numeric.is_integer():
        return None
    return int(numeric)


def _split_csv_attr(value: Any) -> list[str]:
    if value is None:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _read_attrs(dataset: netCDF4.Dataset) -> dict[str, Any]:
    return {
        name: _netcdf_attr_to_json(dataset.getncattr(name))
        for name in dataset.ncattrs()
        if name in TIMELINE_ATTRS or name.startswith("TYWRF_")
    }


def _variable_shape(
    dataset: netCDF4.Dataset,
    variable: str,
    *,
    time_index: int,
) -> tuple[int, ...] | None:
    if variable not in dataset.variables:
        return None
    del time_index
    netcdf_variable = dataset.variables[variable]
    shape = tuple(int(value) for value in netcdf_variable.shape)
    dims = netcdf_variable.dimensions
    if dims and dims[0] == "Time":
        return shape[1:]
    return shape


def _mass_horizontal_shape(
    dataset: netCDF4.Dataset,
    variables: Iterable[str],
    *,
    time_index: int,
) -> tuple[int, int] | None:
    for variable in ("MU", "T", "P", "QVAPOR", *variables):
        shape = _variable_shape(dataset, variable, time_index=time_index)
        if shape is not None and len(shape) >= 2:
            return int(shape[-2]), int(shape[-1])
    if "south_north" in dataset.dimensions and "west_east" in dataset.dimensions:
        return (
            int(len(dataset.dimensions["south_north"])),
            int(len(dataset.dimensions["west_east"])),
        )
    return None


def _movement_from_attrs(attrs: dict[str, Any]) -> dict[str, Any]:
    from_parent_start = _coerce_pair(attrs.get("TYWRF_FROM_PARENT_START"))
    to_parent_start = _coerce_pair(attrs.get("TYWRF_TO_PARENT_START"))
    parent_grid_ratio = _integerish_attr(attrs.get("TYWRF_PARENT_GRID_RATIO"))
    missing = [
        name
        for name, value in (
            ("TYWRF_FROM_PARENT_START", from_parent_start),
            ("TYWRF_TO_PARENT_START", to_parent_start),
            ("TYWRF_PARENT_GRID_RATIO", parent_grid_ratio),
        )
        if value is None
    ]
    if missing:
        return _report_base(
            status="not_available",
            missing_metadata=missing,
            from_parent_start=_pair_to_dict(from_parent_start),
            to_parent_start=_pair_to_dict(to_parent_start),
            parent_grid_ratio=parent_grid_ratio,
            child_delta=None,
            large_movement=None,
            message="cannot reconstruct moving-nest delta from candidate metadata",
        )
    if parent_grid_ratio <= 0:
        return _report_base(
            status="not_available",
            missing_metadata=[],
            from_parent_start=_pair_to_dict(from_parent_start),
            to_parent_start=_pair_to_dict(to_parent_start),
            parent_grid_ratio=parent_grid_ratio,
            child_delta=None,
            large_movement=None,
            message="TYWRF_PARENT_GRID_RATIO must be positive",
        )

    delta_i = (to_parent_start[0] - from_parent_start[0]) * parent_grid_ratio
    delta_j = (to_parent_start[1] - from_parent_start[1]) * parent_grid_ratio
    return _report_base(
        status="available",
        from_parent_start=_pair_to_dict(from_parent_start),
        to_parent_start=_pair_to_dict(to_parent_start),
        parent_grid_ratio=parent_grid_ratio,
        child_delta={"i": delta_i, "j": delta_j},
        manhattan_child_cells=abs(delta_i) + abs(delta_j),
        euclidean_child_cells=math.hypot(delta_i, delta_j),
        large_movement=abs(delta_i) >= LARGE_CHILD_DELTA_I
        or abs(delta_j) >= LARGE_CHILD_DELTA_J,
        message=None,
    )


def _approx_overlap_dimensions(
    variable: str,
    horizontal_shape: tuple[int, int],
    mass_shape: tuple[int, int] | None,
    child_delta: dict[str, int] | None,
) -> dict[str, Any]:
    if child_delta is None:
        return _report_base(status="not_available", message="child_delta unavailable")
    delta_i = int(child_delta["i"])
    delta_j = int(child_delta["j"])
    mask_delta_i = abs(delta_i)
    mask_delta_j = abs(delta_j)
    ny, nx = horizontal_shape
    if variable == "U" and mass_shape == (ny, nx - 1) and delta_i != 0:
        mask_delta_i += 1
    if variable == "V" and mass_shape == (ny - 1, nx) and delta_j != 0:
        mask_delta_j += 1
    return _report_base(
        status="available",
        overlap_dimensions={
            "ny": max(0, ny - mask_delta_j),
            "nx": max(0, nx - mask_delta_i),
        },
        exposed_edge_width={"i": mask_delta_i, "j": mask_delta_j},
    )


def _variable_region_report(
    dataset: netCDF4.Dataset,
    variable: str,
    *,
    movement: dict[str, Any],
    mass_shape: tuple[int, int] | None,
    time_index: int,
) -> dict[str, Any]:
    shape = _variable_shape(dataset, variable, time_index=time_index)
    if shape is None or len(shape) < 2:
        return _report_base(
            status="not_available",
            variable=variable,
            candidate_shape=None,
            message="candidate variable missing or has no horizontal shape",
        )

    horizontal_shape = (int(shape[-2]), int(shape[-1]))
    region = {
        "status": movement.get("status"),
        "from_parent_start": movement.get("from_parent_start"),
        "to_parent_start": movement.get("to_parent_start"),
        "parent_grid_ratio": movement.get("parent_grid_ratio"),
        "child_delta": movement.get("child_delta"),
        "message": movement.get("message"),
    }
    _, mask_info = _field_region_mask(variable, horizontal_shape, mass_shape, region)
    target = _numeric_attr(mask_info.get("target_point_count_2d"))
    overlap = _numeric_attr(mask_info.get("overlap_point_count_2d"))
    if target is not None and overlap is not None and target + overlap > 0.0:
        target_fraction = target / (target + overlap)
    else:
        target_fraction = None
    return _report_base(
        **mask_info,
        variable=variable,
        candidate_shape=list(shape),
        target_region_fraction_2d=target_fraction,
        approximate_geometry=_approx_overlap_dimensions(
            variable,
            horizontal_shape,
            mass_shape,
            movement.get("child_delta"),
        ),
    )


def _static_refresh_summary(attrs: dict[str, Any]) -> dict[str, Any]:
    return _report_base(
        status="available",
        applied=_coerce_bool(attrs.get("TYWRF_STATIC_REFRESH_APPLIED")),
        method=attrs.get("TYWRF_STATIC_REFRESH_METHOD"),
        uses_reference_end=_coerce_bool(
            attrs.get("TYWRF_STATIC_REFRESH_USES_REFERENCE_END")
        ),
        overlap_cells=_numeric_attr(attrs.get("TYWRF_STATIC_REFRESH_OVERLAP_CELLS")),
        exposed_cells=_numeric_attr(attrs.get("TYWRF_STATIC_REFRESH_EXPOSED_CELLS")),
        coordinate_extrapolated_cells=_numeric_attr(
            attrs.get("TYWRF_STATIC_REFRESH_COORD_EXTRAPOLATED_CELLS")
        ),
        hgt_parent_interpolated_cells=_numeric_attr(
            attrs.get("TYWRF_STATIC_REFRESH_HGT_PARENT_INTERPOLATED_CELLS")
        ),
        changed_template_points=_numeric_attr(
            attrs.get("TYWRF_STATIC_REFRESH_CHANGED_TEMPLATE_POINTS")
        ),
    )


def _pressure_refresh_summary(attrs: dict[str, Any]) -> dict[str, Any]:
    return _report_base(
        status="available",
        opt_in=_coerce_bool(attrs.get("TYWRF_PRESSURE_REFRESH_OPT_IN")),
        applied=_coerce_bool(attrs.get("TYWRF_PRESSURE_REFRESH_APPLIED")),
        integration_status=attrs.get("TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS"),
        changed_p_points=_numeric_attr(attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS")),
        changed_pb_points=_numeric_attr(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS")
        ),
        changed_mub_points=_numeric_attr(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS")
        ),
        changed_phb_points=_numeric_attr(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS")
        ),
        refreshed_point_count=_numeric_attr(
            attrs.get("TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT")
        ),
        synced_pb_points=_numeric_attr(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS")
        ),
        synced_mub_points=_numeric_attr(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS")
        ),
        synced_phb_points=_numeric_attr(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS")
        ),
    )


def candidate_timeline(
    candidate_path: Path,
    *,
    variables: Iterable[str],
    time_index: int,
) -> dict[str, Any]:
    variables = tuple(variables)
    with netCDF4.Dataset(candidate_path) as dataset:
        attrs = _read_attrs(dataset)
        movement = _movement_from_attrs(attrs)
        mass_shape = _mass_horizontal_shape(dataset, variables, time_index=time_index)
        variable_regions = {
            variable: _variable_region_report(
                dataset,
                variable,
                movement=movement,
                mass_shape=mass_shape,
                time_index=time_index,
            )
            for variable in variables
        }
        variable_shapes = {
            variable: (
                None
                if (shape := _variable_shape(dataset, variable, time_index=time_index))
                is None
                else list(shape)
            )
            for variable in variables
        }

    state_variables = _split_csv_attr(attrs.get("TYWRF_STATE_VARIABLES")) or list(variables)
    strict_state_variables = [
        variable
        for variable in state_variables
        if variable in STRICT_SELECTED_STATE_VARIABLES
    ]
    parent_interpolated = _split_csv_attr(
        attrs.get("TYWRF_PARENT_INTERPOLATED_STATE_VARIABLES")
    )
    missing_parent_interpolation = sorted(
        set(strict_state_variables) - set(parent_interpolated)
    )
    changed_points = _numeric_attr(attrs.get("TYWRF_SELECTED_FIELD_CHANGED_POINTS"))
    interpolated_points = _numeric_attr(attrs.get("TYWRF_INTERPOLATED_POINTS"))
    count_mismatch = (
        changed_points is not None
        and interpolated_points is not None
        and not math.isclose(changed_points, interpolated_points, rel_tol=0.0, abs_tol=0.0)
    )

    return _report_base(
        status="available",
        candidate=str(candidate_path),
        cycle=_report_base(
            start=attrs.get("TYWRF_CYCLE_START"),
            end=attrs.get("TYWRF_CYCLE_END"),
        ),
        disposition=_report_base(
            gate_candidate=_coerce_bool(attrs.get("TYWRF_GATE_CANDIDATE")),
            integrator_output=_coerce_bool(attrs.get("TYWRF_INTEGRATOR_OUTPUT")),
            candidate_diagnostic_only=_coerce_bool(attrs.get("TYWRF_DIAGNOSTIC_ONLY")),
            validation_gate_only=_coerce_bool(attrs.get("TYWRF_VALIDATION_GATE_ONLY")),
            candidate_kind=attrs.get("TYWRF_CANDIDATE_KIND"),
            candidate_domain=attrs.get("TYWRF_CANDIDATE_DOMAIN"),
            d02_resolution_check=attrs.get("TYWRF_D02_RESOLUTION_CHECK"),
            dx_m=_numeric_attr(attrs.get("DX")),
            dy_m=_numeric_attr(attrs.get("DY")),
        ),
        movement=movement,
        mass_horizontal_shape=(
            None
            if mass_shape is None
            else {"ny": int(mass_shape[0]), "nx": int(mass_shape[1])}
        ),
        variable_shapes=variable_shapes,
        exposed_overlap=variable_regions,
        counts=_report_base(
            selected_field_changed_points=changed_points,
            exposed_exchange_points=_numeric_attr(attrs.get("TYWRF_EXPOSED_EXCHANGE_POINTS")),
            interpolated_points=interpolated_points,
            changed_interpolated_point_count_mismatch=count_mismatch,
        ),
        selected_state_variables=strict_state_variables,
        parent_interpolated_state_variables=parent_interpolated,
        selected_state_variables_without_parent_interpolation=missing_parent_interpolation,
        static_refresh=_static_refresh_summary(attrs),
        pressure_refresh=_pressure_refresh_summary(attrs),
        present_metadata_attrs=sorted(attrs),
    )


def _load_optional_json(path: Path | None, label: str) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    if path is None:
        return None, _report_base(status="not_requested", path=None)
    if not path.exists():
        return None, _report_base(
            status="not_available",
            path=str(path),
            message=f"{label} JSON was requested but does not exist",
        )
    try:
        return json.loads(path.read_text(encoding="utf-8")), _report_base(
            status="available",
            path=str(path),
        )
    except json.JSONDecodeError as exc:
        return None, _report_base(
            status="not_available",
            path=str(path),
            message=f"{label} JSON could not be decoded: {exc}",
        )


def _pipeline_risk_codes(payload: dict[str, Any]) -> list[str]:
    summary_codes = payload.get("summary", {}).get("risk_codes")
    if isinstance(summary_codes, list):
        return [str(code) for code in summary_codes]
    return [str(flag.get("code")) for flag in payload.get("risk_flags", []) if flag.get("code")]


def summarize_pipeline_json(path: Path | None) -> dict[str, Any]:
    payload, status = _load_optional_json(path, "selected-field pipeline")
    if payload is None:
        return status

    movement = payload.get("movement_geometry", {}).get("movement", {})
    spatial = payload.get("spatial_alignment_summary", {})
    evolution = payload.get("evolution_summary", {})
    return _report_base(
        **status,
        audit_status=payload.get("status"),
        audit_candidate_model_pass=payload.get("candidate_model_pass"),
        risk_codes=_pipeline_risk_codes(payload),
        movement_child_delta=movement.get("child_delta"),
        movement_large=movement.get("large_movement"),
        target_region_fraction_flags=[
            code
            for code in _pipeline_risk_codes(payload)
            if "target_region_fraction" in code
        ],
        spatial_modest_shift_improvement=spatial.get("modest_shift_improvement", []),
        spatial_different_best_shifts=spatial.get("different_best_shifts", []),
        evolution_high_normalized_rmse_variables=evolution.get(
            "high_normalized_rmse_variables", []
        ),
        message=(
            "D58 pipeline JSON is summarized for timeline attribution only; "
            "it is not a gate pass or candidate-generation input."
        ),
    )


def _source_line_numbers(text: str, pattern: str) -> list[int]:
    regex = re.compile(pattern)
    return [
        line_number
        for line_number, line in enumerate(text.splitlines(), start=1)
        if regex.search(line)
    ]


def _relative_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _source_file_inventory(path: Path, source_root: Path) -> dict[str, Any]:
    relative = _relative_path(path, source_root)
    if not path.exists():
        return _report_base(
            status="missing",
            path=relative,
            line_count=0,
            hooks={category: {} for category in SOURCE_HOOKS},
            metadata_attrs=[],
        )

    text = path.read_text(encoding="utf-8", errors="replace")
    hooks: dict[str, dict[str, list[int]]] = {}
    for category, symbols in SOURCE_HOOKS.items():
        hooks[category] = {}
        for symbol in symbols:
            if "::" in symbol:
                pattern = re.escape(symbol)
            else:
                pattern = rf"\b{re.escape(symbol)}\b"
            lines = _source_line_numbers(text, pattern)
            if lines:
                hooks[category][symbol] = lines

    return _report_base(
        status="available",
        path=relative,
        line_count=text.count("\n") + (1 if text else 0),
        hooks=hooks,
        metadata_attrs=sorted(set(re.findall(r'"(TYWRF_[A-Z0-9_]+)"', text))),
    )


def _first_hook_line(
    files: list[dict[str, Any]],
    relative_path: str,
    category: str,
    symbol: str,
) -> int | None:
    for entry in files:
        if entry.get("path") != relative_path:
            continue
        lines = entry.get("hooks", {}).get(category, {}).get(symbol, [])
        if lines:
            return int(min(lines))
    return None


def schedule_source_inventory(
    source_root: Path,
    source_files: Iterable[str],
) -> dict[str, Any]:
    files = [
        _source_file_inventory(source_root / relative_path, source_root)
        for relative_path in source_files
    ]
    available = [entry for entry in files if entry["status"] == "available"]
    all_hooks = sorted(
        {
            hook
            for entry in available
            for category in entry["hooks"].values()
            for hook, lines in category.items()
            if lines
        }
    )
    all_attrs = sorted({attr for entry in available for attr in entry["metadata_attrs"]})

    selected_path = "src/tools/selected_field_cycle.cpp"
    schedule_path = "src/dynamics/cycle_schedule.cpp"
    selected_lines = {
        "remap_overlap": _first_hook_line(
            files, selected_path, "selected_field_tool", "remap_child_state_overlap_only"
        ),
        "build_exchange_plan": _first_hook_line(
            files,
            selected_path,
            "selected_field_tool",
            "build_exposed_child_state_exchange_plan",
        ),
        "parent_interpolation": _first_hook_line(
            files,
            selected_path,
            "selected_field_tool",
            "interpolate_parent_to_exposed_child",
        ),
        "static_refresh": _first_hook_line(
            files, selected_path, "selected_field_tool", "refresh_static_fields"
        ),
        "pressure_refresh": _first_hook_line(
            files, selected_path, "selected_field_tool", "apply_pressure_refresh"
        ),
    }
    schedule_lines = {
        "parent_child_interpolation": _first_hook_line(
            files,
            schedule_path,
            "schedule",
            "CycleScheduleCallKind::parent_child_interpolation",
        ),
        "move_check": _first_hook_line(
            files,
            schedule_path,
            "schedule",
            "CycleScheduleCallKind::moving_nest_move_check",
        ),
        "post_move_parent_fill": _first_hook_line(
            files,
            schedule_path,
            "schedule",
            "CycleScheduleCallKind::moving_nest_post_move_parent_fill",
        ),
    }

    def before(lhs: int | None, rhs: int | None) -> bool | None:
        if lhs is None or rhs is None:
            return None
        return lhs < rhs

    aggregate_hooks = _report_base(
        post_move_parent_fill_hook_found=(
            "CycleScheduleCallKind::moving_nest_post_move_parent_fill" in all_hooks
            or "moving_nest_post_move_parent_fill" in all_hooks
        ),
        selected_tool_remap_before_exchange=before(
            selected_lines["remap_overlap"], selected_lines["build_exchange_plan"]
        ),
        selected_tool_exchange_before_interpolation=before(
            selected_lines["build_exchange_plan"], selected_lines["parent_interpolation"]
        ),
        selected_tool_interpolation_before_static_refresh=before(
            selected_lines["parent_interpolation"], selected_lines["static_refresh"]
        ),
        selected_tool_static_before_pressure_refresh=before(
            selected_lines["static_refresh"], selected_lines["pressure_refresh"]
        ),
        schedule_parent_interpolation_before_move_check=before(
            schedule_lines["parent_child_interpolation"], schedule_lines["move_check"]
        ),
        schedule_move_check_before_post_move_parent_fill=before(
            schedule_lines["move_check"], schedule_lines["post_move_parent_fill"]
        ),
        schedule_post_move_parent_fill_after_parent_interpolation=before(
            schedule_lines["parent_child_interpolation"],
            schedule_lines["post_move_parent_fill"],
        ),
        selected_tool_lines=selected_lines,
        cycle_schedule_lines=schedule_lines,
    )

    return _report_base(
        status="available" if available else "not_available",
        source_root=str(source_root),
        files=files,
        matched_hooks=all_hooks,
        matched_metadata_attrs=all_attrs,
        aggregate_hooks=aggregate_hooks,
        message=(
            "Source inventory is static text evidence for schedule/remap hooks; "
            "it does not prove runtime call order."
        ),
    )


def prior_audit_summaries(
    *,
    pipeline_json_path: Path | None,
    spatial_alignment_json_path: Path | None,
    evolution_json_path: Path | None,
) -> dict[str, Any]:
    return _report_base(
        pipeline=summarize_pipeline_json(pipeline_json_path),
        spatial_alignment=summarize_spatial_alignment_json(spatial_alignment_json_path),
        evolution=summarize_evolution_json(evolution_json_path),
        message=(
            "Prior JSON inputs are summarized only for diagnostic attribution and "
            "must not be used to patch, shift, tune, or select candidate fields."
        ),
    )


def _add_flag(
    flags: list[dict[str, Any]],
    code: str,
    severity: str,
    message: str,
    evidence: dict[str, Any] | None = None,
) -> None:
    flags.append(
        _report_base(
            code=code,
            severity=severity,
            message=message,
            evidence=evidence or {},
        )
    )


def _small_shift_records(spatial_summary: dict[str, Any], movement: dict[str, Any]) -> list[dict[str, Any]]:
    if not movement.get("large_movement"):
        return []
    records: list[dict[str, Any]] = []
    variables = spatial_summary.get("variables", {})
    if not isinstance(variables, dict):
        return records
    for variable, entry in variables.items():
        for block_name in ("end_state", "evolution"):
            block = (entry or {}).get(block_name, {})
            best_shift = block.get("best_shift")
            if not isinstance(best_shift, dict):
                continue
            di = _integerish_attr(best_shift.get("di"))
            dj = _integerish_attr(best_shift.get("dj"))
            if di is None or dj is None:
                continue
            if abs(di) <= SMALL_SHIFT_LIMIT and abs(dj) <= SMALL_SHIFT_LIMIT:
                records.append(
                    {
                        "variable": str(variable),
                        "block": block_name,
                        "best_shift": {"di": di, "dj": dj},
                        "baseline_normalized_rmse": block.get(
                            "baseline_normalized_rmse"
                        ),
                        "best_normalized_rmse": block.get("best_normalized_rmse"),
                    }
                )
    return records


def risk_flags(
    timeline: dict[str, Any],
    inventory: dict[str, Any],
    priors: dict[str, Any],
) -> list[dict[str, Any]]:
    flags: list[dict[str, Any]] = []
    movement = timeline.get("movement", {})
    child_delta = movement.get("child_delta")
    aggregate = inventory.get("aggregate_hooks", {})
    if movement.get("large_movement"):
        _add_flag(
            flags,
            "large_movement_interpolation_timing_ambiguity",
            "high",
            "The selected-field candidate contains a large child-grid movement; "
            "post-move parent-fill and interpolation timing require runtime proof.",
            {
                "child_delta": child_delta,
                "schedule_parent_interpolation_before_move_check": aggregate.get(
                    "schedule_parent_interpolation_before_move_check"
                ),
                "schedule_post_move_parent_fill_after_parent_interpolation": aggregate.get(
                    "schedule_post_move_parent_fill_after_parent_interpolation"
                ),
            },
        )

    if not aggregate.get("post_move_parent_fill_hook_found"):
        _add_flag(
            flags,
            "missing_explicit_post_move_parent_fill_evidence",
            "high",
            "Static source inventory did not find an explicit post-move parent-fill hook.",
        )

    missing_parent_interpolated = timeline.get(
        "selected_state_variables_without_parent_interpolation", []
    )
    if missing_parent_interpolated:
        _add_flag(
            flags,
            "selected_fields_not_all_parent_interpolated",
            "high",
            "One or more strict selected state fields are not listed as parent-interpolated.",
            {"variables": missing_parent_interpolated},
        )

    counts = timeline.get("counts", {})
    if counts.get("changed_interpolated_point_count_mismatch"):
        _add_flag(
            flags,
            "changed_interpolated_count_mismatch",
            "medium",
            "Selected-field changed-point count differs from interpolated-point count.",
            {
                "selected_field_changed_points": counts.get(
                    "selected_field_changed_points"
                ),
                "interpolated_points": counts.get("interpolated_points"),
            },
        )

    low_fraction = {
        variable: entry.get("target_region_fraction_2d")
        for variable, entry in timeline.get("exposed_overlap", {}).items()
        if entry.get("target_region_fraction_2d") is not None
        and entry["target_region_fraction_2d"] < 0.5
    }
    if low_fraction:
        _add_flag(
            flags,
            "target_region_fraction_below_half",
            "medium",
            "The exposed target region is less than half the field for one or more variables.",
            {"variables": low_fraction},
        )

    spatial = priors.get("spatial_alignment", {})
    small_shifts = _small_shift_records(spatial, movement)
    if small_shifts:
        _add_flag(
            flags,
            "small_shift_mismatch_from_d57",
            "high",
            "D57 best integer shifts are small relative to the large moving-nest delta.",
            {"small_shift_records": small_shifts, "child_delta": child_delta},
        )

    prior_high_errors = []
    evolution = priors.get("evolution", {})
    prior_high_errors.extend(evolution.get("high_normalized_rmse_variables", []) or [])
    prior_high_errors.extend(spatial.get("high_error_after_best_shift", []) or [])
    pipeline = priors.get("pipeline", {})
    pipeline_risk_codes = pipeline.get("risk_codes", []) or []
    diagnostics_fail = bool(prior_high_errors) or (
        "candidate_claims_gate_or_integrator_but_audits_fail" in pipeline_risk_codes
    )
    disposition = timeline.get("disposition", {})
    if (
        diagnostics_fail
        and (disposition.get("gate_candidate") or disposition.get("integrator_output"))
    ):
        _add_flag(
            flags,
            "candidate_claims_gate_or_integrator_but_diagnostics_fail",
            "high",
            "Candidate metadata claims gate/integrator output while prior diagnostic audits still fail.",
            {
                "gate_candidate": disposition.get("gate_candidate"),
                "integrator_output": disposition.get("integrator_output"),
                "prior_high_errors": prior_high_errors,
                "pipeline_risk_codes": pipeline_risk_codes,
            },
        )

    return flags


def audit_selected_field_timeline(
    candidate_path: Path,
    *,
    pipeline_json_path: Path | None = None,
    spatial_alignment_json_path: Path | None = None,
    evolution_json_path: Path | None = None,
    source_root: Path = PROJECT_ROOT,
    source_files: Iterable[str] = DEFAULT_SOURCE_FILES,
    variables: Iterable[str] = DEFAULT_TIMELINE_VARIABLES,
    time_index: int = -1,
) -> dict[str, Any]:
    timeline = candidate_timeline(
        candidate_path,
        variables=variables,
        time_index=time_index,
    )
    inventory = schedule_source_inventory(source_root, source_files)
    priors = prior_audit_summaries(
        pipeline_json_path=pipeline_json_path,
        spatial_alignment_json_path=spatial_alignment_json_path,
        evolution_json_path=evolution_json_path,
    )
    flags = risk_flags(timeline, inventory, priors)
    return _report_base(
        status="computed_with_flags" if flags else "computed",
        candidate=str(candidate_path),
        candidate_timeline=timeline,
        schedule_source_inventory=inventory,
        prior_audit_summaries=priors,
        risk_flags=flags,
        summary=_report_base(
            risk_flag_count=len(flags),
            risk_codes=[flag["code"] for flag in flags],
            message=(
                "Selected-field timeline audit is diagnostic-only. It reconstructs "
                "moving-nest schedule/remap evidence and never generates, patches, "
                "shifts, tunes, or selects candidate fields."
            ),
        ),
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


def report_to_json(report: dict[str, Any], *, pretty: bool = False) -> str:
    return json.dumps(_strict_json_value(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path, help="Selected-field candidate NetCDF file")
    parser.add_argument(
        "--pipeline-json",
        type=Path,
        help="Optional D58 selected-field pipeline audit JSON",
    )
    parser.add_argument(
        "--spatial-alignment-json",
        type=Path,
        help="Optional D57 selected-field spatial-alignment audit JSON",
    )
    parser.add_argument(
        "--evolution-json",
        type=Path,
        help="Optional D56 selected-field evolution audit JSON",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Repository/source root used for local schedule/remap source inventory",
    )
    parser.add_argument(
        "--source-file",
        action="append",
        default=[],
        help="Relative source file to inventory; repeat to override the default set",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_TIMELINE_VARIABLES),
        help="Candidate variables used for exposed/overlap timeline geometry",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for candidate Time-leading variable shape checks",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_selected_field_timeline(
            args.candidate,
            pipeline_json_path=args.pipeline_json,
            spatial_alignment_json_path=args.spatial_alignment_json,
            evolution_json_path=args.evolution_json,
            source_root=args.source_root,
            source_files=args.source_file or DEFAULT_SOURCE_FILES,
            variables=args.variables,
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
