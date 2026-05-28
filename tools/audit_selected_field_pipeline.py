#!/usr/bin/env python
"""Diagnostic-only audit of the selected-field candidate producer pipeline."""

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

from tools.audit_selected_field_state import (
    DEFAULT_STATE_VARIABLES,
    _coerce_bool,
    _coerce_int,
    _coerce_pair,
    _field_region_mask,
    _infer_target_region,
    _netcdf_attr_to_json,
    _pair_to_dict,
)


CORE_METADATA_ATTRS = (
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
    "I_PARENT_START",
    "J_PARENT_START",
    "DX",
    "DY",
)
PRESSURE_ATTR_PREFIXES = (
    "TYWRF_PRESSURE_",
    "TYWRF_PROVIDER_",
    "TYWRF_BASE_STATE_",
    "TYWRF_WOULD_SYNC_",
    "TYWRF_SYNC_",
    "TYWRF_THERMODYNAMIC_",
    "TYWRF_DRY_RUN_",
)
DEFAULT_SOURCE_FILES = (
    "src/tools/selected_field_cycle.cpp",
    "src/nest/parent_child_interpolation.cpp",
    "src/nest/state_remap.cpp",
    "src/nest/state_exchange.cpp",
    "src/nest/static_fields.cpp",
    "include/tywrf/dynamics/cycle_schedule.hpp",
    "src/dynamics/cycle_schedule.cpp",
)
SOURCE_SYMBOLS = (
    "build_candidate_state",
    "remap_child_state_overlap_only",
    "build_exposed_child_state_exchange_plan",
    "interpolate_parent_to_exposed_child",
    "refresh_static_fields",
    "apply_pressure_refresh",
    "probe_pressure_refresh_provider_readiness",
    "probe_pressure_refresh_dry_run_contract",
    "stamp_gate_metadata",
    "require_d02_resolution",
    "moving_nest_post_move_parent_fill",
    "parent_child_interpolation",
    "two_way_feedback",
)
STRICT_NORMALIZED_RMSE_THRESHOLD = 0.05
AMPLITUDE_NEAR_ONE_RANGE = (0.8, 1.25)
MODEST_SHIFT_REDUCTION_THRESHOLD = 0.25
LARGE_CHILD_DELTA_I = 50
LARGE_CHILD_DELTA_J = 30


def _read_candidate_attrs(dataset: netCDF4.Dataset) -> dict[str, Any]:
    attrs: dict[str, Any] = {}
    for name in dataset.ncattrs():
        if (
            name in CORE_METADATA_ATTRS
            or name.startswith("TYWRF_")
            or name in {"I_PARENT_START", "J_PARENT_START", "DX", "DY", "CEN_LAT", "CEN_LON"}
        ):
            attrs[name] = _netcdf_attr_to_json(dataset.getncattr(name))
    return attrs


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


def _split_csv_attr(value: Any) -> list[str]:
    if value is None:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def candidate_metadata(candidate_path: Path) -> dict[str, Any]:
    with netCDF4.Dataset(candidate_path) as dataset:
        attrs = _read_candidate_attrs(dataset)

    pressure_attrs = {
        key: value
        for key, value in attrs.items()
        if any(key.startswith(prefix) for prefix in PRESSURE_ATTR_PREFIXES)
    }
    return {
        "status": "available",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "candidate": str(candidate_path),
        "attrs": attrs,
        "present_attrs": sorted(attrs),
        "missing_core_attrs": [name for name in CORE_METADATA_ATTRS if name not in attrs],
        "pressure_refresh_attrs": pressure_attrs,
        "disposition": {
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "gate_candidate": _coerce_bool(attrs.get("TYWRF_GATE_CANDIDATE")),
            "candidate_diagnostic_only": _coerce_bool(attrs.get("TYWRF_DIAGNOSTIC_ONLY")),
            "integrator_output": _coerce_bool(attrs.get("TYWRF_INTEGRATOR_OUTPUT")),
            "validation_gate_only": _coerce_bool(attrs.get("TYWRF_VALIDATION_GATE_ONLY")),
            "candidate_kind": attrs.get("TYWRF_CANDIDATE_KIND"),
            "from_parent_start": _pair_to_dict(_coerce_pair(attrs.get("TYWRF_FROM_PARENT_START"))),
            "to_parent_start": _pair_to_dict(_coerce_pair(attrs.get("TYWRF_TO_PARENT_START"))),
            "parent_grid_ratio": _coerce_int(attrs.get("TYWRF_PARENT_GRID_RATIO")),
            "state_variables": _split_csv_attr(attrs.get("TYWRF_STATE_VARIABLES")),
            "parent_interpolated_state_variables": _split_csv_attr(
                attrs.get("TYWRF_PARENT_INTERPOLATED_STATE_VARIABLES")
            ),
            "pressure_refresh_applied": _coerce_bool(
                attrs.get("TYWRF_PRESSURE_REFRESH_APPLIED")
            ),
            "pressure_refresh_readiness_ready": _coerce_bool(
                attrs.get("TYWRF_PRESSURE_REFRESH_READINESS_READY")
            ),
            "pressure_refresh_integration_status": attrs.get(
                "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS"
            ),
        },
    }


def _candidate_variable_shape(
    dataset: netCDF4.Dataset,
    variable: str,
    *,
    time_index: int,
) -> tuple[int, ...] | None:
    if variable not in dataset.variables:
        return None
    del time_index
    netcdf_variable = dataset.variables[variable]
    shape = netcdf_variable.shape
    dims = netcdf_variable.dimensions
    if dims and dims[0] == "Time":
        shape = shape[1:]
    return tuple(int(value) for value in shape)


def movement_geometry(
    candidate_path: Path,
    *,
    variables: Iterable[str],
    time_index: int,
) -> dict[str, Any]:
    with netCDF4.Dataset(candidate_path) as dataset:
        region = {
            **_infer_target_region(dataset),
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
        }
        mass_shape = _candidate_variable_shape(dataset, "MU", time_index=time_index)
        if mass_shape is not None and len(mass_shape) >= 2:
            mass_horizontal_shape = (int(mass_shape[-2]), int(mass_shape[-1]))
        else:
            mass_horizontal_shape = None

        variable_regions: dict[str, Any] = {}
        for variable in variables:
            shape = _candidate_variable_shape(dataset, variable, time_index=time_index)
            if shape is None or len(shape) < 2:
                variable_regions[variable] = {
                    "status": "not_available",
                    "diagnostic_only": True,
                    "candidate_model_pass": "not_applicable",
                    "message": "candidate variable missing or has no horizontal shape",
                }
                continue
            _, mask_info = _field_region_mask(
                variable,
                (int(shape[-2]), int(shape[-1])),
                mass_horizontal_shape,
                region,
            )
            target = _numeric_attr(mask_info.get("target_point_count_2d"))
            overlap = _numeric_attr(mask_info.get("overlap_point_count_2d"))
            if target is not None and overlap is not None and target + overlap > 0:
                target_fraction = target / (target + overlap)
            else:
                target_fraction = None
            variable_regions[variable] = {
                **mask_info,
                "diagnostic_only": True,
                "candidate_model_pass": "not_applicable",
                "candidate_shape": list(shape),
                "target_region_fraction_2d": target_fraction,
            }

    child_delta = region.get("child_delta") or {}
    delta_i = _coerce_int(child_delta.get("i"))
    delta_j = _coerce_int(child_delta.get("j"))
    if delta_i is not None and delta_j is not None:
        movement = {
            "status": "available",
            "child_delta": {"i": delta_i, "j": delta_j},
            "manhattan_child_cells": abs(delta_i) + abs(delta_j),
            "euclidean_child_cells": math.hypot(delta_i, delta_j),
            "large_movement": abs(delta_i) >= LARGE_CHILD_DELTA_I
            or abs(delta_j) >= LARGE_CHILD_DELTA_J,
        }
    else:
        movement = {
            "status": "not_available",
            "child_delta": None,
            "large_movement": None,
        }

    return {
        "status": "available" if region.get("status") == "available" else "not_available",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "region": region,
        "movement": movement,
        "mass_horizontal_shape": (
            None
            if mass_horizontal_shape is None
            else {"ny": mass_horizontal_shape[0], "nx": mass_horizontal_shape[1]}
        ),
        "variable_regions": variable_regions,
    }


def _load_optional_json(path: Path | None, label: str) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    if path is None:
        return None, {
            "status": "not_requested",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "path": None,
        }
    if not path.exists():
        return None, {
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "path": str(path),
            "message": f"{label} JSON was requested but does not exist",
        }
    try:
        return json.loads(path.read_text(encoding="utf-8")), {
            "status": "available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "path": str(path),
        }
    except json.JSONDecodeError as exc:
        return None, {
            "status": "not_available",
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "path": str(path),
            "message": f"{label} JSON could not be decoded: {exc}",
        }


def summarize_evolution_json(path: Path | None) -> dict[str, Any]:
    payload, status = _load_optional_json(path, "evolution")
    if payload is None:
        return status

    variables: dict[str, Any] = {}
    high_error: list[str] = []
    amplitude_near_one_with_high_error: list[str] = []
    target_error_fraction_below_half: list[str] = []
    for entry in payload.get("variables", []):
        variable = str(entry.get("variable"))
        amplitude = entry.get("evolution_amplitude") or {}
        region_split = entry.get("region_split") or {}
        normalized = _numeric_attr(entry.get("normalized_rmse"))
        ratio = _numeric_attr(amplitude.get("amplitude_ratio"))
        capture = _numeric_attr(amplitude.get("capture_fraction"))
        target_error_fraction = _numeric_attr(region_split.get("target_error_fraction"))
        if normalized is not None and normalized > STRICT_NORMALIZED_RMSE_THRESHOLD:
            high_error.append(variable)
        if (
            normalized is not None
            and normalized > STRICT_NORMALIZED_RMSE_THRESHOLD
            and ratio is not None
            and AMPLITUDE_NEAR_ONE_RANGE[0] <= ratio <= AMPLITUDE_NEAR_ONE_RANGE[1]
        ):
            amplitude_near_one_with_high_error.append(variable)
        if target_error_fraction is not None and target_error_fraction < 0.5:
            target_error_fraction_below_half.append(variable)
        variables[variable] = {
            "status": entry.get("status"),
            "normalized_rmse": normalized,
            "rmse": _numeric_attr(entry.get("rmse")),
            "amplitude_ratio": ratio,
            "capture_fraction": capture,
            "target_error_fraction": target_error_fraction,
            "target_region_dominates_global_error": region_split.get(
                "target_region_dominates_global_error"
            ),
        }

    return {
        **status,
        "audit_status": payload.get("status"),
        "audit_candidate_model_pass": payload.get("candidate_model_pass"),
        "summary": payload.get("summary", {}),
        "variables": variables,
        "high_normalized_rmse_variables": high_error,
        "amplitude_near_one_with_high_error": amplitude_near_one_with_high_error,
        "target_error_fraction_below_half": target_error_fraction_below_half,
        "message": (
            "D56 evolution JSON is summarized for diagnostic attribution only; "
            "it is not a validation gate pass."
        ),
    }


def _shift_tuple(value: Any) -> tuple[int | None, int | None] | None:
    if not isinstance(value, dict):
        return None
    return _coerce_int(value.get("di")), _coerce_int(value.get("dj"))


def _alignment_block_summary(block: dict[str, Any]) -> dict[str, Any]:
    global_block = block.get("global") or {}
    baseline = global_block.get("baseline") or {}
    best = global_block.get("best") or {}
    return {
        "status": block.get("status"),
        "best_shift": global_block.get("best_shift"),
        "improved": global_block.get("improved"),
        "baseline_normalized_rmse": _numeric_attr(baseline.get("normalized_rmse")),
        "best_normalized_rmse": _numeric_attr(best.get("normalized_rmse")),
        "normalized_rmse_reduction_fraction": _numeric_attr(
            global_block.get("normalized_rmse_reduction_fraction")
        ),
    }


def summarize_spatial_alignment_json(path: Path | None) -> dict[str, Any]:
    payload, status = _load_optional_json(path, "spatial alignment")
    if payload is None:
        return status

    variables: dict[str, Any] = {}
    modest_improvement: list[str] = []
    different_best_shifts: list[str] = []
    high_error_after_best_shift: list[str] = []
    for entry in payload.get("variables", []):
        variable = str(entry.get("variable"))
        end_state = _alignment_block_summary(entry.get("end_state") or {})
        evolution = _alignment_block_summary(entry.get("evolution") or {})
        end_shift = _shift_tuple(end_state.get("best_shift"))
        evolution_shift = _shift_tuple(evolution.get("best_shift"))
        for label, block in (("end_state", end_state), ("evolution", evolution)):
            reduction = block.get("normalized_rmse_reduction_fraction")
            baseline = block.get("baseline_normalized_rmse")
            if (
                baseline is not None
                and baseline > STRICT_NORMALIZED_RMSE_THRESHOLD
                and reduction is not None
                and reduction <= MODEST_SHIFT_REDUCTION_THRESHOLD
            ):
                modest_improvement.append(f"{variable}:{label}")
            best = block.get("best_normalized_rmse")
            if best is not None and best > STRICT_NORMALIZED_RMSE_THRESHOLD:
                high_error_after_best_shift.append(f"{variable}:{label}")
        if end_shift is not None and evolution_shift is not None and end_shift != evolution_shift:
            different_best_shifts.append(variable)
        variables[variable] = {
            "status": entry.get("status"),
            "end_state": end_state,
            "evolution": evolution,
        }

    return {
        **status,
        "audit_status": payload.get("status"),
        "audit_candidate_model_pass": payload.get("candidate_model_pass"),
        "summary": payload.get("summary", {}),
        "variables": variables,
        "modest_shift_improvement": modest_improvement,
        "different_best_shifts": different_best_shifts,
        "high_error_after_best_shift": high_error_after_best_shift,
        "message": (
            "D57 spatial-alignment JSON is summarized for diagnostic attribution only; "
            "best shifts must not be applied to candidate fields."
        ),
    }


def _source_line_numbers(text: str, pattern: str) -> list[int]:
    regex = re.compile(pattern)
    return [index for index, line in enumerate(text.splitlines(), start=1) if regex.search(line)]


def _source_file_inventory(path: Path, source_root: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "status": "missing",
            "path": str(path.relative_to(source_root) if path.is_relative_to(source_root) else path),
            "matched_symbols": {},
            "matched_attrs": [],
            "matched_phrases": {},
        }
    text = path.read_text(encoding="utf-8", errors="replace")
    matched_symbols = {
        symbol: _source_line_numbers(text, rf"\b{re.escape(symbol)}\b")
        for symbol in SOURCE_SYMBOLS
        if re.search(rf"\b{re.escape(symbol)}\b", text)
    }
    matched_attrs = sorted(set(re.findall(r'"(TYWRF_[A-Z0-9_]+)"', text)))
    matched_phrases = {
        "rejects_reference_end_oracle": _source_line_numbers(
            text, r"rejects end-state/reference-end|Oracle inputs are not accepted"
        ),
        "start_state_only_candidate": _source_line_numbers(
            text, r"from start states only|d01 start-state|d02 start-state"
        ),
        "overlap_then_exposed_parent_fill": _source_line_numbers(
            text, r"overlap_only|exposed.*parent|parent.*exposed"
        ),
        "schedule_or_timing_hook": _source_line_numbers(
            text, r"moving_nest_post_move_parent_fill|snap_to_next_parent_step|parent_child_interpolation"
        ),
    }
    matched_phrases = {key: value for key, value in matched_phrases.items() if value}
    return {
        "status": "available",
        "path": str(path.relative_to(source_root) if path.is_relative_to(source_root) else path),
        "line_count": text.count("\n") + (1 if text else 0),
        "matched_symbols": matched_symbols,
        "matched_attrs": matched_attrs,
        "matched_phrases": matched_phrases,
    }


def source_inventory(source_root: Path, source_files: Iterable[str]) -> dict[str, Any]:
    files = [
        _source_file_inventory(source_root / relative_path, source_root)
        for relative_path in source_files
    ]
    available = [entry for entry in files if entry["status"] == "available"]
    all_symbols = sorted(
        {
            symbol
            for entry in available
            for symbol, lines in entry["matched_symbols"].items()
            if lines
        }
    )
    all_attrs = sorted({attr for entry in available for attr in entry["matched_attrs"]})
    phrases = {
        key
        for entry in available
        for key, lines in entry["matched_phrases"].items()
        if lines
    }
    assumptions = {
        "rejects_reference_end_oracle_inputs": "rejects_reference_end_oracle" in phrases,
        "starts_from_d02_start_state_overlap_remap": "remap_child_state_overlap_only" in all_symbols,
        "exposed_selected_fields_parent_interpolated": (
            "build_exposed_child_state_exchange_plan" in all_symbols
            and "interpolate_parent_to_exposed_child" in all_symbols
        ),
        "static_fields_refreshed_without_reference_end": "refresh_static_fields" in all_symbols,
        "pressure_refresh_is_optional_hook": (
            "probe_pressure_refresh_provider_readiness" in all_symbols
            or "apply_pressure_refresh" in all_symbols
        ),
        "schedule_contains_post_move_parent_fill_hook": (
            "moving_nest_post_move_parent_fill" in all_symbols
        ),
    }
    return {
        "status": "available" if available else "not_available",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "source_root": str(source_root),
        "files": files,
        "matched_symbols": all_symbols,
        "matched_attrs": all_attrs,
        "producer_path_assumptions": assumptions,
        "message": (
            "Source inventory is text-based and identifies likely producer, remap, "
            "pressure-refresh, and timing hooks; it does not prove runtime call order."
        ),
    }


def _add_flag(
    flags: list[dict[str, Any]],
    code: str,
    severity: str,
    message: str,
    evidence: dict[str, Any] | None = None,
) -> None:
    flags.append(
        {
            "code": code,
            "severity": severity,
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "message": message,
            "evidence": evidence or {},
        }
    )


def risk_flags(
    metadata: dict[str, Any],
    geometry: dict[str, Any],
    evolution: dict[str, Any],
    spatial: dict[str, Any],
    inventory: dict[str, Any],
) -> list[dict[str, Any]]:
    flags: list[dict[str, Any]] = []

    low_fraction = {
        variable: info.get("target_region_fraction_2d")
        for variable, info in geometry.get("variable_regions", {}).items()
        if info.get("target_region_fraction_2d") is not None
        and info["target_region_fraction_2d"] < 0.5
    }
    if low_fraction:
        _add_flag(
            flags,
            "target_region_fraction_lt_0_5",
            "medium",
            "The moved/exposed target region is less than half of the field for one or more variables.",
            {"variables": low_fraction},
        )

    target_error_low = evolution.get("target_error_fraction_below_half") or []
    if target_error_low:
        _add_flag(
            flags,
            "target_error_fraction_lt_0_5",
            "medium",
            "Evolution audit reports less than half of global error in the target region.",
            {"variables": target_error_low},
        )

    near_one = evolution.get("amplitude_near_one_with_high_error") or []
    if near_one:
        _add_flag(
            flags,
            "amplitude_near_one_but_rmse_high",
            "high",
            "Evolution amplitude is near WRF while strict normalized RMSE remains high.",
            {"variables": near_one},
        )

    modest = spatial.get("modest_shift_improvement") or []
    if modest:
        _add_flag(
            flags,
            "best_shift_modest_improvement",
            "high",
            "Best integer shifts improve normalized RMSE only modestly.",
            {"variable_blocks": modest},
        )

    differing = spatial.get("different_best_shifts") or []
    if differing:
        _add_flag(
            flags,
            "different_best_shifts_end_vs_evolution",
            "high",
            "End-state and evolution alignment prefer different integer shifts.",
            {"variables": differing},
        )

    disposition = metadata.get("disposition", {})
    audit_high_errors = (evolution.get("high_normalized_rmse_variables") or []) or (
        spatial.get("high_error_after_best_shift") or []
    )
    if (
        (disposition.get("gate_candidate") or disposition.get("integrator_output"))
        and audit_high_errors
    ):
        _add_flag(
            flags,
            "candidate_claims_gate_or_integrator_but_audits_fail",
            "high",
            "Candidate metadata claims gate/integrator output while diagnostic audits still show strict-field errors.",
            {
                "gate_candidate": disposition.get("gate_candidate"),
                "integrator_output": disposition.get("integrator_output"),
                "audit_high_errors": audit_high_errors,
            },
        )

    movement = geometry.get("movement", {})
    child_delta = movement.get("child_delta") or {}
    if movement.get("large_movement"):
        _add_flag(
            flags,
            "large_movement_delta_schedule_remap_risk",
            "high",
            "Movement delta is large; schedule/remap timing may dominate selected-field errors.",
            {"child_delta": child_delta},
        )

    assumptions = inventory.get("producer_path_assumptions", {})
    if not assumptions.get("schedule_contains_post_move_parent_fill_hook"):
        _add_flag(
            flags,
            "schedule_post_move_parent_fill_hook_not_found",
            "medium",
            "Source inventory did not find a post-move parent-fill schedule hook.",
        )
    if not assumptions.get("exposed_selected_fields_parent_interpolated"):
        _add_flag(
            flags,
            "selected_field_parent_interpolation_path_incomplete",
            "medium",
            "Source inventory did not find the full exposed-cell parent interpolation path.",
        )

    return flags


def audit_selected_field_pipeline(
    candidate_path: Path,
    *,
    evolution_json_path: Path | None = None,
    spatial_alignment_json_path: Path | None = None,
    source_root: Path = PROJECT_ROOT,
    source_files: Iterable[str] = DEFAULT_SOURCE_FILES,
    variables: Iterable[str] = DEFAULT_STATE_VARIABLES,
    time_index: int = -1,
) -> dict[str, Any]:
    variables = tuple(variables)
    metadata = candidate_metadata(candidate_path)
    geometry = movement_geometry(candidate_path, variables=variables, time_index=time_index)
    evolution = summarize_evolution_json(evolution_json_path)
    spatial = summarize_spatial_alignment_json(spatial_alignment_json_path)
    inventory = source_inventory(source_root, source_files)
    flags = risk_flags(metadata, geometry, evolution, spatial, inventory)
    return {
        "status": "computed_with_flags" if flags else "computed",
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "candidate": str(candidate_path),
        "candidate_metadata": metadata,
        "movement_geometry": geometry,
        "evolution_summary": evolution,
        "spatial_alignment_summary": spatial,
        "source_inventory": inventory,
        "risk_flags": flags,
        "summary": {
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
            "risk_flag_count": len(flags),
            "risk_codes": [flag["code"] for flag in flags],
            "message": (
                "Selected-field pipeline audit is diagnostic-only. Reference-end-derived "
                "JSON inputs are summarized only for error attribution and must not be "
                "used to generate, patch, shift, tune, or select candidate fields."
            ),
        },
    }


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
        "--evolution-json",
        type=Path,
        help="Optional D56 selected-field evolution audit JSON",
    )
    parser.add_argument(
        "--spatial-alignment-json",
        type=Path,
        help="Optional D57 selected-field spatial-alignment audit JSON",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Repository/source root used for local C++ text inventory",
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
        default=list(DEFAULT_STATE_VARIABLES),
        help="Candidate variables used for movement/exposed-region geometry",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for candidate Time-leading variables",
    )
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_selected_field_pipeline(
            args.candidate,
            evolution_json_path=args.evolution_json,
            spatial_alignment_json_path=args.spatial_alignment_json,
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
