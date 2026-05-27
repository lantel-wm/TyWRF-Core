#!/usr/bin/env python
"""Audit pressure-refresh seed sources without running the C++ probe."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Iterable

import netCDF4


REQUIRED_NAMES = ("P_TOP", "C3F", "C4F", "C3H", "C4H", "ALB")
BASE_STATE_RECONSTRUCTION_REQUIRED_INPUTS = (
    "HT/HGT",
    "P_TOP",
    "C3F",
    "C4F",
    "C3H",
    "C4H",
)
PHB_RECONSTRUCTION_REQUIRED_INPUTS = BASE_STATE_RECONSTRUCTION_REQUIRED_INPUTS
MASS_BASE_STATE_RECONSTRUCTION_SOURCE = "wrf_start_domain_base_state_reconstruction"
FULL_PROVIDER_RECONSTRUCTION_SOURCE = (
    "wrf_start_domain_phb_hypsometric_opt2_reconstruction"
)
FULL_LEVEL_NAMES = ("C3F", "C4F")
MASS_LEVEL_NAMES = ("C3H", "C4H")
TERRAIN_NAMES = ("HT", "HGT")
DIRECT_BASE_STATE_NAMES = ("PB", "T_INIT")
ALB_NAME = "ALB"
PHB_NAME = "PHB"
P_TOP_NAME = "P_TOP"


@dataclass(frozen=True)
class SourceEntry:
    name: str
    domain: str
    path: str
    mass_nx: int
    mass_ny: int
    mass_nz: int
    full_nz: int
    time_index: int = 0
    suitable_for_start_time_truth: bool = False


@dataclass(frozen=True)
class SourceAuditEntry:
    name: str
    domain: str
    path: str
    status: str
    missing_names: list[str]
    p_top_present: bool
    p_top_source: str | None
    alb_available: bool
    direct_alb_source: bool
    direct_phb_present: bool
    direct_phb_shape_valid: bool
    base_state_reconstruction_input_presence: dict[str, bool]
    base_state_reconstruction_inputs_complete: bool
    missing_base_state_reconstruction_inputs: list[str]
    phb_reconstruction_inputs_complete: bool
    missing_phb_reconstruction_inputs: list[str]
    base_state_terrain_source: str | None
    direct_pb_t_init_path_available: bool
    missing_direct_pb_t_init_inputs: list[str]
    base_state_reconstruction_required: bool
    recommended_next_source: str | None
    full_provider_recommended_source: str | None
    diagnostic_only: bool
    can_seed_pressure_refresh: bool
    suitable_for_start_time_truth: bool
    message: str | None = None


@dataclass(frozen=True)
class PressureRefreshSourceAudit:
    status: str
    summary: dict[str, Any]
    entries: list[SourceAuditEntry]


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y"}
    return bool(value)


def _entry_from_mapping(payload: dict[str, Any]) -> SourceEntry:
    try:
        return SourceEntry(
            name=str(payload["name"]),
            domain=str(payload["domain"]),
            path=str(payload["path"]),
            mass_nx=int(payload["mass_nx"]),
            mass_ny=int(payload["mass_ny"]),
            mass_nz=int(payload["mass_nz"]),
            full_nz=int(payload["full_nz"]),
            time_index=int(payload.get("time_index", 0)),
            suitable_for_start_time_truth=_as_bool(
                payload.get(
                    "suitable_for_start_time_truth",
                    payload.get("clean_source", False),
                )
            ),
        )
    except KeyError as exc:
        raise ValueError(f"manifest entry is missing required key: {exc.args[0]}") from exc


def entries_from_manifest(manifest: dict[str, Any]) -> list[SourceEntry]:
    raw_entries = manifest.get("entries")
    if not isinstance(raw_entries, list):
        raise ValueError("manifest must contain an entries list")
    return [_entry_from_mapping(entry) for entry in raw_entries]


def _has_global_attr(dataset: netCDF4.Dataset, name: str) -> bool:
    return name in dataset.ncattrs()


def _p_top_source(dataset: netCDF4.Dataset, time_index: int) -> str | None:
    if _has_global_attr(dataset, P_TOP_NAME):
        return "global_attribute"
    variable = dataset.variables.get(P_TOP_NAME)
    if variable is None:
        return None
    if variable.ndim == 0:
        return "scalar_variable"
    if variable.ndim == 1 and variable.shape[0] > time_index:
        return "time_variable"
    return "invalid_variable"


def _non_time_shape(variable: netCDF4.Variable, time_index: int) -> tuple[int, ...]:
    shape = tuple(int(value) for value in variable.shape)
    if variable.ndim > 0 and variable.dimensions[0] == "Time":
        if shape[0] <= time_index:
            raise ValueError(
                f"{variable.name} has Time length {shape[0]}, cannot read time_index={time_index}"
            )
        return shape[1:]
    return shape


def _check_vector_shape(
    dataset: netCDF4.Dataset,
    name: str,
    expected_length: int,
    time_index: int,
) -> str | None:
    variable = dataset.variables.get(name)
    if variable is None:
        return None
    shape = _non_time_shape(variable, time_index)
    if shape != (expected_length,):
        return f"{name} has shape {shape}, expected ({expected_length},)"
    return None


def _check_alb_shape(dataset: netCDF4.Dataset, entry: SourceEntry) -> str | None:
    variable = dataset.variables.get(ALB_NAME)
    if variable is None:
        return None
    shape = _non_time_shape(variable, entry.time_index)
    expected = (entry.mass_nz, entry.mass_ny, entry.mass_nx)
    if shape != expected:
        return f"ALB has shape {shape}, expected {expected}"
    return None


def _check_phb_shape(dataset: netCDF4.Dataset, entry: SourceEntry) -> str | None:
    variable = dataset.variables.get(PHB_NAME)
    if variable is None:
        return None
    shape = _non_time_shape(variable, entry.time_index)
    expected = (entry.full_nz, entry.mass_ny, entry.mass_nx)
    if shape != expected:
        return f"PHB has shape {shape}, expected {expected}"
    return None


def _check_terrain_source(
    dataset: netCDF4.Dataset,
    entry: SourceEntry,
) -> tuple[str | None, list[str]]:
    errors = []
    expected = (entry.mass_ny, entry.mass_nx)
    for name in TERRAIN_NAMES:
        variable = dataset.variables.get(name)
        if variable is None:
            continue
        shape = _non_time_shape(variable, entry.time_index)
        if shape == expected:
            return name, []
        errors.append(f"{name} has shape {shape}, expected {expected}")
    return None, errors


def _check_mass_3d_shape(
    dataset: netCDF4.Dataset,
    name: str,
    entry: SourceEntry,
) -> str | None:
    variable = dataset.variables.get(name)
    if variable is None:
        return None
    shape = _non_time_shape(variable, entry.time_index)
    expected = (entry.mass_nz, entry.mass_ny, entry.mass_nx)
    if shape != expected:
        return f"{name} has shape {shape}, expected {expected}"
    return None


def _default_reconstruction_presence() -> dict[str, bool]:
    return {
        "HT": False,
        "HGT": False,
        "HT/HGT": False,
        "P_TOP": False,
        "C3F": False,
        "C4F": False,
        "C3H": False,
        "C4H": False,
        "PB": False,
        "T_INIT": False,
        "direct_PB_T_INIT": False,
    }


def _unavailable_reconstruction_inputs() -> list[str]:
    return list(BASE_STATE_RECONSTRUCTION_REQUIRED_INPUTS)


def _missing_reconstruction_inputs(
    *,
    terrain_source: str | None,
    p_top_present: bool,
    vector_inputs_valid: dict[str, bool],
) -> list[str]:
    missing_inputs = []
    if terrain_source is None:
        missing_inputs.append("HT/HGT")
    if not p_top_present:
        missing_inputs.append(P_TOP_NAME)
    for name in (*FULL_LEVEL_NAMES, *MASS_LEVEL_NAMES):
        if not vector_inputs_valid[name]:
            missing_inputs.append(name)
    return missing_inputs


def audit_source_entry(entry: SourceEntry) -> SourceAuditEntry:
    path = Path(entry.path)
    if not path.exists():
        return SourceAuditEntry(
            name=entry.name,
            domain=entry.domain,
            path=str(path),
            status="nonexistent",
            missing_names=[],
            p_top_present=False,
            p_top_source=None,
            alb_available=False,
            direct_alb_source=False,
            direct_phb_present=False,
            direct_phb_shape_valid=False,
            base_state_reconstruction_input_presence=_default_reconstruction_presence(),
            base_state_reconstruction_inputs_complete=False,
            missing_base_state_reconstruction_inputs=_unavailable_reconstruction_inputs(),
            phb_reconstruction_inputs_complete=False,
            missing_phb_reconstruction_inputs=list(PHB_RECONSTRUCTION_REQUIRED_INPUTS),
            base_state_terrain_source=None,
            direct_pb_t_init_path_available=False,
            missing_direct_pb_t_init_inputs=list(DIRECT_BASE_STATE_NAMES),
            base_state_reconstruction_required=False,
            recommended_next_source=None,
            full_provider_recommended_source=None,
            diagnostic_only=True,
            can_seed_pressure_refresh=False,
            suitable_for_start_time_truth=False,
            message="path does not exist",
        )

    try:
        with netCDF4.Dataset(path) as dataset:
            p_top_source = _p_top_source(dataset, entry.time_index)
            p_top_present = p_top_source is not None and p_top_source != "invalid_variable"
            missing_names = [
                name
                for name in REQUIRED_NAMES
                if name != P_TOP_NAME and name not in dataset.variables
            ]
            if not p_top_present:
                missing_names.insert(0, P_TOP_NAME)
            alb_available = ALB_NAME in dataset.variables
            direct_phb_present = PHB_NAME in dataset.variables

            shape_errors = []
            vector_inputs_valid: dict[str, bool] = {}
            for name in FULL_LEVEL_NAMES:
                message = _check_vector_shape(dataset, name, entry.full_nz, entry.time_index)
                if message is not None:
                    shape_errors.append(message)
                vector_inputs_valid[name] = name in dataset.variables and message is None
            for name in MASS_LEVEL_NAMES:
                message = _check_vector_shape(dataset, name, entry.mass_nz, entry.time_index)
                if message is not None:
                    shape_errors.append(message)
                vector_inputs_valid[name] = name in dataset.variables and message is None
            alb_shape_error = _check_alb_shape(dataset, entry)
            if alb_shape_error is not None:
                shape_errors.append(alb_shape_error)
            phb_shape_error = _check_phb_shape(dataset, entry)
            direct_phb_shape_valid = direct_phb_present and phb_shape_error is None

            terrain_source, _terrain_shape_errors = _check_terrain_source(dataset, entry)
            direct_input_shape_errors = {
                name: _check_mass_3d_shape(dataset, name, entry)
                for name in DIRECT_BASE_STATE_NAMES
            }
            direct_input_valid = {
                name: name in dataset.variables and message is None
                for name, message in direct_input_shape_errors.items()
            }
            direct_pb_t_init_path_available = all(
                direct_input_valid[name] for name in DIRECT_BASE_STATE_NAMES
            )
            reconstruction_input_presence = {
                "HT": terrain_source == "HT",
                "HGT": terrain_source == "HGT",
                "HT/HGT": terrain_source is not None,
                "P_TOP": p_top_present,
                **vector_inputs_valid,
                **direct_input_valid,
                "direct_PB_T_INIT": direct_pb_t_init_path_available,
            }
            missing_base_state_reconstruction_inputs = _missing_reconstruction_inputs(
                terrain_source=terrain_source,
                p_top_present=p_top_present,
                vector_inputs_valid=vector_inputs_valid,
            )
            missing_phb_reconstruction_inputs = _missing_reconstruction_inputs(
                terrain_source=terrain_source,
                p_top_present=p_top_present,
                vector_inputs_valid=vector_inputs_valid,
            )
            missing_direct_pb_t_init_inputs = [
                name for name in DIRECT_BASE_STATE_NAMES if not direct_input_valid[name]
            ]
            base_state_reconstruction_inputs_complete = (
                not missing_base_state_reconstruction_inputs
            )
            phb_reconstruction_inputs_complete = not missing_phb_reconstruction_inputs
    except (OSError, RuntimeError, ValueError) as exc:
        return SourceAuditEntry(
            name=entry.name,
            domain=entry.domain,
            path=str(path),
            status="error",
            missing_names=[],
            p_top_present=False,
            p_top_source=None,
            alb_available=False,
            direct_alb_source=False,
            direct_phb_present=False,
            direct_phb_shape_valid=False,
            base_state_reconstruction_input_presence=_default_reconstruction_presence(),
            base_state_reconstruction_inputs_complete=False,
            missing_base_state_reconstruction_inputs=_unavailable_reconstruction_inputs(),
            phb_reconstruction_inputs_complete=False,
            missing_phb_reconstruction_inputs=list(PHB_RECONSTRUCTION_REQUIRED_INPUTS),
            base_state_terrain_source=None,
            direct_pb_t_init_path_available=False,
            missing_direct_pb_t_init_inputs=list(DIRECT_BASE_STATE_NAMES),
            base_state_reconstruction_required=False,
            recommended_next_source=None,
            full_provider_recommended_source=None,
            diagnostic_only=True,
            can_seed_pressure_refresh=False,
            suitable_for_start_time_truth=False,
            message=str(exc),
        )

    if shape_errors:
        return SourceAuditEntry(
            name=entry.name,
            domain=entry.domain,
            path=str(path),
            status="error",
            missing_names=missing_names,
            p_top_present=p_top_present,
            p_top_source=p_top_source,
            alb_available=alb_available,
            direct_alb_source=False,
            direct_phb_present=direct_phb_present,
            direct_phb_shape_valid=direct_phb_shape_valid,
            base_state_reconstruction_input_presence=reconstruction_input_presence,
            base_state_reconstruction_inputs_complete=False,
            missing_base_state_reconstruction_inputs=missing_base_state_reconstruction_inputs,
            phb_reconstruction_inputs_complete=False,
            missing_phb_reconstruction_inputs=missing_phb_reconstruction_inputs,
            base_state_terrain_source=terrain_source,
            direct_pb_t_init_path_available=direct_pb_t_init_path_available,
            missing_direct_pb_t_init_inputs=missing_direct_pb_t_init_inputs,
            base_state_reconstruction_required=False,
            recommended_next_source=None,
            full_provider_recommended_source=None,
            diagnostic_only=True,
            can_seed_pressure_refresh=False,
            suitable_for_start_time_truth=False,
            message="; ".join(shape_errors),
        )

    status = "ok" if not missing_names else "missing"
    can_seed = status == "ok"
    base_state_reconstruction_required = (
        not alb_available
        and base_state_reconstruction_inputs_complete
    )
    recommended_next_source = (
        MASS_BASE_STATE_RECONSTRUCTION_SOURCE
        if base_state_reconstruction_required
        else None
    )
    full_provider_recommended_source = (
        FULL_PROVIDER_RECONSTRUCTION_SOURCE
        if phb_reconstruction_inputs_complete
        else None
    )
    return SourceAuditEntry(
        name=entry.name,
        domain=entry.domain,
        path=str(path),
        status=status,
        missing_names=missing_names,
        p_top_present=p_top_present,
        p_top_source=p_top_source,
        alb_available=alb_available,
        direct_alb_source=alb_available,
        direct_phb_present=direct_phb_present,
        direct_phb_shape_valid=direct_phb_shape_valid,
        base_state_reconstruction_input_presence=reconstruction_input_presence,
        base_state_reconstruction_inputs_complete=base_state_reconstruction_inputs_complete,
        missing_base_state_reconstruction_inputs=missing_base_state_reconstruction_inputs,
        phb_reconstruction_inputs_complete=phb_reconstruction_inputs_complete,
        missing_phb_reconstruction_inputs=missing_phb_reconstruction_inputs,
        base_state_terrain_source=terrain_source,
        direct_pb_t_init_path_available=direct_pb_t_init_path_available,
        missing_direct_pb_t_init_inputs=missing_direct_pb_t_init_inputs,
        base_state_reconstruction_required=base_state_reconstruction_required,
        recommended_next_source=recommended_next_source,
        full_provider_recommended_source=full_provider_recommended_source,
        diagnostic_only=True,
        can_seed_pressure_refresh=can_seed,
        suitable_for_start_time_truth=can_seed and entry.suitable_for_start_time_truth,
    )


def audit_pressure_refresh_sources(entries: Iterable[SourceEntry]) -> PressureRefreshSourceAudit:
    audited_entries = [audit_source_entry(entry) for entry in entries]
    counts = {
        "ok_count": sum(entry.status == "ok" for entry in audited_entries),
        "missing_count": sum(entry.status == "missing" for entry in audited_entries),
        "error_count": sum(entry.status == "error" for entry in audited_entries),
        "nonexistent_count": sum(entry.status == "nonexistent" for entry in audited_entries),
    }
    d02_alb_blockers = [
        entry.name
        for entry in audited_entries
        if entry.domain == "d02" and entry.status == "missing" and ALB_NAME in entry.missing_names
    ]
    base_state_reconstruction_entries = [
        entry.name for entry in audited_entries if entry.base_state_reconstruction_required
    ]
    base_state_reconstruction_inputs_complete_entries = [
        entry.name
        for entry in audited_entries
        if entry.base_state_reconstruction_inputs_complete
    ]
    base_state_reconstruction_missing_inputs_by_entry = {
        entry.name: entry.missing_base_state_reconstruction_inputs
        for entry in audited_entries
        if entry.missing_base_state_reconstruction_inputs
    }
    direct_phb_present_entries = [
        entry.name for entry in audited_entries if entry.direct_phb_present
    ]
    direct_phb_shape_valid_entries = [
        entry.name for entry in audited_entries if entry.direct_phb_shape_valid
    ]
    direct_phb_shape_invalid_entries = [
        entry.name
        for entry in audited_entries
        if entry.direct_phb_present and not entry.direct_phb_shape_valid
    ]
    phb_reconstruction_inputs_complete_entries = [
        entry.name
        for entry in audited_entries
        if entry.phb_reconstruction_inputs_complete
    ]
    phb_reconstruction_missing_inputs_by_entry = {
        entry.name: entry.missing_phb_reconstruction_inputs
        for entry in audited_entries
        if entry.missing_phb_reconstruction_inputs
    }
    full_provider_recommended_entries = [
        entry.name
        for entry in audited_entries
        if entry.full_provider_recommended_source is not None
    ]
    direct_pb_t_init_path_entries = [
        entry.name for entry in audited_entries if entry.direct_pb_t_init_path_available
    ]
    d02_base_state_reconstruction_entries = [
        entry.name
        for entry in audited_entries
        if entry.domain == "d02" and entry.base_state_reconstruction_required
    ]
    failed = bool(
        counts["missing_count"] or counts["error_count"] or counts["nonexistent_count"]
    )
    return PressureRefreshSourceAudit(
        status="failed" if failed else "ok",
        summary={
            "entry_count": len(audited_entries),
            **counts,
            "d02_alb_blocker": bool(d02_alb_blockers),
            "d02_alb_blocker_entries": d02_alb_blockers,
            "base_state_reconstruction_required_inputs": list(
                BASE_STATE_RECONSTRUCTION_REQUIRED_INPUTS
            ),
            "base_state_reconstruction_inputs_complete_count": len(
                base_state_reconstruction_inputs_complete_entries
            ),
            "base_state_reconstruction_inputs_complete_entries": (
                base_state_reconstruction_inputs_complete_entries
            ),
            "base_state_reconstruction_missing_inputs_by_entry": (
                base_state_reconstruction_missing_inputs_by_entry
            ),
            "direct_phb_present_count": len(direct_phb_present_entries),
            "direct_phb_present_entries": direct_phb_present_entries,
            "direct_phb_shape_valid_count": len(direct_phb_shape_valid_entries),
            "direct_phb_shape_valid_entries": direct_phb_shape_valid_entries,
            "direct_phb_shape_invalid_entries": direct_phb_shape_invalid_entries,
            "phb_reconstruction_required_inputs": list(
                PHB_RECONSTRUCTION_REQUIRED_INPUTS
            ),
            "phb_reconstruction_inputs_complete_count": len(
                phb_reconstruction_inputs_complete_entries
            ),
            "phb_reconstruction_inputs_complete_entries": (
                phb_reconstruction_inputs_complete_entries
            ),
            "phb_reconstruction_missing_inputs_by_entry": (
                phb_reconstruction_missing_inputs_by_entry
            ),
            "full_provider_recommended": bool(full_provider_recommended_entries),
            "full_provider_recommended_entries": full_provider_recommended_entries,
            "full_provider_recommended_source": FULL_PROVIDER_RECONSTRUCTION_SOURCE
            if full_provider_recommended_entries
            else None,
            "direct_pb_t_init_path_available_count": len(direct_pb_t_init_path_entries),
            "direct_pb_t_init_path_available_entries": direct_pb_t_init_path_entries,
            "base_state_reconstruction_required": bool(base_state_reconstruction_entries),
            "base_state_reconstruction_entries": base_state_reconstruction_entries,
            "d02_base_state_reconstruction_required": bool(
                d02_base_state_reconstruction_entries
            ),
            "d02_base_state_reconstruction_entries": d02_base_state_reconstruction_entries,
            "recommended_next_source": MASS_BASE_STATE_RECONSTRUCTION_SOURCE
            if d02_base_state_reconstruction_entries
            else None,
        },
        entries=audited_entries,
    )


def report_to_json(report: PressureRefreshSourceAudit, *, pretty: bool = False) -> str:
    return json.dumps(asdict(report), indent=2 if pretty else None)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True, help="JSON manifest to audit")
    parser.add_argument("--output", type=Path, help="Write audit JSON to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        report = audit_pressure_refresh_sources(entries_from_manifest(manifest))
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        parser.error(str(exc))

    text = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if report.status == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
