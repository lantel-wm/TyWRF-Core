#!/usr/bin/env python
"""Audit diagnostic pressure-refresh producer/staging evidence.

This tool is diagnostic-only. It never edits candidate NetCDF files, never
generates candidate fields, and never reports a validation-gate pass.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable

import netCDF4
import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.audit_pressure_refresh_candidate import (  # noqa: E402
    _candidate_metadata,
    _coerce_bool,
    _coerce_int,
    _netcdf_attr_to_json,
)


SOURCE_RELATIVE_PATHS = (
    "src/dynamics/pressure_refresh.cpp",
    "src/dynamics/pressure_refresh_hook.cpp",
    "src/dynamics/pressure_refresh_staging.cpp",
    "include/tywrf/dynamics/pressure_refresh.hpp",
    "include/tywrf/dynamics/pressure_refresh_hook.hpp",
    "include/tywrf/dynamics/pressure_refresh_staging.hpp",
)
CORE_METADATA_ATTRS = (
    "TYWRF_DIAGNOSTIC_ONLY",
    "TYWRF_GATE_CANDIDATE",
    "TYWRF_INTEGRATOR_OUTPUT",
    "TYWRF_VALIDATION_GATE_ONLY",
    "TYWRF_CANDIDATE_KIND",
    "TYWRF_CANDIDATE_MESSAGE",
    "TYWRF_FROM_PARENT_START",
    "TYWRF_TO_PARENT_START",
    "TYWRF_PARENT_GRID_RATIO",
    "TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT",
)
PRESSURE_METADATA_PREFIXES = (
    "TYWRF_PRESSURE_REFRESH_",
    "TYWRF_PRESSURE_COMPUTE_",
    "TYWRF_DRY_RUN_",
    "TYWRF_WOULD_",
    "TYWRF_SYNC_",
)
PRESSURE_HELPER_ATTRS = (
    "TYWRF_PRESSURE_REFRESH_HELPER_NAME",
    "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE",
    "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE",
    "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE",
    "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_SOURCE",
    "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_PROVENANCE",
)
COMPANION_FIELDS = ("PB", "P+PB", "MU", "MUB", "PHB")
WEAK_COMPANION_ABS_TOLERANCE = 100.0
SOURCE_DELTA_ABS_TOLERANCE = 1.0e-6
SOURCE_DELTA_REL_TOLERANCE = 1.0e-5
P_BIAS_ABS_TOLERANCE = 500.0


@dataclass(frozen=True)
class RiskFlag:
    code: str
    severity: str
    message: str
    evidence: dict[str, Any]
    diagnostic_only: bool = True
    candidate_model_pass: str = "not_applicable"


@dataclass(frozen=True)
class ProducerAudit:
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    metadata_contract: dict[str, Any]
    code_producer_inventory: dict[str, Any]
    formula_risk_flags: list[RiskFlag]
    vertical_bias_summary: dict[str, Any]


def _read_relevant_attrs(dataset: netCDF4.Dataset) -> dict[str, Any]:
    attrs: dict[str, Any] = {}
    for name in dataset.ncattrs():
        if name in CORE_METADATA_ATTRS or name.startswith(PRESSURE_METADATA_PREFIXES):
            attrs[name] = _netcdf_attr_to_json(dataset.getncattr(name))
    return attrs


def _present_subset(attrs: dict[str, Any], names: Iterable[str]) -> dict[str, Any]:
    return {name: attrs[name] for name in names if name in attrs}


def _metadata_contract(candidate_path: Path | None) -> dict[str, Any]:
    if candidate_path is None:
        return {
            "status": "not_available",
            "candidate": None,
            "attrs": {},
            "present_attrs": [],
            "helper_source_attrs": {},
            "counts": {},
            "disposition": {
                "diagnostic_only": True,
                "candidate_model_pass": "not_applicable",
                "message": "candidate NetCDF was not supplied",
            },
        }

    with netCDF4.Dataset(candidate_path) as dataset:
        attrs = _read_relevant_attrs(dataset)
        candidate_metadata = _candidate_metadata(dataset)

    counts = {
        "target_column_count": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT")
        ),
        "refreshed_column_count": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT")
        ),
        "refreshed_point_count": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT")
        ),
        "refreshed_p_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS")
        ),
        "changed_p_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS")
        ),
        "changed_pb_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS")
        ),
        "changed_mub_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS")
        ),
        "changed_phb_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS")
        ),
        "synced_pb_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS")
        ),
        "synced_mub_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS")
        ),
        "synced_phb_points": _coerce_int(
            attrs.get("TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS")
        ),
        "would_refresh_p_point_count": _coerce_int(
            attrs.get("TYWRF_WOULD_REFRESH_P_POINT_COUNT")
        ),
    }
    bools = {
        "applied": _coerce_bool(attrs.get("TYWRF_PRESSURE_REFRESH_APPLIED")),
        "compute_called": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED")
        ),
        "provider_ok": _coerce_bool(attrs.get("TYWRF_PRESSURE_REFRESH_PROVIDER_OK")),
        "staging_ok": _coerce_bool(attrs.get("TYWRF_PRESSURE_REFRESH_STAGING_OK")),
        "terrain_override_used": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_TERRAIN_OVERRIDE_USED")
        ),
        "provider_terrain_uses_moved_candidate_hgt": _coerce_bool(
            attrs.get("TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT")
        ),
        "changed_p_matches_refreshed_point_count": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT")
        ),
        "invalid_and_skipped_points_zero": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO")
        ),
        "overlap_halo_untouched": _coerce_bool(
            attrs.get("TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED")
        ),
    }
    helper_source_attrs = _present_subset(attrs, PRESSURE_HELPER_ATTRS)

    return {
        "status": "available",
        "candidate": str(candidate_path),
        "attrs": attrs,
        "present_attrs": sorted(attrs),
        "helper_source_attrs": helper_source_attrs,
        "helper_source_attrs_present": bool(helper_source_attrs),
        "counts": counts,
        "flags": bools,
        "disposition": candidate_metadata.get("disposition", {}),
    }


_FUNCTION_RE = re.compile(
    r"(?P<prefix>(?:\[\[nodiscard\]\]\s*)?(?:constexpr\s+)?(?:template\s*<[^>]+>\s*)?"
    r"(?:[A-Za-z_][\w:<>,\s*&]+\s+)+)"
    r"(?P<name>[A-Za-z_][\w:]*)\s*\([^;{}]*\)\s*(?:noexcept\s*)?(?:const\s*)?\{",
    re.MULTILINE,
)


def _line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def _brace_end(text: str, open_brace_index: int) -> int:
    depth = 0
    for index in range(open_brace_index, len(text)):
        char = text[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return index + 1
    return len(text)


def _function_records(path: Path, relative_path: str, text: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for match in _FUNCTION_RE.finditer(text):
        name = match.group("name")
        open_brace = text.find("{", match.end() - 1)
        if open_brace < 0:
            continue
        body = text[open_brace : _brace_end(text, open_brace)]
        signals = {
            "writes_p": bool(
                re.search(r"\binputs\.p\s*\(", body)
                and re.search(r"\binputs\.p\s*\([^;]+=", body)
            )
            or bool(re.search(r"\bstate_view\.p\s*=", body))
            or bool(re.search(r"\bscratch_view\.p\s*=", body)),
            "computes_perturbation_p": bool(
                "perturbation_pressure_pa" in body
                or "const auto perturbation = total_pressure - pb" in body
            ),
            "calls_pressure_refresh_compute": (
                "refresh_krosa_moving_nest_pressure" in body
                and name != "refresh_krosa_moving_nest_pressure"
            ),
            "stages_pressure_inputs": "PressureRefreshInputs" in body
            or "make_krosa_pressure_refresh_inputs" in name,
            "syncs_base_state_companions": all(
                token in body for token in ("pb", "mub", "phb", "sync_exposed")
            ),
            "uses_moved_candidate_hgt": "moved_candidate_HGT" in body,
        }
        formula_terms = [
            term
            for term in (
                "MU+MUB",
                "PH+PHB",
                "PB subtraction",
                "ALB",
                "C3F/C4F",
                "C3H/C4H",
                "P_TOP",
                "theta_m",
            )
            if _term_present(term, body)
        ]
        if any(signals.values()) or "pressure" in name.lower():
            records.append(
                {
                    "name": name,
                    "line": _line_number(text, match.start()),
                    "signals": signals,
                    "formula_terms": formula_terms,
                    "path": str(path),
                    "relative_path": relative_path,
                }
            )
    return records


def _term_present(term: str, body: str) -> bool:
    if term == "MU+MUB":
        return "inputs.mu" in body and "inputs.mub" in body
    if term == "PH+PHB":
        return "inputs.ph" in body and "inputs.phb" in body
    if term == "PB subtraction":
        return "total_pressure - pb" in body
    if term == "ALB":
        return "inputs.alb" in body or "external_alb" in body
    if term == "C3F/C4F":
        return "c3f" in body and "c4f" in body
    if term == "C3H/C4H":
        return "c3h" in body and "c4h" in body
    if term == "P_TOP":
        return "p_top" in body or "p_top_pa" in body
    if term == "theta_m":
        return "base_potential_temperature_k" in body or "use_theta_m" in body
    return False


def _source_inventory(source_root: Path) -> dict[str, Any]:
    files: list[dict[str, Any]] = []
    all_functions: list[dict[str, Any]] = []
    for relative_path in SOURCE_RELATIVE_PATHS:
        path = source_root / relative_path
        if not path.exists():
            files.append(
                {
                    "relative_path": relative_path,
                    "path": str(path),
                    "status": "missing",
                    "functions": [],
                }
            )
            continue
        text = path.read_text(encoding="utf-8")
        functions = _function_records(path, relative_path, text)
        all_functions.extend(functions)
        files.append(
            {
                "relative_path": relative_path,
                "path": str(path),
                "status": "scanned",
                "functions": functions,
            }
        )

    def matching(signal: str) -> list[dict[str, Any]]:
        return [
            {
                "name": item["name"],
                "relative_path": item["relative_path"],
                "line": item["line"],
            }
            for item in all_functions
            if item["signals"].get(signal)
        ]

    return {
        "status": "scanned" if any(item["status"] == "scanned" for item in files) else "not_available",
        "source_root": str(source_root),
        "files": files,
        "summary": {
            "source_files_scanned": sum(1 for item in files if item["status"] == "scanned"),
            "source_files_missing": sum(1 for item in files if item["status"] == "missing"),
            "p_writer_functions": matching("writes_p"),
            "p_compute_functions": matching("computes_perturbation_p"),
            "refresh_callers": matching("calls_pressure_refresh_compute"),
            "staging_functions": matching("stages_pressure_inputs"),
            "base_state_sync_functions": matching("syncs_base_state_companions"),
            "moved_hgt_mentions": matching("uses_moved_candidate_hgt"),
        },
    }


def _load_vertical_bias(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("vertical-bias JSON must contain an object")
    return payload


def _finite_float(value: Any) -> float | None:
    if isinstance(value, (int, float, np.integer, np.floating)):
        numeric = float(value)
        if math.isfinite(numeric):
            return numeric
    return None


def _close(a: float, b: float) -> bool:
    return abs(a - b) <= max(SOURCE_DELTA_ABS_TOLERANCE, SOURCE_DELTA_REL_TOLERANCE * max(abs(a), abs(b), 1.0))


def _companion_map(attribution: dict[str, Any]) -> dict[str, dict[str, Any]]:
    raw = attribution.get("companion_fields", [])
    if not isinstance(raw, list):
        return {}
    return {
        str(item.get("field")): item
        for item in raw
        if isinstance(item, dict) and item.get("field") is not None
    }


def _vertical_bias_summary(payload: dict[str, Any] | None) -> dict[str, Any]:
    if payload is None:
        return {
            "status": "not_available",
            "source": None,
            "message": "vertical-bias JSON was not supplied",
        }

    attribution = payload.get("worst_level_attribution", [])
    if not isinstance(attribution, list):
        attribution = []
    levels: list[dict[str, Any]] = []
    for item in attribution:
        if not isinstance(item, dict):
            continue
        pressure = item.get("pressure_metrics", {})
        if not isinstance(pressure, dict):
            pressure = {}
        companions = _companion_map(item)
        companion_diffs = {
            field: _finite_float(companions.get(field, {}).get("mean_diff"))
            for field in COMPANION_FIELDS
            if field in companions
        }
        evolution = item.get("source_evolution", {})
        if not isinstance(evolution, dict):
            evolution = {}
        levels.append(
            {
                "level": item.get("level", pressure.get("level")),
                "rank_by_rmse": item.get("rank_by_rmse"),
                "rank_by_abs_bias": item.get("rank_by_abs_bias"),
                "p_mean_diff": _finite_float(pressure.get("mean_diff")),
                "p_rmse": _finite_float(pressure.get("rmse")),
                "p_normalized_rmse": _finite_float(pressure.get("normalized_rmse")),
                "companion_mean_diffs": companion_diffs,
                "source_evolution_delta": _finite_float(
                    evolution.get("mean_candidate_evolution_minus_reference_evolution_p")
                ),
                "source_evolution_status": evolution.get("status"),
            }
        )

    p_biases = [
        abs(item["p_mean_diff"])
        for item in levels
        if _finite_float(item.get("p_mean_diff")) is not None
    ]
    return {
        "status": "available",
        "source": payload.get("candidate_end") or payload.get("candidate"),
        "vertical_bias_status": payload.get("status"),
        "summary": payload.get("summary", {}),
        "worst_levels_by_rmse": payload.get("worst_levels_by_rmse", []),
        "worst_levels_by_abs_bias": payload.get("worst_levels_by_abs_bias", []),
        "levels": levels,
        "max_abs_p_mean_diff": max(p_biases, default=None),
    }


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


def _metadata_risk_flags(metadata: dict[str, Any], flags: list[RiskFlag]) -> None:
    if metadata.get("status") != "available":
        return
    helper_attrs = metadata.get("helper_source_attrs") or {}
    if helper_attrs:
        _add_flag(
            flags,
            "pressure_refresh_helper_metadata_present",
            "info",
            "candidate metadata names an explicit pressure-refresh helper/source path",
            {"helper_source_attrs": helper_attrs},
        )

    typed_flags = metadata.get("flags") or {}
    counts = metadata.get("counts") or {}
    if typed_flags.get("changed_p_matches_refreshed_point_count") is True:
        _add_flag(
            flags,
            "changed_p_points_match_refreshed_count",
            "info",
            "changed P point count matches pressure-refresh refreshed point count",
            {
                "changed_p_points": counts.get("changed_p_points"),
                "refreshed_point_count": counts.get("refreshed_point_count"),
                "refreshed_p_points": counts.get("refreshed_p_points"),
            },
        )

    attrs = metadata.get("attrs") or {}
    terrain_values = [
        value
        for name, value in attrs.items()
        if "TERRAIN" in name or "HGT" in name or name in PRESSURE_HELPER_ATTRS
    ]
    if typed_flags.get("provider_terrain_uses_moved_candidate_hgt") is True or any(
        "moved_candidate_HGT" in str(value) for value in terrain_values
    ):
        _add_flag(
            flags,
            "terrain_provenance_uses_moved_candidate_hgt",
            "warning",
            "pressure-refresh terrain provenance uses moved candidate HGT",
            {
                "provider_terrain_uses_moved_candidate_hgt": typed_flags.get(
                    "provider_terrain_uses_moved_candidate_hgt"
                ),
                "terrain_values": terrain_values,
            },
        )


def _code_risk_flags(inventory: dict[str, Any], flags: list[RiskFlag]) -> None:
    summary = inventory.get("summary") or {}
    writers = summary.get("p_writer_functions") or []
    computes = summary.get("p_compute_functions") or []
    syncs = summary.get("base_state_sync_functions") or []
    if writers and computes:
        _add_flag(
            flags,
            "source_inventory_p_written_by_pressure_refresh",
            "info",
            "source inventory found pressure-refresh code that computes and writes perturbation P",
            {"p_writer_functions": writers, "p_compute_functions": computes},
        )
    if syncs and writers:
        _add_flag(
            flags,
            "source_inventory_companion_sync_then_p_write",
            "info",
            "source inventory found companion-field sync before pressure-refresh P write",
            {"base_state_sync_functions": syncs, "p_writer_functions": writers},
        )


def _vertical_risk_flags(summary: dict[str, Any], flags: list[RiskFlag]) -> None:
    if summary.get("status") != "available":
        return
    for level in summary.get("levels", []):
        if not isinstance(level, dict):
            continue
        level_id = level.get("level")
        p_diff = _finite_float(level.get("p_mean_diff"))
        companion_diffs = level.get("companion_mean_diffs") or {}
        pb_diff = _finite_float(companion_diffs.get("PB"))
        full_diff = _finite_float(companion_diffs.get("P+PB"))
        mub_diff = _finite_float(companion_diffs.get("MUB"))
        phb_diff = _finite_float(companion_diffs.get("PHB"))
        mu_diff = _finite_float(companion_diffs.get("MU"))
        source_delta = _finite_float(level.get("source_evolution_delta"))

        if (
            p_diff is not None
            and abs(p_diff) >= P_BIAS_ABS_TOLERANCE
            and all(
                value is not None and abs(value) <= WEAK_COMPANION_ABS_TOLERANCE
                for value in (pb_diff, mub_diff, phb_diff)
            )
        ):
            _add_flag(
                flags,
                "p_bias_with_weak_base_state_companion_diffs",
                "warning",
                "large perturbation-P bias appears while PB/MUB/PHB differ weakly",
                {
                    "level": level_id,
                    "p_mean_diff": p_diff,
                    "pb_mean_diff": pb_diff,
                    "mub_mean_diff": mub_diff,
                    "phb_mean_diff": phb_diff,
                    "mu_mean_diff": mu_diff,
                },
            )

        if p_diff is not None and full_diff is not None and abs(p_diff) >= P_BIAS_ABS_TOLERANCE and _close(full_diff, p_diff):
            _add_flag(
                flags,
                "full_pressure_inherits_perturbation_p_bias",
                "warning",
                "P+PB target-region bias closely matches perturbation-P bias",
                {
                    "level": level_id,
                    "p_mean_diff": p_diff,
                    "p_plus_pb_mean_diff": full_diff,
                    "pb_mean_diff": pb_diff,
                },
            )

        if p_diff is not None and source_delta is not None and _close(source_delta, p_diff):
            _add_flag(
                flags,
                "source_start_evolution_delta_equals_p_bias",
                "info",
                "source-start evolution delta equals the P mean difference",
                {
                    "level": level_id,
                    "p_mean_diff": p_diff,
                    "source_evolution_delta": source_delta,
                },
            )


def audit_pressure_refresh_producer(
    *,
    candidate: Path | None = None,
    source_root: Path = Path("."),
    vertical_bias_json: Path | None = None,
) -> ProducerAudit:
    metadata = _metadata_contract(candidate)
    inventory = _source_inventory(source_root)
    vertical_payload = _load_vertical_bias(vertical_bias_json)
    vertical_summary = _vertical_bias_summary(vertical_payload)

    flags: list[RiskFlag] = []
    _metadata_risk_flags(metadata, flags)
    _code_risk_flags(inventory, flags)
    _vertical_risk_flags(vertical_summary, flags)

    return ProducerAudit(
        status="computed",
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        metadata_contract=metadata,
        code_producer_inventory=inventory,
        formula_risk_flags=flags,
        vertical_bias_summary=vertical_summary,
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


def report_to_dict(report: ProducerAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(report: ProducerAudit, *, pretty: bool = False) -> str:
    return json.dumps(
        report_to_dict(report),
        indent=2 if pretty else None,
        allow_nan=False,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate",
        type=Path,
        help="Optional candidate NetCDF file to inspect for pressure-refresh metadata",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=Path("."),
        help="Repository/source root containing src/dynamics and include/tywrf/dynamics",
    )
    parser.add_argument(
        "--vertical-bias-json",
        type=Path,
        help="Optional JSON report from audit_pressure_refresh_vertical_bias.py",
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
        report = audit_pressure_refresh_producer(
            candidate=args.candidate,
            source_root=args.source_root,
            vertical_bias_json=args.vertical_bias_json,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))

    text = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
