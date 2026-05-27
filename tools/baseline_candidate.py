#!/usr/bin/env python
"""Generate explicit baseline candidate wrfout files for 6 h validation gates."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import netCDF4
import numpy as np

try:
    from tools.tc_diagnostics import MINIMUM_SLP_VARIABLE_CANDIDATES
    from tools.write_core_wrfout import CORE_WRFOUT_VARIABLES, copy_core_wrfout
except ModuleNotFoundError:
    from tc_diagnostics import MINIMUM_SLP_VARIABLE_CANDIDATES
    from write_core_wrfout import CORE_WRFOUT_VARIABLES, copy_core_wrfout


WRF_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"
BASELINE_CANDIDATE_MODES = ("persistence", "reference_copy")
D02_TARGET_DX_M = 2000.0


@dataclass(frozen=True)
class BaselineCandidateMetadata:
    candidate: str
    source: str
    source_role: str
    domain: str
    start_time: str
    end_time: str
    hours: int
    mode: str
    candidate_kind: str
    reference_copy: bool
    integrator_output: bool
    validation_gate_only: bool
    minimum_slp_gate_ready: bool
    expected_to_meet_thresholds: bool
    copied_slp_candidates: list[str]
    copied_variables: list[str]
    missing_variables: list[str]
    copied_dimensions: list[str]
    times_rewritten: bool
    dx_m: float | None
    dy_m: float | None
    d02_resolution_check: str
    message: str


def parse_wrf_time(value: str) -> datetime:
    for fmt in (WRF_TIME_FORMAT, "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            pass
    raise ValueError(f"time must match YYYY-MM-DD_HH:MM:SS: {value}")


def format_wrf_time(value: str | datetime) -> str:
    if isinstance(value, datetime):
        return value.strftime(WRF_TIME_FORMAT)
    return parse_wrf_time(value).strftime(WRF_TIME_FORMAT)


def normalize_candidate_mode(mode: str) -> str:
    normalized = mode.strip().lower().replace("-", "_")
    if normalized not in BASELINE_CANDIDATE_MODES:
        raise ValueError(
            "baseline candidate mode must be one of: "
            + ", ".join(BASELINE_CANDIDATE_MODES)
        )
    return normalized


def _candidate_source_for_mode(
    mode: str,
    reference_start: Path,
    reference_end: Path,
) -> tuple[Path, str]:
    if mode == "persistence":
        return reference_start, "cycle_start"
    if mode == "reference_copy":
        return reference_end, "reference_end"
    raise ValueError(f"unsupported baseline candidate mode: {mode}")


def _read_float_attr(dataset: netCDF4.Dataset, name: str) -> float | None:
    if name not in dataset.ncattrs():
        return None
    try:
        return float(dataset.getncattr(name))
    except (TypeError, ValueError):
        return None


def _resolution_metadata(source: Path, domain: str) -> tuple[float | None, float | None, str]:
    with netCDF4.Dataset(source) as dataset:
        dx_m = _read_float_attr(dataset, "DX")
        dy_m = _read_float_attr(dataset, "DY")

    if domain != "d02":
        return dx_m, dy_m, "not_applicable"

    if dx_m is None or dy_m is None:
        raise ValueError(f"d02 source must carry DX/DY attributes to enforce 2 km resolution: {source}")
    if not math.isclose(dx_m, D02_TARGET_DX_M, rel_tol=0.0, abs_tol=0.5) or not math.isclose(
        dy_m,
        D02_TARGET_DX_M,
        rel_tol=0.0,
        abs_tol=0.5,
    ):
        raise ValueError(
            f"d02 source resolution must remain 2 km, got DX={dx_m:g} DY={dy_m:g}: {source}"
        )
    return dx_m, dy_m, "d02_2km"


def _write_time_variable(dataset: netCDF4.Dataset, end_time: str) -> bool:
    if "Times" not in dataset.variables:
        return False

    variable = dataset.variables["Times"]
    if len(variable.shape) != 2:
        raise ValueError("Times variable must be a 2D WRF char array")
    if variable.shape[0] != 1:
        raise ValueError(
            f"baseline candidate files must have exactly one Time record, got {variable.shape[0]}"
        )

    date_strlen = int(variable.shape[1])
    if date_strlen < len(end_time):
        raise ValueError(f"Times DateStrLen={date_strlen} cannot store {end_time}")

    padded = end_time.ljust(date_strlen)
    variable[0, :] = np.array(list(padded), dtype="S1")
    return True


def _set_candidate_attrs(
    path: Path,
    metadata: BaselineCandidateMetadata,
) -> None:
    attrs = {
        "TYWRF_CANDIDATE_KIND": metadata.candidate_kind,
        "TYWRF_BASELINE_MODE": metadata.mode,
        "TYWRF_BASELINE_SOURCE_ROLE": metadata.source_role,
        "TYWRF_REFERENCE_COPY": str(metadata.reference_copy).lower(),
        "TYWRF_INTEGRATOR_OUTPUT": str(metadata.integrator_output).lower(),
        "TYWRF_VALIDATION_GATE_ONLY": str(metadata.validation_gate_only).lower(),
        "TYWRF_EXPECTED_TO_MEET_THRESHOLDS": str(metadata.expected_to_meet_thresholds).lower(),
        "TYWRF_CANDIDATE_DOMAIN": metadata.domain,
        "TYWRF_CYCLE_START": metadata.start_time,
        "TYWRF_CYCLE_END": metadata.end_time,
        "TYWRF_CYCLE_HOURS": metadata.hours,
        "TYWRF_CANDIDATE_SOURCE": metadata.source,
        "TYWRF_D02_RESOLUTION_CHECK": metadata.d02_resolution_check,
        "TYWRF_CANDIDATE_MESSAGE": metadata.message,
    }

    with netCDF4.Dataset(path, "a") as dataset:
        for name, value in attrs.items():
            dataset.setncattr(name, value)


def build_baseline_candidate(
    reference_start: Path,
    reference_end: Path,
    candidate: Path,
    *,
    domain: str,
    start_time: str | datetime,
    end_time: str | datetime,
    mode: str,
    variables: Iterable[str] = CORE_WRFOUT_VARIABLES,
    allow_missing: bool = False,
) -> BaselineCandidateMetadata:
    normalized_mode = normalize_candidate_mode(mode)
    start_text = format_wrf_time(start_time)
    end_text = format_wrf_time(end_time)
    hours = int((parse_wrf_time(end_text) - parse_wrf_time(start_text)).total_seconds() // 3600)
    if hours <= 0:
        raise ValueError("baseline candidate cycle length must be positive")

    source, source_role = _candidate_source_for_mode(
        normalized_mode,
        Path(reference_start),
        Path(reference_end),
    )
    if not source.exists():
        raise FileNotFoundError(f"baseline candidate source does not exist: {source}")

    dx_m, dy_m, d02_resolution_check = _resolution_metadata(source, domain)
    copy_summary = copy_core_wrfout(
        source,
        candidate,
        variables=variables,
        allow_missing=allow_missing,
    )

    with netCDF4.Dataset(candidate, "a") as dataset:
        times_rewritten = _write_time_variable(dataset, end_text)

    reference_copy = normalized_mode == "reference_copy"
    copied_slp_candidates = [
        name
        for name in MINIMUM_SLP_VARIABLE_CANDIDATES
        if name in copy_summary.copied_variables
    ]
    minimum_slp_gate_ready = bool(copied_slp_candidates)
    message = (
        "Reference-copy candidate is for validation-gate testing only; it is not a "
        "TyWRF-Core integrator result. It can pass the full gate only when copied "
        "variables include a real SLP/MSLP diagnostic."
        if reference_copy
        else "Persistence candidate copies the cycle-start state to the cycle-end valid time; "
        "threshold failures are expected and quantify the baseline gap."
    )
    metadata = BaselineCandidateMetadata(
        candidate=str(candidate),
        source=str(source),
        source_role=source_role,
        domain=domain,
        start_time=start_text,
        end_time=end_text,
        hours=hours,
        mode=normalized_mode,
        candidate_kind="baseline_candidate",
        reference_copy=reference_copy,
        integrator_output=False,
        validation_gate_only=reference_copy,
        minimum_slp_gate_ready=minimum_slp_gate_ready,
        expected_to_meet_thresholds=reference_copy and minimum_slp_gate_ready,
        copied_slp_candidates=copied_slp_candidates,
        copied_variables=copy_summary.copied_variables,
        missing_variables=copy_summary.missing_variables,
        copied_dimensions=copy_summary.copied_dimensions,
        times_rewritten=times_rewritten,
        dx_m=dx_m,
        dy_m=dy_m,
        d02_resolution_check=d02_resolution_check,
        message=message,
    )
    _set_candidate_attrs(candidate, metadata)
    return metadata


def metadata_to_json(metadata: BaselineCandidateMetadata, *, pretty: bool = False) -> str:
    return json.dumps(asdict(metadata), indent=2 if pretty else None)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-start", required=True, type=Path, help="Cycle-start WRF wrfout")
    parser.add_argument("--reference-end", required=True, type=Path, help="Cycle-end WRF wrfout")
    parser.add_argument("--candidate", required=True, type=Path, help="Output candidate wrfout path")
    parser.add_argument("--domain", required=True, choices=("d01", "d02"), help="WRF domain")
    parser.add_argument("--start", required=True, help="Cycle start time")
    parser.add_argument("--end", required=True, help="Cycle end time")
    parser.add_argument(
        "--mode",
        required=True,
        choices=("persistence", "reference-copy", "reference_copy"),
        help="Baseline candidate generation mode",
    )
    parser.add_argument("--variables", nargs="+", default=list(CORE_WRFOUT_VARIABLES))
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument("--metadata-output", type=Path, help="Write metadata JSON to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    metadata = build_baseline_candidate(
        args.reference_start,
        args.reference_end,
        args.candidate,
        domain=args.domain,
        start_time=args.start,
        end_time=args.end,
        mode=args.mode,
        variables=args.variables,
        allow_missing=args.allow_missing,
    )
    output = metadata_to_json(metadata, pretty=args.pretty)
    if args.metadata_output:
        args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
        args.metadata_output.write_text(output + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
