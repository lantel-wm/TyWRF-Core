#!/usr/bin/env python
"""Generate explicitly marked d02 skeleton cycle candidates."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

import netCDF4

try:
    from tools.baseline_candidate import CORE_WRFOUT_VARIABLES, build_baseline_candidate
    from tools.run_6h_cycle_test import format_wrf_time, parse_wrf_time, wrfout_filename
except ModuleNotFoundError:
    from baseline_candidate import CORE_WRFOUT_VARIABLES, build_baseline_candidate
    from run_6h_cycle_test import format_wrf_time, parse_wrf_time, wrfout_filename


SKELETON_MODES = ("persistence", "identity")
SUPPORTED_DOMAIN = "d02"
DEFAULT_CYCLE_HOURS = 6
SUPPORTED_HOURS = DEFAULT_CYCLE_HOURS


@dataclass(frozen=True)
class SkeletonCycleReport:
    status: str
    mode: str
    skeleton: bool
    not_physical: bool
    integrator_output: bool
    validation_gate_only: bool
    domain: str
    start: str
    end: str
    hours: int
    minutes: int
    source: str
    candidate: str
    reference_end: str
    candidate_kind: str
    copied_variables: list[str]
    missing_variables: list[str]
    copied_dimensions: list[str]
    times_rewritten: bool
    dx_m: float | None
    dy_m: float | None
    d02_resolution_check: str
    suggested_next_step: dict[str, str]
    message: str


def normalize_skeleton_mode(mode: str) -> str:
    normalized = mode.strip().lower().replace("-", "_")
    if normalized not in SKELETON_MODES:
        raise ValueError(
            "skeleton mode must be one of: " + ", ".join(SKELETON_MODES)
        )
    return normalized


def _resolve_cycle_times(
    start: str | datetime,
    *,
    end: str | datetime | None,
    hours: int | None,
    minutes: int | None = None,
) -> tuple[datetime, datetime]:
    start_time = parse_wrf_time(start) if isinstance(start, str) else start
    if isinstance(end, str):
        end_time = parse_wrf_time(end)
    elif isinstance(end, datetime):
        end_time = end
    else:
        end_time = start_time + timedelta(
            minutes=_resolve_duration_minutes(hours=hours, minutes=minutes)
        )
    duration_seconds = int((end_time - start_time).total_seconds())
    if duration_seconds <= 0:
        raise ValueError("skeleton cycle length must be positive")
    if duration_seconds % 60 != 0:
        raise ValueError("skeleton cycle length must align to whole minutes")
    return start_time, end_time


def _resolve_duration_minutes(*, hours: int | None, minutes: int | None) -> int:
    if minutes is not None:
        if minutes <= 0:
            raise ValueError("cycle length minutes must be positive")
        return minutes
    if hours is None or hours <= 0:
        raise ValueError("cycle length must be positive")
    return hours * 60


def _cycle_gate_command(
    reference_dir: Path,
    candidate_dir: Path,
    *,
    start_time: str,
    end_time: str,
    domain: str,
    interval_minutes: int,
) -> str:
    command = (
        "UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python "
        f"tools/cycle_gate.py --reference-dir {reference_dir} --candidate-dir {candidate_dir} "
        f"--start {start_time} --end {end_time} --domain {domain} --pretty"
    )
    if interval_minutes != DEFAULT_CYCLE_HOURS * 60:
        command += f" --interval-minutes {interval_minutes}"
    return command


def _set_skeleton_attrs(path: Path, report: SkeletonCycleReport) -> None:
    attrs = {
        "TYWRF_CANDIDATE_KIND": report.candidate_kind,
        "TYWRF_SKELETON": str(report.skeleton).lower(),
        "TYWRF_SKELETON_MODE": report.mode,
        "TYWRF_SKELETON_NOT_PHYSICAL": str(report.not_physical).lower(),
        "TYWRF_NOT_PHYSICAL": str(report.not_physical).lower(),
        "TYWRF_INTEGRATOR_OUTPUT": str(report.integrator_output).lower(),
        "TYWRF_VALIDATION_GATE_ONLY": str(report.validation_gate_only).lower(),
        "TYWRF_EXPECTED_TO_MEET_THRESHOLDS": "false",
        "TYWRF_CYCLE_START": report.start,
        "TYWRF_CYCLE_END": report.end,
        "TYWRF_CYCLE_HOURS": report.hours,
        "TYWRF_CYCLE_MINUTES": report.minutes,
        "TYWRF_CANDIDATE_SOURCE": report.source,
        "TYWRF_CANDIDATE_MESSAGE": report.message,
    }
    with netCDF4.Dataset(path, "a") as dataset:
        for name, value in attrs.items():
            dataset.setncattr(name, value)


def build_skeleton_cycle_candidate(
    reference_dir: Path,
    candidate_dir: Path,
    *,
    start: str | datetime,
    end: str | datetime | None = None,
    hours: int | None = SUPPORTED_HOURS,
    minutes: int | None = None,
    domain: str = SUPPORTED_DOMAIN,
    mode: str = "persistence",
    variables: Iterable[str] = CORE_WRFOUT_VARIABLES,
    allow_missing: bool = False,
) -> SkeletonCycleReport:
    if domain != SUPPORTED_DOMAIN:
        raise ValueError("skeleton cycle driver currently supports d02 only")

    normalized_mode = normalize_skeleton_mode(mode)
    start_time, end_time = _resolve_cycle_times(
        start,
        end=end,
        hours=hours,
        minutes=minutes,
    )
    start_text = format_wrf_time(start_time)
    end_text = format_wrf_time(end_time)
    duration_minutes = int((end_time - start_time).total_seconds() // 60)

    source = Path(reference_dir) / wrfout_filename(domain, start_time)
    reference_end = Path(reference_dir) / wrfout_filename(domain, end_time)
    candidate = Path(candidate_dir) / wrfout_filename(domain, end_time)

    baseline_metadata = build_baseline_candidate(
        source,
        reference_end,
        candidate,
        domain=domain,
        start_time=start_time,
        end_time=end_time,
        mode="persistence",
        variables=variables,
        allow_missing=allow_missing,
    )

    message = (
        "Skeleton candidate generated from the cycle-start reference state. "
        "It is not physical, not a TyWRF-Core integrator result, and should be "
        "used only to wire the d02 candidate path into the cycle gate."
    )
    report = SkeletonCycleReport(
        status="skeleton_candidate_generated",
        mode=normalized_mode,
        skeleton=True,
        not_physical=True,
        integrator_output=False,
        validation_gate_only=True,
        domain=domain,
        start=start_text,
        end=end_text,
        hours=duration_minutes // 60,
        minutes=duration_minutes,
        source=str(source),
        candidate=str(candidate),
        reference_end=str(reference_end),
        candidate_kind="skeleton_candidate",
        copied_variables=baseline_metadata.copied_variables,
        missing_variables=baseline_metadata.missing_variables,
        copied_dimensions=baseline_metadata.copied_dimensions,
        times_rewritten=baseline_metadata.times_rewritten,
        dx_m=baseline_metadata.dx_m,
        dy_m=baseline_metadata.dy_m,
        d02_resolution_check=baseline_metadata.d02_resolution_check,
        suggested_next_step={
            "tool": "tools/cycle_gate.py",
            "command": _cycle_gate_command(
                Path(reference_dir),
                Path(candidate_dir),
                start_time=start_text,
                end_time=end_text,
                domain=domain,
                interval_minutes=duration_minutes,
            ),
        },
        message=message,
    )
    _set_skeleton_attrs(candidate, report)
    return report


def report_to_json(report: SkeletonCycleReport, *, pretty: bool = False) -> str:
    return json.dumps(asdict(report), indent=2 if pretty else None)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True, type=Path, help="Directory with WRF reference files")
    parser.add_argument("--candidate-dir", required=True, type=Path, help="Directory for skeleton candidates")
    parser.add_argument("--start", required=True, help="Cycle start time, for example 2025-07-26_00:00:00")
    parser.add_argument("--end", help="Cycle end time; defaults to start + 6 h unless --minutes is set")
    parser.add_argument("--hours", type=int, default=SUPPORTED_HOURS, help="Cycle length in hours when --end/--minutes are omitted")
    parser.add_argument("--minutes", type=int, help="Cycle length in minutes when --end is omitted")
    parser.add_argument("--domain", default=SUPPORTED_DOMAIN, choices=(SUPPORTED_DOMAIN,), help="Only d02 is supported")
    parser.add_argument("--mode", default="persistence", choices=SKELETON_MODES, help="Skeleton generation mode")
    parser.add_argument("--variables", nargs="+", default=list(CORE_WRFOUT_VARIABLES))
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument("--output", type=Path, help="Write report JSON to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = build_skeleton_cycle_candidate(
            args.reference_dir,
            args.candidate_dir,
            start=args.start,
            end=args.end,
            hours=args.hours,
            minutes=args.minutes,
            domain=args.domain,
            mode=args.mode,
            variables=args.variables,
            allow_missing=args.allow_missing,
        )
    except (FileNotFoundError, KeyError, ValueError) as exc:
        parser.error(str(exc))

    output = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
