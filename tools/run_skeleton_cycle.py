#!/usr/bin/env python
"""Generate d01+d02 6 h skeleton cycle candidates and a machine-readable report."""

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
    from tools.skeleton_cycle_driver import (
        SUPPORTED_HOURS,
        build_skeleton_cycle_candidate,
        normalize_skeleton_mode,
    )
except ModuleNotFoundError:
    from baseline_candidate import CORE_WRFOUT_VARIABLES, build_baseline_candidate
    from run_6h_cycle_test import format_wrf_time, parse_wrf_time, wrfout_filename
    from skeleton_cycle_driver import (
        SUPPORTED_HOURS,
        build_skeleton_cycle_candidate,
        normalize_skeleton_mode,
    )


DEFAULT_DOMAINS = ("d01", "d02")


@dataclass(frozen=True)
class SkeletonDomainReport:
    domain: str
    status: str
    mode: str
    skeleton: bool
    not_physical: bool
    integrator_output: bool
    validation_gate_only: bool
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
    suggested_gate_command: str
    message: str


@dataclass(frozen=True)
class SkeletonCycleRunReport:
    status: str
    mode: str
    skeleton: bool
    not_physical: bool
    integrator_output: bool
    validation_gate_only: bool
    start: str
    end: str
    hours: int
    reference_dir: str
    candidate_dir: str
    domains: list[SkeletonDomainReport]
    integrator: dict[str, str | bool]
    suggested_gate_commands: dict[str, str]
    message: str


def _resolve_cycle_times(
    start: str | datetime,
    *,
    end: str | datetime | None,
    hours: int,
) -> tuple[datetime, datetime]:
    if hours != SUPPORTED_HOURS:
        raise ValueError("skeleton cycle runner currently supports only 6 h cycles")

    start_time = parse_wrf_time(start) if isinstance(start, str) else start
    end_time = (
        parse_wrf_time(end)
        if isinstance(end, str)
        else end
        if isinstance(end, datetime)
        else start_time + timedelta(hours=SUPPORTED_HOURS)
    )
    if end_time <= start_time or (end_time - start_time) != timedelta(hours=SUPPORTED_HOURS):
        raise ValueError("skeleton cycle runner currently supports exactly one 6 h cycle")
    return start_time, end_time


def _validate_domains(domains: Iterable[str]) -> list[str]:
    requested = list(domains)
    invalid = [domain for domain in requested if domain not in DEFAULT_DOMAINS]
    if invalid:
        raise ValueError(f"unsupported domain(s): {', '.join(invalid)}")
    if not requested:
        raise ValueError("at least one domain is required")
    return requested


def _cycle_gate_command(
    reference_dir: Path,
    candidate_dir: Path,
    *,
    start: str,
    end: str,
    domain: str,
) -> str:
    return (
        "UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python "
        f"tools/cycle_gate.py --reference-dir {reference_dir} --candidate-dir {candidate_dir} "
        f"--start {start} --end {end} --domain {domain} --pretty"
    )


def _set_orchestrator_attrs(
    candidate: Path,
    *,
    domain: str,
    mode: str,
    start: str,
    end: str,
    source: str,
    message: str,
) -> None:
    attrs = {
        "TYWRF_CANDIDATE_KIND": "skeleton_candidate",
        "TYWRF_SKELETON": "true",
        "TYWRF_SKELETON_MODE": mode,
        "TYWRF_SKELETON_ORCHESTRATOR": "tools/run_skeleton_cycle.py",
        "TYWRF_SKELETON_NOT_PHYSICAL": "true",
        "TYWRF_NOT_PHYSICAL": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "false",
        "TYWRF_VALIDATION_GATE_ONLY": "true",
        "TYWRF_EXPECTED_TO_MEET_THRESHOLDS": "false",
        "TYWRF_CANDIDATE_DOMAIN": domain,
        "TYWRF_CYCLE_START": start,
        "TYWRF_CYCLE_END": end,
        "TYWRF_CANDIDATE_SOURCE": source,
        "TYWRF_CANDIDATE_MESSAGE": message,
    }
    with netCDF4.Dataset(candidate, "a") as dataset:
        for name, value in attrs.items():
            dataset.setncattr(name, value)


def _domain_message(domain: str) -> str:
    return (
        f"{domain} skeleton candidate generated from the cycle-start reference state. "
        "It is not physical, not a TyWRF-Core integrator result, and should be "
        "used only to wire the dual-domain 6 h cycle path into validation tools."
    )


def _build_d01_report(
    reference_dir: Path,
    candidate_dir: Path,
    *,
    start_time: datetime,
    end_time: datetime,
    mode: str,
    variables: Iterable[str],
    allow_missing: bool,
) -> SkeletonDomainReport:
    start_text = format_wrf_time(start_time)
    end_text = format_wrf_time(end_time)
    source = reference_dir / wrfout_filename("d01", start_time)
    reference_end = reference_dir / wrfout_filename("d01", end_time)
    candidate = candidate_dir / wrfout_filename("d01", end_time)

    metadata = build_baseline_candidate(
        source,
        reference_end,
        candidate,
        domain="d01",
        start_time=start_time,
        end_time=end_time,
        mode="persistence",
        variables=variables,
        allow_missing=allow_missing,
    )
    message = _domain_message("d01")
    _set_orchestrator_attrs(
        candidate,
        domain="d01",
        mode=mode,
        start=start_text,
        end=end_text,
        source=str(source),
        message=message,
    )
    return SkeletonDomainReport(
        domain="d01",
        status="skeleton_candidate_generated",
        mode=mode,
        skeleton=True,
        not_physical=True,
        integrator_output=False,
        validation_gate_only=True,
        source=str(source),
        candidate=str(candidate),
        reference_end=str(reference_end),
        candidate_kind="skeleton_candidate",
        copied_variables=metadata.copied_variables,
        missing_variables=metadata.missing_variables,
        copied_dimensions=metadata.copied_dimensions,
        times_rewritten=metadata.times_rewritten,
        dx_m=metadata.dx_m,
        dy_m=metadata.dy_m,
        d02_resolution_check=metadata.d02_resolution_check,
        suggested_gate_command=_cycle_gate_command(
            reference_dir,
            candidate_dir,
            start=start_text,
            end=end_text,
            domain="d01",
        ),
        message=message,
    )


def _build_d02_report(
    reference_dir: Path,
    candidate_dir: Path,
    *,
    start_time: datetime,
    end_time: datetime,
    mode: str,
    variables: Iterable[str],
    allow_missing: bool,
) -> SkeletonDomainReport:
    report = build_skeleton_cycle_candidate(
        reference_dir,
        candidate_dir,
        start=start_time,
        end=end_time,
        mode=mode,
        variables=variables,
        allow_missing=allow_missing,
    )
    message = _domain_message("d02")
    _set_orchestrator_attrs(
        Path(report.candidate),
        domain="d02",
        mode=mode,
        start=report.start,
        end=report.end,
        source=report.source,
        message=message,
    )
    return SkeletonDomainReport(
        domain="d02",
        status=report.status,
        mode=report.mode,
        skeleton=report.skeleton,
        not_physical=report.not_physical,
        integrator_output=report.integrator_output,
        validation_gate_only=report.validation_gate_only,
        source=report.source,
        candidate=report.candidate,
        reference_end=report.reference_end,
        candidate_kind=report.candidate_kind,
        copied_variables=report.copied_variables,
        missing_variables=report.missing_variables,
        copied_dimensions=report.copied_dimensions,
        times_rewritten=report.times_rewritten,
        dx_m=report.dx_m,
        dy_m=report.dy_m,
        d02_resolution_check=report.d02_resolution_check,
        suggested_gate_command=_cycle_gate_command(
            reference_dir,
            candidate_dir,
            start=report.start,
            end=report.end,
            domain="d02",
        ),
        message=message,
    )


def build_skeleton_cycle_run(
    reference_dir: Path,
    candidate_dir: Path,
    *,
    start: str | datetime,
    end: str | datetime | None = None,
    hours: int = SUPPORTED_HOURS,
    domains: Iterable[str] = DEFAULT_DOMAINS,
    mode: str = "persistence",
    variables: Iterable[str] = CORE_WRFOUT_VARIABLES,
    allow_missing: bool = False,
) -> SkeletonCycleRunReport:
    normalized_mode = normalize_skeleton_mode(mode)
    requested_domains = _validate_domains(domains)
    requested_variables = tuple(variables)
    start_time, end_time = _resolve_cycle_times(start, end=end, hours=hours)
    start_text = format_wrf_time(start_time)
    end_text = format_wrf_time(end_time)

    domain_reports: list[SkeletonDomainReport] = []
    for domain in requested_domains:
        if domain == "d01":
            domain_reports.append(
                _build_d01_report(
                    reference_dir,
                    candidate_dir,
                    start_time=start_time,
                    end_time=end_time,
                    mode=normalized_mode,
                    variables=requested_variables,
                    allow_missing=allow_missing,
                )
            )
        elif domain == "d02":
            domain_reports.append(
                _build_d02_report(
                    reference_dir,
                    candidate_dir,
                    start_time=start_time,
                    end_time=end_time,
                    mode=normalized_mode,
                    variables=requested_variables,
                    allow_missing=allow_missing,
                )
            )

    return SkeletonCycleRunReport(
        status="skeleton_candidates_generated",
        mode=normalized_mode,
        skeleton=True,
        not_physical=True,
        integrator_output=False,
        validation_gate_only=True,
        start=start_text,
        end=end_text,
        hours=SUPPORTED_HOURS,
        reference_dir=str(reference_dir),
        candidate_dir=str(candidate_dir),
        domains=domain_reports,
        integrator={
            "status": "not_run",
            "integrator_output": False,
            "executable": "pending",
            "message": "Replace this Python skeleton orchestration with the TyWRF-Core executable invocation when available.",
        },
        suggested_gate_commands={
            domain_report.domain: domain_report.suggested_gate_command
            for domain_report in domain_reports
        },
        message=(
            "Generated dual-domain skeleton candidates only. These files are "
            "explicitly not physical and not TyWRF-Core integrator output."
        ),
    )


def report_to_json(report: SkeletonCycleRunReport, *, pretty: bool = False) -> str:
    return json.dumps(asdict(report), indent=2 if pretty else None)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True, type=Path, help="Directory with WRF reference files")
    parser.add_argument("--candidate-dir", required=True, type=Path, help="Directory for skeleton candidates")
    parser.add_argument("--start", required=True, help="Cycle start time, for example 2025-07-26_00:00:00")
    parser.add_argument("--end", help="Cycle end time; defaults to start + 6 h")
    parser.add_argument("--hours", type=int, default=SUPPORTED_HOURS, help="Cycle length; only 6 is supported")
    parser.add_argument(
        "--domain",
        choices=DEFAULT_DOMAINS,
        action="append",
        dest="domains",
        help="Domain to generate; repeat for multiple domains. Defaults to d01 and d02.",
    )
    parser.add_argument("--mode", default="persistence", choices=("persistence", "identity"), help="Skeleton mode")
    parser.add_argument("--variables", nargs="+", default=list(CORE_WRFOUT_VARIABLES))
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument("--output", type=Path, help="Write report JSON to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = build_skeleton_cycle_run(
            args.reference_dir,
            args.candidate_dir,
            start=args.start,
            end=args.end,
            hours=args.hours,
            domains=args.domains if args.domains else DEFAULT_DOMAINS,
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
