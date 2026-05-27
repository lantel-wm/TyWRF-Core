#!/usr/bin/env python
"""Plan 6 h cycles and generate explicit baseline candidate files."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

try:
    from tools.baseline_candidate import (
        CORE_WRFOUT_VARIABLES,
        build_baseline_candidate,
        normalize_candidate_mode,
    )
except ModuleNotFoundError:
    from baseline_candidate import (
        CORE_WRFOUT_VARIABLES,
        build_baseline_candidate,
        normalize_candidate_mode,
    )


WRF_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"
DEFAULT_DOMAINS = ("d01", "d02")


@dataclass(frozen=True)
class DomainCyclePlan:
    domain: str
    reference_start: str
    reference_end: str
    candidate_start: str
    candidate_end: str
    reference_start_exists: bool
    reference_end_exists: bool
    candidate_start_exists: bool
    candidate_end_exists: bool


@dataclass(frozen=True)
class CyclePlan:
    status: str
    start_time: str
    end_time: str
    hours: int
    reference_dir: str
    candidate_dir: str
    domains: list[DomainCyclePlan]
    integrator: dict[str, str]
    diagnostics: dict[str, dict[str, str]]


def parse_wrf_time(value: str) -> datetime:
    """Parse WRF wrfout timestamps and ISO-like alternatives."""

    for fmt in (WRF_TIME_FORMAT, "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            pass
    raise ValueError(f"time must match YYYY-MM-DD_HH:MM:SS: {value}")


def format_wrf_time(value: datetime) -> str:
    return value.strftime(WRF_TIME_FORMAT)


def wrfout_filename(domain: str, valid_time: datetime) -> str:
    return f"wrfout_{domain}_{format_wrf_time(valid_time)}"


def _validate_domains(domains: Iterable[str]) -> list[str]:
    requested = list(domains)
    invalid = [domain for domain in requested if domain not in DEFAULT_DOMAINS]
    if invalid:
        raise ValueError(f"unsupported domain(s): {', '.join(invalid)}")
    if not requested:
        raise ValueError("at least one domain is required")
    return requested


def build_cycle_plan(
    reference_dir: Path,
    candidate_dir: Path,
    start: str | datetime,
    *,
    hours: int = 6,
    domains: Iterable[str] = DEFAULT_DOMAINS,
    dry_run: bool = True,
) -> CyclePlan:
    if hours <= 0:
        raise ValueError("cycle length must be positive")

    start_time = parse_wrf_time(start) if isinstance(start, str) else start
    end_time = start_time + timedelta(hours=hours)
    requested_domains = _validate_domains(domains)

    domain_plans = []
    for domain in requested_domains:
        reference_start = reference_dir / wrfout_filename(domain, start_time)
        reference_end = reference_dir / wrfout_filename(domain, end_time)
        candidate_start = candidate_dir / wrfout_filename(domain, start_time)
        candidate_end = candidate_dir / wrfout_filename(domain, end_time)
        domain_plans.append(
            DomainCyclePlan(
                domain=domain,
                reference_start=str(reference_start),
                reference_end=str(reference_end),
                candidate_start=str(candidate_start),
                candidate_end=str(candidate_end),
                reference_start_exists=reference_start.exists(),
                reference_end_exists=reference_end.exists(),
                candidate_start_exists=candidate_start.exists(),
                candidate_end_exists=candidate_end.exists(),
            )
        )

    return CyclePlan(
        status="dry_run" if dry_run else "not_implemented",
        start_time=format_wrf_time(start_time),
        end_time=format_wrf_time(end_time),
        hours=hours,
        reference_dir=str(reference_dir),
        candidate_dir=str(candidate_dir),
        domains=domain_plans,
        integrator={
            "status": "pending",
            "message": "TyWRF-Core executable invocation is not wired into this planner yet.",
        },
        diagnostics={
            "tc": {
                "status": "pending",
                "message": "TC center, MSLP, Vmax, and rainfall diagnostics are not implemented yet.",
            }
        },
    )


def _baseline_cycle_windows(
    start: str,
    *,
    end: str | None,
    hours: int,
    interval_hours: int,
) -> list[tuple[datetime, datetime]]:
    start_time = parse_wrf_time(start)
    if end is None:
        if hours <= 0:
            raise ValueError("cycle length must be positive")
        return [(start_time, start_time + timedelta(hours=hours))]

    if interval_hours <= 0:
        raise ValueError("interval must be positive")
    end_time = parse_wrf_time(end)
    if end_time <= start_time:
        raise ValueError("end time must be after start time")

    interval = timedelta(hours=interval_hours)
    windows = []
    cycle_start = start_time
    cycle_end = cycle_start + interval
    while cycle_end <= end_time:
        windows.append((cycle_start, cycle_end))
        cycle_start = cycle_end
        cycle_end = cycle_start + interval

    if not windows or windows[-1][1] != end_time:
        raise ValueError("end time must align with start + N * interval")
    return windows


def build_baseline_candidate_report(
    reference_dir: Path,
    candidate_dir: Path,
    *,
    start: str,
    end: str | None,
    hours: int,
    interval_hours: int,
    domains: Iterable[str],
    mode: str,
    variables: Iterable[str] = CORE_WRFOUT_VARIABLES,
    allow_missing: bool = False,
) -> dict[str, object]:
    requested_domains = _validate_domains(domains)
    normalized_mode = normalize_candidate_mode(mode)
    windows = _baseline_cycle_windows(
        start,
        end=end,
        hours=hours,
        interval_hours=interval_hours,
    )

    cycles = []
    candidates = []
    for index, (cycle_start, cycle_end) in enumerate(windows):
        domain_metadata = []
        for domain in requested_domains:
            metadata = build_baseline_candidate(
                reference_dir / wrfout_filename(domain, cycle_start),
                reference_dir / wrfout_filename(domain, cycle_end),
                candidate_dir / wrfout_filename(domain, cycle_end),
                domain=domain,
                start_time=cycle_start,
                end_time=cycle_end,
                mode=normalized_mode,
                variables=variables,
                allow_missing=allow_missing,
            )
            metadata_dict = asdict(metadata)
            domain_metadata.append(metadata_dict)
            candidates.append(metadata_dict)

        cycles.append(
            {
                "index": index,
                "start_time": format_wrf_time(cycle_start),
                "end_time": format_wrf_time(cycle_end),
                "domains": domain_metadata,
            }
        )

    return {
        "status": "baseline_candidate_generated",
        "mode": normalized_mode,
        "reference_dir": str(reference_dir),
        "candidate_dir": str(candidate_dir),
        "start_time": format_wrf_time(windows[0][0]),
        "end_time": format_wrf_time(windows[-1][1]),
        "domains": requested_domains,
        "cycle_count": len(windows),
        "cycles": cycles,
        "candidates": candidates,
        "integrator": {
            "status": "not_run",
            "message": "Generated baseline candidate files only; TyWRF-Core integrator was not invoked.",
        },
        "diagnostics": {
            "tc": {
                "status": "pending",
                "message": "Run compare_wrfout.py or the cycle gate against these candidates to quantify the gap.",
            }
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True, type=Path, help="Directory with WRF reference files")
    parser.add_argument("--candidate-dir", required=True, type=Path, help="Directory for TyWRF-Core output")
    parser.add_argument(
        "--domain",
        choices=DEFAULT_DOMAINS,
        action="append",
        dest="domains",
        help="Domain to plan; repeat for multiple domains. Defaults to d01 and d02.",
    )
    parser.add_argument("--start", required=True, help="Cycle start time, for example 2025-07-26_00:00:00")
    parser.add_argument("--end", help="Final cycle end time for gate mode")
    parser.add_argument("--hours", type=int, default=6, help="Cycle length in hours")
    parser.add_argument("--interval", type=int, default=6, help="Cycle interval in hours for gate mode")
    parser.add_argument("--dry-run", action="store_true", help="Print planned inputs without running the integrator")
    parser.add_argument(
        "--mode",
        choices=("persistence", "reference-copy", "reference_copy"),
        help="Generate baseline candidates in this mode instead of running the gate",
    )
    parser.add_argument("--variables", nargs="+", default=list(CORE_WRFOUT_VARIABLES))
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument("--output", type=Path, help="Write the JSON plan to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.mode and not args.dry_run:
        domains = args.domains if args.domains else ["d02"]
        try:
            report = build_baseline_candidate_report(
                args.reference_dir,
                args.candidate_dir,
                start=args.start,
                end=args.end,
                hours=args.hours,
                interval_hours=args.interval,
                domains=domains,
                mode=args.mode,
                variables=args.variables,
                allow_missing=args.allow_missing,
            )
        except (FileNotFoundError, ValueError) as exc:
            parser.error(str(exc))

        report_json = json.dumps(report, indent=2 if args.pretty else None)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(report_json + "\n", encoding="utf-8")
        print(report_json)
        return 0

    if not args.dry_run:
        try:
            from tools.cycle_gate import evaluate_cycles, report_to_json as gate_report_to_json
        except ModuleNotFoundError:
            from cycle_gate import evaluate_cycles, report_to_json as gate_report_to_json

        domains = args.domains if args.domains else ["d02"]
        if len(domains) != 1:
            parser.error("gate mode accepts exactly one --domain; default is d02")
        try:
            report = evaluate_cycles(
                args.reference_dir,
                args.candidate_dir,
                args.start,
                end=args.end,
                hours=None if args.end else args.hours,
                interval_hours=args.interval,
                domain=domains[0],
            )
        except ValueError as exc:
            parser.error(str(exc))

        report_json = gate_report_to_json(report, pretty=args.pretty)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(report_json + "\n", encoding="utf-8")
        print(report_json)
        return 0 if report.status == "passed" else 1

    domains = args.domains if args.domains else list(DEFAULT_DOMAINS)
    try:
        plan = build_cycle_plan(
            args.reference_dir,
            args.candidate_dir,
            args.start,
            hours=args.hours,
            domains=domains,
            dry_run=args.dry_run,
        )
    except ValueError as exc:
        parser.error(str(exc))

    report_json = json.dumps(asdict(plan), indent=2 if args.pretty else None)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0 if args.dry_run else 2


if __name__ == "__main__":
    raise SystemExit(main())
