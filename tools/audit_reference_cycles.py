#!/usr/bin/env python
"""Audit KROSA d02 WRF reference cycle coverage and validation fields."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

import netCDF4

try:
    from tools.tc_diagnostics import MINIMUM_SLP_VARIABLE_CANDIDATES
except ModuleNotFoundError:
    from tc_diagnostics import MINIMUM_SLP_VARIABLE_CANDIDATES


WRF_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"
DEFAULT_REFERENCE_DIR = Path(
    "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/"
    "output_gfs_analysis/2025wp12/2025072600/WRF"
)
DEFAULT_DOMAIN = "d02"
DEFAULT_START = "2025-07-26_00:00:00"
DEFAULT_END = "2025-08-02_00:00:00"
DEFAULT_INTERVAL_HOURS = 6

TARGET_VARIABLES = (
    "U",
    "V",
    "T",
    "PH",
    "MU",
    "P",
    "QVAPOR",
    "XLAT",
    "XLONG",
    "U10",
    "V10",
)

SLP_CANDIDATE_VARIABLES = MINIMUM_SLP_VARIABLE_CANDIDATES
PRESSURE_PROXY_VARIABLES = ("PSFC",)


@dataclass(frozen=True)
class CycleAudit:
    valid_time: str
    filename: str
    path: str
    exists: bool
    status: str
    missing_variables: list[str]
    available_slp_candidates: list[str]
    available_pressure_proxy_candidates: list[str]
    message: str | None = None


@dataclass(frozen=True)
class ReferenceCycleAudit:
    status: str
    reference_dir: str
    domain: str
    start_time: str
    end_time: str
    interval_hours: int
    cycle_count: int
    missing_files: list[str]
    missing_variables: dict[str, list[str]]
    available_slp_candidates: dict[str, list[str]]
    available_pressure_proxy_candidates: dict[str, list[str]]
    cycles_without_slp_candidates: list[str]
    unreadable_files: dict[str, str]
    cycles: list[CycleAudit]


def parse_wrf_time(value: str) -> datetime:
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


def _cycle_times(start_time: datetime, end_time: datetime, interval_hours: int) -> list[datetime]:
    if interval_hours <= 0:
        raise ValueError("interval-hours must be positive")
    if end_time < start_time:
        raise ValueError("end time must be greater than or equal to start time")

    times = []
    current = start_time
    interval = timedelta(hours=interval_hours)
    while current <= end_time:
        times.append(current)
        current += interval
    return times


def _cycle_status(missing_variables: list[str], slp_candidates: list[str]) -> str:
    if missing_variables and not slp_candidates:
        return "missing_variables_and_slp_candidates"
    if missing_variables:
        return "missing_variables"
    if not slp_candidates:
        return "missing_slp_candidates"
    return "ok"


def audit_reference_cycles(
    reference_dir: Path,
    *,
    domain: str = DEFAULT_DOMAIN,
    start: str | datetime = DEFAULT_START,
    end: str | datetime = DEFAULT_END,
    interval_hours: int = DEFAULT_INTERVAL_HOURS,
    required_variables: Iterable[str] = TARGET_VARIABLES,
    slp_candidates: Iterable[str] = SLP_CANDIDATE_VARIABLES,
    pressure_proxy_candidates: Iterable[str] = PRESSURE_PROXY_VARIABLES,
) -> ReferenceCycleAudit:
    start_time = parse_wrf_time(start) if isinstance(start, str) else start
    end_time = parse_wrf_time(end) if isinstance(end, str) else end
    required = tuple(required_variables)
    slp_candidate_names = tuple(slp_candidates)
    pressure_proxy_names = tuple(pressure_proxy_candidates)

    cycles: list[CycleAudit] = []
    unreadable_files: dict[str, str] = {}

    for valid_time in _cycle_times(start_time, end_time, interval_hours):
        filename = wrfout_filename(domain, valid_time)
        path = reference_dir / filename
        valid_time_text = format_wrf_time(valid_time)
        if not path.exists():
            cycles.append(
                CycleAudit(
                    valid_time=valid_time_text,
                    filename=filename,
                    path=str(path),
                    exists=False,
                    status="missing_file",
                    missing_variables=[],
                    available_slp_candidates=[],
                    available_pressure_proxy_candidates=[],
                )
            )
            continue

        try:
            with netCDF4.Dataset(path) as dataset:
                names = dataset.variables
                missing_variables = [name for name in required if name not in names]
                available_slp_candidates = [
                    name for name in slp_candidate_names if name in names
                ]
                available_pressure_proxy_candidates = [
                    name for name in pressure_proxy_names if name in names
                ]
        except OSError as exc:
            unreadable_files[filename] = str(exc)
            cycles.append(
                CycleAudit(
                    valid_time=valid_time_text,
                    filename=filename,
                    path=str(path),
                    exists=True,
                    status="unreadable",
                    missing_variables=[],
                    available_slp_candidates=[],
                    available_pressure_proxy_candidates=[],
                    message=str(exc),
                )
            )
            continue

        cycles.append(
            CycleAudit(
                valid_time=valid_time_text,
                filename=filename,
                path=str(path),
                exists=True,
                status=_cycle_status(missing_variables, available_slp_candidates),
                missing_variables=missing_variables,
                available_slp_candidates=available_slp_candidates,
                available_pressure_proxy_candidates=available_pressure_proxy_candidates,
            )
        )

    missing_files = [cycle.filename for cycle in cycles if not cycle.exists]
    missing_variables = {
        cycle.filename: cycle.missing_variables
        for cycle in cycles
        if cycle.missing_variables
    }
    available_slp_candidates = {
        cycle.filename: cycle.available_slp_candidates
        for cycle in cycles
        if cycle.exists and cycle.status != "unreadable"
    }
    available_pressure_proxy_candidates = {
        cycle.filename: cycle.available_pressure_proxy_candidates
        for cycle in cycles
        if cycle.exists and cycle.status != "unreadable"
    }
    cycles_without_slp_candidates = [
        cycle.filename
        for cycle in cycles
        if cycle.exists
        and cycle.status != "unreadable"
        and not cycle.available_slp_candidates
    ]
    failed = bool(
        missing_files or missing_variables or cycles_without_slp_candidates or unreadable_files
    )

    return ReferenceCycleAudit(
        status="failed" if failed else "ok",
        reference_dir=str(reference_dir),
        domain=domain,
        start_time=format_wrf_time(start_time),
        end_time=format_wrf_time(end_time),
        interval_hours=interval_hours,
        cycle_count=len(cycles),
        missing_files=missing_files,
        missing_variables=missing_variables,
        available_slp_candidates=available_slp_candidates,
        available_pressure_proxy_candidates=available_pressure_proxy_candidates,
        cycles_without_slp_candidates=cycles_without_slp_candidates,
        unreadable_files=unreadable_files,
        cycles=cycles,
    )


def report_to_json(report: ReferenceCycleAudit, *, pretty: bool = False) -> str:
    return json.dumps(asdict(report), indent=2 if pretty else None)


def report_to_table(report: ReferenceCycleAudit) -> str:
    headers = (
        "valid_time",
        "status",
        "missing_variables",
        "slp_candidates",
        "pressure_proxies",
        "filename",
    )
    rows = [
        (
            cycle.valid_time,
            cycle.status,
            ",".join(cycle.missing_variables) if cycle.missing_variables else "-",
            ",".join(cycle.available_slp_candidates)
            if cycle.available_slp_candidates
            else "-",
            ",".join(cycle.available_pressure_proxy_candidates)
            if cycle.available_pressure_proxy_candidates
            else "-",
            cycle.filename,
        )
        for cycle in report.cycles
    ]
    widths = [
        max(len(header), *(len(row[index]) for row in rows))
        for index, header in enumerate(headers)
    ]
    lines = [
        " | ".join(header.ljust(widths[index]) for index, header in enumerate(headers)),
        "-+-".join("-" * width for width in widths),
    ]
    lines.extend(
        " | ".join(value.ljust(widths[index]) for index, value in enumerate(row))
        for row in rows
    )
    lines.append(
        f"status={report.status} cycle_count={report.cycle_count} "
        f"missing_files={len(report.missing_files)} "
        f"missing_variable_files={len(report.missing_variables)} "
        f"cycles_without_slp_candidates={len(report.cycles_without_slp_candidates)}"
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=DEFAULT_REFERENCE_DIR,
        help="Directory containing KROSA WRF reference wrfout files",
    )
    parser.add_argument("--domain", default=DEFAULT_DOMAIN)
    parser.add_argument("--start", default=DEFAULT_START)
    parser.add_argument("--end", default=DEFAULT_END)
    parser.add_argument("--interval-hours", type=int, default=DEFAULT_INTERVAL_HOURS)
    parser.add_argument(
        "--required-variables",
        nargs="+",
        default=list(TARGET_VARIABLES),
        help="Variables required in every existing wrfout file",
    )
    parser.add_argument(
        "--slp-candidates",
        nargs="+",
        default=list(SLP_CANDIDATE_VARIABLES),
        help="Real sea-level pressure variables to detect; pressure proxies do not satisfy this gate",
    )
    parser.add_argument(
        "--pressure-proxy-candidates",
        nargs="+",
        default=list(PRESSURE_PROXY_VARIABLES),
        help="Pressure proxy variables to record separately, for example PSFC",
    )
    parser.add_argument("--format", choices=("json", "table"), default="json")
    parser.add_argument("--output", type=Path, help="Write the report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_reference_cycles(
            args.reference_dir,
            domain=args.domain,
            start=args.start,
            end=args.end,
            interval_hours=args.interval_hours,
            required_variables=args.required_variables,
            slp_candidates=args.slp_candidates,
            pressure_proxy_candidates=args.pressure_proxy_candidates,
        )
    except ValueError as exc:
        parser.error(str(exc))

    text = (
        report_to_table(report)
        if args.format == "table"
        else report_to_json(report, pretty=args.pretty)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if report.status == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
