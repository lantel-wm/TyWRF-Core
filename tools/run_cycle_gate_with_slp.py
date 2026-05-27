#!/usr/bin/env python
"""Derive SLP for cycle endpoints, then run the d02 cycle gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

try:
    from tools.cycle_gate import (
        DEFAULT_DOMAIN,
        DEFAULT_INTERVAL_HOURS,
        DEFAULT_INTERVAL_MINUTES,
        CycleGateReport,
        cycle_end_times,
        evaluate_cycles,
        format_wrf_time,
        parse_wrf_time,
        report_to_dict as gate_report_to_dict,
        resolve_end_time,
        wrfout_filename,
    )
    from tools.derive_mslp import DeriveMSLPError, write_derived_mslp
except ModuleNotFoundError:
    from cycle_gate import (
        DEFAULT_DOMAIN,
        DEFAULT_INTERVAL_HOURS,
        DEFAULT_INTERVAL_MINUTES,
        CycleGateReport,
        cycle_end_times,
        evaluate_cycles,
        format_wrf_time,
        parse_wrf_time,
        report_to_dict as gate_report_to_dict,
        resolve_end_time,
        wrfout_filename,
    )
    from derive_mslp import DeriveMSLPError, write_derived_mslp


DEFAULT_DERIVED_DIR = Path("build/validation/cycle_gate_with_slp")


@dataclass(frozen=True)
class DerivedCycleFile:
    role: str
    status: str
    cycle_end_time: str
    source: str
    destination: str
    summary: dict[str, object] | None = None
    message: str | None = None


@dataclass(frozen=True)
class CycleGateWithSLPReport:
    status: str
    domain: str
    reference_dir: str
    candidate_dir: str
    derived_reference_dir: str
    derived_candidate_dir: str
    start_time: str
    end_time: str
    interval_hours: int
    interval_minutes: int
    derivations: list[DerivedCycleFile]
    gate_status: str
    gate_report: CycleGateReport


def _strict_json_value(value):
    if isinstance(value, dict):
        return {key: _strict_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_strict_json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _resolve_interval_minutes(
    *,
    interval_hours: int | None,
    interval_minutes: int | None,
) -> int:
    if interval_minutes is not None:
        if interval_minutes <= 0:
            raise ValueError("interval minutes must be positive")
        return interval_minutes
    if interval_hours is None:
        return DEFAULT_INTERVAL_MINUTES
    if interval_hours <= 0:
        raise ValueError("interval must be positive")
    return interval_hours * 60


def _derive_cycle_file(
    *,
    role: str,
    source_dir: Path,
    destination_dir: Path,
    domain: str,
    cycle_end,
) -> DerivedCycleFile:
    source = source_dir / wrfout_filename(domain, cycle_end)
    destination = destination_dir / wrfout_filename(domain, cycle_end)
    try:
        summary = write_derived_mslp(source, destination)
    except (OSError, DeriveMSLPError) as exc:
        return DerivedCycleFile(
            role=role,
            status="failed",
            cycle_end_time=format_wrf_time(cycle_end),
            source=str(source),
            destination=str(destination),
            message=str(exc),
        )

    return DerivedCycleFile(
        role=role,
        status="passed",
        cycle_end_time=format_wrf_time(cycle_end),
        source=str(source),
        destination=str(destination),
        summary=_strict_json_value(asdict(summary)),
    )


def run_cycle_gate_with_slp(
    reference_dir: Path,
    candidate_dir: Path,
    start: str,
    *,
    end: str | None = None,
    hours: int | None = None,
    interval_hours: int | None = DEFAULT_INTERVAL_HOURS,
    interval_minutes: int | None = None,
    domain: str = DEFAULT_DOMAIN,
    derived_dir: Path = DEFAULT_DERIVED_DIR,
    allow_validation_gate_only: bool = False,
) -> CycleGateWithSLPReport:
    """Derive endpoint SLP files and run the cycle gate on the derived copies."""

    if domain not in {"d01", "d02"}:
        raise ValueError(f"unsupported domain: {domain}")

    start_time = parse_wrf_time(start)
    end_time = resolve_end_time(start_time, end=end, hours=hours)
    resolved_interval_minutes = _resolve_interval_minutes(
        interval_hours=interval_hours,
        interval_minutes=interval_minutes,
    )
    cycle_ends = cycle_end_times(
        start_time,
        end_time,
        interval_hours=None,
        interval_minutes=resolved_interval_minutes,
    )

    derived_reference_dir = derived_dir / "reference"
    derived_candidate_dir = derived_dir / "candidate"
    derivations: list[DerivedCycleFile] = []
    for cycle_end in cycle_ends:
        derivations.append(
            _derive_cycle_file(
                role="reference",
                source_dir=reference_dir,
                destination_dir=derived_reference_dir,
                domain=domain,
                cycle_end=cycle_end,
            )
        )
        derivations.append(
            _derive_cycle_file(
                role="candidate",
                source_dir=candidate_dir,
                destination_dir=derived_candidate_dir,
                domain=domain,
                cycle_end=cycle_end,
            )
        )

    gate_report = evaluate_cycles(
        derived_reference_dir,
        derived_candidate_dir,
        start_time,
        end=end_time,
        interval_hours=None,
        interval_minutes=resolved_interval_minutes,
        domain=domain,
        allow_validation_gate_only=allow_validation_gate_only,
    )
    derivations_failed = any(item.status != "passed" for item in derivations)
    status = "failed" if derivations_failed or gate_report.status != "passed" else "passed"

    return CycleGateWithSLPReport(
        status=status,
        domain=domain,
        reference_dir=str(reference_dir),
        candidate_dir=str(candidate_dir),
        derived_reference_dir=str(derived_reference_dir),
        derived_candidate_dir=str(derived_candidate_dir),
        start_time=format_wrf_time(start_time),
        end_time=format_wrf_time(end_time),
        interval_hours=resolved_interval_minutes // 60,
        interval_minutes=resolved_interval_minutes,
        derivations=derivations,
        gate_status=gate_report.status,
        gate_report=gate_report,
    )


def report_to_dict(report: CycleGateWithSLPReport) -> dict[str, object]:
    return _strict_json_value(
        {
            "status": report.status,
            "domain": report.domain,
            "reference_dir": report.reference_dir,
            "candidate_dir": report.candidate_dir,
            "derived_reference_dir": report.derived_reference_dir,
            "derived_candidate_dir": report.derived_candidate_dir,
            "start_time": report.start_time,
            "end_time": report.end_time,
            "interval_hours": report.interval_hours,
            "interval_minutes": report.interval_minutes,
            "derivations": [asdict(item) for item in report.derivations],
            "gate_status": report.gate_status,
            "gate_report": gate_report_to_dict(report.gate_report),
        }
    )


def report_to_json(report: CycleGateWithSLPReport, *, pretty: bool = False) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True, type=Path, help="Directory with WRF reference files")
    parser.add_argument("--candidate-dir", required=True, type=Path, help="Directory with TyWRF-Core output")
    parser.add_argument("--start", required=True, help="First cycle start time, for example 2025-07-26_00:00:00")
    parser.add_argument("--end", help="Final cycle end time, for example 2025-07-26_06:00:00")
    parser.add_argument("--hours", type=int, help="Duration from --start to the final cycle end")
    parser.add_argument("--interval", type=int, default=DEFAULT_INTERVAL_HOURS, help="Cycle interval in hours")
    parser.add_argument(
        "--interval-minutes",
        type=int,
        help="Cycle interval in minutes; overrides --interval for 10 min progressive validation gates",
    )
    parser.add_argument("--domain", choices=("d01", "d02"), default=DEFAULT_DOMAIN, help="Domain to gate")
    parser.add_argument(
        "--derived-dir",
        type=Path,
        default=DEFAULT_DERIVED_DIR,
        help=f"Root directory for SLP-derived reference/candidate copies; default: {DEFAULT_DERIVED_DIR}",
    )
    parser.add_argument(
        "--allow-validation-gate-only",
        action="store_true",
        help=(
            "Allow candidate files marked TYWRF_VALIDATION_GATE_ONLY=true to exercise "
            "the gate implementation; default strict mode fails them as non-integrator candidates."
        ),
    )
    parser.add_argument("--output", type=Path, help="Write the combined JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = run_cycle_gate_with_slp(
            args.reference_dir,
            args.candidate_dir,
            args.start,
            end=args.end,
            hours=args.hours,
            interval_hours=args.interval,
            interval_minutes=args.interval_minutes,
            domain=args.domain,
            derived_dir=args.derived_dir,
            allow_validation_gate_only=args.allow_validation_gate_only,
        )
    except ValueError as exc:
        parser.error(str(exc))

    report_json = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0 if report.status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
