#!/usr/bin/env python
"""Combine KROSA 10 minute moving-nest and field-delta diagnostics."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
from typing import Any, Iterable

try:
    from tools.analyze_cycle_delta import (
        DEFAULT_THRESHOLDS as DEFAULT_DELTA_THRESHOLDS,
        DEFAULT_VARIABLES as DEFAULT_DELTA_VARIABLES,
        analyze_cycle_delta,
        parse_threshold_overrides,
        report_to_dict as delta_report_to_dict,
        resolve_cycle_files,
    )
    from tools.audit_moving_nest import (
        DEFAULT_DOMAIN,
        DEFAULT_PARENT_GRID_RATIO,
        audit_moving_nest,
        report_to_dict as moving_nest_report_to_dict,
        resolve_moving_nest_files,
    )
except ModuleNotFoundError as exc:
    if exc.name != "tools":
        raise
    from analyze_cycle_delta import (
        DEFAULT_THRESHOLDS as DEFAULT_DELTA_THRESHOLDS,
        DEFAULT_VARIABLES as DEFAULT_DELTA_VARIABLES,
        analyze_cycle_delta,
        parse_threshold_overrides,
        report_to_dict as delta_report_to_dict,
        resolve_cycle_files,
    )
    from audit_moving_nest import (
        DEFAULT_DOMAIN,
        DEFAULT_PARENT_GRID_RATIO,
        audit_moving_nest,
        report_to_dict as moving_nest_report_to_dict,
        resolve_moving_nest_files,
    )


DIAGNOSTIC_ONLY_MESSAGE = (
    "diagnostic-only report; no candidate/model pass is created or evaluated"
)


@dataclass(frozen=True)
class ReportFiles:
    start_file: Path
    end_file: Path
    domain: str
    start_time: str
    end_time: str


@dataclass(frozen=True)
class TenMinuteDiagnosticsReport:
    status: str
    diagnostic_only: bool
    candidate_model_pass: str
    disposition: dict[str, Any]
    reference_dir: str
    domain: str
    start_time: str
    end_time: str
    start_file: str
    end_file: str
    summary: dict[str, Any]
    movement_audit: dict[str, Any]
    field_delta: dict[str, Any]


def resolve_report_files(
    reference_dir: Path,
    *,
    domain: str = DEFAULT_DOMAIN,
    start: str,
    end: str,
) -> ReportFiles:
    """Resolve files through both underlying diagnostic APIs and require agreement."""
    movement_start, movement_end, movement_domain, movement_start_time, movement_end_time = (
        resolve_moving_nest_files(
            reference_dir=reference_dir,
            domain=domain,
            start=start,
            end=end,
        )
    )
    delta_start, delta_end, delta_domain, delta_start_time, delta_end_time = resolve_cycle_files(
        reference_dir=reference_dir,
        domain=domain,
        start=start,
        end=end,
    )

    if (movement_start, movement_end) != (delta_start, delta_end):
        raise ValueError(
            "moving-nest and cycle-delta diagnostics resolved different WRF files"
        )
    if (movement_domain, movement_start_time, movement_end_time) != (
        delta_domain,
        delta_start_time,
        delta_end_time,
    ):
        raise ValueError(
            "moving-nest and cycle-delta diagnostics resolved different segment metadata"
        )
    if movement_start_time is None or movement_end_time is None:
        raise ValueError("start and end times are required for the 10 minute report")

    return ReportFiles(
        start_file=movement_start,
        end_file=movement_end,
        domain=movement_domain,
        start_time=movement_start_time,
        end_time=movement_end_time,
    )


def report_10min_diagnostics(
    reference_dir: Path,
    *,
    domain: str = DEFAULT_DOMAIN,
    start: str,
    end: str,
    variables: Iterable[str] = DEFAULT_DELTA_VARIABLES,
    thresholds: dict[str, float] | None = DEFAULT_DELTA_THRESHOLDS,
    parent_grid_ratio: int = DEFAULT_PARENT_GRID_RATIO,
    log_file: Path | None = None,
) -> TenMinuteDiagnosticsReport:
    files = resolve_report_files(
        reference_dir,
        domain=domain,
        start=start,
        end=end,
    )
    movement = audit_moving_nest(
        files.start_file,
        files.end_file,
        domain=files.domain,
        start_time=files.start_time,
        end_time=files.end_time,
        parent_grid_ratio=parent_grid_ratio,
        log_file=log_file,
    )
    delta = analyze_cycle_delta(
        files.start_file,
        files.end_file,
        variables=variables,
        thresholds=thresholds,
        domain=files.domain,
        start_time=files.start_time,
        end_time=files.end_time,
    )

    movement_payload = moving_nest_report_to_dict(movement)
    delta_payload = delta_report_to_dict(delta)
    summary = _combined_summary(movement_payload, delta_payload)
    return TenMinuteDiagnosticsReport(
        status=_combined_status(movement.status, delta.status),
        diagnostic_only=True,
        candidate_model_pass="not_applicable",
        disposition={
            "mode": "diagnostic-only",
            "creates_model_candidate": False,
            "candidate_model_pass": "not_applicable",
            "message": DIAGNOSTIC_ONLY_MESSAGE,
        },
        reference_dir=str(reference_dir),
        domain=files.domain,
        start_time=files.start_time,
        end_time=files.end_time,
        start_file=str(files.start_file),
        end_file=str(files.end_file),
        summary=summary,
        movement_audit=movement_payload,
        field_delta=delta_payload,
    )


def _combined_status(movement_status: str, delta_status: str) -> str:
    if movement_status == "computed" and delta_status == "computed":
        return "computed"
    if "incomplete" in {movement_status, delta_status}:
        return "incomplete"
    return "computed_with_flags"


def _combined_summary(
    movement_payload: dict[str, Any],
    delta_payload: dict[str, Any],
) -> dict[str, Any]:
    movement_summary = movement_payload.get("summary", {})
    delta_summary = delta_payload.get("summary", {})
    resolution = movement_summary.get("resolution", {})
    return {
        "moved": movement_summary.get("moved"),
        "parent_delta": movement_summary.get("parent_delta"),
        "child_cell_delta": movement_summary.get("child_cell_delta"),
        "first_failing_variable": delta_summary.get("first_failing_variable"),
        "strict_threshold_exceeded_count": delta_summary.get(
            "strict_threshold_exceeded"
        ),
        "largest_normalized_deltas": delta_summary.get("largest_normalized_deltas"),
        "d02_2km": resolution.get("d02_2km") if isinstance(resolution, dict) else None,
        "d02_2km_status": (
            resolution.get("status") if isinstance(resolution, dict) else None
        ),
        "d02_resolution": resolution,
        "diagnostic_only": True,
        "candidate_model_pass": "not_applicable",
        "message": DIAGNOSTIC_ONLY_MESSAGE,
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
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def report_to_dict(report: TenMinuteDiagnosticsReport) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(report: TenMinuteDiagnosticsReport, *, pretty: bool = False) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-dir",
        type=Path,
        required=True,
        help="Directory with WRF reference files",
    )
    parser.add_argument("--domain", choices=(DEFAULT_DOMAIN,), default=DEFAULT_DOMAIN)
    parser.add_argument(
        "--start",
        required=True,
        help="Cycle start time, for example 2025-07-26_00:00:00",
    )
    parser.add_argument(
        "--end",
        required=True,
        help="Cycle end time, for example 2025-07-26_00:10:00",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=list(DEFAULT_DELTA_VARIABLES),
        help="Variables to analyze; defaults to strict field gate variables",
    )
    parser.add_argument(
        "--threshold",
        action="append",
        default=[],
        metavar="VARIABLE=VALUE",
        help="Override or add a normalized RMSE threshold; repeat as needed",
    )
    parser.add_argument(
        "--no-thresholds",
        action="store_true",
        help="Report deltas without strict-threshold exceedance flags",
    )
    parser.add_argument(
        "--parent-grid-ratio",
        type=int,
        default=DEFAULT_PARENT_GRID_RATIO,
        help="Parent-to-child horizontal grid ratio; KROSA d02 uses 5",
    )
    parser.add_argument("--log-file", type=Path, help="Optional WRF rsl log to mine for move lines")
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        thresholds = None if args.no_thresholds else parse_threshold_overrides(args.threshold)
        report = report_10min_diagnostics(
            args.reference_dir,
            domain=args.domain,
            start=args.start,
            end=args.end,
            variables=args.variables,
            thresholds=thresholds,
            parent_grid_ratio=args.parent_grid_ratio,
            log_file=args.log_file,
        )
    except ValueError as exc:
        parser.error(str(exc))

    report_json = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
