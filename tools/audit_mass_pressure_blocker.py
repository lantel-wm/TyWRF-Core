#!/usr/bin/env python
"""Diagnostic-only audit for the 00:10 mass/pressure blocker.

This tool reads candidate and WRF reference NetCDF files, compares selected
mass/pressure fields against both reference-start and reference-end states, and
emits evidence for diagnosis only. Reference-end truth is never suitable input
for model updates.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

import netCDF4
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if __package__ in {None, ""} and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


DEFAULT_FIELDS = ("P", "MU", "PSFC")
ALLOWED_FIELDS = ("P", "MU", "PSFC", "PB", "MUB", "PH", "PHB")
UNCHANGED_ABS_TOLERANCE = 0.0
UNCHANGED_NORMALIZED_TOLERANCE = 0.0


class AuditError(RuntimeError):
    """Raised when the audit cannot produce trustworthy diagnostic evidence."""


@dataclass(frozen=True)
class ErrorMetrics:
    rmse: float
    normalized_rmse: float
    max_abs_error: float
    valid_count: int
    total_count: int
    reference_finite_count: int
    candidate_finite_count: int
    sum_squared_error: float


@dataclass(frozen=True)
class MassPressureFieldAudit:
    field: str
    shape: tuple[int, ...]
    rmse: float
    normalized_rmse: float
    max_abs_error: float
    candidate_vs_reference_end: ErrorMetrics
    candidate_vs_reference_start: ErrorMetrics
    reference_end_vs_reference_start: ErrorMetrics
    unchanged_from_reference_start: bool
    changed_in_candidate: bool
    reference_changed_from_start: bool
    candidate_failed_to_change_with_reference: bool


@dataclass(frozen=True)
class MassPressureBlockerAudit:
    candidate: str
    reference_start: str
    reference_end: str
    domain: str
    status: str
    diagnostic_only: bool
    uses_reference_end_truth: bool
    uses_oracle_for_model_update: bool
    advances_00_20: bool
    gate_evidence: bool
    no_gate_pass_claim: bool
    fields: list[MassPressureFieldAudit]
    summary: dict[str, Any]
    diagnosis: dict[str, Any]


def _guard_metadata() -> dict[str, bool]:
    return {
        "diagnostic_only": True,
        "uses_reference_end_truth": True,
        "uses_oracle_for_model_update": False,
        "advances_00_20": False,
        "gate_evidence": False,
        "no_gate_pass_claim": True,
    }


def _require_file(path: Path, role: str) -> None:
    if not path.is_file():
        raise AuditError(f"{role} file is missing: {path}")


def parse_fields(raw_fields: Iterable[str]) -> tuple[str, ...]:
    fields: list[str] = []
    for raw_value in raw_fields:
        for raw_field in str(raw_value).split(","):
            field = raw_field.strip().upper()
            if field:
                fields.append(field)

    if not fields:
        raise AuditError("field list must contain at least one field")

    duplicates = sorted({field for field in fields if fields.count(field) > 1})
    if duplicates:
        raise AuditError(f"duplicate field(s) in field list: {', '.join(duplicates)}")

    invalid = [field for field in fields if field not in ALLOWED_FIELDS]
    if invalid:
        raise AuditError(
            "invalid field(s) for mass/pressure audit: "
            f"{', '.join(invalid)}; allowed fields: {', '.join(ALLOWED_FIELDS)}"
        )

    return tuple(fields)


def _owner_path(owner: str, paths: dict[str, Path]) -> Path:
    return paths[owner]


def _read_required_variable(
    dataset: netCDF4.Dataset,
    *,
    owner: str,
    paths: dict[str, Path],
    field: str,
) -> np.ndarray:
    if field not in dataset.variables:
        raise AuditError(
            f"required variable {field} missing from {owner}: {_owner_path(owner, paths)}"
        )

    try:
        masked = np.ma.asarray(dataset.variables[field][:], dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise AuditError(
            f"required variable {field} in {owner} is not numeric: {exc}"
        ) from exc

    mask = np.ma.getmaskarray(masked)
    values = np.asarray(masked.filled(np.nan), dtype=np.float64)
    if values.size == 0:
        raise AuditError(f"required variable {field} in {owner} has no samples")
    if np.any(mask):
        raise AuditError(f"required variable {field} in {owner} contains masked values")
    if not np.all(np.isfinite(values)):
        raise AuditError(f"required variable {field} in {owner} contains nonfinite values")
    return values


def error_metrics(reference: np.ndarray, candidate: np.ndarray) -> ErrorMetrics:
    if reference.shape != candidate.shape:
        raise AuditError(f"shape mismatch: {reference.shape} != {candidate.shape}")
    if reference.size == 0:
        raise AuditError("cannot compute metrics for an empty field")

    diff = candidate - reference
    sum_squared_error = float(np.sum(diff * diff))
    rmse = float(math.sqrt(sum_squared_error / reference.size))
    scale = float(math.sqrt(np.mean(reference * reference)))
    if scale == 0.0:
        normalized = 0.0 if rmse == 0.0 else math.inf
    else:
        normalized = rmse / scale

    return ErrorMetrics(
        rmse=rmse,
        normalized_rmse=float(normalized),
        max_abs_error=float(np.max(np.abs(diff))),
        valid_count=int(reference.size),
        total_count=int(reference.size),
        reference_finite_count=int(reference.size),
        candidate_finite_count=int(candidate.size),
        sum_squared_error=sum_squared_error,
    )


def _is_unchanged(metrics: ErrorMetrics) -> bool:
    return (
        metrics.max_abs_error <= UNCHANGED_ABS_TOLERANCE
        and metrics.normalized_rmse <= UNCHANGED_NORMALIZED_TOLERANCE
    )


def _audit_field(
    datasets: dict[str, netCDF4.Dataset],
    paths: dict[str, Path],
    field: str,
) -> MassPressureFieldAudit:
    start = _read_required_variable(
        datasets["reference_start"],
        owner="reference_start",
        paths=paths,
        field=field,
    )
    end = _read_required_variable(
        datasets["reference_end"],
        owner="reference_end",
        paths=paths,
        field=field,
    )
    candidate = _read_required_variable(
        datasets["candidate"],
        owner="candidate",
        paths=paths,
        field=field,
    )

    if start.shape != end.shape:
        raise AuditError(
            f"required variable {field} shape mismatch: "
            f"reference_start {start.shape} != reference_end {end.shape}"
        )
    if start.shape != candidate.shape:
        raise AuditError(
            f"required variable {field} shape mismatch: "
            f"reference_start {start.shape} != candidate {candidate.shape}"
        )

    candidate_vs_end = error_metrics(end, candidate)
    candidate_vs_start = error_metrics(start, candidate)
    end_vs_start = error_metrics(start, end)
    unchanged_from_start = _is_unchanged(candidate_vs_start)
    reference_changed = not _is_unchanged(end_vs_start)

    return MassPressureFieldAudit(
        field=field,
        shape=tuple(int(value) for value in start.shape),
        rmse=candidate_vs_end.rmse,
        normalized_rmse=candidate_vs_end.normalized_rmse,
        max_abs_error=candidate_vs_end.max_abs_error,
        candidate_vs_reference_end=candidate_vs_end,
        candidate_vs_reference_start=candidate_vs_start,
        reference_end_vs_reference_start=end_vs_start,
        unchanged_from_reference_start=unchanged_from_start,
        changed_in_candidate=not unchanged_from_start,
        reference_changed_from_start=reference_changed,
        candidate_failed_to_change_with_reference=(
            unchanged_from_start and reference_changed
        ),
    )


def audit_mass_pressure_blocker(
    *,
    candidate_path: Path,
    reference_start_path: Path,
    reference_end_path: Path,
    domain: str,
    fields: Iterable[str] = DEFAULT_FIELDS,
) -> MassPressureBlockerAudit:
    parsed_fields = parse_fields(fields)
    paths = {
        "candidate": candidate_path,
        "reference_start": reference_start_path,
        "reference_end": reference_end_path,
    }
    _require_file(candidate_path, "candidate")
    _require_file(reference_start_path, "reference_start")
    _require_file(reference_end_path, "reference_end")

    try:
        with netCDF4.Dataset(reference_start_path) as reference_start, netCDF4.Dataset(
            reference_end_path
        ) as reference_end, netCDF4.Dataset(candidate_path) as candidate:
            datasets = {
                "reference_start": reference_start,
                "reference_end": reference_end,
                "candidate": candidate,
            }
            field_audits = [
                _audit_field(datasets, paths, field) for field in parsed_fields
            ]
    except OSError as exc:
        raise AuditError(f"failed to open NetCDF input: {exc}") from exc

    unchanged_blockers = [
        audit.field
        for audit in field_audits
        if audit.candidate_failed_to_change_with_reference
    ]
    changed_fields = [audit.field for audit in field_audits if audit.changed_in_candidate]
    unchanged_fields = [
        audit.field for audit in field_audits if audit.unchanged_from_reference_start
    ]
    reference_changed_fields = [
        audit.field for audit in field_audits if audit.reference_changed_from_start
    ]
    summary = {
        "fields": list(parsed_fields),
        "total": len(field_audits),
        "changed_in_candidate": len(changed_fields),
        "unchanged_from_reference_start": len(unchanged_fields),
        "reference_changed_from_start": len(reference_changed_fields),
        "unchanged_despite_reference_change": len(unchanged_blockers),
        "changed_fields": changed_fields,
        "unchanged_fields": unchanged_fields,
        "reference_changed_fields": reference_changed_fields,
        "unchanged_blocker_fields": unchanged_blockers,
        "first_unchanged_mass_pressure_blocker": (
            unchanged_blockers[0] if unchanged_blockers else None
        ),
    }

    return MassPressureBlockerAudit(
        candidate=str(candidate_path),
        reference_start=str(reference_start_path),
        reference_end=str(reference_end_path),
        domain=domain,
        status="computed",
        **_guard_metadata(),
        fields=field_audits,
        summary=summary,
        diagnosis={
            **_guard_metadata(),
            "message": (
                "Mass/pressure blocker audit only: reference-end truth is used "
                "to diagnose the 00:10 failure and must not be used to update "
                "candidate model state or claim a gate pass."
            ),
        },
    )


def _strict_json_value(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _strict_json_value(asdict(value))
    if isinstance(value, dict):
        return {str(key): _strict_json_value(item) for key, item in value.items()}
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


def report_to_dict(report: MassPressureBlockerAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(
    report: MassPressureBlockerAudit,
    *,
    pretty: bool = False,
) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def _failure_payload(
    *,
    candidate_path: Path,
    reference_start_path: Path,
    reference_end_path: Path,
    domain: str,
    message: str,
) -> dict[str, Any]:
    return {
        "candidate": str(candidate_path),
        "reference_start": str(reference_start_path),
        "reference_end": str(reference_end_path),
        "domain": domain,
        "status": "failed",
        **_guard_metadata(),
        "fields": [],
        "summary": {"failed": 1},
        "error": {"message": message},
        "diagnosis": {
            **_guard_metadata(),
            "message": (
                "Audit failed closed; no mass/pressure diagnostic evidence was "
                "accepted from this run."
            ),
        },
    }


def failure_to_json(payload: dict[str, Any], *, pretty: bool = False) -> str:
    return json.dumps(_strict_json_value(payload), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path, help="Candidate 00:10 wrfout")
    parser.add_argument(
        "--reference-start",
        required=True,
        type=Path,
        help="WRF reference start wrfout, normally 00:00",
    )
    parser.add_argument(
        "--reference-end",
        required=True,
        type=Path,
        help="WRF reference end wrfout, normally 00:10",
    )
    parser.add_argument("--domain", required=True, help="Domain label, e.g. d01 or d02")
    parser.add_argument("--output", required=True, type=Path, help="JSON output path")
    parser.add_argument(
        "--fields",
        nargs="+",
        default=list(DEFAULT_FIELDS),
        help=(
            "Mass/pressure fields to audit; accepts space- or comma-separated "
            "names. Default: P MU PSFC"
        ),
    )
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = audit_mass_pressure_blocker(
            candidate_path=args.candidate,
            reference_start_path=args.reference_start,
            reference_end_path=args.reference_end,
            domain=args.domain,
            fields=args.fields,
        )
        payload = report_to_json(report, pretty=args.pretty)
        status = 0
    except AuditError as exc:
        payload = failure_to_json(
            _failure_payload(
                candidate_path=args.candidate,
                reference_start_path=args.reference_start,
                reference_end_path=args.reference_end,
                domain=args.domain,
                message=str(exc),
            ),
            pretty=args.pretty,
        )
        status = 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
