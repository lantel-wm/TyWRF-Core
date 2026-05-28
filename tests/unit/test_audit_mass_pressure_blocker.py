import json
import math
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.audit_mass_pressure_blocker import (
    AuditError,
    audit_mass_pressure_blocker,
    main as audit_main,
    report_to_json,
)


def _ensure_dim(dataset: netCDF4.Dataset, name: str, size: int) -> None:
    if name in dataset.dimensions:
        assert len(dataset.dimensions[name]) == size
        return
    dataset.createDimension(name, size)


def _variable_dims(name: str, data: np.ndarray) -> tuple[str, ...]:
    if data.ndim == 3:
        return ("Time", "bottom_top", "south_north", "west_east")
    if data.ndim == 2:
        return ("Time", "south_north", "west_east")
    return ("Time", *(f"{name}_dim_{axis}" for axis in range(data.ndim)))


def _write_wrfout(path: Path, variables: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        for name, raw_values in variables.items():
            values = np.asarray(raw_values, dtype=np.float64)
            dims = _variable_dims(name, values)
            for dim_name, size in zip(dims[1:], values.shape, strict=True):
                _ensure_dim(dataset, dim_name, size)
            variable = dataset.createVariable(name, "f8", dims)
            variable[0, ...] = values


def _base_fields(value: float) -> dict[str, np.ndarray]:
    return {
        "P": np.full((2, 2, 2), value),
        "MU": np.full((2, 2), value + 10.0),
        "PSFC": np.full((2, 2), value + 1000.0),
    }


def _write_triplet(
    tmp_path: Path,
    *,
    start_fields: dict[str, np.ndarray] | None = None,
    end_fields: dict[str, np.ndarray] | None = None,
    candidate_fields: dict[str, np.ndarray] | None = None,
) -> tuple[Path, Path, Path]:
    reference_start = tmp_path / "reference_start.nc"
    reference_end = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate.nc"
    _write_wrfout(reference_start, start_fields or _base_fields(10.0))
    _write_wrfout(reference_end, end_fields or _base_fields(12.0))
    _write_wrfout(candidate, candidate_fields or _base_fields(11.0))
    return candidate, reference_start, reference_end


def _audit_payload(
    candidate: Path,
    reference_start: Path,
    reference_end: Path,
    *,
    fields: tuple[str, ...] = ("P", "MU", "PSFC"),
) -> dict[str, object]:
    report = audit_mass_pressure_blocker(
        candidate_path=candidate,
        reference_start_path=reference_start,
        reference_end_path=reference_end,
        domain="d02",
        fields=fields,
    )
    return json.loads(report_to_json(report))


def test_success_reports_metrics_and_diagnostic_guards(tmp_path: Path) -> None:
    candidate, reference_start, reference_end = _write_triplet(tmp_path)
    output = tmp_path / "mass_pressure_audit.json"

    exit_code = audit_main(
        [
            "--candidate",
            str(candidate),
            "--reference-start",
            str(reference_start),
            "--reference-end",
            str(reference_end),
            "--domain",
            "d02",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    by_field = {entry["field"]: entry for entry in payload["fields"]}

    assert exit_code == 0
    assert payload["status"] == "computed"
    assert payload["domain"] == "d02"
    assert payload["diagnostic_only"] is True
    assert payload["uses_reference_end_truth"] is True
    assert payload["uses_oracle_for_model_update"] is False
    assert payload["advances_00_20"] is False
    assert payload["gate_evidence"] is False
    assert payload["no_gate_pass_claim"] is True

    p_audit = by_field["P"]
    assert p_audit["shape"] == [1, 2, 2, 2]
    assert math.isclose(p_audit["rmse"], 1.0)
    assert math.isclose(p_audit["normalized_rmse"], 1.0 / 12.0)
    assert p_audit["max_abs_error"] == 1.0
    assert p_audit["changed_in_candidate"] is True
    assert p_audit["unchanged_from_reference_start"] is False
    assert p_audit["reference_changed_from_start"] is True
    assert p_audit["candidate_failed_to_change_with_reference"] is False
    assert p_audit["candidate_vs_reference_end"]["valid_count"] == 8

    assert payload["summary"]["total"] == 3
    assert payload["summary"]["changed_in_candidate"] == 3
    assert payload["summary"]["unchanged_despite_reference_change"] == 0
    assert payload["diagnosis"]["uses_reference_end_truth"] is True
    assert payload["diagnosis"]["gate_evidence"] is False


def test_unchanged_from_start_detection_marks_blocker(tmp_path: Path) -> None:
    start = _base_fields(10.0)
    end = _base_fields(12.0)
    candidate_fields = {name: values.copy() for name, values in start.items()}
    candidate, reference_start, reference_end = _write_triplet(
        tmp_path,
        start_fields=start,
        end_fields=end,
        candidate_fields=candidate_fields,
    )

    payload = _audit_payload(candidate, reference_start, reference_end)
    by_field = {entry["field"]: entry for entry in payload["fields"]}

    assert by_field["P"]["unchanged_from_reference_start"] is True
    assert by_field["P"]["changed_in_candidate"] is False
    assert by_field["P"]["reference_changed_from_start"] is True
    assert by_field["P"]["candidate_failed_to_change_with_reference"] is True
    assert payload["summary"]["unchanged_from_reference_start"] == 3
    assert payload["summary"]["unchanged_despite_reference_change"] == 3
    assert payload["summary"]["first_unchanged_mass_pressure_blocker"] == "P"


def test_missing_required_variable_fails_closed(tmp_path: Path) -> None:
    candidate_fields = _base_fields(11.0)
    del candidate_fields["PSFC"]
    candidate, reference_start, reference_end = _write_triplet(
        tmp_path,
        candidate_fields=candidate_fields,
    )

    with pytest.raises(AuditError, match="required variable PSFC missing from candidate"):
        audit_mass_pressure_blocker(
            candidate_path=candidate,
            reference_start_path=reference_start,
            reference_end_path=reference_end,
            domain="d02",
        )


def test_shape_mismatch_fails_closed(tmp_path: Path) -> None:
    candidate_fields = _base_fields(11.0)
    candidate_fields["P"] = np.full((3, 2, 2), 11.0)
    candidate, reference_start, reference_end = _write_triplet(
        tmp_path,
        candidate_fields=candidate_fields,
    )

    with pytest.raises(AuditError, match="required variable P shape mismatch"):
        audit_mass_pressure_blocker(
            candidate_path=candidate,
            reference_start_path=reference_start,
            reference_end_path=reference_end,
            domain="d02",
            fields=("P",),
        )


def test_nonfinite_values_fail_closed(tmp_path: Path) -> None:
    candidate_fields = _base_fields(11.0)
    candidate_fields["P"] = candidate_fields["P"].copy()
    candidate_fields["P"][0, 0, 0] = np.nan
    candidate, reference_start, reference_end = _write_triplet(
        tmp_path,
        candidate_fields=candidate_fields,
    )

    with pytest.raises(AuditError, match="required variable P in candidate contains nonfinite"):
        audit_mass_pressure_blocker(
            candidate_path=candidate,
            reference_start_path=reference_start,
            reference_end_path=reference_end,
            domain="d02",
        )


def test_invalid_cli_fields_fail_closed_with_json(tmp_path: Path) -> None:
    candidate, reference_start, reference_end = _write_triplet(tmp_path)
    output = tmp_path / "invalid_fields.json"

    exit_code = audit_main(
        [
            "--candidate",
            str(candidate),
            "--reference-start",
            str(reference_start),
            "--reference-end",
            str(reference_end),
            "--domain",
            "d02",
            "--output",
            str(output),
            "--fields",
            "P,U",
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 1
    assert payload["status"] == "failed"
    assert payload["diagnostic_only"] is True
    assert payload["uses_reference_end_truth"] is True
    assert payload["uses_oracle_for_model_update"] is False
    assert payload["gate_evidence"] is False
    assert "invalid field(s)" in payload["error"]["message"]
    assert "U" in payload["error"]["message"]
