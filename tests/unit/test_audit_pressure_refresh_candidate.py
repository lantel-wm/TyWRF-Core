import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.audit_pressure_refresh_candidate import (
    audit_pressure_refresh_candidate,
    main as audit_main,
    report_to_json,
)


def _ensure_dim(dataset: netCDF4.Dataset, name: str, size: int) -> None:
    if name not in dataset.dimensions:
        dataset.createDimension(name, size)


def _variable_dims(name: str, data: np.ndarray) -> tuple[str, ...]:
    if data.ndim == 3:
        vertical = "bottom_top_stag" if name == "PHB" else "bottom_top"
        return ("Time", vertical, "south_north", "west_east")
    if data.ndim == 2:
        return ("Time", "south_north", "west_east")
    return tuple(f"{name}_dim_{axis}" for axis in range(data.ndim))


def _write_wrfout(
    path: Path,
    variables: dict[str, np.ndarray],
    *,
    attrs: dict[str, object] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        for key, value in (attrs or {}).items():
            dataset.setncattr(key, value)

        for name, raw_values in variables.items():
            values = np.asarray(raw_values, dtype=np.float64)
            dims = _variable_dims(name, values)
            if values.ndim >= 2:
                _ensure_dim(dataset, "south_north", values.shape[-2])
                _ensure_dim(dataset, "west_east", values.shape[-1])
            if values.ndim == 3:
                vertical_dim = dims[1]
                _ensure_dim(dataset, vertical_dim, values.shape[0])
            if values.ndim < 2:
                for dim_name, size in zip(dims, values.shape, strict=True):
                    _ensure_dim(dataset, dim_name, size)

            variable = dataset.createVariable(name, "f8", dims)
            if dims and dims[0] == "Time":
                variable[0, ...] = values
            else:
                variable[:] = values


def _base_attrs() -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "false",
        "TYWRF_GATE_CANDIDATE": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
        "TYWRF_FROM_PARENT_START": "10,20",
        "TYWRF_TO_PARENT_START": "11,20",
        "TYWRF_PARENT_GRID_RATIO": 1,
        "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS": "applied_to_candidate",
        "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY": "false",
        "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY": "false",
        "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 2,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT": "2",
        "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT": 4,
        "TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT": "0",
        "TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT": 0,
        "TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS": "false",
        "TYWRF_PRESSURE_REFRESH_TOUCHED_HALO_CELLS": "false",
        "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS": 4,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS": 4,
        "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS": 4,
        "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS": 2,
        "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS": 4,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT": "true",
        "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO": "true",
        "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED": "true",
    }


def test_metadata_extraction_and_global_metrics(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    ref_p = np.full((2, 2, 4), 10.0)
    cand_p = ref_p.copy()
    cand_p[0, 0, 0] += 1.0

    _write_wrfout(reference, {"P": ref_p, "PB": ref_p})
    _write_wrfout(candidate, {"P": cand_p, "PB": ref_p}, attrs=_base_attrs())

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_candidate(
                reference,
                candidate,
                variables=("P", "PB"),
                thresholds={"P": 0.01, "PB": 0.05},
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    disposition = payload["metadata"]["disposition"]
    assert disposition["gate_candidate"] is True
    assert disposition["candidate_diagnostic_only"] is False
    assert disposition["integrator_output"] is True
    assert disposition["pressure_refresh_applied"] is True
    assert disposition["pressure_refresh_integration_status"] == "applied_to_candidate"
    assert disposition["pressure_refresh_target_column_count"] == 2
    assert disposition["pressure_refresh_refreshed_column_count"] == 2
    assert disposition["pressure_refresh_invalid_point_count"] == 0
    assert disposition["pressure_refresh_skipped_point_count"] == 0
    assert disposition["pressure_refresh_touched_overlap_cells"] is False
    assert disposition["pressure_refresh_touched_halo_cells"] is False
    assert (
        disposition["pressure_refresh_changed_p_matches_refreshed_point_count"]
        is True
    )
    assert disposition["pressure_refresh_overlap_halo_untouched"] is True
    assert disposition["pressure_refresh_disposition"]["candidate_model_pass"] == (
        "not_applicable"
    )

    p_result = payload["variables"][0]
    assert p_result["variable"] == "P"
    assert p_result["status"] == "threshold_exceeded"
    assert math.isclose(p_result["rmse"], math.sqrt(1.0 / 16.0))
    assert math.isclose(p_result["normalized_rmse"], math.sqrt(1.0 / 1600.0))
    assert p_result["max_abs_error"] == 1.0
    assert p_result["valid_count"] == 16
    assert p_result["total_count"] == 16
    assert p_result["reference_finite_count"] == 16
    assert p_result["candidate_finite_count"] == 16
    assert payload["summary"]["first_failing_pressure_variable"] == "P"
    assert payload["diagnosis"]["first_failing_pressure_variable"]["variable"] == "P"


def test_region_split_infers_exposed_mask_and_p_error_dominance(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    ref_p = np.full((1, 2, 4), 10.0)
    cand_p = ref_p + 1.0
    cand_p[:, :, 3] = 20.0

    _write_wrfout(reference, {"P": ref_p})
    _write_wrfout(candidate, {"P": cand_p}, attrs=_base_attrs())

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_candidate(
                reference,
                candidate,
                variables=("P",),
                thresholds={"P": 0.05},
            )
        )
    )

    assert payload["region"]["status"] == "available"
    assert payload["region"]["child_delta"] == {"i": 1, "j": 0}
    assert payload["region"]["inferred_target_column_count"] == 2
    assert payload["region"]["target_column_count_matches_metadata"] is True

    split = payload["variables"][0]["region_split"]
    assert split["status"] == "available"
    assert split["target_region"]["valid_count"] == 2
    assert split["non_target_region"]["valid_count"] == 6
    assert math.isclose(split["target_region"]["rmse"], 10.0)
    assert math.isclose(split["non_target_region"]["rmse"], 1.0)
    assert split["target_region_dominates_global_error"] is True
    assert math.isclose(split["target_error_fraction"], 200.0 / 206.0)
    assert payload["summary"]["refreshed_region_p_dominates_global_error"] is True
    assert math.isclose(
        payload["summary"]["p_refreshed_region_error_fraction"],
        200.0 / 206.0,
    )


def test_hidden_seam_metadata_remains_diagnostic_only(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    p = np.ones((1, 2, 4))
    attrs = {
        **_base_attrs(),
        "TYWRF_DIAGNOSTIC_ONLY": "true",
        "TYWRF_GATE_CANDIDATE": "false",
        "TYWRF_INTEGRATOR_OUTPUT": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_pressure_refresh_experimental_apply_v0",
        "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY": "true",
        "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY": "true",
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS": "experimental_apply_test_only",
    }

    _write_wrfout(reference, {"P": p})
    _write_wrfout(candidate, {"P": p}, attrs=attrs)

    payload = json.loads(
        report_to_json(audit_pressure_refresh_candidate(reference, candidate, variables=("P",)))
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    disposition = payload["metadata"]["disposition"]
    assert disposition["hidden_seam_pressure_refresh"] is True
    assert disposition["normal_candidate_pressure_refresh"] is False
    assert (
        disposition["pressure_refresh_disposition"]["status"]
        == "hidden_seam_diagnostic_evidence_only"
    )
    blockers = disposition["candidate_gate_blockers"]
    assert "TYWRF_DIAGNOSTIC_ONLY=true" in blockers
    assert "TYWRF_GATE_CANDIDATE=false" in blockers
    assert "TYWRF_INTEGRATOR_OUTPUT=false" in blockers
    assert "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY=true" in blockers
    assert "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY=true" in blockers
    assert (
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS=experimental_apply_test_only"
        in blockers
    )
    assert payload["diagnosis"]["candidate_model_pass"] == "not_applicable"


def test_missing_optional_fields_are_reported_not_available(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    p = np.full((1, 2, 4), 10.0)

    _write_wrfout(reference, {"P": p})
    _write_wrfout(candidate, {"P": p}, attrs=_base_attrs())

    payload = json.loads(report_to_json(audit_pressure_refresh_candidate(reference, candidate)))
    by_variable = {entry["variable"]: entry for entry in payload["variables"]}

    assert by_variable["P"]["status"] == "ok"
    for name in ("PB", "MUB", "PHB", "SLP"):
        assert by_variable[name]["status"] == "not_available"
        assert by_variable[name]["region_split"]["status"] == "not_available"
        assert "variable missing" in by_variable[name]["message"]
    assert payload["summary"]["not_available"] == 4
    assert payload["summary"]["first_failing_pressure_variable"] is None


def test_region_status_not_available_without_required_metadata(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    p = np.ones((1, 2, 4))

    _write_wrfout(reference, {"P": p})
    _write_wrfout(candidate, {"P": p})

    payload = json.loads(
        report_to_json(audit_pressure_refresh_candidate(reference, candidate, variables=("P",)))
    )

    assert payload["region"]["status"] == "not_available"
    assert payload["variables"][0]["region_split"]["status"] == "not_available"
    assert "TYWRF_FROM_PARENT_START" in payload["region"]["missing_metadata"]
    assert payload["status"] == "computed"


def test_main_writes_json_for_explicit_files(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "audit.json"
    p = np.full((1, 2, 4), 10.0)
    _write_wrfout(reference, {"P": p})
    _write_wrfout(candidate, {"P": p}, attrs=_base_attrs())

    exit_code = audit_main(
        [
            str(reference),
            str(candidate),
            "--variables",
            "P",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["variables"][0]["variable"] == "P"
