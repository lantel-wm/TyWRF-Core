import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.audit_pressure_refresh_formula import (
    audit_pressure_refresh_formula,
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
                _ensure_dim(dataset, dims[1], values.shape[0])
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
        "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 2,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT": 2,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT": 2,
    }


def _by_name(payload: dict[str, object]) -> dict[str, dict[str, object]]:
    return {entry["name"]: entry for entry in payload["comparisons"]}


def test_formula_audit_reports_target_metrics_and_diagnosis(tmp_path: Path) -> None:
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate.nc"
    start = tmp_path / "candidate_start.nc"
    ref_p = np.full((1, 2, 4), 10.0)
    ref_pb = np.full((1, 2, 4), 1000.0)
    cand_p = ref_p + 1.0
    cand_p[:, :, 3] = 20.0
    cand_pb = ref_pb - 1.0
    cand_pb[:, :, 3] = 990.0

    _write_wrfout(reference, {"P": ref_p, "PB": ref_pb})
    _write_wrfout(candidate, {"P": cand_p, "PB": cand_pb}, attrs=_base_attrs())
    _write_wrfout(start, {"P": cand_p, "PB": cand_pb})

    payload = json.loads(
        report_to_json(audit_pressure_refresh_formula(reference, candidate, start))
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["diagnostic_only"] is True
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"
    assert payload["metadata"]["disposition"]["candidate_model_pass"] == (
        "not_applicable"
    )
    assert payload["region"]["status"] == "available"
    assert payload["region"]["child_delta"] == {"i": 1, "j": 0}

    comparisons = _by_name(payload)
    p_compare = comparisons["candidate_p_vs_reference_p"]
    assert p_compare["status"] == "computed"
    assert math.isclose(p_compare["bias_mean"], 26.0 / 8.0)
    split = p_compare["region_split"]
    assert split["target_region"]["valid_count"] == 2
    assert split["non_target_region"]["valid_count"] == 6
    assert math.isclose(split["target_region"]["rmse"], 10.0)
    assert math.isclose(split["non_target_region"]["rmse"], 1.0)
    assert split["target_region_dominates_global_error"] is True
    assert math.isclose(split["target_error_fraction"], 200.0 / 206.0)
    assert math.isclose(split["target_region_bias_mean"], 10.0)

    full_compare = comparisons[
        "candidate_full_pressure_vs_reference_full_pressure"
    ]
    assert full_compare["status"] == "computed"
    assert math.isclose(full_compare["rmse"], 0.0)
    assert comparisons["candidate_p_vs_start_p"]["status"] == "computed"
    assert math.isclose(comparisons["candidate_p_vs_start_p"]["rmse"], 0.0)
    assert comparisons["candidate_full_pressure_vs_reference_p"]["status"] == (
        "computed"
    )
    assert comparisons["candidate_p_vs_reference_full_pressure"]["status"] == (
        "computed"
    )

    diagnosis = payload["diagnosis"]
    assert diagnosis["target_region_p_error_dominates"] is True
    assert math.isclose(diagnosis["p_target_region_error_fraction"], 200.0 / 206.0)
    assert diagnosis["full_pressure_better_than_perturbation"] is True
    assert math.isclose(diagnosis["full_pressure_to_perturbation_rmse_ratio"], 0.0)
    assert diagnosis["candidate_p_close_to_start_p"] is True
    assert math.isclose(diagnosis["candidate_p_bias_mean"], 26.0 / 8.0)
    assert math.isclose(diagnosis["target_region_bias_mean"], 10.0)


def test_missing_optional_pb_and_start_are_not_available(tmp_path: Path) -> None:
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate.nc"
    p = np.full((1, 2, 4), 10.0)
    _write_wrfout(reference, {"P": p})
    _write_wrfout(candidate, {"P": p}, attrs=_base_attrs())

    payload = json.loads(
        report_to_json(audit_pressure_refresh_formula(reference, candidate))
    )
    comparisons = _by_name(payload)

    assert comparisons["candidate_p_vs_reference_p"]["status"] == "computed"
    assert (
        comparisons["candidate_full_pressure_vs_reference_full_pressure"]["status"]
        == "not_available"
    )
    assert "reference_end.PB" in comparisons[
        "candidate_full_pressure_vs_reference_full_pressure"
    ]["missing_inputs"]
    assert "candidate.PB" in comparisons[
        "candidate_full_pressure_vs_reference_full_pressure"
    ]["missing_inputs"]
    assert comparisons["candidate_p_vs_start_p"]["status"] == "not_available"
    assert "candidate_start.P" in comparisons["candidate_p_vs_start_p"][
        "missing_inputs"
    ]
    assert (
        comparisons["candidate_full_pressure_vs_start_full_pressure"]["status"]
        == "not_available"
    )
    assert payload["diagnosis"]["full_pressure_better_than_perturbation"] is None
    assert payload["diagnosis"]["candidate_p_close_to_start_p"] is None


def test_formula_audit_cli_writes_json(tmp_path: Path) -> None:
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "formula_audit.json"
    p = np.full((1, 2, 4), 10.0)
    pb = np.full((1, 2, 4), 1000.0)
    _write_wrfout(reference, {"P": p, "PB": pb})
    _write_wrfout(candidate, {"P": p, "PB": pb}, attrs=_base_attrs())

    exit_code = audit_main(
        [
            "--reference-end",
            str(reference),
            "--candidate",
            str(candidate),
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert _by_name(payload)["candidate_p_vs_reference_p"]["status"] == "computed"
