import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.audit_selected_field_state import (
    audit_selected_field_state,
    main as audit_main,
    report_to_json,
)


def _ensure_dim(dataset: netCDF4.Dataset, name: str, size: int) -> None:
    if name in dataset.dimensions:
        assert len(dataset.dimensions[name]) == size
        return
    dataset.createDimension(name, size)


def _variable_dims(name: str, data: np.ndarray) -> tuple[str, ...]:
    if name == "U" and data.ndim == 3:
        return ("Time", "bottom_top", "south_north", "west_east_stag")
    if name == "U" and data.ndim == 2:
        return ("Time", "south_north", "west_east_stag")
    if name == "V" and data.ndim == 3:
        return ("Time", "bottom_top", "south_north_stag", "west_east")
    if name == "V" and data.ndim == 2:
        return ("Time", "south_north_stag", "west_east")
    if name == "MU" and data.ndim == 2:
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
            for dim_name, size in zip(dims[1:], values.shape, strict=True):
                _ensure_dim(dataset, dim_name, size)
            variable = dataset.createVariable(name, "f8", dims)
            if dims and dims[0] == "Time":
                variable[0, ...] = values
            else:
                variable[:] = values


def _base_attrs(
    *,
    from_start: str = "10,20",
    to_start: str = "11,20",
    ratio: int = 1,
) -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "false",
        "TYWRF_GATE_CANDIDATE": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
        "TYWRF_FROM_PARENT_START": from_start,
        "TYWRF_TO_PARENT_START": to_start,
        "TYWRF_PARENT_GRID_RATIO": ratio,
    }


def test_global_metrics_and_first_failing_state_variable(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    ref_u = np.full((2, 3), 10.0)
    cand_u = ref_u.copy()
    cand_u[0, 0] += 3.0
    cand_u[1, 2] = np.nan
    ref_v = np.full((3, 2), 10.0)
    ref_mu = np.full((2, 2), 100.0)

    _write_wrfout(reference, {"U": ref_u, "V": ref_v, "MU": ref_mu})
    _write_wrfout(
        candidate,
        {"U": cand_u, "V": ref_v, "MU": ref_mu},
        attrs=_base_attrs(),
    )

    payload = json.loads(
        report_to_json(
            audit_selected_field_state(
                reference,
                candidate,
                thresholds={"U": 0.05, "V": 0.05, "MU": 0.05},
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["first_failing_state_variable"] == "U"
    assert payload["diagnosis"]["first_failing_state_variable"]["variable"] == "U"
    assert payload["metadata"]["disposition"]["candidate_model_pass"] == (
        "not_applicable"
    )

    u_result = payload["variables"][0]
    assert u_result["variable"] == "U"
    assert u_result["status"] == "threshold_exceeded"
    assert math.isclose(u_result["rmse"], math.sqrt(9.0 / 5.0))
    assert math.isclose(u_result["normalized_rmse"], math.sqrt(9.0 / 5.0) / 10.0)
    assert u_result["max_abs_error"] == 3.0
    assert u_result["valid_count"] == 5
    assert u_result["total_count"] == 6
    assert u_result["reference_finite_count"] == 6
    assert u_result["candidate_finite_count"] == 5


def test_region_split_handles_mass_and_staggered_uv_shapes(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    attrs = _base_attrs(from_start="10,20", to_start="11,21")

    ref_mu = np.full((3, 4), 100.0)
    ref_u = np.full((1, 3, 5), 10.0)
    ref_v = np.full((1, 4, 4), 10.0)
    cand_mu = ref_mu.copy()
    cand_u = ref_u.copy()
    cand_v = ref_v.copy()

    mu_mask = np.zeros(ref_mu.shape, dtype=bool)
    mu_mask[:, -1:] = True
    mu_mask[-1:, :] = True
    u_mask = np.zeros(ref_u.shape[-2:], dtype=bool)
    u_mask[:, -2:] = True
    u_mask[-1:, :] = True
    v_mask = np.zeros(ref_v.shape[-2:], dtype=bool)
    v_mask[:, -1:] = True
    v_mask[-2:, :] = True
    cand_mu[mu_mask] += 5.0
    cand_u[:, u_mask] += 5.0
    cand_v[:, v_mask] += 5.0

    _write_wrfout(reference, {"U": ref_u, "V": ref_v, "MU": ref_mu})
    _write_wrfout(candidate, {"U": cand_u, "V": cand_v, "MU": cand_mu}, attrs=attrs)

    payload = json.loads(report_to_json(audit_selected_field_state(reference, candidate)))
    assert payload["region"]["status"] == "available"
    assert payload["region"]["child_delta"] == {"i": 1, "j": 1}

    by_variable = {entry["variable"]: entry for entry in payload["variables"]}
    assert by_variable["U"]["region_split"]["status"] == "available"
    assert by_variable["U"]["region_split"]["shape_kind"] == "x_staggered"
    assert by_variable["U"]["region_split"]["target_region"]["valid_count"] == 9
    assert by_variable["U"]["region_split"]["overlap_region"]["valid_count"] == 6
    assert by_variable["U"]["region_split"]["target_region_dominates_global_error"] is True
    assert math.isclose(by_variable["U"]["region_split"]["target_error_fraction"], 1.0)

    assert by_variable["V"]["region_split"]["shape_kind"] == "y_staggered"
    assert by_variable["V"]["region_split"]["target_region"]["valid_count"] == 10
    assert by_variable["V"]["region_split"]["overlap_region"]["valid_count"] == 6

    assert by_variable["MU"]["region_split"]["shape_kind"] == "mass_grid"
    assert by_variable["MU"]["region_split"]["target_region"]["valid_count"] == 6
    assert by_variable["MU"]["region_split"]["overlap_region"]["valid_count"] == 6
    assert payload["diagnosis"]["target_region_dominates_error"] == {
        "U": True,
        "V": True,
        "MU": True,
    }
    assert payload["diagnosis"]["exposed_or_interpolated_cells_dominate_state_error"] is True


def test_start_state_comparison_flags_candidate_close_to_start(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    start = tmp_path / "start.nc"
    start_u = np.full((2, 3), 10.0)
    ref_u = start_u + 5.0
    cand_u = start_u.copy()

    _write_wrfout(reference, {"U": ref_u})
    _write_wrfout(candidate, {"U": cand_u})
    _write_wrfout(start, {"U": start_u})

    payload = json.loads(
        report_to_json(
            audit_selected_field_state(
                reference,
                candidate,
                candidate_start_path=start,
                variables=("U",),
            )
        )
    )

    persistence = payload["variables"][0]["persistence"]
    assert persistence["status"] == "available"
    assert persistence["candidate_close_to_start"] is True
    assert persistence["candidate_vs_start"]["rmse"] == 0.0
    assert persistence["reference_end_vs_start"]["rmse"] == 5.0
    assert (
        persistence["candidate_start_distance_fraction_of_reference_evolution"]
        == 0.0
    )
    assert payload["summary"]["candidate_close_to_start"] == {"U": True}
    assert payload["diagnosis"]["candidate_close_to_start"] == {"U": True}


def test_missing_metadata_reports_region_not_available(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    variables = {
        "U": np.ones((2, 3)),
        "V": np.ones((3, 2)),
        "MU": np.ones((2, 2)),
    }
    _write_wrfout(reference, variables)
    _write_wrfout(candidate, variables)

    payload = json.loads(report_to_json(audit_selected_field_state(reference, candidate)))

    assert payload["region"]["status"] == "not_available"
    assert "TYWRF_FROM_PARENT_START" in payload["region"]["missing_metadata"]
    for entry in payload["variables"]:
        assert entry["region_split"]["status"] == "not_available"
    assert payload["status"] == "computed"


def test_staggered_shape_ambiguity_reports_not_available(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    ref_u = np.ones((2, 3))
    cand_u = ref_u.copy()

    _write_wrfout(reference, {"U": ref_u})
    _write_wrfout(candidate, {"U": cand_u}, attrs=_base_attrs())

    payload = json.loads(
        report_to_json(
            audit_selected_field_state(reference, candidate, variables=("U",))
        )
    )

    assert payload["variables"][0]["status"] == "ok"
    assert payload["variables"][0]["region_split"]["status"] == "not_available"
    assert "mass-grid shape" in payload["variables"][0]["region_split"]["message"]
    assert payload["diagnosis"]["target_region_error_fraction"] == {"U": None}


def test_main_writes_json_for_explicit_files(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    start = tmp_path / "start.nc"
    output = tmp_path / "audit.json"
    ref_u = np.full((2, 3), 11.0)
    cand_u = np.full((2, 3), 10.0)
    start_u = np.full((2, 3), 10.0)

    _write_wrfout(reference, {"U": ref_u})
    _write_wrfout(candidate, {"U": cand_u})
    _write_wrfout(start, {"U": start_u})

    exit_code = audit_main(
        [
            str(reference),
            str(candidate),
            "--candidate-start",
            str(start),
            "--variables",
            "U",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["diagnostic_only"] is True
    assert payload["variables"][0]["variable"] == "U"
    assert payload["variables"][0]["persistence"]["status"] == "available"
