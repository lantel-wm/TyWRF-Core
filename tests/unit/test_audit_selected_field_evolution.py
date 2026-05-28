import json
import math
import os
from pathlib import Path
import subprocess

import netCDF4
import numpy as np

from tools.audit_selected_field_evolution import (
    audit_selected_field_evolution,
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
    return ("Time", *(f"{name}_dim_{axis}" for axis in range(data.ndim)))


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
            variable[0, ...] = values


def _base_attrs(
    *,
    from_start: str = "10,20",
    to_start: str = "11,20",
    ratio: int = 1,
) -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "true",
        "TYWRF_GATE_CANDIDATE": "false",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
        "TYWRF_FROM_PARENT_START": from_start,
        "TYWRF_TO_PARENT_START": to_start,
        "TYWRF_PARENT_GRID_RATIO": ratio,
    }


def test_global_evolution_metrics_and_capture_fraction(tmp_path: Path) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    start = tmp_path / "start.nc"
    start_u = np.full((2, 3), 10.0)
    ref_u = start_u + 3.0
    cand_u = start_u + 2.0
    cand_u[0, 0] = start_u[0, 0] + 5.0
    cand_u[1, 2] = np.nan

    _write_wrfout(reference_end, {"U": ref_u})
    _write_wrfout(candidate_end, {"U": cand_u}, attrs=_base_attrs())
    _write_wrfout(start, {"U": start_u})

    payload = json.loads(
        report_to_json(
            audit_selected_field_evolution(
                reference_end,
                candidate_end,
                start,
                variables=("U",),
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"
    assert payload["status"] == "computed"

    result = payload["variables"][0]
    assert result["variable"] == "U"
    assert result["status"] == "computed"
    assert result["diagnostic_only"] is True
    assert result["candidate_model_pass"] == "not_applicable"
    assert result["valid_count"] == 5
    assert result["total_count"] == 6
    assert math.isclose(result["rmse"], math.sqrt(8.0 / 5.0))
    assert math.isclose(result["normalized_rmse"], math.sqrt(8.0 / 5.0) / 3.0)
    assert result["max_abs_error"] == 2.0

    amplitude = result["evolution_amplitude"]
    assert amplitude["status"] == "available"
    assert math.isclose(amplitude["reference_evolution_rms"], 3.0)
    assert math.isclose(amplitude["candidate_evolution_rms"], math.sqrt(41.0 / 5.0))
    assert math.isclose(amplitude["amplitude_ratio"], math.sqrt(41.0 / 5.0) / 3.0)
    assert math.isclose(amplitude["capture_fraction"], 13.0 / 15.0)
    assert payload["summary"]["capture_fraction"] == {"U": amplitude["capture_fraction"]}


def test_target_overlap_split_handles_staggered_evolution_errors(tmp_path: Path) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    start = tmp_path / "start.nc"
    attrs = _base_attrs(from_start="10,20", to_start="11,21")

    start_mu = np.zeros((3, 4))
    start_u = np.zeros((1, 3, 5))
    start_v = np.zeros((1, 4, 4))
    ref_mu = np.ones_like(start_mu)
    ref_u = np.ones_like(start_u)
    ref_v = np.ones_like(start_v)
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

    _write_wrfout(reference_end, {"U": ref_u, "V": ref_v, "MU": ref_mu})
    _write_wrfout(candidate_end, {"U": cand_u, "V": cand_v, "MU": cand_mu}, attrs=attrs)
    _write_wrfout(start, {"U": start_u, "V": start_v, "MU": start_mu})

    payload = json.loads(
        report_to_json(
            audit_selected_field_evolution(
                reference_end,
                candidate_end,
                start,
                include_level_summary=True,
            )
        )
    )

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

    u_levels = by_variable["U"]["level_summary"]
    assert u_levels["status"] == "available"
    assert u_levels["candidate_model_pass"] == "not_applicable"
    assert u_levels["level_count"] == 1
    assert u_levels["worst_levels"][0]["level"] == 0


def test_missing_start_variable_reports_not_available(tmp_path: Path) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    reference_start = tmp_path / "reference_start.nc"
    candidate_start = tmp_path / "candidate_start.nc"

    _write_wrfout(reference_end, {"U": np.ones((2, 3))})
    _write_wrfout(candidate_end, {"U": np.ones((2, 3))}, attrs=_base_attrs())
    _write_wrfout(reference_start, {"U": np.zeros((2, 3))})
    _write_wrfout(candidate_start, {"MU": np.zeros((2, 3))})

    payload = json.loads(
        report_to_json(
            audit_selected_field_evolution(
                reference_end,
                candidate_end,
                reference_start_path=reference_start,
                candidate_start_path=candidate_start,
                variables=("U",),
            )
        )
    )

    assert payload["status"] == "computed_with_flags"
    result = payload["variables"][0]
    assert result["status"] == "not_available"
    assert result["reference_start_present"] is True
    assert result["candidate_start_present"] is False
    assert "candidate-start" in result["message"]
    assert result["candidate_model_pass"] == "not_applicable"
    assert result["region_split"]["candidate_model_pass"] == "not_applicable"


def test_main_writes_cli_json_output(tmp_path: Path) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    start = tmp_path / "start.nc"
    output = tmp_path / "audit.json"
    start_u = np.full((2, 3), 10.0)

    _write_wrfout(reference_end, {"U": start_u + 2.0})
    _write_wrfout(candidate_end, {"U": start_u + 1.0}, attrs=_base_attrs())
    _write_wrfout(start, {"U": start_u})

    exit_code = audit_main(
        [
            str(reference_end),
            str(candidate_end),
            str(start),
            "--variables",
            "U",
            "--level-summary",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["shared_start_file"] is True
    assert payload["variables"][0]["variable"] == "U"
    assert payload["variables"][0]["level_summary"]["status"] == "not_available"


def test_direct_script_help_invocation() -> None:
    project_root = Path(__file__).resolve().parents[2]
    env = {
        **os.environ,
        "UV_CACHE_DIR": ".uv-cache",
        "UV_PYTHON_INSTALL_DIR": ".uv-python",
    }
    result = subprocess.run(
        [
            "uv",
            "run",
            "python",
            "tools/audit_selected_field_evolution.py",
            "--help",
        ],
        cwd=project_root,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "Audit selected field evolution errors" in result.stdout
    assert "--reference-start" in result.stdout
