import json
import math
import os
from pathlib import Path
import subprocess

import netCDF4
import numpy as np

from tools.audit_selected_field_spatial_alignment import (
    audit_selected_field_spatial_alignment,
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


def test_best_shift_detects_phase_offset_for_end_state_and_evolution(
    tmp_path: Path,
) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    start = tmp_path / "start.nc"
    ref_u = np.arange(30.0).reshape(5, 6)
    start_u = np.zeros_like(ref_u)
    cand_u = np.roll(np.roll(ref_u, 1, axis=-1), -1, axis=-2)

    _write_wrfout(reference_end, {"U": ref_u, "MU": np.zeros((5, 5))})
    _write_wrfout(
        candidate_end,
        {"U": cand_u, "MU": np.zeros((5, 5))},
        attrs=_base_attrs(),
    )
    _write_wrfout(start, {"U": start_u, "MU": np.zeros((5, 5))})

    payload = json.loads(
        report_to_json(
            audit_selected_field_spatial_alignment(
                reference_end,
                candidate_end,
                start,
                variables=("U",),
                max_shift=1,
            )
        )
    )

    result = payload["variables"][0]
    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert result["end_state"]["global"]["best_shift"] == {"di": -1, "dj": 1}
    assert result["evolution"]["global"]["best_shift"] == {"di": -1, "dj": 1}
    assert result["end_state"]["global"]["best"]["rmse"] == 0.0
    assert result["evolution"]["global"]["best"]["rmse"] == 0.0
    assert result["end_state"]["global"]["improved"] is True
    assert "must not be applied" in payload["diagnosis"]["message"]


def test_no_improvement_keeps_zero_shift_best(tmp_path: Path) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    ref_mu = np.full((4, 5), 10.0)
    cand_mu = np.full((4, 5), 12.0)

    _write_wrfout(reference_end, {"MU": ref_mu})
    _write_wrfout(candidate_end, {"MU": cand_mu}, attrs=_base_attrs())

    payload = json.loads(
        report_to_json(
            audit_selected_field_spatial_alignment(
                reference_end,
                candidate_end,
                variables=("MU",),
                max_shift=2,
            )
        )
    )

    global_alignment = payload["variables"][0]["end_state"]["global"]
    assert global_alignment["best_shift"] == {"di": 0, "dj": 0}
    assert global_alignment["improved"] is False
    assert global_alignment["normalized_rmse_improvement_ratio"] == 1.0
    assert global_alignment["shift_count"] == 25
    assert payload["variables"][0]["evolution"]["status"] == "not_requested"


def test_region_split_handles_mass_and_staggered_shapes(tmp_path: Path) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    attrs = _base_attrs(from_start="10,20", to_start="11,21")

    ref_mu = np.ones((3, 4))
    ref_u = np.ones((1, 3, 5))
    ref_v = np.ones((1, 4, 4))
    cand_mu = ref_mu + 1.0
    cand_u = ref_u + 1.0
    cand_v = ref_v + 1.0

    _write_wrfout(reference_end, {"U": ref_u, "V": ref_v, "MU": ref_mu})
    _write_wrfout(candidate_end, {"U": cand_u, "V": cand_v, "MU": cand_mu}, attrs=attrs)

    payload = json.loads(
        report_to_json(
            audit_selected_field_spatial_alignment(
                reference_end,
                candidate_end,
                variables=("U", "V", "MU"),
                max_shift=0,
            )
        )
    )

    by_variable = {entry["variable"]: entry for entry in payload["variables"]}
    u_split = by_variable["U"]["end_state"]["region_split"]
    v_split = by_variable["V"]["end_state"]["region_split"]
    mu_split = by_variable["MU"]["end_state"]["region_split"]

    assert u_split["status"] == "available"
    assert u_split["shape_kind"] == "x_staggered"
    assert u_split["target_region"]["baseline"]["valid_count"] == 9
    assert u_split["overlap_region"]["baseline"]["valid_count"] == 6
    assert v_split["shape_kind"] == "y_staggered"
    assert v_split["target_region"]["baseline"]["valid_count"] == 10
    assert v_split["overlap_region"]["baseline"]["valid_count"] == 6
    assert mu_split["shape_kind"] == "mass_grid"
    assert mu_split["target_region"]["baseline"]["valid_count"] == 6
    assert mu_split["overlap_region"]["baseline"]["valid_count"] == 6


def test_missing_start_reports_end_state_only(tmp_path: Path) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"

    _write_wrfout(reference_end, {"MU": np.ones((2, 3))})
    _write_wrfout(candidate_end, {"MU": np.ones((2, 3))}, attrs=_base_attrs())

    payload = json.loads(
        report_to_json(
            audit_selected_field_spatial_alignment(
                reference_end,
                candidate_end,
                variables=("MU",),
            )
        )
    )

    result = payload["variables"][0]
    assert payload["summary"]["start_provided"] is False
    assert result["end_state"]["status"] == "computed"
    assert result["evolution"]["status"] == "not_requested"
    assert result["evolution"]["candidate_model_pass"] == "not_applicable"


def test_main_writes_cli_json_output(tmp_path: Path) -> None:
    reference_end = tmp_path / "reference_end.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    output = tmp_path / "alignment.json"
    ref_mu = np.arange(12.0).reshape(3, 4)
    cand_mu = np.roll(ref_mu, 1, axis=-1)

    _write_wrfout(reference_end, {"MU": ref_mu})
    _write_wrfout(candidate_end, {"MU": cand_mu}, attrs=_base_attrs())

    exit_code = audit_main(
        [
            str(reference_end),
            str(candidate_end),
            "--variables",
            "MU",
            "--max-shift",
            "1",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["max_shift"] == 1
    assert payload["variables"][0]["end_state"]["global"]["best_shift"] == {
        "di": -1,
        "dj": 0,
    }


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
            "tools/audit_selected_field_spatial_alignment.py",
            "--help",
        ],
        cwd=project_root,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "Diagnostic-only audit" in result.stdout
    assert "--max-shift" in result.stdout
