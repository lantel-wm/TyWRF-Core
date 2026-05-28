import json
import math
from pathlib import Path
import subprocess
import sys

import netCDF4
import numpy as np

from tools.audit_selected_field_wind_source_attribution import (
    audit_selected_field_wind_source_attribution,
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
    to_start: str = "11,21",
    ratio: int = 1,
) -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "false",
        "TYWRF_GATE_CANDIDATE": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
        "TYWRF_CANDIDATE_DOMAIN": "d02",
        "TYWRF_CYCLE_START": "2025-07-26_00:00:00",
        "TYWRF_CYCLE_END": "2025-07-26_00:10:00",
        "DX": 2000.0,
        "DY": 2000.0,
        "TYWRF_D02_RESOLUTION_CHECK": "d02_2km",
        "TYWRF_FROM_PARENT_START": from_start,
        "TYWRF_TO_PARENT_START": to_start,
        "TYWRF_REMAP_FROM_PARENT_START": from_start,
        "TYWRF_REMAP_TO_PARENT_START": to_start,
        "TYWRF_PARENT_GRID_RATIO": ratio,
        "I_PARENT_START": 11,
        "J_PARENT_START": 21,
        "TYWRF_SELECTED_FIELD_TIMELINE_VERSION": "runtime_v0",
        "TYWRF_SELECTED_FIELD_TIMELINE_EVIDENCE_ONLY": "true",
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_COUNT": 2,
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_NAMES": (
            "cycle_start,move_from_to_parent_start"
        ),
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS": (
            "1:cycle_start(time=2025-07-26_00:00:00)|"
            "2:move_from_to_parent_start(from=10_20,to=11_21)"
        ),
    }


def test_staggered_uv_source_metrics_and_explanatory_ratios(tmp_path: Path) -> None:
    d02_start = tmp_path / "d02_start.nc"
    candidate = tmp_path / "candidate.nc"
    reference = tmp_path / "reference.nc"
    start_u = np.full((1, 2, 4), 10.0)
    ref_u = np.full((1, 2, 4), 12.0)
    cand_u = start_u.copy()
    start_v = np.full((1, 3, 3), 20.0)
    ref_v = np.full((1, 3, 3), 23.0)
    cand_v = start_v + 1.0

    _write_wrfout(d02_start, {"U": start_u, "V": start_v, "MU": np.ones((2, 3))})
    _write_wrfout(reference, {"U": ref_u, "V": ref_v, "MU": np.ones((2, 3))})
    _write_wrfout(
        candidate,
        {"U": cand_u, "V": cand_v, "MU": np.ones((2, 3))},
        attrs=_base_attrs(),
    )

    payload = json.loads(
        report_to_json(
            audit_selected_field_wind_source_attribution(
                d02_start,
                candidate,
                reference,
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["gate_evidence"] is False
    assert payload["advances_00_20"] is False
    assert payload["source_pose"]["source_pose_status"] == (
        "raw_start_not_candidate_pose"
    )
    assert payload["source_pose"]["pose_aligned_with_candidate"] is False
    assert payload["summary"]["candidate_vs_d02_start_interpretation"] == (
        "raw_pose_only"
    )
    by_variable = {entry["variable"]: entry for entry in payload["variables"]}

    u_whole = by_variable["U"]["whole_domain"]
    assert u_whole["source_pose_status"] == "raw_start_not_candidate_pose"
    assert u_whole["candidate_vs_d02_start_interpretation"] == "raw_pose_only"
    assert u_whole["candidate_vs_reference"]["rmse"] == 2.0
    assert u_whole["d02_start_vs_reference"]["rmse"] == 2.0
    assert u_whole["candidate_vs_d02_start"]["rmse"] == 0.0
    assert u_whole["candidate_vs_d02_start"]["pose_aligned_with_candidate"] is False
    assert (
        u_whole["candidate_vs_d02_start"][
            "pose_aligned_start_persistence_evidence"
        ]
        is False
    )
    assert u_whole["candidate_vs_reference"]["sum_squared_error"] == 32.0
    assert u_whole["d02_start_vs_reference"]["sum_squared_error"] == 32.0
    assert u_whole["candidate_vs_d02_start"]["sum_squared_error"] == 0.0
    ratios = u_whole["explanatory_ratios"]
    assert ratios["status"] == "raw_pose_only"
    assert ratios["source_pose_status"] == "raw_start_not_candidate_pose"
    assert ratios["pose_aligned_start_persistence_evidence"] is False
    assert ratios["not_shifted_start_persistence_evidence"] is True
    assert ratios["start_state_persistence_fraction"] == 1.0
    assert ratios["candidate_delta_fraction_of_candidate_error"] == 0.0
    assert ratios["candidate_error_cosine_with_start_persistence"] == 1.0

    v_whole = by_variable["V"]["whole_domain"]
    assert v_whole["candidate_vs_reference"]["rmse"] == 2.0
    assert v_whole["d02_start_vs_reference"]["rmse"] == 3.0
    assert v_whole["candidate_vs_d02_start"]["rmse"] == 1.0
    assert math.isclose(
        v_whole["explanatory_ratios"]["start_state_persistence_fraction"],
        9.0 / 4.0,
    )


def test_region_breakdown_reports_whole_exposed_overlap_for_staggered_uv(
    tmp_path: Path,
) -> None:
    d02_start = tmp_path / "d02_start.nc"
    candidate = tmp_path / "candidate.nc"
    reference = tmp_path / "reference.nc"
    ref_u = np.full((1, 3, 5), 12.0)
    ref_v = np.full((1, 4, 4), 22.0)
    start_u = ref_u - 2.0
    start_v = ref_v - 2.0
    cand_u = start_u.copy()
    cand_v = start_v.copy()
    cand_u[:, :, -2:] = ref_u[:, :, -2:]
    cand_u[:, -1:, :] = ref_u[:, -1:, :]
    cand_v[:, :, -1:] = ref_v[:, :, -1:]
    cand_v[:, -2:, :] = ref_v[:, -2:, :]

    _write_wrfout(d02_start, {"U": start_u, "V": start_v, "MU": np.ones((3, 4))})
    _write_wrfout(reference, {"U": ref_u, "V": ref_v, "MU": np.ones((3, 4))})
    _write_wrfout(
        candidate,
        {"U": cand_u, "V": cand_v, "MU": np.ones((3, 4))},
        attrs=_base_attrs(),
    )

    payload = json.loads(
        report_to_json(
            audit_selected_field_wind_source_attribution(
                d02_start,
                candidate,
                reference,
            )
        )
    )
    by_variable = {entry["variable"]: entry for entry in payload["variables"]}

    u_breakdown = by_variable["U"]["region_breakdown"]
    assert u_breakdown["region_breakdown_status"] == "available"
    assert u_breakdown["shape_kind"] == "x_staggered"
    assert u_breakdown["exposed_region"]["candidate_vs_reference"]["valid_count"] == 9
    assert u_breakdown["overlap_region"]["candidate_vs_reference"]["valid_count"] == 6
    assert u_breakdown["overlap_region_dominates_candidate_error"] is True
    assert u_breakdown["exposed_candidate_error_fraction"] == 0.0
    assert (
        u_breakdown["overlap_region"]["explanatory_ratios"][
            "start_state_persistence_fraction"
        ]
        == 1.0
    )

    v_breakdown = by_variable["V"]["region_breakdown"]
    assert v_breakdown["shape_kind"] == "y_staggered"
    assert v_breakdown["exposed_region"]["candidate_vs_reference"]["valid_count"] == 10
    assert v_breakdown["overlap_region"]["candidate_vs_reference"]["valid_count"] == 6
    assert payload["summary"]["region_breakdown_status"] == {
        "U": "available",
        "V": "available",
    }


def test_same_pose_metadata_marks_raw_start_as_pose_aligned_diagnostic(
    tmp_path: Path,
) -> None:
    d02_start = tmp_path / "d02_start.nc"
    candidate = tmp_path / "candidate.nc"
    reference = tmp_path / "reference.nc"
    variables = {
        "U": np.ones((1, 2, 4)),
        "V": np.ones((1, 3, 3)),
        "MU": np.ones((2, 3)),
    }

    _write_wrfout(d02_start, variables)
    _write_wrfout(reference, variables)
    _write_wrfout(
        candidate,
        variables,
        attrs=_base_attrs(from_start="10,20", to_start="10,20"),
    )

    payload = json.loads(
        report_to_json(
            audit_selected_field_wind_source_attribution(
                d02_start,
                candidate,
                reference,
                variables=("U",),
            )
        )
    )

    result = payload["variables"][0]["whole_domain"]
    ratios = result["explanatory_ratios"]
    assert payload["source_pose"]["source_pose_status"] == "same_pose"
    assert payload["source_pose"]["pose_aligned_with_candidate"] is True
    assert payload["summary"]["source_pose_status"] == "same_pose"
    assert result["candidate_vs_d02_start_interpretation"] == (
        "same_pose_start_persistence_diagnostic"
    )
    assert ratios["status"] == "available"
    assert ratios["pose_aligned_start_persistence_evidence"] is True
    assert ratios["not_shifted_start_persistence_evidence"] is False


def test_insufficient_metadata_region_and_pose_fail_closed(tmp_path: Path) -> None:
    d02_start = tmp_path / "d02_start.nc"
    candidate = tmp_path / "candidate.nc"
    reference = tmp_path / "reference.nc"
    start_u = np.full((1, 2, 4), 10.0)
    ref_u = np.full((1, 2, 4), 12.0)
    cand_u = start_u.copy()

    _write_wrfout(d02_start, {"U": start_u, "V": np.ones((1, 3, 3))})
    _write_wrfout(reference, {"U": ref_u, "V": np.ones((1, 3, 3))})
    _write_wrfout(candidate, {"U": cand_u, "V": np.ones((1, 3, 3))})

    payload = json.loads(
        report_to_json(
            audit_selected_field_wind_source_attribution(
                d02_start,
                candidate,
                reference,
                variables=("U",),
            )
        )
    )

    result = payload["variables"][0]
    assert result["status"] == "computed"
    assert result["region_breakdown_status"] == "insufficient_metadata"
    assert result["region_breakdown"]["region_breakdown_status"] == (
        "insufficient_metadata"
    )
    assert result["whole_domain"]["source_pose_status"] == "insufficient_metadata"
    assert result["whole_domain"]["pose_aligned_with_candidate"] is False
    assert result["whole_domain"]["explanatory_ratios"]["status"] == "raw_pose_only"
    assert (
        result["whole_domain"]["explanatory_ratios"][
            "pose_aligned_start_persistence_evidence"
        ]
        is False
    )
    assert payload["source_pose"]["source_pose_status"] == "insufficient_metadata"
    assert "TYWRF_FROM_PARENT_START" in payload["movement_region"]["missing_metadata"]
    assert payload["status"] == "computed_with_flags"


def test_cli_output_includes_metadata_summary_and_parent_context(
    tmp_path: Path,
) -> None:
    d02_start = tmp_path / "d02_start.nc"
    candidate = tmp_path / "candidate.nc"
    reference = tmp_path / "reference.nc"
    parent_start = tmp_path / "parent_start.nc"
    parent_end = tmp_path / "parent_end.nc"
    output = tmp_path / "source_attribution.json"
    variables = {
        "U": np.ones((1, 2, 4)),
        "V": np.ones((1, 3, 3)),
        "MU": np.ones((2, 3)),
    }

    _write_wrfout(d02_start, variables)
    _write_wrfout(reference, variables)
    _write_wrfout(candidate, variables, attrs=_base_attrs())
    _write_wrfout(parent_start, {"U": np.ones((1, 2, 4)), "V": np.ones((1, 3, 3))})
    _write_wrfout(parent_end, {"U": np.ones((1, 2, 4)) + 2.0, "V": np.ones((1, 3, 3))})

    exit_code = audit_main(
        [
            str(d02_start),
            str(candidate),
            str(reference),
            "--d01-parent-start",
            str(parent_start),
            "--d01-parent-end",
            str(parent_end),
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["summary"]["source_time_levels"]["d02_start"]["time_index"] == -1
    assert payload["summary"]["source_time_levels"]["candidate_end"]["time_index"] == -1
    assert payload["summary"]["source_time_levels"]["reference_end"]["time_index"] == -1
    assert payload["summary"]["source_pose_status"] == "raw_start_not_candidate_pose"
    assert payload["summary"]["pose_aligned_with_candidate"] is False
    assert payload["summary"]["candidate_metadata_summary"]["dx"] == 2000.0
    assert payload["summary"]["candidate_metadata_summary"]["dy"] == 2000.0
    assert payload["summary"]["candidate_metadata_summary"]["cycle_end"] == (
        "2025-07-26_00:10:00"
    )
    assert payload["summary"]["candidate_metadata_summary"]["timeline_event_names"] == [
        "cycle_start",
        "move_from_to_parent_start",
    ]
    assert payload["parent_context"]["status"] == "available"
    assert (
        payload["parent_context"]["variables"]["U"]["parent_end_vs_parent_start"][
            "rmse"
        ]
        == 2.0
    )
    assert payload["diagnostic_only"] is True
    assert payload["gate_evidence"] is False
    assert payload["advances_00_20"] is False


def test_help_smoke() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "tools/audit_selected_field_wind_source_attribution.py",
            "--help",
        ],
        check=False,
        cwd=Path(__file__).resolve().parents[2],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "d02 start-state wrfout file" in result.stdout
