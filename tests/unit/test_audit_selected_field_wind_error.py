import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.audit_selected_field_wind_error import (
    audit_selected_field_wind_error,
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


def test_staggered_uv_metrics_include_max_location_and_values(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    ref_u = np.full((1, 2, 4), 10.0)
    ref_v = np.full((1, 3, 3), 20.0)
    cand_u = ref_u.copy()
    cand_v = ref_v.copy()
    cand_u[0, 1, 3] = 16.0
    cand_v[0, 2, 0] = 14.0

    _write_wrfout(reference, {"U": ref_u, "V": ref_v, "MU": np.ones((2, 3))})
    _write_wrfout(
        candidate,
        {"U": cand_u, "V": cand_v, "MU": np.ones((2, 3))},
        attrs=_base_attrs(),
    )

    payload = json.loads(
        report_to_json(audit_selected_field_wind_error(reference, candidate))
    )

    assert payload["diagnostic_only"] is True
    assert payload["gate_evidence"] is False
    assert payload["advances_00_20"] is False

    by_variable = {entry["variable"]: entry for entry in payload["variables"]}
    u_metrics = by_variable["U"]["whole_domain"]
    v_metrics = by_variable["V"]["whole_domain"]
    assert math.isclose(u_metrics["rmse"], math.sqrt(36.0 / 8.0))
    assert math.isclose(u_metrics["normalized_rmse"], math.sqrt(36.0 / 8.0) / 10.0)
    assert u_metrics["max_abs"] == 6.0
    assert u_metrics["valid_count"] == 8
    assert u_metrics["max_abs_location_index"] == [0, 1, 3]
    assert u_metrics["reference_value_at_max_abs"] == 10.0
    assert u_metrics["candidate_value_at_max_abs"] == 16.0

    assert math.isclose(v_metrics["rmse"], math.sqrt(36.0 / 9.0))
    assert math.isclose(v_metrics["normalized_rmse"], math.sqrt(36.0 / 9.0) / 20.0)
    assert v_metrics["max_abs_location_index"] == [0, 2, 0]
    assert v_metrics["reference_value_at_max_abs"] == 20.0
    assert v_metrics["candidate_value_at_max_abs"] == 14.0


def test_region_breakdown_handles_uv_staggered_exposed_and_overlap(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    ref_u = np.ones((1, 3, 5))
    ref_v = np.ones((1, 4, 4))
    cand_u = ref_u.copy()
    cand_v = ref_v.copy()
    cand_u[:, :, -2:] += 2.0
    cand_u[:, -1:, :] += 2.0
    cand_v[:, :, -1:] += 3.0
    cand_v[:, -2:, :] += 3.0

    _write_wrfout(reference, {"U": ref_u, "V": ref_v, "MU": np.ones((3, 4))})
    _write_wrfout(
        candidate,
        {"U": cand_u, "V": cand_v, "MU": np.ones((3, 4))},
        attrs=_base_attrs(from_start="10,20", to_start="11,21"),
    )

    payload = json.loads(
        report_to_json(audit_selected_field_wind_error(reference, candidate))
    )

    by_variable = {entry["variable"]: entry for entry in payload["variables"]}
    u_breakdown = by_variable["U"]["region_breakdown"]
    v_breakdown = by_variable["V"]["region_breakdown"]

    assert u_breakdown["region_breakdown_status"] == "available"
    assert u_breakdown["shape_kind"] == "x_staggered"
    assert u_breakdown["exposed_region"]["valid_count"] == 9
    assert u_breakdown["overlap_region"]["valid_count"] == 6
    assert u_breakdown["exposed_region_dominates_global_error"] is True

    assert v_breakdown["shape_kind"] == "y_staggered"
    assert v_breakdown["exposed_region"]["valid_count"] == 10
    assert v_breakdown["overlap_region"]["valid_count"] == 6
    assert payload["summary"]["region_breakdown_status"] == {
        "U": "available",
        "V": "available",
    }


def test_candidate_metadata_reports_guard_fields_and_timeline(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    variables = {
        "U": np.ones((1, 2, 4)),
        "V": np.ones((1, 3, 3)),
        "MU": np.ones((2, 3)),
    }

    _write_wrfout(reference, variables)
    _write_wrfout(candidate, variables, attrs=_base_attrs())

    payload = json.loads(
        report_to_json(audit_selected_field_wind_error(reference, candidate))
    )
    metadata = payload["metadata"]

    assert metadata["gate_evidence"] is False
    assert metadata["advances_00_20"] is False
    assert metadata["parsed"]["dx"] == 2000.0
    assert metadata["parsed"]["dy"] == 2000.0
    assert metadata["parsed"]["d02_resolution_check"] == "d02_2km"
    assert metadata["parsed"]["from_parent_start"] == {"i": 10, "j": 20}
    assert metadata["parsed"]["to_parent_start"] == {"i": 11, "j": 21}
    assert metadata["parsed"]["remap_from_parent_start"] == {"i": 10, "j": 20}
    assert metadata["parsed"]["remap_to_parent_start"] == {"i": 11, "j": 21}
    assert metadata["parsed"]["timeline_event_count"] == 2
    assert metadata["parsed"]["timeline_event_names"] == [
        "cycle_start",
        "move_from_to_parent_start",
    ]
    assert metadata["parsed"]["timeline_events"][1]["name"] == (
        "move_from_to_parent_start"
    )
    assert metadata["parsed"]["timeline_events"][1]["fields"] == {
        "from": "10_20",
        "to": "11_21",
    }


def test_insufficient_metadata_region_breakdown_fails_closed(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    ref_u = np.full((1, 2, 4), 10.0)
    cand_u = ref_u.copy()
    cand_u[0, 0, 0] = 13.0

    _write_wrfout(reference, {"U": ref_u, "V": np.ones((1, 3, 3))})
    _write_wrfout(candidate, {"U": cand_u, "V": np.ones((1, 3, 3))})

    payload = json.loads(
        report_to_json(
            audit_selected_field_wind_error(reference, candidate, variables=("U",))
        )
    )

    result = payload["variables"][0]
    assert result["status"] == "computed"
    assert result["whole_domain"]["max_abs_location_index"] == [0, 0, 0]
    assert result["region_breakdown_status"] == "insufficient_metadata"
    assert result["region_breakdown"]["region_breakdown_status"] == (
        "insufficient_metadata"
    )
    assert result["region_breakdown"]["whole_domain"]["valid_count"] == 8
    assert "TYWRF_FROM_PARENT_START" in payload["movement_region"]["missing_metadata"]
    assert payload["status"] == "computed_with_flags"


def test_main_writes_cli_json_output(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "wind_error.json"
    variables = {
        "U": np.ones((1, 2, 4)),
        "V": np.ones((1, 3, 3)),
        "MU": np.ones((2, 3)),
    }

    _write_wrfout(reference, variables)
    _write_wrfout(candidate, variables, attrs=_base_attrs())

    exit_code = audit_main(
        [
            str(reference),
            str(candidate),
            "--domain",
            "d02",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["domain"] == "d02"
    assert payload["diagnostic_only"] is True
    assert payload["gate_evidence"] is False
    assert payload["advances_00_20"] is False
    assert [entry["variable"] for entry in payload["variables"]] == ["U", "V"]
