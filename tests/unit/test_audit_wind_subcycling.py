import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.audit_wind_subcycling import (
    audit_wind_subcycling,
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
    if data.ndim == 2:
        return ("Time", f"{name}_south_north", f"{name}_west_east")
    if data.ndim == 3:
        return ("Time", f"{name}_bottom_top", f"{name}_south_north", f"{name}_west_east")
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
                _ensure_dim(dataset, dim_name, int(size))
            variable = dataset.createVariable(name, "f8", dims)
            variable[0, ...] = values


def _wind_attrs() -> dict[str, object]:
    return {
        "TYWRF_WIND_TENDENCY_SUBSTEP_COUNT": 75,
        "TYWRF_WIND_TENDENCY_SUBSTEP_DT_SECONDS": 8.0,
        "TYWRF_WIND_TENDENCY_TOTAL_SECONDS": 600.0,
        "TYWRF_WIND_TENDENCY_ADVECTING_VELOCITY_MODE": "refreshed",
        "TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS": "cross_component",
        "TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION": "average",
        "TYWRF_WIND_TENDENCY_ADVECTION_FORM": "conservative",
    }


def _variables_by_name(payload: dict) -> dict[str, dict]:
    variables = payload["candidates"][0]["variables"]
    return {entry["variable"]: entry for entry in variables}


def test_function_report_contains_global_boundary_vertical_delta_and_metadata(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.nc"
    baseline = tmp_path / "baseline.nc"
    candidate = tmp_path / "candidate.nc"

    ref_u = np.full((2, 3, 4), 10.0)
    ref_v = np.full((2, 4, 3), 20.0)
    cand_u = ref_u + 1.0
    cand_v = ref_v.copy()
    cand_v[:, 0, :] -= 2.0

    _write_wrfout(reference, {"U": ref_u, "V": ref_v})
    _write_wrfout(baseline, {"U": ref_u, "V": ref_v})
    _write_wrfout(candidate, {"U": cand_u, "V": cand_v}, attrs=_wind_attrs())

    payload = json.loads(
        report_to_json(
            audit_wind_subcycling(
                reference,
                {"self_advection": candidate},
                baseline_wrfout=baseline,
                boundary_bands=(0, 1),
            ),
            pretty=True,
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["advances_00_20"] is False
    assert payload["uses_reference_end_truth"] is False
    assert payload["summary"]["candidate_count"] == 1

    metadata = payload["candidates"][0]["metadata"]["global_attrs"]
    assert metadata["TYWRF_WIND_TENDENCY_SUBSTEP_COUNT"] == 75
    assert metadata["TYWRF_WIND_TENDENCY_SUBSTEP_DT_SECONDS"] == 8.0
    assert metadata["TYWRF_WIND_TENDENCY_TOTAL_SECONDS"] == 600.0
    assert metadata["TYWRF_WIND_TENDENCY_ADVECTING_VELOCITY_MODE"] == "refreshed"
    assert metadata["TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS"] == "cross_component"
    assert metadata["TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION"] == "average"
    assert metadata["TYWRF_WIND_TENDENCY_ADVECTION_FORM"] == "conservative"

    by_name = _variables_by_name(payload)
    u = by_name["U"]
    assert u["status"] == "computed"
    assert u["reference_shape"] == [2, 3, 4]
    assert math.isclose(u["global"]["rmse"], 1.0)
    assert math.isclose(u["global"]["normalized_rmse"], 0.1)
    assert u["global"]["max_abs_error"] == 1.0
    assert u["global"]["valid_count"] == 24

    assert u["delta_vs_baseline"]["status"] == "computed"
    assert u["delta_vs_baseline"]["changed_count"] == 24
    assert math.isclose(u["delta_vs_baseline"]["rmse"], 1.0)
    assert math.isclose(u["delta_vs_baseline"]["normalized_rmse"], 0.1)
    assert u["delta_vs_baseline"]["max_abs_diff"] == 1.0
    assert u["delta_vs_baseline"]["differing_count"] == 24
    assert u["delta_vs_baseline"]["max_abs_delta"] == 1.0

    band0 = u["boundary_bands"]["bands"][0]
    assert band0["band"] == 0
    assert band0["boundary"]["count"] == 20
    assert band0["interior"]["count"] == 4
    assert math.isclose(band0["boundary"]["rmse"], 1.0)
    assert math.isclose(band0["interior"]["rmse"], 1.0)

    assert u["vertical_levels"]["status"] == "available"
    assert u["vertical_levels"]["level_count"] == 2
    assert [level["k"] for level in u["vertical_levels"]["levels"]] == [0, 1]
    assert all(
        math.isclose(level["metrics"]["rmse"], 1.0)
        for level in u["vertical_levels"]["levels"]
    )

    v = by_name["V"]
    assert v["status"] == "computed"
    assert v["reference_shape"] == [2, 4, 3]
    assert v["delta_vs_baseline"]["changed_count"] == 6
    assert math.isclose(v["delta_vs_baseline"]["rmse"], 1.0)
    assert math.isclose(v["delta_vs_baseline"]["normalized_rmse"], 0.05)
    assert v["delta_vs_baseline"]["max_abs_diff"] == 2.0
    assert v["delta_vs_baseline"]["differing_count"] == 6
    assert v["boundary_bands"]["horizontal_shape"] == [4, 3]


def test_same_candidate_and_baseline_path_reports_zero_delta(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"

    ref_u = np.full((1, 2, 3), 3.0)
    ref_v = np.full((1, 3, 2), 4.0)
    cand_u = ref_u + 2.0
    cand_v = ref_v - 1.0

    _write_wrfout(reference, {"U": ref_u, "V": ref_v})
    _write_wrfout(candidate, {"U": cand_u, "V": cand_v}, attrs=_wind_attrs())

    payload = audit_wind_subcycling(
        reference,
        {"same": candidate},
        baseline_wrfout=candidate,
    )
    by_name = _variables_by_name(payload)

    for variable in ("U", "V"):
        delta = by_name[variable]["delta_vs_baseline"]
        assert delta["status"] == "computed"
        assert delta["same_path_as_baseline"] is True
        assert delta["changed_count"] == 0
        assert delta["rmse"] == 0.0
        assert delta["normalized_rmse"] == 0.0
        assert delta["max_abs_diff"] == 0.0
        assert delta["differing_count"] == 0
        assert delta["max_abs_delta"] == 0.0


def test_baseline_only_nonwind_is_ignored_and_staggered_uv_shapes_compute(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.nc"
    baseline = tmp_path / "baseline.nc"
    candidate = tmp_path / "candidate.nc"

    ref_u = np.full((2, 3), 4.0)
    ref_v = np.full((3, 2), 5.0)
    cand_u = ref_u.copy()
    cand_v = ref_v + 1.0

    _write_wrfout(reference, {"U": ref_u, "V": ref_v, "T": np.ones((2, 2))})
    _write_wrfout(baseline, {"T": np.full((2, 2), 9.0)})
    _write_wrfout(candidate, {"U": cand_u, "V": cand_v})

    payload = audit_wind_subcycling(
        reference,
        [("candidate", candidate)],
        baseline_wrfout=baseline,
    )
    by_name = _variables_by_name(payload)

    assert sorted(by_name) == ["U", "V"]
    assert by_name["U"]["status"] == "computed"
    assert by_name["V"]["status"] == "computed"
    assert by_name["U"]["reference_shape"] == [2, 3]
    assert by_name["V"]["reference_shape"] == [3, 2]
    assert by_name["U"]["vertical_levels"]["status"] == "no_vertical_dimension"
    assert by_name["V"]["vertical_levels"]["status"] == "no_vertical_dimension"
    assert by_name["U"]["delta_vs_baseline"]["status"] == "missing_baseline_variable"
    assert by_name["V"]["delta_vs_baseline"]["status"] == "missing_baseline_variable"
    assert by_name["U"]["delta_vs_baseline"]["changed_count"] is None
    assert by_name["V"]["delta_vs_baseline"]["max_abs_diff"] is None
    assert by_name["U"]["delta_vs_baseline"]["differing_count"] is None
    assert by_name["V"]["delta_vs_baseline"]["max_abs_delta"] is None


def test_main_writes_pretty_json_output_path(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "audit.json"

    _write_wrfout(
        reference,
        {
            "U": np.ones((1, 2, 3)),
            "V": np.ones((1, 3, 2)),
        },
    )
    _write_wrfout(
        candidate,
        {
            "U": np.ones((1, 2, 3)) + 0.5,
            "V": np.ones((1, 3, 2)),
        },
        attrs=_wind_attrs(),
    )

    exit_code = audit_main(
        [
            "--reference-wrfout",
            str(reference),
            "--candidate",
            f"subcycled={candidate}",
            "--variables",
            "U,V",
            "--boundary-bands",
            "0,1",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    text = output.read_text(encoding="utf-8")
    payload = json.loads(text)

    assert exit_code == 0
    assert text.startswith("{\n  ")
    assert payload["reference_wrfout"] == str(reference)
    assert payload["candidates"][0]["label"] == "subcycled"
    assert payload["variables"] == ["U", "V"]
    assert payload["boundary_bands"] == [0, 1]
    assert payload["candidates"][0]["metadata"]["global_attrs"][
        "TYWRF_WIND_TENDENCY_SUBSTEP_COUNT"
    ] == 75
    assert payload["candidates"][0]["metadata"]["global_attrs"][
        "TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS"
    ] == "cross_component"
