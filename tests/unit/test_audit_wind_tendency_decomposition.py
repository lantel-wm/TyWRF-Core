import json
import math
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.audit_wind_tendency_decomposition import (
    audit_wind_tendency_decomposition,
    main as audit_main,
    report_to_json,
    resolve_child_delta,
)


def _ensure_dim(dataset: netCDF4.Dataset, name: str, size: int) -> None:
    if name in dataset.dimensions:
        assert len(dataset.dimensions[name]) == size
        return
    dataset.createDimension(name, size)


def _variable_dims(name: str, data: np.ndarray) -> tuple[str, ...]:
    if name == "U" and data.ndim == 3:
        return ("Time", "bottom_top", "south_north", "west_east_stag")
    if name == "V" and data.ndim == 3:
        return ("Time", "bottom_top", "south_north_stag", "west_east")
    if data.ndim == 2:
        return ("Time", f"{name}_south_north", f"{name}_west_east")
    if data.ndim == 3:
        return ("Time", f"{name}_bottom_top", f"{name}_south_north", f"{name}_west_east")
    return ("Time", *(f"{name}_dim_{axis}" for axis in range(data.ndim)))


def _write_wrfout(path: Path, variables: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        for name, raw_values in variables.items():
            values = np.asarray(raw_values, dtype=np.float64)
            dims = _variable_dims(name, values)
            for dim_name, size in zip(dims[1:], values.shape, strict=True):
                _ensure_dim(dataset, dim_name, int(size))
            variable = dataset.createVariable(name, "f8", dims)
            variable[0, ...] = values


def _result_by_name(payload: dict) -> dict[str, dict]:
    return {entry["variable"]: entry for entry in payload["results"]}


def _required_cli_args(tmp_path: Path) -> list[str]:
    return [
        "--reference-start-wrfout",
        str(tmp_path / "reference_start.nc"),
        "--reference-end-wrfout",
        str(tmp_path / "reference_end.nc"),
        "--candidate-start-wrfout",
        str(tmp_path / "candidate_start.nc"),
        "--candidate-end-wrfout",
        str(tmp_path / "candidate_end.nc"),
    ]


def _assert_cli_usage_error(args: list[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        audit_main(args)
    assert exc_info.value.code == 2


def test_global_residual_boundary_and_vertical_metrics(tmp_path: Path) -> None:
    reference_start = tmp_path / "reference_start.nc"
    reference_end = tmp_path / "reference_end.nc"
    candidate_start = tmp_path / "candidate_start.nc"
    candidate_end = tmp_path / "candidate_end.nc"

    shape = (2, 3, 4)
    _write_wrfout(reference_start, {"U": np.full(shape, 2.0)})
    _write_wrfout(reference_end, {"U": np.full(shape, 12.0)})
    _write_wrfout(candidate_start, {"U": np.full(shape, 5.0)})
    _write_wrfout(candidate_end, {"U": np.full(shape, 16.0)})

    payload = json.loads(
        report_to_json(
            audit_wind_tendency_decomposition(
                reference_start,
                reference_end,
                candidate_start,
                candidate_end,
                variables=("U",),
                boundary_bands=(0, 1),
            ),
            pretty=True,
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["uses_reference_end_truth"] is True
    assert payload["advances_00_20"] is False
    assert payload["increment_seconds"] == 600.0

    u = _result_by_name(payload)["U"]
    assert u["status"] == "computed"
    assert u["reference_increment_shape"] == [2, 3, 4]
    assert u["candidate_increment_shape"] == [2, 3, 4]
    assert math.isclose(u["global"]["rmse"], 1.0)
    assert math.isclose(u["global"]["normalized_rmse"], 0.1)
    assert math.isclose(u["global"]["max_abs_error"], 1.0)
    assert u["global"]["valid_count"] == 24
    assert math.isclose(u["global"]["sum_squared_residual"], 24.0)
    assert math.isclose(u["reference_increment_stats"]["rms"], 10.0)
    assert math.isclose(u["candidate_increment_stats"]["rms"], 11.0)

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


def test_child_delta_region_split_is_deterministic(tmp_path: Path) -> None:
    reference_start = tmp_path / "reference_start.nc"
    reference_end = tmp_path / "reference_end.nc"
    candidate_start = tmp_path / "candidate_start.nc"
    candidate_end = tmp_path / "candidate_end.nc"

    shape = (2, 4, 5)
    _write_wrfout(reference_start, {"U": np.zeros(shape)})
    _write_wrfout(reference_end, {"U": np.full(shape, 10.0)})
    _write_wrfout(candidate_start, {"U": np.zeros(shape)})
    _write_wrfout(candidate_end, {"U": np.full(shape, 12.0)})

    payload = audit_wind_tendency_decomposition(
        reference_start,
        reference_end,
        candidate_start,
        candidate_end,
        variables=("U",),
        boundary_bands=(0,),
        child_delta=(2, 1),
        child_delta_source="child_delta",
    )
    regions = _result_by_name(payload)["U"]["regions"]
    by_name = regions["regions_by_name"]

    assert regions["status"] == "available"
    assert regions["horizontal_shape"] == [4, 5]
    assert regions["child_delta"] == [2, 1]
    assert regions["child_delta_units"] == "child_grid_cells"
    assert regions["core_boundary_band"] == 0
    assert by_name["overlap"]["metrics"]["selected_count"] == 18
    assert by_name["core"]["metrics"]["selected_count"] == 8
    assert by_name["east_exposed"]["metrics"]["selected_count"] == 12
    assert by_name["north_exposed"]["metrics"]["selected_count"] == 6
    assert by_name["corner"]["metrics"]["selected_count"] == 4
    assert math.isclose(
        by_name["overlap"]["sum_squared_residual_fraction_global"],
        0.45,
    )
    assert math.isclose(
        by_name["east_exposed"]["sum_squared_residual_fraction_global"],
        0.3,
    )
    assert math.isclose(
        by_name["north_exposed"]["sum_squared_residual_fraction_global"],
        0.15,
    )
    assert math.isclose(
        by_name["corner"]["sum_squared_residual_fraction_global"],
        0.1,
    )


def test_parent_start_delta_is_reported_as_unscaled() -> None:
    child_delta, source = resolve_child_delta(
        from_parent_start=(114, 96),
        to_parent_start=(126, 103),
    )

    assert child_delta == (12, 7)
    assert source == "parent_start_difference_unscaled"


def test_main_writes_pretty_json_output_path_with_default_uv(tmp_path: Path) -> None:
    reference_start = tmp_path / "reference_start.nc"
    reference_end = tmp_path / "reference_end.nc"
    candidate_start = tmp_path / "candidate_start.nc"
    candidate_end = tmp_path / "candidate_end.nc"
    output = tmp_path / "decomposition.json"

    _write_wrfout(
        reference_start,
        {
            "U": np.ones((1, 2, 3)),
            "V": np.ones((1, 3, 2)),
        },
    )
    _write_wrfout(
        reference_end,
        {
            "U": np.ones((1, 2, 3)) * 3.0,
            "V": np.ones((1, 3, 2)) * 4.0,
        },
    )
    _write_wrfout(
        candidate_start,
        {
            "U": np.ones((1, 2, 3)) * 5.0,
            "V": np.ones((1, 3, 2)) * 6.0,
        },
    )
    _write_wrfout(
        candidate_end,
        {
            "U": np.ones((1, 2, 3)) * 8.0,
            "V": np.ones((1, 3, 2)) * 10.0,
        },
    )

    exit_code = audit_main(
        [
            "--reference-start-wrfout",
            str(reference_start),
            "--reference-end-wrfout",
            str(reference_end),
            "--candidate-start-wrfout",
            str(candidate_start),
            "--candidate-end-wrfout",
            str(candidate_end),
            "--child-delta",
            "1",
            "0",
            "--boundary-bands",
            "0,1",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    text = output.read_text(encoding="utf-8")
    payload = json.loads(text)
    by_name = _result_by_name(payload)

    assert exit_code == 0
    assert text.startswith("{\n  ")
    assert payload["diagnostic_only"] is True
    assert payload["uses_reference_end_truth"] is True
    assert payload["advances_00_20"] is False
    assert payload["variables"] == ["U", "V"]
    assert payload["boundary_bands"] == [0, 1]
    assert payload["child_delta"] == [1, 0]
    assert payload["child_delta_source"] == "child_delta"
    assert math.isclose(by_name["U"]["global"]["rmse"], 1.0)
    assert math.isclose(by_name["V"]["global"]["rmse"], 1.0)


def test_missing_requested_variable_fails_closed(tmp_path: Path) -> None:
    reference_start = tmp_path / "reference_start.nc"
    reference_end = tmp_path / "reference_end.nc"
    candidate_start = tmp_path / "candidate_start.nc"
    candidate_end = tmp_path / "candidate_end.nc"

    shape = (2, 3, 4)
    _write_wrfout(reference_start, {})
    _write_wrfout(reference_end, {"U": np.ones(shape)})
    _write_wrfout(candidate_start, {"U": np.ones(shape)})
    _write_wrfout(candidate_end, {"U": np.ones(shape)})

    payload = audit_wind_tendency_decomposition(
        reference_start,
        reference_end,
        candidate_start,
        candidate_end,
        variables=("U",),
        boundary_bands=(0,),
    )
    u = _result_by_name(payload)["U"]

    assert payload["status"] == "computed_with_flags"
    assert payload["summary"]["flagged"] == 1
    assert payload["summary"]["computed"] == 0
    assert u["status"] == "missing_variable"
    assert u["source_shapes"]["reference_start"] is None
    assert u["global"] is None
    assert u["reference_increment_stats"] is None
    assert u["candidate_increment_stats"] is None
    assert u["boundary_bands"]["status"] == "not_computed"
    assert u["vertical_levels"]["status"] == "not_computed"
    assert u["regions"]["status"] == "not_computed"


@pytest.mark.parametrize(
    ("reference_start_shape", "reference_end_shape", "candidate_shape"),
    [
        ((2, 3, 4), (2, 3, 5), (2, 3, 4)),
        ((2, 3, 4), (2, 3, 4), (2, 3, 5)),
    ],
)
def test_shape_mismatch_fails_closed_without_global_metrics(
    tmp_path: Path,
    reference_start_shape: tuple[int, int, int],
    reference_end_shape: tuple[int, int, int],
    candidate_shape: tuple[int, int, int],
) -> None:
    reference_start = tmp_path / "reference_start.nc"
    reference_end = tmp_path / "reference_end.nc"
    candidate_start = tmp_path / "candidate_start.nc"
    candidate_end = tmp_path / "candidate_end.nc"

    _write_wrfout(reference_start, {"U": np.zeros(reference_start_shape)})
    _write_wrfout(reference_end, {"U": np.ones(reference_end_shape)})
    _write_wrfout(candidate_start, {"U": np.zeros(candidate_shape)})
    _write_wrfout(candidate_end, {"U": np.ones(candidate_shape)})

    payload = audit_wind_tendency_decomposition(
        reference_start,
        reference_end,
        candidate_start,
        candidate_end,
        variables=("U",),
        boundary_bands=(0,),
    )
    u = _result_by_name(payload)["U"]

    assert payload["status"] == "computed_with_flags"
    assert payload["summary"]["flagged"] == 1
    assert u["status"] == "shape_mismatch"
    assert u["global"] is None
    assert u["boundary_bands"]["status"] == "not_computed"
    assert u["vertical_levels"]["status"] == "not_computed"
    assert u["regions"]["status"] == "not_computed"


def test_cli_rejects_child_delta_with_parent_start_pair(tmp_path: Path) -> None:
    _assert_cli_usage_error(
        [
            *_required_cli_args(tmp_path),
            "--child-delta",
            "1",
            "0",
            "--from-parent-start",
            "114",
            "96",
            "--to-parent-start",
            "126",
            "103",
        ]
    )


@pytest.mark.parametrize(
    "parent_args",
    [
        ["--from-parent-start", "114", "96"],
        ["--to-parent-start", "126", "103"],
    ],
)
def test_cli_rejects_one_sided_parent_start(
    tmp_path: Path,
    parent_args: list[str],
) -> None:
    _assert_cli_usage_error([*_required_cli_args(tmp_path), *parent_args])


@pytest.mark.parametrize("boundary_bands", ["-1", "0,a"])
def test_cli_rejects_invalid_boundary_bands(
    tmp_path: Path,
    boundary_bands: str,
) -> None:
    _assert_cli_usage_error(
        [
            *_required_cli_args(tmp_path),
            "--boundary-bands",
            boundary_bands,
        ]
    )


def test_cli_rejects_empty_variables(tmp_path: Path) -> None:
    _assert_cli_usage_error(
        [
            *_required_cli_args(tmp_path),
            "--variables",
            ",",
        ]
    )
