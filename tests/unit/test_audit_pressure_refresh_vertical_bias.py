import json
import math
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.audit_pressure_refresh_vertical_bias import (
    audit_pressure_refresh_vertical_bias,
    main as audit_main,
    report_to_json,
)


def _ensure_dim(dataset: netCDF4.Dataset, name: str, size: int) -> None:
    if name not in dataset.dimensions:
        dataset.createDimension(name, size)


def _variable_dims(name: str, data: np.ndarray) -> tuple[str, ...]:
    if data.ndim == 3:
        vertical = "bottom_top_stag" if name in {"PH", "PHB"} else "bottom_top"
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
        "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT": 6,
    }


def _by_level(payload: dict[str, object]) -> dict[int, dict[str, object]]:
    return {entry["level"]: entry for entry in payload["per_level_p"]}


def _attribution_by_level(payload: dict[str, object]) -> dict[int, dict[str, object]]:
    return {entry["level"]: entry for entry in payload["worst_level_attribution"]}


def _field_by_name(attribution: dict[str, object]) -> dict[str, dict[str, object]]:
    return {entry["field"]: entry for entry in attribution["companion_fields"]}


def test_worst_level_ranking_and_companion_field_means(tmp_path: Path) -> None:
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate_end.nc"
    ref_p = np.full((3, 2, 4), 100.0)
    cand_p = ref_p.copy()
    cand_p[0, :, 3] += 5.0
    cand_p[1, :, 3] -= 9.0
    cand_p[2, 0, 3] += 7.0
    cand_p[2, 1, 3] -= 7.0

    ref_pb = np.full((3, 2, 4), 900.0)
    cand_pb = ref_pb.copy()
    cand_pb[1, :, 3] += 2.0
    ref_mu = np.full((2, 4), 50.0)
    cand_mu = ref_mu.copy()
    cand_mu[:, 3] += 3.0
    ref_ph = np.stack([np.full((2, 4), value) for value in (0.0, 10.0, 20.0, 30.0)])
    cand_ph = ref_ph + 1.0

    _write_wrfout(reference, {"P": ref_p, "PB": ref_pb, "MU": ref_mu, "PH": ref_ph})
    _write_wrfout(
        candidate,
        {"P": cand_p, "PB": cand_pb, "MU": cand_mu, "PH": cand_ph},
        attrs=_base_attrs(),
    )

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_vertical_bias(
                reference,
                candidate,
                worst_level_count=2,
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"
    assert payload["region"]["status"] == "available"
    assert payload["worst_levels_by_rmse"] == [1, 2]
    assert payload["worst_levels_by_abs_bias"] == [1, 0]

    per_level = _by_level(payload)
    assert math.isclose(per_level[1]["mean_candidate_p"], 91.0)
    assert math.isclose(per_level[1]["mean_reference_p"], 100.0)
    assert math.isclose(per_level[1]["mean_diff"], -9.0)
    assert math.isclose(per_level[1]["rmse"], 9.0)
    assert math.isclose(per_level[1]["normalized_rmse"], 0.09)
    assert per_level[1]["valid_count"] == 2

    attribution = _attribution_by_level(payload)
    assert attribution[1]["rank_by_rmse"] == 1
    assert attribution[1]["rank_by_abs_bias"] == 1
    fields = _field_by_name(attribution[1])
    assert math.isclose(fields["PB"]["mean_diff"], 2.0)
    assert math.isclose(fields["P+PB"]["mean_diff"], -7.0)
    assert math.isclose(fields["MU"]["mean_diff"], 3.0)
    assert fields["PH"]["alignment"] == "staggered_vertical_adjacent_mean"
    assert math.isclose(fields["PH"]["mean_diff"], 1.0)


def test_missing_optional_companion_fields_are_not_available(tmp_path: Path) -> None:
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate_end.nc"
    p = np.full((1, 2, 4), 100.0)
    cand_p = p.copy()
    cand_p[:, :, 3] += 1.0

    _write_wrfout(reference, {"P": p})
    _write_wrfout(candidate, {"P": cand_p}, attrs=_base_attrs())

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_vertical_bias(
                reference,
                candidate,
                worst_level_count=1,
            )
        )
    )
    fields = _field_by_name(payload["worst_level_attribution"][0])

    assert fields["PB"]["status"] == "not_available"
    assert fields["P+PB"]["status"] == "not_available"
    assert fields["MU"]["status"] == "not_available"
    assert fields["MUB"]["status"] == "not_available"
    assert fields["PH"]["status"] == "not_available"
    assert fields["PHB"]["status"] == "not_available"
    assert fields["HGT"]["status"] == "not_available"
    assert fields["PSFC"]["status"] == "not_available"
    assert fields["PB"]["candidate_model_pass"] == "not_applicable"


def test_source_start_evolution_summary(tmp_path: Path) -> None:
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate_end.nc"
    source = tmp_path / "source_start.nc"
    source_p = np.full((2, 2, 4), 100.0)
    ref_p = source_p.copy()
    cand_p = source_p.copy()
    ref_p[0, :, 3] += 4.0
    cand_p[0, :, 3] -= 6.0
    ref_p[1, :, 3] += 1.0
    cand_p[1, :, 3] += 2.0

    _write_wrfout(reference, {"P": ref_p})
    _write_wrfout(candidate, {"P": cand_p}, attrs=_base_attrs())
    _write_wrfout(source, {"P": source_p})

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_vertical_bias(
                reference,
                candidate,
                source,
                worst_level_count=1,
            )
        )
    )

    assert payload["worst_levels_by_rmse"] == [0]
    evolution = payload["worst_level_attribution"][0]["source_evolution"]
    assert evolution["status"] == "computed"
    assert evolution["diagnostic_only"] is True
    assert evolution["candidate_model_pass"] == "not_applicable"
    assert math.isclose(evolution["mean_source_p"], 100.0)
    assert math.isclose(evolution["mean_candidate_minus_source_p"], -6.0)
    assert math.isclose(evolution["mean_reference_minus_source_p"], 4.0)
    assert math.isclose(
        evolution["mean_candidate_evolution_minus_reference_evolution_p"],
        -10.0,
    )
    assert evolution["valid_count"] == 2


def test_cli_json_output_and_help_invocation(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate_end.nc"
    output = tmp_path / "vertical_bias.json"
    p = np.full((1, 2, 4), 100.0)
    cand_p = p.copy()
    cand_p[:, :, 3] += 1.0
    _write_wrfout(reference, {"P": p})
    _write_wrfout(candidate, {"P": cand_p}, attrs=_base_attrs())

    exit_code = audit_main(
        [
            "--reference-end",
            str(reference),
            "--candidate-end",
            str(candidate),
            "--output",
            str(output),
            "--pretty",
            "--worst-level-count",
            "1",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["worst_level_by_rmse"] == 0
    assert payload["per_level_p"][0]["mean_diff"] == 1.0

    with pytest.raises(SystemExit) as excinfo:
        audit_main(["--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--reference-end" in captured.out
    assert "--candidate-end" in captured.out
