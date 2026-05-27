import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.compare_wrfout import STRICT_CORE_VARIABLES, compare_files, main as compare_main
from tools.compare_wrfout import normalized_rmse, report_to_json


def _write_dataset(
    path: Path,
    variables: dict[str, np.ndarray],
    *,
    attrs: dict[str, str] | None = None,
) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        for name, value in (attrs or {}).items():
            dataset.setncattr(name, value)
        for name, values in variables.items():
            data = np.asarray(values, dtype=np.float64)
            dimensions = []
            for axis, size in enumerate(data.shape):
                dimension = f"{name}_dim_{axis}"
                dataset.createDimension(dimension, size)
                dimensions.append(dimension)
            variable = dataset.createVariable(name, "f8", tuple(dimensions))
            variable[:] = data


def _strict_core_fields(value: float = 1.0) -> dict[str, np.ndarray]:
    return {
        name: np.array([value, value + 1.0, value + 2.0])
        for name in STRICT_CORE_VARIABLES
    }


def _write_dataset_with_tc_fields(path: Path) -> None:
    shape = (3, 4)
    latitude = np.array(
        [
            [10.0, 10.0, 10.0, 10.0],
            [11.0, 11.0, 11.0, 11.0],
            [12.0, 12.0, 12.0, 12.0],
        ],
        dtype=np.float64,
    )
    longitude = np.array(
        [
            [120.0, 121.0, 122.0, 123.0],
            [120.0, 121.0, 122.0, 123.0],
            [120.0, 121.0, 122.0, 123.0],
        ],
        dtype=np.float64,
    )
    psfc = np.full(shape, 100000.0)
    psfc[1, 2] = 95000.0
    slp = np.full(shape, 1000.0)
    slp[1, 2] = 950.0
    u10 = np.zeros(shape)
    v10 = np.zeros(shape)
    u10[0, 1] = 30.0
    v10[0, 1] = 40.0
    rainc = np.arange(12, dtype=np.float64).reshape(shape)
    rainnc = np.full(shape, 0.5)

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", 1)
        dataset.createDimension("south_north", shape[0])
        dataset.createDimension("west_east", shape[1])
        for name, values in {
            "U": np.ones(shape),
            "XLAT": latitude,
            "XLONG": longitude,
            "SLP": slp,
            "PSFC": psfc,
            "U10": u10,
            "V10": v10,
            "RAINC": rainc,
            "RAINNC": rainnc,
        }.items():
            variable = dataset.createVariable(name, "f8", ("Time", "south_north", "west_east"))
            if name == "SLP":
                variable.units = "hPa"
            variable[0, :, :] = values


def test_normalized_rmse_identical_arrays() -> None:
    reference = np.array([1.0, 2.0, 3.0])
    rmse, norm, max_abs = normalized_rmse(reference, reference.copy())

    assert rmse == 0.0
    assert norm == 0.0
    assert max_abs == 0.0


def test_normalized_rmse_nonzero_error() -> None:
    reference = np.array([1.0, 2.0, 3.0])
    candidate = np.array([2.0, 2.0, 4.0])
    rmse, norm, max_abs = normalized_rmse(reference, candidate)

    assert math.isclose(rmse, math.sqrt(2.0 / 3.0))
    assert math.isclose(norm, math.sqrt(2.0 / 14.0))
    assert max_abs == 1.0


def test_compare_files_reports_metrics_and_pending_tc_diagnostics(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    _write_dataset(
        reference,
        {
            "U": np.array([1.0, 2.0, 3.0]),
            "RAINC": np.array([0.0, 2.0, 4.0]),
        },
    )
    _write_dataset(
        candidate,
        {
            "U": np.array([1.0, 2.1, 2.9]),
            "RAINC": np.array([0.0, 2.5, 3.5]),
        },
    )

    report = compare_files(reference, candidate, variables=("U", "RAINC"))

    assert report.status == "ok"
    assert report.summary == {"ok": 2, "total": 2, "failed": 0}
    assert report.variables[0].variable == "U"
    assert report.variables[0].status == "ok"
    assert report.variables[0].threshold == 0.05
    assert report.variables[0].valid_count == 3
    assert report.variables[0].total_count == 3
    assert report.variables[1].threshold is None
    assert report.diagnostics["tc"]["status"] == "pending"


def test_compare_files_includes_tc_diagnostics_when_requested(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    _write_dataset_with_tc_fields(reference)
    _write_dataset_with_tc_fields(candidate)

    report = compare_files(reference, candidate, variables=("U",), include_tc_diagnostics=True)

    assert report.status == "ok"
    tc = report.diagnostics["tc"]
    assert tc.status == "ok"
    assert tc.reference.center.j == 1
    assert tc.reference.center.i == 2
    assert tc.reference.minimum_slp_status == "ok"
    assert tc.reference.minimum_slp_hpa == 950.0
    assert tc.reference.mslp_proxy_hpa == 950.0
    assert "PSFC-min proxy" in tc.reference.mslp_proxy_label
    assert tc.candidate.vmax_10m_ms == 50.0
    assert tc.center_error_km == 0.0
    assert tc.mslp_proxy_abs_error_hpa == 0.0


def test_compare_files_fails_status_when_tc_diagnostics_fail(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    _write_dataset_with_tc_fields(reference)
    _write_dataset_with_tc_fields(candidate)

    with netCDF4.Dataset(candidate, "a") as dataset:
        dataset.variables["SLP"][0, :, :] = 1000.0
        dataset.variables["SLP"][0, 2, 3] = 960.0

    report = compare_files(reference, candidate, variables=("U",), include_tc_diagnostics=True)

    assert report.status == "failed"
    assert report.summary["failed"] == 0
    assert report.summary["diagnostics_failed"] == 1
    assert report.diagnostics["tc"].status == "threshold_exceeded"


def test_compare_files_reports_threshold_exceeded(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    _write_dataset(reference, {"U": np.array([1.0, 2.0, 3.0])})
    _write_dataset(candidate, {"U": np.array([2.0, 3.0, 4.0])})

    report = compare_files(reference, candidate, variables=("U",), thresholds={"U": 0.01})

    assert report.status == "failed"
    assert report.summary["threshold_exceeded"] == 1
    assert report.summary["failed"] == 1
    assert report.variables[0].status == "threshold_exceeded"
    assert "exceeds threshold" in report.variables[0].message


def test_compare_files_reports_missing_and_shape_mismatch_statuses(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    _write_dataset(
        reference,
        {
            "U": np.array([1.0, 2.0, 3.0]),
            "V": np.array([4.0, 5.0]),
        },
    )
    _write_dataset(candidate, {"U": np.array([1.0, 2.0])})

    report = compare_files(reference, candidate, variables=("U", "V"))

    assert report.status == "failed"
    assert [result.status for result in report.variables] == ["shape_mismatch", "missing_candidate"]
    assert report.variables[0].reference_shape == (3,)
    assert report.variables[0].candidate_shape == (2,)
    assert report.variables[1].reference_shape == (2,)


def test_compare_files_reports_missing_strict_fields_explicitly(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    candidate_fields = _strict_core_fields()
    del candidate_fields["PH"]
    del candidate_fields["QVAPOR"]

    _write_dataset(reference, _strict_core_fields())
    _write_dataset(candidate, candidate_fields)

    report = compare_files(reference, candidate, variables=STRICT_CORE_VARIABLES)
    results = {result.variable: result for result in report.variables}
    strict = report.diagnostics["strict_fields"]

    assert report.status == "failed"
    assert results["PH"].status == "missing_candidate"
    assert results["PH"].message == "required strict field missing from candidate"
    assert results["QVAPOR"].status == "missing_candidate"
    assert strict["status"] == "failed"
    assert strict["missing_candidate"] == ["PH", "QVAPOR"]
    assert strict["numeric_failures"] == []
    assert strict["first_failure"]["variable"] == "PH"
    assert strict["first_failure"]["failure_kind"] == "missing_field"


def test_compare_files_reports_numeric_first_failure_when_strict_fields_present(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    candidate_fields = _strict_core_fields()
    candidate_fields["U"] = np.array([2.0, 3.0, 4.0])

    _write_dataset(reference, _strict_core_fields())
    _write_dataset(candidate, candidate_fields)

    report = compare_files(reference, candidate, variables=STRICT_CORE_VARIABLES)
    strict = report.diagnostics["strict_fields"]

    assert report.status == "failed"
    assert strict["missing_candidate"] == []
    assert strict["numeric_failures"] == ["U"]
    assert strict["first_failure"]["variable"] == "U"
    assert strict["first_failure"]["status"] == "threshold_exceeded"
    assert strict["first_failure"]["failure_kind"] == "numeric_rmse"
    assert strict["first_failure"]["rmse"] is not None
    assert strict["first_failure"]["max_abs_error"] == 1.0


def test_compare_files_rejects_marked_diagnostic_oracle_candidate(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    fields = _strict_core_fields()
    _write_dataset(reference, fields)
    _write_dataset(
        candidate,
        fields,
        attrs={
            "TYWRF_DIAGNOSTIC_ONLY": "true",
            "TYWRF_CANDIDATE_KIND": "reference_delta_oracle",
        },
    )

    report = compare_files(reference, candidate, variables=STRICT_CORE_VARIABLES)
    metadata = report.diagnostics["candidate_metadata"]

    assert report.status == "failed"
    assert report.summary["failed"] == 0
    assert report.summary["candidate_metadata_failed"] == 1
    assert report.diagnostics["strict_fields"]["status"] == "ok"
    assert metadata["status"] == "failed"
    assert "TYWRF_DIAGNOSTIC_ONLY=true" in metadata["disqualifiers"]
    assert "TYWRF_CANDIDATE_KIND=reference_delta_oracle" in metadata["disqualifiers"]


def test_compare_files_ignores_nonfinite_pairs_and_writes_strict_json(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    _write_dataset(reference, {"U": np.array([1.0, np.nan, 3.0])})
    _write_dataset(candidate, {"U": np.array([2.0, 2.0, 3.0])})

    report = compare_files(reference, candidate, variables=("U",), thresholds=None)
    payload = json.loads(report_to_json(report))

    assert report.status == "ok"
    assert report.variables[0].valid_count == 2
    assert report.variables[0].total_count == 3
    assert payload["variables"][0]["valid_count"] == 2


def test_main_writes_json_report(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "compare.json"
    _write_dataset(reference, {"U": np.array([1.0, 2.0, 3.0])})
    _write_dataset(candidate, {"U": np.array([1.0, 2.0, 3.0])})

    exit_code = compare_main(
        [
            str(reference),
            str(candidate),
            "--variables",
            "U",
            "--output",
            str(output),
            "--pretty",
        ]
    )

    assert exit_code == 0
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["status"] == "ok"
    assert payload["variables"][0]["variable"] == "U"


def test_main_writes_json_report_with_tc_diagnostics_flag(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "compare-tc.json"
    _write_dataset_with_tc_fields(reference)
    _write_dataset_with_tc_fields(candidate)

    exit_code = compare_main(
        [
            str(reference),
            str(candidate),
            "--variables",
            "U",
            "--tc-diagnostics",
            "--output",
            str(output),
            "--pretty",
        ]
    )

    assert exit_code == 0
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["diagnostics"]["tc"]["status"] == "ok"
    assert payload["diagnostics"]["tc"]["reference"]["center"]["latitude"] == 11.0
    assert payload["diagnostics"]["tc"]["reference"]["mslp_proxy_hpa"] == 950.0
    assert payload["diagnostics"]["tc"]["reference"]["rainfall"]["maximum_mm"] == 11.5
