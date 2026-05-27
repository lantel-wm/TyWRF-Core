import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.compare_wrfout import compare_files, main as compare_main
from tools.compare_wrfout import normalized_rmse, report_to_json


def _write_dataset(path: Path, variables: dict[str, np.ndarray]) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        for name, values in variables.items():
            data = np.asarray(values, dtype=np.float64)
            dimensions = []
            for axis, size in enumerate(data.shape):
                dimension = f"{name}_dim_{axis}"
                dataset.createDimension(dimension, size)
                dimensions.append(dimension)
            variable = dataset.createVariable(name, "f8", tuple(dimensions))
            variable[:] = data


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
