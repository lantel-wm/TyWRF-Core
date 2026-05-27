import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.analyze_cycle_delta import (
    analyze_cycle_delta,
    analyze_reference_cycle,
    main as delta_main,
    report_to_json,
    resolve_cycle_files,
)


START = "2025-07-26_00:00:00"
END = "2025-07-26_00:10:00"


def _write_dataset(path: Path, variables: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def test_analyze_cycle_delta_reports_json_for_explicit_files(tmp_path: Path) -> None:
    start = tmp_path / "wrfout_d02_start"
    end = tmp_path / "wrfout_d02_end"
    _write_dataset(start, {"U": np.full((2, 2), 9.0)})
    _write_dataset(end, {"U": np.full((2, 2), 10.0)})

    report = analyze_cycle_delta(start, end, variables=("U",), thresholds={"U": 0.05})
    payload = json.loads(report_to_json(report, pretty=True))

    assert report.status == "computed"
    assert payload["start_file"] == str(start)
    assert payload["end_file"] == str(end)
    assert payload["variables"][0]["variable"] == "U"
    assert payload["variables"][0]["shape"] == [2, 2]
    assert payload["variables"][0]["finite_pair_count"] == 4
    assert payload["variables"][0]["total_pair_count"] == 4
    assert math.isclose(payload["variables"][0]["rmse_persistence"], 1.0)
    assert math.isclose(payload["variables"][0]["normalized_rmse"], 0.1)
    assert math.isclose(payload["variables"][0]["mean_abs_delta"], 1.0)
    assert math.isclose(payload["variables"][0]["max_abs_delta"], 1.0)
    assert payload["variables"][0]["exceeds_threshold"] is True
    assert payload["summary"]["first_failing_variable"] == "U"


def test_threshold_failure_order_uses_variable_order_and_largest_delta_sort(
    tmp_path: Path,
) -> None:
    start = tmp_path / "start.nc"
    end = tmp_path / "end.nc"
    _write_dataset(
        start,
        {
            "U": np.full((2,), 9.0),
            "V": np.full((2,), 2.0),
            "T": np.full((2,), 100.0),
        },
    )
    _write_dataset(
        end,
        {
            "U": np.full((2,), 10.0),
            "V": np.full((2,), 1.0),
            "T": np.full((2,), 99.0),
        },
    )

    report = analyze_cycle_delta(
        start,
        end,
        variables=("U", "V", "T"),
        thresholds={"U": 0.05, "V": 0.05, "T": 0.05},
    )
    payload = json.loads(report_to_json(report))

    assert payload["summary"]["strict_threshold_exceeded"] == 2
    assert payload["summary"]["first_failing_variable"] == "U"
    largest = payload["summary"]["largest_normalized_deltas"]
    assert [item["variable"] for item in largest] == ["V", "U", "T"]
    assert largest[0]["normalized_rmse"] > largest[1]["normalized_rmse"]


def test_nonfinite_pairs_are_ignored_and_json_is_strict(tmp_path: Path) -> None:
    start = tmp_path / "start.nc"
    end = tmp_path / "end.nc"
    _write_dataset(start, {"U": np.array([1.0, np.nan, 5.0, np.inf])})
    _write_dataset(end, {"U": np.array([2.0, 3.0, np.nan, 4.0])})

    report = analyze_cycle_delta(start, end, variables=("U",), thresholds={"U": 0.05})
    payload = json.loads(report_to_json(report))

    assert payload["variables"][0]["finite_pair_count"] == 1
    assert payload["variables"][0]["total_pair_count"] == 4
    assert math.isclose(payload["variables"][0]["rmse_persistence"], 1.0)
    assert math.isclose(payload["variables"][0]["normalized_rmse"], 0.5)


def test_reference_dir_path_resolution_and_report_metadata(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    start = reference_dir / f"wrfout_d02_{START}"
    end = reference_dir / f"wrfout_d02_{END}"
    _write_dataset(start, {"U": np.ones((2, 2))})
    _write_dataset(end, {"U": np.ones((2, 2))})

    resolved_start, resolved_end, domain, start_time, end_time = resolve_cycle_files(
        reference_dir=reference_dir,
        domain="d02",
        start=START,
        end=END,
    )
    report = analyze_reference_cycle(
        reference_dir,
        domain="d02",
        start=START,
        end=END,
        variables=("U",),
    )

    assert resolved_start == start
    assert resolved_end == end
    assert domain == "d02"
    assert start_time == START
    assert end_time == END
    assert report.domain == "d02"
    assert report.start_time == START
    assert report.end_time == END
    assert report.summary["first_failing_variable"] is None


def test_main_writes_reference_dir_json_without_gate_exit_failure(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    output = tmp_path / "delta.json"
    _write_dataset(reference_dir / f"wrfout_d02_{START}", {"U": np.full((2,), 9.0)})
    _write_dataset(reference_dir / f"wrfout_d02_{END}", {"U": np.full((2,), 10.0)})

    exit_code = delta_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--domain",
            "d02",
            "--start",
            START,
            "--end",
            END,
            "--variables",
            "U",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["status"] == "computed"
    assert payload["summary"]["first_failing_variable"] == "U"
    assert payload["variables"][0]["exceeds_threshold"] is True
