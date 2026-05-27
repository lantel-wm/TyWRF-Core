from pathlib import Path

from tools.run_6h_cycle_test import (
    build_baseline_candidate_report,
    build_cycle_plan,
    format_wrf_time,
    main as cycle_main,
    parse_wrf_time,
    wrfout_filename,
)

import netCDF4
import numpy as np


def test_wrf_time_parsing_and_filename_formatting() -> None:
    value = parse_wrf_time("2025-07-26_00:00:00")

    assert format_wrf_time(value) == "2025-07-26_00:00:00"
    assert wrfout_filename("d02", value) == "wrfout_d02_2025-07-26_00:00:00"
    assert parse_wrf_time("2025-07-26T00:00:00") == value


def test_build_cycle_plan_resolves_d01_d02_start_and_end_files(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    reference_dir.mkdir()
    candidate_dir.mkdir()

    (reference_dir / "wrfout_d01_2025-07-26_00:00:00").touch()
    (reference_dir / "wrfout_d01_2025-07-26_06:00:00").touch()
    (candidate_dir / "wrfout_d01_2025-07-26_00:00:00").touch()
    (candidate_dir / "wrfout_d02_2025-07-26_06:00:00").touch()

    plan = build_cycle_plan(
        reference_dir,
        candidate_dir,
        "2025-07-26_00:00:00",
        domains=("d01", "d02"),
    )

    assert plan.status == "dry_run"
    assert plan.start_time == "2025-07-26_00:00:00"
    assert plan.end_time == "2025-07-26_06:00:00"
    assert plan.hours == 6
    assert plan.minutes == 360
    assert [domain.domain for domain in plan.domains] == ["d01", "d02"]
    assert plan.domains[0].reference_start.endswith("wrfout_d01_2025-07-26_00:00:00")
    assert plan.domains[0].candidate_start.endswith("wrfout_d01_2025-07-26_00:00:00")
    assert plan.domains[0].reference_start_exists is True
    assert plan.domains[0].reference_end_exists is True
    assert plan.domains[0].candidate_start_exists is True
    assert plan.domains[0].candidate_end_exists is False
    assert plan.domains[1].reference_end.endswith("wrfout_d02_2025-07-26_06:00:00")
    assert plan.domains[1].candidate_end_exists is True
    assert plan.diagnostics["tc"]["status"] == "pending"


def test_cycle_main_dry_run_writes_json_plan(tmp_path: Path) -> None:
    output = tmp_path / "cycle-plan.json"

    exit_code = cycle_main(
        [
            "--reference-dir",
            str(tmp_path / "reference"),
            "--candidate-dir",
            str(tmp_path / "candidate"),
            "--start",
            "2025-07-26_00:00:00",
            "--dry-run",
            "--output",
            str(output),
            "--pretty",
        ]
    )

    assert exit_code == 0
    text = output.read_text(encoding="utf-8")
    assert '"status": "dry_run"' in text
    assert '"domain": "d01"' in text
    assert '"domain": "d02"' in text


def _make_d02_wrfout(path: Path, valid_time: str, value: float) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.setncattr("DX", 2000.0)
        dataset.setncattr("DY", 2000.0)
        dataset.createDimension("Time", None)
        dataset.createDimension("DateStrLen", 19)
        dataset.createDimension("south_north", 2)
        dataset.createDimension("west_east", 2)

        times = dataset.createVariable("Times", "S1", ("Time", "DateStrLen"))
        times[0, :] = np.array(list(valid_time), dtype="S1")

        psfc = dataset.createVariable("PSFC", "f4", ("Time", "south_north", "west_east"))
        psfc[0, :, :] = np.full((2, 2), value, dtype=np.float32)


def test_build_baseline_candidate_report_generates_d02_multi_cycle_files(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    reference_dir.mkdir()
    for hour, value in ((0, 100000.0), (6, 99000.0), (12, 98000.0)):
        valid_time = f"2025-07-26_{hour:02d}:00:00"
        _make_d02_wrfout(reference_dir / f"wrfout_d02_{valid_time}", valid_time, value)

    report = build_baseline_candidate_report(
        reference_dir,
        candidate_dir,
        start="2025-07-26_00:00:00",
        end="2025-07-26_12:00:00",
        hours=6,
        interval_hours=6,
        domains=("d02",),
        mode="persistence",
        variables=("Times", "PSFC"),
    )

    assert report["status"] == "baseline_candidate_generated"
    assert report["mode"] == "persistence"
    assert report["cycle_count"] == 2
    assert [Path(item["candidate"]).name for item in report["candidates"]] == [
        "wrfout_d02_2025-07-26_06:00:00",
        "wrfout_d02_2025-07-26_12:00:00",
    ]

    with netCDF4.Dataset(candidate_dir / "wrfout_d02_2025-07-26_12:00:00") as dataset:
        assert dataset.getncattr("TYWRF_REFERENCE_COPY") == "false"
        np.testing.assert_array_equal(dataset.variables["PSFC"][:], np.full((1, 2, 2), 99000.0))


def test_build_baseline_candidate_report_supports_10_minute_d02_window(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    reference_dir.mkdir()
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
    )
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:10:00",
        "2025-07-26_00:10:00",
        99900.0,
    )

    report = build_baseline_candidate_report(
        reference_dir,
        candidate_dir,
        start="2025-07-26_00:00:00",
        end=None,
        hours=6,
        minutes=10,
        interval_hours=6,
        domains=("d02",),
        mode="persistence",
        variables=("Times", "PSFC"),
    )

    assert report["cycle_count"] == 1
    assert report["hours"] == 0
    assert report["minutes"] == 10
    assert report["interval_hours"] == 0
    assert report["interval_minutes"] == 10
    candidate = report["candidates"][0]
    assert candidate["candidate"].endswith("wrfout_d02_2025-07-26_00:10:00")
    assert candidate["hours"] == 0
    assert candidate["minutes"] == 10
    assert candidate["integrator_output"] is False
    assert candidate["reference_copy"] is False

    with netCDF4.Dataset(candidate_dir / "wrfout_d02_2025-07-26_00:10:00") as dataset:
        assert dataset.getncattr("TYWRF_CYCLE_HOURS") == 0
        assert dataset.getncattr("TYWRF_CYCLE_MINUTES") == 10
        assert dataset.getncattr("TYWRF_INTEGRATOR_OUTPUT") == "false"


def test_build_baseline_candidate_report_accepts_minute_interval(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    reference_dir.mkdir()
    for minute, value in ((0, 100000.0), (10, 99000.0), (20, 98000.0)):
        valid_time = f"2025-07-26_00:{minute:02d}:00"
        _make_d02_wrfout(reference_dir / f"wrfout_d02_{valid_time}", valid_time, value)

    report = build_baseline_candidate_report(
        reference_dir,
        candidate_dir,
        start="2025-07-26_00:00:00",
        end="2025-07-26_00:20:00",
        hours=6,
        interval_hours=6,
        interval_minutes=10,
        domains=("d02",),
        mode="persistence",
        variables=("Times", "PSFC"),
    )

    assert report["cycle_count"] == 2
    assert [Path(item["candidate"]).name for item in report["candidates"]] == [
        "wrfout_d02_2025-07-26_00:10:00",
        "wrfout_d02_2025-07-26_00:20:00",
    ]


def test_cycle_main_reference_copy_writes_metadata_json(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    output = tmp_path / "report.json"
    reference_dir.mkdir()
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
    )
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_06:00:00",
        "2025-07-26_06:00:00",
        99000.0,
    )

    exit_code = cycle_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--candidate-dir",
            str(candidate_dir),
            "--start",
            "2025-07-26_00:00:00",
            "--domain",
            "d02",
            "--mode",
            "reference-copy",
            "--variables",
            "Times",
            "PSFC",
            "--output",
            str(output),
            "--pretty",
        ]
    )

    assert exit_code == 0
    text = output.read_text(encoding="utf-8")
    assert '"mode": "reference_copy"' in text
    assert '"reference_copy": true' in text
    with netCDF4.Dataset(candidate_dir / "wrfout_d02_2025-07-26_06:00:00") as dataset:
        assert dataset.getncattr("TYWRF_REFERENCE_COPY") == "true"
