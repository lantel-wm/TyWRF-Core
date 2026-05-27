from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.skeleton_cycle_driver import (
    build_skeleton_cycle_candidate,
    main as skeleton_main,
    normalize_skeleton_mode,
)


def _make_d02_wrfout(path: Path, valid_time: str, value: float, *, dx: float = 2000.0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.setncattr("DX", dx)
        dataset.setncattr("DY", dx)
        dataset.createDimension("Time", None)
        dataset.createDimension("DateStrLen", 19)
        dataset.createDimension("south_north", 2)
        dataset.createDimension("west_east", 2)

        times = dataset.createVariable("Times", "S1", ("Time", "DateStrLen"))
        times[0, :] = np.array(list(valid_time), dtype="S1")

        psfc = dataset.createVariable("PSFC", "f4", ("Time", "south_north", "west_east"))
        psfc[0, :, :] = np.full((2, 2), value, dtype=np.float32)


def _read_time(dataset: netCDF4.Dataset) -> str:
    return b"".join(dataset.variables["Times"][0, :]).decode("ascii").strip()


def test_skeleton_driver_generates_d02_cycle_end_candidate(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
    )

    report = build_skeleton_cycle_candidate(
        reference_dir,
        candidate_dir,
        start="2025-07-26_00:00:00",
        mode="persistence",
        variables=("Times", "PSFC"),
    )

    assert report.status == "skeleton_candidate_generated"
    assert report.mode == "persistence"
    assert report.skeleton is True
    assert report.not_physical is True
    assert report.integrator_output is False
    assert report.validation_gate_only is True
    assert report.domain == "d02"
    assert report.start == "2025-07-26_00:00:00"
    assert report.end == "2025-07-26_06:00:00"
    assert report.hours == 6
    assert report.minutes == 360
    assert report.source.endswith("wrfout_d02_2025-07-26_00:00:00")
    assert report.candidate.endswith("wrfout_d02_2025-07-26_06:00:00")
    assert "tools/cycle_gate.py" in report.suggested_next_step["command"]
    assert "--domain d02" in report.suggested_next_step["command"]

    with netCDF4.Dataset(candidate_dir / "wrfout_d02_2025-07-26_06:00:00") as dataset:
        assert _read_time(dataset) == "2025-07-26_06:00:00"
        np.testing.assert_array_equal(dataset.variables["PSFC"][:], np.full((1, 2, 2), 100000.0))
        assert dataset.getncattr("TYWRF_CANDIDATE_KIND") == "skeleton_candidate"
        assert dataset.getncattr("TYWRF_SKELETON") == "true"
        assert dataset.getncattr("TYWRF_SKELETON_MODE") == "persistence"
        assert dataset.getncattr("TYWRF_SKELETON_NOT_PHYSICAL") == "true"
        assert dataset.getncattr("TYWRF_NOT_PHYSICAL") == "true"
        assert dataset.getncattr("TYWRF_INTEGRATOR_OUTPUT") == "false"
        assert dataset.getncattr("TYWRF_VALIDATION_GATE_ONLY") == "true"


def test_identity_mode_is_explicitly_reported_as_skeleton_identity(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
    )

    report = build_skeleton_cycle_candidate(
        reference_dir,
        candidate_dir,
        start="2025-07-26_00:00:00",
        end="2025-07-26_06:00:00",
        mode="identity",
        variables=("Times", "PSFC"),
    )

    assert report.mode == "identity"
    assert report.integrator_output is False
    with netCDF4.Dataset(candidate_dir / "wrfout_d02_2025-07-26_06:00:00") as dataset:
        assert dataset.getncattr("TYWRF_SKELETON_MODE") == "identity"
        assert dataset.getncattr("TYWRF_BASELINE_MODE") == "persistence"


def test_skeleton_driver_generates_10min_cycle_end_candidate(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
    )

    report = build_skeleton_cycle_candidate(
        reference_dir,
        candidate_dir,
        start="2025-07-26_00:00:00",
        end="2025-07-26_00:10:00",
        variables=("Times", "PSFC"),
    )

    assert report.end == "2025-07-26_00:10:00"
    assert report.hours == 0
    assert report.minutes == 10
    assert report.candidate.endswith("wrfout_d02_2025-07-26_00:10:00")
    assert "--end 2025-07-26_00:10:00" in report.suggested_next_step["command"]
    assert "--interval-minutes 10" in report.suggested_next_step["command"]

    with netCDF4.Dataset(candidate_dir / "wrfout_d02_2025-07-26_00:10:00") as dataset:
        assert _read_time(dataset) == "2025-07-26_00:10:00"
        assert dataset.getncattr("TYWRF_CYCLE_HOURS") == 0
        assert dataset.getncattr("TYWRF_CYCLE_MINUTES") == 10
        assert dataset.getncattr("TYWRF_NOT_PHYSICAL") == "true"
        assert dataset.getncattr("TYWRF_INTEGRATOR_OUTPUT") == "false"


def test_skeleton_driver_minutes_option_resolves_cycle_end(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
    )

    report = build_skeleton_cycle_candidate(
        reference_dir,
        candidate_dir,
        start="2025-07-26_00:00:00",
        hours=None,
        minutes=10,
        variables=("Times", "PSFC"),
    )

    assert report.end == "2025-07-26_00:10:00"
    assert report.minutes == 10


def test_skeleton_driver_rejects_non_positive_cycle_length(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="positive"):
        build_skeleton_cycle_candidate(
            tmp_path / "reference",
            tmp_path / "candidate",
            start="2025-07-26_00:00:00",
            end="2025-07-26_00:00:00",
            variables=("Times", "PSFC"),
        )


def test_skeleton_driver_reuses_d02_resolution_check(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
        dx=3000.0,
    )

    with pytest.raises(ValueError, match="d02 source resolution must remain 2 km"):
        build_skeleton_cycle_candidate(
            reference_dir,
            tmp_path / "candidate",
            start="2025-07-26_00:00:00",
            variables=("Times", "PSFC"),
        )


def test_mode_normalization() -> None:
    assert normalize_skeleton_mode("persistence") == "persistence"
    assert normalize_skeleton_mode("identity") == "identity"
    with pytest.raises(ValueError, match="skeleton mode"):
        normalize_skeleton_mode("reference-copy")


def test_skeleton_driver_cli_writes_report_json(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    output = tmp_path / "report.json"
    _make_d02_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
    )

    exit_code = skeleton_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--candidate-dir",
            str(candidate_dir),
            "--start",
            "2025-07-26_00:00:00",
            "--mode",
            "identity",
            "--minutes",
            "10",
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
    assert '"status": "skeleton_candidate_generated"' in text
    assert '"skeleton": true' in text
    assert '"integrator_output": false' in text
    assert '"mode": "identity"' in text
    assert '"end": "2025-07-26_00:10:00"' in text
    assert '"minutes": 10' in text
    assert "tools/cycle_gate.py" in text
