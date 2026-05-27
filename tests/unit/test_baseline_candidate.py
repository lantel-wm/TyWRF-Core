from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.baseline_candidate import build_baseline_candidate, normalize_candidate_mode


def _make_wrfout(path: Path, valid_time: str, value: float, *, dx: float = 2000.0) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.setncattr("DX", dx)
        dataset.setncattr("DY", dx)
        dataset.createDimension("Time", None)
        dataset.createDimension("DateStrLen", 19)
        dataset.createDimension("south_north", 2)
        dataset.createDimension("west_east", 3)

        times = dataset.createVariable("Times", "S1", ("Time", "DateStrLen"))
        times[0, :] = np.array(list(valid_time), dtype="S1")

        psfc = dataset.createVariable("PSFC", "f4", ("Time", "south_north", "west_east"))
        psfc[0, :, :] = np.full((2, 3), value, dtype=np.float32)

        slp = dataset.createVariable("SLP", "f4", ("Time", "south_north", "west_east"))
        slp.units = "hPa"
        slp[0, :, :] = np.full((2, 3), value / 100.0, dtype=np.float32)


def _read_time(dataset: netCDF4.Dataset) -> str:
    return b"".join(dataset.variables["Times"][0, :]).decode("ascii").strip()


def test_persistence_candidate_uses_start_file_and_rewrites_valid_time(tmp_path: Path) -> None:
    start = tmp_path / "wrfout_d02_2025-07-26_00:00:00"
    end = tmp_path / "wrfout_d02_2025-07-26_06:00:00"
    candidate = tmp_path / "candidate" / end.name
    _make_wrfout(start, "2025-07-26_00:00:00", 100000.0)
    _make_wrfout(end, "2025-07-26_06:00:00", 99000.0)

    metadata = build_baseline_candidate(
        start,
        end,
        candidate,
        domain="d02",
        start_time="2025-07-26_00:00:00",
        end_time="2025-07-26_06:00:00",
        mode="persistence",
        variables=("Times", "PSFC"),
    )

    assert metadata.mode == "persistence"
    assert metadata.reference_copy is False
    assert metadata.integrator_output is False
    assert metadata.minimum_slp_gate_ready is False
    assert metadata.expected_to_meet_thresholds is False
    assert metadata.copied_slp_candidates == []
    assert metadata.source == str(start)
    assert metadata.d02_resolution_check == "d02_2km"

    with netCDF4.Dataset(candidate) as dataset:
        assert _read_time(dataset) == "2025-07-26_06:00:00"
        np.testing.assert_array_equal(dataset.variables["PSFC"][:], np.full((1, 2, 3), 100000.0))
        assert dataset.getncattr("TYWRF_BASELINE_MODE") == "persistence"
        assert dataset.getncattr("TYWRF_REFERENCE_COPY") == "false"
        assert dataset.getncattr("TYWRF_INTEGRATOR_OUTPUT") == "false"


def test_reference_copy_candidate_is_explicitly_marked_as_gate_only(tmp_path: Path) -> None:
    start = tmp_path / "wrfout_d02_2025-07-26_00:00:00"
    end = tmp_path / "wrfout_d02_2025-07-26_06:00:00"
    candidate = tmp_path / "candidate" / end.name
    _make_wrfout(start, "2025-07-26_00:00:00", 100000.0)
    _make_wrfout(end, "2025-07-26_06:00:00", 99000.0)

    metadata = build_baseline_candidate(
        start,
        end,
        candidate,
        domain="d02",
        start_time="2025-07-26_00:00:00",
        end_time="2025-07-26_06:00:00",
        mode="reference-copy",
        variables=("Times", "PSFC"),
    )

    assert metadata.mode == "reference_copy"
    assert metadata.reference_copy is True
    assert metadata.validation_gate_only is True
    assert metadata.integrator_output is False
    assert metadata.minimum_slp_gate_ready is False
    assert metadata.expected_to_meet_thresholds is False
    assert metadata.copied_slp_candidates == []
    assert metadata.source == str(end)

    with netCDF4.Dataset(candidate) as dataset:
        assert _read_time(dataset) == "2025-07-26_06:00:00"
        np.testing.assert_array_equal(dataset.variables["PSFC"][:], np.full((1, 2, 3), 99000.0))
        assert dataset.getncattr("TYWRF_BASELINE_MODE") == "reference_copy"
        assert dataset.getncattr("TYWRF_REFERENCE_COPY") == "true"
        assert dataset.getncattr("TYWRF_VALIDATION_GATE_ONLY") == "true"
        assert dataset.getncattr("TYWRF_INTEGRATOR_OUTPUT") == "false"


def test_reference_copy_candidate_marks_gate_ready_when_real_slp_is_copied(
    tmp_path: Path,
) -> None:
    start = tmp_path / "wrfout_d02_2025-07-26_00:00:00"
    end = tmp_path / "wrfout_d02_2025-07-26_06:00:00"
    candidate = tmp_path / "candidate" / end.name
    _make_wrfout(start, "2025-07-26_00:00:00", 100000.0)
    _make_wrfout(end, "2025-07-26_06:00:00", 99000.0)

    metadata = build_baseline_candidate(
        start,
        end,
        candidate,
        domain="d02",
        start_time="2025-07-26_00:00:00",
        end_time="2025-07-26_06:00:00",
        mode="reference-copy",
        variables=("Times", "PSFC", "SLP"),
    )

    assert metadata.minimum_slp_gate_ready is True
    assert metadata.expected_to_meet_thresholds is True
    assert metadata.copied_slp_candidates == ["SLP"]


def test_d02_candidate_rejects_non_2km_source(tmp_path: Path) -> None:
    start = tmp_path / "wrfout_d02_2025-07-26_00:00:00"
    end = tmp_path / "wrfout_d02_2025-07-26_06:00:00"
    _make_wrfout(start, "2025-07-26_00:00:00", 100000.0, dx=3000.0)
    _make_wrfout(end, "2025-07-26_06:00:00", 99000.0)

    with pytest.raises(ValueError, match="d02 source resolution must remain 2 km"):
        build_baseline_candidate(
            start,
            end,
            tmp_path / "candidate" / end.name,
            domain="d02",
            start_time="2025-07-26_00:00:00",
            end_time="2025-07-26_06:00:00",
            mode="persistence",
            variables=("Times", "PSFC"),
        )


def test_mode_normalization() -> None:
    assert normalize_candidate_mode("reference-copy") == "reference_copy"
    assert normalize_candidate_mode("persistence") == "persistence"
    with pytest.raises(ValueError, match="baseline candidate mode"):
        normalize_candidate_mode("integrator")
