import json
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.run_skeleton_cycle import build_skeleton_cycle_run, main as run_skeleton_main


def _make_wrfout(path: Path, valid_time: str, value: float, *, dx: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def _read_time(dataset: netCDF4.Dataset) -> str:
    return b"".join(dataset.variables["Times"][0, :]).decode("ascii").strip()


def _make_dual_domain_start(reference_dir: Path) -> None:
    _make_wrfout(
        reference_dir / "wrfout_d01_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
        dx=10000.0,
    )
    _make_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        99900.0,
        dx=2000.0,
    )


def test_run_skeleton_cycle_generates_d01_and_d02_candidates(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _make_dual_domain_start(reference_dir)

    report = build_skeleton_cycle_run(
        reference_dir,
        candidate_dir,
        start="2025-07-26_00:00:00",
        variables=("Times", "PSFC"),
    )

    assert report.status == "skeleton_candidates_generated"
    assert report.skeleton is True
    assert report.not_physical is True
    assert report.integrator_output is False
    assert report.validation_gate_only is True
    assert report.integrator["status"] == "not_run"
    assert [domain.domain for domain in report.domains] == ["d01", "d02"]
    assert set(report.suggested_gate_commands) == {"d01", "d02"}

    reports_by_domain = {domain.domain: domain for domain in report.domains}
    assert reports_by_domain["d01"].d02_resolution_check == "not_applicable"
    assert reports_by_domain["d02"].d02_resolution_check == "d02_2km"

    for domain, expected_value in (("d01", 100000.0), ("d02", 99900.0)):
        domain_report = reports_by_domain[domain]
        assert domain_report.source.endswith(f"wrfout_{domain}_2025-07-26_00:00:00")
        assert domain_report.candidate.endswith(f"wrfout_{domain}_2025-07-26_06:00:00")
        assert domain_report.reference_end.endswith(f"wrfout_{domain}_2025-07-26_06:00:00")
        assert domain_report.copied_variables == ["Times", "PSFC"]
        assert domain_report.missing_variables == []
        assert f"--domain {domain}" in domain_report.suggested_gate_command
        assert "tools/cycle_gate.py" in domain_report.suggested_gate_command

        with netCDF4.Dataset(candidate_dir / f"wrfout_{domain}_2025-07-26_06:00:00") as dataset:
            assert _read_time(dataset) == "2025-07-26_06:00:00"
            np.testing.assert_array_equal(dataset.variables["PSFC"][:], np.full((1, 2, 3), expected_value))
            assert dataset.getncattr("TYWRF_CANDIDATE_KIND") == "skeleton_candidate"
            assert dataset.getncattr("TYWRF_SKELETON") == "true"
            assert dataset.getncattr("TYWRF_NOT_PHYSICAL") == "true"
            assert dataset.getncattr("TYWRF_INTEGRATOR_OUTPUT") == "false"
            assert dataset.getncattr("TYWRF_VALIDATION_GATE_ONLY") == "true"
            assert dataset.getncattr("TYWRF_CANDIDATE_DOMAIN") == domain
            assert dataset.getncattr("TYWRF_SKELETON_ORCHESTRATOR") == "tools/run_skeleton_cycle.py"


def test_run_skeleton_cycle_d01_only_does_not_require_2km(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    _make_wrfout(
        reference_dir / "wrfout_d01_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
        dx=12000.0,
    )

    report = build_skeleton_cycle_run(
        reference_dir,
        tmp_path / "candidate",
        start="2025-07-26_00:00:00",
        domains=("d01",),
        variables=("Times", "PSFC"),
    )

    assert [domain.domain for domain in report.domains] == ["d01"]
    assert report.domains[0].dx_m == 12000.0
    assert report.domains[0].d02_resolution_check == "not_applicable"


def test_run_skeleton_cycle_reuses_d02_resolution_check(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    _make_wrfout(
        reference_dir / "wrfout_d01_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        100000.0,
        dx=10000.0,
    )
    _make_wrfout(
        reference_dir / "wrfout_d02_2025-07-26_00:00:00",
        "2025-07-26_00:00:00",
        99900.0,
        dx=3000.0,
    )

    with pytest.raises(ValueError, match="d02 source resolution must remain 2 km"):
        build_skeleton_cycle_run(
            reference_dir,
            tmp_path / "candidate",
            start="2025-07-26_00:00:00",
            variables=("Times", "PSFC"),
        )


def test_run_skeleton_cycle_reports_missing_variables_when_allowed(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    _make_dual_domain_start(reference_dir)

    report = build_skeleton_cycle_run(
        reference_dir,
        tmp_path / "candidate",
        start="2025-07-26_00:00:00",
        variables=("Times", "PSFC", "QRAIN"),
        allow_missing=True,
    )

    for domain_report in report.domains:
        assert domain_report.copied_variables == ["Times", "PSFC"]
        assert domain_report.missing_variables == ["QRAIN"]


def test_run_skeleton_cycle_cli_writes_machine_readable_report(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    output = tmp_path / "report.json"
    _make_dual_domain_start(reference_dir)

    exit_code = run_skeleton_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--candidate-dir",
            str(candidate_dir),
            "--start",
            "2025-07-26_00:00:00",
            "--mode",
            "identity",
            "--variables",
            "Times",
            "PSFC",
            "--output",
            str(output),
            "--pretty",
        ]
    )

    assert exit_code == 0
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["status"] == "skeleton_candidates_generated"
    assert report["mode"] == "identity"
    assert report["skeleton"] is True
    assert report["not_physical"] is True
    assert report["integrator_output"] is False
    assert report["integrator"]["integrator_output"] is False
    assert {domain["domain"] for domain in report["domains"]} == {"d01", "d02"}
    assert "d01" in report["suggested_gate_commands"]
    assert "d02" in report["suggested_gate_commands"]
