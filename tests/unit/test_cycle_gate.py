import json
from pathlib import Path

import netCDF4
import numpy as np

from tools.cycle_gate import evaluate_cycles, main as gate_main, report_to_json


CORE_VARIABLES = ("U", "V", "T", "PH", "MU", "P", "QVAPOR")
START = "2025-07-26_00:00:00"
END = "2025-07-26_06:00:00"
END_FILE = "wrfout_d02_2025-07-26_06:00:00"


def _tc_fields(*, center_shift: bool = False, slp_offset_hpa: float = 0.0, vmax_offset: float = 0.0):
    shape = (5, 5)
    latitude = np.repeat(np.arange(10.0, 15.0)[:, None], shape[1], axis=1)
    longitude = np.repeat(np.arange(120.0, 125.0)[None, :], shape[0], axis=0)
    slp = np.full(shape, 1000.0)
    center = (2, 2)
    if center_shift:
        center = (2, 4)
    slp[center] = 950.0 + slp_offset_hpa
    u10 = np.zeros(shape)
    v10 = np.zeros(shape)
    u10[1, 1] = 30.0 + vmax_offset
    v10[1, 1] = 40.0
    return {
        "XLAT": latitude,
        "XLONG": longitude,
        "SLP": slp,
        "U10": u10,
        "V10": v10,
    }


def _write_wrfout(path: Path, *, field_offset: float = 0.0, omit_slp: bool = False, **tc_kwargs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    shape = (5, 5)
    fields = {name: np.full(shape, 10.0 + field_offset) for name in CORE_VARIABLES}
    fields.update(_tc_fields(**tc_kwargs))
    if omit_slp:
        fields.pop("SLP")

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", 1)
        dataset.createDimension("south_north", shape[0])
        dataset.createDimension("west_east", shape[1])
        for name, values in fields.items():
            variable = dataset.createVariable(name, "f8", ("Time", "south_north", "west_east"))
            variable[0, :, :] = values


def _write_pair(tmp_path: Path, **candidate_kwargs) -> tuple[Path, Path]:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _write_wrfout(reference_dir / END_FILE)
    _write_wrfout(candidate_dir / END_FILE, **candidate_kwargs)
    return reference_dir, candidate_dir


def test_cycle_gate_passes_matching_d02_cycle(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(tmp_path)

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    payload = json.loads(report_to_json(report))

    assert report.status == "passed"
    assert report.summary == {"total": 1, "passed": 1, "failed": 0}
    assert payload["domain"] == "d02"
    assert {field["status"] for field in payload["cycles"][0]["fields"]} == {"passed"}
    assert {metric["status"] for metric in payload["cycles"][0]["diagnostics"]} == {"passed"}


def test_cycle_gate_fails_field_normalized_rmse_threshold(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(tmp_path, field_offset=1.0)

    report = evaluate_cycles(reference_dir, candidate_dir, START, hours=6)
    fields = {field.variable: field for field in report.cycles[0].fields}

    assert report.status == "failed"
    assert report.cycles[0].status == "failed"
    assert fields["U"].status == "failed"
    assert fields["U"].source_status == "threshold_exceeded"
    assert fields["U"].normalized_rmse > 0.05


def test_cycle_gate_marks_missing_candidate_file_not_available(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _write_wrfout(reference_dir / END_FILE)
    candidate_dir.mkdir()

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)

    assert report.status == "failed"
    assert "missing candidate file" in report.cycles[0].message
    assert {field.status for field in report.cycles[0].fields} == {"not_available"}
    assert {metric.status for metric in report.cycles[0].diagnostics} == {"not_available"}


def test_cycle_gate_fails_tc_diagnostic_thresholds(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        center_shift=True,
        slp_offset_hpa=6.0,
        vmax_offset=10.0,
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    diagnostics = {metric.name: metric for metric in report.cycles[0].diagnostics}

    assert report.status == "failed"
    assert diagnostics["storm_center"].status == "failed"
    assert diagnostics["minimum_slp"].status == "failed"
    assert diagnostics["vmax10m"].status == "failed"


def test_cycle_gate_marks_missing_slp_diagnostic_not_available(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(tmp_path, omit_slp=True)

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    diagnostics = {metric.name: metric for metric in report.cycles[0].diagnostics}

    assert report.status == "failed"
    assert diagnostics["minimum_slp"].status == "not_available"
    assert "minimum SLP diagnostics failed" in diagnostics["minimum_slp"].message


def test_cycle_gate_cli_writes_json(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(tmp_path)
    output = tmp_path / "gate.json"

    exit_code = gate_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--candidate-dir",
            str(candidate_dir),
            "--start",
            START,
            "--end",
            END,
            "--output",
            str(output),
            "--pretty",
        ]
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert exit_code == 0
    assert payload["status"] == "passed"
    assert payload["cycles"][0]["end_time"] == END
