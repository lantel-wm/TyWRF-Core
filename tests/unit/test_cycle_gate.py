import json
from pathlib import Path

import netCDF4
import numpy as np

from tools.cycle_gate import evaluate_cycles, main as gate_main, report_to_json


CORE_VARIABLES = ("U", "V", "T", "PH", "MU", "P", "QVAPOR")
START = "2025-07-26_00:00:00"
END = "2025-07-26_06:00:00"
END_FILE = "wrfout_d02_2025-07-26_06:00:00"
PRODUCTION_CANDIDATE_ATTRS = {
    "TYWRF_GATE_CANDIDATE": "true",
    "TYWRF_INTEGRATOR_OUTPUT": "true",
    "TYWRF_VALIDATION_GATE_ONLY": "false",
    "TYWRF_CANDIDATE_KIND": "integrator_candidate",
}


def _production_attrs(overrides: dict[str, object] | None = None) -> dict[str, object]:
    attrs = dict(PRODUCTION_CANDIDATE_ATTRS)
    attrs.update(overrides or {})
    return attrs


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


def _write_wrfout(
    path: Path,
    *,
    field_offset: float = 0.0,
    omit_slp: bool = False,
    attrs: dict[str, object] | None = None,
    **tc_kwargs,
) -> None:
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
        for name, value in (attrs or {}).items():
            dataset.setncattr(name, value)
        for name, values in fields.items():
            variable = dataset.createVariable(name, "f8", ("Time", "south_north", "west_east"))
            variable[0, :, :] = values


def _write_pair(tmp_path: Path, **candidate_kwargs) -> tuple[Path, Path]:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    candidate_kwargs = dict(candidate_kwargs)
    candidate_attrs = _production_attrs(candidate_kwargs.pop("attrs", None))
    _write_wrfout(reference_dir / END_FILE)
    _write_wrfout(candidate_dir / END_FILE, attrs=candidate_attrs, **candidate_kwargs)
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


def test_cycle_gate_reports_first_failure_for_10_min_progressive_run(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    end_times = (
        "2025-07-26_00:10:00",
        "2025-07-26_00:20:00",
        "2025-07-26_00:30:00",
        "2025-07-26_00:40:00",
        "2025-07-26_00:50:00",
        "2025-07-26_01:00:00",
    )
    for end_time in end_times:
        _write_wrfout(reference_dir / f"wrfout_d02_{end_time}")
        _write_wrfout(
            candidate_dir / f"wrfout_d02_{end_time}",
            attrs=_production_attrs(),
            field_offset=1.0 if end_time == "2025-07-26_00:20:00" else 0.0,
        )

    report = evaluate_cycles(
        reference_dir,
        candidate_dir,
        START,
        end="2025-07-26_01:00:00",
        interval_minutes=10,
    )
    payload = json.loads(report_to_json(report))

    assert report.status == "failed"
    assert report.summary == {"total": 6, "passed": 5, "failed": 1}
    assert payload["interval_minutes"] == 10
    assert payload["cycles"][1]["end_time"] == "2025-07-26_00:20:00"
    assert payload["first_failure"]["cycle_index"] == 2
    assert payload["first_failure"]["end_time"] == "2025-07-26_00:20:00"
    assert payload["first_failure"]["field"] == "U"
    assert payload["first_failure"]["field_status"] == "failed"
    assert payload["first_failure"]["diagnostic"] is None


def test_cycle_gate_rejects_diagnostic_only_candidate_metadata(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={
            "TYWRF_DIAGNOSTIC_ONLY": "true",
            "TYWRF_CANDIDATE_KIND": "cpp_skeleton_remap_overlap_diagnostic",
        },
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    diagnostics = {metric.name: metric for metric in report.cycles[0].diagnostics}

    assert report.status == "failed"
    assert report.cycles[0].status == "failed"
    assert diagnostics["candidate_metadata"].status == "failed"
    assert "TYWRF_DIAGNOSTIC_ONLY=true" in diagnostics["candidate_metadata"].message
    assert "TYWRF_CANDIDATE_KIND=cpp_skeleton_remap_overlap_diagnostic" in (
        diagnostics["candidate_metadata"].message or ""
    )
    assert {field.status for field in report.cycles[0].fields} == {"passed"}


def test_cycle_gate_rejects_gate_candidate_false_metadata(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={
            "TYWRF_GATE_CANDIDATE": "false",
            "TYWRF_VALIDATION_GATE_ONLY": "false",
            "TYWRF_CANDIDATE_KIND": "closure_artifact",
        },
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert "TYWRF_GATE_CANDIDATE=false" in (metadata.message or "")
    assert "TYWRF_CANDIDATE_KIND=closure_artifact" in (metadata.message or "")


def test_cycle_gate_rejects_explicit_non_integrator_candidate_by_default(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={
            "TYWRF_CANDIDATE_KIND": "baseline_candidate",
            "TYWRF_INTEGRATOR_OUTPUT": "false",
        },
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert "TYWRF_INTEGRATOR_OUTPUT=false" in (metadata.message or "")


def test_cycle_gate_rejects_remap_candidate_kind(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={"TYWRF_CANDIDATE_KIND": "parent_remap_candidate"},
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert "TYWRF_CANDIDATE_KIND=parent_remap_candidate" in (metadata.message or "")


def test_cycle_gate_rejects_oracle_candidate_kind(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={"TYWRF_CANDIDATE_KIND": "reference_delta_oracle_candidate"},
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert "TYWRF_CANDIDATE_KIND=reference_delta_oracle_candidate" in (
        metadata.message or ""
    )


def test_cycle_gate_rejects_missing_positive_candidate_metadata_in_strict_path(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _write_wrfout(reference_dir / END_FILE)
    _write_wrfout(candidate_dir / END_FILE)

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert "TYWRF_GATE_CANDIDATE is not true" in (metadata.message or "")
    assert "TYWRF_INTEGRATOR_OUTPUT is not true" in (metadata.message or "")


def test_cycle_gate_allows_reference_copy_only_when_explicitly_requested(tmp_path: Path) -> None:
    attrs = {
        "TYWRF_CANDIDATE_KIND": "baseline_candidate",
        "TYWRF_REFERENCE_COPY": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "false",
        "TYWRF_VALIDATION_GATE_ONLY": "true",
    }
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _write_wrfout(reference_dir / END_FILE)
    _write_wrfout(candidate_dir / END_FILE, attrs=attrs)

    default_report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    allowed_report = evaluate_cycles(
        reference_dir,
        candidate_dir,
        START,
        end=END,
        allow_validation_gate_only=True,
    )

    default_metadata = {metric.name: metric for metric in default_report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]
    allowed_metadata = {metric.name: metric for metric in allowed_report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]
    assert default_report.status == "failed"
    assert "TYWRF_VALIDATION_GATE_ONLY=true" in (default_metadata.message or "")
    assert allowed_report.status == "passed"
    assert allowed_metadata.status == "passed"


def test_cycle_gate_legacy_flag_does_not_allow_unmarked_candidate(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _write_wrfout(reference_dir / END_FILE)
    _write_wrfout(candidate_dir / END_FILE)

    report = evaluate_cycles(
        reference_dir,
        candidate_dir,
        START,
        end=END,
        allow_validation_gate_only=True,
    )
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert "TYWRF_GATE_CANDIDATE is not true" in (metadata.message or "")


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
