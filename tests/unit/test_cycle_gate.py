import json
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.cycle_gate import evaluate_cycles, main as gate_main, report_to_json


CORE_VARIABLES = ("U", "V", "T", "PH", "MU", "P", "QVAPOR")
START = "2025-07-26_00:00:00"
END = "2025-07-26_06:00:00"
END_0010 = "2025-07-26_00:10:00"
END_0020 = "2025-07-26_00:20:00"
END_FILE = "wrfout_d02_2025-07-26_06:00:00"
END_0010_FILE = "wrfout_d02_2025-07-26_00:10:00"
END_0020_FILE = "wrfout_d02_2025-07-26_00:20:00"
SELECTED_FIELD_INTEGRATOR_KIND = "selected_field_integrator_v0"
SELECTED_FIELD_DIAGNOSTIC_ADAPTER_KIND = "selected_field_diagnostic_adapter_v0"
PRODUCTION_CANDIDATE_ATTRS = {
    "TYWRF_GATE_CANDIDATE": "true",
    "TYWRF_INTEGRATOR_OUTPUT": "true",
    "TYWRF_VALIDATION_GATE_ONLY": "false",
    "TYWRF_CANDIDATE_KIND": "integrator_candidate",
}
WIND_TENDENCY_SUBCYCLING_ATTRS = {
    "TYWRF_WIND_TENDENCY_SUBSTEP_COUNT": 75,
    "TYWRF_WIND_TENDENCY_SUBSTEP_DT_SECONDS": 8,
    "TYWRF_WIND_TENDENCY_TOTAL_SECONDS": 600,
}
WIND_TENDENCY_ADVECTING_VELOCITY_MODE_ATTR = (
    "TYWRF_WIND_TENDENCY_ADVECTING_VELOCITY_MODE"
)
WIND_TENDENCY_ADVECTING_COMPONENTS_ATTR = (
    "TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS"
)
WIND_TENDENCY_ADVECTING_COLLOCATION_ATTR = (
    "TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION"
)


def _wind_tendency_advecting_collocation(advecting_components: str) -> str:
    return "average" if advecting_components == "cross_component" else "same_grid"


def _production_attrs(overrides: dict[str, object] | None = None) -> dict[str, object]:
    attrs = dict(PRODUCTION_CANDIDATE_ATTRS)
    attrs.update(overrides or {})
    return attrs


def _wind_tendency_evidence_attrs(
    overrides: dict[str, object] | None = None,
) -> dict[str, object]:
    overrides = overrides or {}
    attrs = {
        "TYWRF_WIND_TENDENCY_OPT_IN": "true",
        "TYWRF_WIND_TENDENCY_APPLIED": "true",
        "TYWRF_WIND_TENDENCY_SOURCE_KIND": "self_advection",
        "TYWRF_WIND_TENDENCY_GATE_EVIDENCE": "true",
        "TYWRF_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE": "true",
        "TYWRF_WIND_TENDENCY_USES_REFERENCE_END_TRUTH": "false",
        "TYWRF_WIND_TENDENCY_ZERO_OR_IDENTITY_ONLY": "false",
        WIND_TENDENCY_ADVECTING_VELOCITY_MODE_ATTR: "refreshed",
        WIND_TENDENCY_ADVECTING_COMPONENTS_ATTR: "same_component",
        **WIND_TENDENCY_SUBCYCLING_ATTRS,
    }
    attrs.update(overrides)
    if WIND_TENDENCY_ADVECTING_COLLOCATION_ATTR not in overrides:
        attrs[WIND_TENDENCY_ADVECTING_COLLOCATION_ATTR] = (
            _wind_tendency_advecting_collocation(
                str(attrs[WIND_TENDENCY_ADVECTING_COMPONENTS_ATTR])
            )
        )
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


def _write_10min_progressive_pair(
    tmp_path: Path,
    *,
    candidate_attrs: dict[str, object],
) -> tuple[Path, Path]:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    for file_name in (END_0010_FILE, END_0020_FILE):
        _write_wrfout(reference_dir / file_name)
        _write_wrfout(candidate_dir / file_name, attrs=candidate_attrs)
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


def test_cycle_gate_selected_field_integrator_kind_reaches_numeric_checks(
    tmp_path: Path,
) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        field_offset=1.0,
        attrs={"TYWRF_CANDIDATE_KIND": SELECTED_FIELD_INTEGRATOR_KIND},
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    diagnostics = {metric.name: metric for metric in report.cycles[0].diagnostics}
    fields = {field.variable: field for field in report.cycles[0].fields}

    assert report.status == "failed"
    assert diagnostics["candidate_metadata"].status == "passed"
    assert fields["U"].status == "failed"
    assert fields["U"].source_status == "threshold_exceeded"


def test_cycle_gate_selected_field_integrator_kind_requires_positive_metadata(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _write_wrfout(reference_dir / END_FILE)
    _write_wrfout(
        candidate_dir / END_FILE,
        attrs={"TYWRF_CANDIDATE_KIND": SELECTED_FIELD_INTEGRATOR_KIND},
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert "TYWRF_GATE_CANDIDATE is not true" in (metadata.message or "")
    assert "TYWRF_INTEGRATOR_OUTPUT is not true" in (metadata.message or "")


@pytest.mark.parametrize(
    "candidate_kind",
    (
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_oracle",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_diagnostic",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_remap",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_closure",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_helper",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_probe",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_adapter",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_dry_run",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_staging",
        f"{SELECTED_FIELD_INTEGRATOR_KIND}_experimental",
    ),
)
def test_cycle_gate_selected_field_integrator_rejects_artifact_kind_tokens(
    tmp_path: Path,
    candidate_kind: str,
) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={"TYWRF_CANDIDATE_KIND": candidate_kind},
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert f"TYWRF_CANDIDATE_KIND={candidate_kind}" in (metadata.message or "")


@pytest.mark.parametrize(
    ("metadata_name", "metadata_value"),
    (
        ("TYWRF_PRESSURE_REFRESH_HELPER_NAME", "pressure_refresh_helper_v0"),
        ("TYWRF_PRESSURE_COLUMN_PROBE_SOURCE", "pressure_column_probe"),
        ("TYWRF_OUTPUT_ADAPTER", "netcdf_adapter_v0"),
        ("TYWRF_RUN_MODE", "dry_run"),
        ("TYWRF_PRESSURE_REFRESH_STAGING_NAME", "pressure_refresh_staging_v0"),
        (
            "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS",
            "experimental_apply_test_only",
        ),
    ),
)
def test_cycle_gate_rejects_related_artifact_metadata_tokens(
    tmp_path: Path,
    metadata_name: str,
    metadata_value: str,
) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={
            "TYWRF_CANDIDATE_KIND": SELECTED_FIELD_INTEGRATOR_KIND,
            metadata_name: metadata_value,
        },
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert f"{metadata_name}={metadata_value}" in (metadata.message or "")
    assert {field.status for field in report.cycles[0].fields} == {"passed"}


def test_cycle_gate_allows_pressure_refresh_production_readiness_metadata(
    tmp_path: Path,
) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={
            "TYWRF_CANDIDATE_KIND": SELECTED_FIELD_INTEGRATOR_KIND,
            "TYWRF_PRESSURE_REFRESH_OPT_IN": "true",
            "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
            "TYWRF_PRESSURE_REFRESH_SOURCE_SYNC_OK": "true",
            "TYWRF_PRESSURE_REFRESH_PRODUCTION_SOURCE": (
                "krosa_moving_nest_pressure_refresh"
            ),
            "TYWRF_BASE_STATE_SOURCE_SYNC_READINESS_CHECK": "true",
            "TYWRF_BASE_STATE_SOURCE_SYNC_APPLIED": "false",
            "TYWRF_SOURCE_SYNC_PLANNED_PB_POINT_COUNT": 25,
            "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK": "true",
            "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK_CALLED": "true",
            "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK_OK": "true",
            "TYWRF_PRESSURE_REFRESH_PLANNED_P_POINT_COUNT": 25,
            "TYWRF_PRESSURE_READINESS_INVALID_P_POINT_COUNT": 0,
            "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_TARGET_COLUMN_COUNT": 25,
            "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS": (
                ":pressure_refresh_readiness("
                "base_state_source_sync_readiness_check=true,"
                "source_sync_planned_pb_points=25,"
                "readiness_invalid_p_points=0)"
            ),
        },
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "passed"
    assert metadata.status == "passed"


def test_cycle_gate_rejects_selected_field_diagnostic_adapter_ready_counts(
    tmp_path: Path,
) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={
            "TYWRF_CANDIDATE_KIND": SELECTED_FIELD_DIAGNOSTIC_ADAPTER_KIND,
            "TYWRF_SELECTED_FIELD_DIAGNOSTIC_ADAPTER": (
                SELECTED_FIELD_DIAGNOSTIC_ADAPTER_KIND
            ),
            "TYWRF_SELECTED_FIELD_DIAGNOSTIC_ADAPTER_READY": "true",
            "TYWRF_SELECTED_FIELD_DIAGNOSTIC_ADAPTER_READY_COUNT": 25,
            "TYWRF_SELECTED_FIELD_DIAGNOSTIC_ADAPTER_TOTAL_COUNT": 25,
            "TYWRF_SELECTED_FIELD_STAGING_READY": "true",
            "TYWRF_SELECTED_FIELD_STAGING_COLUMN_COUNT": 25,
        },
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    diagnostics = {metric.name: metric for metric in report.cycles[0].diagnostics}
    metadata = diagnostics["candidate_metadata"]

    assert report.status == "failed"
    assert metadata.status == "failed"
    assert f"TYWRF_CANDIDATE_KIND={SELECTED_FIELD_DIAGNOSTIC_ADAPTER_KIND}" in (
        metadata.message or ""
    )
    assert (
        f"TYWRF_SELECTED_FIELD_DIAGNOSTIC_ADAPTER={SELECTED_FIELD_DIAGNOSTIC_ADAPTER_KIND}"
        in (metadata.message or "")
    )
    assert "TYWRF_SELECTED_FIELD_DIAGNOSTIC_ADAPTER_READY=true" in (
        metadata.message or ""
    )
    assert "TYWRF_SELECTED_FIELD_STAGING_READY=true" in (metadata.message or "")
    assert {field.status for field in report.cycles[0].fields} == {"passed"}
    assert {
        metric.status
        for name, metric in diagnostics.items()
        if name != "candidate_metadata"
    } == {"passed"}


def test_cycle_gate_allows_false_related_metadata_flags_for_integrator_candidate(
    tmp_path: Path,
) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs={
            "TYWRF_CANDIDATE_KIND": SELECTED_FIELD_INTEGRATOR_KIND,
            "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY": "false",
            "TYWRF_PRESSURE_COLUMN_PROBE_ENABLED": "false",
        },
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)

    assert report.status == "passed"


@pytest.mark.parametrize(
    ("wind_attrs", "expected_message"),
    (
        (
            {"TYWRF_WIND_TENDENCY_SOURCE_KIND": "zero_tendency"},
            "TYWRF_WIND_TENDENCY_SOURCE_KIND=zero_tendency",
        ),
        (
            {"TYWRF_WIND_TENDENCY_SOURCE_KIND": "identity_update"},
            "TYWRF_WIND_TENDENCY_SOURCE_KIND=identity_update",
        ),
        (
            {"TYWRF_WIND_TENDENCY_SOURCE_KIND": "oracle_reference_end_truth"},
            "TYWRF_WIND_TENDENCY_SOURCE_KIND=oracle_reference_end_truth",
        ),
        (
            {"TYWRF_WIND_TENDENCY_SOURCE_KIND": "reference-end truth"},
            "TYWRF_WIND_TENDENCY_SOURCE_KIND=reference-end truth",
        ),
        (
            {"TYWRF_WIND_TENDENCY_GATE_EVIDENCE": "false"},
            "TYWRF_WIND_TENDENCY_GATE_EVIDENCE=false",
        ),
        (
            {"TYWRF_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE": "false"},
            "TYWRF_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE=false",
        ),
        (
            {"TYWRF_WIND_TENDENCY_USES_REFERENCE_END_TRUTH": "true"},
            "TYWRF_WIND_TENDENCY_USES_REFERENCE_END_TRUTH=true",
        ),
        (
            {"TYWRF_WIND_TENDENCY_ZERO_OR_IDENTITY_ONLY": "true"},
            "TYWRF_WIND_TENDENCY_ZERO_OR_IDENTITY_ONLY=true",
        ),
    ),
)
def test_cycle_gate_rejects_wind_tendency_placeholder_metadata_at_0010(
    tmp_path: Path,
    wind_attrs: dict[str, object],
    expected_message: str,
) -> None:
    attrs = _production_attrs(_wind_tendency_evidence_attrs(wind_attrs))
    reference_dir, candidate_dir = _write_10min_progressive_pair(
        tmp_path,
        candidate_attrs=attrs,
    )

    report = evaluate_cycles(
        reference_dir,
        candidate_dir,
        START,
        end=END_0020,
        interval_minutes=10,
    )
    payload = json.loads(report_to_json(report))
    cycle = report.cycles[0]
    diagnostics = {metric.name: metric for metric in cycle.diagnostics}
    metadata = diagnostics["candidate_metadata"]

    assert report.status == "failed"
    assert report.summary == {"total": 1, "passed": 0, "failed": 1}
    assert len(report.cycles) == 1
    assert cycle.end_time == END_0010
    assert [item["end_time"] for item in payload["cycles"]] == [END_0010]
    assert metadata.status == "failed"
    assert expected_message in (metadata.message or "")
    assert {field.status for field in cycle.fields} == {"passed"}
    assert {
        metric.status
        for name, metric in diagnostics.items()
        if name != "candidate_metadata"
    } == {"passed"}
    assert payload["first_failure"]["cycle_index"] == 1
    assert payload["first_failure"]["end_time"] == END_0010
    assert payload["first_failure"]["field"] is None
    assert payload["first_failure"]["diagnostic"] == "candidate_metadata"


def test_cycle_gate_rejects_wind_tendency_non_evidence_even_when_validation_gate_only_allowed(
    tmp_path: Path,
) -> None:
    attrs = _production_attrs(
        _wind_tendency_evidence_attrs(
            {
                "TYWRF_WIND_TENDENCY_GATE_EVIDENCE": "false",
            }
        )
    )
    attrs.update(
        {
            "TYWRF_VALIDATION_GATE_ONLY": "true",
            "TYWRF_INTEGRATOR_OUTPUT": "false",
        }
    )
    reference_dir, candidate_dir = _write_10min_progressive_pair(
        tmp_path,
        candidate_attrs=attrs,
    )

    report = evaluate_cycles(
        reference_dir,
        candidate_dir,
        START,
        end=END_0020,
        interval_minutes=10,
        allow_validation_gate_only=True,
    )
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "failed"
    assert report.summary == {"total": 1, "passed": 0, "failed": 1}
    assert metadata.status == "failed"
    assert "TYWRF_WIND_TENDENCY_GATE_EVIDENCE=false" in (
        metadata.message or ""
    )
    assert "TYWRF_VALIDATION_GATE_ONLY=true" not in (metadata.message or "")


@pytest.mark.parametrize("advecting_velocity_mode", ("refreshed", "frozen"))
@pytest.mark.parametrize("advecting_components", ("same_component", "cross_component"))
def test_cycle_gate_allows_wind_tendency_evidence_metadata(
    tmp_path: Path,
    advecting_velocity_mode: str,
    advecting_components: str,
) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        attrs=_wind_tendency_evidence_attrs(
            {
                WIND_TENDENCY_ADVECTING_VELOCITY_MODE_ATTR: advecting_velocity_mode,
                WIND_TENDENCY_ADVECTING_COMPONENTS_ATTR: advecting_components,
            }
        ),
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    metadata = {metric.name: metric for metric in report.cycles[0].diagnostics}[
        "candidate_metadata"
    ]

    assert report.status == "passed"
    assert metadata.status == "passed"


@pytest.mark.parametrize("advecting_velocity_mode", ("refreshed", "frozen"))
@pytest.mark.parametrize("advecting_components", ("same_component", "cross_component"))
def test_cycle_gate_self_advection_subcycling_metadata_enters_normal_field_gate(
    tmp_path: Path,
    advecting_velocity_mode: str,
    advecting_components: str,
) -> None:
    reference_dir, candidate_dir = _write_pair(
        tmp_path,
        field_offset=1.0,
        attrs=_wind_tendency_evidence_attrs(
            {
                WIND_TENDENCY_ADVECTING_VELOCITY_MODE_ATTR: advecting_velocity_mode,
                WIND_TENDENCY_ADVECTING_COMPONENTS_ATTR: advecting_components,
            }
        ),
    )

    report = evaluate_cycles(reference_dir, candidate_dir, START, end=END)
    payload = json.loads(report_to_json(report))
    cycle = report.cycles[0]
    diagnostics = {metric.name: metric for metric in cycle.diagnostics}
    fields = {field.variable: field for field in cycle.fields}
    with netCDF4.Dataset(candidate_dir / END_FILE) as dataset:
        for name, expected_value in WIND_TENDENCY_SUBCYCLING_ATTRS.items():
            assert dataset.getncattr(name) == expected_value
        assert (
            dataset.getncattr(WIND_TENDENCY_ADVECTING_VELOCITY_MODE_ATTR)
            == advecting_velocity_mode
        )
        assert (
            dataset.getncattr(WIND_TENDENCY_ADVECTING_COMPONENTS_ATTR)
            == advecting_components
        )
        assert dataset.getncattr(WIND_TENDENCY_ADVECTING_COLLOCATION_ATTR) == (
            _wind_tendency_advecting_collocation(advecting_components)
        )

    assert report.status == "failed"
    assert cycle.status == "failed"
    assert diagnostics["candidate_metadata"].status == "passed"
    assert payload["cycles"][0]["diagnostics"][0]["name"] == "candidate_metadata"
    assert payload["cycles"][0]["diagnostics"][0]["status"] == "passed"
    assert fields["U"].status == "failed"
    assert fields["U"].source_status == "threshold_exceeded"
    assert {
        metric.status
        for name, metric in diagnostics.items()
        if name != "candidate_metadata"
    } == {"passed"}
    assert payload["first_failure"]["field"] == "U"
    assert payload["first_failure"]["diagnostic"] is None


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
