import json
import os
from pathlib import Path
import subprocess

import netCDF4

from tools.audit_selected_field_timeline import (
    audit_selected_field_timeline,
    main as audit_main,
    report_to_json,
)


def _ensure_dim(dataset: netCDF4.Dataset, name: str, size: int | None) -> None:
    if name in dataset.dimensions:
        if size is not None:
            assert len(dataset.dimensions[name]) == size
        return
    dataset.createDimension(name, size)


def _write_candidate(
    path: Path,
    *,
    attrs: dict[str, object],
    mass_shape: tuple[int, int] = (210, 210),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        for key, value in attrs.items():
            dataset.setncattr(key, value)

        mass_ny, mass_nx = mass_shape
        _ensure_dim(dataset, "bottom_top", 2)
        _ensure_dim(dataset, "bottom_top_stag", 3)
        _ensure_dim(dataset, "south_north", mass_ny)
        _ensure_dim(dataset, "west_east", mass_nx)
        _ensure_dim(dataset, "south_north_stag", mass_ny + 1)
        _ensure_dim(dataset, "west_east_stag", mass_nx + 1)
        dataset.createVariable(
            "U",
            "f4",
            ("Time", "bottom_top", "south_north", "west_east_stag"),
        )
        dataset.createVariable(
            "V",
            "f4",
            ("Time", "bottom_top", "south_north_stag", "west_east"),
        )
        for name in ("T", "P", "QVAPOR"):
            dataset.createVariable(
                name,
                "f4",
                ("Time", "bottom_top", "south_north", "west_east"),
            )
        dataset.createVariable(
            "PH",
            "f4",
            ("Time", "bottom_top_stag", "south_north", "west_east"),
        )
        dataset.createVariable("MU", "f4", ("Time", "south_north", "west_east"))


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _base_attrs(
    *,
    from_start: str = "1,1",
    to_start: str = "13,8",
    ratio: int = 5,
    changed: float = 1000.0,
    interpolated: float = 900.0,
) -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "false",
        "TYWRF_GATE_CANDIDATE": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
        "TYWRF_CANDIDATE_DOMAIN": "d02",
        "TYWRF_D02_RESOLUTION_CHECK": "d02_2km",
        "TYWRF_CYCLE_START": "2025-07-26_00:00:00",
        "TYWRF_CYCLE_END": "2025-07-26_00:10:00",
        "TYWRF_FROM_PARENT_START": from_start,
        "TYWRF_TO_PARENT_START": to_start,
        "TYWRF_PARENT_GRID_RATIO": ratio,
        "TYWRF_STATE_VARIABLES": "U,V,T,PH,MU,P,QVAPOR",
        "TYWRF_PARENT_INTERPOLATED_STATE_VARIABLES": "U,V,T,PH,MU,QVAPOR",
        "TYWRF_SELECTED_FIELD_CHANGED_POINTS": changed,
        "TYWRF_EXPOSED_EXCHANGE_POINTS": 900.0,
        "TYWRF_INTERPOLATED_POINTS": interpolated,
        "TYWRF_STATIC_REFRESH_APPLIED": "true",
        "TYWRF_STATIC_REFRESH_METHOD": "overlap_shift_xlat_xlong_extrapolate_hgt_parent_bilinear_v0",
        "TYWRF_STATIC_REFRESH_USES_REFERENCE_END": "false",
        "TYWRF_STATIC_REFRESH_OVERLAP_CELLS": 26250.0,
        "TYWRF_STATIC_REFRESH_EXPOSED_CELLS": 17850.0,
        "TYWRF_PRESSURE_REFRESH_APPLIED": "false",
        "DX": 2000.0,
        "DY": 2000.0,
    }


def _runtime_event_names(*, include_pressure_column_probe: bool = False) -> str:
    names = [
        "cycle_start",
        "move_from_to_parent_start",
        "overlap_remap",
        "exchange_plan_build",
        "parent_interpolation",
        "selected_field_change_summary",
        "static_refresh",
        "pressure_refresh_readiness",
        "pressure_refresh_apply",
        "cycle_end",
        "output_write_preparation",
    ]
    if include_pressure_column_probe:
        names.insert(9, "pressure_column_probe")
    return ",".join(names)


def _runtime_events(
    *,
    include_pressure_column_probe: bool = False,
    exchange_points: int = 900,
    interpolated_points: int = 900,
    static_overlap_cells: int = 26250,
    static_exposed_cells: int = 17850,
    static_hgt_parent_interpolated_cells: int = 17850,
    static_changed_template_points: int = 53000,
) -> str:
    events = [
        "1:cycle_start(cycle_start=2025-07-26_00:00:00,cycle_end=2025-07-26_00:10:00)",
        (
            "2:move_from_to_parent_start(from_i=1,from_j=1,to_i=13,to_j=8,"
            "parent_grid_ratio=5,parent_delta_i=12,parent_delta_j=7,"
            "child_delta_i=60,child_delta_j=35)"
        ),
        (
            "3:overlap_remap(copied_points=1200,copied_fields=7,"
            "child_delta_i=60,child_delta_j=35,needs_parent_fill=true,"
            "needs_derived_pressure_refresh=true)"
        ),
        (
            "4:exchange_plan_build(exchange_points="
            f"{exchange_points},requires_parent_interpolation=true,"
            "modifies_overlap=false,modifies_halo=false)"
        ),
        (
            "5:parent_interpolation(interpolated_points="
            f"{interpolated_points},wrote_overlap=false,wrote_halo=false)"
        ),
        "6:selected_field_change_summary(changed_points=900)",
        (
            "7:static_refresh(overlap_cells="
            f"{static_overlap_cells},exposed_cells={static_exposed_cells},"
            "coord_extrapolated_cells=200,"
            f"hgt_parent_interpolated_cells={static_hgt_parent_interpolated_cells},"
            f"changed_template_points={static_changed_template_points},"
            "uses_reference_end=false)"
        ),
        "8:pressure_refresh_readiness(opt_in=false,ready=not_applicable,status=skipped)",
        "9:pressure_refresh_apply(opt_in=false,applied=false,status=skipped)",
    ]
    if include_pressure_column_probe:
        events.append(
            "10:pressure_column_probe(enabled=true,phase_count=2,"
            "record_count=50,evidence_only=true)"
        )
    cycle_end_index = len(events) + 1
    events.extend(
        [
            (
                f"{cycle_end_index}:cycle_end("
                "cycle_start=2025-07-26_00:00:00,"
                "cycle_end=2025-07-26_00:10:00)"
            ),
            (
                f"{cycle_end_index + 1}:output_write_preparation("
                "time_index=0,variable_count=7,times=2025-07-26_00:10:00,"
                "state_write=pending,metadata_write=pending)"
            ),
        ]
    )
    return "|".join(events)


def _runtime_attrs(
    *,
    include_pressure_column_probe: bool = False,
    **overrides: object,
) -> dict[str, object]:
    attrs: dict[str, object] = {
        "TYWRF_SELECTED_FIELD_TIMELINE_VERSION": "runtime_v0",
        "TYWRF_SELECTED_FIELD_TIMELINE_EVIDENCE_ONLY": "true",
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_COUNT": (
            12 if include_pressure_column_probe else 11
        ),
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_NAMES": _runtime_event_names(
            include_pressure_column_probe=include_pressure_column_probe
        ),
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS": _runtime_events(
            include_pressure_column_probe=include_pressure_column_probe
        ),
        "TYWRF_STATIC_REFRESH_HGT_PARENT_INTERPOLATED_CELLS": 17850.0,
        "TYWRF_STATIC_REFRESH_CHANGED_TEMPLATE_POINTS": 53000.0,
    }
    attrs.update(overrides)
    return attrs


def _write_source_snippets(source_root: Path, *, include_post_fill: bool = True) -> None:
    selected = source_root / "src/tools/selected_field_cycle.cpp"
    selected.parent.mkdir(parents=True, exist_ok=True)
    selected.write_text(
        """
void build_candidate_state() {
  remap_child_state_overlap_only();
  build_exposed_child_state_exchange_plan();
  interpolate_parent_to_exposed_child();
  refresh_static_fields();
  apply_pressure_refresh();
  stamp_gate_metadata();
  "TYWRF_INTERPOLATED_POINTS";
}
""",
        encoding="utf-8",
    )
    schedule = source_root / "src/dynamics/cycle_schedule.cpp"
    schedule.parent.mkdir(parents=True, exist_ok=True)
    post_fill = (
        "  CycleScheduleCallKind::moving_nest_post_move_parent_fill;\n"
        if include_post_fill
        else ""
    )
    schedule.write_text(
        f"""
void schedule() {{
  CycleScheduleCallKind::parent_child_interpolation;
  CycleScheduleCallKind::moving_nest_move_check;
{post_fill}  snap_to_next_parent_step();
}}
""",
        encoding="utf-8",
    )
    header = source_root / "include/tywrf/dynamics/cycle_schedule.hpp"
    header.parent.mkdir(parents=True, exist_ok=True)
    header.write_text(
        """
enum class CycleScheduleCallKind {
  parent_child_interpolation,
  moving_nest_move_check,
  moving_nest_post_move_parent_fill,
};
""",
        encoding="utf-8",
    )
    exchange = source_root / "src/nest/state_exchange.cpp"
    exchange.parent.mkdir(parents=True, exist_ok=True)
    exchange.write_text(
        """
void build_exposed_child_state_exchange_plan() {
  StateExchangeField::u;
  StateExchangeField::v;
  StateExchangeField::mu;
  StateExchangeField::qvapor;
  StateExchangeField::t;
  StateExchangeField::ph;
}
""",
        encoding="utf-8",
    )


def test_runtime_timeline_attrs_are_parsed_and_serialized(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    attrs = {
        **_base_attrs(changed=900.0, interpolated=900.0),
        **_runtime_attrs(),
    }
    _write_candidate(candidate, attrs=attrs)
    _write_source_snippets(source_root)

    payload = json.loads(
        report_to_json(
            audit_selected_field_timeline(candidate, source_root=source_root)
        )
    )

    runtime = payload["candidate_timeline"]["runtime_timeline"]
    assert runtime["diagnostic_only"] is True
    assert runtime["candidate_model_pass"] == "not_applicable"
    assert runtime["status"] == "available"
    assert runtime["version"] == "runtime_v0"
    assert runtime["evidence_only"] is True
    assert runtime["parsed_event_count"] == 11
    assert runtime["event_count_consistent"] is True
    assert runtime["event_names_consistent"] is True
    assert runtime["expected_order_match"] is True
    assert runtime["parsed_events"]["events"][0]["name"] == "cycle_start"
    assert runtime["parsed_events"]["events"][1]["fields"]["child_delta_i"] == "60"
    assert runtime["parsed_events"]["events_by_name"]["static_refresh"]["fields"][
        "hgt_parent_interpolated_cells"
    ] == "17850"
    assert runtime["key_count_parity"]["all_available_counts_consistent"] is True
    assert runtime["key_count_parity"]["mismatch_count"] == 0
    assert "runtime_timeline_event_count_mismatch" not in set(
        payload["summary"]["risk_codes"]
    )
    assert "runtime_timeline_event_order_mismatch" not in set(
        payload["summary"]["risk_codes"]
    )
    assert "runtime_count_parity_mismatch" not in set(payload["summary"]["risk_codes"])


def test_runtime_timeline_pressure_column_probe_order_is_allowed(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    attrs = {
        **_base_attrs(changed=900.0, interpolated=900.0),
        **_runtime_attrs(include_pressure_column_probe=True),
    }
    _write_candidate(candidate, attrs=attrs)
    _write_source_snippets(source_root)

    payload = audit_selected_field_timeline(candidate, source_root=source_root)
    runtime = payload["candidate_timeline"]["runtime_timeline"]
    risk_codes = set(payload["summary"]["risk_codes"])

    assert runtime["diagnostic_only"] is True
    assert runtime["candidate_model_pass"] == "not_applicable"
    assert runtime["evidence_only"] is True
    assert runtime["parsed_event_count"] == 12
    assert runtime["parsed_event_names"][9] == "pressure_column_probe"
    assert runtime["parsed_events"]["events_by_name"]["pressure_column_probe"][
        "fields"
    ]["evidence_only"] == "true"
    assert runtime["event_count_consistent"] is True
    assert runtime["event_names_consistent"] is True
    assert runtime["expected_order_match"] is True
    assert "runtime_timeline_event_count_mismatch" not in risk_codes
    assert "runtime_timeline_event_order_mismatch" not in risk_codes


def test_runtime_timeline_malformed_string_is_reported(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    attrs = {
        **_base_attrs(changed=900.0, interpolated=900.0),
        **_runtime_attrs(
            TYWRF_SELECTED_FIELD_TIMELINE_EVENTS="1:cycle_start(cycle_start=broken"
        ),
    }
    _write_candidate(candidate, attrs=attrs)
    _write_source_snippets(source_root)

    payload = audit_selected_field_timeline(candidate, source_root=source_root)
    runtime = payload["candidate_timeline"]["runtime_timeline"]

    assert runtime["status"] == "malformed"
    assert runtime["parsed_events"]["malformed"] is True
    assert runtime["parsed_events"]["parse_errors"]
    assert "malformed_runtime_timeline_string" in set(payload["summary"]["risk_codes"])


def test_runtime_timeline_count_and_parity_mismatches_are_flagged(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    attrs = {
        **_base_attrs(changed=900.0, interpolated=901.0),
        **_runtime_attrs(TYWRF_SELECTED_FIELD_TIMELINE_EVENT_COUNT=12),
    }
    _write_candidate(candidate, attrs=attrs)
    _write_source_snippets(source_root)

    payload = audit_selected_field_timeline(candidate, source_root=source_root)
    runtime = payload["candidate_timeline"]["runtime_timeline"]
    risk_codes = set(payload["summary"]["risk_codes"])

    assert runtime["event_count_consistent"] is False
    assert runtime["key_count_parity"]["mismatch_count"] == 1
    assert runtime["key_count_parity"]["mismatches"][0]["label"] == "interpolated_points"
    assert "runtime_timeline_event_count_mismatch" in risk_codes
    assert "runtime_count_parity_mismatch" in risk_codes


def test_runtime_timeline_missing_attrs_are_nonfatal_but_flagged(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    _write_candidate(candidate, attrs=_base_attrs(changed=900.0, interpolated=900.0))
    _write_source_snippets(source_root)

    payload = audit_selected_field_timeline(candidate, source_root=source_root)
    runtime = payload["candidate_timeline"]["runtime_timeline"]

    assert runtime["status"] == "not_available"
    assert set(runtime["missing_required_attrs"]) == {
        "TYWRF_SELECTED_FIELD_TIMELINE_VERSION",
        "TYWRF_SELECTED_FIELD_TIMELINE_EVIDENCE_ONLY",
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_COUNT",
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_NAMES",
        "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS",
    }
    assert "missing_runtime_timeline_attrs_on_new_candidate" in set(
        payload["summary"]["risk_codes"]
    )


def test_timeline_audit_flags_large_delta_count_mismatch_and_prior_failures(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    pipeline_json = tmp_path / "pipeline.json"
    spatial_json = tmp_path / "spatial.json"
    evolution_json = tmp_path / "evolution.json"
    _write_candidate(candidate, attrs=_base_attrs())
    _write_source_snippets(source_root)
    _write_json(
        pipeline_json,
        {
            "status": "computed_with_flags",
            "candidate_model_pass": "not_applicable",
            "summary": {
                "risk_codes": [
                    "large_movement_delta_schedule_remap_risk",
                    "candidate_claims_gate_or_integrator_but_audits_fail",
                ]
            },
            "movement_geometry": {
                "movement": {
                    "child_delta": {"i": 60, "j": 35},
                    "large_movement": True,
                }
            },
        },
    )
    _write_json(
        spatial_json,
        {
            "status": "computed",
            "candidate_model_pass": "not_applicable",
            "variables": [
                {
                    "variable": "U",
                    "status": "computed",
                    "end_state": {
                        "global": {
                            "best_shift": {"di": 1, "dj": 0},
                            "baseline": {"normalized_rmse": 0.08},
                            "best": {"normalized_rmse": 0.07},
                            "normalized_rmse_reduction_fraction": 0.125,
                        }
                    },
                    "evolution": {
                        "global": {
                            "best_shift": {"di": 1, "dj": 1},
                            "baseline": {"normalized_rmse": 0.09},
                            "best": {"normalized_rmse": 0.08},
                            "normalized_rmse_reduction_fraction": 0.111,
                        }
                    },
                }
            ],
        },
    )
    _write_json(
        evolution_json,
        {
            "status": "computed",
            "candidate_model_pass": "not_applicable",
            "variables": [
                {
                    "variable": "U",
                    "status": "computed",
                    "normalized_rmse": 0.08,
                    "evolution_amplitude": {"amplitude_ratio": 1.02},
                    "region_split": {"target_error_fraction": 0.4},
                }
            ],
        },
    )

    payload = json.loads(
        report_to_json(
            audit_selected_field_timeline(
                candidate,
                pipeline_json_path=pipeline_json,
                spatial_alignment_json_path=spatial_json,
                evolution_json_path=evolution_json,
                source_root=source_root,
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["candidate_timeline"]["movement"]["child_delta"] == {
        "i": 60,
        "j": 35,
    }
    assert payload["candidate_timeline"]["movement"]["large_movement"] is True
    assert payload["candidate_timeline"]["counts"][
        "changed_interpolated_point_count_mismatch"
    ] is True
    assert payload["candidate_timeline"][
        "selected_state_variables_without_parent_interpolation"
    ] == ["P"]
    assert (
        payload["candidate_timeline"]["exposed_overlap"]["MU"][
            "target_region_fraction_2d"
        ]
        < 0.5
    )
    assert payload["schedule_source_inventory"]["aggregate_hooks"][
        "post_move_parent_fill_hook_found"
    ] is True
    assert payload["schedule_source_inventory"]["aggregate_hooks"][
        "schedule_post_move_parent_fill_after_parent_interpolation"
    ] is True

    risk_codes = set(payload["summary"]["risk_codes"])
    assert "large_movement_interpolation_timing_ambiguity" in risk_codes
    assert "changed_interpolated_count_mismatch" in risk_codes
    assert "selected_fields_not_all_parent_interpolated" in risk_codes
    assert "target_region_fraction_below_half" in risk_codes
    assert "small_shift_mismatch_from_d57" in risk_codes
    assert "candidate_claims_gate_or_integrator_but_diagnostics_fail" in risk_codes


def test_missing_optional_prior_json_paths_are_reported_without_failure(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    _write_candidate(
        candidate,
        attrs=_base_attrs(
            from_start="10,20",
            to_start="11,21",
            ratio=1,
            changed=100.0,
            interpolated=100.0,
        ),
        mass_shape=(20, 20),
    )

    payload = audit_selected_field_timeline(
        candidate,
        pipeline_json_path=tmp_path / "missing_pipeline.json",
        spatial_alignment_json_path=tmp_path / "missing_spatial.json",
        evolution_json_path=tmp_path / "missing_evolution.json",
        source_root=source_root,
    )

    priors = payload["prior_audit_summaries"]
    assert priors["pipeline"]["status"] == "not_available"
    assert priors["spatial_alignment"]["status"] == "not_available"
    assert priors["evolution"]["status"] == "not_available"
    assert "does not exist" in priors["pipeline"]["message"]
    assert payload["schedule_source_inventory"]["status"] == "not_available"

    no_optional = audit_selected_field_timeline(candidate, source_root=source_root)
    assert no_optional["prior_audit_summaries"]["pipeline"]["status"] == "not_requested"
    assert (
        no_optional["prior_audit_summaries"]["spatial_alignment"]["status"]
        == "not_requested"
    )
    assert no_optional["prior_audit_summaries"]["evolution"]["status"] == "not_requested"


def test_source_inventory_flags_missing_post_move_parent_fill(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    _write_candidate(candidate, attrs=_base_attrs())
    _write_source_snippets(source_root, include_post_fill=False)

    payload = audit_selected_field_timeline(candidate, source_root=source_root)

    assert payload["schedule_source_inventory"]["aggregate_hooks"][
        "post_move_parent_fill_hook_found"
    ] is False
    assert "missing_explicit_post_move_parent_fill_evidence" in set(
        payload["summary"]["risk_codes"]
    )


def test_main_writes_cli_json_output(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "timeline.json"
    source_root = tmp_path / "source"
    _write_candidate(candidate, attrs=_base_attrs(changed=900.0, interpolated=900.0))
    _write_source_snippets(source_root)

    exit_code = audit_main(
        [
            str(candidate),
            "--source-root",
            str(source_root),
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["candidate"] == str(candidate)
    assert payload["candidate_timeline"]["cycle"]["end"] == "2025-07-26_00:10:00"
    assert payload["summary"]["diagnostic_only"] is True


def test_direct_script_help_invocation() -> None:
    project_root = Path(__file__).resolve().parents[2]
    env = {
        **os.environ,
        "UV_CACHE_DIR": ".uv-cache",
        "UV_PYTHON_INSTALL_DIR": ".uv-python",
    }
    result = subprocess.run(
        [
            "uv",
            "run",
            "python",
            "tools/audit_selected_field_timeline.py",
            "--help",
        ],
        cwd=project_root,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "Diagnostic-only audit" in result.stdout
    assert "--pipeline-json" in result.stdout
    assert "--spatial-alignment-json" in result.stdout
