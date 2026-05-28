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
