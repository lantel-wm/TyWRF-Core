import json
import os
from pathlib import Path
import subprocess

import netCDF4

from tools.audit_selected_field_pipeline import (
    audit_selected_field_pipeline,
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
    mass_shape: tuple[int, int] = (1000, 1000),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        for key, value in attrs.items():
            dataset.setncattr(key, value)

        mass_ny, mass_nx = mass_shape
        _ensure_dim(dataset, "bottom_top", 2)
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
        dataset.createVariable("MU", "f4", ("Time", "south_north", "west_east"))


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _write_source_snippets(source_root: Path) -> None:
    selected = source_root / "src/tools/selected_field_cycle.cpp"
    selected.parent.mkdir(parents=True, exist_ok=True)
    selected.write_text(
        """
void build_candidate_state() {
  remap_child_state_overlap_only();
  build_exposed_child_state_exchange_plan();
  interpolate_parent_to_exposed_child();
  refresh_static_fields();
  probe_pressure_refresh_provider_readiness();
  probe_pressure_refresh_dry_run_contract();
  apply_pressure_refresh();
  stamp_gate_metadata();
  require_d02_resolution();
  // Oracle inputs are not accepted; selected-field candidate from start states only.
  "TYWRF_GATE_CANDIDATE";
  "TYWRF_INTEGRATOR_OUTPUT";
  "TYWRF_PRESSURE_REFRESH_APPLIED";
}
""",
        encoding="utf-8",
    )
    schedule = source_root / "src/dynamics/cycle_schedule.cpp"
    schedule.parent.mkdir(parents=True, exist_ok=True)
    schedule.write_text(
        """
void schedule() {
  moving_nest_post_move_parent_fill();
  parent_child_interpolation();
  two_way_feedback();
  snap_to_next_parent_step();
}
""",
        encoding="utf-8",
    )
    interpolation = source_root / "src/nest/parent_child_interpolation.cpp"
    interpolation.parent.mkdir(parents=True, exist_ok=True)
    interpolation.write_text(
        """
void interpolate_parent_to_exposed_child() {
  // exposed cells are parent interpolated after overlap remap
}
""",
        encoding="utf-8",
    )


def _base_attrs(
    *,
    from_start: str = "1,1",
    to_start: str = "13,8",
    ratio: int = 5,
) -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "false",
        "TYWRF_GATE_CANDIDATE": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
        "TYWRF_CANDIDATE_DOMAIN": "d02",
        "TYWRF_D02_RESOLUTION_CHECK": "d02_2km",
        "TYWRF_FROM_PARENT_START": from_start,
        "TYWRF_TO_PARENT_START": to_start,
        "TYWRF_PARENT_GRID_RATIO": ratio,
        "TYWRF_STATE_VARIABLES": "U,V,T,PH,MU,P,QVAPOR",
        "TYWRF_PARENT_INTERPOLATED_STATE_VARIABLES": "U,V,T,PH,MU,QVAPOR",
        "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
        "TYWRF_PRESSURE_REFRESH_READINESS_READY": "true",
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS": 123.0,
        "DX": 2000.0,
        "DY": 2000.0,
    }


def test_pipeline_audit_combines_metadata_json_source_and_flags(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    evolution_json = tmp_path / "evolution.json"
    spatial_json = tmp_path / "spatial.json"
    source_root = tmp_path / "source"
    _write_candidate(candidate, attrs=_base_attrs())
    _write_source_snippets(source_root)
    _write_json(
        evolution_json,
        {
            "status": "computed",
            "candidate_model_pass": "not_applicable",
            "summary": {"diagnostic_only": True},
            "variables": [
                {
                    "variable": "U",
                    "status": "computed",
                    "normalized_rmse": 0.08,
                    "rmse": 1.0,
                    "evolution_amplitude": {
                        "amplitude_ratio": 1.02,
                        "capture_fraction": 0.91,
                    },
                    "region_split": {
                        "target_error_fraction": 0.4,
                        "target_region_dominates_global_error": False,
                    },
                }
            ],
        },
    )
    _write_json(
        spatial_json,
        {
            "status": "computed",
            "candidate_model_pass": "not_applicable",
            "summary": {"diagnostic_only": True},
            "variables": [
                {
                    "variable": "U",
                    "status": "computed",
                    "end_state": {
                        "global": {
                            "best_shift": {"di": 1, "dj": 0},
                            "improved": True,
                            "baseline": {"normalized_rmse": 0.08},
                            "best": {"normalized_rmse": 0.07},
                            "normalized_rmse_reduction_fraction": 0.125,
                        }
                    },
                    "evolution": {
                        "global": {
                            "best_shift": {"di": -1, "dj": 1},
                            "improved": True,
                            "baseline": {"normalized_rmse": 0.09},
                            "best": {"normalized_rmse": 0.08},
                            "normalized_rmse_reduction_fraction": 0.111,
                        }
                    },
                }
            ],
        },
    )

    payload = json.loads(
        report_to_json(
            audit_selected_field_pipeline(
                candidate,
                evolution_json_path=evolution_json,
                spatial_alignment_json_path=spatial_json,
                source_root=source_root,
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["candidate_metadata"]["disposition"]["gate_candidate"] is True
    assert (
        payload["candidate_metadata"]["pressure_refresh_attrs"][
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS"
        ]
        == 123.0
    )
    assert payload["movement_geometry"]["movement"]["child_delta"] == {"i": 60, "j": 35}
    assert payload["movement_geometry"]["movement"]["large_movement"] is True
    assert (
        payload["movement_geometry"]["variable_regions"]["MU"][
            "target_region_fraction_2d"
        ]
        < 0.5
    )
    assert (
        payload["evolution_summary"]["variables"]["U"]["amplitude_ratio"]
        == 1.02
    )
    assert payload["spatial_alignment_summary"]["different_best_shifts"] == ["U"]
    assert payload["source_inventory"]["producer_path_assumptions"][
        "exposed_selected_fields_parent_interpolated"
    ] is True
    assert payload["source_inventory"]["producer_path_assumptions"][
        "schedule_contains_post_move_parent_fill_hook"
    ] is True

    risk_codes = set(payload["summary"]["risk_codes"])
    assert "target_region_fraction_lt_0_5" in risk_codes
    assert "target_error_fraction_lt_0_5" in risk_codes
    assert "amplitude_near_one_but_rmse_high" in risk_codes
    assert "best_shift_modest_improvement" in risk_codes
    assert "different_best_shifts_end_vs_evolution" in risk_codes
    assert "candidate_claims_gate_or_integrator_but_audits_fail" in risk_codes
    assert "large_movement_delta_schedule_remap_risk" in risk_codes


def test_missing_optional_json_paths_are_reported_without_failure(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "source"
    _write_candidate(
        candidate,
        attrs=_base_attrs(from_start="10,20", to_start="11,21", ratio=1),
        mass_shape=(10, 10),
    )

    payload = audit_selected_field_pipeline(
        candidate,
        evolution_json_path=tmp_path / "missing_evolution.json",
        spatial_alignment_json_path=tmp_path / "missing_spatial.json",
        source_root=source_root,
    )

    assert payload["evolution_summary"]["status"] == "not_available"
    assert payload["spatial_alignment_summary"]["status"] == "not_available"
    assert "does not exist" in payload["evolution_summary"]["message"]
    assert payload["source_inventory"]["status"] == "not_available"

    no_optional = audit_selected_field_pipeline(candidate, source_root=source_root)
    assert no_optional["evolution_summary"]["status"] == "not_requested"
    assert no_optional["spatial_alignment_summary"]["status"] == "not_requested"


def test_main_writes_cli_json_output(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "pipeline.json"
    source_root = tmp_path / "source"
    _write_candidate(candidate, attrs=_base_attrs(), mass_shape=(20, 20))
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
    assert payload["candidate_metadata"]["disposition"]["parent_grid_ratio"] == 5
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
            "tools/audit_selected_field_pipeline.py",
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
    assert "--evolution-json" in result.stdout
    assert "--spatial-alignment-json" in result.stdout
