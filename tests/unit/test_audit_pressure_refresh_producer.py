import json
import subprocess
import sys
from pathlib import Path

import netCDF4

from tools.audit_pressure_refresh_producer import (
    audit_pressure_refresh_producer,
    main as audit_main,
    report_to_json,
)


def _write_candidate(path: Path, attrs: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        for key, value in attrs.items():
            dataset.setncattr(key, value)


def _base_attrs() -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "false",
        "TYWRF_GATE_CANDIDATE": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
        "TYWRF_FROM_PARENT_START": "10,20",
        "TYWRF_TO_PARENT_START": "11,20",
        "TYWRF_PARENT_GRID_RATIO": 1,
        "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS": "applied_to_candidate",
        "TYWRF_PRESSURE_REFRESH_HELPER_NAME": "refresh_krosa_moving_nest_pressure",
        "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE": "/tmp/wrfinput_d01",
        "TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT": "true",
        "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE": "moved_candidate_HGT",
        "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE": "override:moved_candidate_HGT",
        "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 2,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT": 2,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT": 6,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS": 6,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS": 6,
        "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS": 2,
        "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS": 2,
        "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS": 4,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT": "true",
        "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO": "true",
        "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED": "true",
    }


def _write_source_tree(root: Path) -> None:
    dynamics = root / "src" / "dynamics"
    include = root / "include" / "tywrf" / "dynamics"
    dynamics.mkdir(parents=True)
    include.mkdir(parents=True)
    (dynamics / "pressure_refresh.cpp").write_text(
        """
namespace tywrf::dynamics {
bool compute_krosa_pressure(const PressureRefreshInputs& inputs, float& perturbation_pressure_pa) {
  const auto mu_total = inputs.mu(0, 0) + inputs.mub(0, 0);
  const auto pfu = inputs.c3f.values[1] * mu_total + inputs.c4f.values[1] + inputs.p_top_pa;
  const auto pfd = inputs.c3f.values[0] * mu_total + inputs.c4f.values[0] + inputs.p_top_pa;
  const auto phm = inputs.c3h.values[0] * mu_total + inputs.c4h.values[0] + inputs.p_top_pa;
  const auto phi = inputs.ph(0, 0, 1) + inputs.phb(0, 0, 1);
  const auto alpha = phi / (phm * log(pfd / pfu)) + inputs.alb(0, 0, 0);
  const auto pb = inputs.pb(0, 0, 0);
  const auto total_pressure = 100000.0 * alpha;
  const auto perturbation = total_pressure - pb;
  perturbation_pressure_pa = perturbation;
  return true;
}
PressureRefreshReport refresh_krosa_moving_nest_pressure(const PressureRefreshInputs& inputs) {
  float refreshed_p = 0.0F;
  compute_krosa_pressure(inputs, refreshed_p);
  inputs.p(0, 0, 0) = refreshed_p;
  return {};
}
}
""",
        encoding="utf-8",
    )
    (dynamics / "pressure_refresh_hook.cpp").write_text(
        """
namespace tywrf::dynamics {
KrosaPressureRefreshHookReport apply_krosa_moving_nest_pressure_refresh_hook() {
  sync_exposed_3d(plan.mass, provider_views.pb, scratch_view.pb);
  sync_exposed_2d(plan.surface, provider_views.mub, scratch_view.mub);
  sync_exposed_3d(plan.w_full, provider_views.phb, scratch_view.phb);
  auto staging = make_krosa_pressure_refresh_inputs(scratch_view, provider_views.alb, metadata);
  auto compute = refresh_krosa_moving_nest_pressure(staging.inputs);
  sync_exposed_3d(plan.mass, const_view(scratch_view.p), state_view.p);
  return {};
}
}
""",
        encoding="utf-8",
    )
    (dynamics / "pressure_refresh_staging.cpp").write_text(
        """
namespace tywrf::dynamics {
PressureRefreshStagingResult make_krosa_pressure_refresh_inputs() {
  PressureRefreshInputs inputs{};
  inputs.p = state.p;
  inputs.pb = const_view(state.pb);
  inputs.alb = external_alb;
  return {};
}
}
""",
        encoding="utf-8",
    )
    for name in (
        "pressure_refresh.hpp",
        "pressure_refresh_hook.hpp",
        "pressure_refresh_staging.hpp",
    ):
        (include / name).write_text("// test header\n", encoding="utf-8")


def _write_vertical_bias(path: Path) -> None:
    payload = {
        "status": "computed",
        "candidate_end": "/tmp/candidate.nc",
        "summary": {
            "diagnostic_only": True,
            "candidate_model_pass": "not_applicable",
        },
        "worst_levels_by_rmse": [0],
        "worst_levels_by_abs_bias": [0],
        "worst_level_attribution": [
            {
                "level": 0,
                "rank_by_rmse": 1,
                "rank_by_abs_bias": 1,
                "pressure_metrics": {
                    "level": 0,
                    "mean_diff": -4167.791,
                    "rmse": 4167.791,
                    "normalized_rmse": 0.2,
                },
                "companion_fields": [
                    {"field": "PB", "mean_diff": 1.0},
                    {"field": "P+PB", "mean_diff": -4167.791},
                    {"field": "MU", "mean_diff": 11.625},
                    {"field": "MUB", "mean_diff": 1.0},
                    {"field": "PHB", "mean_diff": 1.0},
                ],
                "source_evolution": {
                    "status": "computed",
                    "mean_candidate_evolution_minus_reference_evolution_p": -4167.791,
                },
            }
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _flag_codes(payload: dict[str, object]) -> set[str]:
    return {item["code"] for item in payload["formula_risk_flags"]}


def test_metadata_code_and_vertical_risk_detection(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "repo"
    vertical = tmp_path / "vertical_bias.json"
    _write_candidate(candidate, _base_attrs())
    _write_source_tree(source_root)
    _write_vertical_bias(vertical)

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_producer(
                candidate=candidate,
                source_root=source_root,
                vertical_bias_json=vertical,
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert set(payload) >= {
        "metadata_contract",
        "code_producer_inventory",
        "formula_risk_flags",
        "vertical_bias_summary",
    }
    assert payload["metadata_contract"]["helper_source_attrs_present"] is True
    assert payload["metadata_contract"]["counts"]["changed_p_points"] == 6
    assert payload["code_producer_inventory"]["summary"]["p_writer_functions"]
    assert payload["code_producer_inventory"]["summary"]["p_compute_functions"]
    assert payload["vertical_bias_summary"]["levels"][0]["p_mean_diff"] == -4167.791

    codes = _flag_codes(payload)
    assert "pressure_refresh_helper_metadata_present" in codes
    assert "changed_p_points_match_refreshed_count" in codes
    assert "terrain_provenance_uses_moved_candidate_hgt" in codes
    assert "source_inventory_p_written_by_pressure_refresh" in codes
    assert "p_bias_with_weak_base_state_companion_diffs" in codes
    assert "full_pressure_inherits_perturbation_p_bias" in codes
    assert "source_start_evolution_delta_equals_p_bias" in codes


def test_missing_optional_inputs_are_reported_not_available(tmp_path: Path) -> None:
    source_root = tmp_path / "empty_repo"
    source_root.mkdir()
    payload = json.loads(
        report_to_json(audit_pressure_refresh_producer(source_root=source_root))
    )

    assert payload["metadata_contract"]["status"] == "not_available"
    assert payload["vertical_bias_summary"]["status"] == "not_available"
    assert payload["code_producer_inventory"]["status"] == "not_available"
    assert payload["formula_risk_flags"] == []


def test_cli_writes_json_and_help_invocation(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "repo"
    output = tmp_path / "producer.json"
    _write_candidate(candidate, _base_attrs())
    _write_source_tree(source_root)

    exit_code = audit_main(
        [
            "--candidate",
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
    assert payload["metadata_contract"]["status"] == "available"
    assert payload["diagnostic_only"] is True

    result = subprocess.run(
        [sys.executable, "tools/audit_pressure_refresh_producer.py", "--help"],
        cwd=Path(__file__).resolve().parents[2],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "--vertical-bias-json" in result.stdout
