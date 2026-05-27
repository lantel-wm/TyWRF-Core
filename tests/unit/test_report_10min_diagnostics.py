import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.report_10min_diagnostics import (
    main as report_main,
    report_10min_diagnostics,
    report_to_json,
    resolve_report_files,
)


START = "2025-07-26_00:00:00"
END = "2025-07-26_00:10:00"


def _write_wrfout(
    path: Path,
    *,
    i_parent_start: int,
    j_parent_start: int,
    u_value: float = 10.0,
    v_value: float = 5.0,
    dx: float = 2000.0,
    dy: float = 2000.0,
    hgt_offset: float = 0.0,
    attrs: dict[str, object] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    y = np.arange(4, dtype=np.float64).reshape(4, 1)
    x = np.arange(5, dtype=np.float64).reshape(1, 5)
    lat = 20.0 + y * 0.1 + x * 0.01
    lon = 140.0 + y * 0.05 + x * 0.02
    hgt = 10.0 + y * 2.0 + x * 3.0 + hgt_offset

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", 1)
        dataset.createDimension("south_north", lat.shape[0])
        dataset.createDimension("west_east", lat.shape[1])
        dataset.setncattr("I_PARENT_START", i_parent_start)
        dataset.setncattr("J_PARENT_START", j_parent_start)
        dataset.setncattr("DX", dx)
        dataset.setncattr("DY", dy)
        for name, value in (attrs or {}).items():
            dataset.setncattr(name, value)
        fields = {
            "XLAT": lat,
            "XLONG": lon,
            "HGT": hgt,
            "U": np.full(lat.shape, u_value, dtype=np.float64),
            "V": np.full(lat.shape, v_value, dtype=np.float64),
        }
        for name, values in fields.items():
            variable = dataset.createVariable(
                name,
                "f8",
                ("Time", "south_north", "west_east"),
            )
            variable[0, :, :] = values


def test_combined_report_summarizes_moved_segment_and_delta_failure(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=10,
        j_parent_start=20,
        u_value=9.0,
        v_value=5.0,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=12,
        j_parent_start=21,
        u_value=10.0,
        v_value=5.0,
        hgt_offset=1.0,
    )

    report = report_10min_diagnostics(
        reference_dir,
        domain="d02",
        start=START,
        end=END,
        variables=("U", "V"),
        thresholds={"U": 0.05, "V": 0.05},
    )
    payload = json.loads(report_to_json(report, pretty=True))

    assert payload["status"] == "computed"
    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["disposition"]["mode"] == "diagnostic-only"
    assert payload["disposition"]["creates_model_candidate"] is False
    assert "no candidate/model pass" in payload["disposition"]["message"]
    assert payload["summary"]["moved"] is True
    assert payload["summary"]["parent_delta"] == {"i": 2, "j": 1}
    assert payload["summary"]["child_cell_delta"] == {"i": 10, "j": 5}
    assert payload["summary"]["first_failing_variable"] == "U"
    assert payload["summary"]["strict_threshold_exceeded_count"] == 1
    assert payload["summary"]["largest_normalized_deltas"][0]["variable"] == "U"
    assert math.isclose(
        payload["summary"]["largest_normalized_deltas"][0]["normalized_rmse"],
        0.1,
    )
    assert payload["summary"]["d02_2km"] is True
    assert payload["summary"]["d02_2km_status"] == "ok"
    assert payload["movement_audit"]["summary"]["hgt_delta"]["status"] == "computed"
    assert payload["field_delta"]["summary"]["computed"] == 2


def test_report_surfaces_parent_fill_candidate_metadata_when_supplied(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    candidate_file = tmp_path / "candidate" / f"wrfout_d02_{END}_parent_fill"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=10,
        j_parent_start=20,
        u_value=9.0,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=12,
        j_parent_start=21,
        u_value=10.0,
    )
    _write_wrfout(
        candidate_file,
        i_parent_start=12,
        j_parent_start=21,
        attrs={
            "TYWRF_DIAGNOSTIC_REMAP_PARENT_FILL": "true",
            "TYWRF_MINIMUM_STATIC_REFRESH_FIELDS": "XLAT,XLONG,HGT",
            "TYWRF_STAGGERED_STATIC_COORDS_STATUS": "pending_unless_emitted_later",
            "TYWRF_P_DERIVED_REFRESH_STATUS": (
                "pending_derive_or_recompute_after_parent_fill_not_direct_wrf_parent_fill"
            ),
            "TYWRF_PRESSURE_REFRESH_REQUIRED": "true",
            "TYWRF_PRESSURE_REFRESH_OPT_IN": "true",
            "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
            "TYWRF_PRESSURE_REFRESH_REQUIREMENT_STATUS": "required_after_parent_fill",
            "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS": "diagnostic_only_not_integrated",
            "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY": "false",
            "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY": "false",
            "TYWRF_PRESSURE_REFRESH_FORMULA_STATUS": "staged_for_later_kernel",
            "TYWRF_PRESSURE_REFRESH_FORMULA_STAGING_NAME": "wrf_pressure_base_plus_perturbation",
            "TYWRF_PRESSURE_REFRESH_REGION_STAGING_NAME": "parent_fill_child_domain",
            "TYWRF_PRESSURE_REFRESH_THERMODYNAMIC_MODE": "dry_reference",
            "TYWRF_PRESSURE_REFRESH_REQUIRED_INPUTS": "PB,P,MU,MUB,PH,PHB,T",
            "TYWRF_PRESSURE_REFRESH_OUTPUT_FIELD": "P",
            "TYWRF_PRESSURE_REFRESH_HELPER_NAME": "refresh_pressure_after_parent_fill",
            "TYWRF_PRESSURE_REFRESH_ALB_SOURCE": "test_alb_provider",
            "TYWRF_PRESSURE_REFRESH_PROVIDER_OK": "true",
            "TYWRF_PRESSURE_REFRESH_STAGING_OK": "true",
            "TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED": "true",
            "TYWRF_PRESSURE_REFRESH_TERRAIN_OVERRIDE_USED": "true",
            "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE": "moved_candidate_HGT",
            "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE": "override:moved_candidate_HGT",
            "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS": 42,
            "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS": "43",
            "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS": 44,
            "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 9,
            "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT": "9",
            "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT": 90,
            "TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT": "0",
            "TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT": 0,
            "TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS": "false",
            "TYWRF_PRESSURE_REFRESH_TOUCHED_HALO_CELLS": "false",
            "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS": "45",
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS": "90",
            "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS": 91,
            "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS": "92",
            "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS": 93,
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT": "true",
            "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO": "true",
            "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED": "true",
            "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE": "candidate_attrs",
            "TYWRF_PRESSURE_REFRESH_METADATA_TIME_INDEX": 0,
            "TYWRF_DIRECT_WRF_END_STATE_ORACLE_STATUS": (
                "diagnostic_only_nonphysical_non_gate"
            ),
            "TYWRF_DIAGNOSTIC_ONLY": "true",
            "TYWRF_GATE_CANDIDATE": "false",
            "TYWRF_INTEGRATOR_OUTPUT": "false",
            "TYWRF_VALIDATION_GATE_ONLY": "false",
            "TYWRF_CANDIDATE_KIND": "pressure_refresh_remap_diagnostic",
        },
    )

    payload = json.loads(
        report_to_json(
            report_10min_diagnostics(
                reference_dir,
                domain="d02",
                start=START,
                end=END,
                variables=("U",),
                thresholds={"U": 0.05},
                candidate_file=candidate_file,
            ),
            pretty=True,
        )
    )

    assert payload["status"] == "computed"
    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["candidate_file"] == str(candidate_file)
    assert payload["candidate_metadata"]["status"] == "available"
    assert payload["candidate_metadata"]["path"] == str(candidate_file)
    assert payload["candidate_metadata"]["gate_candidate"] is False
    assert payload["candidate_metadata"]["candidate_diagnostic_only"] is True
    assert payload["candidate_metadata"]["integrator_output"] is False
    assert payload["candidate_metadata"]["validation_gate_only"] is False
    assert (
        payload["candidate_metadata"]["candidate_kind"]
        == "pressure_refresh_remap_diagnostic"
    )
    assert payload["candidate_metadata"]["candidate_gate_eligible"] is False
    assert payload["candidate_metadata"]["candidate_gate_blockers"] == [
        "TYWRF_DIAGNOSTIC_ONLY=true",
        "TYWRF_GATE_CANDIDATE=false",
        "TYWRF_INTEGRATOR_OUTPUT=false",
        "TYWRF_CANDIDATE_KIND=pressure_refresh_remap_diagnostic",
    ]
    assert payload["candidate_metadata"]["candidate_model_pass"] == "not_applicable"
    assert (
        payload["candidate_metadata"]["attrs"]["TYWRF_DIAGNOSTIC_REMAP_PARENT_FILL"]
        == "true"
    )
    assert (
        payload["candidate_metadata"]["attrs"]["TYWRF_PRESSURE_REFRESH_REQUIRED"]
        == "true"
    )
    assert (
        payload["candidate_metadata"]["attrs"][
            "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS"
        ]
        == "diagnostic_only_not_integrated"
    )
    assert (
        payload["candidate_metadata"]["attrs"]["TYWRF_PRESSURE_REFRESH_OUTPUT_FIELD"]
        == "P"
    )
    assert (
        payload["candidate_metadata"]["attrs"]["TYWRF_PRESSURE_REFRESH_ALB_SOURCE"]
        == "test_alb_provider"
    )
    assert (
        payload["candidate_metadata"]["attrs"][
            "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS"
        ]
        == 42
    )
    assert (
        payload["candidate_metadata"]["attrs"][
            "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT"
        ]
        == 9
    )
    assert (
        payload["candidate_metadata"]["attrs"][
            "TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS"
        ]
        == "false"
    )
    assert (
        payload["candidate_metadata"]["attrs"][
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT"
        ]
        == "true"
    )
    assert payload["parent_fill_metadata"]["status"] == "available"
    assert payload["parent_fill_metadata"]["pressure_refresh_metadata_status"] == "available"
    assert payload["parent_fill_metadata"]["diagnostic_remap_parent_fill"] is True
    assert payload["parent_fill_metadata"]["minimum_static_refresh_fields"] == [
        "XLAT",
        "XLONG",
        "HGT",
    ]
    assert (
        payload["parent_fill_metadata"]["p_derived_refresh_status"]
        == "pending_derive_or_recompute_after_parent_fill_not_direct_wrf_parent_fill"
    )
    assert (
        payload["parent_fill_metadata"]["direct_wrf_end_state_oracle_status"]
        == "diagnostic_only_nonphysical_non_gate"
    )
    assert payload["parent_fill_metadata"]["pressure_refresh_required"] is True
    assert payload["parent_fill_metadata"]["pressure_refresh_opt_in"] is True
    assert payload["parent_fill_metadata"]["pressure_refresh_applied"] is True
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_experimental_apply"]
        is False
    )
    assert (
        payload["parent_fill_metadata"]["experimental_pressure_refresh_apply"]
        is False
    )
    assert (
        payload["parent_fill_metadata"]["normal_candidate_pressure_refresh"]
        is False
    )
    assert payload["parent_fill_metadata"]["hidden_seam_pressure_refresh"] is False
    assert payload["parent_fill_metadata"]["candidate_gate_eligible"] is False
    assert payload["parent_fill_metadata"]["candidate_gate_blockers"] == [
        "TYWRF_DIAGNOSTIC_ONLY=true",
        "TYWRF_GATE_CANDIDATE=false",
        "TYWRF_INTEGRATOR_OUTPUT=false",
        "TYWRF_CANDIDATE_KIND=pressure_refresh_remap_diagnostic",
    ]
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_disposition"]["status"]
        == "diagnostic_helper_evidence_only"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_disposition"][
            "candidate_model_pass"
        ]
        == "not_applicable"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_requirement_status"]
        == "required_after_parent_fill"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_integration_status"]
        == "diagnostic_only_not_integrated"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_formula_status"]
        == "staged_for_later_kernel"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_formula_staging_name"]
        == "wrf_pressure_base_plus_perturbation"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_region_staging_name"]
        == "parent_fill_child_domain"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_thermodynamic_mode"]
        == "dry_reference"
    )
    assert payload["parent_fill_metadata"]["pressure_refresh_required_inputs"] == [
        "PB",
        "P",
        "MU",
        "MUB",
        "PH",
        "PHB",
        "T",
    ]
    assert payload["parent_fill_metadata"]["pressure_refresh_output_field"] == "P"
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_helper_name"]
        == "refresh_pressure_after_parent_fill"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_alb_source"]
        == "test_alb_provider"
    )
    assert payload["parent_fill_metadata"]["pressure_refresh_provider_ok"] is True
    assert payload["parent_fill_metadata"]["pressure_refresh_staging_ok"] is True
    assert payload["parent_fill_metadata"]["pressure_refresh_compute_called"] is True
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_terrain_override_used"]
        is True
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_terrain_source"]
        == "moved_candidate_HGT"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_terrain_provenance"]
        == "override:moved_candidate_HGT"
    )
    assert payload["parent_fill_metadata"]["pressure_refresh_synced_pb_points"] == 42
    assert payload["parent_fill_metadata"]["pressure_refresh_synced_mub_points"] == 43
    assert payload["parent_fill_metadata"]["pressure_refresh_synced_phb_points"] == 44
    assert payload["parent_fill_metadata"]["pressure_refresh_target_column_count"] == 9
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_refreshed_column_count"]
        == 9
    )
    assert payload["parent_fill_metadata"]["pressure_refresh_refreshed_point_count"] == 90
    assert payload["parent_fill_metadata"]["pressure_refresh_skipped_point_count"] == 0
    assert payload["parent_fill_metadata"]["pressure_refresh_invalid_point_count"] == 0
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_touched_overlap_cells"]
        is False
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_touched_halo_cells"]
        is False
    )
    assert payload["parent_fill_metadata"]["pressure_refresh_refreshed_p_points"] == 45
    assert payload["parent_fill_metadata"]["pressure_refresh_changed_p_points"] == 90
    assert payload["parent_fill_metadata"]["pressure_refresh_changed_pb_points"] == 91
    assert payload["parent_fill_metadata"]["pressure_refresh_changed_mub_points"] == 92
    assert payload["parent_fill_metadata"]["pressure_refresh_changed_phb_points"] == 93
    assert (
        payload["parent_fill_metadata"][
            "pressure_refresh_changed_p_matches_refreshed_point_count"
        ]
        is True
    )
    assert (
        payload["parent_fill_metadata"][
            "pressure_refresh_invalid_and_skipped_points_zero"
        ]
        is True
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_overlap_halo_untouched"]
        is True
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_metadata_source"]
        == "candidate_attrs"
    )
    assert payload["parent_fill_metadata"]["pressure_refresh_metadata_time_index"] == 0
    assert isinstance(
        payload["parent_fill_metadata"]["pressure_refresh_provider_ok"], bool
    )
    assert isinstance(
        payload["parent_fill_metadata"]["pressure_refresh_synced_pb_points"], int
    )
    assert isinstance(
        payload["parent_fill_metadata"]["pressure_refresh_metadata_source"], str
    )
    assert payload["summary"]["candidate_file"] == str(candidate_file)
    assert payload["summary"]["candidate_metadata_status"] == "available"
    assert payload["summary"]["candidate_gate_candidate"] is False
    assert payload["summary"]["candidate_gate_eligible"] is False
    assert payload["summary"]["candidate_gate_blockers"] == [
        "TYWRF_DIAGNOSTIC_ONLY=true",
        "TYWRF_GATE_CANDIDATE=false",
        "TYWRF_INTEGRATOR_OUTPUT=false",
        "TYWRF_CANDIDATE_KIND=pressure_refresh_remap_diagnostic",
    ]
    assert payload["summary"]["candidate_diagnostic_only"] is True
    assert payload["summary"]["candidate_integrator_output"] is False
    assert payload["summary"]["candidate_validation_gate_only"] is False
    assert payload["summary"]["candidate_kind"] == "pressure_refresh_remap_diagnostic"
    assert payload["summary"]["parent_fill_metadata_status"] == "available"
    assert payload["summary"]["pressure_refresh_metadata_status"] == "available"
    assert payload["summary"]["diagnostic_remap_parent_fill"] is True
    assert payload["summary"]["minimum_static_refresh_fields"] == [
        "XLAT",
        "XLONG",
        "HGT",
    ]
    assert payload["summary"]["pressure_refresh_required"] is True
    assert payload["summary"]["pressure_refresh_opt_in"] is True
    assert payload["summary"]["pressure_refresh_applied"] is True
    assert (
        payload["summary"]["pressure_refresh_integration_status"]
        == "diagnostic_only_not_integrated"
    )
    assert payload["summary"]["pressure_refresh_experimental_apply"] is False
    assert payload["summary"]["experimental_pressure_refresh_apply"] is False
    assert payload["summary"]["normal_candidate_pressure_refresh"] is False
    assert payload["summary"]["hidden_seam_pressure_refresh"] is False
    assert payload["summary"]["pressure_refresh_required_inputs"] == [
        "PB",
        "P",
        "MU",
        "MUB",
        "PH",
        "PHB",
        "T",
    ]
    assert payload["summary"]["pressure_refresh_output_field"] == "P"
    assert payload["summary"]["pressure_refresh_alb_source"] == "test_alb_provider"
    assert payload["summary"]["pressure_refresh_provider_ok"] is True
    assert payload["summary"]["pressure_refresh_staging_ok"] is True
    assert payload["summary"]["pressure_refresh_compute_called"] is True
    assert payload["summary"]["pressure_refresh_terrain_override_used"] is True
    assert payload["summary"]["pressure_refresh_terrain_source"] == "moved_candidate_HGT"
    assert (
        payload["summary"]["pressure_refresh_terrain_provenance"]
        == "override:moved_candidate_HGT"
    )
    assert payload["summary"]["pressure_refresh_synced_pb_points"] == 42
    assert payload["summary"]["pressure_refresh_synced_mub_points"] == 43
    assert payload["summary"]["pressure_refresh_synced_phb_points"] == 44
    assert payload["summary"]["pressure_refresh_target_column_count"] == 9
    assert payload["summary"]["pressure_refresh_refreshed_column_count"] == 9
    assert payload["summary"]["pressure_refresh_refreshed_point_count"] == 90
    assert payload["summary"]["pressure_refresh_skipped_point_count"] == 0
    assert payload["summary"]["pressure_refresh_invalid_point_count"] == 0
    assert payload["summary"]["pressure_refresh_touched_overlap_cells"] is False
    assert payload["summary"]["pressure_refresh_touched_halo_cells"] is False
    assert payload["summary"]["pressure_refresh_refreshed_p_points"] == 45
    assert payload["summary"]["pressure_refresh_changed_p_points"] == 90
    assert payload["summary"]["pressure_refresh_changed_pb_points"] == 91
    assert payload["summary"]["pressure_refresh_changed_mub_points"] == 92
    assert payload["summary"]["pressure_refresh_changed_phb_points"] == 93
    assert (
        payload["summary"][
            "pressure_refresh_changed_p_matches_refreshed_point_count"
        ]
        is True
    )
    assert (
        payload["summary"]["pressure_refresh_invalid_and_skipped_points_zero"]
        is True
    )
    assert payload["summary"]["pressure_refresh_overlap_halo_untouched"] is True
    assert payload["summary"]["pressure_refresh_metadata_source"] == "candidate_attrs"
    assert payload["summary"]["pressure_refresh_metadata_time_index"] == 0
    assert (
        payload["summary"]["pressure_refresh_disposition"]["status"]
        == "diagnostic_helper_evidence_only"
    )
    assert (
        payload["summary"]["pressure_refresh_disposition"]["pressure_refresh_applied"]
        is True
    )
    assert (
        payload["summary"]["pressure_refresh_disposition"]["candidate_gate_eligible"]
        is False
    )
    assert (
        payload["summary"]["pressure_refresh_disposition"]["candidate_model_pass"]
        == "not_applicable"
    )
    assert isinstance(payload["summary"]["pressure_refresh_compute_called"], bool)
    assert isinstance(payload["summary"]["pressure_refresh_refreshed_p_points"], int)
    assert isinstance(payload["summary"]["pressure_refresh_alb_source"], str)
    assert payload["summary"]["diagnostic_only"] is True
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"
    assert payload["disposition"]["creates_model_candidate"] is False
    assert payload["disposition"]["candidate_gate_eligible"] is False
    assert (
        payload["disposition"]["pressure_refresh_disposition"]["status"]
        == "diagnostic_helper_evidence_only"
    )


def test_report_surfaces_candidate_gate_eligible_metadata(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    candidate_file = tmp_path / "candidate" / f"wrfout_d02_{END}_candidate"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=10,
        j_parent_start=20,
        u_value=9.0,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=12,
        j_parent_start=21,
        u_value=10.0,
    )
    _write_wrfout(
        candidate_file,
        i_parent_start=12,
        j_parent_start=21,
        attrs={
            "TYWRF_DIAGNOSTIC_ONLY": "false",
            "TYWRF_GATE_CANDIDATE": "true",
            "TYWRF_INTEGRATOR_OUTPUT": "true",
            "TYWRF_VALIDATION_GATE_ONLY": "false",
            "TYWRF_CANDIDATE_KIND": "integrator_candidate",
            "TYWRF_PRESSURE_REFRESH_APPLIED": "false",
        },
    )

    payload = json.loads(
        report_to_json(
            report_10min_diagnostics(
                reference_dir,
                domain="d02",
                start=START,
                end=END,
                variables=("U",),
                thresholds={"U": 0.05},
                candidate_file=candidate_file,
            )
        )
    )

    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["candidate_metadata"]["candidate_gate_eligible"] is True
    assert payload["candidate_metadata"]["candidate_gate_blockers"] == []
    assert payload["candidate_metadata"]["candidate_diagnostic_only"] is False
    assert payload["candidate_metadata"]["gate_candidate"] is True
    assert payload["candidate_metadata"]["integrator_output"] is True
    assert payload["candidate_metadata"]["validation_gate_only"] is False
    assert payload["candidate_metadata"]["candidate_kind"] == "integrator_candidate"
    assert payload["summary"]["candidate_gate_eligible"] is True
    assert payload["summary"]["candidate_gate_blockers"] == []
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"
    assert payload["disposition"]["candidate_gate_eligible"] is True
    assert (
        payload["summary"]["pressure_refresh_disposition"]["status"] == "not_applied"
    )


def test_report_distinguishes_normal_pressure_refresh_candidate_metadata(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    candidate_file = tmp_path / "candidate" / f"wrfout_d02_{END}_candidate"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=10,
        j_parent_start=20,
        u_value=9.0,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=12,
        j_parent_start=21,
        u_value=10.0,
    )
    _write_wrfout(
        candidate_file,
        i_parent_start=12,
        j_parent_start=21,
        attrs={
            "TYWRF_DIAGNOSTIC_ONLY": "false",
            "TYWRF_GATE_CANDIDATE": "true",
            "TYWRF_INTEGRATOR_OUTPUT": "true",
            "TYWRF_VALIDATION_GATE_ONLY": "false",
            "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
            "TYWRF_PRESSURE_REFRESH_OPT_IN": "true",
            "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
            "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS": "applied_to_candidate",
            "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY": "false",
            "TYWRF_PRESSURE_REFRESH_TERRAIN_OVERRIDE_USED": "true",
            "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE": "moved_candidate_HGT",
            "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE": "override:moved_candidate_HGT",
            "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 11,
            "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT": 11,
            "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT": 110,
            "TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT": 0,
            "TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT": 0,
            "TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS": "false",
            "TYWRF_PRESSURE_REFRESH_TOUCHED_HALO_CELLS": "false",
            "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS": 110,
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS": 110,
            "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS": 111,
            "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS": 112,
            "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS": 113,
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT": "true",
            "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO": "true",
            "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED": "true",
        },
    )

    payload = json.loads(
        report_to_json(
            report_10min_diagnostics(
                reference_dir,
                domain="d02",
                start=START,
                end=END,
                variables=("U",),
                thresholds={"U": 0.05},
                candidate_file=candidate_file,
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["candidate_metadata"]["candidate_gate_eligible"] is True
    assert payload["candidate_metadata"]["candidate_gate_blockers"] == []
    assert payload["parent_fill_metadata"]["pressure_refresh_metadata_status"] == "available"
    assert payload["parent_fill_metadata"]["normal_candidate_pressure_refresh"] is True
    assert payload["parent_fill_metadata"]["hidden_seam_pressure_refresh"] is False
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_disposition"]["status"]
        == "normal_candidate_metadata_only"
    )
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_disposition"][
            "candidate_model_pass"
        ]
        == "not_applicable"
    )
    assert payload["summary"]["normal_candidate_pressure_refresh"] is True
    assert payload["summary"]["pressure_refresh_metadata_status"] == "available"
    assert payload["summary"]["hidden_seam_pressure_refresh"] is False
    assert payload["summary"]["pressure_refresh_target_column_count"] == 11
    assert payload["summary"]["pressure_refresh_refreshed_column_count"] == 11
    assert payload["summary"]["pressure_refresh_refreshed_point_count"] == 110
    assert payload["summary"]["pressure_refresh_changed_p_points"] == 110
    assert (
        payload["summary"]["pressure_refresh_changed_p_matches_refreshed_point_count"]
        is True
    )
    assert payload["summary"]["pressure_refresh_overlap_halo_untouched"] is True
    assert (
        payload["disposition"]["pressure_refresh_disposition"]["status"]
        == "normal_candidate_metadata_only"
    )


def test_report_rejects_hidden_seam_pressure_refresh_metadata(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    candidate_file = tmp_path / "candidate" / f"wrfout_d02_{END}_hidden_seam"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=10,
        j_parent_start=20,
        u_value=9.0,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=12,
        j_parent_start=21,
        u_value=10.0,
    )
    _write_wrfout(
        candidate_file,
        i_parent_start=12,
        j_parent_start=21,
        attrs={
            "TYWRF_DIAGNOSTIC_ONLY": "true",
            "TYWRF_GATE_CANDIDATE": "false",
            "TYWRF_INTEGRATOR_OUTPUT": "false",
            "TYWRF_VALIDATION_GATE_ONLY": "false",
            "TYWRF_CANDIDATE_KIND": "selected_field_pressure_refresh_experimental_apply_v0",
            "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY": "true",
            "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY": "true",
            "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
            "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS": "experimental_apply_test_only",
            "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 4,
            "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT": 4,
            "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT": 40,
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS": 40,
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT": "true",
            "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO": "true",
            "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED": "true",
        },
    )

    payload = json.loads(
        report_to_json(
            report_10min_diagnostics(
                reference_dir,
                domain="d02",
                start=START,
                end=END,
                variables=("U",),
                thresholds={"U": 0.05},
                candidate_file=candidate_file,
            )
        )
    )

    blockers = payload["candidate_metadata"]["candidate_gate_blockers"]
    assert payload["candidate_metadata"]["candidate_gate_eligible"] is False
    assert "TYWRF_DIAGNOSTIC_ONLY=true" in blockers
    assert "TYWRF_GATE_CANDIDATE=false" in blockers
    assert "TYWRF_INTEGRATOR_OUTPUT=false" in blockers
    assert "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY=true" in blockers
    assert "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY=true" in blockers
    assert (
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS=experimental_apply_test_only"
        in blockers
    )
    assert payload["parent_fill_metadata"]["hidden_seam_pressure_refresh"] is True
    assert payload["parent_fill_metadata"]["pressure_refresh_metadata_status"] == "available"
    assert payload["parent_fill_metadata"]["normal_candidate_pressure_refresh"] is False
    assert (
        payload["parent_fill_metadata"]["pressure_refresh_disposition"]["status"]
        == "hidden_seam_diagnostic_evidence_only"
    )
    assert payload["summary"]["hidden_seam_pressure_refresh"] is True
    assert payload["summary"]["candidate_gate_eligible"] is False
    assert payload["summary"]["pressure_refresh_refreshed_point_count"] == 40
    assert (
        payload["summary"]["pressure_refresh_disposition"]["candidate_model_pass"]
        == "not_applicable"
    )
    assert payload["candidate_model_pass"] == "not_applicable"


def test_report_marks_missing_positive_candidate_metadata_ineligible(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    candidate_file = tmp_path / "candidate" / f"wrfout_d02_{END}_candidate"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=10,
        j_parent_start=20,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=10,
        j_parent_start=20,
    )
    _write_wrfout(candidate_file, i_parent_start=10, j_parent_start=20)

    payload = json.loads(
        report_to_json(
            report_10min_diagnostics(
                reference_dir,
                domain="d02",
                start=START,
                end=END,
                variables=("U",),
                thresholds={"U": 0.05},
                candidate_file=candidate_file,
            )
        )
    )

    assert payload["candidate_metadata"]["candidate_gate_eligible"] is False
    assert payload["candidate_metadata"]["candidate_gate_blockers"] == [
        "TYWRF_GATE_CANDIDATE is not true",
        "TYWRF_INTEGRATOR_OUTPUT is not true",
    ]
    assert payload["summary"]["candidate_gate_eligible"] is False
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"


def test_resolve_report_files_uses_reference_dir_naming(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"

    files = resolve_report_files(
        reference_dir,
        domain="d02",
        start=START,
        end=END,
    )

    assert files.start_file == reference_dir / f"wrfout_d02_{START}"
    assert files.end_file == reference_dir / f"wrfout_d02_{END}"
    assert files.domain == "d02"
    assert files.start_time == START
    assert files.end_time == END


def test_unmoved_report_remains_diagnostic_only_without_candidate_pass(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=7,
        j_parent_start=8,
        u_value=4.0,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=7,
        j_parent_start=8,
        u_value=4.0,
    )

    payload = json.loads(
        report_to_json(
            report_10min_diagnostics(
                reference_dir,
                domain="d02",
                start=START,
                end=END,
                variables=("U",),
                thresholds={"U": 0.05},
            )
        )
    )

    assert payload["summary"]["moved"] is False
    assert payload["summary"]["parent_delta"] == {"i": 0, "j": 0}
    assert payload["summary"]["child_cell_delta"] == {"i": 0, "j": 0}
    assert payload["summary"]["strict_threshold_exceeded_count"] == 0
    assert payload["summary"]["first_failing_variable"] is None
    assert payload["summary"]["diagnostic_only"] is True
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"
    assert "diagnostic-only" in payload["summary"]["message"]


def test_cli_writes_json_and_passes_log_file_to_movement_audit(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    output = tmp_path / "report.json"
    log_file = tmp_path / "rsl.out.0000"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=1,
        j_parent_start=2,
        u_value=3.0,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=3,
        j_parent_start=3,
        u_value=3.0,
    )
    log_file.write_text(
        "\n".join(
            [
                " 2025-07-26_00:00:40 move (rel cd) :            1           1",
                "  moving            2           1           1",
                " 2025-07-26_00:01:20 move (rel cd) :            1           0",
                "  moving            2           1           0",
            ]
        ),
        encoding="utf-8",
    )

    exit_code = report_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--domain",
            "d02",
            "--start",
            START,
            "--end",
            END,
            "--variables",
            "U",
            "--log-file",
            str(log_file),
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["start_file"] == str(reference_dir / f"wrfout_d02_{START}")
    assert payload["end_file"] == str(reference_dir / f"wrfout_d02_{END}")
    assert payload["summary"]["moved"] is True
    assert payload["movement_audit"]["log_events"]["status"] == "available"
    assert payload["movement_audit"]["log_events"]["net_applied_parent_delta"] == {
        "i": 2,
        "j": 1,
    }
    assert (
        payload["movement_audit"]["log_events"]["net_applied_matches_parent_delta"]
        is True
    )


def test_cli_accepts_candidate_file_metadata_without_model_pass(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    output = tmp_path / "report.json"
    candidate_file = tmp_path / "candidate" / f"wrfout_d02_{END}_parent_fill"
    _write_wrfout(
        reference_dir / f"wrfout_d02_{START}",
        i_parent_start=1,
        j_parent_start=2,
        u_value=3.0,
    )
    _write_wrfout(
        reference_dir / f"wrfout_d02_{END}",
        i_parent_start=3,
        j_parent_start=3,
        u_value=3.0,
    )
    _write_wrfout(
        candidate_file,
        i_parent_start=3,
        j_parent_start=3,
        attrs={
            "TYWRF_DIAGNOSTIC_REMAP_PARENT_FILL": "true",
            "TYWRF_MINIMUM_STATIC_REFRESH_FIELDS": "XLAT,XLONG,HGT",
            "TYWRF_STAGGERED_STATIC_COORDS_STATUS": "pending_unless_emitted_later",
            "TYWRF_P_DERIVED_REFRESH_STATUS": "pending_refresh",
            "TYWRF_DIRECT_WRF_END_STATE_ORACLE_STATUS": (
                "diagnostic_only_nonphysical_non_gate"
            ),
            "TYWRF_GATE_CANDIDATE": "false",
        },
    )

    exit_code = report_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--domain",
            "d02",
            "--candidate-file",
            str(candidate_file),
            "--start",
            START,
            "--end",
            END,
            "--variables",
            "U",
            "--output",
            str(output),
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["candidate_file"] == str(candidate_file)
    assert payload["candidate_metadata"]["status"] == "available"
    assert payload["parent_fill_metadata"]["diagnostic_remap_parent_fill"] is True
    assert payload["summary"]["candidate_gate_candidate"] is False
    assert payload["summary"]["diagnostic_only"] is True
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"
