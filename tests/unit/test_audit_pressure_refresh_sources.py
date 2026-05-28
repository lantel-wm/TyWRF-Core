import json
import math
import os
from pathlib import Path
import subprocess

import netCDF4
import numpy as np

from tools.audit_pressure_refresh_sources import (
    PressureDiagnosticRequest,
    SourceEntry,
    audit_pressure_refresh_sources,
    main as audit_main,
    report_to_json,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _write_vector(dataset: netCDF4.Dataset, name: str, length: int) -> None:
    variable = dataset.createVariable(name, "f4", ("Time", f"{name}_dim"))
    variable[0, :] = list(range(length))


def _write_source(
    path: Path,
    *,
    include_alb: bool = True,
    include_phb: bool = True,
    include_terrain: bool = True,
    include_direct_pb_t_init: bool = False,
    bad_alb_shape: bool = False,
    bad_phb_shape: bool = False,
    p_top_as_attr: bool = False,
) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        dataset.createDimension("bottom_top", 2)
        dataset.createDimension("bottom_top_stag", 3)
        dataset.createDimension("south_north", 3)
        dataset.createDimension("west_east", 4)
        dataset.createDimension("west_east_bad", 5)
        dataset.createDimension("C3F_dim", 3)
        dataset.createDimension("C4F_dim", 3)
        dataset.createDimension("C3H_dim", 2)
        dataset.createDimension("C4H_dim", 2)

        if p_top_as_attr:
            dataset.setncattr("P_TOP", 5000.0)
        else:
            p_top = dataset.createVariable("P_TOP", "f4", ("Time",))
            p_top[0] = 5000.0

        _write_vector(dataset, "C3F", 3)
        _write_vector(dataset, "C4F", 3)
        _write_vector(dataset, "C3H", 2)
        _write_vector(dataset, "C4H", 2)

        if include_terrain:
            hgt = dataset.createVariable("HGT", "f4", ("Time", "south_north", "west_east"))
            hgt[0, :, :] = 100.0

        if include_direct_pb_t_init:
            pb = dataset.createVariable(
                "PB",
                "f4",
                ("Time", "bottom_top", "south_north", "west_east"),
            )
            pb[0, :, :, :] = 95000.0
            t_init = dataset.createVariable(
                "T_INIT",
                "f4",
                ("Time", "bottom_top", "south_north", "west_east"),
            )
            t_init[0, :, :, :] = 0.0

        if include_alb:
            west_east_dim = "west_east_bad" if bad_alb_shape else "west_east"
            alb = dataset.createVariable(
                "ALB",
                "f4",
                ("Time", "bottom_top", "south_north", west_east_dim),
            )
            alb[0, :, :, :] = 1.0

        if include_phb:
            west_east_dim = "west_east_bad" if bad_phb_shape else "west_east"
            phb = dataset.createVariable(
                "PHB",
                "f4",
                ("Time", "bottom_top_stag", "south_north", west_east_dim),
            )
            phb[0, :, :, :] = 100.0


def _ensure_dim(dataset: netCDF4.Dataset, name: str, size: int) -> None:
    if name not in dataset.dimensions:
        dataset.createDimension(name, size)


def _write_wrfout(
    path: Path,
    variables: dict[str, np.ndarray],
    *,
    attrs: dict[str, object] | None = None,
) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        for key, value in (attrs or {}).items():
            dataset.setncattr(key, value)
        for name, raw_values in variables.items():
            values = np.asarray(raw_values, dtype=np.float64)
            if values.ndim == 3:
                dims = ("Time", "bottom_top", "south_north", "west_east")
                _ensure_dim(dataset, "bottom_top", values.shape[0])
                _ensure_dim(dataset, "south_north", values.shape[1])
                _ensure_dim(dataset, "west_east", values.shape[2])
            elif values.ndim == 2:
                dims = ("Time", "south_north", "west_east")
                _ensure_dim(dataset, "south_north", values.shape[0])
                _ensure_dim(dataset, "west_east", values.shape[1])
            else:
                dims = tuple(f"{name}_dim_{axis}" for axis in range(values.ndim))
                for dim_name, size in zip(dims, values.shape, strict=True):
                    _ensure_dim(dataset, dim_name, size)
            variable = dataset.createVariable(name, "f8", dims)
            if dims and dims[0] == "Time":
                variable[0, ...] = values
            else:
                variable[:] = values


def _pressure_attrs() -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "false",
        "TYWRF_GATE_CANDIDATE": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "pressure_refresh_candidate",
        "TYWRF_FROM_PARENT_START": "10,20",
        "TYWRF_TO_PARENT_START": "11,20",
        "TYWRF_PARENT_GRID_RATIO": 1,
        "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS": "applied_to_candidate",
        "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 2,
    }


def _entry(name: str, domain: str, path: Path, **overrides: object) -> SourceEntry:
    payload = {
        "name": name,
        "domain": domain,
        "path": str(path),
        "mass_nx": 4,
        "mass_ny": 3,
        "mass_nz": 2,
        "full_nz": 3,
        "time_index": 0,
    }
    payload.update(overrides)
    return SourceEntry(**payload)


def test_audit_reports_complete_source_as_seedable_and_clean_only_when_marked(
    tmp_path: Path,
) -> None:
    source = tmp_path / "wrfinput_d01"
    _write_source(source)

    report = audit_pressure_refresh_sources(
        [_entry("d01_wrfinput", "d01", source, suitable_for_start_time_truth=True)]
    )
    payload = json.loads(report_to_json(report))

    assert payload["status"] == "ok"
    assert payload["summary"]["ok_count"] == 1
    entry = payload["entries"][0]
    assert entry["status"] == "ok"
    assert entry["missing_names"] == []
    assert entry["p_top_present"] is True
    assert entry["p_top_source"] == "time_variable"
    assert entry["alb_available"] is True
    assert entry["direct_alb_source"] is True
    assert entry["direct_phb_present"] is True
    assert entry["direct_phb_shape_valid"] is True
    assert entry["base_state_reconstruction_inputs_complete"] is True
    assert entry["missing_base_state_reconstruction_inputs"] == []
    assert entry["phb_reconstruction_inputs_complete"] is True
    assert entry["missing_phb_reconstruction_inputs"] == []
    assert entry["base_state_terrain_source"] == "HGT"
    assert entry["base_state_reconstruction_input_presence"] == {
        "HT": False,
        "HGT": True,
        "HT/HGT": True,
        "P_TOP": True,
        "C3F": True,
        "C4F": True,
        "C3H": True,
        "C4H": True,
        "PB": False,
        "T_INIT": False,
        "direct_PB_T_INIT": False,
    }
    assert entry["direct_pb_t_init_path_available"] is False
    assert entry["missing_direct_pb_t_init_inputs"] == ["PB", "T_INIT"]
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
    assert (
        entry["full_provider_recommended_source"]
        == "wrf_start_domain_phb_hypsometric_opt2_reconstruction"
    )
    assert entry["diagnostic_only"] is True
    assert entry["can_seed_pressure_refresh"] is True
    assert entry["suitable_for_start_time_truth"] is True


def test_audit_reports_d02_start_time_missing_alb_as_blocker(tmp_path: Path) -> None:
    source = tmp_path / "wrfout_d02_2025-07-26_00:00:00"
    _write_source(source, include_alb=False, p_top_as_attr=True)

    report = audit_pressure_refresh_sources([_entry("d02_start_wrfout", "d02", source)])
    payload = json.loads(report_to_json(report))

    assert payload["status"] == "failed"
    assert payload["summary"]["missing_count"] == 1
    assert payload["summary"]["d02_alb_blocker"] is True
    assert payload["summary"]["d02_alb_blocker_entries"] == ["d02_start_wrfout"]
    assert payload["summary"]["base_state_reconstruction_required_inputs"] == [
        "HT/HGT",
        "P_TOP",
        "C3F",
        "C4F",
        "C3H",
        "C4H",
    ]
    assert payload["summary"]["base_state_reconstruction_inputs_complete_count"] == 1
    assert payload["summary"]["base_state_reconstruction_inputs_complete_entries"] == [
        "d02_start_wrfout"
    ]
    assert payload["summary"]["base_state_reconstruction_missing_inputs_by_entry"] == {}
    assert payload["summary"]["direct_phb_present_count"] == 1
    assert payload["summary"]["direct_phb_present_entries"] == ["d02_start_wrfout"]
    assert payload["summary"]["direct_phb_shape_valid_count"] == 1
    assert payload["summary"]["direct_phb_shape_valid_entries"] == ["d02_start_wrfout"]
    assert payload["summary"]["direct_phb_shape_invalid_entries"] == []
    assert payload["summary"]["phb_reconstruction_required_inputs"] == [
        "HT/HGT",
        "P_TOP",
        "C3F",
        "C4F",
        "C3H",
        "C4H",
    ]
    assert payload["summary"]["phb_reconstruction_inputs_complete_count"] == 1
    assert payload["summary"]["phb_reconstruction_inputs_complete_entries"] == [
        "d02_start_wrfout"
    ]
    assert payload["summary"]["phb_reconstruction_missing_inputs_by_entry"] == {}
    assert payload["summary"]["full_provider_recommended"] is True
    assert payload["summary"]["full_provider_recommended_entries"] == [
        "d02_start_wrfout"
    ]
    assert (
        payload["summary"]["full_provider_recommended_source"]
        == "wrf_start_domain_phb_hypsometric_opt2_reconstruction"
    )
    assert payload["summary"]["direct_pb_t_init_path_available_count"] == 0
    assert payload["summary"]["direct_pb_t_init_path_available_entries"] == []
    assert payload["summary"]["base_state_reconstruction_required"] is True
    assert payload["summary"]["base_state_reconstruction_entries"] == ["d02_start_wrfout"]
    assert payload["summary"]["d02_base_state_reconstruction_required"] is True
    assert payload["summary"]["d02_base_state_reconstruction_entries"] == [
        "d02_start_wrfout"
    ]
    assert (
        payload["summary"]["recommended_next_source"]
        == "wrf_start_domain_base_state_reconstruction"
    )
    entry = payload["entries"][0]
    assert entry["status"] == "missing"
    assert entry["missing_names"] == ["ALB"]
    assert entry["p_top_present"] is True
    assert entry["p_top_source"] == "global_attribute"
    assert entry["alb_available"] is False
    assert entry["direct_alb_source"] is False
    assert entry["direct_phb_present"] is True
    assert entry["direct_phb_shape_valid"] is True
    assert entry["base_state_reconstruction_inputs_complete"] is True
    assert entry["missing_base_state_reconstruction_inputs"] == []
    assert entry["phb_reconstruction_inputs_complete"] is True
    assert entry["missing_phb_reconstruction_inputs"] == []
    assert entry["base_state_terrain_source"] == "HGT"
    assert entry["base_state_reconstruction_required"] is True
    assert entry["recommended_next_source"] == "wrf_start_domain_base_state_reconstruction"
    assert (
        entry["full_provider_recommended_source"]
        == "wrf_start_domain_phb_hypsometric_opt2_reconstruction"
    )
    assert entry["can_seed_pressure_refresh"] is False
    assert entry["suitable_for_start_time_truth"] is False


def test_audit_reports_missing_reconstruction_inputs_for_d02_start_time(
    tmp_path: Path,
) -> None:
    source = tmp_path / "wrfout_d02_2025-07-26_00:00:00"
    _write_source(source, include_alb=False, include_terrain=False, p_top_as_attr=True)

    payload = json.loads(
        report_to_json(audit_pressure_refresh_sources([_entry("d02_start_wrfout", "d02", source)]))
    )

    assert payload["status"] == "failed"
    assert payload["summary"]["d02_alb_blocker"] is True
    assert payload["summary"]["base_state_reconstruction_required"] is False
    assert payload["summary"]["d02_base_state_reconstruction_required"] is False
    assert payload["summary"]["recommended_next_source"] is None
    assert payload["summary"]["base_state_reconstruction_missing_inputs_by_entry"] == {
        "d02_start_wrfout": ["HT/HGT"]
    }
    assert payload["summary"]["phb_reconstruction_missing_inputs_by_entry"] == {
        "d02_start_wrfout": ["HT/HGT"]
    }
    assert payload["summary"]["full_provider_recommended_source"] is None
    entry = payload["entries"][0]
    assert entry["missing_names"] == ["ALB"]
    assert entry["base_state_reconstruction_inputs_complete"] is False
    assert entry["missing_base_state_reconstruction_inputs"] == ["HT/HGT"]
    assert entry["direct_phb_present"] is True
    assert entry["direct_phb_shape_valid"] is True
    assert entry["phb_reconstruction_inputs_complete"] is False
    assert entry["missing_phb_reconstruction_inputs"] == ["HT/HGT"]
    assert entry["base_state_terrain_source"] is None
    assert entry["base_state_reconstruction_input_presence"]["HT/HGT"] is False
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
    assert entry["full_provider_recommended_source"] is None


def test_audit_reports_optional_direct_pb_t_init_path(tmp_path: Path) -> None:
    source = tmp_path / "wrfinput_d01"
    _write_source(source, include_direct_pb_t_init=True)

    payload = json.loads(
        report_to_json(audit_pressure_refresh_sources([_entry("d01_wrfinput", "d01", source)]))
    )

    assert payload["summary"]["direct_pb_t_init_path_available_count"] == 1
    assert payload["summary"]["direct_pb_t_init_path_available_entries"] == [
        "d01_wrfinput"
    ]
    entry = payload["entries"][0]
    assert entry["direct_pb_t_init_path_available"] is True
    assert entry["missing_direct_pb_t_init_inputs"] == []
    assert entry["base_state_reconstruction_input_presence"]["PB"] is True
    assert entry["base_state_reconstruction_input_presence"]["T_INIT"] is True
    assert entry["base_state_reconstruction_input_presence"]["direct_PB_T_INIT"] is True


def test_audit_reports_nonexistent_candidate(tmp_path: Path) -> None:
    missing = tmp_path / "wrfinput_d02"

    payload = json.loads(
        report_to_json(audit_pressure_refresh_sources([_entry("missing_d02_input", "d02", missing)]))
    )

    assert payload["status"] == "failed"
    assert payload["summary"]["nonexistent_count"] == 1
    entry = payload["entries"][0]
    assert entry["status"] == "nonexistent"
    assert entry["missing_names"] == []
    assert entry["direct_alb_source"] is False
    assert entry["direct_phb_present"] is False
    assert entry["direct_phb_shape_valid"] is False
    assert entry["phb_reconstruction_inputs_complete"] is False
    assert entry["missing_phb_reconstruction_inputs"] == [
        "HT/HGT",
        "P_TOP",
        "C3F",
        "C4F",
        "C3H",
        "C4H",
    ]
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
    assert entry["full_provider_recommended_source"] is None
    assert entry["can_seed_pressure_refresh"] is False


def test_audit_reports_bad_shape_as_error(tmp_path: Path) -> None:
    source = tmp_path / "bad_shape.nc"
    _write_source(source, bad_alb_shape=True)

    payload = json.loads(report_to_json(audit_pressure_refresh_sources([_entry("bad", "d01", source)])))

    assert payload["status"] == "failed"
    assert payload["summary"]["error_count"] == 1
    entry = payload["entries"][0]
    assert entry["status"] == "error"
    assert entry["missing_names"] == []
    assert "ALB has shape (2, 3, 5), expected (2, 3, 4)" in entry["message"]
    assert entry["direct_alb_source"] is False
    assert entry["direct_phb_present"] is True
    assert entry["direct_phb_shape_valid"] is True
    assert entry["phb_reconstruction_inputs_complete"] is False
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
    assert entry["full_provider_recommended_source"] is None
    assert entry["can_seed_pressure_refresh"] is False


def test_audit_reports_direct_phb_shape_invalid_without_making_phb_required(
    tmp_path: Path,
) -> None:
    source = tmp_path / "bad_phb_shape.nc"
    _write_source(source, bad_phb_shape=True)

    payload = json.loads(
        report_to_json(audit_pressure_refresh_sources([_entry("bad_phb", "d01", source)]))
    )

    assert payload["status"] == "ok"
    assert payload["summary"]["direct_phb_present_count"] == 1
    assert payload["summary"]["direct_phb_shape_valid_count"] == 0
    assert payload["summary"]["direct_phb_shape_invalid_entries"] == ["bad_phb"]
    assert payload["summary"]["phb_reconstruction_inputs_complete_count"] == 1
    entry = payload["entries"][0]
    assert entry["status"] == "ok"
    assert entry["direct_phb_present"] is True
    assert entry["direct_phb_shape_valid"] is False
    assert entry["phb_reconstruction_inputs_complete"] is True
    assert entry["missing_phb_reconstruction_inputs"] == []
    assert (
        entry["full_provider_recommended_source"]
        == "wrf_start_domain_phb_hypsometric_opt2_reconstruction"
    )


def test_later_restart_is_not_clean_truth_without_explicit_manifest_flag(
    tmp_path: Path,
) -> None:
    source = tmp_path / "wrfrst_d02_later"
    _write_source(source)

    payload = json.loads(
        report_to_json(audit_pressure_refresh_sources([_entry("d02_later_restart", "d02", source)]))
    )

    entry = payload["entries"][0]
    assert entry["status"] == "ok"
    assert entry["direct_alb_source"] is True
    assert entry["direct_phb_present"] is True
    assert entry["direct_phb_shape_valid"] is True
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
    assert (
        entry["full_provider_recommended_source"]
        == "wrf_start_domain_phb_hypsometric_opt2_reconstruction"
    )
    assert entry["can_seed_pressure_refresh"] is True
    assert entry["suitable_for_start_time_truth"] is False


def test_later_restart_direct_alb_does_not_clear_d02_start_time_blocker(
    tmp_path: Path,
) -> None:
    start_source = tmp_path / "wrfout_d02_2025-07-26_00:00:00"
    restart_source = tmp_path / "wrfrst_d02_later"
    _write_source(start_source, include_alb=False, p_top_as_attr=True)
    _write_source(restart_source)

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_sources(
                [
                    _entry("d02_start_wrfout", "d02", start_source),
                    _entry("d02_later_restart", "d02", restart_source),
                ]
            )
        )
    )

    assert payload["status"] == "failed"
    assert payload["summary"]["d02_alb_blocker"] is True
    assert payload["summary"]["d02_alb_blocker_entries"] == ["d02_start_wrfout"]
    assert payload["summary"]["d02_base_state_reconstruction_required"] is True
    entries = {entry["name"]: entry for entry in payload["entries"]}
    assert entries["d02_start_wrfout"]["base_state_reconstruction_required"] is True
    assert entries["d02_start_wrfout"]["suitable_for_start_time_truth"] is False
    assert entries["d02_later_restart"]["direct_alb_source"] is True
    assert entries["d02_later_restart"]["direct_phb_present"] is True
    assert entries["d02_later_restart"]["direct_phb_shape_valid"] is True
    assert entries["d02_later_restart"]["suitable_for_start_time_truth"] is False


def test_pressure_diagnostics_missing_optional_source_reports_not_requested(
    tmp_path: Path,
) -> None:
    source = tmp_path / "wrfinput_d01"
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate.nc"
    _write_source(source)
    reference_p = np.full((2, 2, 4), 100.0)
    candidate_p = reference_p.copy()
    candidate_p[:, :, 3] += 2.0
    _write_wrfout(reference, {"P": reference_p})
    _write_wrfout(candidate, {"P": candidate_p}, attrs=_pressure_attrs())

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_sources(
                [_entry("d01_wrfinput", "d01", source)],
                pressure_diagnostic=PressureDiagnosticRequest(
                    reference_end_path=str(reference),
                    candidate_path=str(candidate),
                ),
            )
        )
    )

    diagnostics = payload["pressure_diagnostics"]
    assert diagnostics["status"] == "computed"
    assert diagnostics["candidate_model_pass"] == "not_applicable"
    assert diagnostics["evolution_from_source_start"]["status"] == "not_requested"
    assert diagnostics["evolution_from_source_start"]["diagnostic_only"] is True
    assert diagnostics["evolution_from_source_start"]["source_start"] is None
    assert diagnostics["target_region_p"]["aggregate"]["valid_count"] == 4
    assert diagnostics["target_region_p"]["per_vertical_level"][0]["valid_count"] == 2


def test_pressure_diagnostics_ranks_worst_vertical_levels(tmp_path: Path) -> None:
    source = tmp_path / "wrfinput_d01"
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate.nc"
    start = tmp_path / "source_start.nc"
    _write_source(source)
    reference_p = np.full((3, 2, 4), 100.0)
    candidate_p = reference_p.copy()
    candidate_p[0, :, 3] += 1.0
    candidate_p[1, :, 3] += 5.0
    candidate_p[2, :, 3] -= 9.0
    source_p = reference_p - 10.0
    _write_wrfout(reference, {"P": reference_p})
    _write_wrfout(candidate, {"P": candidate_p}, attrs=_pressure_attrs())
    _write_wrfout(start, {"P": source_p})

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_sources(
                [_entry("d01_wrfinput", "d01", source)],
                pressure_diagnostic=PressureDiagnosticRequest(
                    reference_end_path=str(reference),
                    candidate_path=str(candidate),
                    source_start_path=str(start),
                    worst_level_count=2,
                ),
            )
        )
    )

    target_p = payload["pressure_diagnostics"]["target_region_p"]
    assert [level["level_index"] for level in target_p["worst_levels_by_rmse"]] == [
        2,
        1,
    ]
    assert [
        level["level_index"] for level in target_p["worst_levels_by_bias_magnitude"]
    ] == [2, 1]
    levels = target_p["per_vertical_level"]
    assert math.isclose(levels[2]["rmse"], 9.0)
    assert math.isclose(levels[2]["mean_bias"], -9.0)
    assert math.isclose(levels[2]["mean_bias_magnitude"], 9.0)

    evolution = payload["pressure_diagnostics"]["evolution_from_source_start"]
    assert evolution["status"] == "available"
    assert evolution["candidate_model_pass"] == "not_applicable"
    assert math.isclose(
        evolution["reference_end_vs_source_start"]["mean_bias"],
        10.0,
    )
    assert len(evolution["per_vertical_level"]) == 3


def test_pressure_diagnostics_disposition_is_diagnostic_only(
    tmp_path: Path,
) -> None:
    source = tmp_path / "wrfinput_d01"
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate.nc"
    _write_source(source)
    p = np.full((1, 2, 4), 100.0)
    _write_wrfout(reference, {"P": p})
    _write_wrfout(candidate, {"P": p}, attrs=_pressure_attrs())

    payload = json.loads(
        report_to_json(
            audit_pressure_refresh_sources(
                [_entry("d01_wrfinput", "d01", source)],
                pressure_diagnostic=PressureDiagnosticRequest(
                    reference_end_path=str(reference),
                    candidate_path=str(candidate),
                ),
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["diagnostic_only"] is True
    assert payload["summary"]["candidate_model_pass"] == "not_applicable"
    assert payload["summary"]["pressure_diagnostics_status"] == "computed"
    assert payload["pressure_diagnostics"]["diagnostic_only"] is True
    assert payload["pressure_diagnostics"]["candidate_model_pass"] == "not_applicable"
    assert (
        payload["pressure_diagnostics"]["metadata"]["disposition"][
            "candidate_model_pass"
        ]
        == "not_applicable"
    )


def test_direct_script_help_invocation_uses_project_import_path() -> None:
    environment = {
        **os.environ,
        "UV_CACHE_DIR": ".uv-cache",
        "UV_PYTHON_INSTALL_DIR": ".uv-python",
    }

    completed = subprocess.run(
        [
            "uv",
            "run",
            "python",
            "tools/audit_pressure_refresh_sources.py",
            "--help",
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0
    assert "--manifest" in completed.stdout
    assert "ModuleNotFoundError" not in completed.stderr


def test_main_reads_manifest_and_writes_json(tmp_path: Path) -> None:
    source = tmp_path / "wrfinput_d01"
    output = tmp_path / "audit.json"
    manifest = tmp_path / "manifest.json"
    _write_source(source)
    manifest.write_text(
        json.dumps(
            {
                "entries": [
                    {
                        "name": "d01_wrfinput",
                        "domain": "d01",
                        "path": str(source),
                        "mass_nx": 4,
                        "mass_ny": 3,
                        "mass_nz": 2,
                        "full_nz": 3,
                        "time_index": 0,
                        "clean_source": True,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    exit_code = audit_main(
        ["--manifest", str(manifest), "--output", str(output), "--pretty"]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["entries"][0]["suitable_for_start_time_truth"] is True
