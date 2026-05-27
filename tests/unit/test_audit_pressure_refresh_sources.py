import json
from pathlib import Path

import netCDF4

from tools.audit_pressure_refresh_sources import (
    SourceEntry,
    audit_pressure_refresh_sources,
    main as audit_main,
    report_to_json,
)


def _write_vector(dataset: netCDF4.Dataset, name: str, length: int) -> None:
    variable = dataset.createVariable(name, "f4", ("Time", f"{name}_dim"))
    variable[0, :] = list(range(length))


def _write_source(
    path: Path,
    *,
    include_alb: bool = True,
    bad_alb_shape: bool = False,
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

        if include_alb:
            west_east_dim = "west_east_bad" if bad_alb_shape else "west_east"
            alb = dataset.createVariable(
                "ALB",
                "f4",
                ("Time", "bottom_top", "south_north", west_east_dim),
            )
            alb[0, :, :, :] = 1.0


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
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
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
    assert entry["base_state_reconstruction_required"] is True
    assert entry["recommended_next_source"] == "wrf_start_domain_base_state_reconstruction"
    assert entry["can_seed_pressure_refresh"] is False
    assert entry["suitable_for_start_time_truth"] is False


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
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
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
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
    assert entry["can_seed_pressure_refresh"] is False


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
    assert entry["base_state_reconstruction_required"] is False
    assert entry["recommended_next_source"] is None
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
    assert entries["d02_later_restart"]["suitable_for_start_time_truth"] is False


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
