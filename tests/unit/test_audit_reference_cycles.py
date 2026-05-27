import json
from pathlib import Path

import netCDF4

from tools.audit_reference_cycles import (
    SLP_CANDIDATE_VARIABLES,
    TARGET_VARIABLES,
    audit_reference_cycles,
    main as audit_main,
    report_to_table,
)


def _write_wrfout(
    path: Path,
    variables: tuple[str, ...],
    *,
    dx: float | None = 10000.0,
    dy: float | None = 10000.0,
) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        if dx is not None:
            dataset.setncattr("DX", dx)
        if dy is not None:
            dataset.setncattr("DY", dy)
        for variable in dict.fromkeys(variables):
            dataset.createVariable(variable, "f4", ("Time",))


def _write_cycle_set(
    directory: Path,
    times: tuple[str, ...],
    *,
    domains: tuple[str, ...] = ("d02",),
    variables: tuple[str, ...] = (*TARGET_VARIABLES, "SLP", "PSFC"),
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for domain in domains:
        for valid_time in times:
            _write_wrfout(directory / f"wrfout_{domain}_{valid_time}", variables)


def test_audit_reports_complete_synthetic_minute_directory_for_both_domains(
    tmp_path: Path,
) -> None:
    _write_cycle_set(
        tmp_path,
        (
            "2025-07-26_00:00:00",
            "2025-07-26_00:10:00",
            "2025-07-26_00:20:00",
        ),
        domains=("d01", "d02"),
    )

    report = audit_reference_cycles(
        tmp_path,
        start="2025-07-26_00:00:00",
        end="2025-07-26_00:20:00",
        interval_minutes=10,
    )

    assert report.status == "ok"
    assert report.domains == ["d01", "d02"]
    assert report.interval_minutes == 10
    assert report.cycle_count == 6
    assert report.missing_files == []
    assert report.missing_variables == {}
    assert report.missing_grid_spacing == {}
    assert all(candidates == ["SLP"] for candidates in report.available_slp_candidates.values())
    assert all(
        candidates == ["PSFC"]
        for candidates in report.available_pressure_proxy_candidates.values()
    )
    assert all(names == ["DX", "DY"] for names in report.available_grid_spacing.values())


def test_audit_reports_missing_minute_files(tmp_path: Path) -> None:
    _write_cycle_set(
        tmp_path,
        (
            "2025-07-26_00:00:00",
            "2025-07-26_00:20:00",
        ),
        domains=("d01", "d02"),
    )

    report = audit_reference_cycles(
        tmp_path,
        start="2025-07-26_00:00:00",
        end="2025-07-26_00:20:00",
        interval_minutes=10,
    )

    assert report.status == "failed"
    assert report.cycle_count == 6
    assert report.missing_files == [
        "d01/wrfout_d01_2025-07-26_00:10:00",
        "d02/wrfout_d02_2025-07-26_00:10:00",
    ]
    assert report.cycles[1].status == "missing_file"


def test_audit_reports_missing_variables(tmp_path: Path) -> None:
    _write_cycle_set(tmp_path, ("2025-07-26_00:00:00",))
    without_v = tuple(variable for variable in TARGET_VARIABLES if variable != "V")
    _write_wrfout(tmp_path / "wrfout_d02_2025-07-26_00:10:00", (*without_v, "SLP", "PSFC"))

    report = audit_reference_cycles(
        tmp_path,
        start="2025-07-26_00:00:00",
        end="2025-07-26_00:10:00",
        interval_minutes=10,
        domain="d02",
    )

    filename = "d02/wrfout_d02_2025-07-26_00:10:00"
    assert report.status == "failed"
    assert report.missing_variables == {filename: ["V"]}
    assert report.cycles[1].status == "missing_variables"


def test_audit_reports_missing_dx_dy(tmp_path: Path) -> None:
    _write_wrfout(
        tmp_path / "wrfout_d02_2025-07-26_00:00:00",
        (*TARGET_VARIABLES, "SLP", "PSFC"),
        dx=None,
    )

    report = audit_reference_cycles(
        tmp_path,
        start="2025-07-26_00:00:00",
        end="2025-07-26_00:00:00",
        domain="d02",
    )

    filename = "d02/wrfout_d02_2025-07-26_00:00:00"
    assert report.status == "failed"
    assert report.missing_grid_spacing == {filename: ["DX"]}
    assert report.available_grid_spacing[filename] == ["DY"]
    assert report.cycles[0].status == "missing_grid_spacing"


def test_audit_identifies_slp_candidates(tmp_path: Path) -> None:
    _write_cycle_set(
        tmp_path,
        ("2025-07-26_00:00:00",),
        variables=(*TARGET_VARIABLES, "SLP", "MSLP", "PSFC"),
    )

    report = audit_reference_cycles(
        tmp_path,
        start="2025-07-26_00:00:00",
        end="2025-07-26_00:00:00",
        domain="d02",
        slp_candidates=SLP_CANDIDATE_VARIABLES,
    )

    filename = "d02/wrfout_d02_2025-07-26_00:00:00"
    assert report.status == "ok"
    assert report.available_slp_candidates[filename] == ["SLP", "MSLP"]
    assert report.available_pressure_proxy_candidates[filename] == ["PSFC"]
    assert report.cycles_without_slp_candidates == []


def test_audit_keeps_psfc_proxy_separate_from_real_slp_candidates(tmp_path: Path) -> None:
    _write_cycle_set(
        tmp_path,
        ("2025-07-26_00:00:00",),
        variables=(*TARGET_VARIABLES, "PSFC"),
    )

    report = audit_reference_cycles(
        tmp_path,
        start="2025-07-26_00:00:00",
        end="2025-07-26_00:00:00",
        domain="d02",
        slp_candidates=SLP_CANDIDATE_VARIABLES,
    )

    filename = "d02/wrfout_d02_2025-07-26_00:00:00"
    assert report.status == "failed"
    assert report.cycles[0].status == "missing_slp_candidates"
    assert report.available_slp_candidates[filename] == []
    assert report.available_pressure_proxy_candidates[filename] == ["PSFC"]
    assert report.cycles_without_slp_candidates == [filename]


def test_audit_main_writes_json_and_table_reports(tmp_path: Path) -> None:
    _write_cycle_set(tmp_path, ("2025-07-26_00:00:00",))
    output = tmp_path / "audit.json"

    exit_code = audit_main(
        [
            "--reference-dir",
            str(tmp_path),
            "--domain",
            "d02",
            "--start",
            "2025-07-26_00:00:00",
            "--end",
            "2025-07-26_00:00:00",
            "--output",
            str(output),
            "--pretty",
        ]
    )

    assert exit_code == 0
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["status"] == "ok"
    assert payload["cycle_count"] == 1

    table = report_to_table(
        audit_reference_cycles(
            tmp_path,
            start="2025-07-26_00:00:00",
            end="2025-07-26_00:00:00",
            domain="d02",
        )
    )
    assert "valid_time" in table
    assert "status=ok cycle_count=1" in table
    assert "pressure_proxies" in table
