from pathlib import Path

from tools.run_6h_cycle_test import (
    build_cycle_plan,
    format_wrf_time,
    main as cycle_main,
    parse_wrf_time,
    wrfout_filename,
)


def test_wrf_time_parsing_and_filename_formatting() -> None:
    value = parse_wrf_time("2025-07-26_00:00:00")

    assert format_wrf_time(value) == "2025-07-26_00:00:00"
    assert wrfout_filename("d02", value) == "wrfout_d02_2025-07-26_00:00:00"
    assert parse_wrf_time("2025-07-26T00:00:00") == value


def test_build_cycle_plan_resolves_d01_d02_start_and_end_files(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    reference_dir.mkdir()
    candidate_dir.mkdir()

    (reference_dir / "wrfout_d01_2025-07-26_00:00:00").touch()
    (reference_dir / "wrfout_d01_2025-07-26_06:00:00").touch()
    (candidate_dir / "wrfout_d01_2025-07-26_00:00:00").touch()
    (candidate_dir / "wrfout_d02_2025-07-26_06:00:00").touch()

    plan = build_cycle_plan(
        reference_dir,
        candidate_dir,
        "2025-07-26_00:00:00",
        domains=("d01", "d02"),
    )

    assert plan.status == "dry_run"
    assert plan.start_time == "2025-07-26_00:00:00"
    assert plan.end_time == "2025-07-26_06:00:00"
    assert [domain.domain for domain in plan.domains] == ["d01", "d02"]
    assert plan.domains[0].reference_start.endswith("wrfout_d01_2025-07-26_00:00:00")
    assert plan.domains[0].candidate_start.endswith("wrfout_d01_2025-07-26_00:00:00")
    assert plan.domains[0].reference_start_exists is True
    assert plan.domains[0].reference_end_exists is True
    assert plan.domains[0].candidate_start_exists is True
    assert plan.domains[0].candidate_end_exists is False
    assert plan.domains[1].reference_end.endswith("wrfout_d02_2025-07-26_06:00:00")
    assert plan.domains[1].candidate_end_exists is True
    assert plan.diagnostics["tc"]["status"] == "pending"


def test_cycle_main_dry_run_writes_json_plan(tmp_path: Path) -> None:
    output = tmp_path / "cycle-plan.json"

    exit_code = cycle_main(
        [
            "--reference-dir",
            str(tmp_path / "reference"),
            "--candidate-dir",
            str(tmp_path / "candidate"),
            "--start",
            "2025-07-26_00:00:00",
            "--dry-run",
            "--output",
            str(output),
            "--pretty",
        ]
    )

    assert exit_code == 0
    text = output.read_text(encoding="utf-8")
    assert '"status": "dry_run"' in text
    assert '"domain": "d01"' in text
    assert '"domain": "d02"' in text
