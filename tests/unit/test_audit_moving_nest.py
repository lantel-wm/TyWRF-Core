import json
import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.audit_moving_nest import (
    KROSA_END,
    KROSA_START,
    audit_moving_nest,
    audit_reference_moving_nest,
    main as audit_main,
    report_to_json,
    resolve_moving_nest_files,
)


def _write_wrfout(
    path: Path,
    *,
    i_parent_start: int,
    j_parent_start: int,
    dx: float = 2000.0,
    dy: float = 2000.0,
    lat_offset: float = 0.0,
    lon_offset: float = 0.0,
    hgt_offset: float = 0.0,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    y = np.arange(4, dtype=np.float64).reshape(4, 1)
    x = np.arange(5, dtype=np.float64).reshape(1, 5)
    lat = 20.0 + y * 0.1 + x * 0.01 + lat_offset
    lon = 140.0 + y * 0.05 + x * 0.02 + lon_offset
    hgt = 10.0 + y * 2.0 + x * 3.0 + hgt_offset

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        dataset.createDimension("south_north", lat.shape[0])
        dataset.createDimension("west_east", lat.shape[1])
        dataset.setncattr("I_PARENT_START", i_parent_start)
        dataset.setncattr("J_PARENT_START", j_parent_start)
        dataset.setncattr("DX", dx)
        dataset.setncattr("DY", dy)
        for name, values in {"XLAT": lat, "XLONG": lon, "HGT": hgt}.items():
            variable = dataset.createVariable(name, "f8", ("Time", "south_north", "west_east"))
            variable[0, :, :] = values


def test_audit_reports_explicit_changed_files_as_diagnostic_json(tmp_path: Path) -> None:
    start = tmp_path / "wrfout_d02_start"
    end = tmp_path / "wrfout_d02_end"
    _write_wrfout(start, i_parent_start=10, j_parent_start=20)
    _write_wrfout(
        end,
        i_parent_start=13,
        j_parent_start=18,
        lat_offset=0.1,
        lon_offset=0.2,
        hgt_offset=4.0,
    )

    report = audit_moving_nest(start, end)
    payload = json.loads(report_to_json(report, pretty=True))

    assert payload["status"] == "computed"
    assert payload["summary"]["start_pose"]["i_parent_start"] == 10
    assert payload["summary"]["end_pose"]["j_parent_start"] == 18
    assert payload["summary"]["parent_delta"] == {"i": 3, "j": -2}
    assert payload["summary"]["child_cell_delta"] == {"i": 15, "j": -10}
    assert payload["summary"]["moved"] is True
    assert payload["summary"]["resolution"]["status"] == "ok"
    assert payload["summary"]["resolution"]["d02_2km"] is True
    assert payload["summary"]["hgt_delta"]["status"] == "computed"
    assert math.isclose(payload["summary"]["hgt_delta"]["rms_m"], 4.0)
    assert math.isclose(payload["summary"]["hgt_delta"]["max_abs_m"], 4.0)
    assert payload["summary"]["drift_km"] > 0.0
    flags = payload["summary"]["validation_flags"]["expected_krosa_first_10min"]
    assert flags["applicable"] is False
    assert flags["matches_expected_krosa_first_10min"] is None


def test_reference_dir_path_resolution_and_krosa_flags(tmp_path: Path) -> None:
    reference_dir = tmp_path / "reference"
    start = reference_dir / f"wrfout_d02_{KROSA_START}"
    end = reference_dir / f"wrfout_d02_{KROSA_END}"
    _write_wrfout(start, i_parent_start=114, j_parent_start=96)
    _write_wrfout(end, i_parent_start=126, j_parent_start=103, hgt_offset=1.0)

    resolved_start, resolved_end, domain, start_time, end_time = resolve_moving_nest_files(
        reference_dir=reference_dir,
        domain="d02",
        start=KROSA_START,
        end=KROSA_END,
    )
    report = audit_reference_moving_nest(reference_dir, start=KROSA_START, end=KROSA_END)
    payload = json.loads(report_to_json(report))

    assert resolved_start == start
    assert resolved_end == end
    assert domain == "d02"
    assert start_time == KROSA_START
    assert end_time == KROSA_END
    assert payload["summary"]["parent_delta"] == {"i": 12, "j": 7}
    assert payload["summary"]["child_cell_delta"] == {"i": 60, "j": 35}
    flags = payload["summary"]["validation_flags"]["expected_krosa_first_10min"]
    assert flags["applicable"] is True
    assert flags["start_pose_matches"] is True
    assert flags["end_pose_matches"] is True
    assert flags["parent_delta_matches"] is True
    assert flags["matches_expected_krosa_first_10min"] is True


def test_resolution_check_marks_non_2km_d02(tmp_path: Path) -> None:
    start = tmp_path / "start.nc"
    end = tmp_path / "end.nc"
    _write_wrfout(start, i_parent_start=1, j_parent_start=1, dx=3000.0, dy=2000.0)
    _write_wrfout(end, i_parent_start=1, j_parent_start=1, dx=3000.0, dy=2000.0)

    report = audit_moving_nest(start, end)
    payload = json.loads(report_to_json(report))

    assert payload["status"] == "computed_with_flags"
    assert payload["summary"]["resolution"]["status"] == "not_2km"
    assert payload["summary"]["resolution"]["d02_2km"] is False
    assert payload["summary"]["resolution"]["start_dx_m"] == 3000.0


def test_unmoved_case_reports_zero_delta_and_zero_hgt_change(tmp_path: Path) -> None:
    start = tmp_path / "start.nc"
    end = tmp_path / "end.nc"
    _write_wrfout(start, i_parent_start=5, j_parent_start=6)
    _write_wrfout(end, i_parent_start=5, j_parent_start=6)

    payload = json.loads(report_to_json(audit_moving_nest(start, end)))

    assert payload["summary"]["parent_delta"] == {"i": 0, "j": 0}
    assert payload["summary"]["child_cell_delta"] == {"i": 0, "j": 0}
    assert payload["summary"]["moved"] is False
    assert math.isclose(payload["summary"]["hgt_delta"]["rms_m"], 0.0)
    assert math.isclose(payload["summary"]["hgt_delta"]["max_abs_m"], 0.0)
    assert math.isclose(payload["summary"]["drift_km"], 0.0)


def test_log_events_parse_first_segment_moves(tmp_path: Path) -> None:
    start = tmp_path / f"wrfout_d02_{KROSA_START}"
    end = tmp_path / f"wrfout_d02_{KROSA_END}"
    log_file = tmp_path / "rsl.out.0000"
    _write_wrfout(start, i_parent_start=114, j_parent_start=96)
    _write_wrfout(end, i_parent_start=116, j_parent_start=97)
    log_file.write_text(
        "\n".join(
            [
                " 2025-07-26_00:00:40 vortex center (in nest x and y):    167.1478       142.9944",
                " 2025-07-26_00:00:40 grid   center (in nest x and y):    105.5000       105.5000",
                " 2025-07-26_00:00:40 disp          :    6.000000       6.000000",
                " 2025-07-26_00:00:40 move (rel cd) :            1           1",
                "  moving            2           1           1",
                " 2025-07-26_00:01:20 move (rel cd) :            1           0",
                "  moving            2           1           0",
                " 2025-07-26_00:10:40 move (rel cd) :            9           9",
                "  moving            2           9           9",
            ]
        ),
        encoding="utf-8",
    )

    payload = json.loads(
        report_to_json(
            audit_moving_nest(
                start,
                end,
                start_time=KROSA_START,
                end_time=KROSA_END,
                log_file=log_file,
            )
        )
    )

    assert payload["log_events"]["status"] == "available"
    assert payload["log_events"]["event_count"] == 2
    assert payload["log_events"]["applied_move_count"] == 2
    assert payload["log_events"]["net_applied_parent_delta"] == {"i": 2, "j": 1}
    assert payload["log_events"]["net_applied_matches_parent_delta"] is True
    assert payload["log_events"]["events"][0]["vortex_center_nest_xy"]["x"] == 167.1478


def test_main_writes_json_for_explicit_files(tmp_path: Path) -> None:
    start = tmp_path / "wrfout_d02_start"
    end = tmp_path / "wrfout_d02_end"
    output = tmp_path / "audit.json"
    _write_wrfout(start, i_parent_start=1, j_parent_start=2)
    _write_wrfout(end, i_parent_start=2, j_parent_start=2)

    exit_code = audit_main(
        [
            "--start-file",
            str(start),
            "--end-file",
            str(end),
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["summary"]["parent_delta"] == {"i": 1, "j": 0}
    assert payload["log_events"]["status"] == "not_requested"
