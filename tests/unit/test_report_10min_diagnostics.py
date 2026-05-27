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
