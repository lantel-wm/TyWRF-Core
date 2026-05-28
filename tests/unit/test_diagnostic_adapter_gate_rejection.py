from __future__ import annotations

import json
from pathlib import Path

import netCDF4
import numpy as np

from tools.compare_wrfout import STRICT_CORE_VARIABLES
from tools.cycle_gate import main as cycle_gate_main


START = "2025-07-26_00:00:00"
END_00_10 = "2025-07-26_00:10:00"
END_FILE = f"wrfout_d02_{END_00_10}"
ADAPTER_KIND = "selected_field_diagnostic_adapter_v0"

D72_DIAGNOSTIC_ADAPTER_REPORT = {
    "status": "selected_field_diagnostic_adapter_reported",
    "candidate_kind": ADAPTER_KIND,
    "diagnostic_only": True,
    "gate_candidate": False,
    "integrator_output": False,
    "diagnostic_adapter_opt_in": True,
    "diagnostic_adapter_ok": True,
    "diagnostic_adapter_diagnostic_only": True,
    "diagnostic_adapter_gate_candidate": False,
    "diagnostic_adapter_integrator_output": False,
    "diagnostic_adapter_writes_netcdf": False,
    "diagnostic_adapter_writes_candidate": False,
    "diagnostic_adapter_called_d68_exchange": True,
    "diagnostic_adapter_called_d69_recompute": True,
    "diagnostic_adapter_integration_status": "staging_report_only_no_gate_no_integrator",
    "diagnostic_adapter_c3h_count": 2,
    "diagnostic_adapter_c4h_count": 2,
    "diagnostic_adapter_exposed_mass_cell_count": 25,
    "diagnostic_adapter_recomputed_point_count": 50,
    "diagnostic_adapter_invalid_point_count": 0,
}


def _bool_attr(value: bool) -> str:
    return "true" if value else "false"


def _attrs_from_adapter_report(report: dict[str, object]) -> dict[str, object]:
    assert report["status"] == "selected_field_diagnostic_adapter_reported"
    assert report["candidate_kind"] == ADAPTER_KIND
    assert report["diagnostic_only"] is True
    assert report["gate_candidate"] is False
    assert report["integrator_output"] is False
    assert report["diagnostic_adapter_diagnostic_only"] is True
    assert report["diagnostic_adapter_gate_candidate"] is False
    assert report["diagnostic_adapter_integrator_output"] is False
    assert report["diagnostic_adapter_writes_candidate"] is False
    assert report["diagnostic_adapter_writes_netcdf"] is False
    assert (
        report["diagnostic_adapter_integration_status"]
        == "staging_report_only_no_gate_no_integrator"
    )

    return {
        "DX": 2000.0,
        "DY": 2000.0,
        "TYWRF_DIAGNOSTIC_ONLY": _bool_attr(report["diagnostic_only"]),
        "TYWRF_GATE_CANDIDATE": _bool_attr(report["gate_candidate"]),
        "TYWRF_INTEGRATOR_OUTPUT": _bool_attr(report["integrator_output"]),
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": report["candidate_kind"],
        "TYWRF_DIAGNOSTIC_ADAPTER_OPT_IN": _bool_attr(
            report["diagnostic_adapter_opt_in"]
        ),
        "TYWRF_DIAGNOSTIC_ADAPTER_OK": _bool_attr(report["diagnostic_adapter_ok"]),
        "TYWRF_DIAGNOSTIC_ADAPTER_DIAGNOSTIC_ONLY": _bool_attr(
            report["diagnostic_adapter_diagnostic_only"]
        ),
        "TYWRF_DIAGNOSTIC_ADAPTER_GATE_CANDIDATE": _bool_attr(
            report["diagnostic_adapter_gate_candidate"]
        ),
        "TYWRF_DIAGNOSTIC_ADAPTER_INTEGRATOR_OUTPUT": _bool_attr(
            report["diagnostic_adapter_integrator_output"]
        ),
        "TYWRF_DIAGNOSTIC_ADAPTER_WRITES_NETCDF": _bool_attr(
            report["diagnostic_adapter_writes_netcdf"]
        ),
        "TYWRF_DIAGNOSTIC_ADAPTER_WRITES_CANDIDATE": _bool_attr(
            report["diagnostic_adapter_writes_candidate"]
        ),
        "TYWRF_DIAGNOSTIC_ADAPTER_CALLED_D68_EXCHANGE": _bool_attr(
            report["diagnostic_adapter_called_d68_exchange"]
        ),
        "TYWRF_DIAGNOSTIC_ADAPTER_CALLED_D69_RECOMPUTE": _bool_attr(
            report["diagnostic_adapter_called_d69_recompute"]
        ),
        "TYWRF_DIAGNOSTIC_ADAPTER_INTEGRATION_STATUS": report[
            "diagnostic_adapter_integration_status"
        ],
        "TYWRF_DIAGNOSTIC_ADAPTER_C3H_COUNT": report[
            "diagnostic_adapter_c3h_count"
        ],
        "TYWRF_DIAGNOSTIC_ADAPTER_C4H_COUNT": report[
            "diagnostic_adapter_c4h_count"
        ],
        "TYWRF_DIAGNOSTIC_ADAPTER_EXPOSED_MASS_CELL_COUNT": report[
            "diagnostic_adapter_exposed_mass_cell_count"
        ],
        "TYWRF_DIAGNOSTIC_ADAPTER_RECOMPUTED_POINT_COUNT": report[
            "diagnostic_adapter_recomputed_point_count"
        ],
        "TYWRF_DIAGNOSTIC_ADAPTER_INVALID_POINT_COUNT": report[
            "diagnostic_adapter_invalid_point_count"
        ],
        "TYWRF_CANDIDATE_MESSAGE": (
            "Diagnostic-only selected-field base-state adapter report; not a "
            "validation gate pass or normal integrator output; staging buffers only."
        ),
    }


def _tc_fields() -> dict[str, np.ndarray]:
    shape = (5, 5)
    latitude = np.repeat(np.arange(10.0, 15.0)[:, None], shape[1], axis=1)
    longitude = np.repeat(np.arange(120.0, 125.0)[None, :], shape[0], axis=0)
    slp = np.full(shape, 1000.0)
    slp[2, 2] = 950.0
    u10 = np.zeros(shape)
    v10 = np.zeros(shape)
    u10[1, 1] = 30.0
    v10[1, 1] = 40.0
    return {
        "XLAT": latitude,
        "XLONG": longitude,
        "SLP": slp,
        "U10": u10,
        "V10": v10,
    }


def _write_wrfout(
    path: Path,
    *,
    attrs: dict[str, object] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    shape = (5, 5)
    fields = {name: np.full(shape, 10.0) for name in STRICT_CORE_VARIABLES}
    fields.update(_tc_fields())

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", 1)
        dataset.createDimension("DateStrLen", 19)
        dataset.createDimension("south_north", shape[0])
        dataset.createDimension("west_east", shape[1])

        for name, value in (attrs or {}).items():
            dataset.setncattr(name, value)

        times = dataset.createVariable("Times", "S1", ("Time", "DateStrLen"))
        times[0, :] = np.array(list(END_00_10), dtype="S1")

        for name, values in fields.items():
            variable = dataset.createVariable(
                name,
                "f8",
                ("Time", "south_north", "west_east"),
            )
            variable[0, :, :] = values


def test_d72_diagnostic_adapter_report_is_rejected_by_strict_0010_gate(
    tmp_path: Path,
    capsys,
) -> None:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _write_wrfout(reference_dir / END_FILE)
    _write_wrfout(
        candidate_dir / END_FILE,
        attrs=_attrs_from_adapter_report(D72_DIAGNOSTIC_ADAPTER_REPORT),
    )

    report_path = tmp_path / "strict_gate.json"
    exit_code = cycle_gate_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--candidate-dir",
            str(candidate_dir),
            "--start",
            START,
            "--end",
            END_00_10,
            "--interval-minutes",
            "10",
            "--output",
            str(report_path),
            "--pretty",
        ]
    )
    capsys.readouterr()

    payload = json.loads(report_path.read_text(encoding="utf-8"))
    cycle = payload["cycles"][0]
    diagnostics = {metric["name"]: metric for metric in cycle["diagnostics"]}
    metadata = diagnostics["candidate_metadata"]

    assert exit_code == 1
    assert payload["status"] == "failed"
    assert payload["summary"] == {"total": 1, "passed": 0, "failed": 1}
    assert payload["interval_minutes"] == 10
    assert cycle["start_time"] == START
    assert cycle["end_time"] == END_00_10
    assert "2025-07-26_00:20:00" not in json.dumps(payload)

    assert metadata["status"] == "failed"
    assert "TYWRF_DIAGNOSTIC_ONLY=true" in metadata["message"]
    assert "TYWRF_GATE_CANDIDATE=false" in metadata["message"]
    assert "TYWRF_INTEGRATOR_OUTPUT=false" in metadata["message"]
    assert f"TYWRF_CANDIDATE_KIND={ADAPTER_KIND}" in metadata["message"]
    assert "TYWRF_DIAGNOSTIC_ADAPTER_OK=true" in metadata["message"]
    assert (
        "TYWRF_DIAGNOSTIC_ADAPTER_INTEGRATION_STATUS="
        "staging_report_only_no_gate_no_integrator"
    ) in metadata["message"]

    assert {field["status"] for field in cycle["fields"]} == {"passed"}
    assert {
        metric["status"]
        for name, metric in diagnostics.items()
        if name != "candidate_metadata"
    } == {"passed"}
    assert payload["first_failure"]["cycle_index"] == 1
    assert payload["first_failure"]["end_time"] == END_00_10
    assert payload["first_failure"]["field"] is None
    assert payload["first_failure"]["diagnostic"] == "candidate_metadata"
