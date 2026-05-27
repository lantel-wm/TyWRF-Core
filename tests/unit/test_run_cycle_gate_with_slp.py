import json
from pathlib import Path

import netCDF4
import numpy as np

from tools.run_cycle_gate_with_slp import (
    main as run_main,
    report_to_dict,
    run_cycle_gate_with_slp,
)


START = "2025-07-26_00:00:00"
END = "2025-07-26_06:00:00"
END2 = "2025-07-26_12:00:00"
END_FILE = "wrfout_d02_2025-07-26_06:00:00"
END_FILE2 = "wrfout_d02_2025-07-26_12:00:00"


def _tc_fields(shape: tuple[int, int]) -> dict[str, np.ndarray]:
    latitude = np.repeat(np.arange(10.0, 10.0 + shape[0])[:, None], shape[1], axis=1)
    longitude = np.repeat(np.arange(120.0, 120.0 + shape[1])[None, :], shape[0], axis=0)
    psfc = np.full(shape, 100000.0)
    psfc[2, 2] = 95000.0
    u10 = np.zeros(shape)
    v10 = np.zeros(shape)
    u10[1, 1] = 30.0
    v10[1, 1] = 40.0
    return {
        "XLAT": latitude,
        "XLONG": longitude,
        "PSFC": psfc,
        "U10": u10,
        "V10": v10,
    }


def _write_wrfout(
    path: Path,
    *,
    omit: set[str] | None = None,
    field_offset: float = 0.0,
) -> None:
    omit = omit or set()
    path.parent.mkdir(parents=True, exist_ok=True)

    horizontal_shape = (5, 5)
    mass_shape = (1, 2, *horizontal_shape)
    staggered_shape = (1, 3, *horizontal_shape)
    surface_shape = (1, *horizontal_shape)

    fields: dict[str, tuple[tuple[str, ...], np.ndarray]] = {
        "U": (("Time", "south_north", "west_east"), np.full(surface_shape, 10.0 + field_offset)),
        "V": (("Time", "south_north", "west_east"), np.full(surface_shape, 11.0 + field_offset)),
        "MU": (("Time", "south_north", "west_east"), np.full(surface_shape, 12.0 + field_offset)),
        "P": (("Time", "bottom_top", "south_north", "west_east"), np.full(mass_shape, 1000.0 + field_offset)),
        "PB": (("Time", "bottom_top", "south_north", "west_east"), np.full(mass_shape, 95000.0)),
        "T": (("Time", "bottom_top", "south_north", "west_east"), np.full(mass_shape, 5.0 + field_offset)),
        "QVAPOR": (("Time", "bottom_top", "south_north", "west_east"), np.full(mass_shape, 0.012)),
        "PH": (("Time", "bottom_top_stag", "south_north", "west_east"), np.zeros(staggered_shape)),
        "PHB": (("Time", "bottom_top_stag", "south_north", "west_east"), np.zeros(staggered_shape)),
    }
    for name, values in _tc_fields(horizontal_shape).items():
        fields[name] = (("Time", "south_north", "west_east"), values[None, :, :])

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", 1)
        dataset.createDimension("bottom_top", 2)
        dataset.createDimension("bottom_top_stag", 3)
        dataset.createDimension("south_north", horizontal_shape[0])
        dataset.createDimension("west_east", horizontal_shape[1])
        for name, (dimensions, values) in fields.items():
            if name in omit:
                continue
            variable = dataset.createVariable(name, "f8", dimensions)
            if name == "PSFC":
                variable.units = "Pa"
            variable[:] = values


def _write_pair(tmp_path: Path, *, omit_candidate: set[str] | None = None) -> tuple[Path, Path]:
    reference_dir = tmp_path / "reference"
    candidate_dir = tmp_path / "candidate"
    _write_wrfout(reference_dir / END_FILE)
    _write_wrfout(candidate_dir / END_FILE, omit=omit_candidate)
    return reference_dir, candidate_dir


def test_run_cycle_gate_with_slp_derives_then_passes_gate(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(tmp_path)
    derived_dir = tmp_path / "derived"

    report = run_cycle_gate_with_slp(
        reference_dir,
        candidate_dir,
        START,
        end=END,
        derived_dir=derived_dir,
    )
    payload = report_to_dict(report)

    assert report.status == "passed"
    assert payload["derived_reference_dir"] == str(derived_dir / "reference")
    assert payload["derived_candidate_dir"] == str(derived_dir / "candidate")
    assert len(payload["derivations"]) == 2
    assert {item["status"] for item in payload["derivations"]} == {"passed"}
    assert payload["gate_status"] == "passed"
    assert payload["gate_report"]["status"] == "passed"
    assert payload["gate_report"]["summary"] == {"total": 1, "passed": 1, "failed": 0}

    with netCDF4.Dataset(derived_dir / "reference" / END_FILE) as dataset:
        assert dataset.getncattr("TYWRF_DERIVED_SLP") == "true"
        assert dataset.variables["SLP"].units == "hPa"
        assert float(np.min(dataset.variables["SLP"][:])) == 950.0


def test_run_cycle_gate_with_slp_supports_multiple_cycle_endpoints(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(tmp_path)
    _write_wrfout(reference_dir / END_FILE2)
    _write_wrfout(candidate_dir / END_FILE2)

    report = run_cycle_gate_with_slp(
        reference_dir,
        candidate_dir,
        START,
        end=END2,
        derived_dir=tmp_path / "derived",
    )
    payload = report_to_dict(report)

    assert report.status == "passed"
    assert len(payload["derivations"]) == 4
    assert payload["gate_report"]["summary"] == {"total": 2, "passed": 2, "failed": 0}


def test_run_cycle_gate_with_slp_fails_derivation_without_psfc_fallback(tmp_path: Path) -> None:
    reference_dir, candidate_dir = _write_pair(tmp_path, omit_candidate={"QVAPOR"})
    derived_dir = tmp_path / "derived"

    report = run_cycle_gate_with_slp(
        reference_dir,
        candidate_dir,
        START,
        end=END,
        derived_dir=derived_dir,
    )
    payload = report_to_dict(report)
    candidate_derivation = next(
        item for item in payload["derivations"] if item["role"] == "candidate"
    )

    assert report.status == "failed"
    assert candidate_derivation["status"] == "failed"
    assert "QVAPOR" in candidate_derivation["message"]
    assert not (derived_dir / "candidate" / END_FILE).exists()
    assert payload["gate_status"] == "failed"
    assert payload["gate_report"]["cycles"][0]["candidate"] == str(
        derived_dir / "candidate" / END_FILE
    )


def test_run_cycle_gate_with_slp_cli_writes_combined_json(tmp_path: Path, capsys) -> None:
    reference_dir, candidate_dir = _write_pair(tmp_path)
    output = tmp_path / "report.json"

    exit_code = run_main(
        [
            "--reference-dir",
            str(reference_dir),
            "--candidate-dir",
            str(candidate_dir),
            "--start",
            START,
            "--end",
            END,
            "--derived-dir",
            str(tmp_path / "derived"),
            "--output",
            str(output),
            "--pretty",
        ]
    )

    assert exit_code == 0
    assert json.loads(capsys.readouterr().out)["status"] == "passed"
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["status"] == "passed"
    assert payload["gate_report"]["cycles"][0]["end_time"] == END
