from pathlib import Path

import netCDF4
import numpy as np

from tools.extract_reference import main as extract_main
from tools.extract_reference import summarize


def _make_reference(path: Path) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        dataset.createDimension("DateStrLen", 19)
        dataset.createDimension("bottom_top", 2)
        dataset.createDimension("south_north", 3)
        dataset.createDimension("west_east_stag", 5)

        times = dataset.createVariable("Times", "S1", ("Time", "DateStrLen"))
        times[0, :] = np.array(list("2025-07-26_00:00:00"), dtype="S1")

        u = dataset.createVariable("U", "f4", ("Time", "bottom_top", "south_north", "west_east_stag"))
        u[:] = np.ones((1, 2, 3, 5), dtype=np.float32)


def test_summarize_reports_dimensions_present_missing_and_variable_metadata(tmp_path: Path) -> None:
    path = tmp_path / "wrfout_d01_2025-07-26_00:00:00"
    _make_reference(path)

    summary = summarize(path, core_variables=("Times", "U", "V"))

    assert summary.status == "missing_core_variables"
    assert summary.dimensions["bottom_top"] == 2
    assert summary.present_core_variables == ["Times", "U"]
    assert summary.missing_core_variables == ["V"]
    assert summary.variables["U"].dimensions == ("Time", "bottom_top", "south_north", "west_east_stag")
    assert summary.variables["U"].shape == (1, 2, 3, 5)
    assert summary.variables["U"].dtype == "float32"


def test_main_writes_json_report(tmp_path: Path) -> None:
    path = tmp_path / "wrfout_d01_2025-07-26_00:00:00"
    output = tmp_path / "summary.json"
    _make_reference(path)

    exit_code = extract_main(
        [str(path), "--core-variables", "Times", "U", "--output", str(output), "--pretty"]
    )

    assert exit_code == 0
    assert '"status": "ok"' in output.read_text(encoding="utf-8")
