from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.write_core_wrfout import CORE_WRFOUT_VARIABLES, copy_core_wrfout


REFERENCE_WRFOUT_D01 = Path(
    "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/"
    "output_gfs_analysis/2025wp12/2025072600/WRF/"
    "wrfout_d01_2025-07-26_00:00:00"
)


def _make_source(path: Path) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.setncattr("TITLE", "unit test source")
        dataset.createDimension("Time", None)
        dataset.createDimension("DateStrLen", 19)
        dataset.createDimension("bottom_top", 2)
        dataset.createDimension("south_north", 3)
        dataset.createDimension("west_east", 4)
        dataset.createDimension("west_east_stag", 5)
        dataset.createDimension("soil_layers_stag", 4)

        times = dataset.createVariable("Times", "S1", ("Time", "DateStrLen"))
        times[0, :] = np.array(list("2025-07-26_00:00:00"), dtype="S1")

        u = dataset.createVariable("U", "f4", ("Time", "bottom_top", "south_north", "west_east_stag"))
        u.setncattr("MemoryOrder", "XYZ")
        u[:] = np.arange(1 * 2 * 3 * 5, dtype=np.float32).reshape(1, 2, 3, 5)

        psfc = dataset.createVariable("PSFC", "f4", ("Time", "south_north", "west_east"))
        psfc[:] = np.full((1, 3, 4), 100000.0, dtype=np.float32)


def test_copy_core_wrfout_copies_requested_variables_and_dimensions(tmp_path: Path) -> None:
    source = tmp_path / "source.nc"
    destination = tmp_path / "core.nc"
    _make_source(source)

    summary = copy_core_wrfout(source, destination, variables=("Times", "U", "PSFC"))

    assert summary.copied_variables == ["Times", "U", "PSFC"]
    assert summary.missing_variables == []
    assert summary.copied_dimensions == [
        "Time",
        "DateStrLen",
        "bottom_top",
        "south_north",
        "west_east",
        "west_east_stag",
        "soil_layers_stag",
    ]

    with netCDF4.Dataset(destination) as dataset:
        assert dataset.getncattr("TITLE") == "unit test source"
        assert dataset.getncattr("TYWRF_CORE_WRFOUT") == "true"
        assert dataset.dimensions["Time"].isunlimited()
        assert len(dataset.dimensions["soil_layers_stag"]) == 4
        assert set(dataset.variables) == {"Times", "U", "PSFC"}
        assert dataset.variables["U"].getncattr("MemoryOrder") == "XYZ"
        np.testing.assert_array_equal(dataset.variables["PSFC"][:], np.full((1, 3, 4), 100000.0))


def test_copy_core_wrfout_requires_variables_by_default(tmp_path: Path) -> None:
    source = tmp_path / "source.nc"
    _make_source(source)

    with pytest.raises(KeyError, match="missing required variables"):
        copy_core_wrfout(source, tmp_path / "core.nc", variables=("Times", "V"))


def test_copy_core_wrfout_can_allow_missing_variables(tmp_path: Path) -> None:
    source = tmp_path / "source.nc"
    destination = tmp_path / "core.nc"
    _make_source(source)

    summary = copy_core_wrfout(source, destination, variables=("Times", "V"), allow_missing=True)

    assert summary.copied_variables == ["Times"]
    assert summary.missing_variables == ["V"]


@pytest.mark.skipif(
    not REFERENCE_WRFOUT_D01.exists(),
    reason="KROSA reference wrfout_d01 is not present",
)
def test_reference_wrfout_has_core_variables_and_subset_copy_preserves_schema(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "core_d01.nc"

    with netCDF4.Dataset(REFERENCE_WRFOUT_D01) as source:
        missing = [name for name in CORE_WRFOUT_VARIABLES if name not in source.variables]
        source_dimensions = list(source.dimensions)
        source_data_model = source.data_model
        source_dx = source.getncattr("DX")

    assert missing == []

    summary = copy_core_wrfout(
        REFERENCE_WRFOUT_D01,
        destination,
        variables=("Times", "XLAT", "RAINNC"),
    )

    assert summary.copied_variables == ["Times", "XLAT", "RAINNC"]
    assert summary.copied_dimensions == source_dimensions

    with netCDF4.Dataset(destination) as dataset:
        assert dataset.data_model == source_data_model
        assert list(dataset.dimensions) == source_dimensions
        assert dataset.dimensions["Time"].isunlimited()
        assert dataset.getncattr("DX") == source_dx
        assert dataset.getncattr("TYWRF_CORE_WRFOUT") == "true"
        assert dataset.variables["XLAT"].dimensions == ("Time", "south_north", "west_east")
        assert dataset.variables["RAINNC"].dimensions == ("Time", "south_north", "west_east")
