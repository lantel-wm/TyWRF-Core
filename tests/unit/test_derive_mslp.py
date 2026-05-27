import json
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tools.derive_mslp import (
    GRAVITY_M_S2,
    KAPPA,
    RD_J_KG_K,
    DeriveMSLPError,
    derive_mslp_hpa,
    main as derive_main,
    write_derived_mslp,
)


def _write_core_fields(
    path: Path,
    *,
    omit: set[str] | None = None,
    qvapor_shape_mismatch: bool = False,
    include_existing_slp: bool = False,
) -> dict[str, np.ndarray]:
    omit = omit or set()
    shape = (1, 2, 3)
    mass_shape = (1, 2, 2, 3)
    staggered_shape = (1, 3, 2, 3)

    psfc = np.full(shape, 95000.0, dtype=np.float64)
    perturbation_pressure = np.full(mass_shape, 250.0, dtype=np.float64)
    base_pressure = np.full(mass_shape, 96000.0, dtype=np.float64)
    perturbation_potential_temperature = np.full(mass_shape, 5.0, dtype=np.float64)
    qvapor = np.full(mass_shape, 0.012, dtype=np.float64)
    surface_height_m = np.array(
        [[[0.0, 50.0, 100.0], [150.0, 200.0, 250.0]]],
        dtype=np.float64,
    )
    geopotential = np.zeros(staggered_shape, dtype=np.float64)
    geopotential[:, 0, :, :] = surface_height_m * GRAVITY_M_S2
    geopotential_base = np.zeros_like(geopotential)

    fields = {
        "PSFC": psfc,
        "P": perturbation_pressure,
        "PB": base_pressure,
        "T": perturbation_potential_temperature,
        "QVAPOR": qvapor,
        "PH": geopotential,
        "PHB": geopotential_base,
    }

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.setncattr("TITLE", "derive mslp unit test")
        dataset.createDimension("Time", None)
        dataset.createDimension("bottom_top", 2)
        dataset.createDimension("bottom_top_stag", 3)
        dataset.createDimension("south_north", 2)
        dataset.createDimension("west_east", 3)
        if qvapor_shape_mismatch:
            dataset.createDimension("west_east_qvapor", 2)

        for name, values in fields.items():
            if name in omit:
                continue
            dimensions = ("Time", "south_north", "west_east")
            if name in {"P", "PB", "T", "QVAPOR"}:
                dimensions = ("Time", "bottom_top", "south_north", "west_east")
            if name in {"PH", "PHB"}:
                dimensions = ("Time", "bottom_top_stag", "south_north", "west_east")
            if name == "QVAPOR" and qvapor_shape_mismatch:
                dimensions = ("Time", "bottom_top", "south_north", "west_east_qvapor")
                values = values[:, :, :, :2]
            variable = dataset.createVariable(name, "f8", dimensions)
            variable[:] = values

        if include_existing_slp:
            slp = dataset.createVariable("SLP", "f8", ("Time", "south_north", "west_east"))
            slp.units = "hPa"
            slp[:, :, :] = 9999.0

    return fields


def _expected_slp_hpa(fields: dict[str, np.ndarray]) -> np.ndarray:
    pressure_pa = fields["P"][:, 0, :, :] + fields["PB"][:, 0, :, :]
    theta_k = fields["T"][:, 0, :, :] + 300.0
    qvapor = fields["QVAPOR"][:, 0, :, :]
    surface_height_m = (fields["PH"][:, 0, :, :] + fields["PHB"][:, 0, :, :]) / GRAVITY_M_S2
    temperature_k = theta_k * np.power(pressure_pa / 100000.0, KAPPA)
    virtual_temperature_k = temperature_k * (1.0 + 0.61 * qvapor)
    return fields["PSFC"] * np.exp(
        GRAVITY_M_S2 * surface_height_m / (RD_J_KG_K * virtual_temperature_k)
    ) / 100.0


def test_derive_mslp_hpa_uses_hypsometric_virtual_temperature(tmp_path: Path) -> None:
    source = tmp_path / "wrfout.nc"
    fields = _write_core_fields(source)

    with netCDF4.Dataset(source) as dataset:
        derived = derive_mslp_hpa(dataset)

    expected = _expected_slp_hpa(fields)
    np.testing.assert_allclose(derived, expected, rtol=2e-7)
    assert derived.shape == fields["PSFC"].shape
    assert not np.allclose(derived, fields["PSFC"] / 100.0)


def test_write_derived_mslp_copies_source_and_replaces_existing_slp(tmp_path: Path) -> None:
    source = tmp_path / "source.nc"
    destination = tmp_path / "with_slp.nc"
    fields = _write_core_fields(source, include_existing_slp=True)

    summary = write_derived_mslp(source, destination)

    assert summary.variable == "SLP"
    assert summary.units == "hPa"
    assert summary.replaced_source_variable is True
    assert summary.shape == fields["PSFC"].shape
    assert summary.finite_count == fields["PSFC"].size
    assert summary.total_count == fields["PSFC"].size

    with netCDF4.Dataset(destination) as dataset:
        assert dataset.getncattr("TITLE") == "derive mslp unit test"
        assert dataset.getncattr("TYWRF_DERIVED_SLP") == "true"
        assert "P" in dataset.variables
        assert "SLP" in dataset.variables
        assert dataset.variables["SLP"].dimensions == dataset.variables["PSFC"].dimensions
        assert dataset.variables["SLP"].units == "hPa"
        assert dataset.variables["SLP"].derived_from == "P,PB,T,QVAPOR,PH,PHB,PSFC"
        np.testing.assert_allclose(
            dataset.variables["SLP"][:],
            _expected_slp_hpa(fields).astype(np.float32),
            rtol=2e-7,
        )
        assert not np.all(dataset.variables["SLP"][:] == 9999.0)


def test_derive_mslp_requires_all_core_source_fields(tmp_path: Path) -> None:
    source = tmp_path / "missing.nc"
    _write_core_fields(source, omit={"QVAPOR"})

    with netCDF4.Dataset(source) as dataset:
        with pytest.raises(DeriveMSLPError, match="missing required variable.*QVAPOR"):
            derive_mslp_hpa(dataset)


def test_derive_mslp_rejects_shape_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "mismatch.nc"
    _write_core_fields(source, qvapor_shape_mismatch=True)

    with netCDF4.Dataset(source) as dataset:
        with pytest.raises(DeriveMSLPError, match="source shapes do not match"):
            derive_mslp_hpa(dataset)


def test_derive_mslp_cli_writes_json_summary(tmp_path: Path, capsys) -> None:
    source = tmp_path / "source.nc"
    destination = tmp_path / "destination.nc"
    _write_core_fields(source)

    exit_code = derive_main([str(source), str(destination), "--pretty"])

    assert exit_code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["variable"] == "SLP"
    assert payload["finite_count"] == 6
    assert destination.exists()
