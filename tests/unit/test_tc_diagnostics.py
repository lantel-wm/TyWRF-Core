import math
from pathlib import Path

import netCDF4
import numpy as np

from tools.tc_diagnostics import compare_tc_diagnostics, diagnose_file


def _write_tc_dataset(
    path: Path,
    *,
    psfc: np.ndarray,
    u10: np.ndarray,
    v10: np.ndarray,
    rainc: np.ndarray,
    rainnc: np.ndarray,
    latitude: np.ndarray | None = None,
    longitude: np.ndarray | None = None,
) -> None:
    y_size, x_size = psfc.shape
    if latitude is None:
        latitude = np.array(
            [
                [10.0, 10.0, 10.0, 10.0],
                [11.0, 11.0, 11.0, 11.0],
                [12.0, 12.0, 12.0, 12.0],
            ],
            dtype=np.float64,
        )
    if longitude is None:
        longitude = np.array(
            [
                [120.0, 121.0, 122.0, 123.0],
                [120.0, 121.0, 122.0, 123.0],
                [120.0, 121.0, 122.0, 123.0],
            ],
            dtype=np.float64,
        )

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", 1)
        dataset.createDimension("south_north", y_size)
        dataset.createDimension("west_east", x_size)

        for name, values in {
            "XLAT": latitude,
            "XLONG": longitude,
            "PSFC": psfc,
            "U10": u10,
            "V10": v10,
            "RAINC": rainc,
            "RAINNC": rainnc,
        }.items():
            variable = dataset.createVariable(name, "f8", ("Time", "south_north", "west_east"))
            variable[0, :, :] = np.asarray(values, dtype=np.float64)


def test_diagnose_file_reports_center_psfc_proxy_vmax_and_rainfall(tmp_path: Path) -> None:
    path = tmp_path / "wrfout_d02_2025-07-26_06:00:00"
    psfc = np.array(
        [
            [100000.0, 99500.0, 99000.0, 99800.0],
            [98500.0, 97000.0, 95000.0, 98000.0],
            [99900.0, 99600.0, 99400.0, 99700.0],
        ],
        dtype=np.float64,
    )
    u10 = np.zeros_like(psfc)
    v10 = np.zeros_like(psfc)
    u10[0, 1] = 30.0
    v10[0, 1] = 40.0
    rainc = np.array(
        [
            [0.0, 1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0, 7.0],
            [8.0, 9.0, 10.0, 11.0],
        ],
        dtype=np.float64,
    )
    rainnc = np.full_like(rainc, 0.5)
    _write_tc_dataset(path, psfc=psfc, u10=u10, v10=v10, rainc=rainc, rainnc=rainnc)

    diagnostics = diagnose_file(path)

    assert diagnostics.status == "ok"
    assert diagnostics.center.j == 1
    assert diagnostics.center.i == 2
    assert diagnostics.center.latitude == 11.0
    assert diagnostics.center.longitude == 122.0
    assert diagnostics.psfc_min_pa == 95000.0
    assert diagnostics.mslp_proxy_hpa == 950.0
    assert "PSFC-min proxy" in diagnostics.mslp_proxy_label
    assert diagnostics.vmax_10m_ms == 50.0
    assert diagnostics.vmax_location.j == 0
    assert diagnostics.vmax_location.i == 1
    assert diagnostics.rainfall.minimum_mm == 0.5
    assert diagnostics.rainfall.maximum_mm == 11.5
    assert diagnostics.rainfall.mean_mm == 6.0
    assert diagnostics.rainfall.total_grid_mm == 72.0
    assert diagnostics.rainfall.valid_count == 12
    assert diagnostics.rainfall.total_count == 12


def test_compare_tc_diagnostics_reports_domain_errors(tmp_path: Path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    shape = (3, 4)
    psfc_reference = np.full(shape, 100000.0)
    psfc_reference[1, 1] = 95000.0
    psfc_candidate = psfc_reference.copy()
    psfc_candidate[1, 1] = 95200.0
    u10_reference = np.zeros(shape)
    v10_reference = np.zeros(shape)
    u10_candidate = np.zeros(shape)
    v10_candidate = np.zeros(shape)
    u10_reference[0, 0] = 10.0
    v10_reference[0, 0] = 10.0
    u10_candidate[0, 0] = 13.0
    v10_candidate[0, 0] = 14.0
    rainc_reference = np.ones(shape)
    rainnc_reference = np.ones(shape)
    rainc_candidate = np.full(shape, 2.0)
    rainnc_candidate = np.ones(shape)
    _write_tc_dataset(
        reference,
        psfc=psfc_reference,
        u10=u10_reference,
        v10=v10_reference,
        rainc=rainc_reference,
        rainnc=rainnc_reference,
    )
    _write_tc_dataset(
        candidate,
        psfc=psfc_candidate,
        u10=u10_candidate,
        v10=v10_candidate,
        rainc=rainc_candidate,
        rainnc=rainnc_candidate,
    )

    comparison = compare_tc_diagnostics(reference, candidate)

    assert comparison.status == "ok"
    assert comparison.center_error_km == 0.0
    assert comparison.mslp_proxy_error_hpa == 2.0
    assert comparison.mslp_proxy_abs_error_hpa == 2.0
    assert math.isclose(comparison.vmax_error_ms, math.sqrt(365.0) - math.sqrt(200.0))
    assert comparison.rainfall_mean_error_mm == 1.0
    assert comparison.rainfall_max_error_mm == 1.0
    assert comparison.rainfall_total_grid_error_mm == 12.0
