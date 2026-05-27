#!/usr/bin/env python
"""Compute baseline tropical-cyclone diagnostics from WRF-compatible output."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import netCDF4
import numpy as np


TC_CENTER_THRESHOLD_KM = 20.0
MINIMUM_SLP_THRESHOLD_HPA = 5.0
MSLP_PROXY_THRESHOLD_HPA = MINIMUM_SLP_THRESHOLD_HPA
VMAX_THRESHOLD_MS = 5.0
MSLP_PROXY_LABEL = "PSFC-min proxy; not WRF sea-level pressure"
MINIMUM_SLP_STATUS_OK = "ok"
MINIMUM_SLP_STATUS_NOT_AVAILABLE = "not_available"
MINIMUM_SLP_VARIABLE_CANDIDATES = (
    "SLP",
    "slp",
    "MSLP",
    "mslp",
    "AFWA_MSLP",
    "afwa_mslp",
    "PMSL",
    "pmsl",
    "PRMSL",
    "prmsl",
    "SEA_LEVEL_PRESSURE",
    "sea_level_pressure",
)


class TCDiagnosticsError(ValueError):
    """Raised when a NetCDF file cannot provide the required TC diagnostics."""


@dataclass(frozen=True)
class GridPointDiagnostic:
    j: int
    i: int
    latitude: float
    longitude: float


@dataclass(frozen=True)
class RainfallSummary:
    minimum_mm: float
    maximum_mm: float
    mean_mm: float
    total_grid_mm: float
    valid_count: int
    total_count: int


@dataclass(frozen=True)
class TCDiagnostics:
    path: str
    status: str
    time_index: int
    center: GridPointDiagnostic
    center_source: str
    minimum_slp_hpa: float | None
    minimum_slp_status: str
    minimum_slp_source: str | None
    minimum_slp_source_units: str | None
    minimum_slp_reason: str | None
    minimum_slp_location: GridPointDiagnostic | None
    psfc_min_pa: float
    psfc_min_location: GridPointDiagnostic
    mslp_proxy_hpa: float
    mslp_proxy_label: str
    vmax_10m_ms: float
    vmax_location: GridPointDiagnostic
    rainfall: RainfallSummary


@dataclass(frozen=True)
class TCDiagnosticsComparison:
    status: str
    reference: TCDiagnostics | None
    candidate: TCDiagnostics | None
    thresholds: dict[str, float]
    center_error_km: float | None = None
    minimum_slp_status: str | None = None
    minimum_slp_error_hpa: float | None = None
    minimum_slp_abs_error_hpa: float | None = None
    minimum_slp_message: str | None = None
    mslp_proxy_error_hpa: float | None = None
    mslp_proxy_abs_error_hpa: float | None = None
    vmax_error_ms: float | None = None
    vmax_abs_error_ms: float | None = None
    rainfall_mean_error_mm: float | None = None
    rainfall_max_error_mm: float | None = None
    rainfall_total_grid_error_mm: float | None = None
    message: str | None = None


@dataclass(frozen=True)
class _MinimumSLPDiagnostic:
    status: str
    hpa: float | None
    source: str | None
    source_units: str | None
    reason: str | None
    location: GridPointDiagnostic | None


def _read_2d(dataset: netCDF4.Dataset, variable: str, time_index: int) -> np.ndarray:
    if variable not in dataset.variables:
        raise TCDiagnosticsError(f"missing required variable {variable}")

    try:
        data = np.ma.asarray(dataset.variables[variable][:], dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise TCDiagnosticsError(f"variable {variable} is not numeric: {exc}") from exc

    try:
        if data.ndim == 2:
            selected = data
        elif data.ndim == 3:
            selected = data[time_index, :, :]
        else:
            raise TCDiagnosticsError(
                f"variable {variable} must be 2D or Time+2D, got shape {data.shape}"
            )
    except IndexError as exc:
        raise TCDiagnosticsError(
            f"time index {time_index} is out of range for variable {variable} with shape {data.shape}"
        ) from exc

    return np.asarray(np.ma.filled(selected, np.nan), dtype=np.float64)


def _require_same_shape(arrays: dict[str, np.ndarray]) -> tuple[int, int]:
    shapes = {name: values.shape for name, values in arrays.items()}
    unique_shapes = set(shapes.values())
    if len(unique_shapes) != 1:
        detail = ", ".join(f"{name}={shape}" for name, shape in sorted(shapes.items()))
        raise TCDiagnosticsError(f"diagnostic variable shapes do not match: {detail}")

    shape = next(iter(unique_shapes))
    if len(shape) != 2:
        raise TCDiagnosticsError(f"diagnostic variables must be 2D after time selection, got {shape}")
    return int(shape[0]), int(shape[1])


def _grid_point(j: int, i: int, latitude: np.ndarray, longitude: np.ndarray) -> GridPointDiagnostic:
    return GridPointDiagnostic(
        j=int(j),
        i=int(i),
        latitude=float(latitude[j, i]),
        longitude=float(longitude[j, i]),
    )


def _finite_pair_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(arrays[0].shape, dtype=bool)
    for values in arrays:
        mask &= np.isfinite(values)
    return mask


def _argmin_with_mask(values: np.ndarray, valid_mask: np.ndarray, variable: str) -> tuple[int, int]:
    if not np.any(valid_mask):
        raise TCDiagnosticsError(f"no finite samples available for {variable}")
    masked = np.where(valid_mask, values, np.inf)
    return tuple(int(index) for index in np.unravel_index(np.argmin(masked), values.shape))


def _argmax_with_mask(values: np.ndarray, valid_mask: np.ndarray, variable: str) -> tuple[int, int]:
    if not np.any(valid_mask):
        raise TCDiagnosticsError(f"no finite samples available for {variable}")
    masked = np.where(valid_mask, values, -np.inf)
    return tuple(int(index) for index in np.unravel_index(np.argmax(masked), values.shape))


def _rainfall_summary(rainc: np.ndarray, rainnc: np.ndarray) -> RainfallSummary:
    _require_same_shape({"RAINC": rainc, "RAINNC": rainnc})
    total = rainc + rainnc
    valid_mask = _finite_pair_mask(rainc, rainnc, total)
    valid_values = total[valid_mask]
    if valid_values.size == 0:
        raise TCDiagnosticsError("no finite samples available for accumulated rainfall")

    return RainfallSummary(
        minimum_mm=float(np.min(valid_values)),
        maximum_mm=float(np.max(valid_values)),
        mean_mm=float(np.mean(valid_values)),
        total_grid_mm=float(np.sum(valid_values)),
        valid_count=int(valid_values.size),
        total_count=int(total.size),
    )


def _find_minimum_slp_variable(dataset: netCDF4.Dataset) -> str | None:
    for candidate in MINIMUM_SLP_VARIABLE_CANDIDATES:
        if candidate in dataset.variables:
            return candidate

    variables_by_lower = {name.lower(): name for name in dataset.variables}
    for candidate in MINIMUM_SLP_VARIABLE_CANDIDATES:
        match = variables_by_lower.get(candidate.lower())
        if match is not None:
            return match
    return None


def _variable_units(dataset: netCDF4.Dataset, variable: str) -> str | None:
    units = getattr(dataset.variables[variable], "units", None)
    if units is None:
        return None
    return str(units)


def _pressure_to_hpa(values: np.ndarray, units: str | None) -> np.ndarray:
    normalized_units = (units or "").strip().lower().replace(" ", "")
    if normalized_units in {"pa", "pascal", "pascals"}:
        return values / 100.0
    if normalized_units in {
        "hpa",
        "hectopascal",
        "hectopascals",
        "mb",
        "mbar",
        "millibar",
        "millibars",
    }:
        return values

    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return values
    if float(np.nanmedian(np.abs(finite_values))) > 2000.0:
        return values / 100.0
    return values


def _minimum_slp_not_available(reason: str) -> _MinimumSLPDiagnostic:
    return _MinimumSLPDiagnostic(
        status=MINIMUM_SLP_STATUS_NOT_AVAILABLE,
        hpa=None,
        source=None,
        source_units=None,
        reason=reason,
        location=None,
    )


def _minimum_slp_diagnostic(
    dataset: netCDF4.Dataset,
    *,
    time_index: int,
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> _MinimumSLPDiagnostic:
    variable = _find_minimum_slp_variable(dataset)
    if variable is None:
        checked = ", ".join(MINIMUM_SLP_VARIABLE_CANDIDATES)
        return _minimum_slp_not_available(
            f"no sea-level pressure variable found; checked {checked}"
        )

    units = _variable_units(dataset, variable)
    try:
        values = _pressure_to_hpa(_read_2d(dataset, variable, time_index), units)
        _require_same_shape({variable: values, "XLAT": latitude, "XLONG": longitude})
        valid_mask = _finite_pair_mask(values, latitude, longitude)
        minimum_j, minimum_i = _argmin_with_mask(values, valid_mask, variable)
    except TCDiagnosticsError as exc:
        return _minimum_slp_not_available(
            f"sea-level pressure variable {variable} is present but unusable: {exc}"
        )

    return _MinimumSLPDiagnostic(
        status=MINIMUM_SLP_STATUS_OK,
        hpa=float(values[minimum_j, minimum_i]),
        source=variable,
        source_units=units,
        reason=None,
        location=_grid_point(minimum_j, minimum_i, latitude, longitude),
    )


def compute_tc_diagnostics(
    dataset: netCDF4.Dataset,
    *,
    path: str = "",
    time_index: int = -1,
) -> TCDiagnostics:
    """Return baseline TC diagnostics for a WRF-compatible domain output file."""

    psfc = _read_2d(dataset, "PSFC", time_index)
    latitude = _read_2d(dataset, "XLAT", time_index)
    longitude = _read_2d(dataset, "XLONG", time_index)
    u10 = _read_2d(dataset, "U10", time_index)
    v10 = _read_2d(dataset, "V10", time_index)
    rainc = _read_2d(dataset, "RAINC", time_index)
    rainnc = _read_2d(dataset, "RAINNC", time_index)

    _require_same_shape({"PSFC": psfc, "XLAT": latitude, "XLONG": longitude})
    _require_same_shape({"U10": u10, "V10": v10, "XLAT": latitude, "XLONG": longitude})

    center_mask = _finite_pair_mask(psfc, latitude, longitude)
    center_j, center_i = _argmin_with_mask(psfc, center_mask, "PSFC")
    psfc_min_location = _grid_point(center_j, center_i, latitude, longitude)
    minimum_slp = _minimum_slp_diagnostic(
        dataset,
        time_index=time_index,
        latitude=latitude,
        longitude=longitude,
    )
    center = minimum_slp.location if minimum_slp.location is not None else psfc_min_location
    center_source = minimum_slp.source if minimum_slp.source is not None else "PSFC"

    speed = np.sqrt(u10 * u10 + v10 * v10)
    wind_mask = _finite_pair_mask(speed, latitude, longitude)
    vmax_j, vmax_i = _argmax_with_mask(speed, wind_mask, "U10/V10 magnitude")

    psfc_min_pa = float(psfc[center_j, center_i])
    return TCDiagnostics(
        path=path,
        status="ok",
        time_index=int(time_index),
        center=center,
        center_source=center_source,
        minimum_slp_hpa=minimum_slp.hpa,
        minimum_slp_status=minimum_slp.status,
        minimum_slp_source=minimum_slp.source,
        minimum_slp_source_units=minimum_slp.source_units,
        minimum_slp_reason=minimum_slp.reason,
        minimum_slp_location=minimum_slp.location,
        psfc_min_pa=psfc_min_pa,
        psfc_min_location=psfc_min_location,
        mslp_proxy_hpa=psfc_min_pa / 100.0,
        mslp_proxy_label=MSLP_PROXY_LABEL,
        vmax_10m_ms=float(speed[vmax_j, vmax_i]),
        vmax_location=_grid_point(vmax_j, vmax_i, latitude, longitude),
        rainfall=_rainfall_summary(rainc, rainnc),
    )


def diagnose_file(path: Path, *, time_index: int = -1) -> TCDiagnostics:
    with netCDF4.Dataset(path) as dataset:
        return compute_tc_diagnostics(dataset, path=str(path), time_index=time_index)


def haversine_km(
    latitude_a: float,
    longitude_a: float,
    latitude_b: float,
    longitude_b: float,
) -> float:
    radius_km = 6371.0
    lat_a = math.radians(latitude_a)
    lat_b = math.radians(latitude_b)
    delta_lat = math.radians(latitude_b - latitude_a)
    delta_lon = math.radians(longitude_b - longitude_a)
    sin_lat = math.sin(delta_lat / 2.0)
    sin_lon = math.sin(delta_lon / 2.0)
    value = sin_lat * sin_lat + math.cos(lat_a) * math.cos(lat_b) * sin_lon * sin_lon
    return float(2.0 * radius_km * math.asin(min(1.0, math.sqrt(value))))


def compare_tc_diagnostics(
    reference_path: Path,
    candidate_path: Path,
    *,
    time_index: int = -1,
) -> TCDiagnosticsComparison:
    thresholds = {
        "center_error_km": TC_CENTER_THRESHOLD_KM,
        "minimum_slp_abs_error_hpa": MINIMUM_SLP_THRESHOLD_HPA,
        "vmax_abs_error_ms": VMAX_THRESHOLD_MS,
    }

    try:
        reference = diagnose_file(reference_path, time_index=time_index)
    except TCDiagnosticsError as exc:
        return TCDiagnosticsComparison(
            status="failed",
            reference=None,
            candidate=None,
            thresholds=thresholds,
            message=f"reference diagnostics failed: {exc}",
        )

    try:
        candidate = diagnose_file(candidate_path, time_index=time_index)
    except TCDiagnosticsError as exc:
        return TCDiagnosticsComparison(
            status="failed",
            reference=reference,
            candidate=None,
            thresholds=thresholds,
            message=f"candidate diagnostics failed: {exc}",
        )

    center_error_km = haversine_km(
        reference.center.latitude,
        reference.center.longitude,
        candidate.center.latitude,
        candidate.center.longitude,
    )
    minimum_slp_status = MINIMUM_SLP_STATUS_NOT_AVAILABLE
    minimum_slp_error_hpa = None
    minimum_slp_abs_error_hpa = None
    minimum_slp_message = None
    if (
        reference.minimum_slp_status == MINIMUM_SLP_STATUS_OK
        and candidate.minimum_slp_status == MINIMUM_SLP_STATUS_OK
        and reference.minimum_slp_hpa is not None
        and candidate.minimum_slp_hpa is not None
    ):
        minimum_slp_status = MINIMUM_SLP_STATUS_OK
        minimum_slp_error_hpa = candidate.minimum_slp_hpa - reference.minimum_slp_hpa
        minimum_slp_abs_error_hpa = abs(minimum_slp_error_hpa)
    else:
        reasons = []
        if reference.minimum_slp_status != MINIMUM_SLP_STATUS_OK:
            reasons.append(f"reference: {reference.minimum_slp_reason}")
        if candidate.minimum_slp_status != MINIMUM_SLP_STATUS_OK:
            reasons.append(f"candidate: {candidate.minimum_slp_reason}")
        minimum_slp_message = "minimum SLP comparison not_available; " + "; ".join(reasons)

    mslp_proxy_error_hpa = candidate.mslp_proxy_hpa - reference.mslp_proxy_hpa
    vmax_error_ms = candidate.vmax_10m_ms - reference.vmax_10m_ms
    rainfall_mean_error_mm = candidate.rainfall.mean_mm - reference.rainfall.mean_mm
    rainfall_max_error_mm = candidate.rainfall.maximum_mm - reference.rainfall.maximum_mm
    rainfall_total_grid_error_mm = (
        candidate.rainfall.total_grid_mm - reference.rainfall.total_grid_mm
    )

    exceeded = []
    if center_error_km > thresholds["center_error_km"]:
        exceeded.append("center_error_km")
    minimum_slp_unavailable = minimum_slp_status != MINIMUM_SLP_STATUS_OK
    if (
        minimum_slp_abs_error_hpa is not None
        and minimum_slp_abs_error_hpa > thresholds["minimum_slp_abs_error_hpa"]
    ):
        exceeded.append("minimum_slp_abs_error_hpa")
    if abs(vmax_error_ms) > thresholds["vmax_abs_error_ms"]:
        exceeded.append("vmax_abs_error_ms")

    if exceeded:
        status = "threshold_exceeded"
    elif minimum_slp_unavailable:
        status = "failed"
    else:
        status = "ok"
    messages = []
    if exceeded:
        messages.append("TC diagnostic thresholds exceeded: " + ", ".join(exceeded))
    if minimum_slp_message is not None:
        messages.append(minimum_slp_message)
    message = " ".join(messages) if messages else None

    return TCDiagnosticsComparison(
        status=status,
        reference=reference,
        candidate=candidate,
        thresholds=thresholds,
        center_error_km=center_error_km,
        minimum_slp_status=minimum_slp_status,
        minimum_slp_error_hpa=minimum_slp_error_hpa,
        minimum_slp_abs_error_hpa=minimum_slp_abs_error_hpa,
        minimum_slp_message=minimum_slp_message,
        mslp_proxy_error_hpa=mslp_proxy_error_hpa,
        mslp_proxy_abs_error_hpa=abs(mslp_proxy_error_hpa),
        vmax_error_ms=vmax_error_ms,
        vmax_abs_error_ms=abs(vmax_error_ms),
        rainfall_mean_error_mm=rainfall_mean_error_mm,
        rainfall_max_error_mm=rainfall_max_error_mm,
        rainfall_total_grid_error_mm=rainfall_total_grid_error_mm,
        message=message,
    )
