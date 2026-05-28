#!/usr/bin/env python
"""Run d02 6 h cycle acceptance gates against WRF reference output."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timedelta
import json
import math
from pathlib import Path

import netCDF4
import numpy as np

try:
    from tools.compare_wrfout import (
        STRICT_CORE_VARIABLES,
        VariableComparison,
        compare_files,
    )
    from tools.tc_diagnostics import MINIMUM_SLP_VARIABLE_CANDIDATES
except ModuleNotFoundError:
    from compare_wrfout import STRICT_CORE_VARIABLES, VariableComparison, compare_files
    from tc_diagnostics import MINIMUM_SLP_VARIABLE_CANDIDATES


WRF_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"
DEFAULT_DOMAIN = "d02"
DEFAULT_INTERVAL_HOURS = 6
DEFAULT_INTERVAL_MINUTES = DEFAULT_INTERVAL_HOURS * 60
GATE_FIELD_THRESHOLD = 0.05
GATE_FIELD_THRESHOLDS = {name: GATE_FIELD_THRESHOLD for name in STRICT_CORE_VARIABLES}
METADATA_GATE_THRESHOLD = 0.0
BOOL_TRUE_VALUES = {"1", "true", "t", "yes", "y", "on"}
BOOL_FALSE_VALUES = {"0", "false", "f", "no", "n", "off"}
BLOCKING_CANDIDATE_KIND_TOKENS = (
    "diagnostic",
    "closure",
    "remap",
    "oracle",
    "helper",
    "probe",
    "adapter",
    "dry_run",
    "staging",
    "experimental",
)
BLOCKING_RELATED_METADATA_TOKENS = (
    "helper",
    "probe",
    "adapter",
    "dry_run",
    "staging",
    "experimental",
)
WIND_TENDENCY_SOURCE_KIND_BLOCKING_TOKENS = (
    "zero",
    "identity",
    "oracle",
    "reference_end",
    "reference_end_truth",
    "reference_truth",
    "end_truth",
    "placeholder",
    "non_evidence",
    "nonevidence",
)
WIND_TENDENCY_FALSE_EVIDENCE_ATTRS = (
    "TYWRF_WIND_TENDENCY_GATE_EVIDENCE",
    "TYWRF_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE",
)
WIND_TENDENCY_TRUE_BLOCKER_ATTRS = (
    "TYWRF_WIND_TENDENCY_USES_REFERENCE_END_TRUTH",
    "TYWRF_WIND_TENDENCY_ZERO_OR_IDENTITY_ONLY",
)
PRESSURE_GRADIENT_WIND_TENDENCY_PREFIX = (
    "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY"
)
PRESSURE_GRADIENT_WIND_TENDENCY_REQUIRED_GROUPS = {
    "switches": (
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_OPT_IN",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_APPLIED",
    ),
    "mode": (
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_MODE",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_SOURCE_KIND",
    ),
    "status": ("TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_STATUS",),
    "evidence": (
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_GATE_EVIDENCE",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE",
    ),
    "counts": (
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_U_UPDATED_COUNT",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_V_UPDATED_COUNT",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_U_SKIPPED_COUNT",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_V_SKIPPED_COUNT",
    ),
    "units": (
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_UNITS",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_UNITS",
    ),
    "alpha": (
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_KIND",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_VALUE",
    ),
    "safety": (
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_DIAGNOSTIC_ONLY",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_USES_ORACLE",
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_USES_REFERENCE_END_TRUTH",
    ),
}
PRESSURE_GRADIENT_WIND_TENDENCY_FALSE_EVIDENCE_ATTRS = (
    "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_GATE_EVIDENCE",
    "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE",
)
PRESSURE_GRADIENT_WIND_TENDENCY_TRUE_BLOCKER_ATTRS = (
    "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_DIAGNOSTIC_ONLY",
    "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_USES_ORACLE",
    "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_USES_REFERENCE_END_TRUTH",
)


@dataclass(frozen=True)
class GateMetric:
    name: str
    status: str
    threshold: float
    value: float | None = None
    units: str | None = None
    message: str | None = None


@dataclass(frozen=True)
class GateField:
    variable: str
    status: str
    threshold: float
    rmse: float | None = None
    normalized_rmse: float | None = None
    max_abs_error: float | None = None
    valid_count: int | None = None
    total_count: int | None = None
    source_status: str | None = None
    message: str | None = None


@dataclass(frozen=True)
class CycleGate:
    status: str
    domain: str
    start_time: str
    end_time: str
    reference: str
    candidate: str
    fields: list[GateField]
    diagnostics: list[GateMetric]
    message: str | None = None


@dataclass(frozen=True)
class CycleGateReport:
    status: str
    domain: str
    reference_dir: str
    candidate_dir: str
    start_time: str
    end_time: str
    interval_hours: int
    interval_minutes: int
    cycles: list[CycleGate]
    summary: dict[str, int]
    first_failure: dict[str, object | None] | None


def parse_wrf_time(value: str) -> datetime:
    for fmt in (WRF_TIME_FORMAT, "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            pass
    raise ValueError(f"time must match YYYY-MM-DD_HH:MM:SS: {value}")


def format_wrf_time(value: datetime) -> str:
    return value.strftime(WRF_TIME_FORMAT)


def wrfout_filename(domain: str, valid_time: datetime) -> str:
    return f"wrfout_{domain}_{format_wrf_time(valid_time)}"


def _resolve_interval_minutes(
    *,
    interval_hours: int | None = DEFAULT_INTERVAL_HOURS,
    interval_minutes: int | None = None,
) -> int:
    if interval_minutes is not None:
        if interval_minutes <= 0:
            raise ValueError("interval minutes must be positive")
        return interval_minutes
    if interval_hours is None:
        raise ValueError("interval is required")
    if interval_hours <= 0:
        raise ValueError("interval must be positive")
    return interval_hours * 60


def cycle_end_times(
    start_time: datetime,
    end_time: datetime,
    interval_hours: int | None = DEFAULT_INTERVAL_HOURS,
    *,
    interval_minutes: int | None = None,
) -> list[datetime]:
    resolved_interval_minutes = _resolve_interval_minutes(
        interval_hours=interval_hours,
        interval_minutes=interval_minutes,
    )
    if resolved_interval_minutes <= 0:
        raise ValueError("interval must be positive")
    if end_time <= start_time:
        raise ValueError("end time must be after start time")

    interval = timedelta(minutes=resolved_interval_minutes)
    current = start_time + interval
    ends: list[datetime] = []
    while current <= end_time:
        ends.append(current)
        current += interval

    if not ends or ends[-1] != end_time:
        raise ValueError("end time must align with start + N * interval")
    return ends


def resolve_end_time(
    start_time: datetime,
    *,
    end: str | datetime | None = None,
    hours: int | None = None,
) -> datetime:
    if end is not None and hours is not None:
        raise ValueError("use either end or hours, not both")
    if end is not None:
        return parse_wrf_time(end) if isinstance(end, str) else end
    duration = DEFAULT_INTERVAL_HOURS if hours is None else hours
    if duration <= 0:
        raise ValueError("hours must be positive")
    return start_time + timedelta(hours=duration)


def _missing_fields(message: str) -> list[GateField]:
    return [
        GateField(
            variable=variable,
            status="not_available",
            threshold=threshold,
            source_status="not_available",
            message=message,
        )
        for variable, threshold in GATE_FIELD_THRESHOLDS.items()
    ]


def _missing_diagnostics(message: str) -> list[GateMetric]:
    return [
        GateMetric(
            name="storm_center",
            status="not_available",
            threshold=20.0,
            units="km",
            message=message,
        ),
        GateMetric(
            name="minimum_slp",
            status="not_available",
            threshold=5.0,
            units="hPa",
            message=message,
        ),
        GateMetric(
            name="vmax10m",
            status="not_available",
            threshold=5.0,
            units="m s-1",
            message=message,
        ),
    ]


def _field_status(comparison: VariableComparison) -> str:
    if comparison.status == "ok":
        return "passed"
    if comparison.status == "threshold_exceeded":
        return "failed"
    return "not_available"


def _gate_field(comparison: VariableComparison) -> GateField:
    return GateField(
        variable=comparison.variable,
        status=_field_status(comparison),
        threshold=GATE_FIELD_THRESHOLDS[comparison.variable],
        rmse=comparison.rmse,
        normalized_rmse=comparison.normalized_rmse,
        max_abs_error=comparison.max_abs_error,
        valid_count=comparison.valid_count,
        total_count=comparison.total_count,
        source_status=comparison.status,
        message=comparison.message,
    )


def _metric_status(value: float | None, threshold: float) -> str:
    if value is None or not math.isfinite(value):
        return "not_available"
    return "passed" if value <= threshold else "failed"


def _read_bool_attr(dataset: netCDF4.Dataset, name: str) -> bool | None:
    if name not in dataset.ncattrs():
        return None
    raw_value = dataset.getncattr(name)
    if isinstance(raw_value, (bool, np.bool_)):
        return bool(raw_value)
    text = str(raw_value).strip().lower()
    if text in BOOL_TRUE_VALUES:
        return True
    if text in BOOL_FALSE_VALUES:
        return False
    return None


def _read_text_attr(dataset: netCDF4.Dataset, name: str) -> str | None:
    if name not in dataset.ncattrs():
        return None
    return str(dataset.getncattr(name)).strip()


def _metadata_text(value: object) -> str:
    return str(value).strip()


def _metadata_token_text(value: object) -> str:
    return _metadata_text(value).lower().replace("-", "_").replace(" ", "_")


def _is_false_metadata_value(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return not bool(value)
    return _metadata_token_text(value) in BOOL_FALSE_VALUES


def _related_metadata_disqualifiers(
    attrs: dict[str, object],
    *,
    excluded_names: set[str] | None = None,
) -> list[str]:
    excluded_names = excluded_names or set()
    disqualifiers: list[str] = []
    for name, value in sorted(attrs.items()):
        if name in excluded_names or not name.startswith("TYWRF_"):
            continue
        if _is_false_metadata_value(value):
            continue
        name_text = _metadata_token_text(name)
        value_text = _metadata_token_text(value)
        if any(
            token in name_text or token in value_text
            for token in BLOCKING_RELATED_METADATA_TOKENS
        ):
            disqualifiers.append(f"{name}={_metadata_text(value)}")
    return disqualifiers


def _wind_tendency_metadata_disqualifiers(attrs: dict[str, object]) -> list[str]:
    disqualifiers: list[str] = []

    source_kind = attrs.get("TYWRF_WIND_TENDENCY_SOURCE_KIND")
    if source_kind is not None:
        source_kind_text = _metadata_token_text(source_kind)
        if any(
            token in source_kind_text
            for token in WIND_TENDENCY_SOURCE_KIND_BLOCKING_TOKENS
        ):
            disqualifiers.append(
                "TYWRF_WIND_TENDENCY_SOURCE_KIND="
                f"{_metadata_text(source_kind)}"
            )

    for name in WIND_TENDENCY_FALSE_EVIDENCE_ATTRS:
        value = attrs.get(name)
        if value is not None and _is_false_metadata_value(value):
            disqualifiers.append(f"{name}={_metadata_text(value)}")

    for name in WIND_TENDENCY_TRUE_BLOCKER_ATTRS:
        value = attrs.get(name)
        if value is not None and _read_metadata_bool(value) is True:
            disqualifiers.append(f"{name}={_metadata_text(value)}")

    return disqualifiers


def _read_metadata_bool(value: object) -> bool | None:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = _metadata_token_text(value)
    if text in BOOL_TRUE_VALUES:
        return True
    if text in BOOL_FALSE_VALUES:
        return False
    return None


def _read_metadata_float(value: object) -> float | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _read_metadata_int(value: object) -> int | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric) or not numeric.is_integer():
        return None
    return int(numeric)


def _pressure_gradient_wind_tendency_metadata_disqualifiers(
    attrs: dict[str, object],
) -> list[str]:
    if not any(
        name.startswith(PRESSURE_GRADIENT_WIND_TENDENCY_PREFIX)
        for name in attrs
    ):
        return []

    disqualifiers: list[str] = []
    for group, names in PRESSURE_GRADIENT_WIND_TENDENCY_REQUIRED_GROUPS.items():
        missing = [name for name in names if name not in attrs]
        if missing:
            disqualifiers.append(
                "pressure-gradient wind tendency missing "
                f"{group} metadata: {', '.join(missing)}"
            )

    mode = attrs.get("TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_MODE")
    if (
        mode is not None
        and _metadata_token_text(mode) != "first_order_constant_alpha"
    ):
        disqualifiers.append(
            "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_MODE="
            f"{_metadata_text(mode)}"
        )

    for name in PRESSURE_GRADIENT_WIND_TENDENCY_REQUIRED_GROUPS["switches"]:
        value = attrs.get(name)
        if value is not None and _read_metadata_bool(value) is not True:
            disqualifiers.append(f"{name}={_metadata_text(value)}")

    source_kind = attrs.get(
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_SOURCE_KIND"
    )
    if source_kind is not None:
        source_kind_text = _metadata_token_text(source_kind)
        if any(
            token in source_kind_text
            for token in WIND_TENDENCY_SOURCE_KIND_BLOCKING_TOKENS
        ):
            disqualifiers.append(
                "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_SOURCE_KIND="
                f"{_metadata_text(source_kind)}"
            )

    units = attrs.get("TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_UNITS")
    if units is not None and _metadata_text(units).lower() != "m s-2":
        disqualifiers.append(
            "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_UNITS="
            f"{_metadata_text(units)}"
        )

    alpha_units = attrs.get(
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_UNITS"
    )
    if (
        alpha_units is not None
        and _metadata_text(alpha_units).lower() != "m3 kg-1"
    ):
        disqualifiers.append(
            "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_UNITS="
            f"{_metadata_text(alpha_units)}"
        )

    alpha_kind = attrs.get(
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_KIND"
    )
    if (
        alpha_kind is not None
        and _metadata_token_text(alpha_kind) != "constant_specific_volume"
    ):
        disqualifiers.append(
            "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_KIND="
            f"{_metadata_text(alpha_kind)}"
        )

    alpha_value = attrs.get(
        "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_VALUE"
    )
    if (
        alpha_value is not None
        and _read_metadata_float(alpha_value) is None
    ):
        disqualifiers.append(
            "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_ALPHA_VALUE="
            f"{_metadata_text(alpha_value)}"
        )

    status = attrs.get("TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_STATUS")
    if status is not None and _metadata_token_text(status) != "ok":
        disqualifiers.append(
            "TYWRF_PRESSURE_GRADIENT_WIND_TENDENCY_STATUS="
            f"{_metadata_text(status)}"
        )

    for name in PRESSURE_GRADIENT_WIND_TENDENCY_REQUIRED_GROUPS["counts"]:
        value = attrs.get(name)
        count = _read_metadata_int(value) if value is not None else None
        if value is not None and (count is None or count < 0):
            disqualifiers.append(f"{name}={_metadata_text(value)}")

    for name in PRESSURE_GRADIENT_WIND_TENDENCY_FALSE_EVIDENCE_ATTRS:
        value = attrs.get(name)
        if value is not None and _read_metadata_bool(value) is not True:
            disqualifiers.append(f"{name}={_metadata_text(value)}")

    for name in PRESSURE_GRADIENT_WIND_TENDENCY_TRUE_BLOCKER_ATTRS:
        value = attrs.get(name)
        if value is not None and _read_metadata_bool(value) is not False:
            disqualifiers.append(f"{name}={_metadata_text(value)}")

    return disqualifiers


def _candidate_metadata_gate(
    candidate_path: Path,
    *,
    allow_validation_gate_only: bool = False,
) -> GateMetric:
    try:
        with netCDF4.Dataset(candidate_path) as dataset:
            attrs = {name: dataset.getncattr(name) for name in dataset.ncattrs()}
            diagnostic_only = _read_bool_attr(dataset, "TYWRF_DIAGNOSTIC_ONLY")
            gate_candidate = _read_bool_attr(dataset, "TYWRF_GATE_CANDIDATE")
            integrator_output = _read_bool_attr(dataset, "TYWRF_INTEGRATOR_OUTPUT")
            validation_gate_only = _read_bool_attr(dataset, "TYWRF_VALIDATION_GATE_ONLY")
            candidate_kind = _read_text_attr(dataset, "TYWRF_CANDIDATE_KIND")
    except OSError as exc:
        return GateMetric(
            name="candidate_metadata",
            status="not_available",
            threshold=METADATA_GATE_THRESHOLD,
            message=f"candidate metadata check failed: {exc}",
        )

    disqualifiers: list[str] = []
    validation_gate_only_allowed = (
        validation_gate_only is True and allow_validation_gate_only
    )
    if diagnostic_only is True:
        disqualifiers.append("TYWRF_DIAGNOSTIC_ONLY=true")
    if gate_candidate is False:
        disqualifiers.append("TYWRF_GATE_CANDIDATE=false")
    elif gate_candidate is not True and not validation_gate_only_allowed:
        disqualifiers.append("TYWRF_GATE_CANDIDATE is not true")
    if validation_gate_only is True and not allow_validation_gate_only:
        disqualifiers.append("TYWRF_VALIDATION_GATE_ONLY=true")
    if integrator_output is False and not validation_gate_only_allowed:
        disqualifiers.append("TYWRF_INTEGRATOR_OUTPUT=false")
    elif integrator_output is not True and not validation_gate_only_allowed:
        disqualifiers.append("TYWRF_INTEGRATOR_OUTPUT is not true")

    kind = _metadata_token_text(candidate_kind) if candidate_kind else ""
    if kind and any(token in kind for token in BLOCKING_CANDIDATE_KIND_TOKENS):
        disqualifiers.append(f"TYWRF_CANDIDATE_KIND={candidate_kind}")
    disqualifiers.extend(_wind_tendency_metadata_disqualifiers(attrs))
    disqualifiers.extend(
        _pressure_gradient_wind_tendency_metadata_disqualifiers(attrs)
    )
    disqualifiers.extend(
        _related_metadata_disqualifiers(
            attrs,
            excluded_names={"TYWRF_CANDIDATE_KIND"},
        )
    )

    if disqualifiers:
        return GateMetric(
            name="candidate_metadata",
            status="failed",
            threshold=METADATA_GATE_THRESHOLD,
            value=1.0,
            message="candidate is not eligible for the default validation gate: "
            + ", ".join(disqualifiers),
        )

    return GateMetric(
        name="candidate_metadata",
        status="passed",
        threshold=METADATA_GATE_THRESHOLD,
        value=0.0,
    )


def _read_2d(dataset: netCDF4.Dataset, variable: str, time_index: int = -1) -> np.ndarray:
    if variable not in dataset.variables:
        raise ValueError(f"missing required variable {variable}")
    data = np.ma.asarray(dataset.variables[variable][:], dtype=np.float64)
    if data.ndim == 2:
        selected = data
    elif data.ndim == 3:
        selected = data[time_index, :, :]
    else:
        raise ValueError(f"variable {variable} must be 2D or Time+2D, got shape {data.shape}")
    return np.asarray(np.ma.filled(selected, np.nan), dtype=np.float64)


def _read_slp_hpa(dataset: netCDF4.Dataset, time_index: int = -1) -> np.ndarray:
    variables_by_lower = {name.lower(): name for name in dataset.variables}
    for candidate in MINIMUM_SLP_VARIABLE_CANDIDATES:
        variable = (
            candidate
            if candidate in dataset.variables
            else variables_by_lower.get(candidate.lower())
        )
        if variable is None:
            continue

        values = _read_2d(dataset, variable, time_index=time_index)
        units = (
            str(getattr(dataset.variables[variable], "units", ""))
            .strip()
            .lower()
            .replace(" ", "")
        )
        if units in {"pa", "pascal", "pascals"}:
            return values / 100.0
        if units in {
            "hpa",
            "hectopascal",
            "hectopascals",
            "mb",
            "mbar",
            "millibar",
            "millibars",
        }:
            return values

        finite = values[np.isfinite(values)]
        if finite.size and float(np.nanmedian(np.abs(finite))) > 2000.0:
            return values / 100.0
        return values
    raise ValueError("missing required real sea-level pressure variable")


def _require_same_shape(arrays: dict[str, np.ndarray]) -> None:
    shapes = {name: values.shape for name, values in arrays.items()}
    if len(set(shapes.values())) != 1:
        detail = ", ".join(f"{name}={shape}" for name, shape in sorted(shapes.items()))
        raise ValueError(f"diagnostic variable shapes do not match: {detail}")


def _minimum_location(values: np.ndarray, latitude: np.ndarray, longitude: np.ndarray) -> tuple[int, int]:
    _require_same_shape({"values": values, "XLAT": latitude, "XLONG": longitude})
    valid = np.isfinite(values) & np.isfinite(latitude) & np.isfinite(longitude)
    if not np.any(valid):
        raise ValueError("no finite samples available for minimum SLP")
    masked = np.where(valid, values, np.inf)
    return tuple(int(index) for index in np.unravel_index(np.argmin(masked), values.shape))


def _maximum_wind(u10: np.ndarray, v10: np.ndarray) -> float:
    _require_same_shape({"U10": u10, "V10": v10})
    speed = np.sqrt(u10 * u10 + v10 * v10)
    finite = speed[np.isfinite(speed)]
    if finite.size == 0:
        raise ValueError("no finite samples available for Vmax10m")
    return float(np.max(finite))


def _haversine_km(
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


def _gate_diagnostics(reference_path: Path, candidate_path: Path) -> list[GateMetric]:
    center_threshold = 20.0
    slp_threshold = 5.0
    vmax_threshold = 5.0

    try:
        with netCDF4.Dataset(reference_path) as reference, netCDF4.Dataset(candidate_path) as candidate:
            ref_slp = _read_slp_hpa(reference)
            cand_slp = _read_slp_hpa(candidate)
            ref_lat = _read_2d(reference, "XLAT")
            ref_lon = _read_2d(reference, "XLONG")
            cand_lat = _read_2d(candidate, "XLAT")
            cand_lon = _read_2d(candidate, "XLONG")
            ref_j, ref_i = _minimum_location(ref_slp, ref_lat, ref_lon)
            cand_j, cand_i = _minimum_location(cand_slp, cand_lat, cand_lon)
            center_error_km = _haversine_km(
                float(ref_lat[ref_j, ref_i]),
                float(ref_lon[ref_j, ref_i]),
                float(cand_lat[cand_j, cand_i]),
                float(cand_lon[cand_j, cand_i]),
            )
            slp_error_hpa = abs(float(cand_slp[cand_j, cand_i]) - float(ref_slp[ref_j, ref_i]))
    except (OSError, ValueError, IndexError) as exc:
        center_metric = GateMetric(
            name="storm_center",
            status="not_available",
            threshold=center_threshold,
            units="km",
            message=f"SLP-based storm center diagnostics failed: {exc}",
        )
        slp_metric = GateMetric(
            name="minimum_slp",
            status="not_available",
            threshold=slp_threshold,
            units="hPa",
            message=f"minimum SLP diagnostics failed: {exc}",
        )
    else:
        center_metric = GateMetric(
            name="storm_center",
            status=_metric_status(center_error_km, center_threshold),
            value=center_error_km,
            threshold=center_threshold,
            units="km",
        )
        slp_metric = GateMetric(
            name="minimum_slp",
            status=_metric_status(slp_error_hpa, slp_threshold),
            value=slp_error_hpa,
            threshold=slp_threshold,
            units="hPa",
        )

    try:
        with netCDF4.Dataset(reference_path) as reference, netCDF4.Dataset(candidate_path) as candidate:
            ref_vmax = _maximum_wind(_read_2d(reference, "U10"), _read_2d(reference, "V10"))
            cand_vmax = _maximum_wind(_read_2d(candidate, "U10"), _read_2d(candidate, "V10"))
            vmax_error_ms = abs(cand_vmax - ref_vmax)
    except (OSError, ValueError, IndexError) as exc:
        vmax_metric = GateMetric(
            name="vmax10m",
            status="not_available",
            threshold=vmax_threshold,
            units="m s-1",
            message=f"Vmax10m diagnostics failed: {exc}",
        )
    else:
        vmax_metric = GateMetric(
            name="vmax10m",
            status=_metric_status(vmax_error_ms, vmax_threshold),
            value=vmax_error_ms,
            threshold=vmax_threshold,
            units="m s-1",
        )

    return [center_metric, slp_metric, vmax_metric]


def evaluate_cycle(
    reference_dir: Path,
    candidate_dir: Path,
    cycle_start: datetime,
    cycle_end: datetime,
    *,
    domain: str = DEFAULT_DOMAIN,
    allow_validation_gate_only: bool = False,
) -> CycleGate:
    reference_path = reference_dir / wrfout_filename(domain, cycle_end)
    candidate_path = candidate_dir / wrfout_filename(domain, cycle_end)

    if not reference_path.exists():
        message = f"missing reference file: {reference_path}"
        return CycleGate(
            status="failed",
            domain=domain,
            start_time=format_wrf_time(cycle_start),
            end_time=format_wrf_time(cycle_end),
            reference=str(reference_path),
            candidate=str(candidate_path),
            fields=_missing_fields(message),
            diagnostics=_missing_diagnostics(message),
            message=message,
        )
    if not candidate_path.exists():
        message = f"missing candidate file: {candidate_path}"
        return CycleGate(
            status="failed",
            domain=domain,
            start_time=format_wrf_time(cycle_start),
            end_time=format_wrf_time(cycle_end),
            reference=str(reference_path),
            candidate=str(candidate_path),
            fields=_missing_fields(message),
            diagnostics=_missing_diagnostics(message),
            message=message,
        )

    metadata_metric = _candidate_metadata_gate(
        candidate_path,
        allow_validation_gate_only=allow_validation_gate_only,
    )
    comparison = compare_files(
        reference_path,
        candidate_path,
        variables=STRICT_CORE_VARIABLES,
        thresholds=GATE_FIELD_THRESHOLDS,
        include_tc_diagnostics=False,
    )
    fields = [_gate_field(variable_comparison) for variable_comparison in comparison.variables]
    diagnostics = [metadata_metric, *_gate_diagnostics(reference_path, candidate_path)]
    failed = any(field.status != "passed" for field in fields) or any(
        metric.status != "passed" for metric in diagnostics
    )

    return CycleGate(
        status="failed" if failed else "passed",
        domain=domain,
        start_time=format_wrf_time(cycle_start),
        end_time=format_wrf_time(cycle_end),
        reference=str(reference_path),
        candidate=str(candidate_path),
        fields=fields,
        diagnostics=diagnostics,
    )


def evaluate_cycles(
    reference_dir: Path,
    candidate_dir: Path,
    start: str | datetime,
    *,
    end: str | datetime | None = None,
    hours: int | None = None,
    interval_hours: int | None = DEFAULT_INTERVAL_HOURS,
    interval_minutes: int | None = None,
    domain: str = DEFAULT_DOMAIN,
    allow_validation_gate_only: bool = False,
) -> CycleGateReport:
    if domain not in {"d01", "d02"}:
        raise ValueError(f"unsupported domain: {domain}")
    start_time = parse_wrf_time(start) if isinstance(start, str) else start
    end_time = resolve_end_time(start_time, end=end, hours=hours)
    resolved_interval_minutes = _resolve_interval_minutes(
        interval_hours=interval_hours,
        interval_minutes=interval_minutes,
    )
    end_times = cycle_end_times(
        start_time,
        end_time,
        interval_hours=None,
        interval_minutes=resolved_interval_minutes,
    )

    cycle_start = start_time
    cycles = []
    for cycle_end in end_times:
        cycle = evaluate_cycle(
            reference_dir,
            candidate_dir,
            cycle_start,
            cycle_end,
            domain=domain,
            allow_validation_gate_only=allow_validation_gate_only,
        )
        cycles.append(cycle)
        if _candidate_metadata_failed(cycle):
            break
        cycle_start = cycle_end

    passed = sum(1 for cycle in cycles if cycle.status == "passed")
    failed = len(cycles) - passed
    first_failure = _first_failure_summary(cycles)
    return CycleGateReport(
        status="passed" if failed == 0 else "failed",
        domain=domain,
        reference_dir=str(reference_dir),
        candidate_dir=str(candidate_dir),
        start_time=format_wrf_time(start_time),
        end_time=format_wrf_time(end_time),
        interval_hours=resolved_interval_minutes // 60,
        interval_minutes=resolved_interval_minutes,
        cycles=cycles,
        summary={"total": len(cycles), "passed": passed, "failed": failed},
        first_failure=first_failure,
    )


def _candidate_metadata_failed(cycle: CycleGate) -> bool:
    return any(
        metric.name == "candidate_metadata" and metric.status != "passed"
        for metric in cycle.diagnostics
    )


def _first_failure_summary(cycles: list[CycleGate]) -> dict[str, object | None] | None:
    for index, cycle in enumerate(cycles, start=1):
        if cycle.status == "passed":
            continue

        failed_field = next((field for field in cycle.fields if field.status != "passed"), None)
        failed_diagnostic = next(
            (metric for metric in cycle.diagnostics if metric.status != "passed"),
            None,
        )
        return {
            "cycle_index": index,
            "start_time": cycle.start_time,
            "end_time": cycle.end_time,
            "field": failed_field.variable if failed_field else None,
            "field_status": failed_field.status if failed_field else None,
            "diagnostic": failed_diagnostic.name if failed_diagnostic else None,
            "diagnostic_status": failed_diagnostic.status if failed_diagnostic else None,
            "message": cycle.message
            or (failed_field.message if failed_field and failed_field.message else None)
            or (
                failed_diagnostic.message
                if failed_diagnostic and failed_diagnostic.message
                else None
            ),
        }
    return None


def _strict_json_value(value):
    if is_dataclass(value) and not isinstance(value, type):
        return _strict_json_value(asdict(value))
    if isinstance(value, dict):
        return {key: _strict_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_strict_json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def report_to_dict(report: CycleGateReport) -> dict:
    return _strict_json_value(report)


def report_to_json(report: CycleGateReport, *, pretty: bool = False) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True, type=Path, help="Directory with WRF reference files")
    parser.add_argument("--candidate-dir", required=True, type=Path, help="Directory with TyWRF-Core output")
    parser.add_argument("--start", required=True, help="First cycle start time, for example 2025-07-26_00:00:00")
    parser.add_argument("--end", help="Final cycle end time, for example 2025-07-26_06:00:00")
    parser.add_argument("--hours", type=int, help="Duration from --start to the final cycle end")
    parser.add_argument("--interval", type=int, default=DEFAULT_INTERVAL_HOURS, help="Cycle interval in hours")
    parser.add_argument(
        "--interval-minutes",
        type=int,
        help="Cycle interval in minutes; overrides --interval for 10 min progressive validation gates",
    )
    parser.add_argument("--domain", choices=("d01", "d02"), default=DEFAULT_DOMAIN, help="Domain to gate")
    parser.add_argument(
        "--allow-validation-gate-only",
        action="store_true",
        help=(
            "Allow files marked TYWRF_VALIDATION_GATE_ONLY=true to exercise the gate "
            "implementation; default strict mode fails them as non-integrator candidates."
        ),
    )
    parser.add_argument("--output", type=Path, help="Write the JSON gate report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = evaluate_cycles(
            args.reference_dir,
            args.candidate_dir,
            args.start,
            end=args.end,
            hours=args.hours,
            interval_hours=args.interval,
            interval_minutes=args.interval_minutes,
            domain=args.domain,
            allow_validation_gate_only=args.allow_validation_gate_only,
        )
    except ValueError as exc:
        parser.error(str(exc))

    report_json = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0 if report.status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
