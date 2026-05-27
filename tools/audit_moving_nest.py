#!/usr/bin/env python
"""Audit d02 moving-nest placement over the first KROSA 10 minute segment."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime
import json
import math
from pathlib import Path
import re
from typing import Any

import netCDF4
import numpy as np


WRF_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"
DEFAULT_DOMAIN = "d02"
DEFAULT_PARENT_GRID_RATIO = 5
EXPECTED_D02_DX_M = 2000.0
EXPECTED_D02_DY_M = 2000.0
RESOLUTION_TOLERANCE_M = 1.0e-6
KROSA_START = "2025-07-26_00:00:00"
KROSA_END = "2025-07-26_00:10:00"
KROSA_EXPECTED_START = {"i": 114, "j": 96}
KROSA_EXPECTED_END = {"i": 126, "j": 103}
TIME_RE = re.compile(r"(\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})")
WRFOUT_RE = re.compile(r"wrfout_(?P<domain>d\d\d)_(?P<time>\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})$")


@dataclass(frozen=True)
class NestPose:
    file: str
    i_parent_start: int | None
    j_parent_start: int | None
    dx_m: float | None
    dy_m: float | None
    grid_shape: tuple[int, int] | None
    center_index: tuple[int, int] | None
    center_lat: float | None
    center_lon: float | None
    mean_lat: float | None
    mean_lon: float | None


@dataclass(frozen=True)
class MovingNestAudit:
    status: str
    start_file: str
    end_file: str
    domain: str
    start_time: str | None
    end_time: str | None
    parent_grid_ratio: int
    summary: dict[str, Any]
    log_events: dict[str, Any]


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


def infer_wrfout_metadata(path: Path) -> tuple[str | None, str | None]:
    match = WRFOUT_RE.search(path.name)
    if not match:
        return None, None
    return match.group("domain"), match.group("time")


def resolve_moving_nest_files(
    *,
    start_file: Path | None = None,
    end_file: Path | None = None,
    reference_dir: Path | None = None,
    domain: str = DEFAULT_DOMAIN,
    start: str | datetime | None = None,
    end: str | datetime | None = None,
) -> tuple[Path, Path, str, str | None, str | None]:
    explicit = start_file is not None or end_file is not None
    named = reference_dir is not None or start is not None or end is not None
    if explicit and named:
        raise ValueError(
            "use either --start-file/--end-file or --reference-dir/--domain/--start/--end"
        )
    if domain != DEFAULT_DOMAIN:
        raise ValueError("moving-nest audit currently supports d02 only")

    if explicit:
        if start_file is None or end_file is None:
            raise ValueError("--start-file and --end-file must be provided together")
        inferred_domain, inferred_start = infer_wrfout_metadata(start_file)
        _, inferred_end = infer_wrfout_metadata(end_file)
        return start_file, end_file, inferred_domain or domain, inferred_start, inferred_end

    if reference_dir is None or start is None or end is None:
        raise ValueError(
            "provide --start-file/--end-file or --reference-dir, --domain, --start, and --end"
        )

    start_time = parse_wrf_time(start) if isinstance(start, str) else start
    end_time = parse_wrf_time(end) if isinstance(end, str) else end
    if end_time <= start_time:
        raise ValueError("end time must be after start time")

    return (
        reference_dir / wrfout_filename(domain, start_time),
        reference_dir / wrfout_filename(domain, end_time),
        domain,
        format_wrf_time(start_time),
        format_wrf_time(end_time),
    )


def _as_json_scalar(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return _as_json_scalar(value.reshape(-1)[0])
        return [_as_json_scalar(item) for item in value.tolist()]
    return value


def _numeric_attr(dataset: netCDF4.Dataset, name: str) -> float | int | None:
    if name not in dataset.ncattrs():
        return None
    value = _as_json_scalar(dataset.getncattr(name))
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if numeric.is_integer():
        return int(numeric)
    return numeric


def _read_2d_variable(dataset: netCDF4.Dataset, name: str) -> np.ndarray | None:
    if name not in dataset.variables:
        return None
    data = np.ma.asarray(dataset.variables[name][:], dtype=np.float64)
    if data.ndim >= 3:
        data = data[0]
    data = np.ma.squeeze(data)
    if data.ndim != 2:
        return None
    return np.asarray(np.ma.filled(data, np.nan), dtype=np.float64)


def _finite_mean(data: np.ndarray) -> float | None:
    finite = np.asarray(data, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None
    return float(np.mean(finite))


def read_pose(path: Path) -> tuple[NestPose, dict[str, np.ndarray | None]]:
    with netCDF4.Dataset(path) as dataset:
        lat = _read_2d_variable(dataset, "XLAT")
        lon = _read_2d_variable(dataset, "XLONG")
        hgt = _read_2d_variable(dataset, "HGT")
        grid_shape = tuple(int(v) for v in lat.shape) if lat is not None else None
        center_index = None
        center_lat = None
        center_lon = None
        if lat is not None and lon is not None and lat.shape == lon.shape:
            j = int(lat.shape[0] // 2)
            i = int(lat.shape[1] // 2)
            center_index = (j, i)
            center_lat = float(lat[j, i]) if np.isfinite(lat[j, i]) else None
            center_lon = float(lon[j, i]) if np.isfinite(lon[j, i]) else None

        pose = NestPose(
            file=str(path),
            i_parent_start=_coerce_int(_numeric_attr(dataset, "I_PARENT_START")),
            j_parent_start=_coerce_int(_numeric_attr(dataset, "J_PARENT_START")),
            dx_m=_coerce_float(_numeric_attr(dataset, "DX")),
            dy_m=_coerce_float(_numeric_attr(dataset, "DY")),
            grid_shape=grid_shape,
            center_index=center_index,
            center_lat=center_lat,
            center_lon=center_lon,
            mean_lat=_finite_mean(lat) if lat is not None else None,
            mean_lon=_finite_mean(lon) if lon is not None else None,
        )
        return pose, {"XLAT": lat, "XLONG": lon, "HGT": hgt}


def _coerce_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _coerce_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _sub_dict(end_value: float | None, start_value: float | None) -> float | None:
    if end_value is None or start_value is None:
        return None
    return float(end_value - start_value)


def _parent_delta(start_pose: NestPose, end_pose: NestPose) -> dict[str, int | None]:
    i_delta = (
        None
        if start_pose.i_parent_start is None or end_pose.i_parent_start is None
        else end_pose.i_parent_start - start_pose.i_parent_start
    )
    j_delta = (
        None
        if start_pose.j_parent_start is None or end_pose.j_parent_start is None
        else end_pose.j_parent_start - start_pose.j_parent_start
    )
    return {"i": i_delta, "j": j_delta}


def _child_cell_delta(
    parent_delta: dict[str, int | None],
    ratio: int,
) -> dict[str, int | None]:
    return {
        "i": None if parent_delta["i"] is None else parent_delta["i"] * ratio,
        "j": None if parent_delta["j"] is None else parent_delta["j"] * ratio,
    }


def _moved(parent_delta: dict[str, int | None]) -> bool | None:
    if parent_delta["i"] is None or parent_delta["j"] is None:
        return None
    return parent_delta["i"] != 0 or parent_delta["j"] != 0


def _haversine_km(
    start_lat: float | None,
    start_lon: float | None,
    end_lat: float | None,
    end_lon: float | None,
) -> float | None:
    if None in (start_lat, start_lon, end_lat, end_lon):
        return None
    radius_km = 6371.0088
    lat1 = math.radians(float(start_lat))
    lon1 = math.radians(float(start_lon))
    lat2 = math.radians(float(end_lat))
    lon2 = math.radians(float(end_lon))
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    return float(2.0 * radius_km * math.asin(math.sqrt(a)))


def _hgt_delta(start_hgt: np.ndarray | None, end_hgt: np.ndarray | None) -> dict[str, Any]:
    if start_hgt is None or end_hgt is None:
        return {
            "status": "missing",
            "rms_m": None,
            "max_abs_m": None,
            "finite_pair_count": 0,
            "total_pair_count": None,
        }
    if start_hgt.shape != end_hgt.shape:
        return {
            "status": "shape_mismatch",
            "start_shape": tuple(int(v) for v in start_hgt.shape),
            "end_shape": tuple(int(v) for v in end_hgt.shape),
            "rms_m": None,
            "max_abs_m": None,
            "finite_pair_count": 0,
            "total_pair_count": int(end_hgt.size),
        }
    finite = np.isfinite(start_hgt) & np.isfinite(end_hgt)
    total = int(end_hgt.size)
    if not np.any(finite):
        return {
            "status": "no_valid_pairs",
            "rms_m": None,
            "max_abs_m": None,
            "finite_pair_count": 0,
            "total_pair_count": total,
        }
    delta = end_hgt[finite] - start_hgt[finite]
    return {
        "status": "computed",
        "rms_m": float(math.sqrt(np.mean(delta * delta))),
        "max_abs_m": float(np.max(np.abs(delta))),
        "finite_pair_count": int(delta.size),
        "total_pair_count": total,
    }


def _resolution_check(start_pose: NestPose, end_pose: NestPose, domain: str) -> dict[str, Any]:
    values = {
        "start_dx_m": start_pose.dx_m,
        "start_dy_m": start_pose.dy_m,
        "end_dx_m": end_pose.dx_m,
        "end_dy_m": end_pose.dy_m,
    }
    if any(value is None for value in values.values()):
        status = "missing"
        is_2km = None
    elif domain != DEFAULT_DOMAIN:
        status = "unsupported_domain"
        is_2km = None
    else:
        checks = (
            abs(float(start_pose.dx_m) - EXPECTED_D02_DX_M) <= RESOLUTION_TOLERANCE_M,
            abs(float(start_pose.dy_m) - EXPECTED_D02_DY_M) <= RESOLUTION_TOLERANCE_M,
            abs(float(end_pose.dx_m) - EXPECTED_D02_DX_M) <= RESOLUTION_TOLERANCE_M,
            abs(float(end_pose.dy_m) - EXPECTED_D02_DY_M) <= RESOLUTION_TOLERANCE_M,
        )
        is_2km = all(checks)
        status = "ok" if is_2km else "not_2km"
    return {
        "status": status,
        "d02_2km": is_2km,
        "expected_dx_m": EXPECTED_D02_DX_M,
        "expected_dy_m": EXPECTED_D02_DY_M,
        "tolerance_m": RESOLUTION_TOLERANCE_M,
        **values,
    }


def _krosa_flags(
    *,
    domain: str,
    start_time: str | None,
    end_time: str | None,
    start_pose: NestPose,
    end_pose: NestPose,
    parent_delta: dict[str, int | None],
) -> dict[str, Any]:
    expected_delta = {
        "i": KROSA_EXPECTED_END["i"] - KROSA_EXPECTED_START["i"],
        "j": KROSA_EXPECTED_END["j"] - KROSA_EXPECTED_START["j"],
    }
    applicable = domain == DEFAULT_DOMAIN and start_time == KROSA_START and end_time == KROSA_END
    if not applicable:
        return {
            "applicable": False,
            "expected_start": KROSA_EXPECTED_START,
            "expected_end": KROSA_EXPECTED_END,
            "expected_parent_delta": expected_delta,
            "start_pose_matches": None,
            "end_pose_matches": None,
            "parent_delta_matches": None,
            "matches_expected_krosa_first_10min": None,
        }

    start_matches = (
        start_pose.i_parent_start == KROSA_EXPECTED_START["i"]
        and start_pose.j_parent_start == KROSA_EXPECTED_START["j"]
    )
    end_matches = (
        end_pose.i_parent_start == KROSA_EXPECTED_END["i"]
        and end_pose.j_parent_start == KROSA_EXPECTED_END["j"]
    )
    delta_matches = parent_delta == expected_delta
    return {
        "applicable": True,
        "expected_start": KROSA_EXPECTED_START,
        "expected_end": KROSA_EXPECTED_END,
        "expected_parent_delta": expected_delta,
        "start_pose_matches": start_matches,
        "end_pose_matches": end_matches,
        "parent_delta_matches": delta_matches,
        "matches_expected_krosa_first_10min": start_matches and end_matches and delta_matches,
    }


def parse_log_events(
    log_file: Path | None,
    *,
    start_time: str | None,
    end_time: str | None,
) -> dict[str, Any]:
    if log_file is None:
        return {"status": "not_requested", "file": None, "events": []}
    if start_time is None or end_time is None:
        return {
            "status": "unavailable",
            "file": str(log_file),
            "events": [],
            "message": "start/end times are unavailable for log filtering",
        }
    if not log_file.exists():
        return {
            "status": "unavailable",
            "file": str(log_file),
            "events": [],
            "message": "log file does not exist",
        }

    start_dt = parse_wrf_time(start_time)
    end_dt = parse_wrf_time(end_time)
    events_by_time: dict[str, dict[str, Any]] = {}
    ordered_times: list[str] = []
    last_move_time: str | None = None

    for line_number, line in enumerate(log_file.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        time_match = TIME_RE.search(line)
        line_time = time_match.group(1) if time_match else None
        if line_time is not None:
            line_dt = parse_wrf_time(line_time)
            if line_dt < start_dt or line_dt > end_dt:
                last_move_time = None
                continue
            event = events_by_time.setdefault(line_time, {"time": line_time, "line_numbers": []})
            if line_time not in ordered_times:
                ordered_times.append(line_time)
            event["line_numbers"].append(line_number)
            if "vortex center (in nest x and y)" in line:
                values = _trailing_floats(line)
                if len(values) >= 2:
                    event["vortex_center_nest_xy"] = {"x": values[-2], "y": values[-1]}
            elif "grid   center (in nest x and y)" in line:
                values = _trailing_floats(line)
                if len(values) >= 2:
                    event["grid_center_nest_xy"] = {"x": values[-2], "y": values[-1]}
            elif "disp" in line:
                values = _trailing_floats(line)
                if len(values) >= 2:
                    event["disp"] = {"i": values[-2], "j": values[-1]}
            elif "move (rel cd)" in line:
                values = _trailing_ints(line)
                if len(values) >= 2:
                    event["move_rel_cd"] = {"i": values[-2], "j": values[-1]}
                    last_move_time = line_time
            continue

        if "moving" in line and last_move_time is not None:
            values = _trailing_ints(line)
            if len(values) >= 3:
                event = events_by_time[last_move_time]
                event.setdefault("line_numbers", []).append(line_number)
                event["applied_move"] = {
                    "domain": values[-3],
                    "i": values[-2],
                    "j": values[-1],
                }
            last_move_time = None

    events = [events_by_time[time] for time in ordered_times]
    applied_moves = [
        event["applied_move"]
        for event in events
        if isinstance(event.get("applied_move"), dict)
    ]
    net_applied = {
        "i": int(sum(move["i"] for move in applied_moves)),
        "j": int(sum(move["j"] for move in applied_moves)),
    }
    return {
        "status": "available" if events else "unavailable",
        "file": str(log_file),
        "start_time": start_time,
        "end_time": end_time,
        "event_count": len(events),
        "applied_move_count": len(applied_moves),
        "net_applied_parent_delta": net_applied,
        "events": events,
        "message": None if events else "no moving-nest log events matched the requested window",
    }


def _trailing_floats(line: str) -> list[float]:
    return [float(value) for value in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", line)]


def _trailing_ints(line: str) -> list[int]:
    return [int(value) for value in re.findall(r"[-+]?\d+", line)]


def audit_moving_nest(
    start_path: Path,
    end_path: Path,
    *,
    domain: str = DEFAULT_DOMAIN,
    start_time: str | None = None,
    end_time: str | None = None,
    parent_grid_ratio: int = DEFAULT_PARENT_GRID_RATIO,
    log_file: Path | None = None,
) -> MovingNestAudit:
    start_pose, start_fields = read_pose(start_path)
    end_pose, end_fields = read_pose(end_path)
    parent_delta = _parent_delta(start_pose, end_pose)
    child_delta = _child_cell_delta(parent_delta, parent_grid_ratio)
    moved = _moved(parent_delta)
    resolution = _resolution_check(start_pose, end_pose, domain)
    hgt_delta = _hgt_delta(start_fields["HGT"], end_fields["HGT"])
    drift_km = _haversine_km(
        start_pose.center_lat,
        start_pose.center_lon,
        end_pose.center_lat,
        end_pose.center_lon,
    )
    log_events = parse_log_events(log_file, start_time=start_time, end_time=end_time)
    if (
        log_events.get("status") == "available"
        and parent_delta["i"] is not None
        and parent_delta["j"] is not None
    ):
        log_delta = log_events.get("net_applied_parent_delta", {})
        log_events["net_applied_matches_parent_delta"] = (
            log_delta.get("i") == parent_delta["i"] and log_delta.get("j") == parent_delta["j"]
        )

    summary = {
        "start_pose": start_pose,
        "end_pose": end_pose,
        "parent_delta": parent_delta,
        "child_cell_delta": child_delta,
        "drift_km": drift_km,
        "moved": moved,
        "resolution": resolution,
        "coordinate_delta": {
            "center_lat": _sub_dict(end_pose.center_lat, start_pose.center_lat),
            "center_lon": _sub_dict(end_pose.center_lon, start_pose.center_lon),
            "mean_lat": _sub_dict(end_pose.mean_lat, start_pose.mean_lat),
            "mean_lon": _sub_dict(end_pose.mean_lon, start_pose.mean_lon),
        },
        "hgt_delta": hgt_delta,
        "validation_flags": {
            "expected_krosa_first_10min": _krosa_flags(
                domain=domain,
                start_time=start_time,
                end_time=end_time,
                start_pose=start_pose,
                end_pose=end_pose,
                parent_delta=parent_delta,
            )
        },
    }
    status = "computed"
    if resolution["status"] in {"missing", "not_2km"} or hgt_delta["status"] != "computed":
        status = "incomplete" if resolution["status"] == "missing" else "computed_with_flags"

    return MovingNestAudit(
        status=status,
        start_file=str(start_path),
        end_file=str(end_path),
        domain=domain,
        start_time=start_time,
        end_time=end_time,
        parent_grid_ratio=parent_grid_ratio,
        summary=summary,
        log_events=log_events,
    )


def audit_reference_moving_nest(
    reference_dir: Path,
    *,
    domain: str = DEFAULT_DOMAIN,
    start: str | datetime = KROSA_START,
    end: str | datetime = KROSA_END,
    parent_grid_ratio: int = DEFAULT_PARENT_GRID_RATIO,
    log_file: Path | None = None,
) -> MovingNestAudit:
    start_path, end_path, resolved_domain, start_time, end_time = resolve_moving_nest_files(
        reference_dir=reference_dir,
        domain=domain,
        start=start,
        end=end,
    )
    return audit_moving_nest(
        start_path,
        end_path,
        domain=resolved_domain,
        start_time=start_time,
        end_time=end_time,
        parent_grid_ratio=parent_grid_ratio,
        log_file=log_file,
    )


def _strict_json_value(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _strict_json_value(asdict(value))
    if isinstance(value, dict):
        return {key: _strict_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_strict_json_value(item) for item in value]
    if isinstance(value, np.generic):
        return _strict_json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def report_to_dict(report: MovingNestAudit) -> dict[str, Any]:
    return _strict_json_value(report)


def report_to_json(report: MovingNestAudit, *, pretty: bool = False) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_argument_group("input selection")
    input_group.add_argument("--start-file", type=Path, help="Cycle start d02 WRF NetCDF file")
    input_group.add_argument("--end-file", type=Path, help="Cycle end d02 WRF NetCDF file")
    input_group.add_argument("--reference-dir", type=Path, help="Directory with WRF reference files")
    input_group.add_argument("--domain", choices=(DEFAULT_DOMAIN,), default=DEFAULT_DOMAIN)
    input_group.add_argument("--start", help="Cycle start time, for example 2025-07-26_00:00:00")
    input_group.add_argument("--end", help="Cycle end time, for example 2025-07-26_00:10:00")
    parser.add_argument(
        "--parent-grid-ratio",
        type=int,
        default=DEFAULT_PARENT_GRID_RATIO,
        help="Parent-to-child horizontal grid ratio; KROSA d02 uses 5",
    )
    parser.add_argument("--log-file", type=Path, help="Optional WRF rsl log to mine for move lines")
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        start_file, end_file, domain, start_time, end_time = resolve_moving_nest_files(
            start_file=args.start_file,
            end_file=args.end_file,
            reference_dir=args.reference_dir,
            domain=args.domain,
            start=args.start,
            end=args.end,
        )
    except ValueError as exc:
        parser.error(str(exc))

    report = audit_moving_nest(
        start_file,
        end_file,
        domain=domain,
        start_time=start_time,
        end_time=end_time,
        parent_grid_ratio=args.parent_grid_ratio,
        log_file=args.log_file,
    )
    report_json = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_json + "\n", encoding="utf-8")
    print(report_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
