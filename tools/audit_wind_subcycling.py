#!/usr/bin/env python
"""Audit U/V wind subcycling candidates against WRF reference output."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Mapping
from contextlib import ExitStack
import json
import math
from pathlib import Path
from typing import Any

import netCDF4
import numpy as np


DEFAULT_VARIABLES = ("U", "V")
DEFAULT_BOUNDARY_BANDS = (0, 1, 2, 5, 10, 20)
WIND_SUBCYCLING_ATTRS = (
    "TYWRF_WIND_TENDENCY_SUBSTEP_COUNT",
    "TYWRF_WIND_TENDENCY_SUBSTEP_DT_SECONDS",
    "TYWRF_WIND_TENDENCY_TOTAL_SECONDS",
    "TYWRF_WIND_TENDENCY_ADVECTING_VELOCITY_MODE",
    "TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS",
    "TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION",
    "TYWRF_WIND_TENDENCY_ADVECTION_FORM",
)


def _attr_to_json(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, (list, tuple)):
        return [_attr_to_json(item) for item in value]
    if hasattr(value, "tolist"):
        return _attr_to_json(value.tolist())
    if hasattr(value, "item"):
        return _attr_to_json(value.item())
    return str(value).strip()


def _json_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _read_variable_data(
    dataset: netCDF4.Dataset,
    name: str,
    *,
    time_index: int,
) -> Any:
    variable = dataset.variables[name]
    if variable.ndim > 0 and variable.dimensions[0] == "Time":
        length = int(variable.shape[0])
        index = time_index if time_index >= 0 else length + time_index
        if index < 0 or index >= length:
            raise IndexError(
                f"{name} has Time length {length}, cannot read time_index={time_index}"
            )
        return variable[index, ...]
    return variable[:]


def _finite_values(data: Any) -> tuple[np.ndarray, np.ndarray]:
    masked = np.ma.masked_invalid(np.ma.asarray(data, dtype=np.float64))
    values = np.asarray(masked.filled(np.nan), dtype=np.float64)
    finite = np.isfinite(values) & ~np.ma.getmaskarray(masked)
    return values, finite


def error_metrics(
    reference: Any,
    candidate: Any,
    *,
    selection_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    ref_values, ref_finite = _finite_values(reference)
    cand_values, cand_finite = _finite_values(candidate)
    if ref_values.shape != cand_values.shape:
        raise ValueError(f"shape mismatch: {ref_values.shape} != {cand_values.shape}")

    if selection_mask is None:
        selected = np.ones(ref_values.shape, dtype=bool)
    else:
        selected = np.asarray(selection_mask, dtype=bool)
        if selected.shape != ref_values.shape:
            selected = np.broadcast_to(selected, ref_values.shape)

    valid = selected & ref_finite & cand_finite
    valid_count = int(np.count_nonzero(valid))
    selected_count = int(np.count_nonzero(selected))
    if valid_count == 0:
        return {
            "rmse": None,
            "normalized_rmse": None,
            "max_abs_error": None,
            "valid_count": 0,
            "count": 0,
            "selected_count": selected_count,
        }

    ref_sample = ref_values[valid]
    diff = cand_values[valid] - ref_sample
    rmse = float(math.sqrt(float(np.mean(diff * diff))))
    scale = float(math.sqrt(float(np.mean(ref_sample * ref_sample))))
    if scale == 0.0:
        normalized = 0.0 if rmse == 0.0 else None
    else:
        normalized = float(rmse / scale)
    return {
        "rmse": rmse,
        "normalized_rmse": normalized,
        "max_abs_error": float(np.max(np.abs(diff))),
        "valid_count": valid_count,
        "count": valid_count,
        "selected_count": selected_count,
    }


def _zero_delta_metrics(candidate: Any) -> dict[str, Any]:
    cand_values, cand_finite = _finite_values(candidate)
    valid_count = int(np.count_nonzero(cand_finite))
    return {
        "status": "computed",
        "changed_count": 0,
        "rmse": 0.0,
        "normalized_rmse": 0.0,
        "max_abs_diff": 0.0,
        "valid_count": valid_count,
        "count": valid_count,
        # D88 compatibility aliases.
        "differing_count": 0,
        "max_abs_delta": 0.0,
    }


def baseline_delta_metrics(candidate: Any, baseline: Any) -> dict[str, Any]:
    cand_values, cand_finite = _finite_values(candidate)
    base_values, base_finite = _finite_values(baseline)
    if cand_values.shape != base_values.shape:
        raise ValueError(f"shape mismatch: {cand_values.shape} != {base_values.shape}")

    valid = cand_finite & base_finite
    valid_count = int(np.count_nonzero(valid))
    if valid_count == 0:
        return {
            "status": "no_valid_samples",
            "changed_count": 0,
            "rmse": None,
            "normalized_rmse": None,
            "max_abs_diff": None,
            "valid_count": 0,
            "count": 0,
            # D88 compatibility aliases.
            "differing_count": 0,
            "max_abs_delta": None,
        }
    delta = cand_values[valid] - base_values[valid]
    rmse = float(math.sqrt(float(np.mean(delta * delta))))
    scale = float(math.sqrt(float(np.mean(base_values[valid] * base_values[valid]))))
    if scale == 0.0:
        normalized = 0.0 if rmse == 0.0 else None
    else:
        normalized = float(rmse / scale)
    max_abs = float(np.max(np.abs(delta)))
    changed_count = int(np.count_nonzero(delta != 0.0))
    return {
        "status": "computed",
        "changed_count": changed_count,
        "rmse": rmse,
        "normalized_rmse": normalized,
        "max_abs_diff": max_abs,
        "valid_count": valid_count,
        "count": valid_count,
        # D88 compatibility aliases.
        "differing_count": changed_count,
        "max_abs_delta": max_abs,
    }


def _horizontal_boundary_distance(shape: tuple[int, ...]) -> np.ndarray | None:
    if len(shape) < 2:
        return None
    ny = int(shape[-2])
    nx = int(shape[-1])
    y = np.arange(ny, dtype=np.int64)[:, None]
    x = np.arange(nx, dtype=np.int64)[None, :]
    return np.minimum(np.minimum(y, x), np.minimum(ny - 1 - y, nx - 1 - x))


def _broadcast_horizontal_mask(mask_2d: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    prefix = (1,) * (len(shape) - 2)
    return np.broadcast_to(mask_2d.reshape(prefix + mask_2d.shape), shape)


def boundary_band_statistics(
    reference: Any,
    candidate: Any,
    *,
    bands: Iterable[int] = DEFAULT_BOUNDARY_BANDS,
) -> dict[str, Any]:
    ref_values = np.asarray(reference)
    shape = tuple(int(value) for value in ref_values.shape)
    distance = _horizontal_boundary_distance(shape)
    if distance is None:
        return {
            "status": "no_horizontal_dimension",
            "horizontal_shape": None,
            "bands": [],
        }

    entries: list[dict[str, Any]] = []
    for band in bands:
        boundary_2d = distance <= int(band)
        interior_2d = ~boundary_2d
        entries.append(
            {
                "band": int(band),
                "boundary": error_metrics(
                    reference,
                    candidate,
                    selection_mask=_broadcast_horizontal_mask(boundary_2d, shape),
                ),
                "interior": error_metrics(
                    reference,
                    candidate,
                    selection_mask=_broadcast_horizontal_mask(interior_2d, shape),
                ),
            }
        )
    return {
        "status": "available",
        "horizontal_shape": [shape[-2], shape[-1]],
        "distance_definition": "min(i, j, nx - 1 - i, ny - 1 - j)",
        "bands": entries,
    }


def vertical_level_statistics(reference: Any, candidate: Any) -> dict[str, Any]:
    ref_values = np.asarray(reference)
    shape = tuple(int(value) for value in ref_values.shape)
    if len(shape) < 3:
        return {
            "status": "no_vertical_dimension",
            "k_axis": None,
            "levels": [
                {
                    "k": None,
                    "metrics": error_metrics(reference, candidate),
                }
            ],
        }

    k_axis = len(shape) - 3
    levels = []
    for k in range(shape[k_axis]):
        levels.append(
            {
                "k": int(k),
                "metrics": error_metrics(
                    np.take(reference, k, axis=k_axis),
                    np.take(candidate, k, axis=k_axis),
                ),
            }
        )
    return {
        "status": "available",
        "k_axis": int(k_axis),
        "level_count": int(shape[k_axis]),
        "levels": levels,
    }


def _candidate_metadata(candidate: netCDF4.Dataset) -> dict[str, Any]:
    attrs = {
        name: _attr_to_json(candidate.getncattr(name))
        for name in WIND_SUBCYCLING_ATTRS
        if name in candidate.ncattrs()
    }
    return {
        "global_attrs": attrs,
        "present_attrs": sorted(attrs),
        "missing_attrs": [name for name in WIND_SUBCYCLING_ATTRS if name not in attrs],
    }


def _shape_of(data: Any | None) -> list[int] | None:
    if data is None:
        return None
    return [int(value) for value in np.asarray(data).shape]


def _same_existing_path(left: Path | None, right: Path | None) -> bool:
    if left is None or right is None:
        return False
    try:
        return left.resolve() == right.resolve()
    except OSError:
        return left.absolute() == right.absolute()


def _read_optional_variable(
    dataset: netCDF4.Dataset,
    variable: str,
    *,
    time_index: int,
) -> tuple[str, Any | None, str | None]:
    if variable not in dataset.variables:
        return "missing", None, "variable missing"
    try:
        return "available", _read_variable_data(dataset, variable, time_index=time_index), None
    except (IndexError, TypeError, ValueError) as exc:
        return "read_error", None, str(exc)


def _delta_vs_baseline(
    candidate_data: Any | None,
    baseline: netCDF4.Dataset | None,
    variable: str,
    *,
    candidate_status: str,
    baseline_is_candidate: bool,
    time_index: int,
) -> dict[str, Any]:
    if baseline is None:
        return {
            "status": "not_requested",
            "changed_count": None,
            "rmse": None,
            "normalized_rmse": None,
            "max_abs_diff": None,
            "differing_count": None,
            "max_abs_delta": None,
        }
    if candidate_status != "available" or candidate_data is None:
        return {
            "status": "candidate_unavailable",
            "changed_count": None,
            "rmse": None,
            "normalized_rmse": None,
            "max_abs_diff": None,
            "differing_count": None,
            "max_abs_delta": None,
        }
    if baseline_is_candidate:
        return {
            "baseline_present": True,
            "candidate_shape": _shape_of(candidate_data),
            "baseline_shape": _shape_of(candidate_data),
            "same_path_as_baseline": True,
            **_zero_delta_metrics(candidate_data),
        }
    baseline_status, baseline_data, message = _read_optional_variable(
        baseline,
        variable,
        time_index=time_index,
    )
    result = {
        "baseline_present": baseline_status == "available",
        "candidate_shape": _shape_of(candidate_data),
        "baseline_shape": _shape_of(baseline_data),
        "same_path_as_baseline": False,
        "changed_count": None,
        "rmse": None,
        "normalized_rmse": None,
        "max_abs_diff": None,
        "differing_count": None,
        "max_abs_delta": None,
    }
    if baseline_status != "available":
        return {
            **result,
            "status": "missing_baseline_variable"
            if baseline_status == "missing"
            else "baseline_read_error",
            "message": message,
        }
    if np.asarray(candidate_data).shape != np.asarray(baseline_data).shape:
        return {
            **result,
            "status": "shape_mismatch",
            "message": (
                f"candidate shape {np.asarray(candidate_data).shape} != "
                f"baseline shape {np.asarray(baseline_data).shape}"
            ),
        }
    try:
        return {
            **result,
            **baseline_delta_metrics(candidate_data, baseline_data),
            "baseline_shape": _shape_of(baseline_data),
        }
    except (TypeError, ValueError) as exc:
        return {
            **result,
            "status": "non_numeric",
            "message": str(exc),
        }


def compare_candidate_variable(
    reference: netCDF4.Dataset,
    candidate: netCDF4.Dataset,
    baseline: netCDF4.Dataset | None,
    variable: str,
    *,
    baseline_is_candidate: bool = False,
    boundary_bands: Iterable[int] = DEFAULT_BOUNDARY_BANDS,
    time_index: int = -1,
) -> dict[str, Any]:
    reference_status, reference_data, reference_message = _read_optional_variable(
        reference,
        variable,
        time_index=time_index,
    )
    candidate_status, candidate_data, candidate_message = _read_optional_variable(
        candidate,
        variable,
        time_index=time_index,
    )

    result: dict[str, Any] = {
        "variable": variable,
        "status": "not_computed",
        "reference_present": reference_status == "available",
        "candidate_present": candidate_status == "available",
        "reference_shape": _shape_of(reference_data),
        "candidate_shape": _shape_of(candidate_data),
        "global": None,
        "reference_metrics": None,
        "boundary_bands": {
            "status": "not_computed",
            "bands": [],
        },
        "vertical_levels": {
            "status": "not_computed",
            "levels": [],
        },
        "delta_vs_baseline": _delta_vs_baseline(
            candidate_data,
            baseline,
            variable,
            candidate_status=candidate_status,
            baseline_is_candidate=baseline_is_candidate,
            time_index=time_index,
        ),
    }
    if reference_status != "available" or candidate_status != "available":
        missing = []
        if reference_status != "available":
            missing.append("reference")
        if candidate_status != "available":
            missing.append("candidate")
        messages = [
            item
            for item in (reference_message, candidate_message)
            if item is not None and item != "variable missing"
        ]
        return {
            **result,
            "status": "missing_variable"
            if "missing" in {reference_status, candidate_status}
            else "read_error",
            "message": "; ".join(messages)
            if messages
            else f"variable missing from {' and '.join(missing)}",
        }

    assert reference_data is not None
    assert candidate_data is not None
    if np.asarray(reference_data).shape != np.asarray(candidate_data).shape:
        return {
            **result,
            "status": "shape_mismatch",
            "message": (
                f"reference shape {np.asarray(reference_data).shape} != "
                f"candidate shape {np.asarray(candidate_data).shape}"
            ),
        }

    try:
        global_metrics = error_metrics(reference_data, candidate_data)
        boundary = boundary_band_statistics(
            reference_data,
            candidate_data,
            bands=boundary_bands,
        )
        vertical = vertical_level_statistics(reference_data, candidate_data)
    except (TypeError, ValueError) as exc:
        return {
            **result,
            "status": "non_numeric",
            "message": f"variable cannot be compared numerically: {exc}",
        }

    status = "computed" if global_metrics["valid_count"] > 0 else "no_valid_samples"
    return {
        **result,
        "status": status,
        "global": global_metrics,
        "reference_metrics": global_metrics,
        "boundary_bands": boundary,
        "vertical_levels": vertical,
        "message": None if status == "computed" else "no finite sample pairs available",
    }


def _normalize_candidates(
    candidates: Mapping[str, Path] | Iterable[tuple[str, Path | str]],
) -> list[tuple[str, Path]]:
    if isinstance(candidates, Mapping):
        items = list(candidates.items())
    else:
        items = list(candidates)
    if not items:
        raise ValueError("at least one --candidate LABEL=PATH is required")

    normalized: list[tuple[str, Path]] = []
    seen: set[str] = set()
    for label, path in items:
        label = str(label).strip()
        if not label:
            raise ValueError("candidate label cannot be empty")
        if label in seen:
            raise ValueError(f"duplicate candidate label: {label}")
        seen.add(label)
        normalized.append((label, Path(path)))
    return normalized


def audit_wind_subcycling(
    reference_wrfout: Path,
    candidates: Mapping[str, Path] | Iterable[tuple[str, Path | str]],
    *,
    baseline_wrfout: Path | None = None,
    variables: Iterable[str] = DEFAULT_VARIABLES,
    boundary_bands: Iterable[int] = DEFAULT_BOUNDARY_BANDS,
    time_index: int = -1,
) -> dict[str, Any]:
    candidate_paths = _normalize_candidates(candidates)
    variables = tuple(str(variable).strip() for variable in variables if str(variable).strip())
    boundary_bands = tuple(int(band) for band in boundary_bands)

    with ExitStack() as stack:
        reference = stack.enter_context(netCDF4.Dataset(reference_wrfout))
        baseline = (
            stack.enter_context(netCDF4.Dataset(baseline_wrfout))
            if baseline_wrfout is not None
            else None
        )
        candidate_reports = []
        for label, path in candidate_paths:
            with netCDF4.Dataset(path) as candidate:
                variable_reports = [
                    compare_candidate_variable(
                        reference,
                        candidate,
                        baseline,
                        variable,
                        baseline_is_candidate=_same_existing_path(
                            path,
                            baseline_wrfout,
                        ),
                        boundary_bands=boundary_bands,
                        time_index=time_index,
                    )
                    for variable in variables
                ]
                candidate_reports.append(
                    {
                        "label": label,
                        "path": str(path),
                        "metadata": _candidate_metadata(candidate),
                        "variables": variable_reports,
                    }
                )

    computed = sum(
        1
        for candidate in candidate_reports
        for variable in candidate["variables"]
        if variable["status"] == "computed"
    )
    total = len(candidate_reports) * len(variables)
    return {
        "status": "computed" if computed == total else "computed_with_flags",
        "diagnostic_only": True,
        "advances_00_20": False,
        "uses_reference_end_truth": False,
        "reference_wrfout": str(reference_wrfout),
        "baseline_wrfout": None if baseline_wrfout is None else str(baseline_wrfout),
        "variables": list(variables),
        "boundary_bands": list(boundary_bands),
        "time_index": int(time_index),
        "candidates": candidate_reports,
        "summary": {
            "candidate_count": len(candidate_reports),
            "variable_count": len(variables),
            "candidate_variable_count": total,
            "computed": computed,
            "flagged": total - computed,
        },
    }


def report_to_dict(report: dict[str, Any]) -> dict[str, Any]:
    return _json_value(report)


def report_to_json(report: dict[str, Any], *, pretty: bool = False) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


def parse_labeled_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"candidate must use LABEL=PATH: {value}")
    label, raw_path = value.split("=", 1)
    label = label.strip()
    raw_path = raw_path.strip()
    if not label or not raw_path:
        raise ValueError(f"candidate must use LABEL=PATH: {value}")
    return label, Path(raw_path)


def parse_csv(value: str, *, name: str) -> tuple[str, ...]:
    items = tuple(item.strip() for item in value.split(",") if item.strip())
    if not items:
        raise ValueError(f"{name} cannot be empty")
    return items


def parse_boundary_bands(value: str) -> tuple[int, ...]:
    raw_items = parse_csv(value, name="boundary bands")
    bands = []
    for raw in raw_items:
        try:
            band = int(raw)
        except ValueError as exc:
            raise ValueError(f"boundary band must be an integer: {raw}") from exc
        if band < 0:
            raise ValueError(f"boundary band must be non-negative: {raw}")
        bands.append(band)
    return tuple(bands)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-wrfout",
        type=Path,
        required=True,
        help="WRF reference wrfout used as the error target",
    )
    parser.add_argument(
        "--baseline-wrfout",
        type=Path,
        help="Optional baseline candidate wrfout used for delta comparison",
    )
    parser.add_argument(
        "--candidate",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="Candidate wrfout to audit; repeat for multiple candidates",
    )
    parser.add_argument(
        "--variables",
        default=",".join(DEFAULT_VARIABLES),
        help="Comma-separated variables to audit; default: U,V",
    )
    parser.add_argument(
        "--boundary-bands",
        default=",".join(str(value) for value in DEFAULT_BOUNDARY_BANDS),
        help="Comma-separated horizontal boundary bands; default: 0,1,2,5,10,20",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time index used for Time-leading variables; default: -1",
    )
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON")
    parser.add_argument("--output", type=Path, help="Optional JSON output path")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        variables = parse_csv(args.variables, name="variables")
        boundary_bands = parse_boundary_bands(args.boundary_bands)
        candidates = [parse_labeled_path(value) for value in args.candidate]
        report = audit_wind_subcycling(
            args.reference_wrfout,
            candidates,
            baseline_wrfout=args.baseline_wrfout,
            variables=variables,
            boundary_bands=boundary_bands,
            time_index=args.time_index,
        )
    except ValueError as exc:
        parser.error(str(exc))

    payload = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
