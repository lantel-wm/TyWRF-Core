#!/usr/bin/env python
"""Decompose U/V implied wind tendency residuals for diagnostic audits.

The moving-nest region split uses child-grid cells.  Pass ``--child-delta I J``
when the nest displacement is already expressed on the child output grid.  When
``--from-parent-start`` and ``--to-parent-start`` are used, the reported delta
is the unscaled parent-start index difference; use ``--child-delta`` for
parent-ratio-scaled D90/D91 region widths.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable
from contextlib import ExitStack
import json
import math
from pathlib import Path
from typing import Any

import netCDF4
import numpy as np


DEFAULT_VARIABLES = ("U", "V")
DEFAULT_BOUNDARY_BANDS = (0, 1, 2, 5, 10, 20, 40)
IMPLIED_INCREMENT_SECONDS = 600.0


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
    time_index: int = -1,
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


def _shape_of(data: Any | None) -> list[int] | None:
    if data is None:
        return None
    return [int(value) for value in np.asarray(data).shape]


def error_metrics(
    reference_increment: Any,
    candidate_increment: Any,
    *,
    selection_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    ref_values, ref_finite = _finite_values(reference_increment)
    cand_values, cand_finite = _finite_values(candidate_increment)
    if ref_values.shape != cand_values.shape:
        raise ValueError(f"shape mismatch: {ref_values.shape} != {cand_values.shape}")

    if selection_mask is None:
        selected = np.ones(ref_values.shape, dtype=bool)
    else:
        selected = np.asarray(selection_mask, dtype=bool)
        if selected.shape != ref_values.shape:
            selected = np.broadcast_to(selected, ref_values.shape)

    valid = selected & ref_finite & cand_finite
    selected_count = int(np.count_nonzero(selected))
    valid_count = int(np.count_nonzero(valid))
    if valid_count == 0:
        return {
            "rmse": None,
            "normalized_rmse": None,
            "max_abs_error": None,
            "valid_count": 0,
            "count": 0,
            "selected_count": selected_count,
            "sum_squared_residual": None,
        }

    ref_sample = ref_values[valid]
    residual = cand_values[valid] - ref_sample
    rmse = float(math.sqrt(float(np.mean(residual * residual))))
    scale = float(math.sqrt(float(np.mean(ref_sample * ref_sample))))
    if scale == 0.0:
        normalized = 0.0 if rmse == 0.0 else None
    else:
        normalized = float(rmse / scale)
    return {
        "rmse": rmse,
        "normalized_rmse": normalized,
        "max_abs_error": float(np.max(np.abs(residual))),
        "valid_count": valid_count,
        "count": valid_count,
        "selected_count": selected_count,
        "sum_squared_residual": float(np.sum(residual * residual)),
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
    reference_increment: Any,
    candidate_increment: Any,
    *,
    bands: Iterable[int] = DEFAULT_BOUNDARY_BANDS,
) -> dict[str, Any]:
    ref_values = np.asarray(reference_increment)
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
                    reference_increment,
                    candidate_increment,
                    selection_mask=_broadcast_horizontal_mask(boundary_2d, shape),
                ),
                "interior": error_metrics(
                    reference_increment,
                    candidate_increment,
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


def vertical_level_statistics(reference_increment: Any, candidate_increment: Any) -> dict[str, Any]:
    ref_values = np.asarray(reference_increment)
    shape = tuple(int(value) for value in ref_values.shape)
    if len(shape) < 3:
        return {
            "status": "no_vertical_dimension",
            "k_axis": None,
            "levels": [
                {
                    "k": None,
                    "metrics": error_metrics(reference_increment, candidate_increment),
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
                    np.take(reference_increment, k, axis=k_axis),
                    np.take(candidate_increment, k, axis=k_axis),
                ),
            }
        )
    return {
        "status": "available",
        "k_axis": int(k_axis),
        "level_count": int(shape[k_axis]),
        "levels": levels,
    }


def _clamped_width(delta: int, limit: int) -> int:
    return min(abs(int(delta)), int(limit))


def _exposed_mask_1d(delta: int, limit: int) -> np.ndarray:
    width = _clamped_width(delta, limit)
    mask = np.zeros(limit, dtype=bool)
    if width == 0:
        return mask
    if delta >= 0:
        mask[limit - width :] = True
    else:
        mask[:width] = True
    return mask


def region_masks_2d(
    horizontal_shape: tuple[int, int],
    child_delta: tuple[int, int],
    *,
    core_boundary_band: int,
) -> dict[str, np.ndarray]:
    ny, nx = (int(horizontal_shape[0]), int(horizontal_shape[1]))
    delta_i, delta_j = (int(child_delta[0]), int(child_delta[1]))
    exposed_i = _exposed_mask_1d(delta_i, nx)[None, :]
    exposed_j = _exposed_mask_1d(delta_j, ny)[:, None]
    exposed_i_2d = np.broadcast_to(exposed_i, (ny, nx))
    exposed_j_2d = np.broadcast_to(exposed_j, (ny, nx))

    overlap = ~(exposed_i_2d | exposed_j_2d)
    corner = exposed_i_2d & exposed_j_2d
    distance = _horizontal_boundary_distance((ny, nx))
    assert distance is not None
    core = overlap & (distance > int(core_boundary_band))

    i_name = "east_exposed" if delta_i >= 0 else "west_exposed"
    j_name = "north_exposed" if delta_j >= 0 else "south_exposed"
    return {
        "overlap": overlap,
        "core": core,
        i_name: exposed_i_2d & ~exposed_j_2d,
        j_name: exposed_j_2d & ~exposed_i_2d,
        "corner": corner,
    }


def region_decomposition(
    reference_increment: Any,
    candidate_increment: Any,
    *,
    child_delta: tuple[int, int] | None,
    child_delta_source: str | None,
    core_boundary_band: int,
) -> dict[str, Any]:
    ref_values = np.asarray(reference_increment)
    shape = tuple(int(value) for value in ref_values.shape)
    if child_delta is None:
        return {
            "status": "not_requested",
            "regions": [],
            "regions_by_name": {},
        }
    if len(shape) < 2:
        return {
            "status": "no_horizontal_dimension",
            "regions": [],
            "regions_by_name": {},
        }

    global_metrics = error_metrics(reference_increment, candidate_increment)
    global_sse = global_metrics["sum_squared_residual"]
    masks = region_masks_2d(
        (shape[-2], shape[-1]),
        child_delta,
        core_boundary_band=core_boundary_band,
    )

    entries = []
    by_name = {}
    for name, mask_2d in masks.items():
        mask = _broadcast_horizontal_mask(mask_2d, shape)
        metrics = error_metrics(
            reference_increment,
            candidate_increment,
            selection_mask=mask,
        )
        sse = metrics["sum_squared_residual"]
        if global_sse in (None, 0.0) or sse is None:
            fraction = None
        else:
            fraction = float(sse / global_sse)
        entry = {
            "name": name,
            "metrics": metrics,
            "sum_squared_residual_fraction_global": fraction,
        }
        entries.append(entry)
        by_name[name] = entry

    delta_i, delta_j = child_delta
    return {
        "status": "available",
        "horizontal_shape": [shape[-2], shape[-1]],
        "child_delta": [int(delta_i), int(delta_j)],
        "child_delta_units": "child_grid_cells"
        if child_delta_source == "child_delta"
        else "parent_start_index_cells_unscaled",
        "child_delta_source": child_delta_source,
        "positive_i_exposed_side": "east",
        "positive_j_exposed_side": "north",
        "core_boundary_band": int(core_boundary_band),
        "region_set_semantics": (
            "overlap, exposed side strips, and corner are exclusive; core is a "
            "subset of overlap with boundary distance greater than core_boundary_band"
        ),
        "regions": entries,
        "regions_by_name": by_name,
    }


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


def _increment(
    start: Any,
    end: Any,
) -> np.ndarray:
    start_values = np.asarray(np.ma.asarray(start, dtype=np.float64).filled(np.nan))
    end_values = np.asarray(np.ma.asarray(end, dtype=np.float64).filled(np.nan))
    if start_values.shape != end_values.shape:
        raise ValueError(f"start shape {start_values.shape} != end shape {end_values.shape}")
    return end_values - start_values


def _increment_stats(increment: Any) -> dict[str, Any]:
    values, finite = _finite_values(increment)
    valid_count = int(np.count_nonzero(finite))
    if valid_count == 0:
        return {
            "valid_count": 0,
            "rms": None,
            "mean": None,
            "mean_abs": None,
            "max_abs": None,
        }
    sample = values[finite]
    return {
        "valid_count": valid_count,
        "rms": float(math.sqrt(float(np.mean(sample * sample)))),
        "mean": float(np.mean(sample)),
        "mean_abs": float(np.mean(np.abs(sample))),
        "max_abs": float(np.max(np.abs(sample))),
    }


def compare_variable_tendency_decomposition(
    reference_start: netCDF4.Dataset,
    reference_end: netCDF4.Dataset,
    candidate_start: netCDF4.Dataset,
    candidate_end: netCDF4.Dataset,
    variable: str,
    *,
    boundary_bands: Iterable[int] = DEFAULT_BOUNDARY_BANDS,
    child_delta: tuple[int, int] | None = None,
    child_delta_source: str | None = None,
    time_index: int = -1,
) -> dict[str, Any]:
    sources = {
        "reference_start": _read_optional_variable(
            reference_start,
            variable,
            time_index=time_index,
        ),
        "reference_end": _read_optional_variable(
            reference_end,
            variable,
            time_index=time_index,
        ),
        "candidate_start": _read_optional_variable(
            candidate_start,
            variable,
            time_index=time_index,
        ),
        "candidate_end": _read_optional_variable(
            candidate_end,
            variable,
            time_index=time_index,
        ),
    }
    shapes = {name: _shape_of(data) for name, (_, data, _) in sources.items()}
    missing = [
        name
        for name, (status, _, _) in sources.items()
        if status != "available"
    ]
    result: dict[str, Any] = {
        "variable": variable,
        "status": "not_computed",
        "diagnostic_only": True,
        "uses_reference_end_truth": True,
        "advances_00_20": False,
        "increment_seconds": IMPLIED_INCREMENT_SECONDS,
        "residual_definition": (
            "(candidate_end - candidate_start) - "
            "(reference_end - reference_start)"
        ),
        "source_shapes": shapes,
        "global": None,
        "reference_increment_stats": None,
        "candidate_increment_stats": None,
        "boundary_bands": {
            "status": "not_computed",
            "bands": [],
        },
        "vertical_levels": {
            "status": "not_computed",
            "levels": [],
        },
        "regions": {
            "status": "not_computed",
            "regions": [],
            "regions_by_name": {},
        },
    }
    if missing:
        messages = [
            f"{name}: {message}"
            for name, (_, _, message) in sources.items()
            if message is not None and name in missing
        ]
        return {
            **result,
            "status": "missing_variable"
            if any(sources[name][0] == "missing" for name in missing)
            else "read_error",
            "message": "; ".join(messages) if messages else "variable missing",
        }

    data = {name: value for name, (_, value, _) in sources.items()}
    assert all(value is not None for value in data.values())
    try:
        reference_increment = _increment(data["reference_start"], data["reference_end"])
        candidate_increment = _increment(data["candidate_start"], data["candidate_end"])
    except (TypeError, ValueError) as exc:
        return {
            **result,
            "status": "shape_mismatch" if "shape" in str(exc) else "non_numeric",
            "message": str(exc),
        }

    if reference_increment.shape != candidate_increment.shape:
        return {
            **result,
            "status": "shape_mismatch",
            "reference_increment_shape": _shape_of(reference_increment),
            "candidate_increment_shape": _shape_of(candidate_increment),
            "message": (
                f"reference increment shape {reference_increment.shape} != "
                f"candidate increment shape {candidate_increment.shape}"
            ),
        }

    try:
        global_metrics = error_metrics(reference_increment, candidate_increment)
        boundary = boundary_band_statistics(
            reference_increment,
            candidate_increment,
            bands=boundary_bands,
        )
        vertical = vertical_level_statistics(reference_increment, candidate_increment)
        regions = region_decomposition(
            reference_increment,
            candidate_increment,
            child_delta=child_delta,
            child_delta_source=child_delta_source,
            core_boundary_band=max(int(band) for band in boundary_bands),
        )
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
        "reference_increment_shape": _shape_of(reference_increment),
        "candidate_increment_shape": _shape_of(candidate_increment),
        "global": global_metrics,
        "reference_increment_stats": _increment_stats(reference_increment),
        "candidate_increment_stats": _increment_stats(candidate_increment),
        "boundary_bands": boundary,
        "vertical_levels": vertical,
        "regions": regions,
        "message": None if status == "computed" else "no finite sample pairs available",
    }


def resolve_child_delta(
    *,
    child_delta: tuple[int, int] | None = None,
    from_parent_start: tuple[int, int] | None = None,
    to_parent_start: tuple[int, int] | None = None,
) -> tuple[tuple[int, int] | None, str | None]:
    if child_delta is not None and (
        from_parent_start is not None or to_parent_start is not None
    ):
        raise ValueError("use --child-delta or --from-parent-start/--to-parent-start, not both")
    if child_delta is not None:
        return (int(child_delta[0]), int(child_delta[1])), "child_delta"
    if from_parent_start is None and to_parent_start is None:
        return None, None
    if from_parent_start is None or to_parent_start is None:
        raise ValueError("--from-parent-start and --to-parent-start must be provided together")
    return (
        int(to_parent_start[0]) - int(from_parent_start[0]),
        int(to_parent_start[1]) - int(from_parent_start[1]),
    ), "parent_start_difference_unscaled"


def audit_wind_tendency_decomposition(
    reference_start_wrfout: Path,
    reference_end_wrfout: Path,
    candidate_start_wrfout: Path,
    candidate_end_wrfout: Path,
    *,
    variables: Iterable[str] = DEFAULT_VARIABLES,
    boundary_bands: Iterable[int] = DEFAULT_BOUNDARY_BANDS,
    child_delta: tuple[int, int] | None = None,
    child_delta_source: str | None = None,
    time_index: int = -1,
) -> dict[str, Any]:
    variables = tuple(str(variable).strip() for variable in variables if str(variable).strip())
    if not variables:
        raise ValueError("variables cannot be empty")
    boundary_bands = tuple(int(band) for band in boundary_bands)
    if not boundary_bands:
        raise ValueError("boundary bands cannot be empty")

    with ExitStack() as stack:
        reference_start = stack.enter_context(netCDF4.Dataset(reference_start_wrfout))
        reference_end = stack.enter_context(netCDF4.Dataset(reference_end_wrfout))
        candidate_start = stack.enter_context(netCDF4.Dataset(candidate_start_wrfout))
        candidate_end = stack.enter_context(netCDF4.Dataset(candidate_end_wrfout))
        results = [
            compare_variable_tendency_decomposition(
                reference_start,
                reference_end,
                candidate_start,
                candidate_end,
                variable,
                boundary_bands=boundary_bands,
                child_delta=child_delta,
                child_delta_source=child_delta_source,
                time_index=time_index,
            )
            for variable in variables
        ]

    computed = sum(1 for item in results if item["status"] == "computed")
    return {
        "status": "computed" if computed == len(results) else "computed_with_flags",
        "diagnostic_only": True,
        "uses_reference_end_truth": True,
        "advances_00_20": False,
        "reference_start_wrfout": str(reference_start_wrfout),
        "reference_end_wrfout": str(reference_end_wrfout),
        "candidate_start_wrfout": str(candidate_start_wrfout),
        "candidate_end_wrfout": str(candidate_end_wrfout),
        "variables": list(variables),
        "increment_seconds": IMPLIED_INCREMENT_SECONDS,
        "boundary_bands": list(boundary_bands),
        "child_delta": None if child_delta is None else [int(child_delta[0]), int(child_delta[1])],
        "child_delta_source": child_delta_source,
        "child_delta_semantics": (
            "child_delta is [i, j] in child-grid cells; from/to parent starts are "
            "reported as unscaled parent-start index differences"
        ),
        "time_index": int(time_index),
        "results": results,
        "summary": {
            "diagnostic_only": True,
            "uses_reference_end_truth": True,
            "advances_00_20": False,
            "variable_count": len(results),
            "computed": computed,
            "flagged": len(results) - computed,
        },
    }


def report_to_dict(report: dict[str, Any]) -> dict[str, Any]:
    return _json_value(report)


def report_to_json(report: dict[str, Any], *, pretty: bool = False) -> str:
    return json.dumps(report_to_dict(report), indent=2 if pretty else None, allow_nan=False)


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


def _pair(values: list[int] | None) -> tuple[int, int] | None:
    if values is None:
        return None
    return int(values[0]), int(values[1])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Read-only D92 audit: compare 600 s implied U/V increments and "
            "decompose candidate-minus-reference residuals."
        )
    )
    parser.add_argument("--reference-start-wrfout", required=True, type=Path)
    parser.add_argument("--reference-end-wrfout", required=True, type=Path)
    parser.add_argument("--candidate-start-wrfout", required=True, type=Path)
    parser.add_argument("--candidate-end-wrfout", required=True, type=Path)
    parser.add_argument(
        "--variables",
        default=",".join(DEFAULT_VARIABLES),
        help="Comma-separated variables to audit; default: U,V",
    )
    parser.add_argument(
        "--from-parent-start",
        nargs=2,
        type=int,
        metavar=("I", "J"),
        help="Original parent-start indices; paired with --to-parent-start.",
    )
    parser.add_argument(
        "--to-parent-start",
        nargs=2,
        type=int,
        metavar=("I", "J"),
        help="New parent-start indices; difference is unscaled in the report.",
    )
    parser.add_argument(
        "--child-delta",
        nargs=2,
        type=int,
        metavar=("I", "J"),
        help="Moving-nest displacement in child-grid cells for region decomposition.",
    )
    parser.add_argument(
        "--boundary-bands",
        default=",".join(str(value) for value in DEFAULT_BOUNDARY_BANDS),
        help="Comma-separated horizontal boundary bands; default: 0,1,2,5,10,20,40",
    )
    parser.add_argument("--output", type=Path, help="Write JSON report to this path")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        variables = parse_csv(args.variables, name="variables")
        boundary_bands = parse_boundary_bands(args.boundary_bands)
        child_delta, child_delta_source = resolve_child_delta(
            child_delta=_pair(args.child_delta),
            from_parent_start=_pair(args.from_parent_start),
            to_parent_start=_pair(args.to_parent_start),
        )
        report = audit_wind_tendency_decomposition(
            args.reference_start_wrfout,
            args.reference_end_wrfout,
            args.candidate_start_wrfout,
            args.candidate_end_wrfout,
            variables=variables,
            boundary_bands=boundary_bands,
            child_delta=child_delta,
            child_delta_source=child_delta_source,
        )
    except (IndexError, OSError, ValueError) as exc:
        parser.error(str(exc))

    text = report_to_json(report, pretty=args.pretty)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
