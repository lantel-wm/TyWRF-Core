#!/usr/bin/env python
"""Derive a sea-level pressure diagnostic from WRF-compatible core fields."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import netCDF4
import numpy as np


REQUIRED_VARIABLES = ("P", "PB", "T", "QVAPOR", "PH", "PHB", "PSFC")
SLP_VARIABLE_NAME = "SLP"
DERIVED_SLP_METHOD = "hypsometric_lowest_mass_level_virtual_temperature"
P0_PA = 100000.0
RD_J_KG_K = 287.05
CP_J_KG_K = 1004.0
KAPPA = RD_J_KG_K / CP_J_KG_K
GRAVITY_M_S2 = 9.80665


class DeriveMSLPError(ValueError):
    """Raised when the source file cannot provide derived MSLP."""


@dataclass(frozen=True)
class DerivedMSLPSummary:
    source: str
    destination: str
    variable: str
    units: str
    method: str
    required_variables: list[str]
    shape: tuple[int, ...]
    finite_count: int
    total_count: int
    minimum_hpa: float
    maximum_hpa: float
    replaced_source_variable: bool


def _copy_attributes(source, destination) -> None:
    attrs = {name: source.getncattr(name) for name in source.ncattrs() if name != "_FillValue"}
    destination.setncatts(attrs)


def _read_numeric(dataset: netCDF4.Dataset, variable_name: str) -> np.ndarray:
    variable = dataset.variables[variable_name]
    try:
        data = np.ma.asarray(variable[:], dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise DeriveMSLPError(f"variable {variable_name} is not numeric: {exc}") from exc
    return np.asarray(np.ma.filled(data, np.nan), dtype=np.float64)


def _require_variables(dataset: netCDF4.Dataset) -> None:
    missing = [name for name in REQUIRED_VARIABLES if name not in dataset.variables]
    if missing:
        raise DeriveMSLPError(
            "missing required variable(s) for derived SLP: " + ", ".join(missing)
        )


def _require_safe_output_variable(variable_name: str) -> None:
    source_collisions = {
        required.lower() for required in REQUIRED_VARIABLES
    }
    if variable_name.lower() in source_collisions:
        raise DeriveMSLPError(
            f"derived SLP output variable {variable_name!r} would replace a required "
            "source field; choose a diagnostic name such as SLP"
        )


def _bottom_level(
    dataset: netCDF4.Dataset,
    variable_name: str,
    vertical_dimension: str,
) -> np.ndarray:
    variable = dataset.variables[variable_name]
    if vertical_dimension not in variable.dimensions:
        raise DeriveMSLPError(
            f"variable {variable_name} must include vertical dimension {vertical_dimension}; "
            f"got dimensions {variable.dimensions}"
        )

    vertical_axis = variable.dimensions.index(vertical_dimension)
    if variable.shape[vertical_axis] < 1:
        raise DeriveMSLPError(f"variable {variable_name} has no bottom level")

    data = _read_numeric(dataset, variable_name)
    return np.take(data, indices=0, axis=vertical_axis)


def _require_same_shape(arrays: dict[str, np.ndarray]) -> tuple[int, ...]:
    shapes = {name: values.shape for name, values in arrays.items()}
    unique_shapes = set(shapes.values())
    if len(unique_shapes) != 1:
        detail = ", ".join(f"{name}={shape}" for name, shape in sorted(shapes.items()))
        raise DeriveMSLPError(f"derived SLP source shapes do not match: {detail}")
    return next(iter(unique_shapes))


def derive_mslp_hpa(dataset: netCDF4.Dataset) -> np.ndarray:
    """Return derived sea-level pressure in hPa for every PSFC grid point.

    The approximation uses the lowest WRF mass level to estimate virtual
    temperature, the bottom geopotential staggered level for terrain height,
    then hydrostatically reduces surface pressure to sea level:

    ``SLP = PSFC * exp(g * z_surface / (Rd * Tv_lowest))``.
    """

    _require_variables(dataset)

    psfc_pa = _read_numeric(dataset, "PSFC")
    pressure_pa = _bottom_level(dataset, "P", "bottom_top") + _bottom_level(
        dataset, "PB", "bottom_top"
    )
    theta_k = _bottom_level(dataset, "T", "bottom_top") + 300.0
    qvapor_kg_kg = _bottom_level(dataset, "QVAPOR", "bottom_top")
    geopotential_m2_s2 = _bottom_level(dataset, "PH", "bottom_top_stag") + _bottom_level(
        dataset, "PHB", "bottom_top_stag"
    )
    surface_height_m = geopotential_m2_s2 / GRAVITY_M_S2

    shape = _require_same_shape(
        {
            "PSFC": psfc_pa,
            "P+PB bottom": pressure_pa,
            "T bottom": theta_k,
            "QVAPOR bottom": qvapor_kg_kg,
            "PH+PHB surface": surface_height_m,
        }
    )

    qvapor_nonnegative = np.maximum(qvapor_kg_kg, 0.0)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        temperature_k = theta_k * np.power(pressure_pa / P0_PA, KAPPA)
        virtual_temperature_k = temperature_k * (1.0 + 0.61 * qvapor_nonnegative)
        valid = (
            np.isfinite(psfc_pa)
            & np.isfinite(pressure_pa)
            & np.isfinite(theta_k)
            & np.isfinite(qvapor_kg_kg)
            & np.isfinite(surface_height_m)
            & np.isfinite(virtual_temperature_k)
            & (psfc_pa > 0.0)
            & (pressure_pa > 0.0)
            & (theta_k > 0.0)
            & (virtual_temperature_k > 0.0)
        )
        slp_hpa = (
            psfc_pa
            * np.exp(GRAVITY_M_S2 * surface_height_m / (RD_J_KG_K * virtual_temperature_k))
        ) / 100.0

    slp_hpa = np.where(valid, slp_hpa, np.nan).astype(np.float64, copy=False)
    finite_count = int(np.count_nonzero(np.isfinite(slp_hpa)))
    if finite_count == 0:
        raise DeriveMSLPError(f"derived SLP has no finite samples for shape {shape}")
    return slp_hpa


def _copy_dataset_with_derived_slp(
    source: netCDF4.Dataset,
    destination: netCDF4.Dataset,
    slp_hpa: np.ndarray,
    *,
    variable_name: str,
) -> None:
    _copy_attributes(source, destination)
    destination.setncattr("TYWRF_DERIVED_SLP", "true")
    destination.setncattr("TYWRF_DERIVED_SLP_DIAGNOSTIC", "true")
    destination.setncattr("TYWRF_DERIVED_SLP_IS_PSFC_PROXY", "false")
    destination.setncattr("TYWRF_DERIVED_SLP_METHOD", DERIVED_SLP_METHOD)
    destination.setncattr("TYWRF_DERIVED_SLP_REQUIRED_VARIABLES", ",".join(REQUIRED_VARIABLES))

    for name, dimension in source.dimensions.items():
        size = None if dimension.isunlimited() else len(dimension)
        destination.createDimension(name, size)

    for name, source_variable in source.variables.items():
        if name == variable_name:
            continue
        fill_value = getattr(source_variable, "_FillValue", None)
        destination_variable = destination.createVariable(
            name,
            source_variable.datatype,
            source_variable.dimensions,
            fill_value=fill_value,
        )
        _copy_attributes(source_variable, destination_variable)
        destination_variable[:] = source_variable[:]

    psfc_variable = source.variables["PSFC"]
    slp_variable = destination.createVariable(
        variable_name,
        "f4",
        psfc_variable.dimensions,
        fill_value=np.float32(np.nan),
    )
    slp_variable.long_name = (
        "sea level pressure diagnostic derived by TyWRF-Core validation tooling"
    )
    slp_variable.description = (
        "Hypsometric reduction from PSFC using lowest mass-level virtual temperature "
        "from P, PB, T, QVAPOR and surface geopotential height from PH, PHB; "
        "this is not a copied WRF end-state SLP field and not a PSFC proxy renamed as SLP"
    )
    slp_variable.units = "hPa"
    slp_variable.FieldType = 104
    slp_variable.MemoryOrder = getattr(psfc_variable, "MemoryOrder", "XY")
    slp_variable.coordinates = getattr(psfc_variable, "coordinates", "XLONG XLAT")
    slp_variable.derived_from = ",".join(REQUIRED_VARIABLES)
    slp_variable.method = DERIVED_SLP_METHOD
    slp_variable.TYWRF_DERIVED_DIAGNOSTIC = "true"
    slp_variable.TYWRF_IS_PSFC_PROXY = "false"
    slp_variable.TYWRF_SOURCE_KIND = "hypsometric_diagnostic_from_core_fields"
    slp_variable[:] = slp_hpa.astype(np.float32)


def write_derived_mslp(
    source_path: Path,
    destination_path: Path,
    *,
    variable_name: str = SLP_VARIABLE_NAME,
) -> DerivedMSLPSummary:
    """Copy a WRF-compatible file and add a derived ``SLP`` diagnostic."""

    if source_path.resolve() == destination_path.resolve():
        raise DeriveMSLPError("source and destination must be different paths")
    _require_safe_output_variable(variable_name)

    destination_path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(source_path) as source:
        slp_hpa = derive_mslp_hpa(source)
        replaced = variable_name in source.variables
        finite_values = slp_hpa[np.isfinite(slp_hpa)]

        with netCDF4.Dataset(destination_path, "w", format=source.data_model) as destination:
            _copy_dataset_with_derived_slp(
                source,
                destination,
                slp_hpa,
                variable_name=variable_name,
            )

    return DerivedMSLPSummary(
        source=str(source_path),
        destination=str(destination_path),
        variable=variable_name,
        units="hPa",
        method="hypsometric_lowest_mass_level_virtual_temperature",
        required_variables=list(REQUIRED_VARIABLES),
        shape=tuple(int(size) for size in slp_hpa.shape),
        finite_count=int(finite_values.size),
        total_count=int(slp_hpa.size),
        minimum_hpa=float(np.min(finite_values)),
        maximum_hpa=float(np.max(finite_values)),
        replaced_source_variable=replaced,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Source WRF-compatible NetCDF file")
    parser.add_argument("destination", type=Path, help="Destination NetCDF file with derived SLP")
    parser.add_argument(
        "--variable-name",
        default=SLP_VARIABLE_NAME,
        help="Output pressure diagnostic variable name; default: SLP",
    )
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        summary = write_derived_mslp(
            args.source,
            args.destination,
            variable_name=args.variable_name,
        )
    except (OSError, DeriveMSLPError) as exc:
        parser.error(str(exc))

    print(json.dumps(asdict(summary), indent=2 if args.pretty else None))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
