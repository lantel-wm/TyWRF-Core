# WRF Compatibility

TyWRF-Core v1 reads WRF standard inputs from the current PGWRF/KROSA case and
writes a reduced, WRF-compatible core `wrfout` set.

## Required Inputs

- `namelist.input`
- `wrfinput_d01`
- `wrfbdy_d01`
- `wrffdda_d01`

The v1 integrator reuses WPS and `real.exe`; it does not reimplement WRF
pre-processing.

## Minimum Output Variables

- Coordinates/static/time: `Times`, `XLAT`, `XLONG`, `HGT`
- Dynamics/core state: `U`, `V`, `W`, `PH`, `PHB`, `T`, `MU`, `MUB`, `P`, `PB`
- Moisture/hydrometeors: `QVAPOR`, `QCLOUD`, `QRAIN`, `QICE`, `QSNOW`,
  `QGRAUP`, `QNICE`, `QNRAIN`
- Near-surface diagnostics: `PSFC`, `U10`, `V10`, `T2`, `Q2`
- Precipitation: `RAINC`, `RAINNC`

## Reference Schema Snapshot

Sampled on 2026-05-27 from:

```text
/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF
```

`wrfout_d01_2025-07-26_00:00:00`:

- `west_east = 265`
- `south_north = 429`
- `bottom_top = 59`
- `bottom_top_stag = 60`
- `west_east_stag = 266`
- `south_north_stag = 430`

`wrfout_d02_2025-07-26_00:00:00`:

- `west_east = 210`
- `south_north = 210`
- `bottom_top = 59`
- `bottom_top_stag = 60`
- `west_east_stag = 211`
- `south_north_stag = 211`

`wrfbdy_d01`:

- `Time = 28`
- `bdy_width = 5`

`wrffdda_d01`:

- `Time = 28`
- nudging fields include U/V/T/PH/MU old/new fields;
- Q nudging fields are named `Q_NDG_OLD` and `Q_NDG_NEW`.

## Phase 2 Baseline Interfaces

- C++ NetCDF schema reading lives in `tywrf::io::read_netcdf_schema`.
- Current KROSA input validation is exposed through:
  - `validate_wrfinput_d01_schema`
  - `validate_wrfbdy_d01_schema`
  - `validate_wrffdda_d01_schema`
- The validators check the d01 KROSA dimension sizes, unlimited `Time`
  status, required variable presence, NetCDF primitive type, and variable
  dimension order. `wrfbdy_d01` validation includes U/V/W/PH/T/MU plus
  QVAPOR, QCLOUD, QRAIN, QICE, QSNOW, QGRAUP, QNICE, and QNRAIN boundary
  values and tendencies.
- Reduced core `wrfout` files can be generated for compatibility smoke tests
  with `tools/write_core_wrfout.py`.
- The reduced writer preserves the source NetCDF data model, all source
  dimensions, global attributes, selected variable attributes, and the minimum
  v1 core variables by default.
