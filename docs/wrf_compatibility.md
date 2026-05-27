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

## KROSA Nesting And FDDA Constants

The v1 nesting interface is fixed to the current KROSA namelist subset:

- d01/d02 horizontal spacing: `10000 m` / `2000 m`;
- WRF namelist extents: d01 `e_we/e_sn = 266/430`, d02 `211/211`;
- NetCDF mass-grid extents: d01 `265 x 429`, d02 `210 x 210`;
- d02 parent start: `i_parent_start = 114`, `j_parent_start = 96`;
- `parent_grid_ratio = 5`, `parent_time_step_ratio = 5`;
- moving nest settings: `vortex_interval = 15,15`,
  `max_vortex_speed = 30,30`, `corral_dist = 10,30`,
  `track_level = 70000 Pa`, `time_to_move = 0,0`.

Spectral nudging is fixed to the KROSA `wrffdda_d01` contract:

- `grid_fdda = 2`;
- `gfdda_inname = "wrffdda_d<domain>"`;
- `gfdda_interval_m = 360`;
- `guv = 0.0003`;
- `xwavenum = 2`, `ywavenum = 4`;
- old/new fields: `U_NDG_OLD/U_NDG_NEW`, `V_NDG_OLD/V_NDG_NEW`,
  `T_NDG_OLD/T_NDG_NEW`, `Q_NDG_OLD/Q_NDG_NEW`, `PH_NDG_OLD/PH_NDG_NEW`,
  and `MU_NDG_OLD/MU_NDG_NEW`.

The current C++ nest module validates these constants and exchange contracts
only. Parent-child interpolation and two-way feedback return explicit
`not_implemented` statuses until the numerical kernels are added.

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
- C++ state-layout loading lives in `tywrf::io::load_wrf_state`.
- `tywrf::io::derive_grid_from_wrf_file` derives `Grid` from WRF mass and
  staggered dimensions:
  - `west_east`, `south_north`, `bottom_top`
  - `west_east_stag = west_east + 1`
  - `south_north_stag = south_north + 1`
  - `bottom_top_stag = bottom_top + 1`
- The loader reads WRF `Time,k,j,i` variables into the canonical TyWRF layout
  `((j*nz)+k)*nx+i`, offset into the field's active region when halos are
  present. The current tested baseline covers representative WRF-shaped fields
  `U`, `T`, and `MU`; the API has dispatch entries for the remaining v1 core
  state fields and reports missing/type/shape errors through
  `tywrf::io::WrfStateIoError`.
- `tywrf::io::write_wrf_state` creates a minimal WRF-shaped NetCDF file for
  State-backed fields. It defines the WRF mass/staggered dimensions and writes
  selected variables from canonical TyWRF layout back to WRF `Time,k,j,i`
  ordering. With no explicit selection it writes all currently State-backed
  core fields returned by `tywrf::io::wrf_state_writable_field_names`.
- The writer can now combine State-backed variables with selected metadata and
  static fields copied from a WRF template file. Set
  `WrfStateWriteOptions::template_path` to a compatible wrfout/wrfinput-style
  file, then request `Times`, `XLAT`, `XLONG`, or `HGT` alongside State-backed
  variables. `Times` can be overwritten with
  `WrfStateWriteOptions::times_value` so cycle-end candidates do not retain the
  template timestamp. `template_time_index` controls which template record is
  copied for time-dependent static fields.
- `Times`, `XLAT`, `XLONG`, and `HGT` are still not generated from `State`;
  `tywrf::io::wrf_state_writer_missing_field_names` reports them as fields that
  require a template source. Requesting them without
  `WrfStateWriteOptions::template_path` fails with `WrfStateIoError` instead of
  silently producing incomplete compatibility output.
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
