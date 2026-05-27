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
and now stages selected moving-nest remap helpers. Full parent-child
interpolation and two-way feedback still return explicit `not_implemented`
statuses until the numerical kernels are added.

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
- C++ forcing-file reads live in `tywrf::io::KrosaForcingReader`. This is a
  minimum KROSA-only interface for `wrfbdy_d01` and `wrffdda_d01`, not a
  general WRF forcing abstraction. It reports the `Time` length, `Times`
  strings, required-variable presence, variable dimensions/shapes, and reads a
  selected float variable at one `Time` index into a contiguous
  `std::vector<float>`.
- The forcing reader keeps NetCDF I/O outside hot kernels. Boundary and
  spectral-nudging kernels should receive staged POD views later; this
  interface only provides metadata and staging buffers for the current
  KROSA-compatible integrator path.
- The current required KROSA forcing lists cover `Times`; `wrfbdy_d01` U/V/W/PH,
  T, MU, QVAPOR, QCLOUD, QRAIN, QICE, QSNOW, QGRAUP, QNICE, and QNRAIN boundary
  values/tendencies; and `wrffdda_d01` old/new U/V/T/Q/PH/MU nudging fields.
- Reduced core `wrfout` files can be generated for compatibility smoke tests
  with `tools/write_core_wrfout.py`.
- The reduced writer preserves the source NetCDF data model, all source
  dimensions, global attributes, selected variable attributes, and the minimum
  v1 core variables by default.

## KROSA d02 Diagnostic Remap Smoke

Agent R12 ran a real KROSA d02 smoke on 2026-05-27 for the C++
`tywrf_skeleton_cycle --diagnostic-remap-overlap` path. This was a
diagnostic-only moving-nest overlap exercise, not a normal validation pass.

The existing built executable was used:

```bash
build/tywrf_skeleton_cycle \
  --domain d02 \
  --state /home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF/wrfout_d02_2025-07-26_00:00:00 \
  --template /home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF/wrfout_d02_2025-07-26_00:00:00 \
  --output /tmp/tywrf_r12_krosa_d02_remap_overlap_wrfout_d02_2025-07-26_06:00:00 \
  --times 2025-07-26_06:00:00 \
  --diagnostic-remap-overlap \
  --from-parent-start 114,96 \
  --to-parent-start 125,105 \
  --pretty
```

Observed report metadata:

- `status = diagnostic_remap_overlap_generated`
- `candidate_kind = cpp_skeleton_remap_overlap_diagnostic`
- `domain = d02`, with `dx_m = 2000` and `dy_m = 2000`
- `times_source = --times:2025-07-26_06:00:00`
- `diagnostic_remap_overlap = true`
- `diagnostic_only = true`
- `gate_candidate = false`
- `integrator_output = false`
- `not_physical = true`
- move: `from_parent_start = 114,96` to `to_parent_start = 125,105`
- `remap_copied_field_count = 25`
- `remap_copied_point_count = 24468580`
- `needs_parent_fill = true`
- `exposed_cells_filled_by_parent = false`
- `unfilled_exposed_cells = nan_sentinel_pending_parent_fill`

Direct NetCDF metadata inspection of the generated candidate confirmed:

- `TYWRF_CANDIDATE_KIND = cpp_skeleton_remap_overlap_diagnostic`
- `TYWRF_DIAGNOSTIC_REMAP_OVERLAP = true`
- `TYWRF_DIAGNOSTIC_ONLY = true`
- `TYWRF_GATE_CANDIDATE = false`
- `TYWRF_NEEDS_PARENT_FILL = true`
- `TYWRF_REMAP_EXPOSED_CELLS_FILLED_BY_PARENT = false`
- `TYWRF_UNFILLED_EXPOSED_CELLS = nan_sentinel_pending_parent_fill`
- `TYWRF_REMAP_FROM_PARENT_START = 114,96`
- `TYWRF_REMAP_TO_PARENT_START = 125,105`
- `TYWRF_REMAP_COPIED_FIELD_COUNT = 25`
- `TYWRF_REMAP_COPIED_POINT_COUNT = 24468580`
- output `Times = 2025-07-26_06:00:00`

This smoke confirms that the C++ skeleton can read the real KROSA d02 00 h
state, write the reduced WRF-compatible d02 candidate at the 06 h timestamp,
and apply the diagnostic overlap copy for a meaningful moving-nest shift.
Newly exposed d02 cells remain intentionally marked with NaN until parent-fill
interpolation is implemented, so this output must not be used as a physical
6 h candidate or as a validation-gate pass.

## Moving-Nest Parent Fill And Pressure Refresh

TyWRF now stages the WRF-style moving-nest start sequence as:

1. overlap remap from the old child state;
2. direct parent-fill for WRF direct-fill fields in newly exposed cells;
3. post-`start_domain` time-level sync for exposed cells;
4. a required derived/recomputed perturbation-pressure refresh.

The direct parent-fill field set intentionally excludes `P`. It fills exposed
cells for `U`, `V`, `W`, `PH`, `PHB`, `T`, `PB`, water species, `MU`, `MUB`,
`PSFC`, `U10`, `V10`, `T2`, `Q2`, `RAINC`, and `RAINNC`; exposed `P` remains
pending with `TYWRF_P_DERIVED_REFRESH_STATUS =
pending_derive_or_recompute_after_parent_fill_not_direct_wrf_parent_fill`.

`time_level_sync` is a standalone helper for the WRF post-`start_domain`
exposed-cell `_2 -> _1` copy. It takes POD field views for `U`, `V`, `T`, `W`,
`PH`, `MU`, and optional `TKE`; it is not a `State` expansion and does not
refresh pressure.

The KROSA pressure-refresh helper now exists, and I/O staging exists for the
pressure-refresh constants. Skeleton metadata reports the refresh as required
and not applied, but the compute path is not invoked from the skeleton/remap
sequence yet. Until that connection exists, pressure-refresh metadata is a
blocker marker, not evidence that exposed-cell `P` was produced.

The pressure-refresh consumer still needs `P_TOP`, `C3F`, `C4F`, `C3H`,
`C4H`, and `ALB`. `ALB` is WRF Registry state variable `alb`: inverse base
density on `ikj`, not the surface albedo field `ALBEDO`; this naming constraint
is retained. Its flags include input/restart/interp/smooth handling (`irdus`)
but no history-output flag, so `wrfout` files normally omitting `ALB` is
expected.

The narrow `PB + T_INIT -> ALB` helper is safe only when `PB` and `T_INIT`
already exist for the same domain and valid time. It is not the d02 start-time
source of truth and must not read later restart-file `ALB` as validation truth.
Later restart-file `ALB` may be used only as a probe or smoke reference for
shape and range checks.

Round D27 targets the complete KROSA non-restart `start_domain` base-state
provider:

```text
HGT/P_TOP/C3F/C4F/C3H/C4H/WRF base-atmosphere constants
  -> PB/T_INIT/MUB/ALB/PHB
```

Round D30 defines the skeleton hook helper boundary as:

```text
parent-fill/remap -> provider -> base-state sync -> staging -> pressure refresh compute -> report
```

This hook may synchronize provider-derived `PB`, `MUB`, and `PHB` only into
newly exposed child cells. Existing overlap cells and all halo cells remain
owned by the prior remap/halo-update contracts and must not be overwritten by
the provider sync.

Provider `ALB` is an external pressure-refresh staging input only. It is not a
`State` field and must not be written back into `State`. Provider `T_INIT` is
the WRF base-state initial temperature used by the provider path; it is not the
WRF perturbation potential temperature variable `T` and must never be copied
into `State::t`.

Calling pressure-refresh compute from the helper only means the skeleton has a
producer for exposed-cell pressure inputs. It is not by itself validated
integrator output. Until a real 10 min d02 validation gate passes, hook
metadata and reports must keep `gate_candidate=false` and
`integrator_output=false`.

Later restart `ALB` or `PHB` remains limited to probe/smoke comparisons for
shape, range, and formula checks. It must not be promoted to start-time truth
or used to bypass the provider and recompute path.

The provider reconstructs `PB`, `T_INIT`, and `MUB` from terrain, hybrid
coordinate coefficients, model top, and the WRF base-atmosphere constants, then
derives `ALB` from that reconstructed same-domain base state. For KROSA
`hypsometric_opt=2`, `PHB` is a log-linear full-level reconstruction, not a
mass-level field. With WRF-style 1-based vertical indexing:

```text
PHB(i,1,j) = HGT(i,j) * g
pfu = C3F(k  ) * MUB(i,j) + C4F(k  ) + P_TOP
pfd = C3F(k-1) * MUB(i,j) + C4F(k-1) + P_TOP
phm = C3H(k-1) * MUB(i,j) + C4H(k-1) + P_TOP
PHB(i,k,j) = PHB(i,k-1,j) + ALB(i,k-1,j) * phm * log(pfd / pfu)
```

The `PHB` full-level provider must use `ALB` and `MUB` from the same domain and
same valid time. Mixing a full-level `PHB` reconstruction with mass-level
`ALB`/`MUB` from another time or a later restart is not compatible with the
KROSA start-time contract. Even after this provider is complete, the
pressure-refresh compute path is still not invoked from the skeleton/remap
sequence; staging fields, diagnostic outputs, and provider probes remain
non-gate artifacts until a real pressure refresh is wired and validated.

Later restart-file `ALB` or `PHB` may be used only for probe/smoke checks such
as shape, finite range, and formula sanity. It must not be promoted to d02
start-time validation truth.

Diagnostic parent-fill candidates remain marked `not_physical`,
`diagnostic_only`, and `gate_candidate = false`. They are non-physical,
non-gate artifacts even after direct parent-fill and time-level sync, and they
must not be used as validation passes without a real pressure-refresh compute
path validated against WRF.

Round D28 narrows the real-file base-state provider probe and I/O reader
boundary. The probe exists to verify and provide the KROSA start-domain
base-state pipeline:

```text
HGT/P_TOP/C3F/C4F/C3H/C4H -> PB/T_INIT/MUB/ALB/PHB
```

Its report, probe files, and staging metadata must continue to identify the
artifact as diagnostic, not a gate candidate, and not integrator output. A
successful provider probe confirms that the required real-file constants and
derived base-state fields are readable and internally consistent; it does not
mean the pressure-refresh compute path has been invoked from skeleton/remap,
and it does not count as a 10 min validation-gate pass.

Later restart-file `ALB` or `PHB` remains limited to smoke/range/probe checks.
Those fields must not be used as d02 start-time truth, even when they agree
with the provider within expected ranges.

Round D29 records the adapter/staging bridge boundary. The
`base_state_reconstruction_provider` source may feed pressure-refresh staging
with same-domain, same-valid-time `ALB`, `PB`, `MUB`, and `PHB` buffers. This
is a source and staging contract only: successful provider adapter/staging
bridge execution does not mean perturbation-pressure refresh was computed.
Any report or output with `pressure_refresh_applied = false` remains ineligible
for a validation-gate pass.

Round D30 keeps the same skeleton/remap hook helper ordering:

```text
parent-fill/remap -> provider -> base-state sync -> staging -> pressure refresh compute -> report
```

The hook must keep provider, exposed-cell base-state sync, staging, and
pressure-refresh compute status separate in reports. Later restart-file `ALB`
or `PHB` remains restricted to probe/smoke use and must not become start-time
truth for the pressure-refresh chain.

## Physics Bridge Compatibility Notes

P6 audited the current PGWRF/WRF tree for the v1 physics bridge strategy. The
primary KROSA WRF path is:

```text
/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/model/WRFV4.6.1
```

`main` and `phys` are symlinks into the shared WP11/WP12 bundle; the bundle
contains the full `dyn_em`, `frame`, `Registry`, and `configure.wrf` tree used
to confirm WRF call order and build settings.

The compatible bridge should initially wrap WRF mediation-layer driver
subroutines, not leaf physics kernels. The fixed KROSA suite maps to:

- Thompson: `module_microphysics_driver::microphysics_driver`, selecting
  `module_mp_thompson::mp_gt_driver`;
- RRTMG: `module_radiation_driver::radiation_driver`, selecting
  `rrtmg_lwrad` and `rrtmg_swrad`;
- YSU: `module_pbl_driver::pbl_driver`, selecting `module_bl_ysu::ysu`;
- MM5 surface layer, Noah LSM, slab ocean, and TC flux:
  `module_surface_driver::surface_driver`, selecting `SFCLAY`, `lsm`, and
  `OCEAN_DRIVER` for the KROSA options;
- KF: `module_cumulus_driver::cumulus_driver` on d01 only. d02 uses
  `cu_physics = 0` and should skip cumulus.

The WRF source call order observed for `dyn_em` is radiation, surface, PBL,
cumulus, then microphysics. TyWRF's schedule-driven skeleton can attach its
physics call point to this sequence, but real calls require a wider staging
contract than the current no-op `TywrfPhysicsStaging` v1.

Additional real-physics staging must include WRF physics-grid derived fields
(`u_phy`, `v_phy`, `th_phy`, `p_hyd`, `p_hyd_w`, `pi_phy`, `rho`, `dz8w`,
`z`, `z_at_w`), static/surface fields (`XLAT`, `XLONG`, `XLAND`, `XICE`,
`TSK`, `SST`, `ALBEDO`, `EMISS`, land-use and soil-category fields), physics
tendencies, radiation accumulators, soil/snow state, and slab-ocean state. See
`bindings/wrf_physics/README.md` for the detailed entrypoint manifest and the
recommended minimal ABI spike.
