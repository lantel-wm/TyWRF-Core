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

Round D31 narrows the opt-in skeleton/remap diagnostic path. When
`tywrf_skeleton_cycle` is invoked in the diagnostic parent-fill mode, it may
call the D30 hook and run:

```text
parent-fill/remap -> provider -> base-state sync -> staging -> pressure refresh compute -> report
```

This is still a diagnostic helper path, not validated TyWRF integrator output.
If the report says `pressure_refresh_applied = true`, that means only that the
pressure-refresh helper compute was called and modified exposed-cell `P` in the
generated artifact. It does not imply that the d02 10 min gate passed, that the
full integrator ran, or that the file is physically accepted.

The CLI report and NetCDF metadata for this opt-in path must keep:

```text
TYWRF_GATE_CANDIDATE = false
TYWRF_INTEGRATOR_OUTPUT = false
TYWRF_VALIDATION_GATE_ONLY = false
```

Direct `ALB` is not required from a history-output file. If direct `ALB` is
missing but same-domain base-state reconstruction inputs are available, the
path may use the provider to stage `ALB` for pressure refresh. Later restart
`ALB` or `PHB` remains limited to probe/smoke comparisons and must not become
start-time truth.

When pressure refresh is requested, any `--variables` selection and the output
field list must include the fields required by the pressure-refresh path,
including the exposed `P` target and the provider/staging fields used to
refresh it (`PB`, `MUB`, `PHB`, and staged provider `ALB`, plus the KROSA
constants `P_TOP`, `C3F`, `C4F`, `C3H`, and `C4H`). A narrowed variable list
that omits required inputs must fail early rather than refreshing an
uninitialized or partially staged state.

## Round D33-D35 Compatibility Boundary

Round D33 separates diagnostic accounting from WRF-compatible integrator output.
Strict gate eligibility now requires positive candidate metadata, not merely the
absence of a diagnostic marker. A gate candidate must identify itself as real
TyWRF integrator output with `TYWRF_GATE_CANDIDATE = true` and
`TYWRF_INTEGRATOR_OUTPUT = true`; artifacts marked as diagnostic, validation
plumbing, oracle/reference-copy, or non-physical remain incompatible with a
strict gate pass.

Candidate kinds that describe oracle behavior are rejected by compatibility
policy. This includes reference copies, WRF end-state copies, WRF end-state
delta products, diagnostic remap outputs, diagnostic pressure-refresh outputs,
diagnostic closures, and any other artifact whose fields are copied from or
corrected by a later WRF truth file. Such files may be useful to test I/O or
report formatting, but they must not be documented as WRF-compatible physical
cycle output.

`report_10min_diagnostics` gate-eligibility data is descriptive context only.
It can explain which metadata would block or allow a strict gate attempt, but
it does not itself produce a pass/fail acceptance. The strict d02 gate remains
the normal validation comparison against the `00:10` WRF reference.

`state_exchange` currently records intended nesting exchange accounting only.
While `performed_interpolation = false`, it means no parent-to-child numerical
interpolation has been applied. Field lists, cell counts, exchange regions, and
planned source/target ranges in that block are planning metadata, not evidence
of a physical d02 parent-fill or a WRF-compatible nesting result.

Round D34 should begin selected-field parent-to-child interpolation as a narrow
compatibility slice. The allowed source set is the cycle start state plus
d01-derived parent fields and KROSA constants available at the same valid time.
The path must not read later d02 restart truth, must not use WRF end-state
delta fields, and must not copy or blend WRF cycle-end output into the child
candidate.

The first selected-field interpolation reports should name exactly which
fields were interpolated, which exposed cells were filled, and whether pressure
refresh still remains pending. They must preserve d02 `DX/DY = 2000 m`, avoid
best-track nudging, and avoid any gate-pass claim until the strict
`2025-07-26_00:10:00` d02 gate passes with a positive-metadata integrator
candidate.

Round D35 defines `selected_field_cycle` as the first intended non-oracle
selected-field candidate path. Its compatibility boundary is the cycle-start
d01/d02 state only, followed by moving-nest remap,
`StateExchangePlan`, and `parent_child_interpolation` for the selected fields.
It must not read WRF reference-end files, later restart truth, end-state
deltas, or any reference-end-derived correction. The output must include every
strict d02 gate field required by the validation plan, and those fields must be
finite so failures are numerical RMSE/diagnostic failures rather than
metadata, missing-field, or sentinel-field failures.

`selected_field_cycle` is not WRF-exact physics. Its first real outputs are
expected to fail the strict `00:10` gate numerically while still exercising the
non-oracle integrator path. It may set the positive candidate metadata only
when all of the following are true:

```text
TYWRF_GATE_CANDIDATE = true
TYWRF_INTEGRATOR_OUTPUT = true
TYWRF_VALIDATION_GATE_ONLY = false
TYWRF_CANDIDATE_KIND = selected_field_integrator_v0
```

Those attributes are valid only if no end-state/reference-end inputs were
used, d02 remains `2 km`, all strict fields are finite, and at least one
selected field changed from the cycle-start child state through the remap /
exchange / interpolation path. If any condition is not met, the file must not
claim positive gate-candidate metadata and the report should name the blocking
condition. `tendency_cycle`, persistence/skeleton, diagnostic remap,
diagnostic pressure-refresh, and diagnostic closure paths remain explicitly
non-gate compatibility artifacts.

Round D35 real KROSA selected-field smoke produced the first positive-metadata
`selected_field_cycle` candidate from `2025-07-26_00:00:00` start states only.
The candidate used the non-oracle selected-field path, preserved d02
`DX/DY = 2000 m`, and passed the strict gate metadata checks:

```text
TYWRF_GATE_CANDIDATE = true
TYWRF_INTEGRATOR_OUTPUT = true
TYWRF_VALIDATION_GATE_ONLY = false
TYWRF_CANDIDATE_KIND = selected_field_integrator_v0
```

This smoke is not a gate pass. The strict RMSE fields had full finite coverage,
so the first failure is numeric: `U` normalized RMSE `0.117875`. `V`, `MU`, and
`P` also fail, with `P` normalized RMSE `0.907405` as the largest current
blocker. `T`, `PH`, and `QVAPOR` pass their strict field thresholds. Remaining
output/diagnostic gaps are that accepted sea-level pressure is unavailable,
`U10` is absent in the v0 smoke output, and Vmax/MSLP diagnostics therefore
cannot yet close the full strict d02 gate.

Round D36 kept the same compatibility boundary and improved output
completeness; it did not produce a strict-gate pass. Commit `7d8d37c`
(`Preserve selected-field diagnostic outputs`) makes `selected_field_cycle`
write the strict d02 fields plus available cycle-start-backed `PB`, `PHB`,
`MUB`, `PSFC`, `U10`, `V10`, `T2`, `Q2`, `RAINC`, and `RAINNC` by default.
This is plumbing/output completeness work so validation tools can see the
candidate shape and near-surface diagnostic variables instead of failing early
on missing files or missing variables.

The D36 real KROSA gate still fails. The strict first failure remains `U`
normalized RMSE `0.117875`; `V` normalized RMSE `0.134244`, `MU` normalized
RMSE `0.133382`, and `P` normalized RMSE `0.907405` also fail. `T`, `PH`, and
`QVAPOR` pass. TC diagnostics report storm-center error `150.622 km` failed,
minimum SLP error `0.364 hPa` passed after derived SLP staging, and Vmax error
`0.769 m s-1` passed. The `U10`/`V10` field comparisons are still failed
because these variables are preserved cycle-start fields, not outputs from a
real surface-diagnostics producer.

Positive selected-field metadata remains exclusive to the non-oracle
`selected_field_cycle` candidate path. It must not be copied onto or used to
legitimize diagnostic closures, diagnostic pressure-refresh artifacts,
diagnostic remap outputs, reference-copy/oracle outputs, derived SLP/rainfall
side products, or report-only files. Metadata or helper telemetry attached to
those artifacts is compatibility context only and must not override the
candidate file's own gate metadata.

Round D37 pressure-refresh work, if enabled, is opt-in and must consume only
cycle-start or provider-backed inputs. Missing, shape-incompatible, non-finite,
or reference-end-derived pressure-refresh inputs must abort the candidate
generation or diagnostic request; they must not silently fall back to
start/end deltas, later restart truth, `PSFC` pressure shortcuts, or helper
telemetry that marks the candidate as successful. Helper reports may describe
provider, staging, and compute status, but they cannot cover or rewrite
candidate NetCDF metadata.

The D37 real KROSA opt-in smoke confirms the plumbing but not the numerics.
`selected_field_cycle --pressure-refresh` generated a positive-metadata
candidate from d01/d02 `00:00` start state and d02 `00:00` template metadata,
with provider, staging, and compute all reporting success and `1,053,150`
changed `P` points. The strict gate still failed, and `P` worsened from the
D36 normalized RMSE `0.907405` to `5.806413`. `U`, `V`, and `MU` remained
failed, while `T`, `PH`, and `QVAPOR` remained passed. This opt-in path must
therefore stay off by default until provider/base-state/static consistency is
fixed.

The remaining compatibility blockers are real numerical producers, not
metadata gaps: moving-nest static/coordinate handling for the shifted d02
domain, pressure refresh as a real candidate producer, and surface diagnostics
for `U10`, `V10`, `T2`, `Q2`, and related fields. D37 and later work must not
use `00:10` reference-end static or coordinate truth, best-track nudging, an
SLP shortcut, or a `U10`/`V10` proxy shortcut to make those blockers appear
closed.

## Round D38 Compatibility Boundary

D37 commit `aba2469` added selected-field `--pressure-refresh` opt-in plumbing.
It is safe wiring, not a numerical improvement and not a default behavior. The
real KROSA opt-in run reported provider, staging, and compute success and
changed `1,053,150` `P` points, but the strict gate still failed: `P`
normalized RMSE worsened from D36 `0.907405` to `5.806413`; `U` `0.117875`,
`V` `0.134244`, and `MU` `0.133382` still failed; `T`, `PH`, and `QVAPOR`
still passed; storm-center error remained `150.622 km`. Compatibility notes,
CLI help, and reports must not describe `--pressure-refresh` as a gate pass or
enable it by default.

D38 work is split across three producer boundaries:

- Static refresh: G37 found stale shifted-d02 static coordinates
  (`XLAT` RMSE `0.6913 deg`, `XLONG` RMSE `1.2463 deg`, same-index
  displacement `152.6 km`). The fix must derive shifted d02 `XLAT`, `XLONG`,
  and `HGT` from the cycle-start moving-nest pose shift and parent/static
  interpolation or reconstruction path. It must not copy `00:10` reference-end
  `XLAT`, `XLONG`, `HGT`, or other shifted-domain truth into the candidate.
- Pressure negative diagnosis: the D37 refresh changed `P` but made validation
  worse, so the next pressure work is diagnosis of provider/static/vertical
  coordinate consistency, exposed-vs-overlap cell handling, and
  `PB`/`MUB`/`PHB`/`ALB` alignment. It remains opt-in until the same
  non-oracle candidate path improves or passes the real `00:10` strict gate.
- Surface ABI v2: preserved `U10`, `V10`, `T2`, `Q2`, `RAINC`, and `RAINNC`
  keep output files complete, but preserved values are not real producers. The
  compatible path must be a real ABI v2 surface/physics bridge, starting with
  SFCLAY/physics-wrapper-owned diagnostics, not metadata relabeling or
  start-state preservation.

The first D38 static-refresh smoke confirms the static side of that boundary.
The selected-field candidate no longer leaves d02 shifted-domain static fields
as the `00:00` template copy: compared with the `00:10` KROSA d02 reference,
`XLAT` RMSE is `0.000287 deg`, `XLONG` RMSE is `0.000195 deg`, and `HGT` RMSE
is `0.541 m` with maximum absolute `HGT` error `20.789 m`. The strict gate
still fails, but storm-center error drops from `150.622 km` to `43.483 km`.
This validates the stale-coordinate fix while keeping the remaining state and
pressure errors visible.

## Round D39 Compatibility Boundary

D38 commit `517ad28` (`Refresh selected-field moving nest statics`) is an
effective shifted-d02 coordinate/static fix, not a validation pass. The current
real KROSA selected-field candidate preserves d02 `DX/DY = 2000 m` and improves
static fields to `XLAT` RMSE `0.000287 deg`, `XLONG` RMSE `0.000195 deg`, and
`HGT` RMSE `0.541 m` with maximum absolute `HGT` error `20.789 m`. Storm-center
distance improves from `150.622 km` to `43.483 km`, but the `00:10` strict gate
still fails on storm-center distance and on `U` `0.117875`, `V` `0.134244`,
`MU` `0.133382`, and `P` `0.907405` normalized RMSE. `T`, `PH`, and `QVAPOR`
passing their thresholds does not close the moving-nest producer gap for the
failed fields.

`--pressure-refresh` is currently not a compatibility improvement and must stay
guarded/off by default. D37 showed that the opt-in path worsened `P` normalized
RMSE to `5.806413`; P38 pressure diagnostics also found same-source wrfout
self-consistency error around 10% for the current refresh formula/staging. A
candidate-facing pressure refresh therefore needs an explicit readiness guard
and self-consistency probe before it can run as a producer. Failing that guard
must abort or leave the candidate on the non-refresh path; it must not silently
enable refresh, relabel helper telemetry as success, or claim improved
compatibility.

The exposed-cell `T`/`PH`/`P` same-pose issue is a producer/interpolation
problem, not a reference-truth problem. Analysis may compare exposed cells at a
common reconstructed child pose to understand which field family is inconsistent,
but the production path must use cycle-start fields, d01-derived parent data,
same-time KROSA constants, and documented interpolation/recompute rules. It
must not read or blend `00:10` WRF d02 truth for `T`, `PH`, `P`, statics, or
storm position.

The surface ABI v2 work is a scaffold boundary only until a real SFCLAY or
physics-wrapper surface producer is invoked. ABI definitions, staging buffers,
or report fields may describe planned `U10`, `V10`, `T2`, `Q2`, `RAINC`, and
`RAINNC` ownership, but they must not claim executed WRF physics, executed
surface diagnostics, or a closed surface producer while the values remain
cycle-start-preserved or scaffold-only.

## Round D40 Compatibility Boundary

D39 commit `0ea9a69` (`Guard pressure refresh and scaffold physics ABI v2`)
keeps pressure refresh behind a readiness guard. In the current real KROSA
smoke, `--pressure-refresh` aborts under that guard and does not generate an
output file. This is the intended compatibility behavior: the guard is a safety
constraint that prevents an unready pressure path from contaminating a
gate-eligible candidate. It is not a pressure producer and must not be reported
as evidence that exposed-cell `P` has been refreshed.

The ABI v2 sidecar scaffold and fixture are compatibility scaffolding only. The
current `_ex` wrapper path reports `wrapper_unavailable` and
`executed_physics = 0`; therefore sidecar fixture success does not mean SFCLAY,
surface diagnostics, precipitation, or any WRF physics package executed. ABI v2
metadata may document argument shape, ownership, and fixture plumbing, but it
must not be counted as a producer for `U10`, `V10`, `T2`, `Q2`, `RAINC`, or
`RAINNC` until the wrapper actually runs and writes those fields.

D40 should keep the next parent-fill slice narrow: exposed-cell `T` and `PH`
may be produced by non-oracle parent-child interpolation from cycle-start d02
state, d01-derived parent fields, same-time KROSA constants, and documented
moving-nest geometry. It must not use `00:10` WRF d02 truth as source data,
delta, bias correction, or storm-position guidance. `P` must not be direct
parent-interpolated into a production field; it remains owned by a derived or
recomputed pressure producer once readiness checks pass.

`PB`, `PHB`, and `MUB` have a read-only ownership boundary for this slice. They
may be consumed as same-domain, same-valid-time base-state/provider inputs for
`T`/`PH` interpolation context and later pressure-refresh staging, but D40 must
not overwrite them through selected-field interpolation or use later restart
truth to repair them. Any future write path for these base-state fields needs a
separate provider/sync design that preserves overlap, halo, and moving-nest
ownership rules.

## Round D41/D42 Compatibility Boundary

D41 commit `fe01be1` (`Allow moved terrain in base-state provider`) adds
`KrosaBaseStateProviderTerrainOverride`. This lets the base-state provider
reconstruct its buffers from an explicit moved candidate `HGT` view instead of
only metadata terrain. That is a provider/probe capability, not a pressure
producer: using moved terrain in the provider does not by itself write `P`,
`PB`, `PHB`, `MUB`, or any other candidate field.

D42 pressure-refresh work must keep that boundary read-only. The intended
selected-field `--pressure-refresh` readiness probe may inject moved candidate
`HGT` into the provider and verify that base-state buffers can be reconstructed,
but the probe result must not be synced back into the candidate NetCDF. Even
when the provider terrain source is `moved_candidate_HGT`,
`thermodynamic_base_state_consistency_ready` remains `false`, so
`--pressure-refresh` must still fail closed and produce no output.

The default selected-field candidate continues to preserve `P`, `PB`, `PHB`,
and `MUB` from start-state ownership. `P` must not be added to parent
interpolation; only the documented non-oracle parent-fill variables may move
through that path. D41/D42 also must not use `00:10` reference-end truth,
restart `PHB`/`ALB` in place of the provider, end-state deltas, or probe and
diagnostic pressure-refresh artifacts as compatibility evidence or gate passes.
The current `00:10` gate is still failed on `U`, `V`, `MU`, `P`, and
storm-center distance.

## Round D43 Hook Terrain-Override Boundary

D43 may extend the pressure-refresh hook/API so the hook can ask the
base-state provider to use a moved-terrain override. This is diagnostic hook
capability only. It does not make the default selected-field candidate a
pressure producer, and it does not change the D42 opt-in readiness guard.

A hook-level unit test may exercise the terrain override and refresh exposed
`P` through the diagnostic pressure-refresh hook. That proves the hook can pass
override terrain into provider-backed staging for an exposed-cell smoke path; it
does not prove that `selected_field_cycle --pressure-refresh` is ready, and it
must not be reported as a real `00:10` gate pass.

The default selected-field candidate and the D42 opt-in guard still keep
`thermodynamic_base_state_consistency_ready=false`. The opt-in path must fail
closed, produce no output file, and avoid writing `P`, `PB`, `PHB`, or `MUB` to
the candidate. `P` remains forbidden from parent interpolation. D43 also must
not use `00:10` reference-end truth, restart `PHB`/`ALB` as a provider
substitute, or probe/diagnostic artifacts as gate evidence. The current real
`00:10` gate remains failed on `U`, `V`, `MU`, `P`, and storm-center distance.

## Round D44 Selected-Field Pressure-Refresh Readiness Contract

D41 through D43 complete only the capability chain needed to inspect the next
pressure-refresh step: moved `HGT` provider override, selected-field moved-HGT
provider probe, and hook terrain override. This chain is not yet a selected-field
pressure producer. The default selected-field candidate must continue to
preserve `P`, `PB`, `PHB`, and `MUB` from start-state ownership.

The next stage must define a controlled readiness contract before
`selected_field_cycle --pressure-refresh` can be allowed to write a candidate.
At minimum, the opt-in path must report provider terrain source/provenance,
base-state reconstruction status, planned base-state and pressure sync counts,
overlap and halo safety for the moved nest, and pressure compute/requested/
applied status. Missing, inconsistent, non-finite, or unsafe readiness data must
keep the path fail-closed.

While the opt-in path is not ready, `--pressure-refresh` must abort before
writing output and must not modify `P`, `PB`, `PHB`, or `MUB`. A hook-level
diagnostic smoke may refresh exposed `P`, but that is hook evidence only; it is
not selected-field readiness, not a gate-eligible candidate, and not a real
`00:10` pass. Shortcuts remain forbidden: `00:10` reference-end truth, restart
`PHB`/`ALB` as a provider substitute, direct parent interpolation of `P`, and
diagnostic artifacts as gate evidence. The current real `00:10` gate remains
failed on `U`, `V`, `MU`, `P`, and storm-center distance.

## Round D45/D46 Dry-Run Contract Boundary

D45 adds a pressure-refresh hook dry-run report for the base-state sync
contract. The dry-run path may reconstruct provider-backed buffers and report
would-sync counts for `PB`, `MUB`, and `PHB`, plus overlap/halo write safety and
no-write/no-compute status. This is contract evidence only: dry-run must keep
`base_state_sync_applied=false`, must not call pressure compute, and must not
modify `P`, `PB`, `PHB`, `MUB`, or any candidate output field.

D46 selected-field work may connect that hook dry-run report into the
`selected_field_cycle --pressure-refresh` opt-in guard so the abort message and
metadata expose the contract status. The path must still fail closed while
`thermodynamic_base_state_consistency_ready=false`. Even if the dry-run contract
is internally consistent and reports zero overlap/halo writes, the tool must not
generate an output file, must not write `P`/`PB`/`PHB`/`MUB`, and must not be
classified as a compatibility gate pass.

Dry-run fields and hook diagnostic smoke remain readiness evidence only. They
do not replace the strict `2025-07-26_00:10:00` selected-field gate and cannot
be used as validation evidence for the real integrator. The current real
`00:10` gate remains failed on `U`, `V`, `MU`, `P`, and storm-center distance.
Shortcuts remain forbidden: `00:10` reference-end truth, restart `PHB`/`ALB` as
a provider substitute, direct parent interpolation of `P`, and diagnostic
artifacts as gate evidence.

## Round D47 Pressure Compute Scratch Dry-Run Boundary

D47 extends only the pressure-refresh hook evidence boundary. The target is a
scratch-buffer pressure compute dry-run inside the hook layer: it may run the
pressure compute path against provider-backed scratch `P` storage and report
would-refresh `P` point counts, invalid `P` point counts, and an explicit safety
status. This is compute-readiness evidence, not candidate production.

The scratch dry-run must not write candidate `P`, `PB`, `PHB`, or `MUB`, must
not generate selected-field output, and must not report
`pressure_refresh_applied=true`. Any report should keep the distinction between
base-state sync dry-run, scratch pressure compute dry-run, and real pressure
refresh application.

The selected-field tool layer is not connected in this round. The D46
`selected_field_cycle --pressure-refresh` opt-in guard therefore remains
fail-closed while `thermodynamic_base_state_consistency_ready=false`, and the
default selected-field candidate still preserves start-state-owned pressure and
base-state fields. Scratch compute evidence is not a compatibility gate pass:
the current real `2025-07-26_00:10:00` gate remains failed on `U`, `V`, `MU`,
`P`, and storm-center distance.

Shortcuts remain forbidden: `00:10` reference-end truth, restart `PHB`/`ALB` as
a provider substitute, direct parent interpolation of `P`, diagnostic artifact
gate evidence, and treating scratch dry-run telemetry as selected-field output
or integrator validation evidence.

## Round D48 Selected-Field Scratch Telemetry Boundary

D48 may expose the D47 hook scratch pressure compute dry-run telemetry through
the selected-field `--pressure-refresh` opt-in guard. The guard report may name
whether the dry-run was requested, called, and ok, plus the would-refresh `P`
point count and invalid `P` point count. These fields are readiness evidence
only; they are not selected-field output values and do not show that pressure
refresh has been applied to the candidate.

While `thermodynamic_base_state_consistency_ready=false`, the opt-in path must
still exit nonzero and produce no output file. The scratch dry-run must not
write candidate `P`, `PB`, `PHB`, or `MUB`, and must not set
`pressure_refresh_applied=true`. The default selected-field candidate numerical
path remains unchanged and continues to preserve start-state-owned pressure and
base-state fields.

Scratch telemetry is not a compatibility gate pass. The real
`2025-07-26_00:10:00` selected-field gate remains failed on `U`, `V`, `MU`,
`P`, and storm-center distance. Shortcuts remain forbidden: `00:10`
reference-end truth, restart `PHB`/`ALB` as a provider substitute, direct
parent interpolation of `P`, diagnostic artifact gate evidence, and treating
scratch dry-run telemetry as gate evidence.

## Round D49 Future Apply Terrain-Parity Boundary

D49 preparation is limited to semantic parity between the future selected-field
pressure-refresh apply path and the existing dry-run path. If a later readiness
change makes the apply path reachable, that path must use the same moved
candidate `HGT` terrain override that the dry-run provider probe uses, rather
than silently falling back to metadata terrain.

This is an apply-before-enable alignment only. It does not mean
`selected_field_cycle --pressure-refresh` has generated output, refreshed
candidate `P`, or become gate-eligible. While
`thermodynamic_base_state_consistency_ready=false`, the opt-in path must still
fail closed with a nonzero exit and no output file. It must not write `P`,
`PB`, `PHB`, or `MUB`, and it must not report
`pressure_refresh_applied=true`.

The default selected-field candidate numerical path remains unchanged and
continues to preserve start-state-owned pressure and base-state fields. Even
when future apply-path terrain parity is ready, scratch telemetry and an
unreachable guarded apply path are not compatibility gate evidence. The real
`2025-07-26_00:10:00` selected-field gate remains failed on `U`, `V`, `MU`,
`P`, and storm-center distance. Shortcuts remain forbidden: `00:10`
reference-end truth, restart `PHB`/`ALB` as a provider substitute, direct
parent interpolation of `P`, diagnostic artifact gate evidence, and treating
scratch telemetry or unreachable apply plumbing as a gate pass.

## Round D50 Pressure-Refresh Apply-Readiness Blockers

D50 documents the remaining blocker list before selected-field pressure-refresh
apply may become reachable. Current readiness evidence is limited to four
pre-apply capabilities: the moved-`HGT` provider probe, the base-state sync
dry-run contract, the scratch pressure compute dry-run, and D49 future apply
moved-`HGT` terrain parity. Together these show that provider, dry-run, scratch
compute, and future apply terrain sources can be inspected consistently; they
do not show candidate mutation safety or validation success.

D50 must not open `thermodynamic_base_state_consistency_ready`. While that flag
remains `false`, `selected_field_cycle --pressure-refresh` must still fail
closed with a nonzero exit, no output file, no candidate writes to `P`, `PB`,
`PHB`, or `MUB`, and no `pressure_refresh_applied=true`. The default
selected-field candidate numerical path remains unchanged and continues to
preserve start-state-owned pressure and base-state fields.

The minimum safety conditions before apply can be enabled are:

- candidate mutation audit proving exactly which candidate fields may change;
- moved-`HGT` apply proof showing the reachable apply path uses the same moved
  candidate terrain override as the provider probe and dry-run path;
- overlap/halo no-write proof showing pressure/base-state writes stay out of
  remapped overlap cells that must be preserved and all halo cells;
- pre/post apply report consistency showing requested, computed, applied, and
  changed counts agree before any output is written;
- a real `2025-07-26_00:10:00` validation plan that compares the produced d02
  file against the 10 min WRF reference and stops at the first failing field.

Scratch telemetry and an unreachable guarded apply path are not compatibility
gate evidence, even when their reports are internally consistent. Shortcuts
remain forbidden: `00:10` reference-end truth, restart `PHB`/`ALB` as a provider
substitute, direct parent interpolation of `P`, diagnostic artifact gate
evidence, and treating scratch telemetry or unreachable apply plumbing as a
gate pass.

## Round D51 Hook Postcondition And Controlled Apply Seam Boundary

D51 still must not open the selected-field pressure-refresh guard by default.
For normal use, `thermodynamic_base_state_consistency_ready=false` remains the
contract: `selected_field_cycle --pressure-refresh` must fail closed with a
nonzero exit, produce no output file, make no candidate writes, and report no
`pressure_refresh_applied=true`. The default selected-field candidate numerical
path remains unchanged.

The hook postcondition hardening target is narrow. If a controlled path reaches
real pressure-refresh apply, a successful apply may change only `P`, `PB`,
`MUB`, and `PHB`. It must prove zero writes to remapped overlap cells and all
halo cells. Invalid compute, non-finite compute, inconsistent counts, or any
partial compute/write sequence cannot be reported as a successful apply.

The controlled selected-field seam is a test boundary only. Its purpose is to
prove that the tool-level apply call consumes the moved candidate `HGT` terrain
used by the selected-field candidate, not metadata terrain from the template or
source file. Seam and experiment outputs are not gate-eligible, are not
compatibility passes, and must not be promoted into normal CLI behavior.

D51 keeps the same shortcut bans as D50: no `00:10` reference-end truth, no
restart `PHB`/`ALB` provider substitute, no direct parent interpolation of `P`,
and no diagnostic artifact gate evidence. A postcondition test or controlled
seam proof is readiness evidence only until a real selected-field
`--pressure-refresh` candidate is produced by the normal guarded path and
passes the progressive `00:10` gate.

## Round D52 Transactional Apply And Report Parity Boundary

D52 still must not open the selected-field `--pressure-refresh` guard by
default. The normal `thermodynamic_base_state_consistency_ready=false` contract
remains fail-closed: nonzero exit, no output file, no candidate writes, no
`pressure_refresh_applied=true`, and no change to the default selected-field
candidate numerical path.

The hook apply target is transactional. A pressure refresh must first compute
into scratch buffers and evaluate the full postcondition there. Only after the
postcondition succeeds may the hook write allowed fields into the candidate. If
scratch compute fails, produces invalid values, touches unsafe regions, or has
inconsistent counts, the candidate must remain completely unchanged, including
`P`, `PB`, `MUB`, `PHB`, remapped overlap cells, halo cells, and unrelated
fields.

The selected-field experimental apply seam remains report-parity evidence only.
If explicitly enabled, its JSON/report and metadata must record the pressure
refresh target, refreshed, skipped, and invalid counts; overlap and halo safety
flags; changed counts for candidate fields; and consistency flags comparing
target/refreshed/changed counts. These records do not make the seam output a
normal candidate.

Hidden seam output remains diagnostic-only, non-gate, and not normal integrator
output. It must carry negative gate/integrator metadata such as
`TYWRF_DIAGNOSTIC_ONLY=true`, `TYWRF_GATE_CANDIDATE=false`, and
`TYWRF_INTEGRATOR_OUTPUT=false`, and must not be promoted into default CLI
behavior or compatibility-pass evidence.

D52 keeps the same shortcut bans as D51: no `00:10` reference-end truth, no
restart `PHB`/`ALB` provider substitute, no direct parent interpolation of `P`,
and no diagnostic artifact gate evidence. Transactional hook tests, report
parity, and hidden seam outputs are readiness evidence only until the normal
guarded selected-field path produces a real candidate and passes the progressive
`00:10` gate.

## Round D53 Strict Readiness-Derived Normal Apply Boundary

D53 may replace the hard-coded
`thermodynamic_base_state_consistency_ready=false` guard with a strictly derived
readiness decision, but only from normal-path provider, dry-run, and
transactional apply evidence. Missing, inconsistent, diagnostic-only, or
hidden-seam evidence must keep `selected_field_cycle --pressure-refresh`
fail-closed with no output file, no candidate writes, and no
`pressure_refresh_applied=true`.

A normal selected-field `--pressure-refresh` candidate may be written only when
all of the following are true:

- moving-nest static refresh is start-state based and does not copy or blend
  `00:10` reference-end static truth;
- the pressure provider uses the selected candidate terrain source
  `moved_candidate_HGT`, not template metadata terrain, restart terrain, or a
  diagnostic substitute;
- the base-state sync dry-run and pressure compute dry-run both pass for the
  same domain, valid time, terrain source, and field selection;
- overlap, halo, invalid-point, skipped-point, and count-consistency reports
  all show no unsafe or unresolved pressure-refresh work;
- the real apply is transactional: scratch compute and postcondition checks
  succeed before candidate writes, and the only fields eligible to change are
  `P`, `PB`, `MUB`, and `PHB`.

These conditions make a produced normal-path artifact candidate-eligible, not a
compatibility pass. The hidden `--experimental-pressure-refresh-apply` seam
remains diagnostic-only, non-gate, and non-integrator output. Probe, dry-run,
diagnostic, hidden-seam, and report-parity artifacts remain invalid as gate
evidence. The first real validation after D53 must still be the progressive
`2025-07-26_00:10:00` gate; if it fails, reports must name the first failed
field or diagnostic and must not advance to `00:20`.

## Round D54 Post-D53 Compatibility State

D53 promoted the normal `selected_field_cycle --pressure-refresh` path from a
fail-closed guard to gate-eligible output only under strict readiness:
start-state static refresh, `moved_candidate_HGT` provider terrain, passing
base-state and pressure dry-runs, clean overlap/halo/invalid/skipped checks,
and successful transactional apply for the same candidate. This is a
compatibility eligibility state, not WRF numerical compatibility.

The D53 real KROSA d02 `2025-07-26_00:10:00` gate failed after candidate
metadata passed. The first failed field was `U` normalized RMSE `0.117875`;
`V` `0.134244`, `MU` `0.133382`, and pressure-refresh `P` `6.33495` also
failed. `T` `0.013121`, `PH` `0.017410`, and `QVAPOR` `0.017017` passed.
Storm-center distance remained failed at `43.483 km`, while minimum SLP error
`0.364 hPa` and Vmax error `0.769 m s-1` passed.

D54 compatibility work must keep the path stopped at `00:10`. The next
candidate-facing work is diagnosis of the `U` first failure and the numerical
`P` pressure-refresh error; no `00:20` progression is compatible until the
normal strict `00:10` gate passes. Passing apply contracts, moved-terrain
parity, or report-count parity only proves readiness plumbing. It does not
show WRF-compatible pressure until `P` is compared numerically against the
`00:10` WRF reference and meets the gate threshold.

Audit reports are compatibility diagnostics only. They may describe metadata,
reference coverage, provenance, changed fields, or why the gate failed, but
they must not produce candidates, repair a candidate into gate eligibility, or
substitute for strict-gate evidence. The hidden
`--experimental-pressure-refresh-apply` seam remains diagnostic-only,
non-gate, and non-integrator output.

## Round D55 Diagnostic Direction

D54 audit of the pressure-refresh candidate shows that the `P` failure is
concentrated in the refreshed target region, not in the untouched field
background. Global `P` normalized RMSE is `6.334952`; target-region `P`
normalized RMSE is `10.33536`; non-target `P` normalized RMSE is `0.338162`;
and the target region accounts for `0.998219` of total `P` squared error.
`PB`, `MUB`, and `PHB` are close enough to pass their audit checks, so D55
should treat the pressure problem as a perturbation-pressure refresh semantics
issue before changing base-state ownership.

The D54 derived-SLP audit can pass with normalized RMSE `0.003248` while the
strict gate still fails `U`, `V`, `MU`, `P`, and storm-center distance. This is
not contradictory: derived SLP is an additional diagnostic field, not a
replacement for the strict `P` field comparison or the TC center requirement.
A passing derived-SLP audit must not be documented or reported as overriding
the failed `P` gate, the failed storm-center diagnostic, or the stop-at-`00:10`
rule.

D55 compatibility work should therefore remain diagnostic-only and focus on two
audits:

- pressure-refresh `P` formula/source/vertical/region semantics, including the
  exact fields and constants used to refresh the target region;
- selected-field `U`, `V`, and `MU` state-error splitting, separating remapped
  overlap, newly exposed target cells, and untouched/non-target cells.

These audits may read WRF reference output only to compute diagnostic
differences, regional masks, and formula/source attribution. Reference output
must not be used to generate candidate values, choose a candidate correction,
tune a refresh formula into the output, move the storm, or promote a diagnostic
artifact into a gate candidate. The hidden pressure-refresh seam remains
diagnostic-only, non-gate, and non-integrator output, and audit tools remain
observation-only rather than candidate producers.

## Round D56 Diagnostic Continuation

D55 integration commit `cc39cdb` preserves the same compatibility state after
full validation: CMake build passed, CTest passed `29/29`, and pytest passed
`141/141`. The integrated outcome is still a failed strict `00:10` gate, not a
compatibility pass.

The pressure formula/source audit remains diagnostic-only. The refreshed
target-region perturbation pressure dominates the `P` failure: global `P`
normalized RMSE is `6.334952`, target-region `P` normalized RMSE is
`10.335362`, non-target `P` normalized RMSE is `0.338162`, target error
fraction is `0.998219`, and target bias is about `-805 Pa`. Full-pressure
normalized RMSE may look small because full pressure has a much larger
reference scale, but the raw RMSE is nearly unchanged, so full pressure cannot
replace strict perturbation-pressure `P` gate evidence. The candidate `P` is
not close to start `P`, so the failure is not a simple persistence baseline.

The selected-field state audit also remains diagnostic-only. The strict gate's
first failed variable is still `U`; the selected-field failures remain `U`
normalized RMSE `0.117875`, `V` `0.134244`, and `MU` `0.133382`. Target-region
error fractions below `0.5`, together with candidate fields that are not close
to the start state, point to a broader selected-field and dynamics gap rather
than an isolated exposed-cell interpolation problem.

D56 compatibility work should continue with observation-only diagnostics:

- a vertical/source-level target-region pressure audit for the refreshed `P`
  source terms, levels, constants, and target mask;
- a selected-field evolution audit comparing candidate evolution against WRF
  reference evolution from the shared start state.

These audits may attribute errors but must not generate candidates, patch
candidate fields, tune formula output from reference-end truth, or promote
diagnostic/probe/hidden-seam artifacts into gate evidence. The validation
sequence remains stopped at `2025-07-26_00:10:00`; no `00:20` progression is
compatible until the strict `00:10` gate passes.

## Round D57 Diagnostic Continuation

D56 commit `985a38f` (`Audit pressure vertical and state evolution`) preserves
the same compatibility state after full validation: CMake build passed, CTest
passed `29/29`, pytest passed `150/150`, and push to `origin/main` succeeded.
This is codebase validation, not a WRF-compatible numerical pass. The strict
d02 `2025-07-26_00:10:00` gate remains failed, and no `00:20` progression is
compatible until that endpoint passes.

The pressure source audit top-level status may remain failed because real
start `wrfout` source entries still miss `ALB`. That status describes
source-seed audit semantics for the pressure diagnostic inputs; it is not a
model pass/fail override and must not be used to promote or demote a candidate
outside the normal strict gate.

The D56 pressure diagnostic subreport keeps the pressure-refresh failure
localized to low-level target-region perturbation pressure. Target-region `P`
RMSE is `1560.573 Pa`, normalized RMSE is `10.335362`, and mean bias is
`-805.059 Pa`. The worst levels by RMSE and bias are low-level mass levels
`0..4`; level `0` has RMSE `4170.444 Pa` and bias `-4167.791 Pa`. The
source/start evolution fraction is `6.57739`, so the compatibility issue
currently looks like a low-level target-region pressure formula/source problem,
not a metadata eligibility problem.

The selected-field evolution audit also remains diagnostic-only. Evolution
normalized RMSE is `U = 0.156800`, `V = 0.160006`, and `MU = 0.197589`.
Target-region fractions are `U = 0.414932`, `V = 0.362148`, and
`MU = 0.329146`; amplitude ratios are close to `1`, and capture fractions are
high but imperfect. This points to a spatial, phase, or timing mismatch in the
selected-field evolution rather than an amplitude-only failure or an error
confined only to the target region.

D57 compatibility work should therefore add observation-only diagnostics:

- a pressure vertical-bias companion-field attribution audit for the low-level
  `P` error, keeping perturbation-pressure `P` as the strict gate field;
- a selected-field spatial alignment/shift diagnostic audit for `U`, `V`, and
  `MU` end/evolution errors.

Both diagnostics may compare with WRF reference output for attribution, but
they must not generate candidates, patch fields, tune formulas with
reference-end truth, or describe diagnostic/audit/probe/hidden-seam artifacts
as gate passes.

## Round D58 Local Audit Direction

D57 commit `ee2fff4` (`Audit pressure bias and spatial alignment`) preserves
the same compatibility state after full validation: CMake build passed, CTest
passed `29/29`, pytest passed `160/160`, and push to `origin/main` succeeded.
This is codebase validation, not a WRF-compatible numerical pass. The strict
d02 `2025-07-26_00:10:00` gate remains failed, and no `00:20` progression is
compatible until that endpoint passes.

The D57 pressure vertical-bias audit keeps the error attribution on
perturbation-pressure `P`. The worst levels remain low mass levels `0..4`;
level `0` has `P` mean difference `-4167.791 Pa` and RMSE `4170.444 Pa`.
Companion base-state fields do not explain that bias: `PB` and `MUB`
companion differences are about `+1 Pa`, `MU` is about `+11.625`, full
pressure `P + PB` inherits the negative `P` bias, and `PSFC` is about
`-376.164 Pa`. The current compatibility conclusion is therefore a
perturbation-`P` producer or staging issue, not a base-state companion mismatch.

The D57 selected-field spatial alignment audit also remains diagnostic-only.
Small horizontal shifts improve global normalized RMSE only modestly: `U` end
`2.46%` and evolution `2.64%`, `V` end `3.40%` and evolution `3.79%`, and
`MU` end `7.57%` and evolution `7.28%`. The selected-field gap is therefore
not explained by a simple small horizontal shift.

D58 compatibility work should stay observation-only and audit local evidence:

- pressure producer/staging C++ and metadata paths that explain how
  perturbation `P` is sourced, staged, and written;
- selected-field pipeline, timing, and remap C++/metadata paths, cross-checked
  against the D56/D57 JSON summaries.

These audits may inspect code, metadata, and prior JSON reports for error
attribution. They must not patch candidate fields, generate candidate files,
tune formulas from reference-end truth, or describe diagnostic/audit/probe/
hidden-seam outputs as gate passes.

## Round D59 Formula and Timeline Audit Direction

D58 commit `5d9adb5` (`Audit pressure and selected-field producers`) preserves
the same compatibility state after full validation: CMake build passed, CTest
passed `29/29`, pytest passed `167/167`, and push to `origin/main` succeeded.
This is codebase validation, not a WRF-compatible numerical pass. The strict
d02 `2025-07-26_00:10:00` gate remains failed, and no `00:20` progression is
compatible until that endpoint passes.

The D58 pressure producer audit
`build/validation/r58_pressure_refresh_producer_audit.json` remains
diagnostic-only. It localizes the perturbation-pressure producer inventory
around `compute_krosa_pressure` (`src/dynamics/pressure_refresh.cpp:275`),
`refresh_krosa_moving_nest_pressure`
(`src/dynamics/pressure_refresh.cpp:357`), and
`apply_krosa_moving_nest_pressure_refresh_hook`
(`src/dynamics/pressure_refresh_hook.cpp:309`). The worst level remains mass
level `0`, and the maximum absolute `P` mean difference is about
`4167.790767 Pa`. Compatibility work should therefore inspect the pressure
formula inputs, constants, masks, and staging handoff before any producer
change is attempted.

The D58 selected-field pipeline audit
`build/validation/r58_selected_field_pipeline_audit.json` is also
diagnostic-only. It reported seven risk flags: target fraction below half,
target error fraction below half, amplitude near one while RMSE remains high,
modest best-shift improvement, end/evolution shift mismatch, gate/integrator
claims while audits fail, and a large movement-delta schedule/remap risk. The
movement child delta is `{i:60,j:35}`, and target fractions are about
`U = 0.407583`, `V = 0.407583`, and `MU = 0.404762`.

D59 compatibility work should remain observation-only:

- audit pressure formula inputs, constants, masks, target staging, and local
  producer metadata without writing candidate fields;
- audit the selected-field schedule/remap timeline around movement, target
  mask exposure, interpolation, pressure refresh, and output metadata.

These audits may inspect code, metadata, local inputs, and prior JSON reports
for attribution. They must not generate candidates, patch candidate fields,
tune formulas from reference-end truth, use diagnostic/probe/hidden-seam output
as gate evidence, lower d02 below `2 km`, or introduce best-track nudging.

## Round D60 Formula-Term and Runtime Timeline Direction

D59 commit `957c3d6` (`Audit pressure inputs and selected-field timeline`)
preserves the same compatibility state after full validation: CMake build
passed, CTest passed `29/29`, pytest passed `175/175`, and push to
`origin/main` succeeded. This remains codebase validation, not a
WRF-compatible numerical pass. The strict d02 `2025-07-26_00:10:00` gate
remains failed, and no `00:20` progression is compatible until that endpoint
passes.

The D59 pressure formula/input audit
`build/validation/r59_pressure_formula_inputs_audit.json` remains
diagnostic-only. It reported status `computed_with_flags` with six risk flags.
The target-region normalized RMSE for `P` is `10.335362`, with target mean
difference `-805.058863 Pa`. The low-level pressure error remains concentrated
at mass level `0`, where `P` mean difference is `-4167.790767 Pa` and RMSE is
`4170.444038 Pa`; `P + PB` at the same level inherits almost the same bias with
mean difference `-4166.786671 Pa`. Companion fields are not large enough to
explain this by themselves: target normalized RMSE is `8.76e-05` for `PB`,
`0.000878` for `MU + MUB`, `0.000805` for `PH + PHB`, `0.012115` for `T`, and
`0.016169` for `QVAPOR`. `HGT` target normalized RMSE is `0.495197`, but its
mean difference is only `-0.085612 m`, so this remains diagnostic evidence
rather than a candidate-repair path.

The D59 selected-field timeline audit
`build/validation/r59_selected_field_timeline_audit.json` is also
diagnostic-only. It reported status `computed_with_flags` with six risk flags.
The movement child delta is `{i:60,j:35}` and is classified as a large
movement. `P` is not listed as parent-interpolated, the changed/interpolated
count mismatch flag is true, and target fractions remain about `0.407583` for
`U` and `V`, and `0.404762` for `T`, `PH`, `MU`, `P`, and `QVAPOR`.

D60 compatibility work should remain observation/telemetry-only:

- add a per-column pressure formula-term probe for worst low-level target
  columns, isolating terms such as `theta/T`, `PH + PHB`, `MU + MUB`,
  coefficients, `ALB`, and `PB` subtraction without patching fields;
- emit selected-field runtime schedule/timeline metadata around movement,
  parent fill, exposed-cell interpolation, static refresh, pressure refresh,
  and output timestamps/counts.

D60 telemetry may explain producer timing and formula-term attribution, but it
must not patch fields, tune formulas from reference-end truth, generate
candidate files from oracle data, use diagnostic/probe/hidden-seam output as
gate evidence, lower d02 below `2 km`, or introduce best-track nudging.

## Round D61 Runtime Observation Direction

D60 commit `7903cbf` (`Add pressure column probe and timeline metadata`)
preserves the same compatibility state after full validation: CMake build
passed, CTest passed `29/29`, pytest passed `178/178`, and push to
`origin/main` succeeded. This remains codebase validation, not a
WRF-compatible numerical pass. The strict d02 `2025-07-26_00:10:00` gate
remains failed, and no `00:20` progression is compatible until that endpoint
passes.

The D60 normal `--pressure-refresh` KROSA candidate output is
`build/validation/r60_pressure_normal/wrfout_d02_2025-07-26_00:10:00`. It emits
runtime timeline metadata with 11 ordered events:
`cycle_start`, `move_from_to_parent_start`, `overlap_remap`,
`exchange_plan_build`, `parent_interpolation`,
`selected_field_change_summary`, `static_refresh`,
`pressure_refresh_readiness`, `pressure_refresh_apply`, `cycle_end`, and
`output_write_preparation`. The metadata records child movement delta `(60,35)`
and pressure-refresh refreshed, synced, and changed counts. The added telemetry
does not change normal candidate numerics: `U`, `V`, `T`, `PH`, `MU`, `P`,
`QVAPOR`, `PB`, `PHB`, `MUB`, and `HGT` match the D54 normal candidate with
maximum absolute difference `0.0`.

The D60 pressure column probe
`build/validation/r60_pressure_formula_column_probe_audit.json` remains
diagnostic-only. It selected worst columns `(160,49)`, `(160,50)`, `(160,51)`,
`(161,49)`, and `(161,50)`. For rank 1 at level `0`, candidate `P` is
`-4172.158691 Pa`, reference `P` is `514.484375 Pa`, and the difference is
`-4686.643066 Pa`; the mean absolute low-level `P` error across levels `0..4`
is `4391.227588 Pa`. Unique risk codes include
`low_level_column_p_error_large`,
`column_candidate_source_delta_matches_p_error`,
`coefficient_terms_missing_from_candidate_netcdf`, and
`formula_input_json_prior_risks_present`.

The D60 selected-field timeline audit
`build/validation/r60_selected_field_timeline_audit.json` is also
diagnostic-only. It confirms that runtime timeline attributes are present, but
it still reports the D59 risks until the audit parser fully structures the
event strings. The strict `00:10` gate report
`build/validation/r60_pressure_normal_gate.json` still fails: first failing
field `U` has normalized RMSE `0.117875`, with `V = 0.134244`,
`MU = 0.133382`, `P = 6.334952`, and storm-center error `43.482716 km`.

D61 compatibility work should remain telemetry/audit-only:

- add optional same-column runtime observation metadata for pressure formula
  attribution, without patching fields or tuning formulas from reference-end
  truth;
- parse selected-field runtime timeline attributes structurally, so movement,
  remap, static refresh, pressure-refresh, and output-preparation events can be
  audited without becoming gate evidence.

D61 telemetry may explain when and where a mismatch enters the local runtime
path, but it must not patch fields, tune formulas from reference-end truth,
generate candidate files from oracle data, convert audit/probe/hidden-seam
outputs into gate evidence, lower d02 below `2 km`, or introduce best-track
nudging.

## Round D62 Pressure Formula Observation Direction

D61 local commit `ab68c51` (`Add runtime pressure probes and timeline audit`)
preserves the same compatibility state after local full validation: CMake build
passed, CTest passed `29/29`, pytest passed `183/183`, and `git diff --check`
passed. This commit is local only at D62 start. Push to `origin/main` is
blocked by GitHub SSH/network timeouts on both normal SSH port `22` and
SSH-over-443; `origin/main` remains at D60 commit `7903cbf`. Do not report D61
as pushed until a later successful push is verified.

The D61 real KROSA smoke output is
`build/validation/r61_pressure_normal/wrfout_d02_2025-07-26_00:10:00`. Its
runtime timeline metadata contains 12 ordered events, with
`pressure_column_probe` inserted between `pressure_refresh_apply` and
`cycle_end`. The pressure-column probe observes columns `(160,49)`,
`(160,50)`, `(160,51)`, `(161,49)`, and `(161,50)` at levels `0..4`, phases
`post_static_refresh` and `post_pressure_refresh`, for a record count of `50`.
The probe still marks `ALB`, `C3F`, `C4F`, `C3H`, `C4H`, `P_TOP`, and
`theta_m` as unavailable in the candidate NetCDF. D61 telemetry does not change
accepted numerical fields relative to D60: `U`, `V`, `T`, `PH`, `MU`, `P`,
`QVAPOR`, `PB`, `PHB`, `MUB`, and `HGT` all have maximum absolute difference
`0.0`.

The strict d02 `2025-07-26_00:10:00` gate still fails for D61. The first
failing field remains `U` with normalized RMSE `0.117875`; related failures
include `V = 0.134244`, `MU = 0.133382`, `P = 6.334952`, and storm-center
error `43.482716 km`. Validation must not advance to `00:20`.

D62 compatibility work should remain diagnostic-only:

- add pressure-core formula observation for runtime POD values such as `ALB`,
  `C3F`, `C4F`, `C3H`, `C4H`, `P_TOP`, `theta`, and related pressure-refresh
  intermediates at requested columns/levels;
- add a Python pressure-column runtime probe audit that parses the D61/D62
  NetCDF attributes, structures before/after deltas, and flags strong negative
  post-refresh perturbation pressure and missing formula terms.

D62 observation may expose internal terms that were unavailable in D61, but it
must not patch candidate fields, tune formulas from reference-end truth,
generate oracle candidates, treat audit/probe/hidden-seam output as gate
evidence, advance to `00:20`, lower d02 below `2 km`, or introduce best-track
nudging.

## Round D63 Formula Observation Integration Direction

D62 is now integrated and synchronized. Commit `655df5b`
(`Add pressure formula observation audit`) has been pushed to `origin/main`,
and the same successful push also delivered the previously blocked D61 commit
`ab68c51` (`Add runtime pressure probes and timeline audit`). The D62 full
validation passed CTest `29/29`, pytest `188/188`, and `git diff --check`.

The D62 runtime probe audit
`build/validation/r62_pressure_column_runtime_probe_audit.json` is
diagnostic-only. It parsed `50/50` pressure-column records consistently. The
active risk flags are `post_pressure_refresh_p_negative`,
`large_p_drop_magnitude`, and `formula_terms_unavailable`. The minimum
post-refresh `P` is `-4174.87305 Pa`, and the largest pressure drop is at
`(161,49,k=0)` with `delta_P = -4257.17773 Pa`.

D63 compatibility work should connect the D62 pressure-core formula
observation to the selected-field `--pressure-column-probe` path. The intended
runtime output is compact NetCDF/JSON metadata for the requested
columns/levels, including the core formula observation terms produced by the
pressure refresh. The paired audit work should parse the new formula
observation attributes and associate them with the existing before/after
pressure-column deltas.

D63 remains diagnostic-only. It must not patch candidate fields, tune formulas
from reference-end truth, generate oracle candidates, promote audit/probe/
hidden-seam output to gate evidence, advance to `00:20`, reduce d02 below
`2 km`, or introduce best-track nudging.

D63 is now complete and synchronized. Commit `3e1f6a7`
(`Add selected-field formula observation telemetry`) has been pushed to
`origin/main`. Full code validation passed CTest `29/29`, pytest `195/195`,
and `git diff --check`.

The D63 real KROSA candidate is
`build/validation/r63_pressure_formula_observation/wrfout_d02_2025-07-26_00:10:00`,
with audit report
`build/validation/r63_pressure_formula_observation/runtime_probe_audit.json`.
It is telemetry/evidence-only, not a compatibility or validation-gate pass.
The audit shows formula observation is present, with `25` formula records, all
recorded and valid, and `0` pressure mismatches against the post-refresh probe
`P`. Therefore the negative perturbation pressure is now traced to the recorded
formula/base-pressure terms themselves, not to a NetCDF encoding problem or an
audit/probe mismatch.

The strict d02 `00:10` gate still fails for the D63 candidate. The first
failing field is `U` normalized RMSE `0.117875`; additional failures are `V`
`0.134244`, `MU` `0.133382`, `P` `6.334952`, and storm-center error
`43.482716 km`. D63 candidate fields `U`, `V`, `T`, `PH`, `MU`, `P`,
`QVAPOR`, `PB`, `PHB`, `MUB`, and `HGT` are numerically identical to D61, with
maximum absolute difference `0.0`.

D64 is complete and synchronized. Commit `2db9279`
(`Add pressure budget runtime audit`) has been pushed to `origin/main`, and
full validation passed CTest `29/29`, pytest `196/196`, and
`git diff --check`. This is diagnostic-only completion; it is not a strict d02
gate pass.

The real audit
`build/validation/r64_pressure_budget_runtime_probe_audit.json` on the D63
candidate has `status = computed_with_flags`. It records `25` pressure budget
records, `25` `total_pressure < PB` records, `25` records where the large drop
is explained by formula/base-pressure subtraction, and pressure mismatch count
`0`. The first tracked point `(i=160,j=49,k=0)` reports
`total_pressure = 95539.724 Pa`, `PB = 99711.8828 Pa`,
`total_pressure_minus_pb = -4172.1588 Pa`, and
`probe_delta_p = -4252.3540025 Pa`.

The source audit confirms that both TyWRF and WRF runtime use
`P = total_pressure - PB`. Compatibility work must not remove the `PB`
subtraction to hide negative perturbation pressure. The next pressure
diagnostics should instead focus on `PH + PHB`, `theta`, `MU + MUB`, and
base-state staging consistency.

D65 compatibility direction is diagnostic-only: add a formula sensitivity
diagnostic and audit WRF moving-nest/base-state call order before any formula
correction. D65 must not patch formulas, generate reference-end or oracle
candidates, use audit/probe/hidden-seam output as gate evidence, advance to
`00:20`, reduce d02 below `2 km`, or introduce best-track nudging.

D65 is complete and synchronized, but remains diagnostic-only rather than a
compatibility gate pass. Commit `07908b1`
(`Add pressure formula sensitivity audit`) has been pushed to `origin/main`,
with full validation passing CTest `29/29`, pytest `197/197`, and
`git diff --check`.

The real sensitivity audit
`build/validation/r65_pressure_formula_sensitivity_audit.json` records `25`
sensitivity records. The `PB - total_pressure` gap ranges from `3593.9011` to
`4174.8731 Pa`, with mean `3904.1623 Pa`. The approximate fractional
total-pressure increase needed ranges from `0.0389068` to `0.0436995`, with
mean `0.0414497`; all `25` records require a large total-pressure increase.
This explains the scale of the pressure mismatch but does not authorize a
formula patch or validation shortcut. The strict d02 `00:10` gate still fails,
so no `00:20` progression is compatible.

B65 concluded that WRF's broad base-state interpolation path includes
`PHB`, `MUB`, `PB`, `ALB`, `T_INIT`, and `HT`, while TyWRF's selected-field path
currently interpolates only `U`, `V`, `T`, `PH`, `MU`, and `QVAPOR`. D66
compatibility direction is therefore to expose a d02 base-field provenance
contract and audit for newly exposed cells, including WRF generated
interpolation-mask semantics. This is a contract/audit step only: it must not
open the selected-field numerical path yet, generate reference-end or oracle
candidates, treat audit/probe/hidden-seam output as gate evidence, reduce d02
below `2 km`, or introduce best-track nudging.

At D68 start, D66 commit `09d6ba2`
(`Add moving-nest base-state exchange contract`) and D67 commit `5a78fbb`
(`Add base-state exchange action diagnostics`) have both been pushed.
`main` and `origin/main` are synchronized at `5a78fbb`. The strict d02
`2025-07-26_00:10:00` gate still fails; no D66 or D67 diagnostic report is a
compatibility gate pass.

The D66 contract is diagnostic-only. It reports active selected fields
`U`, `V`, `MU`, `QVAPOR`, `T`, and `PH`, while the WRF broad base-state
candidate set is `PHB`, `MUB`, `PB`, `ALB`, `T_INIT`, `HT`, and `HGT`.
`PHB`, `MUB`, and `PB` are State-backed candidates; `ALB`, `T_INIT`, `HT`,
and `HGT` are static/provider-backed candidates. No base-state candidate is
marked as selected-field interpolated, `P` is excluded from that base-state
candidate list, `diagnostic_only = true`, and
`enables_selected_field_numerics = false`.

D66 also adds exposed-child exchange planning for the selected fields. The
plan describes newly exposed active child strips, their stagger, point counts,
and whether future parent interpolation is required. It does not interpolate,
does not mutate fields, and reports `performed_interpolation = false`,
`modifies_overlap = false`, and `modifies_halo = false`. D66 validation before
the local commit passed CTest `29/29` and pytest `197/197`.

B66 refined the WRF compatibility target: WRF's generated moving-nest
`imask_nostag` semantics apply to `PHB`, `T_INIT`, `MUB`, `ALB`, `PB`, and
`HT`, and WRF `start_domain` recompute rules govern which generated/base
fields are recomputed after nest movement. TyWRF must therefore distinguish
overlap remap, exposed-mask parent/base interpolation, static/provider
generation, and `start_domain` recompute provenance instead of treating these
fields as ordinary selected-field numerical outputs.

D67 completed an opt-in diagnostic moving-nest base-field provenance/action
report and prepared exposed-mask regression coverage. The report makes
base-field action and source state explicit for exposed d02 cells, and D68 is
the first helper/regression slice to pin WRF-style unstaggered exposed-mask
geometry. The action vocabulary is diagnostic-only:
`interpolate_exposed_cells` for `MUB`,
`recompute_from_mub_after_interpolation` for `PB`, `ALB`, and `T_INIT`,
`preserve_interpolated_when_rebalance_zero` for `PHB`, and
`static_height_input` for `HT`/`HGT`.

D68 compatibility direction is a diagnostic-only exposed base-state exchange
helper and regression slice. The helper/test may expose WRF-style
`PHB`/`MUB`/`HT` exposed interpolation semantics and mark
`PB`/`T_INIT`/`ALB` for recompute, but these are diagnostic/helper semantics,
not selected-field numerical production and not gate evidence. D68 must not
wire the helper into the production selected-field numerical path, generate
reference-end or oracle candidates, treat diagnostic/probe/helper output as a
gate pass, advance validation to `00:20`, reduce d02 below `2 km`, or
introduce best-track nudging. D70 or later may consider hook or diagnostic
wiring for WRF-style exposed base-state policy after this boundary is reviewed.

D68 is now complete, fully validated, pushed, and synchronized. Commit
`c8a83a2` (`Add exposed base-state exchange diagnostics`) records the exposed
base-state exchange diagnostic helper/regression slice. It remains
diagnostic-only: it does not connect to the production selected-field path and
does not make the strict d02 `2025-07-26_00:10:00` gate pass.

D69 is now complete, fully validated, pushed, and synchronized. Commit
`e32ccc9` (`Add exposed MUB base-state recompute provider`) records a
provider/test-only API that recomputes `PB`, `T_INIT`, and staged provider
`ALB` from already exposed-interpolated `MUB` for newly exposed d02 cells. The
D69 API is still not connected to the production selected-field path, is not a
candidate writer, and is not gate evidence. It must not regenerate `MUB` from
`HGT`, rebuild or synchronize `PHB`, write provider `T_INIT` into `State::t`,
or write `ALB` into `State`. `ALB` remains an external staging value for
downstream pressure-refresh consumers only. The strict d02
`2025-07-26_00:10:00` gate still fails.

D70 is now complete, fully validated, pushed, and synchronized. Commit
`5eb6485` (`Add exposed base-state adapter diagnostics`) records the
diagnostic-only adapter/report that composes the D68 exposed base-state
exchange helper with the D69 exposed-`MUB` recompute API. It remains outside
the normal `selected_field_cycle` path: it does not write a production
candidate, does not connect the pressure-refresh gate path, and does not change
production selected-field numerics. Its outputs remain report/dry-run/staging
artifacts with `diagnostic_only = true`, `gate_candidate = false`, and
`integrator_output = false`. The strict d02 `2025-07-26_00:10:00` gate still
fails.

D71 compatibility scope is metadata guard work first. It must not connect the
D70 adapter into production `selected_field_cycle`, even as an opt-in
candidate path, until the guard prevents diagnostic artifacts from masquerading
as normal gate candidates. C++ NetCDF attributes and stdout JSON should be
derived from one unified candidate-disposition decision, rather than separate
boolean logic that can drift between file metadata and reports.

The disposition guard must reject pseudo-positive artifacts by role and kind,
not only by the presence of positive-looking booleans. `helper`, `probe`,
`adapter`, `dry_run`, `staging`, `experimental`, diagnostic, oracle, and
reference-copy outputs remain non-gate even if a report or NetCDF file carries
`gate_candidate = true` or `integrator_output = true` by mistake. D71 and later
compatibility work must not use reference-end truth, oracle/reference-copy
candidates, diagnostic/probe/helper outputs, adapter/dry-run/staging reports,
or restart substitutes as compatibility evidence. It must not advance
validation to `00:20`, reduce d02 below `2 km`, or introduce best-track
nudging.

D71 is now complete, fully validated, pushed, and synchronized at commit
`9b08d09` (`Guard diagnostic adapter gate metadata`). The guard keeps
diagnostic adapter artifacts classified by role/kind before any positive
candidate booleans are considered, so pseudo-positive adapter metadata cannot
be promoted into a normal compatibility candidate. `main` and `origin/main`
are synchronized at this commit, and the strict d02
`2025-07-26_00:10:00` gate still fails.

D72 may wire an opt-in diagnostic adapter report path for the exposed
base-state adapter. That path is limited to reporting and diagnostic
staging/audit metadata. It must not be connected to the normal
pressure-refresh candidate path, must not write or relabel a production
`selected_field_cycle` candidate, and must not change strict-gate
eligibility. Any D72 adapter artifact must remain metadata non-gate:

```text
diagnostic_only = true
gate_candidate = false
integrator_output = false
```

The strict d02 gate must reject D72 diagnostic adapter output by default. Such
output is not a model pass, not pressure-refresh gate evidence, and not
validated integrator output. D72 must not use reference-end truth, oracle or
reference-copy products, diagnostic/probe/helper output, adapter reports, or
restart substitutes as compatibility evidence; it must not advance validation
to `00:20`, reduce d02 below `2 km`, or introduce best-track nudging.

D72 is now complete, verified, pushed, and synchronized at commit `616c6c9`
(`Add diagnostic adapter selected-field report path`). It adds the opt-in
selected-field diagnostic adapter report path behind
`--diagnostic-base-state-adapter-report`. This is staging/report-only wiring
for the exposed base-state adapter: it must not write or relabel a production
`selected_field_cycle` candidate, must not connect the normal pressure-refresh
candidate path, and must not change selected-field gate eligibility.

Any D72 adapter report or artifact remains non-gate compatibility context:

```text
diagnostic_only = true
gate_candidate = false
integrator_output = false
```

D73 is now complete, verified, pushed, and synchronized at commit `b73704f`.
It kept the opt-in adapter smoke/audit boundary diagnostic-only: the strict d02
gate rejects adapter, helper, staging, dry-run, probe, diagnostic, oracle, and
reference-copy roles/kinds before interpreting field RMSE or TC diagnostics.
Even if an adapter-side report shows passing field thresholds, storm-center,
MSLP, or Vmax values, the artifact is still not a gate candidate and cannot
count as the strict `2025-07-26_00:10:00` pass. The current strict gate remains
failed, so validation must stop at `00:10` and must not advance to `00:20`.

D74 is now complete, verified, pushed, and synchronized at commit `b1d2de2`
(`Add base-state source staging provider`). It added the non-oracle
child-shaped `BaseStateSourceStagingProvider` as a diagnostic/test-only source
contract intended to replace `source==child` staging in a later hidden adapter
wiring round. It is not connected to `selected_field_cycle`, does not write
normal candidates, does not mark anything as gate-eligible, and does not relax
strict-gate metadata or disposition checks. Any D74 report, fixture, or staged
buffer remains non-gate compatibility context only:

```text
diagnostic_only = true
gate_candidate = false
integrator_output = false
```

D75 is now complete, verified, pushed, and synchronized at commit `d5a1f99`
(`Wire diagnostic adapter source staging provider`). It wired the D74 provider
only into the hidden diagnostic adapter source-staging side. That path remains
report/staging-only: it does not write or relabel normal candidates, refresh
candidate `P`, change the normal pressure-refresh path, or relax strict gate
metadata. Provider-backed source provenance and staged/count metadata explain
which source populated hidden scratch buffers, but that metadata is itself
disqualifying context and cannot be used as strict-gate evidence.

D76 is complete, verified, and pushed at commit `f7082c1`. It added
source-vs-child delta/provenance planning for hidden base-state adapter source
staging only. Diagnostic deltas, zero-difference summaries, source/child
provenance, and staged/exposed/masked counts for the hidden adapter scratch path
remain non-gate metadata. The presence of any source-child delta or provenance
metadata is disqualifying for strict compatibility evidence, even when it shows
zero difference or useful counts. D76 did not write normal candidates, refresh
candidate `P`, change the normal pressure-refresh path, relax strict gate
metadata, advance validation to `00:20`, reduce d02 below `2 km`, or introduce
best-track nudging. The strict d02 `2025-07-26_00:10:00` gate remains failed.

D77 is complete, verified, pushed, and synchronized at commit `c7e46a7`. It
kept compatibility scope inside the hidden diagnostic provider-derived
base-state source path. The hidden adapter source now uses provider-derived
`PB`, `T_INIT`, `MUB`, `ALB`, and `PHB` plus `output_static` `HGT`/`HT`, with
provider-source metadata and strict-gate regression coverage. This does not
change normal selected-field candidate eligibility or strict-gate eligibility.

Provider-derived source metadata, source-staging metadata, and source-child
delta summaries remain diagnostic-only and disqualifying context. They cannot
be used as strict gate evidence, cannot prove a `2025-07-26_00:10:00` pass,
cannot borrow positive selected-field metadata, and cannot make a normal
candidate compatible by report wording. Any `PHB` reconstruction in this
provider-derived hidden source path remains diagnostic evidence about provider
staging only, not WRF rebalance semantics and not proof that candidate `PHB`,
`P`, or pressure refresh is physically produced.

D78 is complete, verified, pushed, and synchronized at commit `d7ecd11`. It
added the diagnostic adapter provider-source audit CLI for provider-source,
source-staging, and source-child-delta metadata. The audit path is
diagnostic-only and fail-closed: it may summarize hidden diagnostic attributes
and identify the largest deltas for debugging, but it must not write candidate
fields, refresh candidate `P`, change pressure-refresh behavior, relax strict
gate metadata, generate oracle candidates, advance validation to `00:20`,
reduce d02 below `2 km`, or introduce best-track nudging. D78 audit output is
not strict-gate evidence and is not proof that the failed `00:10` endpoint
passed.

D79 may add a normal-path selected-field exposed base-state producer, but its
source boundary is narrower than the D77 hidden diagnostic adapter path. The
producer may write only legitimate candidate base-state fields derived from the
current cycle-start state, d01-derived parent data, and same-time KROSA
metadata/constants. It must not use D77 hidden diagnostic adapter NetCDF
attributes, reference-end truth, WRF end-state deltas, a direct `P` shortcut,
gate-metadata relaxation, a `PSFC`-as-`SLP` proxy, or promotion of diagnostic
artifacts into compatibility evidence.

D79 must not directly patch perturbation `P`; pressure refresh remains a
separately guarded producer. It also must not claim the strict
`2025-07-26_00:10:00` gate passed or proceed to `00:20`; the current first
failure remains `U` normalized RMSE `0.117875` until a real strict-gate output
proves otherwise. Any `PHB` reconstruction in this normal-path base-state
producer remains limited by the documented provider semantics and must not be
described as WRF rebalance equivalence.

D79 is complete, passed focused and full validation, and was pushed at commit
`3d49d9c` (`Add normal selected-field base-state producer`). That commit records
the normal selected-field base-state producer as candidate-path compatibility
work, not as proof that the field gate passed. The real KROSA d02
`2025-07-26_00:10:00` gate still fails, with first failed field `U` normalized
RMSE `0.11787539215928292`; validation must not advance to `00:20`.

D80 may refresh only the normal pressure-refresh production metadata naming:
legacy helper/dry_run/staging blocker labels may move to gate-safe
production/readiness/source-sync names when they describe the same normal
candidate path. This is a metadata eligibility clarification only. It must not
relax the strict metadata guard and must not allow diagnostic, oracle, helper,
probe, adapter, staging, dry_run, or experimental artifacts to pass as
WRF-compatible integrator output. The metadata repair must keep the real
numerical failures visible: failed `U`, `V`, `MU`, `P` RMSE and the failed TC
center diagnostic remain blockers until a real strict `00:10` gate passes.

D80 is complete, passed focused and full validation, passed the real `00:10`
metadata check with `candidate_metadata` accepted, and was pushed at commit
`765bd06` (`Make pressure-refresh production metadata gate-safe`). This commit
records a metadata-eligibility repair for the normal pressure-refresh
production path only. It does not prove numerical compatibility: the real d02
`2025-07-26_00:10:00` gate still fails overall, with first failed field `U`
normalized RMSE `0.11787539215928292`. Validation must remain stopped at
`00:10` and must not advance to `00:20`.

D81 wind-failure audit scope is read-only compatibility diagnosis. It may read
the reference and candidate `U`/`V` fields and decompose the observed wind
errors by region, level, staggered layout, remap/exposed-cell mask, and
metadata provenance. Its output must remain diagnostic-only with
`gate_evidence=false`: it must not write candidate files or candidate fields,
change thresholds, relax metadata rejection, or be described as proof that the
`00:10` gate passed. D81 must not handle `P`, `MU`, physics, or best-track
nudging, and d02 must remain `2 km`.

D81 is complete at commit `3b17b6b`
(`Add selected-field wind error audit`). It passed focused validation, full
validation, the real D80 wind-audit smoke, and was pushed. The audit conclusion
is that the D80 real `00:10` selected-field candidate has domain-wide `U`/`V`
error, not a metadata-only problem. `candidate_metadata` passes, but the field
gate remains failed; the first failed field is still `U` with normalized RMSE
`0.11787539215928292`. Compatibility work must remain stopped at `00:10` and
must not advance to `00:20`.

D82 may address only parent-child interpolation `U`/`V` stagger-aware
coordinate semantics and directly related tests or metadata. The allowed
compatibility boundary is coordinate interpretation for staggered wind fields:
`U` must be handled on the west-east staggered coordinate, `V` on the
south-north staggered coordinate, and any parent-child coordinate metadata must
make that distinction explicit. D82 must not handle `P`, `MU`, physics,
best-track nudging, oracle/reference-copy candidates, or any d02 resolution
change below `2 km`.

Coordinate metadata, selected-field wind audits, and stagger-coordinate reports
are compatibility diagnostics only. They are not gate evidence and must not be
reported as an accepted `00:10` field pass. The only evidence that D82 improved
compatibility is a later real strict `2025-07-26_00:10:00` d02 gate report with
positive integrator metadata and improved field metrics.

D82 is complete, passed focused validation, full validation, and the real
`00:10` gate smoke, and was pushed at commit `0a192d4`
(`Use WRF-style staggered interpolation coordinates`). The real
`2025-07-26_00:10:00` d02 run now passes candidate metadata, and the
staggered-coordinate change gives finite `U`/`V` improvement, but the strict
field gate still fails numerically. The first failed field remains `U` with
normalized RMSE `0.11359509344276145`; `V` normalized RMSE is
`0.12937874064226143`. `P` normalized RMSE is slightly worse at
`6.354144377000247`, so compatibility notes must keep `P` regression risk
visible rather than describing D82 as a pass.

D83 documentation and audit scope is limited to `U`/`V` source/time-level
attribution at the failed `00:10` endpoint. Any D83 output is diagnostic-only
and must carry `gate_evidence=false`. It may explain which source and time
level contributed to `U`/`V` differences, but it must not write a candidate,
alter thresholds, repair metadata, or serve as pass proof. D83 must not handle
`P`, `MU`, physics, or best-track nudging, must not lower d02 below `2 km`, and
must not advance validation to `00:20`.

D83 is complete, passed focused validation, full validation, and the real D82
source-attribution smoke, and was pushed at commit `d8b8453`
(`Add selected-field wind source attribution audit`). Its compatibility result
is diagnostic attribution only. The moving-nest raw d02 start comparison is
classified as raw-pose diagnostics, not shifted-start persistence evidence and
not a validation-gate result. The strict real d02 `2025-07-26_00:10:00` gate
still fails, so validation remains stopped at `00:10` and must not advance to
`00:20`.

D84 may add only a `U`/`V` wind tendency core skeleton and directly related
tests. This is scaffolding for wind-tendency accounting, not a declaration of
real-gate improvement. D84 must not handle `P`, `MU`, physics, best-track
nudging, oracle/reference-copy candidates, or any d02 resolution change below
`2 km`. Any wind-tendency metadata, source/core report, or skeleton output is
compatibility context only and is not gate evidence; only a later real strict
`2025-07-26_00:10:00` gate run with positive integrator metadata and passing
metrics can prove compatibility improvement.

D84 is complete, validated, and pushed at commit `44c7fff`
(`Add wind tendency core skeleton`). The commit adds the wind tendency core
skeleton and reports only skeleton/accounting state. Its wind tendency report
is not validation-gate evidence; the core report explicitly fixes
`validation_gate_evidence=false` so reports cannot be mistaken for a strict
gate result. The prior real strict d02 `2025-07-26_00:10:00` gate still fails,
validation remains stopped at `00:10`, and d02 remains `2 km`.

D85 may add only default-off selected-field `U`/`V` wind tendency opt-in
plumbing. The opt-in path may expose wiring, source selection metadata, and
strict rejection behavior for selected wind tendencies, but it must not become
default behavior and must not claim numerical validation improvement. Zero or
identity placeholder tendency sources are allowed only as strict-gate-negative
plumbing checks: they must be rejected as non-evidence and must not receive
gate credit. D85 must not use any reference-end or oracle source, must not
touch `P`, `MU`, `PB`, or `PHB`, must not invoke pressure refresh, physics, or
best-track nudging, must not change d02 below `2 km`, and must not advance or
report validation at `00:20`.

D85 is complete, validated, and pushed at commit `18bf109`
(`Add selected-field wind tendency opt-in plumbing`). The default selected-field
wind tendency source remains `none`; that default emits no wind tendency
metadata and makes no `U`/`V` tendency update. The `zero` and `identity`
sources are placeholder wiring modes only. They are non-evidence, must be
strict-gate rejected, and must not be described as wind compatibility progress.

D86 may add only a selected-field `U`/`V` non-oracle wind tendency source based
on self-advection/prognostic candidate state. If that source writes positive
candidate metadata, the metadata may make the artifact gate-eligible, but it is
not pass proof by itself. The only acceptable proof is a real strict
`2025-07-26_00:10:00` d02 field and diagnostic gate success from that
non-oracle candidate. D86 must not use reference-end or oracle sources, must
not touch `P`, `MU`, `PB`, or `PHB`, must not invoke pressure refresh,
physics, or best-track nudging, must not change d02 below `2 km`, and must not
advance or report validation at `00:20`.

D86 is complete, validated, and pushed at commit `83162d3`
(`Add self-advection wind tendency source`). The `self_advection` source is a
non-oracle selected-field `U`/`V` tendency source and its real KROSA output
passes the `candidate_metadata` eligibility check, but that metadata result is
not a field pass. The real strict d02 `2025-07-26_00:10:00` gate still fails
first on `U`, so compatibility work remains stopped at `00:10` and must not
advance to `00:20`.

D87 may add only selected-field `U`/`V` `self_advection` subcycling. The
subcycling path must remain non-oracle and may use only candidate/start-state
wind tendency inputs already allowed for the selected-field path. It must not
use reference-end truth, WRF end-state deltas, oracle/reference-copy data, or
reference-end-derived tendencies. It must not touch `P`, `MU`, `PB`, or `PHB`,
must not run pressure refresh, must not invoke physics or best-track nudging,
and must not change d02 resolution below `2 km`.

D87 subcycling metadata may make an artifact gate-eligible only as
`candidate_metadata` from the normal non-oracle candidate path. Metadata,
substep counts, source reports, or tendency reports are not pass proof. The
only compatibility proof remains a real strict `2025-07-26_00:10:00` d02 run
whose fields and TC diagnostics pass; until then validation must not proceed to
`00:20`.

D87 is complete, verified, and pushed at commit `8466e7c`
(`Add self-advection wind subcycling`). It added
`--wind-tendency-substeps N` for the selected-field `self_advection` wind
tendency source. The d02 wind tendency substep is fixed at `8 s`; `N=75`
therefore represents the complete `2025-07-26_00:00:00 ->
2025-07-26_00:10:00` 600 second endpoint.

The real KROSA d02 `00:10` run with `N=75` passed `candidate_metadata`, but the
strict gate failed. The first failed field is still `U`, with normalized RMSE
`0.13578703428452885`; `V` normalized RMSE is `0.15517830284022266`. The SLP
diagnostic gate failed storm-center distance at `43.4827156063485 km`, while
minimum SLP error `0.36407470703125 hPa` and Vmax10m error
`0.7686817572680305 m s-1` passed. Compared with default wind tendency source
`none`, D87 `self_advection` subcycling changed only `U` and `V`; non-wind
fields were unchanged.

D88 compatibility scope is limited to wind subcycling `U`/`V` error
localization and audit tooling. D88 must not treat the D87 metadata pass as a
validation pass, must not advance or report `00:20`, and must not use
reference-end/oracle sources. It must not handle `P`, `MU`, `PB`, or `PHB`,
must not invoke pressure refresh, physics, or best-track nudging, and must not
change d02 resolution below `2 km`.

D88 is complete, verified, and pushed at commit `a785388`
(`Add wind subcycling audit tool`). It added
`tools/audit_wind_subcycling.py` and
`tests/unit/test_audit_wind_subcycling.py` for selected-field wind subcycling
audit coverage: `U`/`V` global metrics, boundary-band splits, vertical-level
peaks, baseline delta against the default wind-tendency path, and subcycling
metadata checks.

The D88 real audit reproduced the D87 KROSA d02 `N=75` wrfout result. The
subcycling metadata is `SUBSTEP_COUNT=75`, `SUBSTEP_DT_SECONDS=8`, and
`TOTAL_SECONDS=600`. The reproduced global normalized RMSE is `U=0.135787`
and `V=0.155178`; deep-interior `band=40` RMSE is `U=2.097462` and
`V=2.463738`; the vertical peak is zero-based `U k=44` with RMSE `2.292039`
and zero-based `V k=7` with RMSE `2.760373`. These are localization/audit
facts, not compatibility-pass evidence.

D89 may only document or run selected-field `U`/`V` `self_advection` A/B work:
frozen versus refreshed advecting velocity, plus read-only cross-component
feasibility checks. The default selected-field behavior must remain the D87/D88
refreshed advecting-velocity path. D89 must not treat `candidate_metadata` or
subcycling metadata as a validation pass; the strict
`2025-07-26_00:10:00` d02 gate still fails and validation must not advance to
`00:20`. D89 must not use reference-end/oracle sources, must not handle `P`,
`MU`, `PB`, or `PHB`, must not invoke pressure refresh, physics, or best-track
nudging, and must not change d02 resolution below `2 km`.

D89 is complete, verified, and pushed at commit `d7e54e9`
(`Add wind advecting velocity mode`). It added
`--wind-tendency-advecting-velocity refreshed|frozen` for selected-field
`self_advection`; the default remains `refreshed`, preserving the D87/D88
behavior.

The real KROSA d02 `00:10` A/B results are compatibility evidence for a failed
gate, not a pass. With refreshed advecting velocity and `N=75`,
`candidate_metadata` passed, but the gate failed with `U` normalized RMSE
`0.13578703428452885` and `V` normalized RMSE
`0.15517830284022266`. With frozen advecting velocity and `N=75`,
`candidate_metadata` also passed, but the gate failed with `U` normalized RMSE
`0.13553969712614714` and `V` normalized RMSE
`0.15933552812814994`. The frozen SLP gate failed storm-center distance at
`43.4827156063485 km`; minimum SLP error `0.36407470703125 hPa` passed, and
Vmax10m error `0.7686817572680305 m s-1` passed. Frozen advecting velocity is
therefore not a fix; the next wind compatibility step is cross-component A/B.

D90 is complete, verified, and pushed at commit `34e2213`
(`Add cross-component wind advecting mode`). It added
`--wind-tendency-advecting-components same_component|cross_component` for
selected-field `U`/`V` `self_advection`. The default `same_component` mode
preserves D89 exactly. D90 records
`TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS=same_component|cross_component`;
`TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION=same_grid` is used for
`same_component`, and `TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION=average` is
used for `cross_component`.

The D90 real KROSA d02 `00:10` results remain failed-gate evidence, not a
validation pass. Same-component/refreshed reported `U` normalized RMSE
`0.13578703428452885` and `V` normalized RMSE `0.15517830284022266`.
Same-component/frozen reported `U` normalized RMSE `0.13553969712614714` and
`V` normalized RMSE `0.15933552812814994`. Cross-component/refreshed reported
`U` normalized RMSE `0.1271909426315702` and `V` normalized RMSE
`0.13542621105084265`. Cross-component/frozen reported `U` normalized RMSE
`0.1262784410537921` and `V` normalized RMSE `0.1343751747688221`. The derived
SLP diagnostic still failed storm-center distance at `43.4827156063485 km`;
minimum SLP error and Vmax10m error passed. The strict `00:10` gate remains
failed and validation must not advance to `00:20`.

D91 is complete, verified, and pushed at commit `2a2de69`
(`Add wind advection form sensitivity`). Its scope remained residual
localization, audit metadata/delta improvements, and default-off
advection-form sensitivity. D91 must not be read as a validation pass, must not
advance to `00:20`, must not use reference-end/oracle data, must not reduce d02
below `2 km`, must not introduce best-track nudging, and must not couple this
wind sensitivity work to pressure refresh or physics.

D91 adds `--wind-tendency-advection-form centered|upwind` for selected-field
`self_advection` and keeps `centered` as the default. The centered
cross/frozen KROSA `00:10` candidate is value-identical to D90 cross/frozen for
all compared fields while writing `TYWRF_WIND_TENDENCY_ADVECTION_FORM=centered`.
The strict gate remains failed: centered cross/frozen reports `U` normalized
RMSE `0.1262784410537921` and `V` normalized RMSE `0.1343751747688221`;
cross/frozen/upwind reports `U` normalized RMSE `0.16812017042527477` and `V`
normalized RMSE `0.17567799996961442`; cross/refreshed/upwind reports `U`
normalized RMSE `0.16408002949211245` and `V` normalized RMSE
`0.17946073108932759`. This upwind sensitivity degrades the current
selected-field wind-only result, so the next compatibility step should move
toward tendency decomposition and pressure-gradient/mass-coupled wind
experiments rather than declaring an advection-form closure.

The D91 wind audit extension also records advection-form metadata, including
`TYWRF_WIND_TENDENCY_ADVECTION_FORM`, and reports explicit baseline-delta
statistics for candidate-versus-baseline comparisons. These diagnostics support
residual localization only; they do not change the strict gate result.

D92 scope is tendency-decomposition audit plus a standalone, default-off
pressure-gradient wind-tendency staging skeleton. D92 must not claim a gate
pass, must not advance to `00:20`, must not use reference-end/oracle data for
model updates, must not reduce d02 below `2 km`, must not introduce best-track
nudging, and must not couple the staging skeleton to physics or pressure
refresh in this round.

The D92 wind tendency decomposition audit contract is diagnostic-only: reports
must set `diagnostic_only=true`, `uses_reference_end_truth=true`, and
`advances_00_20=false`. The residual is
`(candidate_end-candidate_start)-(reference_end-reference_start)`, so
reference-end truth is allowed only for post-run audit decomposition, never for
model updates. `--child-delta` is in child-grid cells; the
`--from-parent-start/--to-parent-start` difference is an unscaled parent index
delta, so D90/D91 parent-ratio-scaled decompositions must pass the scaled width
directly with `--child-delta`. This audit is not gate evidence and cannot
justify `00:20` progression.

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
