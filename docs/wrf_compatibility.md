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
