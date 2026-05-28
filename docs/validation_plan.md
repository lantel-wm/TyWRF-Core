# Validation Plan

Validation starts with progressive 10 min gates against the 1 h KROSA reference
run. Longer 3 h or 6 h cycle tests are deferred until the 10 min sequence is
stable.

## Cycle Procedure

1. Start from the matching WRF reference initial state.
2. Integrate the next TyWRF-Core segment, starting with `00:00 -> 00:10`.
3. Write WRF-compatible core `wrfout` files.
4. Compare against the WRF reference output for the same valid time.
5. Stop at the first failed 10 min endpoint and fix that field/kernel before
   advancing to the next endpoint.

## d02 10 min Strict Gate

The current acceptance gate is intentionally narrower than the full validation
plan. For each d02 cycle-end `wrfout`, starting at
`2025-07-26_00:10:00`, it requires:

- `U`, `V`, `T`, `PH`, `MU`, `P`, `QVAPOR`: normalized RMSE <= `5%`
- d02 storm center error <= `20 km`
- minimum SLP error <= `5 hPa`
- Vmax10m error <= `5 m s-1`

Missing candidate files and missing TC pressure diagnostics are hard failures
reported as `not_available`; they are never treated as passing.
The gate requires a real `SLP`, `MSLP`, `AFWA_MSLP`, `PMSL`, `PRMSL`, or
`SEA_LEVEL_PRESSURE` variable for the storm-center and minimum-SLP diagnostics.
It does not pass the minimum-SLP gate from a `PSFC` proxy.

Longer cycle checks use the same gate interface once the 10 min sequence is
stable:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/cycle_gate.py \
  --reference-dir /path/to/reference/WRF \
  --candidate-dir /path/to/tywrf/output \
  --start 2025-07-26_00:00:00 \
  --end 2025-07-26_06:00:00 \
  --pretty
```

`--hours 6` can be used instead of `--end`. The default domain is `d02`; use
`--interval` to gate multiple cycle endpoints between `--start` and `--end`.
The JSON report contains each cycle, each field gate, each diagnostic gate, and
the top-level `passed` or `failed` status.

For the current 10 min KROSA smoke, use the 1 h / 10 min reference directory
and force the interval in minutes:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/cycle_gate.py \
  --reference-dir /home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF_1h_10min_20260527_172838 \
  --candidate-dir /path/to/tywrf/output \
  --start 2025-07-26_00:00:00 \
  --end 2025-07-26_00:10:00 \
  --domain d02 \
  --interval-minutes 10 \
  --pretty
```

`tywrf_skeleton_cycle` can generate a d02 persistence/skeleton candidate for
this smoke by reading the `00:00` d02 state/template and writing
`wrfout_d02_2025-07-26_00:10:00` with `--times 2025-07-26_00:10:00`. The output
is explicitly marked `not_physical`, `integrator_output=false`, and
`validation_gate_only=true`. Default strict `cycle_gate.py` behavior must reject
that candidate at `00:10`; this is a wiring smoke and is not a passing
integrator result. The skeleton metadata/report records `--interval-minutes 10`
when `--cycle-start` and `--cycle-end` imply a 10 min segment. d02 `DX` and `DY`
must remain `2000 m`.

Diagnostic moving-nest parent-fill outputs and 10-minute diagnostic reports are
also non-physical and non-gate, even when exposed cells have direct parent-fill
values. The current staged sequence is overlap remap, direct parent-fill fields,
post-start-domain time-level sync, then required `P` refresh. Because direct
parent-fill intentionally excludes exposed-cell `P`, and pressure-refresh
compute is not invoked from skeleton/remap yet, these reports must not relax the
d02 strict gate or count as a validation pass.

## Broader Field Thresholds

These thresholds remain the broader validation target for later reports:

- `U`, `V`, `T`, `PH`, `MU`, `P`, `QVAPOR`: normalized RMSE <= `5%`
- `W` and hydrometeors: normalized RMSE <= `10%`
- `PSFC`, `U10`, `V10`, `T2`, `Q2`: normalized RMSE <= `10%`

## Typhoon Diagnostics

- d02 storm center error <= `20 km`
- minimum SLP error <= `5 hPa` when a sea-level pressure field is available
- Vmax error <= `5 m s-1`

The TC diagnostics prefer a real sea-level pressure field when present. The
current accepted candidate names are `SLP`, `slp`, `MSLP`, `mslp`,
`AFWA_MSLP`, `afwa_mslp`, `PMSL`, `pmsl`, `PRMSL`, `prmsl`,
`SEA_LEVEL_PRESSURE`, and `sea_level_pressure`. If one is found, diagnostics
compute `minimum_slp_hpa` from that field and use the minimum SLP grid point for
the reported TC center.

If no accepted sea-level pressure field is present, `minimum_slp_status` is
`not_available`, `minimum_slp_hpa` is null, and the report includes a reason.
The minimum `PSFC` value remains available only as `mslp_proxy_hpa` metadata
with a proxy label; it is not counted as satisfying the minimum SLP objective.

`tools/derive_mslp.py` can add a real diagnostic `SLP` field to a
WRF-compatible core file before running the gate. It copies the source NetCDF
file to a new destination and writes `SLP` in hPa using a hypsometric reduction:
the lowest mass level supplies `P + PB`, `T + 300 K`, and `QVAPOR` for virtual
temperature, the bottom staggered geopotential level from `PH + PHB` supplies
surface height, and `PSFC` supplies surface pressure. This is an approximation
for validation plumbing, but it is not a `PSFC` proxy and it fails clearly if
any required input field is missing or shape-incompatible.

Because the derivation depends on finite `P + PB`, validation tools must not
hide sentinel or pending exposed-cell `P`. A parent-fill candidate with pending
pressure refresh should fail derivation or gate checks until the pressure
producer is fixed. The KROSA pressure-refresh helper and I/O staging exist, and
skeleton metadata reports the refresh as required/not applied, but the compute
path is not invoked from skeleton/remap yet. Reconstructed staging fields,
diagnostic parent-fill, and oracle/reference-copy outputs must not be counted
as gate passes until pressure-refresh compute is truly wired and the 10 min gate
passes.

The pressure-refresh consumer requires `P_TOP`, `C3F`, `C4F`, `C3H`, `C4H`,
and `ALB`. The `ALB` in that source set is WRF Registry `alb`, the inverse base
density, not surface `ALBEDO`; this constraint is retained. It is
input/restart/interp/smooth state without history output, so `wrfout` missing
`ALB` is expected.

The narrow `PB + T_INIT -> ALB` helper is safe only for same-domain, same-time
inputs that already contain `PB` and `T_INIT`. The complete KROSA start-domain
base-state provider targeted in round D27 is:

```text
HGT/P_TOP/C3F/C4F/C3H/C4H/WRF base-atmosphere constants
  -> PB/T_INIT/MUB/ALB/PHB
```

The provider reconstructs `PB`, `T_INIT`, and `MUB`, derives `ALB` from that
same-domain base state, and reconstructs `PHB` on full levels for KROSA
`hypsometric_opt=2`:

```text
PHB(i,1,j) = HGT(i,j) * g
pfu = C3F(k  ) * MUB(i,j) + C4F(k  ) + P_TOP
pfd = C3F(k-1) * MUB(i,j) + C4F(k-1) + P_TOP
phm = C3H(k-1) * MUB(i,j) + C4H(k-1) + P_TOP
PHB(i,k,j) = PHB(i,k-1,j) + ALB(i,k-1,j) * phm * log(pfd / pfu)
```

This full-level `PHB` reconstruction must be paired with mass-level `ALB` and
column `MUB` from the same domain and same valid time. Later restart `ALB` or
`PHB` may be used only as a probe/smoke sample and must not serve as d02
start-time validation truth.

Completing the provider is not evidence that pressure refresh is complete. The
pressure-refresh compute path remains unwired from skeleton/remap, so staging,
diagnostic parent-fill, provider-probe, and oracle/reference-copy outputs must
not be counted as gate passes until pressure-refresh compute is truly wired and
the 10 min gate passes.

Round D28 adds a real-file provider probe and I/O reader boundary for the same
base-state pipeline:

```text
HGT/P_TOP/C3F/C4F/C3H/C4H -> PB/T_INIT/MUB/ALB/PHB
```

The probe/report/staging metadata must remain explicit:

```text
diagnostic = true
gate_candidate = false
integrator_output = false
```

The provider may prove that KROSA real-file constants can drive the base-state
reader and reconstruction, but it is still outside the progressive 10 min gate.
Because pressure-refresh compute is not connected to skeleton/remap yet,
provider success must not relax strict d02 gate criteria or be reported as a
10 min pass. Later restart-file `ALB` or `PHB` may be sampled only for
smoke/range probes and must not become d02 start-time truth.

Round D29 adds the adapter/staging bridge rule for pressure-refresh inputs. The
`base_state_reconstruction_provider` may act as the same-domain/time source for
staged `ALB`, `PB`, `MUB`, and `PHB` buffers consumed by pressure refresh. That
bridge is still not the pressure-refresh compute step. Validation reports must
keep `pressure_refresh_applied` as the deciding marker: if it is `false`, the
candidate cannot pass the d02 gate, even when provider, adapter, and staging
all succeeded.

Round D30 defines the hook helper boundary and must preserve this sequence:

```text
parent-fill/remap -> provider -> base-state sync -> staging -> pressure refresh compute -> report
```

The provider may synchronize `PB`, `MUB`, and `PHB` only to exposed child cells.
The remapped overlap region and halo cells must remain unchanged by this sync.
Provider `ALB` is only an external staging input for pressure refresh and must
not be written into `State`; provider `T_INIT` is the WRF base-state initial
temperature, not perturbation `T`, and must not be written into `State::t`.

Reports should expose parent-fill/remap, provider, base-state sync, staging,
pressure-refresh compute, and final report status separately. A helper-invoked
pressure-refresh compute is still only a diagnostic/skeleton milestone until
the real 10 min d02 gate passes. Before that pass, metadata and reports must
keep `gate_candidate=false` and `integrator_output=false`, even if provider,
staging, and compute all ran.

Later restart-file `ALB` or `PHB` remains usable only for probe/smoke
comparisons and must not be used as start-time truth.

Round D31 defines the opt-in skeleton/remap pressure-refresh diagnostic path.
`tywrf_skeleton_cycle` may call the D30 hook only through an explicit
diagnostic parent-fill mode. That path may produce a report with
`pressure_refresh_applied = true`, but this flag means only that the diagnostic
helper compute was invoked and updated the exposed `P` field in that artifact.
It is not a d02 10 min gate pass and is not validated integrator output.

For this opt-in diagnostic path, CLI/report and NetCDF metadata must remain:

```text
TYWRF_GATE_CANDIDATE = false
TYWRF_INTEGRATOR_OUTPUT = false
TYWRF_VALIDATION_GATE_ONLY = false
```

The default gate must continue to reject or ignore these artifacts unless a
future, explicitly documented diagnostic mode is requested. Even then, its
metrics belong in a diagnostic block and must not satisfy the strict d02 gate.

If direct `ALB` is missing from the source file, but the same-domain
base-state reconstruction inputs are ready, the provider may supply staged
`ALB` for the pressure-refresh helper. Later restart-file `ALB` or `PHB` may
still be sampled only for probe/smoke comparisons; neither field can be used as
start-time truth.

When an opt-in pressure-refresh diagnostic is requested, the selected
`--variables` set and output fields must include the pressure-refresh required
fields. At minimum this covers the exposed `P` target plus the
provider/staging fields required to compute it (`PB`, `MUB`, `PHB`, staged
provider `ALB`) and KROSA constants (`P_TOP`, `C3F`, `C4F`, `C3H`, `C4H`).
Validation tooling should treat missing required fields as a hard diagnostic
failure instead of allowing refresh from uninitialized state.

Round D32 records the expected judgment for the real KROSA opt-in
pressure-refresh diagnostic smoke. A report with
`pressure_refresh_applied = true` provides only diagnostic helper evidence:
the provider, staging, and pressure-refresh compute path ran on real KROSA d02
fields and updated the exposed diagnostic `P` target. It is not a d02 10 min
gate pass, not proof of validated integrator output, and not a substitute for
the progressive `00:10` strict gate.

The strict d02 gate must continue to require real candidate metadata. Files
marked with any of the following metadata are not gate-eligible and must be
rejected or reported as failed by default gate tooling:

```text
TYWRF_GATE_CANDIDATE = false
TYWRF_INTEGRATOR_OUTPUT = false
TYWRF_VALIDATION_GATE_ONLY = false
```

The R32 smoke is therefore limited to two checks: it confirms that the
provider/staging/compute path can execute on real KROSA d02 fields, and it
confirms that the non-gate defense remains active when the artifact is passed
near validation tooling. It must not use later restart truth, must not use
later restart `ALB` or `PHB` as start-time truth, and must preserve d02 at
2 km.

Round D33 tightens the strict-gate metadata boundary. A file is gate-eligible
only when it carries positive candidate metadata from the TyWRF integrator path,
including `TYWRF_GATE_CANDIDATE = true` and
`TYWRF_INTEGRATOR_OUTPUT = true`. Metadata that only says an artifact is near a
gate, validation-only, diagnostic, or eligible for reporting is not enough.
Default strict-gate tooling must reject candidate kinds that identify oracle,
reference-copy, end-state-copy, diagnostic-remap, diagnostic-closure,
diagnostic-pressure-refresh, or other non-integrator artifacts, even when their
numeric fields match the WRF reference.

`report_10min_diagnostics` may report a gate-eligibility block to explain why a
file is or is not suitable for the strict `00:10` gate, but that block is
diagnostic context only. It must not be described as a gate pass and must not
override `cycle_gate.py` strict-gate results. If its input is an oracle,
reference-copy, diagnostic, or metadata-negative candidate, the correct
interpretation remains non-gate diagnostic evidence.

Round D33 also leaves `state_exchange` as planning and accounting metadata
only. A report may describe intended parent/child exchange coverage, selected
fields, exposed-cell counts, or planned interpolation regions, but while
`performed_interpolation = false` it is not evidence that parent-to-child
interpolation ran, not a physical d02 update, and not a validation pass.

Round D34 should start the selected-field parent-to-child interpolation path
from start-state and d01-derived inputs only. The first implementation should
target an explicit selected field set for newly exposed d02 cells and any
required staging fields, then report the fields and cells actually
interpolated. It must not read later d02 restart truth, must not use WRF
cycle-end differences as a delta oracle, and must not copy or blend WRF
end-state fields to improve validation metrics.

Selected-field interpolation output remains a candidate under development until
the strict `2025-07-26_00:10:00` d02 validation gate passes through the normal
gate toolchain. Before that pass, reports may say that interpolation was
performed for selected fields, but they must not claim a gate pass, a physical
10 min forecast acceptance, or full moving-nest compatibility. d02 resolution
must remain `2 km`, and best-track nudging remains out of scope.

Round D35 introduces the intended `selected_field_cycle` boundary for the first
non-oracle selected-field candidate. The path may use only d01/d02
cycle-start inputs, same-time KROSA constants, moving-nest remap,
`StateExchangePlan`, and `parent_child_interpolation`. It must not use WRF
reference-end files, later restart truth, end-state deltas, or any field
derived from a reference-end state. It must write all strict d02 gate fields so
the normal gate can fail on numerical error rather than metadata, missing
variables, or non-finite sentinels.

The expected first outcome is a gate-eligible candidate that fails strict
thresholds numerically, not a claimed pass. Gate tooling may consider
`selected_field_cycle` eligible only when the output carries:

```text
TYWRF_GATE_CANDIDATE = true
TYWRF_INTEGRATOR_OUTPUT = true
TYWRF_VALIDATION_GATE_ONLY = false
TYWRF_CANDIDATE_KIND = selected_field_integrator_v0
```

Those values are allowed only if no end-state/reference-end inputs were used,
d02 `DX/DY` remains `2000 m`, every strict field is finite, and selected fields
changed from the start-state through the remap / exchange / interpolation path.
If any condition is missing, the artifact must be reported as not gate-eligible
before RMSE thresholds are interpreted. `tendency_cycle`, skeleton/persistence,
diagnostic remap, diagnostic pressure-refresh, oracle, reference-copy, and
diagnostic-closure paths remain non-gate even if they produce comparable
fields or diagnostic reports.

Round D35 real KROSA selected-field smoke reached the intended next validation
state: a positive-metadata, non-oracle candidate that fails on real numbers
rather than gate plumbing. `selected_field_cycle` generated the candidate from
`2025-07-26_00:00:00` d01/d02 start states only. Strict metadata passed with
`TYWRF_GATE_CANDIDATE = true`, `TYWRF_INTEGRATOR_OUTPUT = true`,
`TYWRF_VALIDATION_GATE_ONLY = false`, and
`TYWRF_CANDIDATE_KIND = selected_field_integrator_v0`.

The strict field arrays had full finite coverage. The first failing field is
`U`, with normalized RMSE `0.117875` against the `00:10` KROSA d02 reference.
`V`, `MU`, and `P` also fail strict thresholds; `P` normalized RMSE
`0.907405` is the largest blocker. `T`, `PH`, and `QVAPOR` pass their strict
thresholds. This remains a failed gate and must be reported as such.

The remaining diagnostic/output gaps after D35 were separate from the field
RMSE result: accepted `SLP`/MSLP was unavailable, `U10` was absent in the v0
smoke output, and the full MSLP/storm-center/Vmax portions of the strict d02
gate were incomplete. The lack of SLP must not be hidden by a `PSFC` proxy, and
missing or start-state-only 10 m winds must remain hard output/producer gaps
for Vmax10m diagnostics.

Round D36 validation status is a plumbing/output completeness improvement, not
a gate pass. Commit `7d8d37c` makes `selected_field_cycle` default output cover
the strict fields plus available cycle-start-backed `PB`, `PHB`, `MUB`,
`PSFC`, `U10`, `V10`, `T2`, `Q2`, `RAINC`, and `RAINNC`. Metadata and derived
SLP/rainfall side artifacts no longer block report generation by themselves,
but the real KROSA gate still fails numerically and diagnostically.

The D36 strict field result at `2025-07-26_00:10:00` is:

- failed first field: `U`, normalized RMSE `0.117875`;
- additional failed strict fields: `V` `0.134244`, `MU` `0.133382`, and `P`
  `0.907405`;
- passed strict fields: `T`, `PH`, and `QVAPOR`.

The D36 TC diagnostic result is:

- storm center error `150.622 km`, failed against the `20 km` threshold;
- minimum SLP error `0.364 hPa`, passed against the `5 hPa` threshold;
- Vmax error `0.769 m s-1`, passed against the `5 m s-1` threshold.

`U10` and `V10` comparisons still fail as field comparisons because those
variables are preserved cycle-start values, not outputs from a real surface
diagnostics producer. This distinction must remain visible in validation
reports: a passing Vmax diagnostic from preserved fields does not close the
surface-diagnostics producer blocker.

Positive selected-field candidate metadata may be used only for the
non-oracle `selected_field_cycle` candidate file. It cannot be transferred to
diagnostic closure files, diagnostic pressure-refresh files, remap-only
outputs, derived SLP copies, derived rainfall copies, report JSON, oracle
outputs, reference-copy outputs, or any artifact using reference-end truth.
Default gate tooling must continue to evaluate the candidate's own metadata
and candidate kind; helper telemetry cannot override or repair missing,
negative, or incompatible candidate metadata.

Round D37 pressure-refresh validation work is opt-in and must be
cycle-start/provider-backed. If any required pressure-refresh input is missing,
shape-incompatible, non-finite, from `00:10` reference-end truth, or otherwise
not provider-backed for the same domain and cycle-start context, the run must
abort. It must not fall back to a diagnostic closure, a `PSFC`/SLP proxy, WRF
end-state deltas, later restart truth, or metadata-only success.

The first D37 real KROSA opt-in smoke is a negative numerical result. The
candidate generated successfully with `--pressure-refresh`, provider/staging
and pressure compute reported success, and `1,053,150` `P` points changed.
However, the `00:10` strict gate still failed and `P` normalized RMSE worsened
to `5.806413`; D36 without opt-in was `0.907405`. `U`, `V`, and `MU` remained
failed, `T`, `PH`, and `QVAPOR` remained passed, storm-center error remained
`150.622 km`, and MSLP/Vmax still passed. Treat this as safe opt-in plumbing
plus a provider/static consistency diagnostic, not as a validation improvement.

Moving-nest static/coordinate production and surface diagnostics are now
explicit validation blockers. Validation must not use `00:10` reference-end
`XLAT`, `XLONG`, `HGT`, or other static/coordinate fields as oracle truth for
the shifted d02 candidate, and must not use SLP or `U10` shortcuts to bypass
the missing real producers.

Round D38 validation must keep three boundaries separate:

- Static refresh validation may compare candidate `XLAT`, `XLONG`, and `HGT`
  against the `00:10` WRF reference to measure error, but the producer must use
  the cycle-start moving-nest pose shift plus parent/static interpolation or
  reconstruction. G37 measured stale-coordinate symptoms of `XLAT` RMSE
  `0.6913 deg`, `XLONG` RMSE `1.2463 deg`, and same-index displacement
  `152.6 km`; copying reference-end statics would hide that bug and is an
  oracle/reference-end violation.
- Pressure-refresh validation must report D37 as a negative opt-in result:
  provider/staging/compute succeeded and changed `1,053,150` `P` points, but
  `P` normalized RMSE worsened to `5.806413`. The default selected-field gate
  must not enable `--pressure-refresh`, and a refreshed candidate cannot be
  called improved or passing until the strict gate metrics improve through the
  same non-oracle candidate path.
- Surface diagnostics validation must track provenance. Preserved
  cycle-start-backed `U10`, `V10`, `T2`, `Q2`, `RAINC`, and `RAINNC` are useful
  output placeholders but they are not producers. A surface field comparison or
  Vmax result closes only when the fields come from the ABI v2
  SFCLAY/physics-wrapper path or another documented real surface producer.

The first D38 real KROSA static-refresh smoke is a failed gate with improved
coordinates. Static fields compare well against the `00:10` reference:
`XLAT` RMSE `0.000287 deg`, `XLONG` RMSE `0.000195 deg`, and `HGT` RMSE
`0.541 m` with maximum absolute `HGT` error `20.789 m`. The storm-center
diagnostic improves from `150.622 km` to `43.483 km`, but remains above the
`20 km` threshold. Strict field failures remain `U`, `V`, `MU`, and `P`; this
is not a pass and should stop at the `00:10` gate.

Round D39 validation keeps the D38 static-refresh result in the candidate
evidence column but not in the pass column. Commit `517ad28` fixed stale
shifted-d02 static coordinates well enough to reduce the storm-center error,
yet the strict `00:10` gate remains failed with `U` normalized RMSE `0.117875`,
`V` `0.134244`, `MU` `0.133382`, `P` `0.907405`, and storm-center distance
`43.483 km`. `T`, `PH`, and `QVAPOR` passing does not allow validation to
advance to `00:20`; the first failed endpoint is still `00:10`.

`--pressure-refresh` must be treated as guarded experimental plumbing. It is
not a default validation mode and not an improvement claim. D37 worsened `P`
normalized RMSE to `5.806413`, and the P38 same-source wrfout probe showed
about 10% pressure-refresh self-consistency error. Validation tooling and
reports should require a readiness guard that checks required inputs, same-pose
base-state/static consistency, finite staging, and self-consistency before
refresh is allowed to affect a gate-eligible candidate. If the guard fails, the
run should abort or report pressure refresh as skipped; it must not fall back
to reference-end truth, an oracle correction, or a silently refreshed default.

Exposed-cell same-pose analysis for `T`, `PH`, and `P` is diagnostic only
unless it produces values from non-oracle inputs. A valid future producer may
use cycle-start d02 remap state, d01-derived parent fields, same-time KROSA
constants, and documented parent-child interpolation or recompute formulas.
The validation path must not use `00:10` WRF d02 truth as an exposed-cell
source, delta, bias correction, or storm-position guide.

Surface ABI v2 validation must distinguish scaffold readiness from executed
physics. A report may say that ABI v2 staging or SFCLAY call boundaries are
present, but it cannot count `U10`, `V10`, `T2`, `Q2`, `RAINC`, or `RAINNC` as
real surface-producer output until the ABI path actually invokes the documented
surface/physics wrapper and writes those fields from that producer. Preserved
cycle-start values and scaffold-only outputs remain producer gaps for field and
Vmax validation.

Round D40 validation records the D39 guard and scaffold results. Commit
`0ea9a69` makes `--pressure-refresh` fail closed: the real KROSA smoke aborts
on the readiness guard and produces no output. This abort is a successful
safety check, not a pressure-refresh result. Validation reports must keep
`P` as unresolved when the guard blocks refresh, and they must not reinterpret
the guard as a pressure producer, a skipped-pass result, or a reason to advance
past the failed `00:10` endpoint.

D40 parent interpolation should target exposed-cell `T` and `PH` only. A valid
candidate may compare those fields against the `00:10` reference after
generation, but generation must use non-oracle inputs: cycle-start d02 remap
state, d01-derived parent fields, same-time KROSA constants, and documented
parent-child interpolation rules. It must not use `00:10` WRF d02 truth,
reference-end deltas, bias corrections, or storm-center information as an
input. `P` remains outside direct parent interpolation and cannot be promoted
to a production field until a derived/recomputed pressure path passes its
readiness and numerical checks.

`PB`, `PHB`, and `MUB` are read-only validation inputs for this D40 slice. They
may support interpolation context and pressure-refresh staging as
same-domain/same-time base-state fields, but selected-field interpolation must
not overwrite them, and validation must not repair them from later restart or
reference-end truth. If a future provider writes these fields, its report must
separate provider ownership from interpolation ownership.

The ABI v2 sidecar fixture is a validation fixture, not executed physics. The
current `_ex` sidecar path reports `wrapper_unavailable` with
`executed_physics = 0`; therefore fixture success can validate ABI shape,
sidecar metadata, and failure reporting only. It must not satisfy surface field
comparisons or Vmax provenance until the ABI v2 wrapper actually executes the
documented surface/physics producer.

Round D41 validation records `fe01be1` as provider readiness plumbing only.
`KrosaBaseStateProviderTerrainOverride` allows the provider to consume a moved
candidate `HGT` view and report terrain provenance such as
`moved_candidate_HGT`, but provider reconstruction is not pressure-refresh
production and does not make a candidate gate-eligible. It must not be counted
as written `P`, `PB`, `PHB`, `MUB`, `PH`, `T`, or surface diagnostics.

Round D42 should use the moved-terrain provider path only as a read-only
`--pressure-refresh` readiness probe. The probe may verify that moved
candidate `HGT` can reconstruct provider base-state buffers, but it must not
sync those buffers into the candidate, call the pressure hook as a producer, or
write a candidate output when readiness is incomplete. Even with provider
terrain source `moved_candidate_HGT`, `thermodynamic_base_state_consistency_ready`
remains `false`; validation must therefore expect fail-closed behavior and no
output file from `--pressure-refresh`.

The default selected-field validation candidate still preserves `P`, `PB`,
`PHB`, and `MUB` from the start-state source. `P` is not an allowed
parent-interpolation variable and must stay outside direct parent fill until a
derived or recomputed pressure producer passes readiness and numeric checks.
Validation must reject any shortcut using `00:10` reference-end truth, restart
`PHB`/`ALB` as a provider substitute, reference-end deltas, or probe/diagnostic
pressure-refresh artifacts as gate evidence. The strict `00:10` endpoint has
not passed; current failures remain `U`, `V`, `MU`, `P`, and storm-center
distance.

Round D43 validation scope is limited to pressure-refresh hook/API terrain
override support. The hook may accept a moved-terrain override and use it when
reconstructing provider-backed pressure-refresh staging, but that is a
diagnostic hook capability, not the default selected-field producer.

A D43 hook unit test that refreshes exposed `P` with override terrain is only a
hook smoke. It must be reported separately from
`selected_field_cycle --pressure-refresh`, and it must not be used to claim
that the selected-field readiness guard passed or that the real
`2025-07-26_00:10:00` d02 gate passed.

The default selected-field candidate and D42 opt-in guard remain unchanged:
`thermodynamic_base_state_consistency_ready=false`, fail closed, no output file,
and no writes to `P`, `PB`, `PHB`, or `MUB`. `P` is still excluded from parent
interpolation. Validation must continue rejecting `00:10` reference-end truth,
restart `PHB`/`ALB` provider substitutes, and probe or diagnostic
pressure-refresh artifacts as gate evidence. The current strict `00:10` gate is
still failed on `U`, `V`, `MU`, `P`, and storm-center distance.

Round D44 defines the next selected-field pressure-refresh readiness contract.
D41-D43 provide moved `HGT` provider override, a selected-field moved-HGT
provider probe, and hook terrain override, but validation must still treat those
as readiness evidence only. Before `selected_field_cycle --pressure-refresh`
may write any candidate output, its report must explicitly show the provider
terrain source/provenance, base-state reconstruction success, planned sync
counts for base-state and exposed pressure fields, overlap/halo safety for the
moved d02 window, and pressure compute/requested/applied status.

If any readiness item is missing, false, unsafe, or inconsistent,
`--pressure-refresh` must fail closed with no output file. The default path
continues to preserve `P`, `PB`, `PHB`, and `MUB`; validation must not accept
direct parent interpolation of `P`, restart `PHB`/`ALB` as a provider
substitute, `00:10` reference-end truth, or diagnostic artifacts as gate
evidence. A hook diagnostic smoke that refreshes exposed `P` remains useful
hook-level evidence only and must not be counted as selected-field readiness or
a strict `00:10` gate pass. The real `00:10` gate still fails on `U`, `V`,
`MU`, `P`, and storm-center distance.

Round D45 adds hook-level dry-run reporting for the pressure-refresh base-state
sync contract. Validation may inspect the would-sync counts for `PB`, `MUB`,
and `PHB`, overlap/halo write counts, and no-write/no-compute flags, but these
fields are readiness evidence rather than candidate evidence. A valid dry-run
must leave the candidate state unchanged, report that base-state sync was not
applied, and avoid pressure compute.

Round D46 should surface that dry-run contract inside the selected-field
`--pressure-refresh` opt-in guard. The expected validation result is still a
fail-closed abort with no output file while
`thermodynamic_base_state_consistency_ready=false`. Even if the dry-run contract
is ok, validation must reject any output or report that writes `P`, `PB`,
`PHB`, or `MUB`, and must not count the dry-run, hook diagnostic smoke, or any
diagnostic artifact as a strict gate pass.

The strict `2025-07-26_00:10:00` d02 gate remains failed on `U`, `V`, `MU`,
`P`, and storm-center distance. Validation must continue rejecting
reference-end truth, restart `PHB`/`ALB` provider substitutes, direct parent
interpolation of `P`, and diagnostic artifact gate evidence.

Round D47 may add hook-layer scratch-buffer pressure compute dry-run reporting.
Validation may inspect would-refresh `P` point counts, invalid `P` point
counts, and the scratch compute safety status, but only as readiness telemetry.
The scratch path must not write candidate `P`, `PB`, `PHB`, or `MUB`, must not
produce a selected-field output file, and must not set
`pressure_refresh_applied=true`.

Selected-field validation is unchanged in D47: the tool-layer opt-in guard is
not connected to scratch compute and must still fail closed under the D46
`thermodynamic_base_state_consistency_ready=false` rule. A successful hook
scratch dry-run is not a strict d02 gate pass, not selected-field output, and
not evidence that the real `2025-07-26_00:10:00` failures on `U`, `V`, `MU`,
`P`, and storm-center distance are resolved. Validation must continue rejecting
reference-end truth, restart `PHB`/`ALB` provider substitutes, direct parent
interpolation of `P`, and diagnostic artifacts or scratch telemetry as gate
evidence.

Round D48 may surface that hook scratch telemetry in the selected-field
`--pressure-refresh` opt-in guard. Validation reports may inspect dry-run
requested/called/ok fields, would-refresh `P` point counts, and invalid `P`
point counts, but only as readiness evidence. They are not proof that candidate
`P`, `PB`, `PHB`, or `MUB` changed and must not be treated as a strict gate
pass.

The D48 opt-in path must still fail closed while
`thermodynamic_base_state_consistency_ready=false`: nonzero exit, no output
file, no writes to candidate `P`, `PB`, `PHB`, or `MUB`, and no
`pressure_refresh_applied=true`. The default selected-field candidate numerical
path remains unchanged. The real `2025-07-26_00:10:00` gate is still failed on
`U`, `V`, `MU`, `P`, and storm-center distance, and validation must continue
rejecting reference-end truth, restart `PHB`/`ALB` provider substitutes, direct
parent interpolation of `P`, diagnostic artifact gate evidence, and scratch
telemetry as gate evidence.

Round D49 validation records only terrain-source parity for a future
selected-field pressure-refresh apply path. The dry-run guard already uses the
moved candidate `HGT` terrain override; D49 preparation requires any later
reachable apply path to use that same moved-terrain provider input so dry-run
and apply semantics do not diverge before apply is enabled.

This parity is not a validation result. It does not mean
`selected_field_cycle --pressure-refresh` produced an output file or refreshed
candidate pressure. While `thermodynamic_base_state_consistency_ready=false`,
validation must expect fail-closed behavior: nonzero exit, no output, no
candidate writes to `P`, `PB`, `PHB`, or `MUB`, and no
`pressure_refresh_applied=true`. The default selected-field candidate numerical
path remains unchanged.

Even after future apply-path terrain parity is in place, validation must not
count scratch telemetry or an unreachable guarded apply path as a gate pass. The
real `2025-07-26_00:10:00` selected-field gate remains failed on `U`, `V`,
`MU`, `P`, and storm-center distance. Validation must continue rejecting
`00:10` reference-end truth, restart `PHB`/`ALB` provider substitutes, direct
parent interpolation of `P`, diagnostic artifact gate evidence, and scratch or
unreachable-apply telemetry presented as gate evidence.

Round D50 does not change selected-field validation status. The readiness
evidence now available is the moved-`HGT` provider probe, the base-state sync
dry-run contract, the scratch pressure compute dry-run, and future apply
moved-`HGT` terrain parity. These are necessary apply-readiness signals, but
they are not a gate pass and do not prove that candidate `P` can be refreshed
safely.

Validation must keep `thermodynamic_base_state_consistency_ready=false` for
D50. With that guard closed, `selected_field_cycle --pressure-refresh` is still
expected to exit nonzero, produce no output file, leave candidate `P`, `PB`,
`PHB`, and `MUB` unchanged, and report no pressure refresh application. The
default selected-field candidate numerical path remains the validation baseline
for this slice.

Before any later change opens apply, validation requires the following minimum
safety proofs:

- candidate mutation audit with an explicit allowlist of fields that may change
  and a denylist check for fields that must remain untouched;
- moved-`HGT` apply proof that provider probe, dry-run, scratch compute, and
  reachable apply all use the same moved candidate terrain source/provenance;
- overlap/halo no-write proof for remapped overlap cells and all halo cells;
- pre/post apply report consistency for requested, would-refresh, applied,
  invalid, changed, and written point counts;
- a real `2025-07-26_00:10:00` validation run against the WRF 1 h / 10 min
  reference, with RMSE, normalized RMSE, max error, TC diagnostics, and the
  first failed field reported before any later time is attempted.

The 00:10 validation sequence after apply is opened must start from the same
cycle-start WRF state, generate a d02 selected-field `--pressure-refresh`
candidate without reference-end inputs, run the normal comparison tools against
the `00:10` WRF reference, and stop immediately if `U`, `V`, `T`, `PH`, `MU`,
`P`, `QVAPOR`, surface diagnostics, hydrometeors, or TC diagnostics fail their
thresholds. It must not continue to `00:20` until `00:10` passes.

Validation must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
provider substitutes, direct parent interpolation of `P`, diagnostic artifact
gate evidence, scratch telemetry, and unreachable guarded apply plumbing as
gate evidence.

Round D51 adds safety validation for the hook and a controlled selected-field
apply seam, but it does not change selected-field gate status. The default
`--pressure-refresh` path must still fail closed while
`thermodynamic_base_state_consistency_ready=false`: nonzero exit, no output
file, no writes to the candidate, no `pressure_refresh_applied=true`, and no
change to the default selected-field candidate numerical path.

Hook postcondition tests should verify that a real apply, when reached only
through an explicit controlled path, can be accepted only if the changed-field
set is limited to `P`, `PB`, `MUB`, and `PHB`. They must also verify zero
overlap/halo writes and reject invalid compute, non-finite compute,
inconsistent requested/computed/applied/changed counts, or any partial compute
as a failure rather than a successful pressure refresh.

The selected-field controlled seam is allowed only to prove tool-level terrain
provenance: apply must use the moved candidate `HGT`, not metadata terrain.
The seam output, any experimental candidate, and any hook-only report remain
non-gate evidence. They cannot satisfy the strict d02 gate, cannot be reported
as a `00:10` pass, and cannot replace the future normal guarded validation run.

Validation must continue rejecting `00:10` reference-end truth, restart
`PHB`/`ALB` provider substitutes, direct parent interpolation of `P`, diagnostic
artifact gate evidence, scratch telemetry, unreachable guarded apply plumbing,
and controlled-seam outputs presented as gate evidence.

Round D52 validation covers transactional apply semantics and selected-field
experimental report parity without opening the default pressure-refresh guard.
The normal `selected_field_cycle --pressure-refresh` path still runs under
`thermodynamic_base_state_consistency_ready=false` and must fail closed with
nonzero exit, no output file, no candidate writes, no
`pressure_refresh_applied=true`, and no change to the default selected-field
candidate numerical path.

Hook validation for D52 must require transactional behavior: pressure refresh
computes into scratch buffers first, checks the full postcondition, and writes
candidate fields only after that postcondition succeeds. Any failed scratch
compute, invalid value, overlap/halo safety failure, count mismatch, or partial
compute result must leave the candidate completely unchanged.

The selected-field experimental apply seam may be validated only for report
parity. Its report should record pressure-refresh target, refreshed, skipped,
and invalid counts; overlap and halo flags; changed counts for `P`, `PB`,
`MUB`, `PHB`, and unchanged-field audits where available; plus consistency
flags showing whether target/refreshed/changed counts agree. These fields are
not strict-gate metrics and must not be interpreted as a `00:10` pass.

Hidden seam output remains diagnostic-only, not gate-eligible, and not normal
integrator output. Validation must require negative gate/integrator metadata
such as `TYWRF_DIAGNOSTIC_ONLY=true`, `TYWRF_GATE_CANDIDATE=false`, and
`TYWRF_INTEGRATOR_OUTPUT=false`. It must still reject `00:10` reference-end
truth, restart `PHB`/`ALB` provider substitutes, direct parent interpolation of
`P`, diagnostic artifact gate evidence, and hidden seam output presented as
gate evidence.

Round D53 validation may promote the normal selected-field
`--pressure-refresh` path only when readiness is derived from strict
same-artifact evidence rather than from a hard-coded
`thermodynamic_base_state_consistency_ready=false` placeholder. The required
evidence is the normal provider report, base-state sync dry-run, pressure
compute dry-run, and transactional real-apply report for the same candidate.

A normal `--pressure-refresh` output is gate-eligible only if the report shows:

- start-state-based static refresh, with no `00:10` reference-end or oracle
  static fields;
- provider terrain source `moved_candidate_HGT`;
- passing base-state sync dry-run and pressure compute dry-run;
- zero unsafe overlap writes, zero halo writes, zero invalid pressure points,
  zero skipped pressure target points, and consistent requested/refreshed/
  applied/changed counts;
- successful transactional real apply after scratch postcondition checks, with
  candidate writes limited to `P`, `PB`, `MUB`, and `PHB`.

Gate eligibility is not validation success. The hidden
`--experimental-pressure-refresh-apply` seam remains diagnostic-only,
non-gate, and non-integrator output, and diagnostic, probe, dry-run, hidden
seam, report-parity, oracle, and reference-copy artifacts remain invalid as
gate evidence. The first real validation after D53 must be the strict
`2025-07-26_00:10:00` d02 gate against the 1 h / 10 min WRF reference. If that
gate fails, the report must identify the first failed field or TC diagnostic
and the sequence must stop; no `00:20` validation is allowed until `00:10`
passes.

Round D54 records the post-D53 validation state. D53 moved the normal
`selected_field_cycle --pressure-refresh` path from fail-closed to
gate-eligible only when strict same-artifact readiness and transactional apply
evidence pass. That promotion means the candidate may enter the strict gate; it
does not mean the gate passed.

The D53 real KROSA d02 `2025-07-26_00:10:00` gate failed after candidate
metadata passed. The first failed field was `U` with normalized RMSE
`0.117875`. Additional strict field failures were `V` `0.134244`, `MU`
`0.133382`, and `P` `6.33495`. Passing strict fields were `T` `0.013121`,
`PH` `0.017410`, and `QVAPOR` `0.017017`. TC diagnostics still failed storm
center distance at `43.483 km`; minimum SLP error `0.364 hPa` and Vmax error
`0.769 m s-1` passed.

D54 validation work must therefore diagnose the `U` first failure and the
pressure-refresh numerical `P` error before any `00:20` progression. Passing
pressure-refresh apply contracts, terrain parity, mutation audits, or report
count parity is not enough: candidate `P` must show numerical compatibility
against the `00:10` WRF reference through the normal strict gate.

Audit tools and diagnostic reports may explain reference coverage, metadata
eligibility, candidate provenance, and failed metrics, but they are
diagnostic-only. They must not generate candidates, repair candidate metadata,
or replace `cycle_gate.py` evidence. The hidden
`--experimental-pressure-refresh-apply` seam remains diagnostic-only,
non-gate, and non-integrator output.

Round D55 validation direction follows the D54 audit split. The pressure
failure is dominated by the refreshed target region: global `P` normalized RMSE
is `6.334952`, target-region `P` normalized RMSE is `10.33536`, non-target
`P` normalized RMSE is `0.338162`, and the target error fraction is
`0.998219`. `PB`, `MUB`, and `PHB` pass or remain close in the same audit, so
the next validation diagnostics should inspect pressure-refresh
formula/source/vertical/region semantics before treating base-state fields as
failed producers.

The same D54 audit found derived SLP can pass with normalized RMSE `0.003248`.
That does not alter strict-gate status. The `2025-07-26_00:10:00` gate remains
failed on `U`, `V`, `MU`, `P`, and storm-center distance, and validation must
not advance to `00:20`. A passing derived-SLP audit may be reported only as a
diagnostic consistency result; it cannot satisfy the `P` field threshold,
repair storm-center distance, or replace the normal d02 strict gate.

D55 audit reports should add two diagnostic-only views: a pressure-refresh
`P` formula/source/vertical/region audit, and a selected-field `U`/`V`/`MU`
state-error split across remapped overlap, refreshed/exposed target cells, and
non-target cells. These reports may use WRF reference output only as the
comparison target for errors, masks, and attribution. They must not write
candidate fields, copy reference values into a candidate, tune candidate
generation from reference-end truth, or mark hidden-seam output as gate
evidence.

Round D56 starts from integrated D55 commit `cc39cdb`, where full validation of
the codebase passed (`cmake --build`, CTest `29/29`, pytest `141/141`) but the
scientific validation state did not advance. The `2025-07-26_00:10:00` strict
d02 gate remains failed, and no `00:20` validation may be attempted until the
strict `00:10` gate passes.

The D55 pressure formula/source audit is diagnostic-only. It attributes the
strict `P` failure to refreshed target-region perturbation pressure: global
`P` normalized RMSE `6.334952`, target-region `P` normalized RMSE `10.335362`,
non-target `P` normalized RMSE `0.338162`, target error fraction `0.998219`,
and target bias about `-805 Pa`. Full-pressure normalized RMSE is not a
replacement metric for the strict `P` gate: its larger pressure scale can make
normalized error look small while raw RMSE remains nearly unchanged. Candidate
`P` is also not close to start `P`, so persistence cannot explain the audit.

The D55 selected-field state audit is likewise observation-only. The first
failing strict variable remains `U`; `U`, `V`, and `MU` normalized RMSE remain
`0.117875`, `0.134244`, and `0.133382`. Target-region error fractions below
`0.5`, plus candidate fields that are not close to the start state, indicate a
broader selected-field/dynamics evolution gap instead of an isolated
exposed-cell interpolation error.

D56 validation diagnostics should add two non-producing reports:

- a vertical/source-level target-region pressure audit that keeps
  perturbation-pressure `P` as the strict gate field;
- a selected-field evolution audit that compares candidate evolution with WRF
  reference evolution from the same `2025-07-26_00:00:00` start state.

Both reports are observation-only. They may compute diagnostic differences,
regional attribution, and source-level summaries, but they must not generate or
patch candidates, cannot count as gate passes, and cannot describe audit/probe/
hidden-seam output as strict-gate evidence.

D57 starts from D56 commit `985a38f`, where full code validation passed
(`cmake --build`, CTest `29/29`, pytest `150/150`) and the branch was pushed to
`origin/main`. This does not change scientific validation status: the strict
d02 `2025-07-26_00:10:00` gate remains failed, and validation must not advance
to `00:20`.

The D56 pressure source audit may keep a failed top-level status because real
start `wrfout` source entries still do not contain `ALB`. That is a
source-seed audit result for diagnostic readiness, not a strict-gate override
and not an integrator pass/fail substitute.

The pressure vertical diagnostic subreport computed target-region `P` RMSE
`1560.573 Pa`, normalized RMSE `10.335362`, and mean bias `-805.059 Pa`. The
largest level errors are at low-level mass levels `0..4`; level `0` has RMSE
`4170.444 Pa` and bias `-4167.791 Pa`. The source/start evolution fraction is
`6.57739`. Validation should therefore treat current pressure-refresh failure
as a low-level target-region perturbation-pressure formula/source issue until a
non-oracle producer proves otherwise.

The selected-field evolution audit computed evolution normalized RMSE
`U = 0.156800`, `V = 0.160006`, and `MU = 0.197589`. Target-region fractions
were `U = 0.414932`, `V = 0.362148`, and `MU = 0.329146`. Amplitude ratios are
close to `1`, while capture fractions are high but imperfect. This points to a
spatial/phase/timing mismatch in `U`, `V`, and `MU` evolution rather than an
amplitude-only error or a target-region-only failure.

D57 validation diagnostics should add two more observation-only audits:

- a pressure vertical-bias companion-field attribution audit for the low-level
  `P` error;
- a selected-field spatial alignment/shift diagnostic audit for `U`, `V`, and
  `MU` end/evolution error.

Both audits may use WRF reference output for diagnostic comparison and error
attribution only. They must not patch candidates, generate candidate fields,
advance validation to `00:20`, or describe diagnostic/audit/probe/hidden-seam
artifacts as gate passes.

D58 starts from D57 commit `ee2fff4`, where full code validation passed
(`cmake --build`, CTest `29/29`, pytest `160/160`) and the branch was pushed to
`origin/main`. This does not change scientific validation status: the strict
d02 `2025-07-26_00:10:00` gate remains failed, and validation must not advance
to `00:20`.

The D57 pressure vertical-bias audit computed the worst perturbation-pressure
error at low mass levels `0..4`. At level `0`, `P` mean difference is
`-4167.791 Pa` and RMSE is `4170.444 Pa`. Companion fields do not carry a
matching negative bias: `PB` and `MUB` companion differences are about `+1 Pa`,
`MU` is about `+11.625`, full pressure `P + PB` inherits the negative `P`
bias, and `PSFC` is about `-376.164 Pa`. Validation should therefore treat the
current pressure problem as a perturbation-`P` producer/staging issue rather
than a base-state companion mismatch.

The D57 selected-field spatial alignment audit computed only modest gains from
small horizontal shifts. Global normalized RMSE improvements were `U` end
`2.46%` and evolution `2.64%`, `V` end `3.40%` and evolution `3.79%`, and
`MU` end `7.57%` and evolution `7.28%`. This keeps the selected-field failure
outside a simple small-shift explanation.

D58 validation diagnostics should therefore remain non-producing audits:

- a local C++/metadata audit of the perturbation `P` producer and staging path;
- a local C++/metadata and D56/D57 JSON-summary audit of selected-field
  pipeline, timing, and remap behavior.

Both D58 audits are observation-only. They cannot patch or generate
candidates, cannot use reference-end truth to tune output, cannot advance the
progressive sequence to `00:20`, and cannot describe audit/probe/hidden-seam
artifacts as gate passes.

D59 starts from D58 commit `5d9adb5`, where full code validation passed
(`cmake --build`, CTest `29/29`, pytest `167/167`) and the branch was pushed
to `origin/main`. This does not change scientific validation status: the
strict d02 `2025-07-26_00:10:00` gate remains failed, and validation must not
advance to `00:20`.

The D58 pressure producer audit output is
`build/validation/r58_pressure_refresh_producer_audit.json`. It is
diagnostic-only and identifies the local producer inventory around
`compute_krosa_pressure` (`src/dynamics/pressure_refresh.cpp:275`),
`refresh_krosa_moving_nest_pressure`
(`src/dynamics/pressure_refresh.cpp:357`), and
`apply_krosa_moving_nest_pressure_refresh_hook`
(`src/dynamics/pressure_refresh_hook.cpp:309`). The worst level is still mass
level `0`, with maximum absolute `P` mean difference about `4167.790767 Pa`.
Validation should treat this as pressure formula/input and staging evidence,
not as a replacement metric or a gate result.

The D58 selected-field pipeline audit output is
`build/validation/r58_selected_field_pipeline_audit.json`. It is
diagnostic-only and reports seven risk flags: target fraction below half,
target error fraction below half, amplitude near one but RMSE high, modest
best-shift improvement, end/evolution shift mismatch, gate/integrator claims
while audits fail, and large movement-delta schedule/remap risk. The movement
child delta is `{i:60,j:35}`. Target fractions are about `U = 0.407583`,
`V = 0.407583`, and `MU = 0.404762`.

D59 validation diagnostics should remain non-producing audits:

- a pressure formula/input audit for the local producer path, source fields,
  constants, masks, and staging metadata;
- a selected-field schedule/remap timeline audit for movement ordering,
  exposure masks, interpolation, pressure refresh, and output metadata.

Both D59 audits are observation-only. They cannot generate candidates, patch
fields, tune formulas from reference-end truth, promote audit/probe/hidden-seam
outputs to gate evidence, advance to `00:20`, reduce d02 below `2 km`, or add
best-track nudging.

D60 starts from D59 commit `957c3d6`, where full code validation passed
(`cmake --build`, CTest `29/29`, pytest `175/175`) and the branch was pushed
to `origin/main`. This does not change scientific validation status: the
strict d02 `2025-07-26_00:10:00` gate remains failed, and validation must not
advance to `00:20`.

The D59 pressure formula/input audit output is
`build/validation/r59_pressure_formula_inputs_audit.json`. It is
diagnostic-only and reported status `computed_with_flags` with six risk flags.
For `P`, the target-region normalized RMSE is `10.335362` and target mean
difference is `-805.058863 Pa`. Low mass level `0` remains the sharpest
failure, with `P` mean difference `-4167.790767 Pa` and RMSE `4170.444038 Pa`.
`P + PB` at level `0` has mean difference `-4166.786671 Pa`, while companion
target normalized RMSE values are much smaller for `PB` (`8.76e-05`),
`MU + MUB` (`0.000878`), `PH + PHB` (`0.000805`), `T` (`0.012115`), and
`QVAPOR` (`0.016169`). `HGT` target normalized RMSE is `0.495197`, but mean
difference is only `-0.085612 m`. Validation should therefore treat this as
pressure formula-term and input attribution evidence, not as a gate result or
candidate-generation recipe.

The D59 selected-field timeline audit output is
`build/validation/r59_selected_field_timeline_audit.json`. It is
diagnostic-only and reported status `computed_with_flags` with six risk flags.
The movement child delta is `{i:60,j:35}`, large movement is true, `P` is not
listed as parent-interpolated, and the changed/interpolated count mismatch flag
is true. Target fractions remain about `0.407583` for `U` and `V`, and
`0.404762` for `T`, `PH`, `MU`, `P`, and `QVAPOR`.

D60 validation diagnostics should remain telemetry-only:

- a pressure per-column formula-term probe for worst low-level target columns,
  with term-level attribution but no field patching;
- selected-field runtime schedule/timeline metadata emission around move
  checks, post-move parent fill, exposed-cell interpolation, static refresh,
  pressure refresh, and output counts/timestamps.

Both D60 diagnostics are observation-only. They cannot patch fields, tune
formulas from reference-end truth, generate candidate files from oracle data,
promote audit/probe/hidden-seam outputs to gate evidence, advance to `00:20`,
reduce d02 below `2 km`, or add best-track nudging.

D61 starts from D60 commit `7903cbf`, where full code validation passed
(`cmake --build`, CTest `29/29`, pytest `178/178`) and the branch was pushed
to `origin/main`. This does not change scientific validation status: the
strict d02 `2025-07-26_00:10:00` gate remains failed, and validation must not
advance to `00:20`.

The D60 normal `--pressure-refresh` candidate output is
`build/validation/r60_pressure_normal/wrfout_d02_2025-07-26_00:10:00`. Runtime
timeline metadata is present with 11 ordered events:
`cycle_start`, `move_from_to_parent_start`, `overlap_remap`,
`exchange_plan_build`, `parent_interpolation`,
`selected_field_change_summary`, `static_refresh`,
`pressure_refresh_readiness`, `pressure_refresh_apply`, `cycle_end`, and
`output_write_preparation`. The metadata records child movement delta `(60,35)`
and pressure-refresh refreshed, synced, and changed counts. It is telemetry
only: core numerical fields match the D54 normal candidate exactly, with max
absolute difference `0.0` for `U`, `V`, `T`, `PH`, `MU`, `P`, `QVAPOR`, `PB`,
`PHB`, `MUB`, and `HGT`.

The D60 pressure column probe output is
`build/validation/r60_pressure_formula_column_probe_audit.json`. It is
diagnostic-only and selected worst columns `(160,49)`, `(160,50)`, `(160,51)`,
`(161,49)`, and `(161,50)`. Rank 1 level `0` `P` is `-4172.158691 Pa` in the
candidate and `514.484375 Pa` in the reference, for a difference of
`-4686.643066 Pa`; levels `0..4` have mean absolute low-level `P` error
`4391.227588 Pa`. Unique risk codes include
`low_level_column_p_error_large`,
`column_candidate_source_delta_matches_p_error`,
`coefficient_terms_missing_from_candidate_netcdf`, and
`formula_input_json_prior_risks_present`.

The D60 selected-field timeline audit output is
`build/validation/r60_selected_field_timeline_audit.json`. It confirms runtime
timeline attributes are present, but it still reports D59 risks until event
strings are parsed into structured event records. The D60 strict gate report
`build/validation/r60_pressure_normal_gate.json` still fails. The first failing
field is `U` with normalized RMSE `0.117875`; related failures include
`V = 0.134244`, `MU = 0.133382`, `P = 6.334952`, and storm-center error
`43.482716 km`.

D61 validation diagnostics should remain telemetry/audit-only:

- optional same-column runtime observation metadata for pressure attribution,
  without field patches or reference-end formula tuning;
- structural parsing of runtime timeline attributes for selected-field
  movement, remap, static refresh, pressure refresh, and output-preparation
  ordering/counts.

Both D61 tasks are observation-only. They cannot patch fields, tune formulas
from reference-end truth, generate oracle candidates, promote
audit/probe/hidden-seam outputs to gate evidence, advance to `00:20`, reduce
d02 below `2 km`, or add best-track nudging.

D62 starts from local commit `ab68c51`
(`Add runtime pressure probes and timeline audit`). That commit passed local
full validation (`cmake --build`, CTest `29/29`, pytest `183/183`, and
`git diff --check`), but it has not been pushed: both GitHub SSH port `22` and
SSH-over-443 push attempts timed out. `origin/main` therefore remains at
`7903cbf`, and D61 must not be described as remote-synchronized until push is
verified.

The D61 real KROSA smoke output is
`build/validation/r61_pressure_normal/wrfout_d02_2025-07-26_00:10:00`. Runtime
timeline parsing reports `parsed_event_count = 12`; the
`pressure_column_probe` event appears after `pressure_refresh_apply` and before
`cycle_end`, with `expected_order_match = true`, event-count consistency true,
event-name consistency true, and runtime count parity mismatch count `0`. The
probe covers columns `(160,49)`, `(160,50)`, `(160,51)`, `(161,49)`, and
`(161,50)`, levels `0..4`, phases `post_static_refresh` and
`post_pressure_refresh`, and `50` records. The D61 candidate remains
numerically identical to the D60 candidate for `U`, `V`, `T`, `PH`, `MU`, `P`,
`QVAPOR`, `PB`, `PHB`, `MUB`, and `HGT`, with max absolute difference `0.0`.

The D61 strict gate report
`build/validation/r61_pressure_normal_gate.json` still fails at `00:10`. The
first failing field is `U` with normalized RMSE `0.117875`; related failures
are `V = 0.134244`, `MU = 0.133382`, `P = 6.334952`, and storm-center error
`43.482716 km`. This is the controlling validation status. Do not run or report
`00:20` progression until this endpoint passes through the normal strict gate.

D62 validation diagnostics should remain diagnostic-only:

- pressure-core formula observation should expose runtime POD terms such as
  `ALB`, `C3F`, `C4F`, `C3H`, `C4H`, `P_TOP`, `theta`, and pressure-refresh
  intermediates at selected columns/levels;
- a Python pressure-column runtime probe audit should parse the NetCDF probe
  attributes, report before/after pressure deltas, and flag missing formula
  terms or strong negative post-refresh `P`.

Neither D62 diagnostic can generate candidates, patch fields, tune formulas
from reference-end truth, promote probe/audit/hidden-seam output to gate
evidence, advance validation to `00:20`, reduce d02 below `2 km`, or add
best-track nudging.

D62 completion status supersedes the earlier push blocker. Commit `655df5b`
(`Add pressure formula observation audit`) is pushed to `origin/main`, and that
push also delivered D61 commit `ab68c51`
(`Add runtime pressure probes and timeline audit`). The D62 full validation
passed CTest `29/29`, pytest `188/188`, and `git diff --check`.

The D62 runtime probe audit output is
`build/validation/r62_pressure_column_runtime_probe_audit.json`. It parsed
`50/50` records consistently and remains diagnostic-only. The reported risks
are `post_pressure_refresh_p_negative`, `large_p_drop_magnitude`, and
`formula_terms_unavailable`; the minimum post-refresh `P` is
`-4174.87305 Pa`, and the largest drop is `(161,49,k=0)` with
`delta_P = -4257.17773 Pa`.

D63 validation work should wire the D62 core formula observation into the
selected-field `--pressure-column-probe` NetCDF/JSON output and extend the
runtime probe audit to parse those formula observation attributes. These
records may explain the pressure-refresh path at selected columns, but they
remain diagnostic-only. They cannot produce gate candidates, patch fields, tune
formulas from reference-end truth, promote audit/probe/hidden-seam output to
gate evidence, advance validation to `00:20`, reduce d02 below `2 km`, or add
best-track nudging.

D63 is now complete. Commit `3e1f6a7`
(`Add selected-field formula observation telemetry`) is pushed to
`origin/main`, and full validation passed CTest `29/29`, pytest `195/195`, and
`git diff --check`. This is code and telemetry validation only, not a strict
d02 gate pass.

The D63 real candidate
`build/validation/r63_pressure_formula_observation/wrfout_d02_2025-07-26_00:10:00`
and audit
`build/validation/r63_pressure_formula_observation/runtime_probe_audit.json`
confirm that formula observation is present in the selected-field pressure
probe path. The audit reports `25` formula records, all recorded and valid,
and `0` pressure mismatches against the post-refresh probe `P`. This closes the
NetCDF/audit-mismatch suspicion for negative `P`: the negative values are now
evidence of formula/base-pressure term behavior, not serialization or parser
disagreement.

The D63 strict d02 `2025-07-26_00:10:00` gate still fails and remains the
controlling validation status. The first failing field is `U` normalized RMSE
`0.117875`; related failures are `V = 0.134244`, `MU = 0.133382`, `P =
6.334952`, and storm-center error `43.482716 km`. Candidate fields `U`, `V`,
`T`, `PH`, `MU`, `P`, `QVAPOR`, `PB`, `PHB`, `MUB`, and `HGT` are numerically
identical to D61, with max absolute difference `0.0`.

D64 is complete and synchronized. Commit `2db9279`
(`Add pressure budget runtime audit`) is pushed to `origin/main`, and full
validation passed CTest `29/29`, pytest `196/196`, and `git diff --check`.
This is code/audit validation only, not a strict d02 gate pass.

The D64 real audit
`build/validation/r64_pressure_budget_runtime_probe_audit.json` on the D63
candidate has `status = computed_with_flags`. It reports `25` pressure budget
records, `25` `total_pressure < PB` records, `25` records where the large
pressure drop is explained by formula/base-pressure subtraction, and pressure
mismatch count `0`. The first tracked point `(i=160,j=49,k=0)` has
`total_pressure = 95539.724 Pa`, `PB = 99711.8828 Pa`,
`total_pressure_minus_pb = -4172.1588 Pa`, and
`probe_delta_p = -4252.3540025 Pa`.

The audit conclusion is that TyWRF and WRF runtime both use
`P = total_pressure - PB`, so validation work must not remove `PB` subtraction
as a shortcut. The strict d02 `2025-07-26_00:10:00` gate still fails. D65
validation diagnostics should examine formula sensitivity and WRF
moving-nest/base-state call order, with emphasis on `PH + PHB`, `theta`,
`MU + MUB`, and base-state staging.

D65 remains diagnostic-only. It cannot generate candidates, patch formulas or
fields, tune formulas from reference-end truth, promote audit/probe/
hidden-seam output to gate evidence, advance validation to `00:20`, reduce d02
below `2 km`, or add best-track nudging.

D65 is complete, synchronized, and still not a gate pass. Commit `07908b1`
(`Add pressure formula sensitivity audit`) is pushed to `origin/main`; the D65
full validation passed CTest `29/29`, pytest `197/197`, and
`git diff --check`.

The D65 real audit
`build/validation/r65_pressure_formula_sensitivity_audit.json` is
diagnostic-only. It reports `25` sensitivity records, a `PB - total_pressure`
gap range of `3593.9011` to `4174.8731 Pa` with mean `3904.1623 Pa`, and an
approximate fractional total-pressure increase needed of `0.0389068` to
`0.0436995` with mean `0.0414497`. All `25` records require a large
total-pressure increase. These numbers are pressure-attribution evidence, not
permission to tune the pressure formula, patch fields, or accept an
audit-derived candidate.

B65's WRF call-path audit found that the broad WRF base-state interpolation
coverage includes `PHB`, `MUB`, `PB`, `ALB`, `T_INIT`, and `HT`; the TyWRF
selected-field path currently interpolates only `U`, `V`, `T`, `PH`, `MU`, and
`QVAPOR`. D66 validation direction is to expose a d02 base-field provenance
contract and audit, plus WRF generated interpolation mask semantics, so reports
can state whether newly exposed d02 base fields came from overlap remap,
parent/base-state interpolation, generated-field staging, or remained missing.
This must not open the selected-field numerical path yet.

The strict d02 `2025-07-26_00:10:00` gate remains failed. Validation must not
advance to `00:20`, must not generate reference-end/oracle candidates, must not
promote audit/probe/hidden-seam output to gate evidence, must not reduce d02
below `2 km`, and must not add best-track nudging.

At D68 start, D66 commit `09d6ba2`
(`Add moving-nest base-state exchange contract`) and D67 commit `5a78fbb`
(`Add base-state exchange action diagnostics`) have both been pushed.
`main` and `origin/main` are synchronized at `5a78fbb`. D66 validation before
the local commit passed CTest `29/29` and pytest `197/197`; D67 does not
change the strict gate result.

The D66 artifact is a diagnostic contract, not gate evidence. It records the
current selected-field exchange set as `U`, `V`, `MU`, `QVAPOR`, `T`, and
`PH`; records the WRF broad base-state candidate set as `PHB`, `MUB`, `PB`,
`ALB`, `T_INIT`, `HT`, and `HGT`; excludes `P` from that base-state candidate
set; and marks the contract with `diagnostic_only = true` and
`enables_selected_field_numerics = false`. Its exposed-child exchange plan
counts exposed strips and point totals for future interpolation, but reports
`performed_interpolation = false`, `modifies_overlap = false`, and
`modifies_halo = false`.

B66 found that WRF generated unstaggered-mask behavior, represented by
`imask_nostag`, is part of the expected semantics for `PHB`, `T_INIT`, `MUB`,
`ALB`, `PB`, and `HT`. B66 also found that WRF `start_domain` recompute rules
must be represented separately from selected-field interpolation. Validation
reports must therefore separate exposed-mask provenance/action facts from
validated numerical production.

D67 completed an opt-in diagnostic moving-nest base-field provenance/action
report and prepared the exposed-mask regression direction. These diagnostics
describe whether exposed d02 base fields came from overlap remap, parent/base
interpolation, static/provider generation, recompute staging, or remained
missing. The diagnostic action vocabulary includes
`interpolate_exposed_cells`, `recompute_from_mub_after_interpolation`,
`preserve_interpolated_when_rebalance_zero`, and `static_height_input`, but
these labels are report facts only.

D68 validation direction is a diagnostic-only exposed base-state exchange
helper and regression. The helper/test may expose `PHB`/`MUB`/`HT` exposed
interpolation semantics and `PB`/`T_INIT`/`ALB` recompute marks for newly
exposed d02 cells, but those facts are helper/report semantics only. They must
not be connected to the production selected-field numerical path, must not be
treated as evidence for a gate pass, and must not produce reference-end truth
or oracle candidates. The strict d02 `2025-07-26_00:10:00` gate remains failed;
validation must not advance to `00:20`, must not reduce d02 below `2 km`, and
must not add best-track nudging. D70 or later may consider hook or diagnostic
wiring for WRF-style exposed base-state policy; D68 keeps the migration as
documented future work only.

D68 is now complete and synchronized. Commit `c8a83a2`
(`Add exposed base-state exchange diagnostics`) has completed full validation
and has been pushed to `origin/main`, leaving `main` and `origin/main`
synchronized. This is still diagnostic/helper validation only. It does not
connect to the production selected-field path and does not change the strict
d02 `2025-07-26_00:10:00` gate result, which remains failed.

D69 is now complete, fully validated, pushed, and synchronized at commit
`e32ccc9` (`Add exposed MUB base-state recompute provider`). Its validation
coverage is provider/test-only: it verifies recomputing `PB`, `T_INIT`, and
staged provider `ALB` from already exposed-interpolated `MUB`, but it still
does not route that recompute result into the production selected-field path,
write a gate candidate, or report a gate pass. The recompute API must not
regenerate `MUB` from `HGT`, rebuild or sync `PHB`, write `T_INIT` into
`State::t`, or write `ALB` into `State`. The strict d02
`2025-07-26_00:10:00` gate remains failed.

D70 is complete, fully validated, pushed, and synchronized at commit
`5eb6485` (`Add exposed base-state adapter diagnostics`). It added only a
diagnostic adapter/report for the D68 exposed base-state exchange helper and
the D69 exposed-`MUB` recompute API. It did not produce a production
candidate, did not connect the normal `selected_field_cycle`
pressure-refresh gate path, and did not alter production selected-field
numerics. Any D70 report remains non-gate with `diagnostic_only = true`,
`gate_candidate = false`, and `integrator_output = false`. The strict d02
`2025-07-26_00:10:00` gate remains failed.

D71 validation scope is a metadata guard before any adapter connection. D71
must not connect the D70 adapter into the production `selected_field_cycle`
normal path. The guard should make the C++ NetCDF metadata and stdout JSON use
the same candidate-disposition decision so the file attributes and report
fields cannot disagree about gate eligibility.

The default Python strict gate must reject pseudo-positive metadata from
helper, probe, adapter, dry-run, staging, experimental, diagnostic, oracle, and
reference-copy artifacts. That rejection must win even if a malformed artifact
claims positive `gate_candidate` or `integrator_output` fields. Diagnostic,
oracle, probe, helper, adapter, dry-run, staging, and experimental outputs
remain excluded from strict-gate evidence. Validation must not use
reference-end truth or oracle candidates, must not advance to `00:20`, must not
lower d02 below `2 km`, and must not introduce best-track nudging.

D71 is now complete, fully validated, pushed, and synchronized at commit
`9b08d09` (`Guard diagnostic adapter gate metadata`). Its result is a metadata
guard, not a numerical model improvement: `main` and `origin/main` are
synchronized, and the strict d02 `2025-07-26_00:10:00` gate still fails. The
guard requires strict-gate tooling to reject helper/probe/adapter/dry-run/
staging/experimental/diagnostic/oracle/reference-copy roles before reporting
any pass/fail interpretation from field metrics.

D72 may add an opt-in diagnostic adapter report path only. It must remain
separate from the normal pressure-refresh candidate path and cannot be used as
the pressure-refresh producer for a strict d02 candidate. Every D72 diagnostic
adapter report or artifact must carry non-gate metadata:

```text
diagnostic_only = true
gate_candidate = false
integrator_output = false
```

Default strict gate behavior must reject D72 diagnostic adapter output. A D72
diagnostic report can explain adapter wiring state, but it is not a model
pass, cannot satisfy `00:10`, and cannot permit validation to advance to
`00:20`. It must not use reference-end truth, oracle/reference-copy output,
diagnostic/probe/helper artifacts, adapter telemetry, or restart substitutes as
gate evidence; d02 remains `2 km`, and best-track nudging remains out of scope.

D72 is now complete, fully validated, pushed, and synchronized at commit
`616c6c9` (`Add diagnostic adapter selected-field report path`). The new
`--diagnostic-base-state-adapter-report` path is opt-in diagnostic reporting
only. It may expose adapter staging/audit state for the selected-field
exposed-base-state work, but it is not a pressure-refresh producer, not normal
integrator output, and not a strict-gate candidate. Its report and any artifact
metadata must remain:

```text
diagnostic_only = true
gate_candidate = false
integrator_output = false
```

D73 is now complete, verified, pushed, and synchronized at commit `b73704f`.
It verified that the opt-in adapter smoke/audit path reports expected
adapter/helper/staging metadata and that strict-gate rejection remains active.
The rejection is metadata-first: default strict-gate tooling must reject
adapter, helper, staging, dry-run, probe, diagnostic, oracle, and reference-copy
roles/kinds before interpreting RMSE fields or TC diagnostics. Passing-looking
field RMSE, storm-center, MSLP, or Vmax numbers from one of those artifacts must
still be reported as non-gate diagnostic evidence, not as a
`2025-07-26_00:10:00` validation pass. The strict `00:10` gate still fails, so
validation must not proceed to `00:20`.

D74 is complete, verified, pushed, and synchronized at commit `b1d2de2`
(`Add base-state source staging provider`). It staged only the non-oracle
child-shaped `BaseStateSourceStagingProvider` for diagnostic/test use. Its
purpose is to replace `source==child` staging in a later hidden adapter wiring
round, not to change validation behavior in that provider round. D74 did not
connect the provider to `selected_field_cycle`, write or relabel normal
candidates, set positive gate metadata, or weaken the strict
metadata/disposition checks. Default strict-gate tooling must continue to treat
D74 provider reports, fixtures, and staged outputs as non-gate artifacts:

```text
diagnostic_only = true
gate_candidate = false
integrator_output = false
```

D75 is complete, verified, pushed, and synchronized at commit `d5a1f99`
(`Wire diagnostic adapter source staging provider`). Its validation scope was
hidden diagnostic adapter source staging only. The D74 provider may back hidden
adapter scratch-source staging, and reports may expose provider-backed source
provenance plus staged/exposed/masked count metadata. That metadata remains
disqualifying and cannot be used as gate evidence. D75 did not write or relabel
normal candidates, refresh candidate `P`, change the normal pressure-refresh
path, mark a report as integrator output, or relax strict metadata rejection.

D76 is complete, verified, and pushed at commit `f7082c1`. Its validation scope
was diagnostic source-vs-child delta/provenance planning for the hidden
base-state adapter source-staging path. Default strict-gate tooling must treat
any such source-child delta/provenance block as disqualifying diagnostic context
before field RMSE or TC diagnostics are interpreted. This is true even if the
reported source-child differences are all zero, or if the staged/exposed/masked
counts are useful for debugging. Source-child provenance does not prove
integrator output, does not satisfy the pressure-refresh contract, and cannot be
used as evidence that the strict d02 gate passed. D76 did not write normal
candidates, refresh candidate `P`, change the normal pressure-refresh path,
relax strict metadata rejection, or proceed to `00:20`. The strict d02
`2025-07-26_00:10:00` gate still fails; d02 stays `2 km`, and best-track
nudging remains prohibited.

D77 is complete, verified, pushed, and synchronized at commit `c7e46a7`. Its
validation scope remained hidden diagnostic provider-derived base-state source
staging only. The hidden adapter source now uses provider-derived `PB`,
`T_INIT`, `MUB`, `ALB`, and `PHB` plus `output_static` `HGT`/`HT`, with
provider-source metadata and strict-gate regression coverage. This does not
change normal selected-field candidate eligibility or strict-gate eligibility.

Default strict-gate tooling must keep provider-derived source metadata,
source-staging metadata, and source-child delta summaries as diagnostic-only
disqualifying context before interpreting RMSE or TC diagnostics. These
attributes cannot prove integrator output, cannot satisfy the pressure-refresh
contract, cannot refresh candidate `P`, cannot prove that the strict
`2025-07-26_00:10:00` gate passed, and cannot permit validation to advance to
`00:20`.

Any `PHB` reconstruction observed in D77 provider-derived source staging is
diagnostic evidence only. It is not WRF rebalance semantics, not a validated
candidate `PHB` producer, and not a substitute for a normal pressure-refresh or
base-state producer.

D78 is complete, verified, pushed, and synchronized at commit `d7ecd11`. It
added the diagnostic adapter provider-source audit CLI for provider-source,
source-staging, and source-child-delta attributes. The audit is diagnostic-only
and fail-closed: it may summarize those attributes and identify the largest
deltas for debugging only. It must not write candidate fields, refresh candidate `P`,
change pressure-refresh behavior, relax strict metadata rejection, produce gate
candidates, proceed to `00:20`, reduce d02 below `2 km`, or introduce
best-track nudging. D78 audit output must not be used as strict gate evidence
or as proof that the failed `00:10` endpoint passed.

D79 may add a normal selected-field exposed base-state producer for legitimate
candidate base-state fields. Validation may accept only sources from the
current cycle-start state, d01-derived parent data, and same-time KROSA
metadata/constants. The D79 path must not use D77 hidden diagnostic adapter
NetCDF attributes, reference-end truth, WRF end-state deltas, a direct `P`
shortcut, gate metadata relaxation, a `PSFC`-as-`SLP` proxy, or diagnostic
artifact promotion.

D79 must not directly patch perturbation `P`; pressure refresh remains
separately guarded and must keep its own readiness/reporting checks. D79 also
must not claim a strict `2025-07-26_00:10:00` pass or allow validation to
advance to `00:20`. The current first failure remains `U` normalized RMSE
`0.117875` until real strict-gate output proves otherwise. Any `PHB`
reconstruction remains limited by documented provider semantics and must not be
reported as WRF rebalance equivalence.

D79 is complete, passed focused and full validation, and was pushed at commit
`3d49d9c` (`Add normal selected-field base-state producer`). The D79 validation
milestone verifies the normal selected-field base-state producer wiring, but it
does not change the progressive gate state. The real KROSA d02
`2025-07-26_00:10:00` gate still fails; the first failed field is `U` with
normalized RMSE `0.11787539215928292`. Validation must stop at `00:10` and
must not run or report `00:20` as an advanced gate.

D80 may update only normal pressure-refresh metadata eligibility names. The
allowed migration is from helper/dry_run/staging blocker wording to gate-safe
production/readiness/source-sync wording for the normal production path. This
does not relax strict metadata rejection. Default strict-gate tooling must
still reject diagnostic, oracle, helper, probe, adapter, staging, dry_run, and
experimental artifacts before interpreting RMSE or TC diagnostics.

The D80 metadata refresh is not a field validation pass. Reports must continue
to expose the unresolved `U`, `V`, `MU`, and `P` RMSE failures and the failed
TC center diagnostic. Metadata eligibility fixes cannot hide those failures,
cannot convert helper or diagnostic artifacts into candidates, and cannot allow
progression beyond the failed `2025-07-26_00:10:00` endpoint.

D80 is complete at commit `765bd06`
(`Make pressure-refresh production metadata gate-safe`). It passed focused
validation, full validation, the real `00:10` metadata check, and was pushed.
The current real `00:10` gate state is therefore split deliberately:
`candidate_metadata` passes, but the overall gate still fails numerically. The
first failed field remains `U`, with normalized RMSE
`0.11787539215928292`; validation must stop at
`2025-07-26_00:10:00` and must not proceed to `00:20`.

D81 may run only a U/V wind-failure audit. The audit may read reference and
candidate `U`/`V` at the failed `00:10` endpoint and break down the error by
diagnostic categories, but every report or artifact from that audit must be
diagnostic-only and carry `gate_evidence=false`. D81 must not write a
candidate, modify thresholds, change gate metadata, count as `00:10` pass
proof, or advance the validation sequence. It also must not handle `P`, `MU`,
physics, or best-track nudging, and d02 must remain at `2 km`.

D81 is complete at commit `3b17b6b`
(`Add selected-field wind error audit`). It passed focused validation, full
validation, the real D80 wind-audit smoke, and was pushed. The audit result is
that D80's real `00:10` `U`/`V` error is domain-wide. The metadata portion
passes, but the field gate does not; the first failed field remains `U` with
normalized RMSE `0.11787539215928292`. The progressive sequence remains stopped
at `2025-07-26_00:10:00`; no `00:20` run or report may be treated as the next
gate until `00:10` passes.

D82 validation scope is limited to parent-child interpolation `U`/`V`
stagger-coordinate semantics and directly related tests or metadata. It may
clarify or test that `U` interpolation uses west-east staggered coordinates and
`V` interpolation uses south-north staggered coordinates. It must not work on
`P`, `MU`, physics, best-track nudging, threshold changes, oracle candidates,
or any d02 resolution change below `2 km`.

Any D82 coordinate metadata, coordinate audit, or wind-audit report is not gate
evidence. Validation status can change only after a real strict d02
`2025-07-26_00:10:00` gate run reports positive integrator metadata and field
metrics that prove improvement. Until that output exists, the field gate is
failed and validation must not advance to `00:20`.

D82 is complete at commit `0a192d4`
(`Use WRF-style staggered interpolation coordinates`). It passed focused
validation, full validation, the real `00:10` gate smoke, and was pushed. The
real `2025-07-26_00:10:00` d02 gate smoke confirms candidate metadata passes,
and `U`/`V` improved finitely after WRF-style staggered interpolation
coordinates, but the strict gate still fails. The first failed field is `U`
with normalized RMSE `0.11359509344276145`; `V` normalized RMSE is
`0.12937874064226143`; `P` normalized RMSE is slightly worse at
`6.354144377000247`. That `P` change is a regression risk to keep visible, not
a reason to adjust thresholds or broaden the next worker scope.

D83 may only audit `U`/`V` source/time-level attribution for the failed
`2025-07-26_00:10:00` endpoint. Its report must be diagnostic-only and
`gate_evidence=false`; it cannot write a candidate, update thresholds, alter
gate metadata, or be cited as proof that `00:10` passed. D83 must not work on
`P`, `MU`, physics, best-track nudging, or d02 resolution, and it must not run
or report `00:20` as the next validation gate.

Example:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/derive_mslp.py \
  /path/to/wrfout_d02_2025-07-26_06:00:00 \
  /path/to/with-slp/wrfout_d02_2025-07-26_06:00:00 \
  --pretty
```

`tools/run_cycle_gate_with_slp.py` wraps the common validation plumbing: for
each cycle-end file, it derives `SLP` for both reference and candidate into
separate output directories, then runs `tools/cycle_gate.py` against those
derived copies. The combined JSON report records the derived directories,
per-file derivation summaries, and the nested gate report. A derivation failure
is a hard overall failure; the wrapper does not fall back to a `PSFC` proxy.

Example:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/run_cycle_gate_with_slp.py \
  --reference-dir /path/to/reference/WRF \
  --candidate-dir /path/to/tywrf/output \
  --start 2025-07-26_00:00:00 \
  --end 2025-07-26_06:00:00 \
  --derived-dir build/validation/cycle_gate_with_slp \
  --pretty
```

The checked KROSA d02 reference files under
`/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF`
do not contain the accepted SLP candidates. Across 29 `wrfout_d02_*` files, the
pressure-related matches were `P`, `PB`, `PSFC`, and non-equivalent `SEAICE`.

## Opt-In Diagnostic Closures

Diagnostic closures are separate from normal d02 cycle-gate acceptance. A file
modified by a closure is a diagnostic candidate, not physical TyWRF-Core
integrator output, even when it keeps WRF-compatible variable names.

`pressure_mu_closure` is the current proposed closure. Its design contract is
tracked in `docs/diagnostic_closures.md`. It may be used only to study local
pressure and dry column mass consistency from cycle-start fields. It must not
use a WRF reference end-state delta, copy cycle-end reference pressure fields,
patch validation metrics, or hide modified pressure fields inside the normal
skeleton, integrator, or writer path.

The minimum allowed closure inputs are cycle-start `MU`, `P`, `PB`, `PH`,
`PHB`, `T`, `QVAPOR`, and `PSFC`. Future implementations may use a hydrostatic
column consistency approximation or a local mass-pressure closure, but the
result remains a `not_physical` / `diagnostic_closure` artifact.

Closure output must carry machine-readable metadata, including:

```text
diagnostic_closure = "pressure_mu_closure"
diagnostic_candidate = "true"
not_physical = "true"
excluded_from_default_gate = "true"
uses_reference_end_delta = "false"
```

The only opt-in closure metrics currently in scope are `MU` normalized RMSE,
`P` normalized RMSE, derived SLP minimum error, and `PSFC`/`SLP` range sanity
checks. The closure does not claim or improve `U`, `V`, `QVAPOR`, storm-center
distance, Vmax, or rainfall metrics. If a report locates the minimum of a
closure-derived `SLP` field, it must label that point as a closure SLP minimum,
not as the accepted TC storm center.

Default validation tools must either exclude closure-marked files from normal
pass/fail gates or require an explicit diagnostic-closure mode that reports
closure metrics in a separate block. Closure-derived `SLP` does not satisfy the
default minimum-SLP or storm-center gate.

## Python Validation Tools

`tools/audit_reference_cycles.py` checks KROSA reference coverage before cycle
comparisons. Its default audit targets the 1 h / 10 min validation reference
directory, both `d01` and `d02`, and the required validation times from
`2025-07-26_00:00:00` through `2025-07-26_01:00:00`. It verifies the v1 core
output variables, including `Times`, checks `DX`/`DY` as grid-spacing metadata
or variables, reports real SLP candidates from the accepted SLP/MSLP list, and
records `PSFC` separately as a pressure proxy. The report includes
`missing_files`, `missing_variables`, `missing_grid_spacing`,
`available_grid_spacing`, `available_slp_candidates`,
`available_pressure_proxy_candidates`, and `cycle_count`; use `--format table`
for a compact terminal summary. `--interval-minutes` can audit alternate minute
cadences, while `--interval-hours` remains available for older 6 h coverage
checks. This tool is diagnostic-only: it audits reference/candidate coverage and
does not produce gate candidates or supply gate-pass evidence.

`tools/baseline_candidate.py` generates explicitly marked baseline candidate
`wrfout` files without running a TyWRF-Core integrator. `persistence` copies
the cycle-start reference file to the cycle-end candidate path and rewrites
`Times` to the cycle-end valid time; this is expected to fail strict thresholds
and quantifies the no-integration baseline gap. `reference-copy` copies the
cycle-end reference file and marks both JSON metadata and NetCDF attributes with
`reference_copy=true`, `integrator_output=false`, and
`validation_gate_only=true`; it must be used only to verify the gate plumbing.
WRF end-state/oracle output is never accepted as evidence of a TyWRF gate pass.
Its `expected_to_meet_thresholds` flag is true only when a real SLP/MSLP
diagnostic was copied, because a reference copy without SLP still fails the
full d02 cycle gate.
For d02, the tool rejects sources whose `DX`/`DY` attributes are not 2 km.

`tools/extract_reference.py` inspects WRF NetCDF files and reports dimensions,
present core variables, missing core variables, and metadata for each present
core variable. It exits nonzero when required core variables are missing unless
the expected variable list is narrowed with `--core-variables`.

`tools/compare_wrfout.py` compares selected numeric fields between a WRF
reference file and a TyWRF-Core candidate file. The JSON report contains:

- per-variable status: `ok`, `threshold_exceeded`, `missing_reference`,
  `missing_candidate`, `shape_mismatch`, `non_numeric`, or `no_valid_samples`;
- RMSE, normalized RMSE, max absolute error, threshold, sample counts, and
  reference/candidate shapes when metrics are available;
- summary counts and top-level pass/fail status.

Default normalized RMSE thresholds follow the field thresholds above. Use
`--no-thresholds` to report metrics without threshold failures, or repeat
`--threshold VARIABLE=VALUE` to override individual thresholds.

Use `--tc-diagnostics` to add an opt-in TC diagnostics block while preserving
the default RMSE-only comparison behavior. The TC block reports minimum SLP when
available, a separately labeled `PSFC` proxy when minimum SLP is not available,
10 m wind maximum from `sqrt(U10^2 + V10^2)`, accumulated rainfall summaries
from `RAINC + RAINNC`, and reference-candidate errors for center distance,
minimum SLP, proxy pressure metadata, Vmax, and rainfall summary metrics. Use
`--diagnostic-time-index` when a file has multiple time records; the default is
the last record.

`tools/run_6h_cycle_test.py --dry-run` resolves the expected
`wrfout_d01_YYYY-MM-DD_HH:MM:SS` and `wrfout_d02_YYYY-MM-DD_HH:MM:SS` files for
the cycle start and end times without invoking the integrator. The dry-run plan
reports expected WRF reference and TyWRF-Core candidate files for both cycle
endpoints and both domains by default.

Example:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/run_6h_cycle_test.py \
  --reference-dir /path/to/reference/WRF \
  --candidate-dir /path/to/tywrf/output \
  --start 2025-07-26_00:00:00 \
  --dry-run \
  --pretty
```

To generate a d02 persistence baseline candidate for one cycle:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/run_6h_cycle_test.py \
  --reference-dir /path/to/reference/WRF \
  --candidate-dir /path/to/baseline-candidates \
  --start 2025-07-26_00:00:00 \
  --domain d02 \
  --mode persistence \
  --pretty
```

Use `--mode reference-copy` only to validate the gate itself. With `--end`, the
same command generates consecutive d02 cycle-end candidates at `--interval`
hour spacing while preserving the d02 2 km resolution check.

TC diagnostics remain explicitly pending in JSON reports unless
`--tc-diagnostics` is requested.
