# Diagnostic Closures

Diagnostic closures are opt-in validation artifacts. They can be used to study
whether a reduced candidate file has internally consistent diagnostic pressure
fields, but they are not TyWRF-Core physics or dynamics results.

No diagnostic closure may be hidden inside the normal skeleton, dynamics loop,
physics bridge, nesting path, or WRF-compatible writer default path. A closure
output must be visibly separate from normal integrator output and must carry
metadata that prevents it from being mistaken for a physical forecast.

## pressure_mu_closure

`pressure_mu_closure` is a proposed diagnostic candidate for pressure and dry
column mass consistency studies. It is not a forecast correction and is not an
accepted substitute for missing dynamics, physics, lateral-boundary evolution,
spectral nudging, moving-nest feedback, or surface diagnostics.

The closure may produce candidate `MU`, `P`, `PSFC`, or derived `SLP` fields for
diagnostic reports only. These fields must be treated as `not_physical` and
`diagnostic_closure` data even when their names match WRF core variable names.

## Moving-Nest Pressure Refresh Boundary

Moving-nest parent-fill staging is not a pressure closure. Its direct
parent-fill fields intentionally exclude exposed-cell `P`; that field must be
derived or recomputed by a real pressure producer before a physical gate can use
it. Diagnostic parent-fill files and 10-minute reports may surface the pending
refresh metadata, but they must remain non-physical and non-gate artifacts.

The KROSA pressure-refresh helper and I/O staging are present, and skeleton
metadata reports pressure refresh as required and not applied. The compute path
is not invoked from skeleton/remap yet. Therefore parent-fill outputs that lack
applied pressure refresh remain diagnostic only; they must not hide sentinel
`P`, rename `PSFC` as pressure, patch validation metrics, or satisfy a normal
gate.

The pressure-refresh consumer needs `P_TOP`, `C3F`, `C4F`, `C3H`, `C4H`, and
`ALB`. `ALB` here is WRF Registry `alb`, the inverse base density on `ikj`,
not the surface albedo variable `ALBEDO`; this constraint is retained. Registry
flags give it input/restart/interp/smooth behavior (`irdus`) but no history
output, so missing `ALB` in `wrfout` is expected.

The narrow `PB + T_INIT -> ALB` helper is valid only when same-domain,
same-time `PB` and `T_INIT` are already available. It is a helper, not the
complete KROSA d02 start-time provider. Round D27 targets the full base-state
provider:

```text
HGT/P_TOP/C3F/C4F/C3H/C4H/WRF base-atmosphere constants
  -> PB/T_INIT/MUB/ALB/PHB
```

The provider reconstructs `PB`, `T_INIT`, and `MUB`, derives `ALB` from that
same-domain base state, and reconstructs `PHB` on full levels for KROSA
`hypsometric_opt=2` with the WRF log-linear relation:

```text
PHB(i,1,j) = HGT(i,j) * g
pfu = C3F(k  ) * MUB(i,j) + C4F(k  ) + P_TOP
pfd = C3F(k-1) * MUB(i,j) + C4F(k-1) + P_TOP
phm = C3H(k-1) * MUB(i,j) + C4H(k-1) + P_TOP
PHB(i,k,j) = PHB(i,k-1,j) + ALB(i,k-1,j) * phm * log(pfd / pfu)
```

This `PHB` reconstruction must be paired with mass-level `ALB` and column
`MUB` from the same domain and same valid time. Later restart-file `ALB` or
`PHB` can only support probe/smoke checks and must not be promoted to d02
start-time validation truth.

Completing the provider does not complete pressure refresh. The
pressure-refresh compute path remains unwired from skeleton/remap, so
parent-fill, staging, provider-probe, and closure outputs are still
non-physical and not gate-eligible.

The round D28 real-file provider probe and I/O reader are diagnostic boundary
checks for the base-state pipeline:

```text
HGT/P_TOP/C3F/C4F/C3H/C4H -> PB/T_INIT/MUB/ALB/PHB
```

Their purpose is to verify that KROSA real-file constants can feed the
provider and that the reconstructed `PB`, `T_INIT`, `MUB`, `ALB`, and `PHB`
are staged for later pressure-refresh consumers. Probe outputs, JSON reports,
and staging metadata must stay marked as diagnostic, not gate candidates, and
not integrator output. Provider success is not a pressure-refresh application,
because the compute path is still not wired into skeleton/remap, and it is not
a 10 min validation pass.

Later restart-file `ALB` or `PHB` remains acceptable only for smoke/range probe
comparisons. It must not be treated as d02 start-time truth or used to bypass
the provider/recompute path.

Round D29 defines the adapter/staging bridge boundary. The
`base_state_reconstruction_provider` may supply same-domain, same-valid-time
`ALB`, `PB`, `MUB`, and `PHB` buffers to pressure-refresh staging. That bridge
does not apply pressure refresh, does not modify exposed-cell `P`, and does not
make a diagnostic parent-fill or closure artifact gate-eligible. If
`pressure_refresh_applied = false`, the artifact must remain non-physical and
excluded from normal gates, regardless of provider or staging success.

Round D30 defines a skeleton/remap hook helper only with this explicit ordering:

```text
parent-fill/remap -> provider -> base-state sync -> staging -> pressure refresh compute -> report
```

The base-state sync step may copy provider `PB`, `MUB`, and `PHB` only into
newly exposed child cells. It must not rewrite overlap cells that came from
remap, and it must not touch halo cells governed by halo-update contracts.

Provider `ALB` is an external staging input for pressure-refresh compute, not a
normal `State` field. Provider `T_INIT` is base-state initial temperature, not
perturbation `T`; it must not be copied into `State::t` or treated as a closure
replacement for prognostic temperature.

The report step must preserve distinct parent-fill/remap, provider,
base-state-sync, staging, and pressure-compute markers. A helper-invoked compute
does not make the artifact validated integrator output. Until a real 10 min d02
gate passes, diagnostic hook reports must keep `gate_candidate=false` and
`integrator_output=false`.

Later restart `ALB` or `PHB` may still be used only for probe/smoke checks and
must never substitute for start-time truth.

Round D31 defines the opt-in skeleton/remap pressure-refresh diagnostic path.
`tywrf_skeleton_cycle` may call the D30 hook only when the diagnostic
parent-fill path is explicitly requested. This path is not a normal closure,
not a normal skeleton output, and not validated integrator output.

If `pressure_refresh_applied = true`, the only allowed interpretation is that
the diagnostic pressure-refresh helper compute was invoked and modified the
artifact's exposed `P` values. The flag does not prove that the full
parent-child remap, feedback, dynamics, physics, or 10 min validation gate
passed.

The CLI report and NetCDF attributes for this path must keep:

```text
TYWRF_GATE_CANDIDATE = false
TYWRF_INTEGRATOR_OUTPUT = false
TYWRF_VALIDATION_GATE_ONLY = false
```

Direct `ALB` remains optional only in the sense that WRF history files usually
omit it. If direct `ALB` is missing and the same-domain base-state
reconstruction inputs are available, the provider may stage `ALB` for the
pressure-refresh helper. Later restart `ALB` or `PHB` may be used only for
probe/smoke checks and must never substitute for start-time truth.

Any `--variables` narrowing or output-field selection used with this opt-in
path must retain the pressure-refresh required fields: the exposed `P` target,
the provider/staging fields needed by the helper (`PB`, `MUB`, `PHB`, staged
provider `ALB`), and the KROSA vertical/base-state constants (`P_TOP`, `C3F`,
`C4F`, `C3H`, `C4H`). Omitting one of these must fail before compute starts so
the helper cannot refresh uninitialized or partially staged state.

Round D32 records the closure-side interpretation of the real KROSA opt-in
pressure-refresh diagnostic smoke. `pressure_refresh_applied = true` is only
evidence that the diagnostic helper path invoked provider, staging, and
pressure-refresh compute on real KROSA d02 fields. It means the helper updated
the exposed diagnostic `P` values in that artifact; it does not mean the file
is physical, integrator-produced, or accepted by the 10 min d02 gate.

Strict gate tooling must keep requiring real candidate metadata. Any artifact
with the following metadata values remains non-gate and should be rejected or
reported as failed when submitted to the default strict gate:

```text
TYWRF_GATE_CANDIDATE = false
TYWRF_INTEGRATOR_OUTPUT = false
TYWRF_VALIDATION_GATE_ONLY = false
```

The R32 smoke is useful only because it exercises the provider/staging/compute
path against real KROSA d02 fields while verifying that the non-gate metadata
defense still prevents accidental validation credit. It must not use later
restart truth, must not promote later restart `ALB` or `PHB` to start-time
truth, must not lower d02 resolution below 2 km, and must not relax the
best-track-nudging prohibition.

Round D33 extends the closure boundary to gate metadata. The strict d02 gate now
requires positive integrator candidate metadata, such as
`TYWRF_GATE_CANDIDATE = true` and `TYWRF_INTEGRATOR_OUTPUT = true`. A closure,
oracle, diagnostic pressure-refresh, diagnostic remap, reference-copy, or WRF
end-state-delta artifact remains rejected even if a report can compute RMSE or
TC diagnostic numbers from it. Negative or missing gate metadata is a failure,
and diagnostic-only metadata cannot be reinterpreted as a pass.

`report_10min_diagnostics` may include gate-eligibility context for these
artifacts, but that context is only an audit trail. It does not satisfy the
strict gate, does not turn closure-derived `SLP` into an accepted TC diagnostic,
and does not permit an oracle candidate kind to pass.

`state_exchange` reports are likewise outside closure acceptance. With
`performed_interpolation = false`, they describe planned or counted
parent/child exchange only. They do not show that interpolation updated child
fields, do not close pressure or mass diagnostics, and do not make any closure
or diagnostic artifact gate-eligible.

Round D34 selected-field parent-to-child interpolation must remain outside the
diagnostic-closure mechanism. Its allowed inputs are start-state fields,
d01-derived parent data, and same-time KROSA constants. It must not use later
restart truth, WRF end-state fields, or WRF end-state deltas as an oracle. Until
the normal strict `00:10` validation gate passes with positive integrator
metadata, any selected-field interpolation report is implementation evidence
only and must not be labeled as a gate pass.

Round D35 keeps `selected_field_cycle` outside diagnostic closures while
allowing it to become the first non-oracle selected-field integrator candidate.
It may use d01/d02 start-state inputs only, then apply moving-nest remap,
`StateExchangePlan`, and `parent_child_interpolation` to produce all strict d02
gate fields. It is not WRF-exact physics and is expected at first to fail
strict thresholds numerically, not because a closure, oracle, metadata shortcut,
missing field, or non-finite sentinel hid the real error.

The following positive metadata belongs only to that non-oracle selected-field
integrator path:

```text
TYWRF_GATE_CANDIDATE = true
TYWRF_INTEGRATOR_OUTPUT = true
TYWRF_VALIDATION_GATE_ONLY = false
TYWRF_CANDIDATE_KIND = selected_field_integrator_v0
```

Diagnostic closures must not emit or retrofit those values. They are valid
only when no end-state/reference-end input was used, d02 remains `2 km`, every
strict field is finite, and selected fields changed from the cycle-start state
through remap / exchange / interpolation. `tendency_cycle`, skeleton,
diagnostic remap, diagnostic pressure-refresh, closure-modified, oracle, and
reference-copy artifacts remain explicitly non-gate and cannot be described as
selected-field gate candidates.

Round D35 real KROSA selected-field smoke is the counterexample that keeps the
closure boundary clear. `selected_field_cycle` generated a positive-metadata
candidate from `2025-07-26_00:00:00` start states only, and the strict metadata
gate passed, but the result still failed numerically. The first failure was
`U` normalized RMSE `0.117875`; `V`, `MU`, and `P` also failed, with `P`
normalized RMSE `0.907405` as the largest current blocker. `T`, `PH`, and
`QVAPOR` passed thresholds, and all strict RMSE fields had full finite
coverage.

This result must not trigger a diagnostic closure to patch the failed fields.
The remaining gaps after D35 were normal output/diagnostic work: accepted
`SLP` was unavailable, `U10` was absent in the v0 smoke, and `P` needed a real
producer rather than a closure or oracle correction. Round D36 evaluated that
boundary by preserving additional selected-field outputs, not by creating a
closure.

D36 output completeness is not a closure success and not a gate pass. Commit
`7d8d37c` makes `selected_field_cycle` write strict fields plus available
cycle-start-backed `PB`, `PHB`, `MUB`, `PSFC`, `U10`, `V10`, `T2`, `Q2`,
`RAINC`, and `RAINNC` by default. Metadata and derived SLP/rainfall side
artifacts no longer block report plumbing by themselves, but the actual KROSA
gate still fails: `U` normalized RMSE `0.117875`, `V` `0.134244`, `MU`
`0.133382`, and `P` `0.907405` fail; `T`, `PH`, and `QVAPOR` pass. TC
diagnostics still fail storm-center distance at `150.622 km`; MSLP error
`0.364 hPa` and Vmax error `0.769 m s-1` pass, but `U10`/`V10` field
comparisons fail because they are preserved start-state values, not real
surface-diagnostics output.

This result must not trigger a diagnostic closure to patch failed fields or
producer gaps. The positive selected-field metadata belongs only to the
non-oracle candidate generated by `selected_field_cycle`; closures, diagnostic
pressure-refresh outputs, remap-only outputs, derived SLP/rainfall copies,
report JSON, oracle/reference-copy products, and metadata-only helper outputs
must not emit, borrow, or retrofit:

```text
TYWRF_GATE_CANDIDATE = true
TYWRF_INTEGRATOR_OUTPUT = true
TYWRF_VALIDATION_GATE_ONLY = false
TYWRF_CANDIDATE_KIND = selected_field_integrator_v0
```

Round D37 pressure-refresh work remains outside diagnostic closures unless it
is explicitly requested as a diagnostic artifact, and even then it cannot make
a normal gate pass. A candidate-facing pressure-refresh path must use only
cycle-start/provider-backed same-domain inputs. Missing required inputs,
shape mismatches, non-finite staging, or any need for `00:10` reference-end
truth must abort rather than fall back to a closure, reference-end delta,
later restart truth, `PSFC` pressure shortcut, SLP proxy, or helper telemetry
that masks the candidate metadata.

The D37 opt-in smoke demonstrates why this boundary matters. The provider-backed
hook can be invoked inside `selected_field_cycle` and can write a candidate
without using reference-end truth, but the first real KROSA run worsened `P`
normalized RMSE from `0.907405` to `5.806413`. This failure must not be repaired
with a diagnostic closure or by reclassifying helper telemetry as physical
success; the next work is provider/base-state/static consistency, not a
metadata or closure shortcut.

Moving-nest static/coordinate handling and surface diagnostics are real
producer blockers, not closure opportunities. A closure must not copy
`00:10` reference-end `XLAT`, `XLONG`, `HGT`, or shifted-domain static truth
into a candidate, must not use best-track nudging, and must not create SLP or
`U10`/`V10` proxy shortcuts to satisfy storm-center, pressure, or Vmax gates.

Round D38 keeps those blockers outside closure scope:

- Static refresh is a moving-nest/static-field producer problem. A valid
  implementation may use the cycle-start d02 pose, the computed nest shift, and
  parent/static interpolation or reconstruction. It must not copy
  reference-end `XLAT`, `XLONG`, `HGT`, or any other `00:10` static truth into
  the candidate, even if that would reduce the G37 stale-coordinate errors.
- Pressure negative diagnosis must treat D37 as failed numerical evidence, not
  as a closure trigger. The opt-in path changed `1,053,150` `P` points and
  worsened `P` normalized RMSE from `0.907405` to `5.806413`; a closure must
  not patch that result, relabel helper telemetry as physical success, or make
  `--pressure-refresh` default.
- Surface ABI v2 must provide real diagnostics through SFCLAY or a physics
  wrapper. Preserved `U10`, `V10`, `T2`, `Q2`, `RAINC`, and `RAINNC` values are
  placeholders/provenance markers only. A closure must not treat them as
  producers or use them to satisfy surface field or Vmax gates.

The first D38 static-refresh result is deliberately not a closure. It improves
candidate coordinates through start-state pose reconstruction: `XLAT` RMSE
`0.000287 deg`, `XLONG` RMSE `0.000195 deg`, `HGT` RMSE `0.541 m`, and
storm-center error `43.483 km`. It does not copy reference-end static truth and
does not make the `00:10` gate pass; `U`, `V`, `MU`, `P`, and storm-center
thresholds still fail.

Round D39 keeps the post-D38 blockers out of closure scope:

- D38 static refresh, committed as `517ad28`, is a real coordinate/static
  producer improvement but not a gate pass. The improved `XLAT`, `XLONG`, and
  `HGT` metrics and the storm-center reduction to `43.483 km` must not be used
  to justify a closure patch for the still-failed `U`, `V`, `MU`, `P`, or
  storm-center thresholds.
- `--pressure-refresh` cannot be promoted into a closure or default fix. D37
  worsened `P` normalized RMSE to `5.806413`, and P38 showed around 10%
  same-source wrfout self-consistency error for the current refresh
  formula/staging. Any candidate-facing refresh needs an explicit readiness
  guard and self-consistency probe. A failed guard must stop or skip refresh,
  not route through pressure_mu_closure, reference-end deltas, later restart
  truth, or a metadata-only success marker.
- Exposed-cell `T`/`PH`/`P` same-pose analysis may diagnose whether the
  parent-fill, interpolation, base-state, or pressure path is inconsistent, but
  a closure must not provide the missing producer. The accepted fix must use
  non-oracle parent-child interpolation or recompute logic from cycle-start and
  same-time inputs. It must not copy, blend, or bias-correct with `00:10` WRF
  d02 truth.
- Surface ABI v2 is a scaffold until the SFCLAY or physics-wrapper producer
  actually runs. Closures must not treat ABI v2 structs, staging buffers,
  interface reports, or preserved cycle-start `U10`, `V10`, `T2`, `Q2`,
  `RAINC`, and `RAINNC` as executed physics or as a surface producer.

Round D40 keeps the D39 safety and ABI results outside closure scope:

- The `--pressure-refresh` readiness guard added in commit `0ea9a69` currently
  aborts the real KROSA smoke and produces no output. That abort is a safety
  constraint, not a pressure producer. A closure must not reinterpret a guarded
  abort as refreshed `P`, must not fill the missing output from
  `pressure_mu_closure`, and must not bypass the guard with reference-end
  pressure, `PSFC`, SLP, or later restart truth.
- Exposed-cell `T` and `PH` parent interpolation is an integrator producer
  slice, not a diagnostic closure. It may use non-oracle parent-child
  interpolation from cycle-start d02 state, d01-derived parent fields, and
  same-time KROSA constants. It must not use `00:10` WRF d02 truth, end-state
  deltas, bias corrections, or closure outputs as source data. `P` must remain
  outside direct parent interpolation and cannot be closure-patched into a
  production field.
- `PB`, `PHB`, and `MUB` are read-only base-state/provider-owned inputs for
  this slice. Closures must not overwrite them, repair them from later restart
  files, or treat a selected-field interpolation run as an owner of those
  fields. A future provider/sync path must document its ownership separately.
- The ABI v2 sidecar fixture and `_ex` scaffold are not executed physics. The
  current path reports `wrapper_unavailable` and `executed_physics = 0`, so a
  closure must not use sidecar fixture success to satisfy `U10`, `V10`, `T2`,
  `Q2`, `RAINC`, `RAINNC`, Vmax, or surface-producer validation.

Round D41/D42 keeps moved-terrain provider work outside closure scope:

- D41 commit `fe01be1` adds `KrosaBaseStateProviderTerrainOverride`, allowing
  the provider to reconstruct base-state buffers from a moved candidate `HGT`
  view. This is provider/probe plumbing only. A diagnostic closure must not
  treat provider reconstruction as pressure-refresh production or as a writer
  of `P`, `PB`, `PHB`, or `MUB`.
- D42 `--pressure-refresh` readiness work may inject moved candidate `HGT` into
  the provider as a read-only probe. The reconstructed buffers must not be
  synced into the candidate, used to overwrite base-state fields, or emitted as
  a diagnostic pressure-refresh artifact for gate credit.
- Even when provider terrain source is `moved_candidate_HGT`,
  `thermodynamic_base_state_consistency_ready=false`; `--pressure-refresh` must
  fail closed and not generate output. A closure must not fill that missing
  output from restart `PHB`/`ALB`, `pressure_mu_closure`, `00:10`
  reference-end truth, or end-state deltas.
- The default selected-field candidate must continue preserving
  `P`/`PB`/`PHB`/`MUB` start-state ownership. `P` is forbidden from direct
  parent interpolation. The current `00:10` gate is still failed on `U`, `V`,
  `MU`, `P`, and storm-center distance; probe or diagnostic pressure-refresh
  artifacts must not be reported as gate passes.

Round D43 keeps pressure-refresh hook terrain overrides outside closure scope:

- The pressure-refresh hook/API may let the provider consume a moved-terrain
  override. This is diagnostic hook capability only, not the default
  selected-field pressure producer and not a closure input.
- A hook unit test may call the terrain-override path and refresh exposed `P`;
  that demonstrates hook wiring only. It does not make
  `selected_field_cycle --pressure-refresh` readiness pass, and it does not
  count as a real `00:10` gate pass.
- The default selected-field candidate and D42 opt-in guard still require
  `thermodynamic_base_state_consistency_ready=false`, fail closed, no output
  file, and no writes to `P`, `PB`, `PHB`, or `MUB`. A closure must not fill in
  that missing output or reinterpret the hook smoke as candidate production.
- `P` remains forbidden from parent interpolation. Closures must still reject
  `00:10` reference-end truth, restart `PHB`/`ALB` as a provider substitute,
  and probe or diagnostic pressure-refresh artifacts as gate evidence. The
  current `00:10` gate remains failed on `U`, `V`, `MU`, `P`, and storm-center
  distance.

Round D44 defines the selected-field pressure-refresh readiness boundary:

- D41-D43 provide moved `HGT` provider override, selected-field moved-HGT
  provider probe, and hook terrain override. These are readiness capabilities,
  not closure producers and not selected-field gate evidence.
- The next stage must report a controlled readiness contract before
  `--pressure-refresh` can write output: provider terrain source/provenance,
  base-state reconstruction status, planned base-state and pressure sync
  counts, moved-nest overlap/halo safety, and pressure compute/requested/
  applied status.
- Until that contract is complete and positive, `--pressure-refresh` must fail
  closed with no output file. The default selected-field candidate must preserve
  `P`/`PB`/`PHB`/`MUB`, and a closure must not fill or reinterpret those fields
  from diagnostic artifacts.
- A hook diagnostic smoke may refresh exposed `P`, but it is hook-level evidence
  only. It must not be promoted into selected-field readiness or a real `00:10`
  gate pass.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  as a provider substitute, direct parent interpolation of `P`, and diagnostic
  artifact gate evidence. The current `00:10` gate remains failed on `U`, `V`,
  `MU`, `P`, and storm-center distance.

Round D45/D46 keeps hook dry-run reporting outside closure scope:

- D45 hook dry-run reporting may expose would-sync counts for `PB`, `MUB`, and
  `PHB`, overlap/halo write safety, `base_state_sync_applied=false`, and
  no pressure-compute status. Those fields describe a no-write/no-compute
  contract only; they are not closure outputs and must not change candidate
  state.
- D46 selected-field `--pressure-refresh` work may report that dry-run contract
  in the opt-in guard, but the guard must continue to fail closed while
  `thermodynamic_base_state_consistency_ready=false`. A closure must not create
  the missing output file or write `P`, `PB`, `PHB`, or `MUB` on behalf of that
  aborted path.
- Even a successful dry-run contract is readiness evidence, not a gate pass.
  Hook diagnostic smoke is also hook evidence only. Neither may be used to
  claim that the strict `2025-07-26_00:10:00` d02 gate passed.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  as a provider substitute, direct parent interpolation of `P`, and diagnostic
  artifact gate evidence. The current `00:10` gate remains failed on `U`, `V`,
  `MU`, `P`, and storm-center distance.

Round D47 keeps scratch pressure compute dry-run reporting outside closure
scope:

- The hook layer may run pressure compute against scratch `P` buffers only to
  report would-refresh `P` point counts, invalid `P` point counts, and a safety
  status. These values are telemetry, not closure-modified fields.
- The scratch dry-run must not write candidate `P`, `PB`, `PHB`, or `MUB`, must
  not generate selected-field output, and must not set
  `pressure_refresh_applied=true`.
- The selected-field tool layer is not connected in D47. The D46
  `selected_field_cycle --pressure-refresh` opt-in guard remains fail-closed
  while `thermodynamic_base_state_consistency_ready=false`.
- Scratch compute evidence is not a gate pass and must not be promoted into
  closure evidence. The real `2025-07-26_00:10:00` gate remains failed on `U`,
  `V`, `MU`, `P`, and storm-center distance.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  provider substitutes, direct parent interpolation of `P`, diagnostic artifact
  gate evidence, and any diagnostic artifact presented as selected-field
  integrator output.

Round D48 keeps selected-field scratch telemetry outside closure scope:

- The selected-field `--pressure-refresh` opt-in guard may report scratch
  pressure compute dry-run telemetry such as dry-run requested/called/ok,
  would-refresh `P` point counts, and invalid `P` point counts. These fields
  are readiness telemetry only, not closure outputs and not candidate values.
- While `thermodynamic_base_state_consistency_ready=false`, the opt-in guard
  must still fail closed with a nonzero exit and no output file. A closure must
  not create the missing output or reinterpret the abort report as a produced
  candidate.
- The scratch dry-run must not write candidate `P`, `PB`, `PHB`, or `MUB`, and
  must not set `pressure_refresh_applied=true`. The default selected-field
  numerical path remains unchanged.
- Scratch telemetry is not a gate pass. The real
  `2025-07-26_00:10:00` gate remains failed on `U`, `V`, `MU`, `P`, and
  storm-center distance.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  provider substitutes, direct parent interpolation of `P`, diagnostic artifact
  gate evidence, and scratch telemetry presented as gate evidence.

Round D49 keeps future apply-path terrain parity outside closure scope:

- D49 preparation only aligns the future selected-field pressure-refresh apply
  path with the dry-run path so both use the same moved candidate `HGT` terrain
  override when reconstructing provider-backed pressure-refresh staging.
- This terrain parity is not a closure output, not pressure-refresh output, and
  not evidence that `selected_field_cycle --pressure-refresh` produced a
  candidate file.
- While `thermodynamic_base_state_consistency_ready=false`, the opt-in path
  must still fail closed with a nonzero exit, no output file, no candidate
  writes to `P`, `PB`, `PHB`, or `MUB`, and no
  `pressure_refresh_applied=true`.
- The default selected-field numerical path remains unchanged and must not be
  closure-patched.
- Even if future apply-path terrain parity is ready, closures must not treat
  scratch telemetry or an unreachable guarded apply path as a gate pass. The
  real `2025-07-26_00:10:00` gate remains failed on `U`, `V`, `MU`, `P`, and
  storm-center distance.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  provider substitutes, direct parent interpolation of `P`, diagnostic artifact
  gate evidence, and scratch or unreachable-apply telemetry presented as gate
  evidence.

Round D50 keeps selected-field pressure-refresh apply readiness outside closure
scope:

- The current positive evidence is limited to the moved-`HGT` provider probe,
  base-state sync dry-run, scratch pressure compute dry-run, and future apply
  moved-`HGT` terrain parity.
- Those signals are readiness evidence only. They do not produce a closure, do
  not refresh candidate `P`, and do not make a selected-field candidate
  gate-eligible.
- D50 must not open `thermodynamic_base_state_consistency_ready`; while it is
  `false`, `--pressure-refresh` must still fail closed with nonzero exit, no
  output file, no candidate writes to `P`, `PB`, `PHB`, or `MUB`, and no
  `pressure_refresh_applied=true`.
- The default selected-field candidate numerical path remains unchanged and
  must not be closure-patched.
- Before a future apply path can be enabled, closure boundaries require a
  candidate mutation audit, moved-`HGT` apply proof, overlap/halo no-write
  proof, pre/post apply report consistency, and a real
  `2025-07-26_00:10:00` validation plan that keeps closure artifacts out of the
  gate.
- Scratch telemetry and unreachable apply plumbing must not be treated as a
  closure result or as gate evidence.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  provider substitutes, direct parent interpolation of `P`, diagnostic artifact
  gate evidence, and scratch or unreachable-apply telemetry presented as gate
  evidence.

Round D51 keeps hook postcondition and selected-field controlled seam evidence
outside closure scope:

- D51 does not open the default selected-field `--pressure-refresh` guard. The
  normal path must still fail closed with nonzero exit, no output file, no
  candidate writes, no `pressure_refresh_applied=true`, and no change to the
  default selected-field numerical path.
- Hook postcondition hardening is an apply safety check, not a closure. If real
  apply is reached in a controlled path, only `P`, `PB`, `MUB`, and `PHB` may
  change; overlap cells and halo cells must have zero writes; invalid,
  non-finite, inconsistent, or partial compute/write results must fail instead
  of being reported as successful apply.
- The controlled selected-field seam may prove only that tool-level apply uses
  moved candidate `HGT` rather than metadata terrain. A closure must not convert
  that seam, its reports, or experimental output into a pressure producer or
  gate evidence.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  provider substitutes, direct parent interpolation of `P`, diagnostic artifact
  gate evidence, scratch telemetry, unreachable apply plumbing, and controlled
  seam outputs presented as gate evidence.

Round D52 keeps transactional apply and selected-field report parity outside
closure scope:

- D52 still does not open the default selected-field `--pressure-refresh`
  guard. The normal path must fail closed with nonzero exit, no output file, no
  candidate writes, no `pressure_refresh_applied=true`, and no change to the
  default selected-field numerical path.
- Hook transactional apply is not a closure. The hook must compute pressure
  refresh into scratch buffers first, accept the postcondition there, and only
  then write allowed candidate fields. If the postcondition fails, the candidate
  must remain completely unchanged.
- Selected-field experimental apply report parity is audit evidence only. It
  may record target/refreshed/skipped/invalid counts, overlap and halo flags,
  changed counts, and consistency flags, but those fields do not produce a
  closure and do not make the output gate-eligible.
- Hidden seam output remains diagnostic-only, non-gate, and not normal
  integrator output. A closure must not convert it, its metadata, or its report
  into pressure production or gate evidence.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  provider substitutes, direct parent interpolation of `P`, diagnostic artifact
  gate evidence, and any transactional/report-parity or hidden-seam artifact
  presented as a gate pass.

Round D53 keeps strict normal pressure-refresh promotion outside closure scope:

- D53 may replace the hard-coded
  `thermodynamic_base_state_consistency_ready=false` guard only with a derived
  readiness decision from normal-path provider, dry-run, and transactional
  apply evidence. A closure, probe, hidden seam, or report-only artifact cannot
  supply that readiness.
- A normal selected-field `--pressure-refresh` candidate may be written only
  when static refresh is start-state based, provider terrain is
  `moved_candidate_HGT`, base-state sync dry-run passes, pressure compute
  dry-run passes, overlap/halo/invalid/skipped point checks are clean, and the
  transactional real apply succeeds.
- The hidden `--experimental-pressure-refresh-apply` seam remains
  diagnostic-only, non-gate, and non-integrator output. A closure must not
  convert it, its metadata, or its report into pressure production,
  gate-eligible candidate metadata, or validation evidence.
- Candidate eligibility is not a real KROSA gate pass. After any D53 normal
  candidate is produced, the first validation must be the strict
  `2025-07-26_00:10:00` gate. If it fails, reports must name the first failed
  field or diagnostic and must not proceed to `00:20`.
- Closures must still reject `00:10` reference-end truth, restart `PHB`/`ALB`
  provider substitutes, direct parent interpolation of `P`, diagnostic/probe/
  hidden artifacts, and any readiness telemetry presented as gate evidence.

Round D54 records the post-D53 failure boundary:

- D53 moved normal `selected_field_cycle --pressure-refresh` from fail-closed
  to gate-eligible only under strict same-artifact readiness and transactional
  apply. That is not a closure success and not a validation pass.
- The D53 real KROSA d02 `2025-07-26_00:10:00` gate failed after metadata
  passed. First failure was `U` normalized RMSE `0.117875`; `V` `0.134244`,
  `MU` `0.133382`, and `P` `6.33495` also failed. `T` `0.013121`, `PH`
  `0.017410`, and `QVAPOR` `0.017017` passed. Storm-center error
  `43.483 km` failed; minimum SLP error `0.364 hPa` and Vmax error
  `0.769 m s-1` passed.
- D54 must diagnose the `U` first failure and the pressure-refresh numerical
  `P` error before any `00:20` progression. A closure must not patch `U`, `V`,
  `MU`, `P`, or storm center, and must not reinterpret pressure-refresh apply
  contract/parity as WRF-compatible `P`.
- Audit tools and diagnostic reports are observation-only. They may explain
  metadata eligibility, provenance, changed fields, or failures, but they must
  not produce candidates, repair metadata, or count as gate evidence.
- The hidden `--experimental-pressure-refresh-apply` seam remains
  diagnostic-only, non-gate, and non-integrator output. Its reports and
  metadata cannot be converted into closure evidence, candidate evidence, or a
  `00:10` gate pass.

Round D55 keeps the D54 audit conclusions outside closure scope:

- The `P` error is concentrated in the refreshed target region. Global `P`
  normalized RMSE is `6.334952`, target-region `P` normalized RMSE is
  `10.33536`, non-target `P` normalized RMSE is `0.338162`, and the target
  error fraction is `0.998219`.
- `PB`, `MUB`, and `PHB` pass or remain close in the audit, so a closure must
  not rewrite those base-state fields to hide a perturbation-pressure refresh
  semantics problem.
- Derived SLP can pass with normalized RMSE `0.003248`, but this does not
  override the failed strict `P` field, the failed `U`/`V`/`MU` fields, or the
  failed storm-center diagnostic.
- D55 diagnostics may audit pressure-refresh formula/source/vertical/region
  semantics and may split selected-field `U`, `V`, and `MU` errors by
  remapped-overlap, target/exposed, and non-target regions.
- Those audits are observation-only. They may compare against WRF reference
  output to explain errors, but they must not produce candidate fields, patch a
  candidate, tune a formula using reference-end truth, or convert hidden-seam
  output into gate evidence.

Round D56 keeps the integrated D55 outcome outside closure scope:

- D55 commit `cc39cdb` passed code validation (`cmake --build`, CTest `29/29`,
  pytest `141/141`), but that is not a scientific validation pass. The strict
  d02 `2025-07-26_00:10:00` gate remains failed and validation must not advance
  to `00:20`.
- The pressure formula/source audit remains diagnostic-only. The refreshed
  target-region perturbation-pressure error is still the dominant `P` issue:
  global `P` normalized RMSE `6.334952`, target-region `10.335362`,
  non-target `0.338162`, target error fraction `0.998219`, and target bias
  about `-805 Pa`.
- Full-pressure normalized RMSE can appear small because the full-pressure
  reference scale is much larger, but raw RMSE is nearly unchanged. A closure
  must not replace the strict perturbation-pressure `P` gate with a
  full-pressure metric.
- Candidate `P` is not close to start `P`; a closure must not relabel that as
  persistence or use it to justify a pressure patch.
- The selected-field state audit remains diagnostic-only. The first failing
  variable is `U`, with selected-field normalized RMSE values `U = 0.117875`,
  `V = 0.134244`, and `MU = 0.133382`. Target-region fractions below `0.5`
  and candidate fields not close to start fields point to a broader
  selected-field/dynamics gap, not an isolated exposed-cell interpolation hole.
- D56 may add a vertical/source-level target-region pressure audit and a
  selected-field evolution audit against WRF reference evolution from the start
  state. These audits may observe and attribute errors only; they must not
  generate closure candidates, patch fields, tune formulas from reference-end
  truth, or convert diagnostic/probe/hidden-seam output into gate evidence.

Round D57 keeps the D56 audit results outside closure scope:

- D56 commit `985a38f` passed code validation (`cmake --build`, CTest `29/29`,
  pytest `150/150`) and pushed to `origin/main`, but that is still not a
  scientific validation pass. The strict d02 `2025-07-26_00:10:00` gate remains
  failed and validation must not advance to `00:20`.
- The pressure source audit can keep a failed top-level status because real
  start `wrfout` source entries still miss `ALB`. A closure must not interpret
  that source-seed audit status as a model pass/fail override or use it to
  bypass the normal strict gate.
- The pressure diagnostic subreport shows target-region `P` RMSE
  `1560.573 Pa`, normalized RMSE `10.335362`, and mean bias `-805.059 Pa`.
  The worst levels by RMSE and bias are low-level mass levels `0..4`, with
  level `0` RMSE `4170.444 Pa` and bias `-4167.791 Pa`; the source/start
  evolution fraction is `6.57739`. This points to a low-level target-region
  pressure formula/source issue, not an invitation to patch `P` with a closure.
- The selected-field evolution audit reports evolution normalized RMSE
  `U = 0.156800`, `V = 0.160006`, and `MU = 0.197589`, with target-region
  fractions `U = 0.414932`, `V = 0.362148`, and `MU = 0.329146`. Amplitude
  ratios are close to `1` and capture fractions are high but imperfect, so a
  closure must not treat the failure as amplitude-only or target-region-only.
- D57 may add a pressure vertical-bias companion-field attribution audit and a
  selected-field spatial alignment/shift diagnostic audit. Both are
  observation-only; they cannot patch fields, generate candidates, tune
  formulas from reference-end truth, or turn diagnostic/audit/probe/hidden-seam
  artifacts into gate evidence.

Round D58 keeps the D57 audit results outside closure scope:

- D57 commit `ee2fff4` passed code validation (`cmake --build`, CTest `29/29`,
  pytest `160/160`) and pushed to `origin/main`, but that is still not a
  scientific validation pass. The strict d02 `2025-07-26_00:10:00` gate remains
  failed and validation must not advance to `00:20`.
- The pressure vertical-bias audit shows the low-level `P` error is not
  explained by base-state companion fields. The worst mass levels are `0..4`;
  level `0` has `P` mean difference `-4167.791 Pa` and RMSE `4170.444 Pa`.
  `PB` and `MUB` companion differences are about `+1 Pa`, `MU` is about
  `+11.625`, `P + PB` inherits the negative `P` bias, and `PSFC` is about
  `-376.164 Pa`. A closure must not patch that perturbation-`P`
  producer/staging problem or relabel it as base-state agreement.
- The selected-field spatial alignment audit found only modest normalized RMSE
  gains from small shifts: `U` end `2.46%` and evolution `2.64%`, `V` end
  `3.40%` and evolution `3.79%`, and `MU` end `7.57%` and evolution `7.28%`.
  A closure must not treat the selected-field failures as solved by a simple
  horizontal shift or use a shift diagnostic to move candidate fields.
- D58 may audit local C++/metadata pressure producer and staging paths, plus
  selected-field pipeline, timing, and remap paths using local metadata and the
  D56/D57 JSON summaries. Those audits are observation-only: they cannot
  generate closure candidates, patch fields, tune formulas from reference-end
  truth, or convert diagnostic/audit/probe/hidden-seam artifacts into gate
  evidence.

Round D59 keeps the D58 audit results outside closure scope:

- D58 commit `5d9adb5` passed code validation (`cmake --build`, CTest `29/29`,
  pytest `167/167`) and pushed to `origin/main`, but that is still not a
  scientific validation pass. The strict d02 `2025-07-26_00:10:00` gate remains
  failed and validation must not advance to `00:20`.
- The pressure producer audit
  `build/validation/r58_pressure_refresh_producer_audit.json` is
  diagnostic-only. It identifies local producer inventory around
  `compute_krosa_pressure` (`src/dynamics/pressure_refresh.cpp:275`),
  `refresh_krosa_moving_nest_pressure`
  (`src/dynamics/pressure_refresh.cpp:357`), and
  `apply_krosa_moving_nest_pressure_refresh_hook`
  (`src/dynamics/pressure_refresh_hook.cpp:309`). The worst level remains
  mass level `0`, with maximum absolute `P` mean difference about
  `4167.790767 Pa`. A closure must not patch `P`, relabel this as a gate
  result, or tune a formula from reference-end truth.
- The selected-field pipeline audit
  `build/validation/r58_selected_field_pipeline_audit.json` is also
  diagnostic-only. It reports seven risk flags: target fraction below half,
  target error fraction below half, amplitude near one while RMSE remains
  high, modest best-shift improvement, end/evolution shift mismatch,
  gate/integrator claims while audits fail, and large movement-delta
  schedule/remap risk. The movement child delta is `{i:60,j:35}`, with target
  fractions about `U = 0.407583`, `V = 0.407583`, and `MU = 0.404762`. A
  closure must not treat these timeline/remap diagnostics as permission to
  shift, blend, or patch candidate fields.
- D59 may add a pressure formula/input audit and a selected-field
  schedule/remap timeline audit. Both are observation-only; they cannot
  generate closure candidates, patch fields, tune formulas from reference-end
  truth, convert diagnostic/probe/hidden-seam outputs into gate evidence,
  reduce d02 below `2 km`, or introduce best-track nudging.

Round D60 keeps the D59 audit results outside closure scope:

- D59 commit `957c3d6` passed code validation (`cmake --build`, CTest `29/29`,
  pytest `175/175`) and pushed to `origin/main`, but that is still not a
  scientific validation pass. The strict d02 `2025-07-26_00:10:00` gate remains
  failed and validation must not advance to `00:20`.
- The pressure formula/input audit
  `build/validation/r59_pressure_formula_inputs_audit.json` is
  diagnostic-only, with status `computed_with_flags` and six risk flags. It
  reports `P` target-region normalized RMSE `10.335362`, target mean
  difference `-805.058863 Pa`, and low level `0` `P` mean difference
  `-4167.790767 Pa` with RMSE `4170.444038 Pa`. `P + PB` at low level `0`
  inherits the bias with mean difference `-4166.786671 Pa`, while companion
  fields remain much smaller: `PB` target normalized RMSE `8.76e-05`,
  `MU + MUB` `0.000878`, `PH + PHB` `0.000805`, `T` `0.012115`, and `QVAPOR`
  `0.016169`. `HGT` target normalized RMSE is `0.495197`, but its mean
  difference is only `-0.085612 m`. A closure must not patch pressure, tune a
  formula, or generate candidates from these diagnostic differences.
- The selected-field timeline audit
  `build/validation/r59_selected_field_timeline_audit.json` is also
  diagnostic-only, with status `computed_with_flags` and six risk flags. It
  reports movement child delta `{i:60,j:35}`, large movement true, `P` not
  listed as parent-interpolated, changed/interpolated count mismatch true, and
  target fractions about `0.407583` for `U` and `V`, and `0.404762` for `T`,
  `PH`, `MU`, `P`, and `QVAPOR`. A closure must not use this telemetry to
  shift, blend, patch, or relabel candidate fields as passing.
- D60 may add pressure per-column formula-term probes and selected-field
  runtime schedule/timeline metadata emission. Both are observation-only; they
  cannot generate closure candidates, patch fields, tune formulas from
  reference-end truth, convert diagnostic/probe/hidden-seam outputs into gate
  evidence, reduce d02 below `2 km`, or introduce best-track nudging.

Round D61 keeps the D60 telemetry results outside closure scope:

- D60 commit `7903cbf` passed full code validation (`cmake --build`, CTest
  `29/29`, pytest `178/178`) and pushed to `origin/main`, but that is still
  not a scientific validation pass. The strict d02
  `2025-07-26_00:10:00` gate remains failed and validation must not advance to
  `00:20`.
- The normal D60 `--pressure-refresh` candidate output
  `build/validation/r60_pressure_normal/wrfout_d02_2025-07-26_00:10:00`
  emitted runtime timeline metadata with 11 events, from `cycle_start` through
  `output_write_preparation`, including movement delta `(60,35)` and
  pressure-refresh refreshed, synced, and changed counts. It did not modify
  accepted numerical fields relative to the D54 normal candidate: max absolute
  difference is `0.0` for `U`, `V`, `T`, `PH`, `MU`, `P`, `QVAPOR`, `PB`,
  `PHB`, `MUB`, and `HGT`. A closure must not treat metadata emission as a
  repair or as gate evidence.
- The pressure column probe
  `build/validation/r60_pressure_formula_column_probe_audit.json` is
  diagnostic-only. It selected worst columns `(160,49)`, `(160,50)`,
  `(160,51)`, `(161,49)`, and `(161,50)`. Rank 1 level `0` `P` differs by
  `-4686.643066 Pa` (`-4172.158691 Pa` candidate versus `514.484375 Pa`
  reference), and levels `0..4` have mean absolute low-level `P` error
  `4391.227588 Pa`. A closure must not patch those columns, tune formula terms
  from reference-end truth, or use the probe as an accepted pressure producer.
- The selected-field timeline audit
  `build/validation/r60_selected_field_timeline_audit.json` confirms timeline
  attributes are present but still reports D59 risks until the parser fully
  structures event strings. The strict gate report
  `build/validation/r60_pressure_normal_gate.json` still fails with first
  failing field `U` normalized RMSE `0.117875`, plus `V = 0.134244`,
  `MU = 0.133382`, `P = 6.334952`, and storm-center error `43.482716 km`. A
  closure must not relabel these failing metrics as accepted through
  telemetry, hidden-seam output, or diagnostic probes.
- D61 may add optional same-column runtime observation metadata and structural
  timeline parsing. Both are telemetry/audit-only; they cannot generate closure
  candidates, patch fields, tune formulas from reference-end truth, convert
  diagnostic/probe/hidden-seam outputs into gate evidence, reduce d02 below
  `2 km`, or introduce best-track nudging.

Round D62 keeps the D61 runtime probe and timeline audit outside closure scope:

- D61 local commit `ab68c51`
  (`Add runtime pressure probes and timeline audit`) passed local full code
  validation (`cmake --build`, CTest `29/29`, pytest `183/183`, and
  `git diff --check`), but it is not pushed. GitHub SSH port `22` and
  SSH-over-443 push attempts timed out, so `origin/main` remains at D60 commit
  `7903cbf`. A closure document or report must not imply D61 is synchronized to
  origin until push succeeds.
- The D61 real KROSA smoke generated
  `build/validation/r61_pressure_normal/wrfout_d02_2025-07-26_00:10:00` with
  12 runtime timeline events. `pressure_column_probe` is ordered between
  `pressure_refresh_apply` and `cycle_end`. The probe covers five columns
  `(160,49)`, `(160,50)`, `(160,51)`, `(161,49)`, and `(161,50)`, five levels
  `0..4`, two phases `post_static_refresh` and `post_pressure_refresh`, and
  `50` records. It observes the mismatch location but does not close pressure
  or validate the integrator.
- D61 and D60 candidates are numerically identical for `U`, `V`, `T`, `PH`,
  `MU`, `P`, `QVAPOR`, `PB`, `PHB`, `MUB`, and `HGT`, with max absolute
  difference `0.0`. Telemetry therefore must not be described as a field
  repair.
- The strict `00:10` gate still fails: first failing field `U` normalized RMSE
  `0.117875`, with `V = 0.134244`, `MU = 0.133382`, `P = 6.334952`, and
  storm-center error `43.482716 km`. A closure must not relabel this as
  accepted or advance validation to `00:20`.
- D62 may add pressure-core formula observation for runtime POD terms
  (`ALB`, `C3F`, `C4F`, `C3H`, `C4H`, `P_TOP`, `theta`, and related
  pressure-refresh intermediates) and a Python pressure-column runtime probe
  audit. Both are diagnostic-only. They cannot create closure candidates,
  patch `P`, tune pressure from reference-end truth, promote probe/audit/
  hidden-seam output to gate evidence, reduce d02 below `2 km`, or introduce
  best-track nudging.

Round D62 is now synchronized to the remote repository:

- Commit `655df5b` (`Add pressure formula observation audit`) is pushed to
  `origin/main`, and the same push also delivered the previously blocked D61
  commit `ab68c51`.
- D62 full validation passed CTest `29/29`, pytest `188/188`, and
  `git diff --check`.
- The D62 runtime probe audit
  `build/validation/r62_pressure_column_runtime_probe_audit.json` parsed
  `50/50` records consistently. It reports
  `post_pressure_refresh_p_negative`, `large_p_drop_magnitude`, and
  `formula_terms_unavailable`; the minimum post-refresh `P` is
  `-4174.87305 Pa`, and the largest drop is `(161,49,k=0)` with
  `delta_P = -4257.17773 Pa`.

Round D63 may connect the D62 core formula observation to the selected-field
`--pressure-column-probe` NetCDF/JSON output and extend the audit parser for
the new formula observation attributes. This is still outside closure scope:
formula probes and audits cannot become accepted pressure closures, cannot
patch fields, cannot tune output from reference-end truth, cannot use
reference-end oracle candidate generation, cannot promote audit/probe/
hidden-seam output to gate evidence, cannot advance validation to `00:20`,
cannot reduce d02 below `2 km`, and cannot introduce best-track nudging.

Round D63 is complete but remains outside closure scope:

- Commit `3e1f6a7` (`Add selected-field formula observation telemetry`) is
  pushed to `origin/main`.
- Full code validation passed CTest `29/29`, pytest `195/195`, and
  `git diff --check`.
- The real candidate
  `build/validation/r63_pressure_formula_observation/wrfout_d02_2025-07-26_00:10:00`
  and audit
  `build/validation/r63_pressure_formula_observation/runtime_probe_audit.json`
  report formula observation present, `25` formula records, all recorded and
  valid, and `0` pressure mismatches against post-refresh probe `P`.
- Negative post-refresh `P` is therefore traced to the observed formula and
  base-pressure terms, not to a NetCDF writer mismatch, audit parser mismatch,
  or hidden probe disagreement.
- The strict `00:10` gate still fails: first failing field `U` normalized RMSE
  `0.117875`, with `V = 0.134244`, `MU = 0.133382`, `P = 6.334952`, and
  storm-center error `43.482716 km`.
- D63 candidate fields `U`, `V`, `T`, `PH`, `MU`, `P`, `QVAPOR`, `PB`, `PHB`,
  `MUB`, and `HGT` are numerically identical to D61, with max absolute
  difference `0.0`.

Round D64 is complete but remains outside closure scope:

- Commit `2db9279` (`Add pressure budget runtime audit`) is pushed to
  `origin/main`.
- Full validation passed CTest `29/29`, pytest `196/196`, and
  `git diff --check`.
- This is diagnostic-only completion, not a strict d02 gate pass.
- The real audit
  `build/validation/r64_pressure_budget_runtime_probe_audit.json` on the D63
  candidate has `status = computed_with_flags`.
- The audit reports `25` pressure budget records, `25`
  `total_pressure < PB` records, `25` records where the large drop is
  explained by formula/base-pressure subtraction, and pressure mismatch count
  `0`.
- The first tracked point `(i=160,j=49,k=0)` reports
  `total_pressure = 95539.724 Pa`, `PB = 99711.8828 Pa`,
  `total_pressure_minus_pb = -4172.1588 Pa`, and
  `probe_delta_p = -4252.3540025 Pa`.
- The source audit says TyWRF and WRF runtime both use
  `P = total_pressure - PB`, so a closure must not remove `PB` subtraction or
  rewrite pressure to hide this diagnostic result.

Round D65 may add a formula sensitivity diagnostic and a WRF
moving-nest/base-state call-order audit. These diagnostics should focus on
`PH + PHB`, `theta`, `MU + MUB`, and base-state staging consistency before any
formula correction. They remain outside closure acceptance: they cannot patch
pressure, rewrite formula terms, tune from reference-end truth, generate
reference-end or oracle candidates, convert audit/probe/hidden-seam output into
gate evidence, advance validation to `00:20`, reduce d02 below `2 km`, or
introduce best-track nudging.

Round D65 is complete but remains outside closure scope:

- Commit `07908b1` (`Add pressure formula sensitivity audit`) is pushed to
  `origin/main`.
- Full validation passed CTest `29/29`, pytest `197/197`, and
  `git diff --check`.
- This is diagnostic-only completion, not a strict d02 gate pass.
- The real audit
  `build/validation/r65_pressure_formula_sensitivity_audit.json` reports
  `25` sensitivity records.
- The `PB - total_pressure` gap ranges from `3593.9011` to `4174.8731 Pa`,
  with mean `3904.1623 Pa`.
- The approximate fractional total-pressure increase needed ranges from
  `0.0389068` to `0.0436995`, with mean `0.0414497`.
- All `25` records require a large total-pressure increase.

A closure must not treat the D65 sensitivity result as authority to rescale
total pressure, remove `PB` subtraction, tune formula terms from the
`00:10` reference, or write a patched gate candidate. The strict d02 `00:10`
gate remains failed, so validation must not advance to `00:20`.

Round D66 remains a provenance and mask-semantics audit direction, not a
closure or selected-field numerical opening. B65 concluded that WRF's broad
base-state interpolation includes `PHB`, `MUB`, `PB`, `ALB`, `T_INIT`, and
`HT`, while TyWRF's selected-field path currently interpolates only `U`, `V`,
`T`, `PH`, `MU`, and `QVAPOR`. D66 may document and audit the exposed d02
base-field provenance contract and WRF generated interpolation mask semantics,
but it must not generate reference-end or oracle candidates, convert
audit/probe/hidden-seam output into gate evidence, reduce d02 below `2 km`, or
introduce best-track nudging.

At D68 start, D66 commit `09d6ba2`
(`Add moving-nest base-state exchange contract`) and D67 commit `5a78fbb`
(`Add base-state exchange action diagnostics`) have both been pushed.
`main` and `origin/main` are synchronized at `5a78fbb`. D66 validation before
the local commit passed CTest `29/29` and pytest `197/197`; D67 did not make
the strict d02 `00:10` gate pass.

The D66 contract remains outside diagnostic-closure authority. It records
selected fields `U`, `V`, `MU`, `QVAPOR`, `T`, and `PH`; base-state candidates
`PHB`, `MUB`, `PB`, `ALB`, `T_INIT`, `HT`, and `HGT`; no selected-field
interpolation for those base-state candidates; and no `P` base-state candidate.
`PHB`, `MUB`, and `PB` are State-backed candidates, while `ALB`, `T_INIT`,
`HT`, and `HGT` are static/provider-backed candidates. The contract is marked
`diagnostic_only = true` and `enables_selected_field_numerics = false`.
Exposed-child exchange planning counts future interpolation work without
writing fields: `performed_interpolation = false`, `modifies_overlap = false`,
and `modifies_halo = false`.

B66 found that WRF generated unstaggered moving-nest mask semantics
(`imask_nostag`) apply to `PHB`, `T_INIT`, `MUB`, `ALB`, `PB`, and `HT`, and
that WRF `start_domain` recompute rules must remain distinct from closure
patches or selected-field numerical interpolation. A closure must not use those
facts to overwrite base fields, hide missing exposed-mask actions, or repair a
candidate with reference-end truth.

D67 closure boundary: the opt-in diagnostic moving-nest base-field
provenance/action report is an audit artifact only. It exposes base-field
source/action categories for newly exposed d02 cells and prepares WRF-style
mask geometry checks. Action labels such as
`interpolate_exposed_cells`, `recompute_from_mub_after_interpolation`,
`preserve_interpolated_when_rebalance_zero`, and `static_height_input` are
diagnostic report labels, not closure permissions.

D68 may add only a diagnostic exposed base-state exchange helper/regression.
`PHB`/`MUB`/`HT` exposed interpolation and `PB`/`T_INIT`/`ALB` recompute marks
are helper semantics, not closure authority, not production selected-field
numerics, and not gate evidence. A closure must not use them to patch exposed
base fields, hide missing parent-fill work, or convert a diagnostic artifact
into a candidate. D68 must not use reference-end truth, generate oracle
candidates, promote diagnostic/probe/helper output to gate evidence, advance
validation to `00:20`, reduce d02 below `2 km`, or introduce best-track
nudging. Only D69 or later should consider migrating `state_remap`/parent-fill
semantics into a WRF-style exposed base-state policy, and that migration must
remain outside diagnostic-closure shortcuts.

## Hard Prohibitions

The following schemes are forbidden:

- using a WRF reference end-state delta, such as
  `reference_end - reference_start`, to update candidate `MU`, `P`, `PSFC`,
  `SLP`, or any other field;
- treating WRF end-state/oracle output, reference-copy files, or output without
  applied pressure refresh as a validation-gate pass;
- copying, blending, nudging toward, bias-correcting against, or otherwise
  patching metrics with WRF cycle-end reference output;
- copying `00:10` reference-end `XLAT`, `XLONG`, `HGT`, or shifted-domain
  static truth to implement moving-nest static refresh;
- using reference output to move the TC center, adjust the minimum SLP, or tune
  pressure ranges;
- silently replacing normal skeleton or integrator fields with closure fields;
- treating preserved cycle-start `U10`, `V10`, `T2`, `Q2`, `RAINC`, or
  `RAINNC` as a real surface-diagnostics producer;
- allowing a closure-modified file to satisfy the normal validation gate
  without an explicit diagnostic-closure mode;
- emitting closure fields without metadata that marks them as non-physical
  diagnostic artifacts.

## Allowed Inputs

The minimum allowed input set is the cycle-start state only:

- `MU`
- `P`
- `PB`
- `PH`
- `PHB`
- `T`
- `QVAPOR`
- `PSFC`

Optional metadata such as domain id, valid time, grid spacing, and dimension
names may be read to preserve WRF-compatible shape and reporting context. The
closure may also use fixed KROSA constants already encoded in project metadata,
such as the 50 hPa model top, but it must not read a cycle-end WRF reference
file or any field derived from one.

If an implementation starts from a TyWRF candidate file, the candidate file is
allowed only as the container whose metadata and non-pressure fields are copied
forward. The pressure closure itself must remain computable from the cycle-start
fields above plus explicit constants. If candidate pressure tendencies become an
allowed input later, that change must be documented here before implementation.

## Mathematical Boundary

`pressure_mu_closure` may use one of two local column approximations.

### Hydrostatic Column Consistency

For each horizontal column, define:

```text
theta = T + 300 K
p_start = P + PB
z_w = (PH + PHB) / g
```

The closure may solve or approximate a one-dimensional pressure profile that is
consistent with hydrostatic balance in that column:

```text
dPhi ~= -alpha dp
alpha ~= R_d * T_v / p
T_v ~= theta * (p / p0)^(R_d / c_p) * (1 + 0.61 * QVAPOR)
```

The profile may be anchored by cycle-start `PSFC`, the cycle-start pressure
profile, and the inferred top pressure. The resulting perturbation pressure is:

```text
P_closed = p_closed - PB_start
```

`MU_closed` may be derived only as the dry column pressure-thickness quantity
implied by the same column approximation. Without the full WRF eta-coordinate
mass state and real tendencies, this is a consistency diagnostic, not a WRF
mass update.

### Local Mass-Pressure Closure

As a simpler fallback, the closure may preserve the cycle-start vertical
pressure shape and apply a local column correction that makes the reported
`MU`, `P`, and `PSFC` mutually consistent under a documented dry-column pressure
thickness relation. The correction must be local to each column and must not use
horizontal storm displacement, best-track data, cycle-end reference fields, or
reference-derived deltas.

## Limitations

The closure has no wind, moisture, microphysics, radiation, PBL, surface-layer,
land-surface, slab-ocean, lateral-boundary, spectral-nudging, or nesting
dynamics. It cannot forecast storm motion or intensity. It may reduce internal
pressure inconsistencies in a diagnostic candidate, but any agreement with a
WRF cycle-end field is not evidence that the integrator reproduced WRF physics.

The closure does not promise improvements for:

- `U`
- `V`
- `QVAPOR`
- hydrometeors
- accumulated rainfall
- 10 m wind
- storm center
- Vmax

If an opt-in report computes a pressure-minimum location from closure `SLP`,
that location must be labeled as a closure SLP minimum, not as the accepted TC
storm center.

## Observable Metrics

An opt-in diagnostic report may observe:

- `MU` normalized RMSE against the same valid-time WRF reference;
- `P` normalized RMSE against the same valid-time WRF reference;
- derived `SLP` minimum error, when both reference and closure diagnostic SLP
  are present in the diagnostic report;
- `PSFC` range and `SLP` range sanity checks, including min/max and finite
  sample counts.

The report must not claim pass/fail credit for `U`, `V`, `QVAPOR`, storm-center
distance, Vmax, or rainfall from this closure.

## Metadata Requirements

Every closure output file must include global metadata equivalent to:

```text
diagnostic_closure = "pressure_mu_closure"
diagnostic_candidate = "true"
not_physical = "true"
validation_gate_only = "true"
excluded_from_default_gate = "true"
uses_reference_end_delta = "false"
closure_allowed_inputs = "cycle_start:MU,P,PB,PH,PHB,T,QVAPOR,PSFC"
closure_mode = "hydrostatic_column" or "local_mass_pressure"
closure_source_valid_time = "<cycle-start time>"
closure_output_valid_time = "<reported valid time>"
```

Each variable written or modified by the closure must carry variable metadata:

```text
diagnostic_closure = "pressure_mu_closure"
not_physical = "true"
closure_role = "modified" or "derived"
closure_source = "cycle_start"
```

If the output also copies untouched non-pressure variables from an input
candidate container, those variables must not be marked with
`closure_role = "modified"`. The file-level metadata still marks the whole file
as a diagnostic candidate so validation tools do not treat it as normal
integrator output.

## Future API Draft

A future Python-only implementation should be opt-in and explicit, for example:

```text
apply_pressure_mu_closure(
    cycle_start_path,
    output_path,
    valid_time,
    domain,
    mode="hydrostatic_column",
    container_path=None,
    write_fields=("MU", "P", "PSFC", "SLP"),
)
```

Expected behavior:

- reject missing or shape-incompatible required inputs;
- reject cycle-end WRF reference paths and any `reference_end_delta` argument;
- never overwrite the normal integrator output in place by default;
- write the metadata listed above before the file can be used by validation
  tooling;
- preserve d02 2 km metadata when operating on d02 files;
- report touched variables and input source times in machine-readable JSON.

## Future Test Plan

The first implementation should add tests for:

- metadata presence on the file and every modified or derived variable;
- hard rejection of missing required inputs and shape mismatches;
- hard rejection of reference-end files, reference-copy candidates, and any
  reference-end delta option;
- deterministic output from the same cycle-start inputs, independent of which
  WRF cycle-end reference file is available nearby;
- synthetic hydrostatic columns with finite, positive, vertically monotonic
  pressure;
- local mass-pressure closure sanity for `MU`, `P`, and `PSFC` consistency;
- validation reporting that keeps closure metrics in an opt-in diagnostic block
  and prevents them from satisfying the normal d02 cycle gate.
