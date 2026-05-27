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
