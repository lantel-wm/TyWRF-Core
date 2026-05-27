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

## Hard Prohibitions

The following schemes are forbidden:

- using a WRF reference end-state delta, such as
  `reference_end - reference_start`, to update candidate `MU`, `P`, `PSFC`,
  `SLP`, or any other field;
- treating WRF end-state/oracle output, reference-copy files, or output without
  applied pressure refresh as a validation-gate pass;
- copying, blending, nudging toward, bias-correcting against, or otherwise
  patching metrics with WRF cycle-end reference output;
- using reference output to move the TC center, adjust the minimum SLP, or tune
  pressure ranges;
- silently replacing normal skeleton or integrator fields with closure fields;
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
