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

## Hard Prohibitions

The following schemes are forbidden:

- using a WRF reference end-state delta, such as
  `reference_end - reference_start`, to update candidate `MU`, `P`, `PSFC`,
  `SLP`, or any other field;
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
