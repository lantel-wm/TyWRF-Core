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
checks.

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
