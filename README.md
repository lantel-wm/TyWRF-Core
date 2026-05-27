# TyWRF-Core

TyWRF-Core is a CPU-first, CUDA-ready WRF-compatible typhoon integrator for the
current PGWRF/KROSA configuration subset.

The v1 target is deliberately narrow:

- reuse WPS and `real.exe`;
- read the current WRF standard inputs;
- write WRF-compatible core `wrfout` fields;
- validate progressive 10 min gates against the existing WRF reference before
  longer cycle tests;
- defer best-track nudging, arbitrary namelists, full WRF replacement behavior,
  and direct CUDA implementation.

## Current Status

Phase 1 Bootstrap is active. The repository now contains the CMake skeleton,
documentation skeleton, a tiny C++ core library, CTest smoke test, and initial
Python validation tool entry points.

## Build

```bash
cmake -S . -B build -DTYWRF_ENABLE_OPENMP=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## C++ Tools

`tywrf_skeleton_cycle` is the current C++ executable shell for a d02 10 min
candidate. It reads a cycle-start WRF state, copies `XLAT/XLONG/HGT` from the
static template, writes selected core fields, and marks the result as
`skeleton=true`, `not_physical=true`, and `integrator_output=false`. It enforces
d02 `DX/DY=2000m`; it does not run dynamics or physics.

For d02 moving-nest persistence smoke tests, pass the cycle-start wrfout as the
static `--template` and set output `Times` with `--times`. This keeps the 00 h
state paired with the 00 h nest coordinates while labeling the candidate at the
cycle end. Only use cycle-end template coordinates when the state has already
been remapped to the end nest pose by a real integrator.

```bash
./build/tywrf_skeleton_cycle \
  --state /path/to/wrfout_d02_2025-07-26_00:00:00 \
  --template /path/to/wrfout_d02_2025-07-26_00:00:00 \
  --output /path/to/candidate/wrfout_d02_2025-07-26_00:10:00 \
  --cycle-start 2025-07-26_00:00:00 \
  --cycle-end 2025-07-26_00:10:00 \
  --times 2025-07-26_00:10:00 \
  --pretty
```

The output NetCDF metadata and JSON report identify the state source, static
coordinate source, `Times` source, and whether static coordinates came from the
same source file as the state. When `--cycle-start` and `--cycle-end` describe
a 10 min candidate, the metadata/report also records `--interval-minutes 10` as
the gate interval option. The default strict gate still rejects this output
because it is marked `TYWRF_VALIDATION_GATE_ONLY=true` and
`TYWRF_INTEGRATOR_OUTPUT=false`; that failure is expected for the current
skeleton smoke and must not be reported as a physical pass.

An opt-in diagnostic remap path is available for wiring the moving-nest overlap
remap into this C++ skeleton:

```bash
./build/tywrf_skeleton_cycle \
  --state /path/to/wrfout_d02_2025-07-26_00:00:00 \
  --template /path/to/wrfout_d02_2025-07-26_00:00:00 \
  --output /path/to/diagnostic/remap_overlap_wrfout_d02 \
  --cycle-start 2025-07-26_00:00:00 \
  --cycle-end 2025-07-26_00:10:00 \
  --times 2025-07-26_00:10:00 \
  --diagnostic-remap-overlap \
  --from-parent-start 114,96 \
  --to-parent-start 115,96 \
  --pretty
```

This path only copies d02 active cells covered by the old/new moving-nest
overlap. Newly exposed cells are left as NaN until a future parent-fill step is
implemented. Its NetCDF metadata and JSON report mark
`TYWRF_DIAGNOSTIC_REMAP_OVERLAP=true`, `TYWRF_NEEDS_PARENT_FILL`, copied
field/point counts, from/to parent starts, and `TYWRF_NOT_PHYSICAL=true`. It is
a diagnostic candidate only, not a physical integrator result and not a
validation-gate candidate.

## Python Tools

Use the project-local uv environment:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv sync --extra dev
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python -m pytest
```

Useful entry points:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/extract_reference.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/audit_reference_cycles.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/compare_wrfout.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/cycle_gate.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/baseline_candidate.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/skeleton_cycle_driver.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/run_skeleton_cycle.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/write_core_wrfout.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/run_6h_cycle_test.py --help
```

`tools/baseline_candidate.py` and `tools/run_6h_cycle_test.py --mode` generate
explicit baseline candidate files only. `persistence` copies a cycle-start
`wrfout` to the cycle-end filename and rewrites `Times`; it is expected to fail
strict thresholds and is used to quantify the current gap. `reference-copy`
copies the cycle-end reference file and marks JSON/file metadata as
`reference_copy`; it is only for checking the validation gate itself and is not
a TyWRF-Core integrator result. Its metadata marks `minimum_slp_gate_ready`
only when a real SLP/MSLP diagnostic was copied.

`tools/skeleton_cycle_driver.py` is the current d02 skeleton candidate
entry point. It copies the cycle-start reference state to the cycle-end
candidate file in `persistence` or `identity` mode, then marks both JSON and
NetCDF metadata with `skeleton=true`, `not_physical=true`, and
`integrator_output=false`; it also marks `validation_gate_only=true`. The report
includes the suggested `cycle_gate.py` command to run next. It defaults to a
6 h end time for existing cycle tests, but also accepts explicit `--end` times
or `--minutes 10` for progressive 10 min validation wiring; this output is
wiring infrastructure, not a physical TyWRF-Core result.

`tools/run_skeleton_cycle.py` is the dual-domain wrapper for the same skeleton
path. It generates d01 and d02 candidate files, keeps d02 on the existing
`DX/DY=2000m` check, does not require d01 to be 2 km, and writes a
machine-readable report with per-domain source, candidate, reference-end,
copied/missing variables, cycle minutes, and suggested `cycle_gate.py` commands.
Its JSON and NetCDF metadata also mark `skeleton=true`, `not_physical=true`, and
`integrator_output=false`; this is a shell for later C++ executable wiring, not
an integrator result.

The current acceptance gate is `tools/cycle_gate.py`. It defaults to d02 and
checks each cycle-end `wrfout` for `U,V,T,PH,MU,P,QVAPOR` normalized RMSE <=
`0.05`, storm-center error <= `20 km`, minimum SLP error <= `5 hPa`, and
Vmax10m error <= `5 m s-1`. Missing candidate files or missing real SLP
diagnostics fail the gate. For the KROSA progressive validation smoke, run the
first endpoint with `--start 2025-07-26_00:00:00 --end 2025-07-26_00:10:00
--interval-minutes 10`; a C++ skeleton/persistence candidate is expected to
fail strict mode at `00:10`.

TC diagnostics report `minimum_slp_hpa` only from a real sea-level pressure
field such as `SLP`, `MSLP`, `AFWA_MSLP`, `PMSL`, or `PRMSL`. When no accepted
SLP field is present, `minimum_slp_status` is `not_available`; the minimum
`PSFC` value is kept only as labeled proxy metadata and is not treated as
meeting the minimum SLP validation objective.

## Reference Case

- Storm: `2025WP12 KROSA`
- Initial time: `2025-07-26 00 UTC`
- Primary validation length: `1 h`
- Primary validation output interval: `10 min`
- Full reference length: `168 h`
- Full reference output interval: `6 h`
- d01: `10 km`
- d02: `2 km`, moving, two-way nested
- Vertical levels: `60` full eta levels, `59` mass levels

See `AGENTS.md`, `PLAN.md`, and `docs/` for project scope and implementation
rules.
