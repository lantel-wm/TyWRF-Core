# TyWRF-Core

TyWRF-Core is a CPU-first, CUDA-ready WRF-compatible typhoon integrator for the
current PGWRF/KROSA configuration subset.

The v1 target is deliberately narrow:

- reuse WPS and `real.exe`;
- read the current WRF standard inputs;
- write WRF-compatible core `wrfout` fields;
- validate 6 h cycle tests against the existing WRF reference;
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

`tools/skeleton_cycle_driver.py` is the current d02 6 h skeleton candidate
entry point. It copies the cycle-start reference state to the cycle-end
candidate file in `persistence` or `identity` mode, then marks both JSON and
NetCDF metadata with `skeleton=true`, `not_physical=true`, and
`integrator_output=false`; it also marks `validation_gate_only=true`. The report
includes the suggested `cycle_gate.py` command to run next; this output is
wiring infrastructure, not a physical TyWRF-Core result.

`tools/run_skeleton_cycle.py` is the dual-domain wrapper for the same skeleton
path. It generates d01 and d02 6 h candidate files, keeps d02 on the existing
`DX/DY=2000m` check, does not require d01 to be 2 km, and writes a
machine-readable report with per-domain source, candidate, reference-end,
copied/missing variables, and suggested `cycle_gate.py` commands. Its JSON and
NetCDF metadata also mark `skeleton=true`, `not_physical=true`, and
`integrator_output=false`; this is a shell for later C++ executable wiring, not
an integrator result.

The current acceptance gate is `tools/cycle_gate.py`. It defaults to d02 and
checks each 6 h cycle-end `wrfout` for `U,V,T,PH,MU,P,QVAPOR` normalized RMSE <=
`0.05`, storm-center error <= `20 km`, minimum SLP error <= `5 hPa`, and
Vmax10m error <= `5 m s-1`. Missing candidate files or missing real SLP
diagnostics fail the gate.

TC diagnostics report `minimum_slp_hpa` only from a real sea-level pressure
field such as `SLP`, `MSLP`, `AFWA_MSLP`, `PMSL`, or `PRMSL`. When no accepted
SLP field is present, `minimum_slp_status` is `not_available`; the minimum
`PSFC` value is kept only as labeled proxy metadata and is not treated as
meeting the minimum SLP validation objective.

## Reference Case

- Storm: `2025WP12 KROSA`
- Initial time: `2025-07-26 00 UTC`
- Forecast length: `168 h`
- Output interval: `6 h`
- d01: `10 km`
- d02: `2 km`, moving, two-way nested
- Vertical levels: `60` full eta levels, `59` mass levels

See `AGENTS.md`, `PLAN.md`, and `docs/` for project scope and implementation
rules.
