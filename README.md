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
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/compare_wrfout.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/write_core_wrfout.py --help
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/run_6h_cycle_test.py --help
```

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
