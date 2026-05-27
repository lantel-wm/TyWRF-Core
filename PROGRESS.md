# TyWRF-Core Progress

## 2026-05-27

Current phase: Phase 2-4 baseline integration after Phase 1 Bootstrap

Completed work:

- Read `AGENTS.md` and `PLAN.md`.
- Confirmed the local project had no usable implementation files yet.
- Checked the KROSA reference data paths and sampled NetCDF dimensions.
- Added the first CMake project skeleton and C++ smoke test.
- Added documentation skeleton and initial Python validation tool entry points.
- Initialized Git repository and set the initial branch to `main`.
- Switched execution model to the PLAN.md main-agent/sub-agent workflow at the
  user's request.
- Spawned disjoint sub-agent workstreams:
  - Agent A: WRF NetCDF I/O schema and reduced core `wrfout` writer.
  - Agent B: CUDA-ready `FieldView`/`Grid`/`State` layout baseline.
  - Agent C: dynamics loop skeleton and zero/identity call-point sequencing.
  - Agent F: validation/RMSE/regression Python tooling.
- Agent C returned a completed dynamics skeleton subtask for main-agent review.
- Agent F returned a completed validation tooling subtask for main-agent review.
- Agent A returned a completed WRF NetCDF I/O/schema subtask for main-agent
  review.
- Agent B returned a completed CUDA-ready layout subtask for main-agent review.
- Main-agent integration review completed for the first A/B/C/F subagent round.
- Unified validation passed for the integrated Phase 2-4 baseline.
- Created verified Git milestone commit `48a8415`:
  `Bootstrap TyWRF-Core baseline`.
- Updated `AGENTS.md` with the persistent collaboration constraints requested
  by the user: proactive subagent dispatch each round, near-real-time
  `PROGRESS.md` updates, ongoing progress reports, and verified milestone Git
  commits.

Reference observations:

- `wrfout_d01_2025-07-26_00:00:00` mass grid is `265 x 429`, with `59`
  mass levels and `60` staggered vertical levels.
- `wrfout_d02_2025-07-26_00:00:00` mass grid is `210 x 210`, with `59`
  mass levels and `60` staggered vertical levels.
- `wrfbdy_d01` has `28` times and `bdy_width = 5`.
- `wrffdda_d01` has `28` times and includes U/V/T/PH/MU old/new nudging
  fields plus Q nudging fields named `Q_NDG_OLD` and `Q_NDG_NEW`.

Changed files:

- `CMakeLists.txt`
- `.gitignore`
- `AGENTS.md`
- `README.md`
- `PROGRESS.md`
- `docs/architecture.md`
- `docs/cuda_ready_layout.md`
- `docs/validation_plan.md`
- `docs/wrf_compatibility.md`
- `include/tywrf/version.hpp`
- `include/tywrf/field_view.hpp`
- `include/tywrf/grid.hpp`
- `include/tywrf/state.hpp`
- `include/tywrf/dynamics/dynamics_loop.hpp`
- `src/core/version.cpp`
- `src/dynamics/dynamics_loop.cpp`
- `include/tywrf/io/netcdf_schema.hpp`
- `src/io/netcdf_schema.cpp`
- `tests/unit/dynamics_loop_test.cpp`
- `tests/unit/field_layout_test.cpp`
- `tests/unit/netcdf_schema_test.cpp`
- `tests/unit/smoke_test.cpp`
- `tests/unit/test_compare_wrfout.py`
- `tests/unit/test_extract_reference.py`
- `tests/unit/test_run_6h_cycle_test.py`
- `tests/unit/test_write_core_wrfout.py`
- `tools/compare_wrfout.py`
- `tools/extract_reference.py`
- `tools/__init__.py`
- `tools/run_6h_cycle_test.py`
- `tools/write_core_wrfout.py`
- directory placeholders under `src/`, `bindings/`, `tests/`, and
  `third_party/`

Sub-agent status:

- Agent A: completed, pending main-agent integration review.
  - Reported changed files: `src/io/netcdf_schema.cpp`,
    `tests/unit/netcdf_schema_test.cpp`, `tools/write_core_wrfout.py`,
    `tests/unit/test_write_core_wrfout.py`, `docs/wrf_compatibility.md`.
  - Reported tests: `cmake --build build -j`,
    `ctest --test-dir build --output-on-failure`,
    `./build/tywrf_netcdf_schema_test`,
    `UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python -m pytest tests/unit/test_write_core_wrfout.py`,
    and full `tests/unit` pytest.
  - Reported risks: schema is intentionally KROSA-specific; Python writer is a
    compatibility copier, not the future integrator writer; reduced core
    `wrfout` files can still be large.
- Agent B: completed, pending main-agent integration review.
  - Reported changed files: `include/tywrf/field_view.hpp`,
    `include/tywrf/grid.hpp`, `include/tywrf/state.hpp`,
    `tests/unit/field_layout_test.cpp`, `docs/cuda_ready_layout.md`,
    `CMakeLists.txt`.
  - Reported tests: `cmake -S . -B build`, `cmake --build build -j`,
    and `ctest --test-dir build --output-on-failure`.
  - Reported risks: layout is a baseline contract only; no dynamics/physics
    kernels consume `StateView` yet; live NetCDF loading has not been wired to
    field views; halo policy is uniform for now.
- Agent C: completed, pending main-agent integration review.
  - Reported changed files: `include/tywrf/dynamics/dynamics_loop.hpp`,
    `src/dynamics/dynamics_loop.cpp`, `tests/unit/dynamics_loop_test.cpp`,
    `docs/architecture.md`, `CMakeLists.txt`.
  - Reported tests: `cmake --build build -j`,
    `ctest --test-dir build -R tywrf_dynamics_loop_test --output-on-failure`,
    and full `ctest --test-dir build --output-on-failure`.
  - Reported risks: zero/identity dynamics only; moving nest relocation not
    implemented; call order not yet validated against WRF driver internals.
- Agent F: completed, pending main-agent integration review.
  - Reported changed files: `tools/compare_wrfout.py`,
    `tools/extract_reference.py`, `tools/run_6h_cycle_test.py`,
    `tests/unit/test_compare_wrfout.py`,
    `tests/unit/test_extract_reference.py`,
    `tests/unit/test_run_6h_cycle_test.py`, `docs/validation_plan.md`.
  - Reported tests:
    `UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python -m pytest tests/unit/test_compare_wrfout.py tests/unit/test_extract_reference.py tests/unit/test_run_6h_cycle_test.py`.
  - Reported result: 12 Python tests passed.
  - Reported risks: TC diagnostics still pending; 6 h cycle tool remains
    dry-run only; tests use synthetic files rather than the full KROSA dataset.

Git status policy:

- Commit after each main-agent-reviewed and verified integration milestone.
- Do not commit sub-agent intermediate work before unified review and tests.
- Latest verified milestone commit: `48a8415`.

Test results:

- `cmake -S . -B build -DTYWRF_ENABLE_OPENMP=ON` passed.
- `cmake --build build -j` passed.
- `ctest --test-dir build --output-on-failure` passed: 4/4 C++ tests:
  smoke, NetCDF schema, dynamics loop, and field layout.
- `UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python -m pytest`
  passed: 16/16 Python tests.
- `./build/tywrf_netcdf_schema_test` validated the actual KROSA
  `wrfinput_d01`, `wrfbdy_d01`, and `wrffdda_d01` schemas.
- `tools/extract_reference.py` successfully summarized the d01 initial
  reference `wrfout`.
- `tools/extract_reference.py` successfully summarized both d01 and d02 initial
  reference `wrfout` files with no missing core variables.
- `tools/write_core_wrfout.py` successfully wrote a reduced d02 core
  `wrfout` from the real KROSA reference file with all 29 minimum variables.
- `tools/compare_wrfout.py` successfully compared the real d02 reference
  `wrfout` against the reduced core copy for the default 22 comparison
  variables: all zero RMSE.
- `tools/run_6h_cycle_test.py --dry-run` resolved d01/d02 6 h reference start
  and end files for `2025-07-26_00:00:00` to `2025-07-26_06:00:00`.

Blockers:

- None for the current integrated baseline.

Next step:

- Proactively dispatch the next subagent round from the remaining plan after
  this integration checkpoint.
