# TyWRF-Core Agent Context

This file is the persistent context for Codex agents working on **TyWRF-Core**.
It should be read before planning or editing the project.

## Project Identity

- Project name: `TyWRF-Core`
- Target project root: `/home/zzy/Projects/TyWRF-Core`
- The Codex CLI that will execute this project is expected to run directly in the target project environment.
- The current repository `/home/zzy/Projects/tc_sim` contains reference scripts, diagnostics, PGWRF experiments, and plotting utilities.

## Hard Safety Rules

- Do not reduce d02 horizontal resolution; d02 must remain `2 km` for this project.
- Do not start by implementing best-track nudging. It is a later feature, after the compatible integrator and validation framework exist.
- Do not treat the project as a full WRF replacement in v1. v1 targets the current PGWRF/KROSA configuration subset.

## Role Model

The active Codex instance should act as the **main agent**.

Responsibilities:

- Own architecture decisions.
- Split work into independent subtasks.
- Use subagents when available.
- Keep subagent write scopes disjoint.
- When spawning subagents, do not enable fast mode and do not select low
  reasoning; use the inherited/default model and at least normal reasoning
  unless the user explicitly changes this constraint.
- Use a `600` second wait window for subagent rounds when waiting is needed.
- Close completed or no-longer-needed subagents promptly after integrating
  their results.
- Review and integrate all subagent outputs.
- Maintain `PROGRESS.md` with completed work, in-progress work, blockers, and next steps.
- Report only the integrated state to the user.
- Do not directly take over implementation work when subagents are available
  and the work can be split cleanly.
- After each subagent round completes, proactively derive the next work items
  from the full project goal, split them into disjoint scopes, and dispatch the
  next set of subagents without waiting for the user to restate the task.
- Keep reporting current integrated progress to the user while work is ongoing.
- Keep `PROGRESS.md` updated in near real time when subagents start, finish,
  report risks, or when integration/test/Git status changes.
- Write `PROGRESS.md` in Chinese. It is a local progress log, should remain in
  `.gitignore`, and should not be committed.
- Keep `PROGRESS.md` in ascending round/time order, with the newest integrated
  state always appended at the bottom. If reversed order, insertion into older
  rounds, or old rounds drifting to the bottom is detected, fix the ordering
  before continuing development.
- Commit Git changes at verified integration milestones. Do not commit raw
  subagent intermediate work before main-agent review and unified validation.
- Prefer small, coherent commits that match completed phases or reviewed
  integration slices.
- After each successful commit, push the committed branch to `origin`.
- If `origin` is missing or push fails, record the push blocker in
  `PROGRESS.md` and report it to the user immediately.

If subagents are unavailable, execute the same workstreams sequentially.

Suggested subagent workstreams:

- Agent A: WRF NetCDF I/O schema and reference extraction.
- Agent B: CUDA-ready field layout, halo, indexer, and kernel views.
- Agent C: dynamics loop skeleton and CPU/OpenMP kernel dispatch.
- Agent D: WRF Fortran physics bridge feasibility and staging buffers.
- Agent E: moving nest, parent-child interpolation, feedback, and spectral nudging interfaces.
- Agent F: regression tests, RMSE metrics, and TC diagnostics.

## Current WRF Baseline

Primary reference case:

- Storm: `2025WP12 KROSA`
- Initial time: `2025-07-26 00 UTC`
- Primary validation length: `1 h`
- Primary validation output frequency: `10 min`
- Full reference length: `168 h`
- Full reference output frequency: `6 h`
- Model: PGWRF / WRF v4.6.1
- Forcing: GFS analysis as lateral boundary and spectral nudging source.

Important reference paths:

- PGWRF root: `/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF`
- Primary 1 h / 10 min WRF validation output directory:
  `/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF_1h_10min_20260527_172838`
- Full 168 h / 6 h WRF reference output directory:
  `/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF`
- Logs:
  `/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/logs`
- Existing diagnostics:
  `/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/diagnostics`

Reference input files:

- `namelist.input`
- `wrfinput_d01`
- `wrfbdy_d01`
- `wrffdda_d01`

Reference output files:

- Primary validation files: `wrfout_d01_YYYY-MM-DD_HH:MM:00`, `wrfout_d02_YYYY-MM-DD_HH:MM:00`
- Required validation times: `00:00`, `00:10`, `00:20`, `00:30`, `00:40`, `00:50`, `01:00`
- Full reference files: `wrfout_d01_YYYY-MM-DD_HH:00:00`, `wrfout_d02_YYYY-MM-DD_HH:00:00`

Current namelist highlights:

- Domains: d01 `10 km`, d02 `2 km`, two-way nested moving nest.
- Grid sizes in namelist: d01 `266 x 430`, d02 `211 x 211`.
- NetCDF mass-grid dimensions observed in wrfout: d01 `265 x 429`, d02 `210 x 210`.
- Vertical levels: `60` full eta levels, `59` mass levels.
- Model top: `50 hPa`.
- Time step: d01 `40 s`, d02 `8 s`, `parent_time_step_ratio = 5`.
- Boundary interval: `21600 s`.
- History interval: `360 min`.
- MPI baseline: `WRF_NP=44`, `nproc_x=4`, `nproc_y=11`.

Current physics:

- Microphysics: Thompson, `mp_physics = 8`.
- Cumulus: Kain-Fritsch on d01, explicit convection on d02, `cu_physics = 1, 0`.
- Longwave/shortwave radiation: RRTMG, `ra_lw_physics = 4`, `ra_sw_physics = 4`.
- PBL: YSU, `bl_pbl_physics = 1`.
- Surface layer: MM5, `sf_sfclay_physics = 1`.
- Land surface: Noah LSM, `sf_surface_physics = 2`.
- Slab ocean and TC flux are enabled: `sf_ocean_physics = 1`, `isftcflx = 2`.

Spectral nudging:

- `grid_fdda = 2`
- `gfdda_inname = "wrffdda_d<domain>"`
- `gfdda_interval_m = 360`
- `guv = 0.0003`
- `xwavenum = 2`
- `ywavenum = 4`
- `wrffdda_d01` contains U/V/T/Q/PH/MU old/new nudging fields.

Moving nest:

- `vortex_interval = 15, 15`
- `max_vortex_speed = 30, 30`
- `corral_dist = 10, 30`
- `track_level = 70000`
- `time_to_move = 0., 0.`

Timing baseline from the completed 168 h reference run:

- d01 timing: median about `5.30 s/step`, mean about `6.41 s/step`.
- d02 timing: median about `0.61 s/step`, mean about `0.70 s/step`.
- Approximate combined cost per d01 step: `d01 + 5*d02 ~= 8.35 s`.
- d01 is still a major bottleneck even though d02 is high resolution.

## TyWRF-Core Goal

Build a CPU-first, CUDA-ready WRF-compatible typhoon integrator.

First-stage goal:

- Support the current PGWRF/KROSA configuration subset.
- Reuse WPS and `real.exe`; do not reimplement geogrid/ungrib/metgrid/real in v1.
- Read WRF standard inputs.
- Write WRF-compatible core `wrfout`.
- Achieve small RMSE against WRF through progressive 10 min validation gates.
- Only after compatibility is established, optimize CPU performance and migrate hot kernels to CUDA.

This project does not require bitwise identity with WRF.

## Current Development Direction

- Starting with the next development round after the D94 mass/pressure blocker
  closure, prioritize the system integrator instead of additional isolated
  sensitivity probes.
- Split the work into four coordinated tracks:
  - integrator orchestration and domain/time-step scheduling;
  - forcing ingestion, boundary forcing, and spectral nudging staging;
  - WRF physics bridge staging and call-order scaffolding;
  - moving nested domain, parent-child interpolation, and feedback.
- Use WRF v4.6.1/PGWRF behavior as the compatibility reference, but do not copy
  monolithic Fortran control flow directly into hot paths.
- Decompose modules around CUDA-ready kernel boundaries: POD views, explicit
  staging buffers, explicit halo/boundary operations, no NetCDF/logging/hidden
  allocation inside kernels, and CPU/OpenMP implementations that can later map
  to CUDA kernels.
- Keep early integrator slices default-off or non-gate until their metadata,
  non-oracle data dependencies, and `00:10` validation behavior are verified.

## v1 Functional Boundary

Required v1 inputs:

- `namelist.input`
- `wrfinput_d01`
- `wrfbdy_d01`
- `wrffdda_d01`

Required v1 outputs:

- WRF-compatible NetCDF files for d01 and d02 at 10 min validation checkpoints.
- Later full-run output can use 6 h intervals, but v1 validation must support 10 min output.
- The files do not need all 220 WRF variables.
- They must contain the core fields needed for validation and existing diagnostics.

Minimum output variables:

- Coordinates/static/time: `Times`, `XLAT`, `XLONG`, `HGT`
- Dynamics/core state: `U`, `V`, `W`, `PH`, `PHB`, `T`, `MU`, `MUB`, `P`, `PB`
- Moisture/hydrometeors: `QVAPOR`, `QCLOUD`, `QRAIN`, `QICE`, `QSNOW`, `QGRAUP`, `QNICE`, `QNRAIN`
- Near-surface diagnostics: `PSFC`, `U10`, `V10`, `T2`, `Q2`
- Precipitation: `RAINC`, `RAINNC`

Required v1 model features:

- d01 lateral boundary update from `wrfbdy_d01`
- d01 spectral nudging from `wrffdda_d01`
- d02 moving nested domain
- parent-child interpolation
- two-way nesting feedback
- current physics timing and call order as closely as practical

Out of scope for v1:

- Arbitrary WRF namelists
- Full WRF 220-variable output
- Best-track nudging
- Alternative physics suites
- Direct CUDA implementation
- Full WRF replacement

## Technical Stack

Preferred stack:

- C++20 for the main integrator framework.
- CMake for builds.
- OpenMP for CPU parallel reference execution.
- NetCDF-C / NetCDF-Fortran for WRF-compatible I/O.
- Fortran wrappers around WRF physics via `ISO_C_BINDING`.
- Python/pytest tools for reference extraction, RMSE comparison, and typhoon diagnostics.

The CPU implementation must be written as an engineering reference for later CUDA migration, not as a throwaway prototype.

## Python Environment

Use the project-local `uv` environment for Python tools and tests.

Environment policy:

- Keep the virtual environment in `/home/zzy/Projects/TyWRF-Core/.venv`.
- Keep uv caches inside the project with `UV_CACHE_DIR=.uv-cache`.
- Keep uv-managed Python installs inside the project with `UV_PYTHON_INSTALL_DIR=.uv-python`.
- Use Python `3.11` for this project; avoid using the system `python3` directly.
- Runtime and validation dependencies are declared in `pyproject.toml`.
- Development/test dependencies are installed with the `dev` extra.

Create or refresh the environment from the project root:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv venv --python 3.11
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv sync --extra dev
```

Run Python commands through uv from the project root:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python -m pytest
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/compare_wrfout.py --help
```

For quick interactive checks:

```bash
source .venv/bin/activate
python -c "import netCDF4, numpy, xarray; print(netCDF4.__version__, numpy.__version__, xarray.__version__)"
```

## CUDA-Ready Design Rules

Core state layout:

- Use structure-of-arrays.
- Each physical variable owns a contiguous buffer.
- No nested `std::vector` in hot paths.
- No arrays of structs for grid fields.
- Use explicit halo cells.
- Keep `i` as the contiguous dimension.

Canonical flattened 3D indexing:

```text
idx = ((j * nz) + k) * nx + i
```

Kernel API rules:

- Kernels receive POD views only: pointer, shape, stride, halo.
- No NetCDF I/O inside kernels.
- No logging inside kernels.
- No virtual dispatch inside kernels.
- No hidden allocation inside kernels.
- CPU scalar, CPU OpenMP, and future CUDA kernels should share the same view types.

Loop rules:

- Make `i` the innermost loop for contiguous access.
- Prefer tiled loops over large monolithic loops.
- Keep boundary/halo updates explicit.
- Keep staging buffers explicit when calling WRF Fortran physics.

## Validation Standard

Validation is based on progressive 10 min gates using the 1 h WRF reference run. Do not start with a 6 h validation gate; it is too slow for early debugging.

Primary reference directory:

```text
/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF_1h_10min_20260527_172838
```

Procedure:

1. Start TyWRF-Core from the same WRF reference initial state at `2025-07-26_00:00:00`.
2. Integrate to `00:10`.
3. Compare TyWRF-Core output to WRF reference output at `00:10`.
4. If the `00:10` RMSE gate fails, stop immediately and fix the first failing field or kernel. Do not continue to `00:20`.
5. If `00:10` passes, continue to `00:20` and compare again.
6. Continue the same gated sequence through `00:30`, `00:40`, `00:50`, and `01:00`.
7. Only after all 10 min gates pass should longer 3 h or 6 h cycle tests be attempted.

Default field thresholds:

- `U`, `V`, `T`, `PH`, `MU`, `P`, `QVAPOR`: normalized RMSE <= `5%`
- `W` and hydrometeors: normalized RMSE <= `10%`
- `PSFC`, `U10`, `V10`, `T2`, `Q2`: normalized RMSE <= `10%`

Default typhoon diagnostics:

- d02 storm center error <= `20 km`
- MSLP error <= `5 hPa`
- Vmax error <= `5 m s-1`

Validation tools should report:

- RMSE and normalized RMSE per variable/domain/time.
- Max absolute error per variable.
- TC center, MSLP, Vmax, and accumulated rainfall summary.
- Timing per kernel and per 10 min segment.
- The first failed gate and first failed variable.

## Suggested Implementation Phases

Phase 1: Bootstrap

- Create project skeleton.
- Add `README.md`, `PROGRESS.md`, `docs/`, `src/`, `include/`, `tests/`, `tools/`, `bindings/`.
- Add CMake skeleton and a passing smoke test.

Phase 2: WRF I/O Baseline

- Read `namelist.input`, `wrfinput_d01`, `wrfbdy_d01`, `wrffdda_d01`.
- Write WRF-compatible core NetCDF output.
- Confirm existing Python diagnostics can read the output.

Phase 3: Data Layout Freeze

- Implement `FieldView`, `Grid`, `State`, halo/indexer abstractions.
- Add tests for indexing, halo shape, contiguous access, and OpenMP-safe iteration.

Phase 4: Dynamics Skeleton

- Implement d01/d02 time loop, d02 subcycling, boundary/nudging/nest call points.
- Use zero or identity tendencies first to validate time sequencing and I/O.

Phase 5: Physics Bridge

- Build WRF Fortran physics wrapper.
- Define staging buffers.
- Call current physics suite through wrappers.

Phase 6: Progressive 10 min Validation

- Run the first 10 min segment.
- Compare to WRF reference at `00:10`.
- Stop and fix immediately if the 10 min gate fails.
- Advance to 20, 30, 40, 50, and 60 min only after the previous gate passes.
- Record RMSE, TC diagnostics, and timing.

Phase 7: Optimization Preparation

- Identify hot kernels.
- Improve CPU cache behavior and OpenMP schedule.
- Prepare CUDA migration plan for selected kernels.

## Reporting Rules

Keep `PROGRESS.md` current.

Append each new integrated status update to the bottom of `PROGRESS.md`; do not
insert current updates above older rounds.

Each progress update should include:

- Current phase
- Completed work
- Changed files
- Test results
- Blockers
- Next step

Do not report subagent raw outputs directly to the user. The main agent should synthesize and report integrated status only.

## Useful Commands

Start from project root:

```bash
cd /home/zzy/Projects/TyWRF-Core
```

Current PGWRF reference root:

```bash
cd /home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF
```

Primary 1 h / 10 min WRF validation output directory:

```bash
cd /home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF_1h_10min_20260527_172838
```

Full 168 h / 6 h WRF reference output directory:

```bash
cd /home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF
```

Use the existing Python environment for NetCDF inspection when useful:

```bash
/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/.venv-wrfdiag/bin/python
```

## Default First Prompt for Codex CLI

```text
You are the main agent for TyWRF-Core. Read AGENTS.md first. Then begin Phase 1 Bootstrap only.

Important:
- Keep v1 scoped to the current PGWRF/KROSA configuration.
- Use CPU-first, CUDA-ready design.
- Do not implement dynamics or physics yet in Phase 1.
- Create project skeleton, CMake smoke test, README, docs skeleton, and PROGRESS.md.
```
