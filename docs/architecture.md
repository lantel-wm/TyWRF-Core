# Architecture

TyWRF-Core is structured around a small C++20 integrator core with explicit
boundaries between state layout, dynamics kernels, physics staging, nesting,
NetCDF I/O, and validation tools.

## v1 Scope

v1 targets only the current PGWRF/KROSA setup described in `AGENTS.md`:

- d01 at `10 km`;
- d02 at `2 km`;
- moving two-way nest;
- 60 full eta levels and 59 mass levels;
- lateral boundary input from `wrfbdy_d01`;
- spectral nudging input from `wrffdda_d01`;
- WRF physics reached through staged Fortran wrappers.

Out-of-scope items include arbitrary namelists, full WRF output coverage,
best-track nudging, and direct CUDA kernels.

## Module Boundaries

- `include/tywrf/`: public C++ headers and POD view types.
- `src/core/`: build/version utilities and eventually shared integrator
  orchestration.
- `src/dynamics/`: time loop and dynamics kernel dispatch.
- `src/io/`: WRF NetCDF readers and writers.
- `src/nest/`: moving nest, interpolation, and feedback interfaces.
- `src/physics_bridge/`: C++ side of WRF physics staging.
- `bindings/wrf_physics/`: Fortran/C ABI wrappers around WRF physics.
- `tools/`: Python reference extraction, RMSE comparison, and 6 h cycle
  helpers.

Phase 1 intentionally keeps these modules skeletal. Dynamics and physics
implementation begins only after I/O and data layout are pinned down.

## Phase 4 Dynamics Skeleton

The first dynamics implementation is an orchestration skeleton, not a numerical
core. `tywrf::dynamics::DynamicsLoopRunner` owns only domain and time-step
descriptors plus call-point sequencing for the fixed KROSA v1 configuration:

- d01 parent: `10 km`, `40 s`, lateral boundary update enabled, spectral
  nudging enabled;
- d02 moving nest: `2 km`, `8 s`;
- `parent_time_step_ratio = 5`;
- one 6 h segment is `540` d01 steps and `2700` d02 substeps.

Each d01 step emits call points in this order:

1. d01 lateral boundary update;
2. d01 spectral nudging;
3. d01 zero dynamics tendency;
4. d01 physics call point;
5. five d02 subcycles, each with parent-child interpolation, zero dynamics
   tendency, and physics call point;
6. d02-to-d01 feedback;
7. d01 and d02 history output when the step lands on the 6 h history interval.

All dynamics tendencies are currently zero/identity placeholders. The event sink
is callback-based and carries no NetCDF, logging, or state-layout dependency, so
Agent B's field/state interfaces and later physics/nesting work can attach at
these call points without changing the tested time sequencing.
