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

## Nest And Nudging Interfaces

`tywrf::nest` defines the v1 interface baseline for KROSA moving-nest and
spectral-nudging call points. The interface is intentionally contract-first:

- `ParentChildDescriptor` pins d01/d02 spacing, WRF namelist extents, mass-grid
  extents, `parent_grid_ratio = 5`, and `parent_time_step_ratio = 5`.
- `ParentChildPosition` stores the WRF namelist parent-start indices for d02
  (`i_parent_start = 114`, `j_parent_start = 96`) with an explicit index base.
- `MovingNestConfig` records the KROSA vortex interval, maximum vortex speed,
  corral distance, 700 hPa tracking level, and `time_to_move` values for each
  domain.
- `SpectralNudgingConfig` records `grid_fdda = 2`, `wrffdda_d<domain>`, the
  6 h FDDA interval, `guv = 0.0003`, `xwavenum = 2`, `ywavenum = 4`, and the
  required old/new wrffdda field names for U, V, T, Q, PH, and MU.

Numerical parent-child interpolation and two-way feedback are still explicit
`not_implemented` stubs. They validate descriptors, position bounds, timing,
and child-substep contracts before reporting that status, so the dynamics loop
can wire to stable call contracts without pretending the numerics exist.

## Phase 4 Dynamics Skeleton

The first dynamics implementation is an orchestration skeleton, not a numerical
core. `tywrf::dynamics::DynamicsLoopRunner` owns only domain and time-step
descriptors for the fixed KROSA v1 configuration and consumes the
`CycleSchedule` contract for boundary, nudging, moving-nest, interpolation,
feedback, and history call-point sequencing:

- d01 parent: `10 km`, `40 s`, lateral boundary update enabled, spectral
  nudging enabled;
- d02 moving nest: `2 km`, `8 s`;
- `parent_time_step_ratio = 5`;
- one 6 h segment is `540` d01 steps and `2700` d02 substeps.

The runner traverses every schedule call, then inserts current no-op numerical
placeholders at the compatible call points. Segment-start and segment-end
boundary/nudging input refreshes, inclusive moving-nest position updates, and
history output are visible in the callback event stream. Each d01 step emits
these effective call points:

1. d01 lateral boundary update;
2. d01 spectral nudging;
3. d01 zero dynamics tendency;
4. d01 physics call point;
5. five d02 subcycles, each with parent-child interpolation, zero dynamics
   tendency, and physics call point;
6. d02-to-d01 feedback;
7. d01 and d02 history output when the step lands on the 6 h history interval.

All dynamics tendencies are currently zero/identity placeholders. The event sink
is callback-based and carries no NetCDF, logging, or state-layout dependency.
It also records the originating schedule sequence index and nominal schedule
time so later state, physics, and nesting work can attach without changing the
tested 6 h schedule contract.

`tywrf::dynamics::CycleSchedule` freezes the v1 KROSA 6 h
boundary/nudging/nesting contract separately from the executable loop. It is a
plan object only: no dynamics, physics, interpolation, feedback, NetCDF, or
logging runs from this contract. The default schedule asserts:

- `540` d01 parent steps at `40 s`;
- `2700` d02 child substeps at `8 s` with `parent_time_step_ratio = 5`;
- d01 boundary and spectral nudging source refreshes at both `0 s` and
  `21600 s`, bracketing the 6 h segment;
- one d01 boundary-update call point and one d01 spectral-nudging call point on
  every parent step;
- one parent-child interpolation call point per d02 substep;
- one two-way feedback call point after each d01 step's d02 subcycles;
- moving-nest position update call points every nominal `15 min`, including
  both endpoints. Nominal times that do not land on the `40 s` parent-step grid
  are conservatively snapped to the next parent-step boundary, so the first
  nominal `900 s` check is scheduled at `920 s`;
- one d01 and one d02 history-output call point at `21600 s`.

## Phase 5 Physics Bridge ABI Baseline

The initial physics bridge is an ABI and staging contract only. It does not
execute Thompson, RRTMG, YSU, MM5 surface layer, Noah LSM, slab ocean, TC flux,
or Kain-Fritsch physics.

`bindings/wrf_physics/wrf_physics_bridge.h` defines the C-compatible boundary
that a future Fortran `ISO_C_BINDING` wrapper must mirror:

- `TywrfPhysicsGridMetadata` carries domain, KROSA grid spacing, time step,
  active mass dimensions, and full-level count.
- `TywrfPhysicsSuiteConfig` freezes the current KROSA physics option set:
  Thompson `8`, RRTMG `4/4`, YSU `1`, MM5 surface layer `1`, Noah `2`, slab
  ocean `1`, TC flux `2`, and KF cumulus only on d01.
- `TywrfPhysicsField2D` and `TywrfPhysicsField3D` carry POD field views:
  pointer, dimensions, strides, halo widths, and element size.
- `TywrfPhysicsStaging` groups the current required core state fields for a
  single domain physics call point.
- `tywrf_wrf_physics_step` is the stable C entry point. The current
  implementation validates staging and returns `stub_validated` for a valid
  no-op call with `executed_physics = 0`.

The C++ side in `include/tywrf/physics_bridge/staging.hpp` builds staging views
from `Grid` and mutable `StateView` without copying field data. Validation is
strict about the canonical field layout: `i` contiguous, `stride_k = nx`,
`stride_j = nx * nz` for 3D fields, and `stride_j = nx` for 2D fields. Required
field shapes follow the existing state layout: U/V staggered horizontally,
W/PH/PHB on full eta levels, mass fields on mass levels, and near-surface fields
on the mass horizontal grid.

### WRF Physics Entrypoint Strategy

The current PGWRF tree exposes WRF v4.6.1 through symlinks under:

```text
/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/model/WRFV4.6.1
```

The physics entrypoint strategy is to wrap WRF mediation-layer driver
subroutines first, not the individual scheme kernels. The drivers already
encode WRF option constants, selected optional arguments, timing accumulators,
and the call shape used by `dyn_em` for the fixed KROSA suite:

- `module_radiation_driver::radiation_driver` selects RRTMG longwave and
  shortwave for `ra_lw_physics = 4` and `ra_sw_physics = 4`.
- `module_surface_driver::surface_driver` dispatches MM5 surface layer,
  Noah LSM, and slab ocean through the selected KROSA surface options.
- `module_pbl_driver::pbl_driver` selects YSU for `bl_pbl_physics = 1`.
- `module_cumulus_driver::cumulus_driver` selects KF only on d01; d02 keeps
  `cu_physics = 0` and returns immediately.
- `module_microphysics_driver::microphysics_driver` selects Thompson for
  `mp_physics = 8`.

The WRF `dyn_em` call order observed in the source is:

1. `module_first_rk_step_part1.F`: radiation;
2. `module_first_rk_step_part1.F`: surface driver, including surface layer,
   Noah LSM, and slab ocean;
3. `module_first_rk_step_part1.F`: PBL;
4. `module_first_rk_step_part1.F`: cumulus;
5. `solve_em.F`: microphysics.

The next bridge should therefore add a Fortran `ISO_C_BINDING` wrapper that
exports the existing `tywrf_wrf_physics_step` symbol, reconstructs WRF-style
Fortran arrays from C staging, and calls only these fixed-suite drivers. It
should also add a one-time initialization boundary that mirrors the relevant
parts of `module_physics_init::phy_init`, including RRTMG table initialization,
`sfclayinit`, Noah table initialization, KF `kfinit`, and Thompson
`thompson_init`.

The current ABI v1 is enough for the validated no-op bridge, but not enough for
real WRF physics. Real calls need an ABI v2 sidecar or extension for:

- static and surface fields: `XLAT`, `XLONG`, `XLAND`, `XICE`, `HGT`, `TSK`,
  `SST`, `ALBEDO`, `EMISS`, `VEGFRA`, `IVGTYP`, `ISLTYP`;
- derived dynamics fields on WRF physics grids: `u_phy`, `v_phy`, `th_phy`,
  `t_phy`, `p_hyd`, `p_hyd_w`, `pi_phy`, `rho`, `dz8w`, `z`, `z_at_w`;
- physics tendencies and accumulators: `rth*ten`, `rq*ten`, `ru*ten`,
  `rv*ten`, `RAINCV`, `SNOWNC`, `GRAUPELNC`, `STEPRA`, `STEPBL`, `STEPCU`,
  `NCA`, `W0AVG`;
- radiation state: `GLW`, `GSW`, `SWDOWN`, `CLDFRA`, cloud fractions, RRTMG
  flux arrays, effective radii, ozone/aerosol placeholders for the disabled
  KROSA options;
- surface and soil state: `SMOIS`, `SH2O`, `TSLB`, `SNOW`, `SNOWH`, `CANWAT`,
  `TMN`, `ZNT`, `UST`, `PBLH`, `HFX`, `QFX`, `LH`;
- slab-ocean state: `TML`, `T0ML`, `HML`, `H0ML`, `HUML`, `HVML`, `TMOML`
  and KROSA ocean timing scalars.

Fortran should not directly consume arbitrary strided C views for the first
real bridge. The safer first implementation is an explicit staging copy:
TyWRF canonical buffers are packed into WRF-order arrays
`real :: field(ims:ime,kms:kme,jms:jme)` and `real :: field2d(ims:ime,jms:jme)`,
the WRF driver mutates those arrays, and changed fields are unpacked back to
TyWRF state. That preserves CUDA-ready core layout while avoiding fragile
Fortran pointer remapping across halo offsets.
