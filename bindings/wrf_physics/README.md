# WRF Physics Bridge Entrypoint Audit

This manifest records the P6 source audit for the TyWRF-Core v1 physics bridge.
It is a design and staging plan only; no real WRF physics calls are implemented
here.

## Audited Tree

Primary path:

```text
/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/model/WRFV4.6.1
```

In the KROSA tree, `main` and `phys` are symlinks into:

```text
/home/zzy/Projects/tc_sim/pgwrf_2025wp11_wp12_gfs_analysis_bundle/PGWRF/model/WRFV4.6.1
```

The full bundle tree also contains `dyn_em`, `frame`, `Registry`, and
`configure.wrf`, so the audit used the symlinked KROSA `phys/main` paths for
scheme files and the bundle path for `dyn_em` call sites and build metadata.

## KROSA Physics Suite

Observed in the active `namelist.input`:

| Option | d01 | d02 | WRF target |
| --- | ---: | ---: | --- |
| `mp_physics` | 8 | 8 | Thompson |
| `cu_physics` | 1 | 0 | KF on d01, disabled on d02 |
| `ra_lw_physics` | 4 | 4 | RRTMG LW |
| `ra_sw_physics` | 4 | 4 | RRTMG SW |
| `bl_pbl_physics` | 1 | 1 | YSU |
| `sf_sfclay_physics` | 1 | 1 | MM5 surface layer |
| `sf_surface_physics` | 2 | 2 | Noah LSM |
| `sf_ocean_physics` | 1 | 1 | slab ocean / OML |
| `isftcflx` | 2 | 2 | TC flux option inside `SFCLAY` |

## Candidate Entrypoints

Use mediation-layer drivers as the first real wrapper boundary. Leaf schemes
remain useful for smoke spikes, but they lose WRF driver gating and optional
argument checks.

| Suite piece | Driver candidate | Leaf routine selected in KROSA |
| --- | --- | --- |
| Radiation | `module_radiation_driver::radiation_driver` | `module_ra_rrtmg_lw::rrtmg_lwrad`, `module_ra_rrtmg_sw::rrtmg_swrad` |
| Surface | `module_surface_driver::surface_driver` | `module_sf_sfclay::SFCLAY`, `module_sf_noahdrv::lsm`, `module_sf_ocean_driver::OCEAN_DRIVER` |
| PBL | `module_pbl_driver::pbl_driver` | `module_bl_ysu::ysu` |
| Cumulus | `module_cumulus_driver::cumulus_driver` | `module_cu_kf::KFCPS` on d01 only |
| Microphysics | `module_microphysics_driver::microphysics_driver` | `module_mp_thompson::mp_gt_driver` |

Initialization candidates:

| Need | Candidate |
| --- | --- |
| WRF-wide physics setup | `module_physics_init::phy_init` |
| Radiation tables | `rrtmg_lwinit`, `rrtmg_swinit` through `module_physics_init` |
| Surface layer lookup | `module_sf_sfclay::sfclayinit` |
| Noah tables and state | `module_physics_init::landuse_init`, Noah `lsminit` path |
| KF tendencies | `module_cu_kf::kfinit` |
| Thompson tables and aerosols | `module_mp_thompson::thompson_init` |

The observed `dyn_em` order is radiation, surface, PBL, cumulus, then
microphysics. Radiation, surface, PBL, and cumulus are called from
`dyn_em/module_first_rk_step_part1.F`; microphysics is called from
`dyn_em/solve_em.F`.

## ABI Boundary

The exported C symbol should remain:

```c
int32_t tywrf_wrf_physics_step(
    const TywrfPhysicsStaging* staging,
    TywrfPhysicsDiagnostics* diagnostics);
```

The Fortran side should provide:

```fortran
subroutine tywrf_wrf_physics_step(staging, diagnostics) bind(C, name="tywrf_wrf_physics_step")
```

The wrapper should `use iso_c_binding` and WRF modules, convert C pointers into
temporary WRF-order arrays, call the fixed KROSA driver sequence, then unpack
mutated fields back into TyWRF buffers.

Do not pass haloed or strided C views directly to WRF drivers in the first real
implementation. WRF driver arguments are assumed-shape arrays with explicit
Fortran lower bounds such as `(ims:ime,kms:kme,jms:jme)` and
`(ims:ime,jms:jme)`. A pack/call/unpack layer is the lowest-risk bridge and
keeps the TyWRF core layout independent from WRF's Fortran ABI.

## Current ABI Gap

`TywrfPhysicsStaging` v1 only carries core prognostic and near-surface fields.
That is enough for the current no-op validation stub. It is not enough for real
WRF physics drivers.

The header now keeps that v1 struct and `tywrf_wrf_physics_step` behavior
frozen, and adds an ABI v2 sidecar scaffold for the first SFCLAY leaf smoke.
The new `_ex` entry point accepts the unchanged v1 staging pointer plus a chain
of v2 blocks beginning with `TywrfPhysicsBlockHeader`. The implemented v2
validator requires driver context, derived state, static/mask, and SFCLAY
surface blocks and can distinguish a complete sidecar from missing required
field pointers. `_ex` still reports `wrapper_unavailable` with
`executed_physics = 0`; it is not a real WRF physics producer.

The C++ `tywrf::physics_bridge::SidecarFixtureV2` helper is scaffold only. It
owns finite dummy ABI v2 buffers and block headers so unit tests and the next
SFCLAY leaf smoke can construct a complete sidecar chain from existing v1
staging without hand-written fixture wiring. Its provenance is explicitly
`scaffold/finite dummy sidecar; not executed physics`; the values are validator
fixtures, not WRF outputs. The first real physics wrapper remains the future
SFCLAY leaf smoke, not this helper.

ABI v2 needs either a new struct or an extension block for these staging groups:

The detailed field inventory, capability grouping, and backward-compatible ABI
strategy now live in `ABI_V2_MANIFEST.md`; this README keeps only the audit
summary.

| Group | Required fields |
| --- | --- |
| Derived physics state | `u_phy`, `v_phy`, `th_phy`, `t_phy`, `p_hyd`, `p_hyd_w`, `pi_phy`, `rho`, `dz8w`, `z`, `z_at_w` |
| Static and mask fields | `XLAT`, `XLONG`, `XLAND`, `XICE`, `HGT`, `LU_INDEX`, `IVGTYP`, `ISLTYP`, `VEGFRA` |
| Surface diagnostics | `TSK`, `SST`, `ALBEDO`, `EMISS`, `ZNT`, `UST`, `PBLH`, `HFX`, `QFX`, `LH`, `QSFC`, `WSPD`, `BR`, `FM`, `FHH` |
| Soil and snow | `SMOIS`, `SH2O`, `TSLB`, `SNOW`, `SNOWH`, `SNOWC`, `CANWAT`, `TMN`, `SNOALB`, `MAVAIL` |
| Tendencies | `RTHRATEN`, `RTHRATENLW`, `RTHRATENSW`, `RTHBLTEN`, `RQVBLTEN`, `RUBLTEN`, `RVBLTEN`, `RTHCUTEN`, `RQVCUTEN` and active hydrometeor tendency arrays |
| Accumulators | `RAINCV`, `RAINNCV`, `SNOWNC`, `SNOWNCV`, `GRAUPELNC`, `GRAUPELNCV`, `NCA`, `W0AVG`, `STEPRA`, `STEPBL`, `STEPCU` |
| Radiation extras | `GLW`, `GSW`, `SWDOWN`, `CLDFRA`, `CFRACH`, `CFRACL`, `CFRACM`, `RE_CLOUD`, `RE_ICE`, `RE_SNOW`, RRTMG flux arrays |
| Slab ocean | `TML`, `T0ML`, `HML`, `H0ML`, `HUML`, `HVML`, `TMOML`, `OML_GAMMA`, `OML_RELAXATION_TIME` |

Several fields can be initialized from `wrfinput` or copied from WRF templates
for the first spike, but they still need explicit ownership in the bridge.

## Dependency Risks

- The WRF build is compiler-specific: `configure.wrf` uses Intel oneAPI
  `ifx/icx`, `DMPARALLEL=1`, `RWORDSIZE=4`, `WRF_CHEM=0`, and
  `ARCH_LOCAL=-DNONSTANDARD_SYSTEM_FUNC -DWRF_USE_CLM`. The wrapper must be
  compiled with the same Fortran compiler family and WRF module files.
- The compiled library is available as `main/libwrflib.a`; linking a wrapper
  against it also requires WRF external libraries from `configure.wrf`,
  including NetCDF, HDF5, RSL_LITE, and WRF I/O archives.
- WRF initialization and selected schemes read runtime data files from the WRF
  run directory: RRTMG data, Thompson aerosol lookup data, and Noah land-use
  tables are the first expected file dependencies.
- WRF error paths use `wrf_error_fatal`, which can abort rather than return a
  status. The first spike should run as a separate executable or subprocess
  until failure handling is understood.
- Several modules keep `SAVE` or module state, including RRTMG tables,
  Thompson lookup tables, and SFCLAY lookup arrays. Treat WRF physics as
  process-global and non-reentrant until proven otherwise. Do not call it
  concurrently from multiple TyWRF domains or OpenMP threads.
- The selected drivers expect WRF tile/domain/memory indices and optional
  arguments to be internally consistent. Missing an optional argument for the
  active scheme often reaches `wrf_error_fatal`.
- TyWRF U/V are staggered core fields. WRF physics uses mass-grid `u_phy` and
  `v_phy`; interpolation to mass points must happen before staging.

## Next Minimal Spike

1. Add a standalone Fortran wrapper source outside the WRF tree, for example
   `bindings/wrf_physics/tywrf_wrf_physics_wrapper.F90`.
2. Include `wrf_physics_bridge.h` through matching `BIND(C)` Fortran derived
   types or a generated ABI module. Keep the exported C symbol unchanged.
3. Build the wrapper with the WRF `configure.wrf` compiler, flags, module
   include directories, and `main/libwrflib.a`. Do this in TyWRF-Core build
   output, not inside the WRF source tree.
4. Implement `tywrf_wrf_physics_initialize` as a private spike entrypoint or
   test-only symbol that calls the fixed KROSA initialization path.
5. Pack one small d01 tile from reference `wrfinput`/`wrfout` into WRF-order
   staging arrays and call only `module_sf_sfclay::SFCLAY`. Verify the wrapper
   links, runs, and mutates `U10`, `V10`, `T2`, `Q2`, `HFX`, and `QFX`.
6. Expand the spike in this order: surface driver with Noah disabled for a
   minimal smoke, full `surface_driver`, YSU `pbl_driver`, Thompson
   `microphysics_driver`, RRTMG radiation, then d01 KF.
7. Only after the wrapper survives one real driver should TyWRF-Core promote
   ABI v2 fields into the C++ staging API.

## Audit Commands

Representative commands used for this audit:

```bash
rg -n "SUBROUTINE\\s+(microphysics_driver|radiation_driver|pbl_driver|surface_driver|cumulus_driver)" \
  /home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/model/WRFV4.6.1/phys

rg -n "CALL\\s+(radiation_driver|surface_driver|pbl_driver|cumulus_driver|microphysics_driver)" \
  /home/zzy/Projects/tc_sim/pgwrf_2025wp11_wp12_gfs_analysis_bundle/PGWRF/model/WRFV4.6.1/dyn_em

rg -n "rrtmg_lwinit|rrtmg_swinit|sfclayinit|kfinit|thompson_init|phy_init" \
  /home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/model/WRFV4.6.1/phys

rg -n "^SFC|^SCC|^FC|^DM_FC|^RWORDSIZE|^WRF_CHEM|^ARCH_LOCAL|^LIBWRFLIB" \
  /home/zzy/Projects/tc_sim/pgwrf_2025wp11_wp12_gfs_analysis_bundle/PGWRF/model/WRFV4.6.1/configure.wrf
```
