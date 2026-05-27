# WRF Physics ABI v2 Manifest

This manifest turns the P7 read-only physics staging research into an input for
future C header layout work and wrapper tests. It is documentation only. It does
not change the v1 C ABI, CMake targets, tests, or any Fortran wrapper source.

## Scope

ABI v2 is intended to carry the fields that real WRF v4.6.1 KROSA physics
drivers need beyond the current no-op `TywrfPhysicsStaging` v1 contract.

The target suite remains the fixed PGWRF/KROSA subset:

- Thompson microphysics, `mp_physics = 8`;
- Kain-Fritsch cumulus on d01 only, `cu_physics = 1, 0`;
- RRTMG longwave and shortwave radiation, `ra_lw_physics = 4`,
  `ra_sw_physics = 4`;
- YSU PBL, `bl_pbl_physics = 1`;
- MM5 surface layer, Noah LSM, slab ocean, and TC flux,
  `sf_sfclay_physics = 1`, `sf_surface_physics = 2`,
  `sf_ocean_physics = 1`, `isftcflx = 2`.

The manifest uses these source categories:

- `State`: owned by TyWRF core state or the existing v1 staging fields.
- `wrfinput`: read from `wrfinput_d01` or the matching moving-nest input state.
- `wrfout template`: copied or initialized from a WRF reference/template output
  until TyWRF owns the producer.
- `derived`: computed during pack/staging from State, static fields, grid
  metadata, or thermodynamic formulas.
- `persistent`: owned across physics calls, either per domain or
  process-global.

The manifest uses these shape categories:

- `2D`: mass-horizontal field.
- `3D`: mass-level field.
- `full-level`: W-level/full eta field.
- `scalar`: single value or short fixed option vector.
- `open set`: scheme-dependent repeated fields that should be represented by
  counted arrays or capability-specific blocks, not by an unbounded flat struct.

## Capability Groups

Future headers should expose capability bits or equivalent feature flags for
these groups. Names below are proposed identifiers, not current ABI symbols.

| Capability group | Purpose |
| --- | --- |
| `driver_context` | Domain, time, grid, tile, memory, and suite metadata. |
| `derived_state` | WRF physics-grid state derived from TyWRF prognostic fields. |
| `static_mask` | Geolocation, terrain, map factors, land/sea/ice masks, and land-use IDs. |
| `surface_state` | Surface diagnostics and exchange fields used by SFCLAY, Noah, YSU, and radiation. |
| `soil_snow` | Noah soil, liquid water, skin, canopy, snow, and runoff state. |
| `tendencies` | WRF physics tendency arrays consumed or produced by active schemes. |
| `accumulators` | Precipitation, convection, timing, and scheme-step accumulators. |
| `radiation_rrtmg` | RRTMG inputs, diagnostics, cloud fractions, effective radii, and fluxes. |
| `slab_ocean` | Slab ocean mixed-layer state and timing/scalar controls. |
| `process_global` | WRF module initialization, lookup tables, runtime data, and non-reentrant state. |

## Field Inventory

`v2 required` means required for the full KROSA driver wrapper target. `Can
defer` means the field group can be omitted from the first header/layout test or
from a narrow leaf smoke as long as the corresponding capability bit is absent.

| Layer | Fields | Source | Shape category | v2 required | Can defer | Capability group |
| --- | --- | --- | --- | --- | --- | --- |
| Driver context | `domain_id`, `mass_nx`, `mass_ny`, `mass_nz`, `full_nz`, `dx_m`, `dy_m`, `dt_s`, `step_index`, `start_seconds`, `end_seconds` | State | scalar | yes | no | `driver_context` |
| Driver context | WRF memory/domain/tile bounds: `ids:ide`, `jds:jde`, `kds:kde`, `ims:ime`, `jms:jme`, `kms:kme`, `ips:ipe`, `jps:jpe`, `kps:kpe`, tile starts/ends | derived | scalar | yes | no | `driver_context` |
| Driver context | suite options: `mp_physics`, `cu_physics`, `ra_lw_physics`, `ra_sw_physics`, `bl_pbl_physics`, `sf_sfclay_physics`, `sf_surface_physics`, `sf_ocean_physics`, `isftcflx`, `num_moist_species`, `num_soil_layers` | State | scalar | yes | no | `driver_context` |
| Driver context | physics clocks and gates: `itimestep`, `gmt`, `julian`, `xtime`, radiation/PBL/cumulus intervals, current history-cycle seconds | derived | scalar | yes | no | `driver_context` |
| Driver context | constants and config scalars used by drivers, including `g`, `cp`, `rd`, `rv`, `p1000mb`, `r_d`, and KROSA namelist switches not already in v1 | derived | scalar | yes | no | `driver_context` |
| Derived physics state | `u_phy`, `v_phy` mass-grid winds packed from staggered TyWRF `U` and `V` | derived | 3D | yes | only for non-wind leaf smoke | `derived_state` |
| Derived physics state | `th_phy`, `t_phy`, `p_hyd`, `pi_phy`, `rho`, `alt`, `moist`, active hydrometeor mixing ratios on mass levels | State, derived | 3D | yes | no for real driver calls | `derived_state` |
| Derived physics state | `p_hyd_w`, `dz8w`, `z`, `z_at_w`, geopotential-derived height and layer thickness fields | State, derived | full-level | yes | no for PBL/radiation | `derived_state` |
| Derived physics state | edge or helper fields required by selected drivers, such as lowest-model-level pressure/temperature and reciprocal density work arrays | derived | 2D, 3D, open set | yes | yes until a driver proves the exact need | `derived_state` |
| Static and mask | `XLAT`, `XLONG`, `HGT`, `XLAND`, `LANDMASK`, `XICE`, `SEAICE`, `LU_INDEX`, `IVGTYP`, `ISLTYP`, `VEGFRA`, `LAI` | wrfinput, wrfout template | 2D | yes | no for surface/radiation | `static_mask` |
| Static and mask | map and rotation fields: `MAPFAC_M`, `MAPFAC_U`, `MAPFAC_V`, `MAPFAC_MX`, `MAPFAC_MY`, `F`, `E`, `SINALPHA`, `COSALPHA` | wrfinput, wrfout template | 2D | yes | yes until a wrapped driver consumes them | `static_mask` |
| Static and mask | soil and land category counts, land-use table IDs, dominant category metadata | wrfinput, derived | scalar | yes | yes for SFCLAY-only leaf smoke | `static_mask` |
| Surface state | `TSK`, `SST`, `ALBEDO`, `EMISS`, `ZNT`, `UST`, `MOL`, `PBLH`, `PSFC`, `U10`, `V10`, `T2`, `Q2`, `TH2` | State, wrfinput, wrfout template, persistent | 2D | yes | no for surface/PBL | `surface_state` |
| Surface state | surface flux and exchange fields: `HFX`, `QFX`, `LH`, `GRDFLX`, `QSFC`, `WSPD`, `BR`, `FM`, `FHH`, `FLHC`, `FLQC`, `CHS`, `CHS2`, `CQS2` | wrfout template, derived, persistent | 2D | yes | no for SFCLAY/Noah/YSU | `surface_state` |
| Surface state | optional diagnostic surfaces used by active WRF branches, including roughness, wetness, emissivity, and surface category work arrays | wrfout template, derived | 2D, open set | yes when capability present | yes | `surface_state` |
| Soil, snow, canopy | `SMOIS`, `SH2O`, `TSLB` across Noah soil layers | wrfinput, persistent | 3D | yes for Noah | yes for SFCLAY-only leaf smoke | `soil_snow` |
| Soil, snow, canopy | `SNOW`, `SNOWH`, `SNOWC`, `CANWAT`, `TMN`, `SNOALB`, `MAVAIL`, `SFROFF`, `UDROFF` | wrfinput, wrfout template, persistent | 2D | yes for Noah | yes for SFCLAY-only leaf smoke | `soil_snow` |
| Soil, snow, canopy | Noah table selections and per-category defaults needed by `landuse_init`/Noah paths | persistent | scalar, open set | yes for Noah | yes until Noah driver smoke | `soil_snow` |
| Tendencies | heat tendencies: `RTHRATEN`, `RTHRATENLW`, `RTHRATENSW`, `RTHBLTEN`, `RTHCUTEN` | wrfout template, derived, persistent | 3D | yes | no for full driver sequence | `tendencies` |
| Tendencies | moisture tendencies: `RQVBLTEN`, `RQVCUTEN`, `RQCCUTEN`, `RQRCUTEN`, `RQICUTEN`, `RQSCUTEN`, `RQGCUTEN` and active Thompson/KF hydrometeor tendencies | wrfout template, derived, persistent | 3D, open set | yes | d02 cumulus entries can defer | `tendencies` |
| Tendencies | momentum tendencies: `RUBLTEN`, `RVBLTEN`, `RUCUTEN`, `RVCUTEN` | wrfout template, derived, persistent | 3D | yes | d02 cumulus entries can defer | `tendencies` |
| Tendencies | tendency work arrays local to a wrapped driver that do not need long-term ownership | derived | open set | yes when driver requires | yes | `tendencies` |
| Accumulators | v1 precipitation state: `RAINC`, `RAINNC` | State | 2D | yes | no | `accumulators` |
| Accumulators | precipitation increments and species totals: `RAINCV`, `RAINNCV`, `SNOWNC`, `SNOWNCV`, `GRAUPELNC`, `GRAUPELNCV`, optional hail totals if active | wrfout template, persistent | 2D | yes | no for Thompson/KF | `accumulators` |
| Accumulators | convection and vertical-motion accumulators: `NCA`, `W0AVG`, KF trigger/state counters | wrfout template, persistent | 2D, open set | yes for d01 KF | yes for d02 and non-cumulus smoke | `accumulators` |
| Accumulators | scheme step counters and intervals: `STEPRA`, `STEPBL`, `STEPCU`, microphysics/radiation/PBL elapsed counters | persistent, derived | scalar | yes | no for driver sequence tests | `accumulators` |
| Radiation extras | radiation diagnostics: `GLW`, `GSW`, `SWDOWN`, `OLR`, `XLAT`, solar angle/work scalars visible to radiation | wrfout template, derived, persistent | 2D, scalar | yes for RRTMG | yes for non-radiation smoke | `radiation_rrtmg` |
| Radiation extras | cloud fields: `CLDFRA`, `CFRACH`, `CFRACL`, `CFRACM`, cloud liquid/ice/rain/snow/graupel inputs used by RRTMG | State, wrfout template, derived | 3D, 2D | yes for RRTMG | yes for non-radiation smoke | `radiation_rrtmg` |
| Radiation extras | effective radii and optical inputs: `RE_CLOUD`, `RE_ICE`, `RE_SNOW`, ozone/aerosol placeholders for disabled KROSA options | wrfout template, derived, persistent | 3D, open set | yes when radiation capability present | yes | `radiation_rrtmg` |
| Radiation extras | RRTMG flux arrays and heating-rate work fields on layer interfaces | wrfout template, persistent | full-level, open set | yes for RRTMG | yes for non-radiation smoke | `radiation_rrtmg` |
| Slab ocean | mixed-layer state: `TML`, `T0ML`, `HML`, `H0ML`, `HUML`, `HVML`, `TMOML` | wrfinput, wrfout template, persistent | 2D | yes when `sf_ocean_physics = 1` | yes for SFCLAY-only leaf smoke | `slab_ocean` |
| Slab ocean | controls: `OML_GAMMA`, `OML_RELAXATION_TIME`, slab-ocean update counters and interval scalars | wrfinput, derived, persistent | scalar | yes when `sf_ocean_physics = 1` | yes for SFCLAY-only leaf smoke | `slab_ocean` |
| Persistence and process-global state | WRF physics initialization status: `phy_init`, RRTMG `rrtmg_lwinit`/`rrtmg_swinit`, `sfclayinit`, Noah land-use tables, KF `kfinit`, Thompson lookup and aerosol tables | persistent | scalar, open set | yes for real wrapper | yes until wrapper implementation | `process_global` |
| Persistence and process-global state | runtime data paths and loaded table versions for RRTMG, Thompson, Noah, and surface layer lookup data | persistent | scalar, open set | yes for reproducible wrapper tests | yes until wrapper implementation | `process_global` |
| Persistence and process-global state | process-wide non-reentrant state guard, per-domain serialization, initialization reference counts, and last-error capture if WRF abort paths are isolated | persistent | scalar | yes before threaded/domain-concurrent calls | yes for single-process leaf smoke | `process_global` |

## Backward-Compatible ABI Strategy

Freeze `TywrfPhysicsStaging` v1 exactly as it is today. It remains the ABI for
the current no-op validation bridge and for any caller compiled against
`TYWRF_PHYSICS_ABI_VERSION == 1`.

ABI v2 should be additive. Acceptable implementation patterns are:

1. Keep `tywrf_wrf_physics_step(const TywrfPhysicsStaging*, ...)` as the v1
   entrypoint and add a new `_ex` symbol, for example
   `tywrf_wrf_physics_step_ex(const TywrfPhysicsStaging* v1,
   const TywrfPhysicsBlockHeader* blocks, TywrfPhysicsDiagnostics* diagnostics)`.
2. Keep v1 as the fixed prefix and pass v2 data as sidecar blocks. Each block
   starts with a common header:

   ```c
   typedef struct TywrfPhysicsBlockHeader {
     uint32_t struct_size;
     uint32_t abi_version;
     uint64_t capabilities;
     const void* next;
   } TywrfPhysicsBlockHeader;
   ```

3. Alternatively add a new `TywrfPhysicsStagingV2` plus a new entrypoint. The
   v2 struct must still begin with `struct_size`, `abi_version`, `capabilities`,
   and `next` so older validators can reject unknown layouts cleanly and newer
   validators can skip unknown extension blocks.

The header layout test should verify:

- v1 size, offsets, enum values, and entrypoint declaration remain unchanged;
- v2 block headers have stable alignment in C and C++;
- `struct_size` rejects short blocks and tolerates longer compatible blocks;
- `abi_version` is checked before any v2 field dereference;
- `capabilities` controls required block validation;
- `next` supports ordered extension scanning without requiring one monolithic
  v2 struct.

## Wrapper And Test Sequencing

This round intentionally does not add a Fortran wrapper.

Recommended next steps:

1. Add header-only ABI v2 layout declarations plus C/C++ layout tests from this
   manifest.
2. Add an independent SFCLAY leaf smoke after the layout is frozen. It should
   pack one small d01 tile from reference `wrfinput`/`wrfout` fields, call only
   the MM5 `SFCLAY` leaf, and verify mutation of `U10`, `V10`, `T2`, `Q2`,
   `HFX`, and `QFX`.
3. Promote broader driver wrappers only after the SFCLAY smoke proves compiler,
   WRF module, runtime-data, and pack/unpack assumptions.

## Risks

- Several WRF drivers require optional arguments whose exact active set is
  controlled by scheme options and compile-time flags. Missing one can call
  `wrf_error_fatal` instead of returning an error.
- WRF physics modules own process-global state and should be treated as
  non-reentrant until proven otherwise. Do not call the wrapper concurrently
  from multiple domains or OpenMP threads.
- Many fields can be bootstrapped from `wrfout template` files for the first
  smoke, but v2 must still record final ownership so template copying does not
  become hidden model state.
- Staggered TyWRF winds are not WRF physics winds. `u_phy` and `v_phy` must be
  mass-grid derived fields with explicit pack/unpack ownership.
