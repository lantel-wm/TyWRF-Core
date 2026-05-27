#ifndef TYWRF_WRF_PHYSICS_BRIDGE_H
#define TYWRF_WRF_PHYSICS_BRIDGE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define TYWRF_PHYSICS_ABI_VERSION_V1 1
#define TYWRF_PHYSICS_ABI_VERSION_V2 2
#define TYWRF_PHYSICS_ABI_VERSION TYWRF_PHYSICS_ABI_VERSION_V1

#define TYWRF_PHYSICS_CAPABILITY_NONE UINT64_C(0)
#define TYWRF_PHYSICS_CAPABILITY_DRIVER_CONTEXT (UINT64_C(1) << 0)
#define TYWRF_PHYSICS_CAPABILITY_DERIVED_STATE (UINT64_C(1) << 1)
#define TYWRF_PHYSICS_CAPABILITY_STATIC_MASK (UINT64_C(1) << 2)
#define TYWRF_PHYSICS_CAPABILITY_SURFACE_STATE (UINT64_C(1) << 3)
#define TYWRF_PHYSICS_CAPABILITY_SOIL_SNOW (UINT64_C(1) << 4)
#define TYWRF_PHYSICS_CAPABILITY_TENDENCIES (UINT64_C(1) << 5)
#define TYWRF_PHYSICS_CAPABILITY_ACCUMULATORS (UINT64_C(1) << 6)
#define TYWRF_PHYSICS_CAPABILITY_RADIATION_RRTMG (UINT64_C(1) << 7)
#define TYWRF_PHYSICS_CAPABILITY_SLAB_OCEAN (UINT64_C(1) << 8)
#define TYWRF_PHYSICS_CAPABILITY_PROCESS_GLOBAL (UINT64_C(1) << 9)

typedef struct TywrfPhysicsBlockHeader {
  uint32_t struct_size;
  uint32_t abi_version;
  uint64_t capabilities;
  const void* next;
} TywrfPhysicsBlockHeader;

#define TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE \
  ((uint32_t)sizeof(TywrfPhysicsBlockHeader))

static inline int tywrf_physics_block_header_has_min_size(
    const TywrfPhysicsBlockHeader* header,
    const uint32_t min_size) {
  return header != NULL && header->struct_size >= min_size;
}

static inline int tywrf_physics_block_header_has_abi(
    const TywrfPhysicsBlockHeader* header,
    const uint32_t abi_version) {
  return header != NULL && header->abi_version == abi_version;
}

typedef enum TywrfPhysicsStatus {
  TYWRF_PHYSICS_STATUS_STUB_VALIDATED = 0,
  TYWRF_PHYSICS_STATUS_NULL_ARGUMENT = 1,
  TYWRF_PHYSICS_STATUS_UNSUPPORTED_ABI = 2,
  TYWRF_PHYSICS_STATUS_UNSUPPORTED_DOMAIN = 3,
  TYWRF_PHYSICS_STATUS_UNSUPPORTED_SUITE = 4,
  TYWRF_PHYSICS_STATUS_INVALID_TIMING = 5,
  TYWRF_PHYSICS_STATUS_INVALID_DIMENSIONS = 6,
  TYWRF_PHYSICS_STATUS_INVALID_STRIDES = 7,
  TYWRF_PHYSICS_STATUS_INVALID_ELEMENT_SIZE = 8,
  TYWRF_PHYSICS_STATUS_MISSING_REQUIRED_FIELD = 9,
  TYWRF_PHYSICS_STATUS_WRAPPER_UNAVAILABLE = 10
} TywrfPhysicsStatus;

typedef enum TywrfPhysicsDomain {
  TYWRF_PHYSICS_DOMAIN_D01 = 1,
  TYWRF_PHYSICS_DOMAIN_D02 = 2
} TywrfPhysicsDomain;

typedef enum TywrfPhysicsFieldId {
  TYWRF_PHYSICS_FIELD_NONE = 0,
  TYWRF_PHYSICS_FIELD_U = 1,
  TYWRF_PHYSICS_FIELD_V = 2,
  TYWRF_PHYSICS_FIELD_W = 3,
  TYWRF_PHYSICS_FIELD_PH = 4,
  TYWRF_PHYSICS_FIELD_PHB = 5,
  TYWRF_PHYSICS_FIELD_T = 6,
  TYWRF_PHYSICS_FIELD_P = 7,
  TYWRF_PHYSICS_FIELD_PB = 8,
  TYWRF_PHYSICS_FIELD_QVAPOR = 9,
  TYWRF_PHYSICS_FIELD_QCLOUD = 10,
  TYWRF_PHYSICS_FIELD_QRAIN = 11,
  TYWRF_PHYSICS_FIELD_QICE = 12,
  TYWRF_PHYSICS_FIELD_QSNOW = 13,
  TYWRF_PHYSICS_FIELD_QGRAUP = 14,
  TYWRF_PHYSICS_FIELD_QNICE = 15,
  TYWRF_PHYSICS_FIELD_QNRAIN = 16,
  TYWRF_PHYSICS_FIELD_MU = 17,
  TYWRF_PHYSICS_FIELD_MUB = 18,
  TYWRF_PHYSICS_FIELD_PSFC = 19,
  TYWRF_PHYSICS_FIELD_U10 = 20,
  TYWRF_PHYSICS_FIELD_V10 = 21,
  TYWRF_PHYSICS_FIELD_T2 = 22,
  TYWRF_PHYSICS_FIELD_Q2 = 23,
  TYWRF_PHYSICS_FIELD_RAINC = 24,
  TYWRF_PHYSICS_FIELD_RAINNC = 25,
  TYWRF_PHYSICS_FIELD_TSK = 26,
  TYWRF_PHYSICS_FIELD_ZNT = 27,
  TYWRF_PHYSICS_FIELD_UST = 28,
  TYWRF_PHYSICS_FIELD_PBLH = 29,
  TYWRF_PHYSICS_FIELD_HFX = 30,
  TYWRF_PHYSICS_FIELD_QFX = 31,
  TYWRF_PHYSICS_FIELD_LH = 32,
  TYWRF_PHYSICS_FIELD_TH2 = 33,
  TYWRF_PHYSICS_FIELD_XLAND = 34,
  TYWRF_PHYSICS_FIELD_LAKEMASK = 35,
  TYWRF_PHYSICS_FIELD_LU_INDEX = 36,
  TYWRF_PHYSICS_FIELD_HGT = 37,
  TYWRF_PHYSICS_FIELD_XLAT = 38,
  TYWRF_PHYSICS_FIELD_XLONG = 39,
  TYWRF_PHYSICS_FIELD_U_PHY = 40,
  TYWRF_PHYSICS_FIELD_V_PHY = 41,
  TYWRF_PHYSICS_FIELD_T_PHY = 42,
  TYWRF_PHYSICS_FIELD_QV_CURR = 43,
  TYWRF_PHYSICS_FIELD_P_PHY = 44,
  TYWRF_PHYSICS_FIELD_DZ8W = 45
} TywrfPhysicsFieldId;

typedef struct TywrfPhysicsField2D {
  void* data;
  int32_t nx;
  int32_t ny;
  int32_t stride_i;
  int32_t stride_j;
  int32_t halo_i_lower;
  int32_t halo_i_upper;
  int32_t halo_j_lower;
  int32_t halo_j_upper;
  int32_t element_bytes;
} TywrfPhysicsField2D;

typedef struct TywrfPhysicsField3D {
  void* data;
  int32_t nx;
  int32_t ny;
  int32_t nz;
  int32_t stride_i;
  int32_t stride_k;
  int32_t stride_j;
  int32_t halo_i_lower;
  int32_t halo_i_upper;
  int32_t halo_j_lower;
  int32_t halo_j_upper;
  int32_t halo_k_lower;
  int32_t halo_k_upper;
  int32_t element_bytes;
} TywrfPhysicsField3D;

typedef struct TywrfPhysicsGridMetadata {
  int32_t domain_id;
  int32_t mass_nx;
  int32_t mass_ny;
  int32_t mass_nz;
  int32_t full_nz;
  double dx_m;
  double dy_m;
  double dt_s;
  int64_t step_index;
  double start_seconds;
  double end_seconds;
} TywrfPhysicsGridMetadata;

typedef struct TywrfPhysicsSuiteConfig {
  int32_t mp_physics;
  int32_t cu_physics;
  int32_t ra_lw_physics;
  int32_t ra_sw_physics;
  int32_t bl_pbl_physics;
  int32_t sf_sfclay_physics;
  int32_t sf_surface_physics;
  int32_t sf_ocean_physics;
  int32_t isftcflx;
  int32_t num_moist_species;
  int32_t num_soil_layers;
} TywrfPhysicsSuiteConfig;

typedef struct TywrfPhysicsStaging {
  int32_t abi_version;
  TywrfPhysicsGridMetadata grid;
  TywrfPhysicsSuiteConfig suite;

  TywrfPhysicsField3D u;
  TywrfPhysicsField3D v;
  TywrfPhysicsField3D w;
  TywrfPhysicsField3D ph;
  TywrfPhysicsField3D phb;
  TywrfPhysicsField3D t;
  TywrfPhysicsField3D p;
  TywrfPhysicsField3D pb;
  TywrfPhysicsField3D qvapor;
  TywrfPhysicsField3D qcloud;
  TywrfPhysicsField3D qrain;
  TywrfPhysicsField3D qice;
  TywrfPhysicsField3D qsnow;
  TywrfPhysicsField3D qgraup;
  TywrfPhysicsField3D qnice;
  TywrfPhysicsField3D qnrain;

  TywrfPhysicsField2D mu;
  TywrfPhysicsField2D mub;
  TywrfPhysicsField2D psfc;
  TywrfPhysicsField2D u10;
  TywrfPhysicsField2D v10;
  TywrfPhysicsField2D t2;
  TywrfPhysicsField2D q2;
  TywrfPhysicsField2D rainc;
  TywrfPhysicsField2D rainnc;
} TywrfPhysicsStaging;

typedef struct TywrfPhysicsDiagnostics {
  int32_t status;
  int32_t failing_field;
  int32_t validated_field_count;
  int32_t executed_physics;
} TywrfPhysicsDiagnostics;

typedef struct TywrfPhysicsConstantsV2 {
  double g;
  double cp;
  double r_d;
  double r_v;
  double p1000mb;
  double t0;
} TywrfPhysicsConstantsV2;

typedef struct TywrfPhysicsDriverContextV2 {
  TywrfPhysicsBlockHeader header;
  int32_t domain_id;
  int32_t mass_nx;
  int32_t mass_ny;
  int32_t mass_nz;
  int32_t full_nz;
  int32_t ids;
  int32_t ide;
  int32_t jds;
  int32_t jde;
  int32_t kds;
  int32_t kde;
  int32_t ims;
  int32_t ime;
  int32_t jms;
  int32_t jme;
  int32_t kms;
  int32_t kme;
  int32_t its;
  int32_t ite;
  int32_t jts;
  int32_t jte;
  int32_t kts;
  int32_t kte;
  int32_t enable_sfclay;
  int32_t enable_surface_driver;
  int32_t enable_pbl_driver;
  int32_t enable_cumulus_driver;
  int32_t enable_microphysics_driver;
  int32_t enable_radiation_driver;
  int64_t step_index;
  double dx_m;
  double dy_m;
  double dt_s;
  double start_seconds;
  double end_seconds;
  double xtime_minutes;
  TywrfPhysicsSuiteConfig suite;
  const TywrfPhysicsConstantsV2* constants;
} TywrfPhysicsDriverContextV2;

typedef struct TywrfPhysicsDerivedStateV2 {
  TywrfPhysicsBlockHeader header;
  TywrfPhysicsField3D u_phy;
  TywrfPhysicsField3D v_phy;
  TywrfPhysicsField3D t_phy;
  TywrfPhysicsField3D qv_curr;
  TywrfPhysicsField3D p_phy;
  TywrfPhysicsField3D dz8w;
} TywrfPhysicsDerivedStateV2;

typedef struct TywrfPhysicsStaticMaskV2 {
  TywrfPhysicsBlockHeader header;
  TywrfPhysicsField2D xland;
  TywrfPhysicsField2D lakemask;
  TywrfPhysicsField2D lu_index;
  TywrfPhysicsField2D hgt;
  TywrfPhysicsField2D xlat;
  TywrfPhysicsField2D xlong;
} TywrfPhysicsStaticMaskV2;

typedef struct TywrfPhysicsSfclaySurfaceV2 {
  TywrfPhysicsBlockHeader header;
  TywrfPhysicsField2D psfc;
  TywrfPhysicsField2D tsk;
  TywrfPhysicsField2D znt;
  TywrfPhysicsField2D ust;
  TywrfPhysicsField2D pblh;
  TywrfPhysicsField2D hfx;
  TywrfPhysicsField2D qfx;
  TywrfPhysicsField2D lh;
  TywrfPhysicsField2D u10;
  TywrfPhysicsField2D v10;
  TywrfPhysicsField2D t2;
  TywrfPhysicsField2D q2;
  TywrfPhysicsField2D th2;
  TywrfPhysicsField2D scratch_1;
  TywrfPhysicsField2D scratch_2;
} TywrfPhysicsSfclaySurfaceV2;

#define TYWRF_PHYSICS_CONSTANTS_V2_STRUCT_SIZE \
  ((uint32_t)sizeof(TywrfPhysicsConstantsV2))
#define TYWRF_PHYSICS_DRIVER_CONTEXT_V2_STRUCT_SIZE \
  ((uint32_t)sizeof(TywrfPhysicsDriverContextV2))
#define TYWRF_PHYSICS_DERIVED_STATE_V2_STRUCT_SIZE \
  ((uint32_t)sizeof(TywrfPhysicsDerivedStateV2))
#define TYWRF_PHYSICS_STATIC_MASK_V2_STRUCT_SIZE \
  ((uint32_t)sizeof(TywrfPhysicsStaticMaskV2))
#define TYWRF_PHYSICS_SFCLAY_SURFACE_V2_STRUCT_SIZE \
  ((uint32_t)sizeof(TywrfPhysicsSfclaySurfaceV2))

/*
 * Stable C entry point for the future WRF physics wrapper.
 *
 * Fortran ISO_C_BINDING should mirror the structs above with BIND(C), represent
 * field pointers as TYPE(C_PTR), and export this exact symbol with
 * BIND(C, NAME="tywrf_wrf_physics_step"). The current C++ implementation is a
 * deterministic no-op stub: it validates metadata and field views, sets
 * diagnostics.executed_physics to zero, and returns
 * TYWRF_PHYSICS_STATUS_STUB_VALIDATED for a valid staging request.
 */
int32_t tywrf_wrf_physics_step(
    const TywrfPhysicsStaging* staging,
    TywrfPhysicsDiagnostics* diagnostics);

/*
 * ABI v2 sidecar validation and future extended entry point.
 *
 * The sidecar chain is additive to the frozen v1 staging struct. Every v2 block
 * starts with TywrfPhysicsBlockHeader and is linked through header.next. The
 * current implementation validates a minimal SFCLAY-ready sidecar set but does
 * not call WRF. A valid tywrf_wrf_physics_step_ex request returns
 * TYWRF_PHYSICS_STATUS_WRAPPER_UNAVAILABLE with diagnostics.executed_physics
 * set to zero.
 */
int32_t tywrf_physics_validate_sidecar_v2(
    const TywrfPhysicsStaging* staging,
    const TywrfPhysicsBlockHeader* sidecars,
    TywrfPhysicsDiagnostics* diagnostics);

int32_t tywrf_wrf_physics_step_ex(
    const TywrfPhysicsStaging* staging,
    const TywrfPhysicsBlockHeader* sidecars,
    TywrfPhysicsDiagnostics* diagnostics);

const char* tywrf_physics_status_name(int32_t status);

#ifdef __cplusplus
}
#endif

#endif
