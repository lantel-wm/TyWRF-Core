#include "tywrf/grid.hpp"
#include "tywrf/physics_bridge/staging.hpp"
#include "tywrf/state.hpp"
#include "wrf_physics_bridge.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <type_traits>

namespace {

int failures = 0;

void expect(const bool condition, const std::string_view message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

void expect_capability_bit(
    const std::uint64_t bit,
    const std::string_view label,
    std::uint64_t& seen_bits) {
  expect(bit != 0U, std::string(label) + " capability is nonzero");
  expect((bit & (bit - 1U)) == 0U, std::string(label) + " capability is a single bit");
  expect((seen_bits & bit) == 0U, std::string(label) + " capability does not overlap");
  seen_bits |= bit;
}

tywrf::Grid make_test_grid() {
  return tywrf::Grid({
      .mass_nx = 4,
      .mass_ny = 3,
      .mass_nz = 2,
      .full_nz = 3,
      .halo = {1, 1, 1, 1, 1, 0},
  });
}

TywrfPhysicsDiagnostics make_dirty_diagnostics() {
  return {
      .status = -1,
      .failing_field = -1,
      .validated_field_count = -1,
      .executed_physics = -1,
  };
}

void expect_status(
    const std::int32_t status,
    const TywrfPhysicsDiagnostics& diagnostics,
    const std::int32_t expected_status,
    const std::int32_t expected_field,
    const std::string_view label) {
  expect(status == expected_status, std::string(label) + " return status");
  expect(diagnostics.status == expected_status, std::string(label) + " diagnostic status");
  expect(diagnostics.failing_field == expected_field, std::string(label) + " failing field");
  expect(diagnostics.executed_physics == 0, std::string(label) + " remains a no-op stub");
}

TywrfPhysicsBlockHeader make_v2_header(
    const std::uint32_t struct_size,
    const std::uint64_t capability,
    const TywrfPhysicsBlockHeader* next) {
  return {
      .struct_size = struct_size,
      .abi_version = TYWRF_PHYSICS_ABI_VERSION_V2,
      .capabilities = capability,
      .next = next,
  };
}

}  // namespace

int main() {
  static_assert(TYWRF_PHYSICS_ABI_VERSION == TYWRF_PHYSICS_ABI_VERSION_V1);
  static_assert(TYWRF_PHYSICS_ABI_VERSION_V1 == 1);
  static_assert(TYWRF_PHYSICS_ABI_VERSION_V2 == 2);
  static_assert(TYWRF_PHYSICS_STATUS_STUB_VALIDATED == 0);
  static_assert(TYWRF_PHYSICS_STATUS_WRAPPER_UNAVAILABLE == 10);
  static_assert(TYWRF_PHYSICS_FIELD_RAINNC == 25);
  static_assert(TYWRF_PHYSICS_FIELD_TSK == 26);

  static_assert(std::is_standard_layout_v<TywrfPhysicsBlockHeader>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsBlockHeader>);
  static_assert(offsetof(TywrfPhysicsBlockHeader, struct_size) == 0);
  static_assert(offsetof(TywrfPhysicsBlockHeader, abi_version) == sizeof(std::uint32_t));
  static_assert(offsetof(TywrfPhysicsBlockHeader, capabilities) == 2 * sizeof(std::uint32_t));
  static_assert(offsetof(TywrfPhysicsBlockHeader, next) >=
                offsetof(TywrfPhysicsBlockHeader, capabilities) + sizeof(std::uint64_t));
  static_assert(sizeof(TywrfPhysicsBlockHeader) >=
                offsetof(TywrfPhysicsBlockHeader, next) + sizeof(void*));
  static_assert(alignof(TywrfPhysicsBlockHeader) >= alignof(std::uint64_t));
  static_assert(alignof(TywrfPhysicsBlockHeader) >= alignof(void*));
  static_assert(TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE == sizeof(TywrfPhysicsBlockHeader));

  static_assert(std::is_standard_layout_v<TywrfPhysicsField2D>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsField2D>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsField3D>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsField3D>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsGridMetadata>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsGridMetadata>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsSuiteConfig>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsSuiteConfig>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsStaging>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsStaging>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsDiagnostics>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsDiagnostics>);
  static_assert(sizeof(TywrfPhysicsField2D) == 48);
  static_assert(sizeof(TywrfPhysicsField3D) == 64);
  static_assert(sizeof(TywrfPhysicsGridMetadata) == 72);
  static_assert(sizeof(TywrfPhysicsSuiteConfig) == 44);
  static_assert(sizeof(TywrfPhysicsStaging) == 1584);
  static_assert(sizeof(TywrfPhysicsDiagnostics) == 16);
  static_assert(offsetof(TywrfPhysicsStaging, abi_version) == 0);
  static_assert(offsetof(TywrfPhysicsStaging, grid) == 8);

  static_assert(std::is_standard_layout_v<TywrfPhysicsConstantsV2>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsConstantsV2>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsDriverContextV2>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsDriverContextV2>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsDerivedStateV2>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsDerivedStateV2>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsStaticMaskV2>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsStaticMaskV2>);
  static_assert(std::is_standard_layout_v<TywrfPhysicsSfclaySurfaceV2>);
  static_assert(std::is_trivially_copyable_v<TywrfPhysicsSfclaySurfaceV2>);
  static_assert(offsetof(TywrfPhysicsDriverContextV2, header) == 0);
  static_assert(offsetof(TywrfPhysicsDerivedStateV2, header) == 0);
  static_assert(offsetof(TywrfPhysicsStaticMaskV2, header) == 0);
  static_assert(offsetof(TywrfPhysicsSfclaySurfaceV2, header) == 0);
  static_assert(TYWRF_PHYSICS_CONSTANTS_V2_STRUCT_SIZE == sizeof(TywrfPhysicsConstantsV2));
  static_assert(
      TYWRF_PHYSICS_DRIVER_CONTEXT_V2_STRUCT_SIZE ==
      sizeof(TywrfPhysicsDriverContextV2));
  static_assert(
      TYWRF_PHYSICS_DERIVED_STATE_V2_STRUCT_SIZE == sizeof(TywrfPhysicsDerivedStateV2));
  static_assert(TYWRF_PHYSICS_STATIC_MASK_V2_STRUCT_SIZE == sizeof(TywrfPhysicsStaticMaskV2));
  static_assert(
      TYWRF_PHYSICS_SFCLAY_SURFACE_V2_STRUCT_SIZE ==
      sizeof(TywrfPhysicsSfclaySurfaceV2));

  expect(TYWRF_PHYSICS_CAPABILITY_NONE == 0U, "empty v2 capability set is zero");

  std::uint64_t seen_capabilities = 0;
  expect_capability_bit(
      TYWRF_PHYSICS_CAPABILITY_DRIVER_CONTEXT, "driver_context", seen_capabilities);
  expect_capability_bit(
      TYWRF_PHYSICS_CAPABILITY_DERIVED_STATE, "derived_state", seen_capabilities);
  expect_capability_bit(
      TYWRF_PHYSICS_CAPABILITY_STATIC_MASK, "static_mask", seen_capabilities);
  expect_capability_bit(
      TYWRF_PHYSICS_CAPABILITY_SURFACE_STATE, "surface_state", seen_capabilities);
  expect_capability_bit(TYWRF_PHYSICS_CAPABILITY_SOIL_SNOW, "soil_snow", seen_capabilities);
  expect_capability_bit(TYWRF_PHYSICS_CAPABILITY_TENDENCIES, "tendencies", seen_capabilities);
  expect_capability_bit(
      TYWRF_PHYSICS_CAPABILITY_ACCUMULATORS, "accumulators", seen_capabilities);
  expect_capability_bit(
      TYWRF_PHYSICS_CAPABILITY_RADIATION_RRTMG, "radiation_rrtmg", seen_capabilities);
  expect_capability_bit(TYWRF_PHYSICS_CAPABILITY_SLAB_OCEAN, "slab_ocean", seen_capabilities);
  expect_capability_bit(
      TYWRF_PHYSICS_CAPABILITY_PROCESS_GLOBAL, "process_global", seen_capabilities);

  TywrfPhysicsBlockHeader v2_header{
      .struct_size = TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE,
      .abi_version = TYWRF_PHYSICS_ABI_VERSION_V2,
      .capabilities = TYWRF_PHYSICS_CAPABILITY_DRIVER_CONTEXT |
                      TYWRF_PHYSICS_CAPABILITY_DERIVED_STATE |
                      TYWRF_PHYSICS_CAPABILITY_STATIC_MASK,
      .next = nullptr,
  };
  expect(
      v2_header.struct_size == sizeof(TywrfPhysicsBlockHeader),
      "v2 block header struct_size matches layout");
  expect(v2_header.abi_version == 2U, "v2 block header advertises ABI version 2");
  expect(
      (v2_header.capabilities & TYWRF_PHYSICS_CAPABILITY_DRIVER_CONTEXT) != 0U,
      "v2 block header carries requested capabilities");
  expect(
      (v2_header.capabilities & TYWRF_PHYSICS_CAPABILITY_PROCESS_GLOBAL) == 0U,
      "v2 block header omits absent capabilities");
  expect(v2_header.next == nullptr, "v2 block header next may terminate the sidecar chain");
  expect(
      tywrf_physics_block_header_has_min_size(
          &v2_header, TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE) != 0,
      "v2 block header accepts its minimum size");
  expect(
      tywrf_physics_block_header_has_abi(&v2_header, TYWRF_PHYSICS_ABI_VERSION_V2) != 0,
      "v2 block header ABI predicate accepts v2");

  auto larger_v2_header = v2_header;
  larger_v2_header.struct_size += 16U;
  expect(
      tywrf_physics_block_header_has_min_size(
          &larger_v2_header, TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE) != 0,
      "v2 block header tolerates larger compatible blocks");

  auto short_v2_header = v2_header;
  short_v2_header.struct_size = TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE - 1U;
  expect(
      tywrf_physics_block_header_has_min_size(
          &short_v2_header, TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE) == 0,
      "v2 block header rejects short blocks");

  auto v1_header = v2_header;
  v1_header.abi_version = TYWRF_PHYSICS_ABI_VERSION_V1;
  expect(
      tywrf_physics_block_header_has_abi(&v1_header, TYWRF_PHYSICS_ABI_VERSION_V2) == 0,
      "v2 block header ABI predicate rejects v1");
  expect(
      tywrf_physics_block_header_has_min_size(
          nullptr, TYWRF_PHYSICS_BLOCK_HEADER_STRUCT_SIZE) == 0,
      "v2 block header size predicate rejects null");
  expect(
      tywrf_physics_block_header_has_abi(nullptr, TYWRF_PHYSICS_ABI_VERSION_V2) == 0,
      "v2 block header ABI predicate rejects null");

  const auto grid = make_test_grid();
  tywrf::State<double> state(grid);
  auto state_view = state.view();
  state_view.t(state.t.layout().i_begin(), state.t.layout().j_begin(), state.t.layout().k_begin()) =
      301.25;

  auto d01_staging = tywrf::physics_bridge::make_krosa_staging(
      tywrf::physics_bridge::Domain::d01, grid, state_view, 7, 280.0);
  expect(d01_staging.abi_version == TYWRF_PHYSICS_ABI_VERSION, "ABI version is staged");
  expect(d01_staging.grid.domain_id == TYWRF_PHYSICS_DOMAIN_D01, "d01 domain id is staged");
  expect(d01_staging.grid.dx_m == 10'000.0, "d01 dx is staged");
  expect(d01_staging.grid.dt_s == 40.0, "d01 dt is staged");
  expect(d01_staging.suite.mp_physics == 8, "Thompson option is staged");
  expect(d01_staging.suite.cu_physics == 1, "d01 KF option is staged");
  expect(d01_staging.suite.ra_lw_physics == 4, "RRTMG longwave option is staged");
  expect(d01_staging.suite.ra_sw_physics == 4, "RRTMG shortwave option is staged");
  expect(d01_staging.suite.bl_pbl_physics == 1, "YSU option is staged");
  expect(d01_staging.suite.sf_sfclay_physics == 1, "MM5 surface layer option is staged");
  expect(d01_staging.suite.sf_surface_physics == 2, "Noah LSM option is staged");
  expect(d01_staging.suite.sf_ocean_physics == 1, "slab ocean option is staged");
  expect(d01_staging.suite.isftcflx == 2, "TC flux option is staged");
  expect(d01_staging.t.data == state.t.data(), "state T pointer is staged");
  expect(d01_staging.u.nx - d01_staging.u.halo_i_lower - d01_staging.u.halo_i_upper ==
             grid.config().mass_nx + 1,
         "U active nx is staggered");
  expect(d01_staging.v.ny - d01_staging.v.halo_j_lower - d01_staging.v.halo_j_upper ==
             grid.config().mass_ny + 1,
         "V active ny is staggered");
  expect(d01_staging.w.nz - d01_staging.w.halo_k_lower - d01_staging.w.halo_k_upper ==
             grid.config().full_nz,
         "W active nz uses full levels");

  auto diagnostics = make_dirty_diagnostics();
  auto status = tywrf_wrf_physics_step(&d01_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_STUB_VALIDATED,
      TYWRF_PHYSICS_FIELD_NONE, "valid d01 staging");
  expect(diagnostics.validated_field_count == 25, "valid d01 validates all required fields");
  expect(state_view.t(
             state.t.layout().i_begin(), state.t.layout().j_begin(), state.t.layout().k_begin()) ==
             301.25,
         "stub does not mutate staged state");

  auto d02_staging = tywrf::physics_bridge::make_krosa_staging(
      tywrf::physics_bridge::Domain::d02, grid, state_view, 3, 24.0);
  expect(d02_staging.grid.domain_id == TYWRF_PHYSICS_DOMAIN_D02, "d02 domain id is staged");
  expect(d02_staging.grid.dx_m == 2'000.0, "d02 dx is staged");
  expect(d02_staging.grid.dt_s == 8.0, "d02 dt is staged");
  expect(d02_staging.suite.cu_physics == 0, "d02 explicit convection option is staged");

  diagnostics = make_dirty_diagnostics();
  const auto wrapper_status = tywrf::physics_bridge::run_stub_bridge(d02_staging, &diagnostics);
  expect(wrapper_status == tywrf::physics_bridge::Status::stub_validated,
         "C++ stub wrapper returns typed status");
  expect(diagnostics.validated_field_count == 25, "valid d02 validates all required fields");
  expect(tywrf::physics_bridge::status_name(wrapper_status) == "stub_validated",
         "typed status name is available");
  expect(std::string_view(tywrf_physics_status_name(TYWRF_PHYSICS_STATUS_UNSUPPORTED_SUITE)) ==
             "unsupported_suite",
         "C status name is available");
  expect(std::string_view(tywrf_physics_status_name(TYWRF_PHYSICS_STATUS_WRAPPER_UNAVAILABLE)) ==
             "wrapper_unavailable",
         "v2 wrapper unavailable status name is available");

  tywrf::physics_bridge::SidecarFixtureV2 sidecar(d01_staging);
  expect(
      tywrf::physics_bridge::SidecarFixtureV2::provenance() ==
          "scaffold/finite dummy sidecar; not executed physics",
      "ABI v2 helper declares fixture provenance");
  expect(sidecar.sidecars() == &sidecar.driver_context().header, "sidecar chain starts at driver");
  expect(sidecar.driver_context().constants != nullptr, "sidecar helper stages constants");
  expect(sidecar.driver_context().enable_sfclay == 1, "sidecar helper enables SFCLAY leaf only");
  expect(sidecar.derived_state().u_phy.data != nullptr, "sidecar helper stages derived fields");
  expect(sidecar.static_mask().xlat.data != nullptr, "sidecar helper stages optional XLAT");
  expect(sidecar.sfclay_surface().u10.data != nullptr, "sidecar helper stages finite U10 fixture");

  diagnostics = make_dirty_diagnostics();
  status = tywrf_physics_validate_sidecar_v2(&d01_staging, sidecar.sidecars(), &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_STUB_VALIDATED,
      TYWRF_PHYSICS_FIELD_NONE, "complete ABI v2 sidecar validation");
  expect(diagnostics.validated_field_count == 27, "complete ABI v2 sidecar validates fields");

  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step_ex(&d01_staging, sidecar.sidecars(), &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_WRAPPER_UNAVAILABLE,
      TYWRF_PHYSICS_FIELD_NONE, "ABI v2 extended entrypoint");
  expect(diagnostics.validated_field_count == 27, "ABI v2 extended entrypoint validates fields");

  diagnostics = make_dirty_diagnostics();
  const auto typed_v2_status =
      tywrf::physics_bridge::run_stub_bridge_ex(d01_staging, sidecar.sidecars(), &diagnostics);
  expect(typed_v2_status == tywrf::physics_bridge::Status::wrapper_unavailable,
         "C++ ABI v2 stub wrapper reports unavailable wrapper");
  expect(tywrf::physics_bridge::status_name(typed_v2_status) == "wrapper_unavailable",
         "typed ABI v2 status name is available");

  const auto saved_tsk = sidecar.sfclay_surface().tsk;
  sidecar.sfclay_surface().tsk.data = nullptr;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_physics_validate_sidecar_v2(&d01_staging, sidecar.sidecars(), &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_MISSING_REQUIRED_FIELD,
      TYWRF_PHYSICS_FIELD_TSK, "ABI v2 sidecar missing required TSK");
  sidecar.sfclay_surface().tsk = saved_tsk;

  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(nullptr, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_NULL_ARGUMENT,
      TYWRF_PHYSICS_FIELD_NONE, "null staging");
  expect(diagnostics.validated_field_count == 0, "null staging validates no fields");

  auto bad_staging = d01_staging;
  bad_staging.abi_version = 0;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(&bad_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_UNSUPPORTED_ABI,
      TYWRF_PHYSICS_FIELD_NONE, "unsupported ABI");

  bad_staging = d01_staging;
  bad_staging.grid.domain_id = 99;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(&bad_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_UNSUPPORTED_DOMAIN,
      TYWRF_PHYSICS_FIELD_NONE, "unsupported domain");

  bad_staging = d01_staging;
  bad_staging.suite.mp_physics = 6;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(&bad_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_UNSUPPORTED_SUITE,
      TYWRF_PHYSICS_FIELD_NONE, "unsupported suite");

  bad_staging = d01_staging;
  bad_staging.grid.end_seconds += 1.0;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(&bad_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_INVALID_TIMING,
      TYWRF_PHYSICS_FIELD_NONE, "invalid timing");

  bad_staging = d01_staging;
  bad_staging.u.nx += 1;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(&bad_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_INVALID_DIMENSIONS,
      TYWRF_PHYSICS_FIELD_U, "invalid U dimensions");

  bad_staging = d01_staging;
  bad_staging.qvapor.stride_k += 1;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(&bad_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_INVALID_STRIDES,
      TYWRF_PHYSICS_FIELD_QVAPOR, "invalid QVAPOR strides");

  bad_staging = d01_staging;
  bad_staging.t.element_bytes = 16;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(&bad_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_INVALID_ELEMENT_SIZE,
      TYWRF_PHYSICS_FIELD_T, "invalid T element size");

  bad_staging = d01_staging;
  bad_staging.t.data = nullptr;
  diagnostics = make_dirty_diagnostics();
  status = tywrf_wrf_physics_step(&bad_staging, &diagnostics);
  expect_status(
      status, diagnostics, TYWRF_PHYSICS_STATUS_MISSING_REQUIRED_FIELD,
      TYWRF_PHYSICS_FIELD_T, "missing T field");
  expect(diagnostics.validated_field_count == 5, "missing T reports prior validated fields");

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated WRF physics bridge ABI/staging stub status handling\n";
  return 0;
}
