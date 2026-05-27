#include "tywrf/grid.hpp"
#include "tywrf/physics_bridge/staging.hpp"
#include "tywrf/state.hpp"
#include "wrf_physics_bridge.h"

#include <cassert>
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

}  // namespace

int main() {
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
