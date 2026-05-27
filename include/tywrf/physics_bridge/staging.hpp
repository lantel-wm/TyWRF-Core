#pragma once

#include "tywrf/grid.hpp"
#include "tywrf/state.hpp"
#include "wrf_physics_bridge.h"

#include <cstdint>
#include <string_view>
#include <type_traits>

namespace tywrf::physics_bridge {

enum class Domain : std::int32_t {
  d01 = TYWRF_PHYSICS_DOMAIN_D01,
  d02 = TYWRF_PHYSICS_DOMAIN_D02,
};

enum class Status : std::int32_t {
  stub_validated = TYWRF_PHYSICS_STATUS_STUB_VALIDATED,
  null_argument = TYWRF_PHYSICS_STATUS_NULL_ARGUMENT,
  unsupported_abi = TYWRF_PHYSICS_STATUS_UNSUPPORTED_ABI,
  unsupported_domain = TYWRF_PHYSICS_STATUS_UNSUPPORTED_DOMAIN,
  unsupported_suite = TYWRF_PHYSICS_STATUS_UNSUPPORTED_SUITE,
  invalid_timing = TYWRF_PHYSICS_STATUS_INVALID_TIMING,
  invalid_dimensions = TYWRF_PHYSICS_STATUS_INVALID_DIMENSIONS,
  invalid_strides = TYWRF_PHYSICS_STATUS_INVALID_STRIDES,
  invalid_element_size = TYWRF_PHYSICS_STATUS_INVALID_ELEMENT_SIZE,
  missing_required_field = TYWRF_PHYSICS_STATUS_MISSING_REQUIRED_FIELD,
  wrapper_unavailable = TYWRF_PHYSICS_STATUS_WRAPPER_UNAVAILABLE,
};

[[nodiscard]] constexpr TywrfPhysicsSuiteConfig make_krosa_suite_config(
    const Domain domain) noexcept {
  return {
      8,
      domain == Domain::d01 ? 1 : 0,
      4,
      4,
      1,
      1,
      2,
      1,
      2,
      8,
      4,
  };
}

[[nodiscard]] constexpr std::int32_t krosa_grid_spacing_m(const Domain domain) noexcept {
  return domain == Domain::d01 ? 10'000 : 2'000;
}

[[nodiscard]] constexpr std::int32_t krosa_time_step_seconds(const Domain domain) noexcept {
  return domain == Domain::d01 ? 40 : 8;
}

template <typename Real>
[[nodiscard]] constexpr TywrfPhysicsField2D make_physics_field(
    const FieldView2D<Real> view) noexcept {
  static_assert(!std::is_const_v<Real>, "physics staging requires mutable fields");
  static_assert(
      std::is_same_v<Real, float> || std::is_same_v<Real, double>,
      "physics staging supports float or double fields");

  return {
      static_cast<void*>(view.data),
      view.nx,
      view.ny,
      view.stride_i,
      view.stride_j,
      view.halo.i_lower,
      view.halo.i_upper,
      view.halo.j_lower,
      view.halo.j_upper,
      static_cast<std::int32_t>(sizeof(Real)),
  };
}

template <typename Real>
[[nodiscard]] constexpr TywrfPhysicsField3D make_physics_field(
    const FieldView3D<Real> view) noexcept {
  static_assert(!std::is_const_v<Real>, "physics staging requires mutable fields");
  static_assert(
      std::is_same_v<Real, float> || std::is_same_v<Real, double>,
      "physics staging supports float or double fields");

  return {
      static_cast<void*>(view.data),
      view.nx,
      view.ny,
      view.nz,
      view.stride_i,
      view.stride_k,
      view.stride_j,
      view.halo.i_lower,
      view.halo.i_upper,
      view.halo.j_lower,
      view.halo.j_upper,
      view.halo.k_lower,
      view.halo.k_upper,
      static_cast<std::int32_t>(sizeof(Real)),
  };
}

template <typename Real>
[[nodiscard]] TywrfPhysicsStaging make_krosa_staging(
    const Domain domain,
    const Grid& grid,
    const StateView<Real> state,
    const std::int64_t step_index,
    const double start_seconds) noexcept {
  static_assert(!std::is_const_v<Real>, "physics staging requires mutable state");

  const auto mass_shape = grid.mass_shape();
  const auto full_shape = grid.w_shape();
  const auto dt = static_cast<double>(krosa_time_step_seconds(domain));
  const auto dx = static_cast<double>(krosa_grid_spacing_m(domain));

  return {
      TYWRF_PHYSICS_ABI_VERSION,
      {
          static_cast<std::int32_t>(domain),
          mass_shape.nx,
          mass_shape.ny,
          mass_shape.nz,
          full_shape.nz,
          dx,
          dx,
          dt,
          step_index,
          start_seconds,
          start_seconds + dt,
      },
      make_krosa_suite_config(domain),
      make_physics_field(state.u),
      make_physics_field(state.v),
      make_physics_field(state.w),
      make_physics_field(state.ph),
      make_physics_field(state.phb),
      make_physics_field(state.t),
      make_physics_field(state.p),
      make_physics_field(state.pb),
      make_physics_field(state.qvapor),
      make_physics_field(state.qcloud),
      make_physics_field(state.qrain),
      make_physics_field(state.qice),
      make_physics_field(state.qsnow),
      make_physics_field(state.qgraup),
      make_physics_field(state.qnice),
      make_physics_field(state.qnrain),
      make_physics_field(state.mu),
      make_physics_field(state.mub),
      make_physics_field(state.psfc),
      make_physics_field(state.u10),
      make_physics_field(state.v10),
      make_physics_field(state.t2),
      make_physics_field(state.q2),
      make_physics_field(state.rainc),
      make_physics_field(state.rainnc),
  };
}

[[nodiscard]] Status validate_staging(
    const TywrfPhysicsStaging& staging,
    TywrfPhysicsDiagnostics* diagnostics = nullptr) noexcept;

[[nodiscard]] Status validate_sidecar_v2(
    const TywrfPhysicsStaging& staging,
    const TywrfPhysicsBlockHeader* sidecars,
    TywrfPhysicsDiagnostics* diagnostics = nullptr) noexcept;

[[nodiscard]] Status run_stub_bridge(
    const TywrfPhysicsStaging& staging,
    TywrfPhysicsDiagnostics* diagnostics = nullptr) noexcept;

[[nodiscard]] Status run_stub_bridge_ex(
    const TywrfPhysicsStaging& staging,
    const TywrfPhysicsBlockHeader* sidecars,
    TywrfPhysicsDiagnostics* diagnostics = nullptr) noexcept;

[[nodiscard]] std::string_view status_name(Status status) noexcept;

}  // namespace tywrf::physics_bridge
