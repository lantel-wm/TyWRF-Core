#include "tywrf/dynamics/mass_continuity_tendency.hpp"

#include "tywrf/field_view.hpp"
#include "tywrf/state.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
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

bool nearly_equal(const double lhs, const double rhs) {
  return std::abs(lhs - rhs) < 1.0e-12;
}

constexpr tywrf::ActiveShape2D mu_shape() {
  return {4, 3};
}

constexpr tywrf::Halo2D halo() {
  return {1, 1, 1, 1};
}

tywrf::FieldStorage2D<double> make_mu_storage() {
  return tywrf::FieldStorage2D<double>(mu_shape(), halo());
}

tywrf::FieldStorage2D<double> make_full_u_flux_storage() {
  return tywrf::FieldStorage2D<double>({mu_shape().nx + 1, mu_shape().ny}, halo());
}

tywrf::FieldStorage2D<double> make_full_v_flux_storage() {
  return tywrf::FieldStorage2D<double>({mu_shape().nx, mu_shape().ny + 1}, halo());
}

tywrf::FieldStorage2D<double> make_high_boundary_partial_u_flux_storage() {
  return tywrf::FieldStorage2D<double>(mu_shape(), halo());
}

tywrf::FieldStorage2D<double> make_high_boundary_partial_v_flux_storage() {
  return tywrf::FieldStorage2D<double>(mu_shape(), halo());
}

template <typename Storage>
void fill_storage(Storage& storage, const double value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

template <typename Storage>
void fill_u_flux_gradient(Storage& storage) {
  const auto view = storage.view();
  const auto active_nx = view.nx - view.halo.i_lower - view.halo.i_upper;
  const auto active_ny = view.ny - view.halo.j_lower - view.halo.j_upper;
  for (std::int32_t active_j = 0; active_j < active_ny; ++active_j) {
    const auto j = view.halo.j_lower + active_j;
    for (std::int32_t active_i = 0; active_i < active_nx; ++active_i) {
      const auto i = view.halo.i_lower + active_i;
      view(i, j) = 10.0 + 3.0 * static_cast<double>(active_i) +
                   0.25 * static_cast<double>(active_j);
    }
  }
}

template <typename Storage>
void fill_v_flux_gradient(Storage& storage) {
  const auto view = storage.view();
  const auto active_nx = view.nx - view.halo.i_lower - view.halo.i_upper;
  const auto active_ny = view.ny - view.halo.j_lower - view.halo.j_upper;
  for (std::int32_t active_j = 0; active_j < active_ny; ++active_j) {
    const auto j = view.halo.j_lower + active_j;
    for (std::int32_t active_i = 0; active_i < active_nx; ++active_i) {
      const auto i = view.halo.i_lower + active_i;
      view(i, j) = 20.0 + 5.0 * static_cast<double>(active_j) +
                   0.125 * static_cast<double>(active_i);
    }
  }
}

[[nodiscard]] std::int64_t active_points(const tywrf::FieldLayout2D layout) {
  return static_cast<std::int64_t>(layout.active_nx()) *
         static_cast<std::int64_t>(layout.active_ny());
}

[[nodiscard]] bool active(
    const tywrf::FieldLayout2D layout,
    const std::int32_t i,
    const std::int32_t j) {
  return i >= layout.i_begin() && i < layout.i_end() && j >= layout.j_begin() &&
         j < layout.j_end();
}

[[nodiscard]] bool high_boundary_partial_cell(
    const tywrf::FieldLayout2D layout,
    const std::int32_t i,
    const std::int32_t j) {
  if (!active(layout, i, j)) {
    return false;
  }
  const auto active_i = i - layout.i_begin();
  const auto active_j = j - layout.j_begin();
  return active_i == layout.active_nx() - 1 ||
         active_j == layout.active_ny() - 1;
}

template <typename Real>
[[nodiscard]] tywrf::dynamics::MassContinuityTendencyViews<Real> views(
    tywrf::FieldStorage2D<Real>& mu_target,
    const tywrf::FieldStorage2D<Real>& u_face_mass_flux,
    const tywrf::FieldStorage2D<Real>& v_face_mass_flux) {
  return {
      mu_target.view(),
      u_face_mass_flux.view(),
      v_face_mass_flux.view(),
  };
}

[[nodiscard]] tywrf::dynamics::MassContinuityTendencyConfig<double>
enabled_config() {
  return {
      .dt_seconds = 2.0,
      .dx_m = 3.0,
      .dy_m = 5.0,
      .enable_mass_continuity = true,
      .diagnostic_only = true,
      .gate_candidate = false,
      .validation_gate_evidence = false,
  };
}

void test_pod_contract_and_disabled_noop_counts() {
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::MassContinuityTendencyConfig<double>>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::MassContinuityTendencyConfig<double>>);
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::MassContinuityTendencyViews<double>>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::MassContinuityTendencyViews<double>>);
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::MassContinuityTendencyReport>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::MassContinuityTendencyReport>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::MassContinuityTendencyStatus>);

  auto mu = make_mu_storage();
  auto u_flux = make_full_u_flux_storage();
  auto v_flux = make_full_v_flux_storage();
  fill_storage(mu, 100.0);
  fill_u_flux_gradient(u_flux);
  fill_v_flux_gradient(v_flux);

  const auto report =
      tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
          views(mu, u_flux, v_flux),
          {.dt_seconds = 2.0, .dx_m = 3.0, .dy_m = 5.0});

  expect(report.status == tywrf::dynamics::MassContinuityTendencyStatus::ok,
         "disabled status is ok with valid views");
  expect(report.active_mu_points == active_points(mu.layout()),
         "disabled active MU count");
  expect(report.active_points == report.active_mu_points,
         "disabled active total count");
  expect(report.active_u_flux_points == active_points(u_flux.layout()),
         "disabled active U-flux count");
  expect(report.active_v_flux_points == active_points(v_flux.layout()),
         "disabled active V-flux count");
  expect(report.updated_mu_points == 0, "disabled update count");
  expect(report.skipped_mu_points == report.active_mu_points,
         "disabled skipped count");
  expect(report.skipped_disabled_points == report.active_mu_points,
         "disabled skipped-disabled count");
  expect(!report.mass_continuity_enabled, "disabled flag retained");
  expect(!report.halo_writes, "disabled has no halo writes");

  for (std::size_t idx = 0; idx < mu.size(); ++idx) {
    expect(nearly_equal(mu.data()[idx], 100.0), "disabled leaves MU unchanged");
  }
}

void test_simple_divergence_update() {
  auto mu = make_mu_storage();
  auto u_flux = make_full_u_flux_storage();
  auto v_flux = make_full_v_flux_storage();
  fill_storage(mu, 100.0);
  fill_u_flux_gradient(u_flux);
  fill_v_flux_gradient(v_flux);

  const auto report =
      tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
          views(mu, u_flux, v_flux), enabled_config());

  expect(report.status == tywrf::dynamics::MassContinuityTendencyStatus::ok,
         "divergence update status");
  expect(report.updated_mu_points == active_points(mu.layout()),
         "divergence updates all cells with full face fluxes");
  expect(report.skipped_mu_points == 0, "divergence skips no full-face cells");
  expect(report.skipped_boundary_points == 0,
         "divergence has no boundary skips with full faces");

  const auto view = mu.view();
  for (std::int32_t j = 0; j < view.ny; ++j) {
    for (std::int32_t i = 0; i < view.nx; ++i) {
      const auto expected = active(mu.layout(), i, j) ? 96.0 : 100.0;
      expect(nearly_equal(view(i, j), expected),
             active(mu.layout(), i, j) ? "active MU divergence update"
                                       : "MU halo preserved");
    }
  }
}

void test_zero_divergence_preserves_mu() {
  auto mu = make_mu_storage();
  auto u_flux = make_full_u_flux_storage();
  auto v_flux = make_full_v_flux_storage();
  fill_storage(mu, 42.0);
  fill_storage(u_flux, 11.0);
  fill_storage(v_flux, 13.0);

  const auto report =
      tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
          views(mu, u_flux, v_flux), enabled_config());

  expect(report.status == tywrf::dynamics::MassContinuityTendencyStatus::ok,
         "zero divergence status");
  expect(report.updated_mu_points == active_points(mu.layout()),
         "zero divergence update count");
  for (std::size_t idx = 0; idx < mu.size(); ++idx) {
    expect(nearly_equal(mu.data()[idx], 42.0), "zero divergence preserves MU");
  }
}

void test_boundary_skip_counts_for_partial_face_fluxes() {
  auto mu = make_mu_storage();
  auto u_flux = make_high_boundary_partial_u_flux_storage();
  auto v_flux = make_high_boundary_partial_v_flux_storage();
  fill_storage(mu, 100.0);
  fill_u_flux_gradient(u_flux);
  fill_v_flux_gradient(v_flux);

  const auto report =
      tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
          views(mu, u_flux, v_flux), enabled_config());
  const auto expected_updates =
      static_cast<std::int64_t>(mu.layout().active_nx() - 1) *
      static_cast<std::int64_t>(mu.layout().active_ny() - 1);
  const auto expected_skips = active_points(mu.layout()) - expected_updates;

  expect(report.status == tywrf::dynamics::MassContinuityTendencyStatus::ok,
         "partial boundary flux status");
  expect(report.updated_mu_points == expected_updates,
         "partial boundary update count");
  expect(report.skipped_mu_points == expected_skips,
         "partial boundary skipped count");
  expect(report.skipped_boundary_points == expected_skips,
         "partial boundary skipped-boundary count");

  const auto view = mu.view();
  for (std::int32_t j = 0; j < view.ny; ++j) {
    for (std::int32_t i = 0; i < view.nx; ++i) {
      const auto expected = active(mu.layout(), i, j) &&
                                    !high_boundary_partial_cell(mu.layout(), i, j)
                                ? 96.0
                                : 100.0;
      expect(nearly_equal(view(i, j), expected),
             high_boundary_partial_cell(mu.layout(), i, j)
                 ? "partial boundary MU skipped"
                 : "partial boundary MU update or halo preserve");
    }
  }
}

void test_invalid_layout_and_config_fail_closed() {
  auto mu = make_mu_storage();
  auto u_flux = make_full_u_flux_storage();
  auto v_flux = make_full_v_flux_storage();
  fill_storage(mu, 31.0);
  fill_u_flux_gradient(u_flux);
  fill_v_flux_gradient(v_flux);

  auto invalid_views = views(mu, u_flux, v_flux);
  invalid_views.mu_target.data = nullptr;
  auto report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      invalid_views, enabled_config());
  expect(report.status == tywrf::dynamics::MassContinuityTendencyStatus::null_target,
         "null MU target status");

  invalid_views = views(mu, u_flux, v_flux);
  invalid_views.u_face_mass_flux.data = nullptr;
  report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      invalid_views, enabled_config());
  expect(report.status == tywrf::dynamics::MassContinuityTendencyStatus::null_flux,
         "null flux status");

  invalid_views = views(mu, u_flux, v_flux);
  invalid_views.mu_target.stride_i = 2;
  report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      invalid_views, enabled_config());
  expect(report.status ==
             tywrf::dynamics::MassContinuityTendencyStatus::invalid_layout,
         "invalid stride status");

  auto config = enabled_config();
  config.dx_m = 0.0;
  report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      views(mu, u_flux, v_flux), config);
  expect(report.status ==
             tywrf::dynamics::MassContinuityTendencyStatus::invalid_config,
         "non-positive dx status");

  config = enabled_config();
  config.dt_seconds = -1.0;
  report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      views(mu, u_flux, v_flux), config);
  expect(report.status ==
             tywrf::dynamics::MassContinuityTendencyStatus::invalid_config,
         "negative dt status");

  config = enabled_config();
  config.gate_candidate = true;
  report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      views(mu, u_flux, v_flux), config);
  expect(report.status ==
             tywrf::dynamics::MassContinuityTendencyStatus::invalid_config,
         "diagnostic gate-candidate contradiction status");

  for (std::size_t idx = 0; idx < mu.size(); ++idx) {
    expect(nearly_equal(mu.data()[idx], 31.0),
           "invalid contract leaves MU unchanged");
  }
}

void test_nonfinite_values_fail_closed_without_writes() {
  auto mu = make_mu_storage();
  auto u_flux = make_full_u_flux_storage();
  auto v_flux = make_full_v_flux_storage();
  fill_storage(mu, 50.0);
  fill_u_flux_gradient(u_flux);
  fill_v_flux_gradient(v_flux);

  auto u_view = u_flux.view();
  u_view(u_view.halo.i_lower + 1, u_view.halo.j_lower) =
      std::numeric_limits<double>::quiet_NaN();
  auto report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      views(mu, u_flux, v_flux), enabled_config());
  expect(report.status ==
             tywrf::dynamics::MassContinuityTendencyStatus::invalid_flux_value,
         "nonfinite flux status");
  expect(report.updated_mu_points == 0, "nonfinite flux update count");
  for (std::size_t idx = 0; idx < mu.size(); ++idx) {
    expect(nearly_equal(mu.data()[idx], 50.0), "nonfinite flux leaves MU");
  }

  fill_u_flux_gradient(u_flux);
  auto mu_view = mu.view();
  mu_view(mu_view.halo.i_lower, mu_view.halo.j_lower) =
      std::numeric_limits<double>::infinity();
  report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      views(mu, u_flux, v_flux), enabled_config());
  expect(report.status ==
             tywrf::dynamics::MassContinuityTendencyStatus::invalid_mu_value,
         "nonfinite MU status");
  expect(report.updated_mu_points == 0, "nonfinite MU update count");
  expect(std::isinf(mu_view(mu_view.halo.i_lower, mu_view.halo.j_lower)),
         "nonfinite MU value not overwritten");
  for (std::int32_t j = 0; j < mu_view.ny; ++j) {
    for (std::int32_t i = 0; i < mu_view.nx; ++i) {
      if (i == mu_view.halo.i_lower && j == mu_view.halo.j_lower) {
        continue;
      }
      expect(nearly_equal(mu_view(i, j), 50.0), "nonfinite MU leaves others");
    }
  }
}

void test_shape_mismatch_fails_closed() {
  auto mu = make_mu_storage();
  auto u_flux = tywrf::FieldStorage2D<double>(
      tywrf::ActiveShape2D{mu_shape().nx + 2, mu_shape().ny}, halo());
  auto v_flux = make_full_v_flux_storage();
  fill_storage(mu, 80.0);
  fill_storage(u_flux, 1.0);
  fill_storage(v_flux, 1.0);

  auto report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      views(mu, u_flux, v_flux), enabled_config());
  expect(report.status ==
             tywrf::dynamics::MassContinuityTendencyStatus::
                 incompatible_u_flux_shape,
         "U flux shape mismatch status");

  auto good_u_flux = make_full_u_flux_storage();
  auto bad_v_flux = tywrf::FieldStorage2D<double>(
      tywrf::ActiveShape2D{mu_shape().nx + 1, mu_shape().ny + 1}, halo());
  fill_storage(good_u_flux, 1.0);
  fill_storage(bad_v_flux, 1.0);
  report = tywrf::dynamics::apply_dry_air_mass_continuity_mu_tendency(
      views(mu, good_u_flux, bad_v_flux), enabled_config());
  expect(report.status ==
             tywrf::dynamics::MassContinuityTendencyStatus::
                 incompatible_v_flux_shape,
         "V flux shape mismatch status");

  for (std::size_t idx = 0; idx < mu.size(); ++idx) {
    expect(nearly_equal(mu.data()[idx], 80.0), "shape mismatch leaves MU");
  }
}

}  // namespace

int main() {
  test_pod_contract_and_disabled_noop_counts();
  test_simple_divergence_update();
  test_zero_divergence_preserves_mu();
  test_boundary_skip_counts_for_partial_face_fluxes();
  test_invalid_layout_and_config_fail_closed();
  test_nonfinite_values_fail_closed_without_writes();
  test_shape_mismatch_fails_closed();

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated standalone dry-air mass continuity MU tendency\n";
  return 0;
}
