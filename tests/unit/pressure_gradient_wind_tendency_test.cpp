#include "tywrf/dynamics/pressure_gradient_wind_tendency.hpp"

#include "tywrf/grid.hpp"
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

tywrf::Grid make_grid(const tywrf::Halo3D halo = {1, 1, 1, 1, 0, 0}) {
  return tywrf::Grid({
      .mass_nx = 3,
      .mass_ny = 2,
      .mass_nz = 2,
      .full_nz = 3,
      .halo = halo,
  });
}

template <typename Storage>
void fill_storage(Storage& storage, const double value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

template <typename Storage>
void fill_pressure_gradient(Storage& perturbation, Storage& base) {
  const auto p_view = perturbation.view();
  const auto pb_view = base.view();
  for (std::int32_t j = 0; j < p_view.ny; ++j) {
    const auto active_j = j - p_view.halo.j_lower;
    for (std::int32_t k = 0; k < p_view.nz; ++k) {
      const auto active_k = k - p_view.halo.k_lower;
      for (std::int32_t i = 0; i < p_view.nx; ++i) {
        const auto active_i = i - p_view.halo.i_lower;
        p_view(i, j, k) =
            8.0 * static_cast<double>(active_i) +
            9.0 * static_cast<double>(active_j) +
            0.25 * static_cast<double>(active_k);
        pb_view(i, j, k) = 1'000.0 + 0.5 * static_cast<double>(active_k);
      }
    }
  }
}

[[nodiscard]] bool active(
    const tywrf::FieldLayout3D layout,
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k) {
  return i >= layout.i_begin() && i < layout.i_end() && j >= layout.j_begin() &&
         j < layout.j_end() && k >= layout.k_begin() && k < layout.k_end();
}

[[nodiscard]] std::int64_t active_points(const tywrf::FieldLayout3D layout) {
  return static_cast<std::int64_t>(layout.active_nx()) *
         static_cast<std::int64_t>(layout.active_ny()) *
         static_cast<std::int64_t>(layout.active_nz());
}

[[nodiscard]] std::int64_t interior_u_points(const tywrf::Grid& grid) {
  const auto mass = grid.mass_layout();
  return static_cast<std::int64_t>(mass.active_nx() - 1) *
         static_cast<std::int64_t>(mass.active_ny()) *
         static_cast<std::int64_t>(mass.active_nz());
}

[[nodiscard]] std::int64_t interior_v_points(const tywrf::Grid& grid) {
  const auto mass = grid.mass_layout();
  return static_cast<std::int64_t>(mass.active_nx()) *
         static_cast<std::int64_t>(mass.active_ny() - 1) *
         static_cast<std::int64_t>(mass.active_nz());
}

[[nodiscard]] bool internal_u_face(
    const tywrf::FieldLayout3D layout,
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k) {
  if (!active(layout, i, j, k)) {
    return false;
  }
  const auto active_i = i - layout.i_begin();
  return active_i > 0 && active_i < layout.active_nx() - 1;
}

[[nodiscard]] bool internal_v_face(
    const tywrf::FieldLayout3D layout,
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k) {
  if (!active(layout, i, j, k)) {
    return false;
  }
  const auto active_j = j - layout.j_begin();
  return active_j > 0 && active_j < layout.active_ny() - 1;
}

template <typename Real>
[[nodiscard]] tywrf::dynamics::PressureGradientWindTendencyViews<Real> views(
    tywrf::FieldStorage3D<Real>& u_target,
    tywrf::FieldStorage3D<Real>& v_target,
    const tywrf::FieldStorage3D<Real>& perturbation_pressure,
    const tywrf::FieldStorage3D<Real>& base_pressure) {
  return {
      u_target.view(),
      v_target.view(),
      perturbation_pressure.view(),
      base_pressure.view(),
  };
}

[[nodiscard]] tywrf::dynamics::PressureGradientWindTendencyConfig<double>
enabled_config() {
  return {
      .dt_seconds = 2.0,
      .dx_m = 4.0,
      .dy_m = 3.0,
      .constant_specific_volume_m3_per_kg = 0.5,
      .pressure_units = tywrf::dynamics::PressureGradientPressureUnits::pascal,
      .specific_volume_units =
          tywrf::dynamics::PressureGradientSpecificVolumeUnits::
              cubic_meter_per_kilogram,
      .enable_pressure_gradient = true,
      .diagnostic_only = true,
      .gate_candidate = false,
      .validation_gate_evidence = false,
  };
}

void test_pod_contract_and_disabled_counts() {
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::PressureGradientWindTendencyConfig<double>>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureGradientWindTendencyConfig<double>>);
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::PressureGradientWindTendencyViews<double>>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureGradientWindTendencyViews<double>>);
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::PressureGradientWindTendencyReport>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureGradientWindTendencyReport>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureGradientWindTendencyStatus>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureGradientWindTendencyForm>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureGradientPressureUnits>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureGradientSpecificVolumeUnits>);

  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> p(grid.mass_layout());
  tywrf::FieldStorage3D<double> pb(grid.mass_layout());
  fill_storage(u_target, 10.0);
  fill_storage(v_target, 20.0);
  fill_pressure_gradient(p, pb);

  const auto report =
      tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
          views(u_target, v_target, p, pb),
          {.dt_seconds = 2.0, .dx_m = 4.0, .dy_m = 3.0});

  expect(report.status == tywrf::dynamics::PressureGradientWindTendencyStatus::ok,
         "disabled status is ok with valid views");
  expect(report.active_u_points == active_points(grid.u_layout()),
         "disabled U active count");
  expect(report.active_v_points == active_points(grid.v_layout()),
         "disabled V active count");
  expect(report.active_points == report.active_u_points + report.active_v_points,
         "disabled total active count");
  expect(report.updated_u_points == 0, "disabled U update count");
  expect(report.updated_v_points == 0, "disabled V update count");
  expect(report.skipped_u_points == 0, "disabled U skip count");
  expect(report.skipped_v_points == 0, "disabled V skip count");
  expect(!report.pressure_gradient_enabled, "disabled flag retained");
  expect(report.used_base_pressure, "base pressure validated");
  expect(!report.halo_writes, "disabled report has no halo writes");
}

void test_active_only_writes_and_known_gradient() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> p(grid.mass_layout());
  tywrf::FieldStorage3D<double> pb(grid.mass_layout());
  fill_storage(u_target, 100.0);
  fill_storage(v_target, 200.0);
  fill_pressure_gradient(p, pb);

  const auto report =
      tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
          views(u_target, v_target, p, pb), enabled_config());

  expect(report.status == tywrf::dynamics::PressureGradientWindTendencyStatus::ok,
         "known-gradient status");
  expect(report.updated_u_points == interior_u_points(grid),
         "known-gradient interior U update count");
  expect(report.updated_v_points == interior_v_points(grid),
         "known-gradient interior V update count");
  expect(report.skipped_u_points ==
             active_points(grid.u_layout()) - interior_u_points(grid),
         "known-gradient skipped U boundary count");
  expect(report.skipped_v_points ==
             active_points(grid.v_layout()) - interior_v_points(grid),
         "known-gradient skipped V boundary count");

  const auto u_view = u_target.view();
  for (std::int32_t j = 0; j < u_view.ny; ++j) {
    for (std::int32_t k = 0; k < u_view.nz; ++k) {
      for (std::int32_t i = 0; i < u_view.nx; ++i) {
        const auto expected =
            internal_u_face(grid.u_layout(), i, j, k) ? 98.0 : 100.0;
        expect(nearly_equal(u_view(i, j, k), expected),
               internal_u_face(grid.u_layout(), i, j, k)
                   ? "U interior gradient update"
                   : "U boundary or halo preserved");
      }
    }
  }

  const auto v_view = v_target.view();
  for (std::int32_t j = 0; j < v_view.ny; ++j) {
    for (std::int32_t k = 0; k < v_view.nz; ++k) {
      for (std::int32_t i = 0; i < v_view.nx; ++i) {
        const auto expected =
            internal_v_face(grid.v_layout(), i, j, k) ? 197.0 : 200.0;
        expect(nearly_equal(v_view(i, j, k), expected),
               internal_v_face(grid.v_layout(), i, j, k)
                   ? "V interior gradient update"
                   : "V boundary or halo preserved");
      }
    }
  }
}

void test_zero_halo_grid_updates_internal_faces() {
  const auto grid = make_grid({0, 0, 0, 0, 0, 0});
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> p(grid.mass_layout());
  tywrf::FieldStorage3D<double> pb(grid.mass_layout());
  fill_storage(u_target, 10.0);
  fill_storage(v_target, 20.0);
  fill_pressure_gradient(p, pb);

  const auto report =
      tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
          views(u_target, v_target, p, pb), enabled_config());

  expect(report.status == tywrf::dynamics::PressureGradientWindTendencyStatus::ok,
         "zero-halo pressure layout status");
  expect(report.updated_u_points == interior_u_points(grid),
         "zero-halo interior U update count");
  expect(report.updated_v_points == interior_v_points(grid),
         "zero-halo interior V update count");
  expect(report.skipped_u_points ==
             active_points(grid.u_layout()) - interior_u_points(grid),
         "zero-halo skipped U boundary count");
  expect(report.skipped_v_points ==
             active_points(grid.v_layout()) - interior_v_points(grid),
         "zero-halo skipped V boundary count");
}

void test_zero_dt_is_noop_after_validation() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> p(grid.mass_layout());
  tywrf::FieldStorage3D<double> pb(grid.mass_layout());
  fill_storage(u_target, 10.0);
  fill_storage(v_target, 20.0);
  fill_pressure_gradient(p, pb);

  auto config = enabled_config();
  config.dt_seconds = 0.0;
  const auto report =
      tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
          views(u_target, v_target, p, pb), config);

  expect(report.status == tywrf::dynamics::PressureGradientWindTendencyStatus::ok,
         "zero dt status");
  expect(report.updated_u_points == 0, "zero dt U update count");
  expect(report.updated_v_points == 0, "zero dt V update count");
  expect(report.skipped_u_points == 0, "zero dt U skip count");
  expect(report.skipped_v_points == 0, "zero dt V skip count");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], 10.0), "zero dt leaves U");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], 20.0), "zero dt leaves V");
  }
}

void test_invalid_pressure_values_fail_closed() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> p(grid.mass_layout());
  tywrf::FieldStorage3D<double> pb(grid.mass_layout());
  fill_storage(u_target, 31.0);
  fill_storage(v_target, 32.0);
  fill_pressure_gradient(p, pb);

  auto p_view = p.view();
  const auto pb_view = pb.view();
  const auto invalid_i = p_view.halo.i_lower;
  const auto invalid_j = p_view.halo.j_lower;
  const auto invalid_k = p_view.halo.k_lower;
  p_view(invalid_i, invalid_j, invalid_k) =
      -pb_view(invalid_i, invalid_j, invalid_k);

  auto report =
      tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
          views(u_target, v_target, p, pb), enabled_config());
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::
                 invalid_pressure_value,
         "non-positive pressure status");
  expect(report.updated_u_points == 0, "invalid pressure U update count");
  expect(report.updated_v_points == 0, "invalid pressure V update count");
  expect(report.skipped_u_points == 0, "invalid pressure U skip count");
  expect(report.skipped_v_points == 0, "invalid pressure V skip count");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], 31.0),
           "non-positive pressure leaves U");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], 32.0),
           "non-positive pressure leaves V");
  }

  fill_pressure_gradient(p, pb);
  p_view = p.view();
  p_view(invalid_i + 1, invalid_j, invalid_k) =
      std::numeric_limits<double>::quiet_NaN();
  report = tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
      views(u_target, v_target, p, pb), enabled_config());
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::
                 invalid_pressure_value,
         "non-finite pressure status");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], 31.0),
           "non-finite pressure leaves U");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], 32.0),
           "non-finite pressure leaves V");
  }
}

void test_invalid_contracts_fail_closed() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> p(grid.mass_layout());
  tywrf::FieldStorage3D<double> pb(grid.mass_layout());
  fill_storage(u_target, 11.0);
  fill_storage(v_target, 12.0);
  fill_pressure_gradient(p, pb);

  auto invalid_views = views(u_target, v_target, p, pb);
  invalid_views.u_target.data = nullptr;
  auto report =
      tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
          invalid_views, enabled_config());
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::null_target,
         "null target status");

  invalid_views = views(u_target, v_target, p, pb);
  invalid_views.base_pressure.data = nullptr;
  report = tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
      invalid_views, enabled_config());
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::null_pressure,
         "null pressure status");

  invalid_views = views(u_target, v_target, p, pb);
  invalid_views.v_target.stride_i = 2;
  report = tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
      invalid_views, enabled_config());
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::invalid_layout,
         "invalid target layout status");

  tywrf::FieldStorage3D<double> mismatched_pb(grid.w_layout());
  fill_storage(mismatched_pb, 1.0);
  report = tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
      views(u_target, v_target, p, mismatched_pb), enabled_config());
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::
                 mismatched_pressure_layout,
         "mismatched pressure layout status");

  tywrf::FieldStorage3D<double> bad_u(grid.mass_layout());
  fill_storage(bad_u, 13.0);
  report = tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
      views(bad_u, v_target, p, pb), enabled_config());
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::
                 mismatched_target_layout,
         "mismatched target layout status");

  auto config = enabled_config();
  config.dx_m = 0.0;
  report = tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
      views(u_target, v_target, p, pb), config);
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::invalid_config,
         "non-positive dx status");

  config = enabled_config();
  config.pressure_units =
      tywrf::dynamics::PressureGradientPressureUnits::unspecified;
  report = tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
      views(u_target, v_target, p, pb), config);
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::unsupported_config,
         "missing pressure units status");

  config = enabled_config();
  config.form = tywrf::dynamics::PressureGradientWindTendencyForm::wrf_exact;
  report = tywrf::dynamics::apply_horizontal_pressure_gradient_wind_tendency(
      views(u_target, v_target, p, pb), config);
  expect(report.status ==
             tywrf::dynamics::PressureGradientWindTendencyStatus::not_implemented,
         "WRF-exact form not implemented status");

  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], 11.0), "invalid contracts leave U");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], 12.0), "invalid contracts leave V");
  }
}

}  // namespace

int main() {
  test_pod_contract_and_disabled_counts();
  test_active_only_writes_and_known_gradient();
  test_zero_halo_grid_updates_internal_faces();
  test_zero_dt_is_noop_after_validation();
  test_invalid_pressure_values_fail_closed();
  test_invalid_contracts_fail_closed();

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated standalone pressure-gradient wind tendency skeleton\n";
  return 0;
}
