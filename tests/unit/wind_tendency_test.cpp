#include "tywrf/dynamics/wind_tendency.hpp"

#include "tywrf/grid.hpp"
#include "tywrf/state.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
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

tywrf::Grid make_grid() {
  return tywrf::Grid({
      .mass_nx = 4,
      .mass_ny = 3,
      .mass_nz = 2,
      .full_nz = 3,
      .halo = {1, 1, 1, 1, 0, 0},
  });
}

template <typename Storage>
void fill_storage(Storage& storage, const double value) {
  for (std::size_t idx = 0; idx < storage.size(); ++idx) {
    storage.data()[idx] = value;
  }
}

template <typename Storage>
void fill_by_index(Storage& storage) {
  for (std::size_t idx = 0; idx < storage.size(); ++idx) {
    storage.data()[idx] = static_cast<double>(idx) + 0.25;
  }
}

template <typename Storage>
void fill_x_gradient(Storage& storage, const double intercept, const double slope) {
  const auto view = storage.view();
  for (std::int32_t j = 0; j < view.ny; ++j) {
    for (std::int32_t k = 0; k < view.nz; ++k) {
      for (std::int32_t i = 0; i < view.nx; ++i) {
        view(i, j, k) = intercept + slope * static_cast<double>(i);
      }
    }
  }
}

template <typename Storage>
void fill_y_gradient(Storage& storage, const double intercept, const double slope) {
  const auto view = storage.view();
  for (std::int32_t j = 0; j < view.ny; ++j) {
    for (std::int32_t k = 0; k < view.nz; ++k) {
      for (std::int32_t i = 0; i < view.nx; ++i) {
        view(i, j, k) = intercept + slope * static_cast<double>(j);
      }
    }
  }
}

template <typename Storage>
void fill_quadratic_xy(Storage& storage) {
  const auto view = storage.view();
  for (std::int32_t j = 0; j < view.ny; ++j) {
    for (std::int32_t k = 0; k < view.nz; ++k) {
      for (std::int32_t i = 0; i < view.nx; ++i) {
        view(i, j, k) =
            static_cast<double>(i * i) + 10.0 * static_cast<double>(j) +
            0.5 * static_cast<double>(k);
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

template <typename Real>
[[nodiscard]] tywrf::dynamics::WindTendencyViews<Real> views(
    tywrf::FieldStorage3D<Real>& u_target,
    tywrf::FieldStorage3D<Real>& v_target,
    const tywrf::FieldStorage3D<Real>& u_source,
    const tywrf::FieldStorage3D<Real>& v_source,
    const tywrf::FieldStorage3D<Real>& u_wind_x,
    const tywrf::FieldStorage3D<Real>& u_wind_y,
    const tywrf::FieldStorage3D<Real>& v_wind_x,
    const tywrf::FieldStorage3D<Real>& v_wind_y) {
  return {
      {u_target.view(), u_source.view(), u_wind_x.view(), u_wind_y.view()},
      {v_target.view(), v_source.view(), v_wind_x.view(), v_wind_y.view()}};
}

void test_reports_pod_flags_and_active_counts() {
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::WindTendencyConfig<double>>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::WindTendencyConfig<double>>);
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::WindTendencyAdvectionForm>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::WindTendencyAdvectionForm>);
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::WindTendencyFieldViews<double>>);
  static_assert(
      std::is_trivially_copyable_v<
          tywrf::dynamics::WindTendencyFieldViews<double>>);
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::WindTendencyViews<double>>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::WindTendencyViews<double>>);
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::WindTendencyReport>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::WindTendencyReport>);

  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> u_source(grid.u_layout());
  tywrf::FieldStorage3D<double> v_source(grid.v_layout());
  tywrf::FieldStorage3D<double> u_wind_x(grid.u_layout());
  tywrf::FieldStorage3D<double> u_wind_y(grid.u_layout());
  tywrf::FieldStorage3D<double> v_wind_x(grid.v_layout());
  tywrf::FieldStorage3D<double> v_wind_y(grid.v_layout());
  fill_storage(u_target, 1.0);
  fill_storage(v_target, 2.0);
  fill_storage(u_source, 1.0);
  fill_storage(v_source, 2.0);
  fill_storage(u_wind_x, 0.0);
  fill_storage(u_wind_y, 0.0);
  fill_storage(v_wind_x, 0.0);
  fill_storage(v_wind_y, 0.0);

  const auto report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y, v_wind_x,
            v_wind_y),
      tywrf::dynamics::WindTendencyConfig<double>{
          .dt_seconds = 4.0,
          .dx_m = 2.0,
          .dy_m = 3.0,
          .enable_horizontal_advection = true,
          .diagnostic_only = true,
          .gate_candidate = true,
          .validation_gate_evidence = true,
      });

  expect(report.status == tywrf::dynamics::WindTendencyStatus::ok, "report status");
  expect(report.active_u_points == active_points(grid.u_layout()), "U active count");
  expect(report.active_v_points == active_points(grid.v_layout()), "V active count");
  expect(report.active_points == report.active_u_points + report.active_v_points,
         "total active count");
  expect(report.updated_u_points == report.active_u_points, "U updated count");
  expect(report.updated_v_points == report.active_v_points, "V updated count");
  expect(!report.halo_writes, "report marks no halo writes");
  expect(report.diagnostic_only, "diagnostic flag retained");
  expect(report.gate_candidate, "gate candidate flag retained");
  expect(!report.validation_gate_evidence, "core report is not validation evidence");
  expect(report.horizontal_advection_enabled, "advection flag retained");
}

void test_active_only_writes_and_halo_preservation() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> u_source(grid.u_layout());
  tywrf::FieldStorage3D<double> v_source(grid.v_layout());
  tywrf::FieldStorage3D<double> u_wind_x(grid.u_layout());
  tywrf::FieldStorage3D<double> u_wind_y(grid.u_layout());
  tywrf::FieldStorage3D<double> v_wind_x(grid.v_layout());
  tywrf::FieldStorage3D<double> v_wind_y(grid.v_layout());
  fill_storage(u_target, 100.0);
  fill_storage(v_target, 200.0);
  fill_x_gradient(u_source, 0.0, 2.0);
  fill_y_gradient(v_source, 0.0, 3.0);
  fill_storage(u_wind_x, 3.0);
  fill_storage(u_wind_y, 0.0);
  fill_storage(v_wind_x, 0.0);
  fill_storage(v_wind_y, 2.0);

  const auto report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y, v_wind_x,
            v_wind_y),
      {.dt_seconds = 4.0, .dx_m = 2.0, .dy_m = 3.0});

  expect(report.status == tywrf::dynamics::WindTendencyStatus::ok,
         "active-only status");
  const auto u_view = u_target.view();
  const auto v_view = v_target.view();
  for (std::int32_t j = 0; j < u_view.ny; ++j) {
    for (std::int32_t k = 0; k < u_view.nz; ++k) {
      for (std::int32_t i = 0; i < u_view.nx; ++i) {
        const auto expected = active(grid.u_layout(), i, j, k) ? 88.0 : 100.0;
        expect(nearly_equal(u_view(i, j, k), expected),
               active(grid.u_layout(), i, j, k) ? "U active write"
                                                 : "U halo preserved");
      }
    }
  }
  for (std::int32_t j = 0; j < v_view.ny; ++j) {
    for (std::int32_t k = 0; k < v_view.nz; ++k) {
      for (std::int32_t i = 0; i < v_view.nx; ++i) {
        const auto expected = active(grid.v_layout(), i, j, k) ? 192.0 : 200.0;
        expect(nearly_equal(v_view(i, j, k), expected),
               active(grid.v_layout(), i, j, k) ? "V active write"
                                                 : "V halo preserved");
      }
    }
  }
}

void test_i_contiguous_expectation() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  const auto u_view = u_target.view();
  const auto v_view = v_target.view();

  expect(u_view.stride_i == 1, "U stride_i is contiguous");
  expect(v_view.stride_i == 1, "V stride_i is contiguous");
  expect(u_view.index(
             grid.u_layout().i_begin() + 1,
             grid.u_layout().j_begin(),
             grid.u_layout().k_begin()) ==
             u_view.index(
                 grid.u_layout().i_begin(),
                 grid.u_layout().j_begin(),
                 grid.u_layout().k_begin()) +
                 1,
         "U i index increments by one");
  expect(v_view.index(
             grid.v_layout().i_begin() + 1,
             grid.v_layout().j_begin(),
             grid.v_layout().k_begin()) ==
             v_view.index(
                 grid.v_layout().i_begin(),
                 grid.v_layout().j_begin(),
                 grid.v_layout().k_begin()) +
                 1,
         "V i index increments by one");
}

void test_zero_gradient_and_zero_wind_noop() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> u_source(grid.u_layout());
  tywrf::FieldStorage3D<double> v_source(grid.v_layout());
  tywrf::FieldStorage3D<double> u_wind_x(grid.u_layout());
  tywrf::FieldStorage3D<double> u_wind_y(grid.u_layout());
  tywrf::FieldStorage3D<double> v_wind_x(grid.v_layout());
  tywrf::FieldStorage3D<double> v_wind_y(grid.v_layout());
  fill_by_index(u_target);
  fill_by_index(v_target);
  fill_storage(u_source, 9.0);
  fill_storage(v_source, 8.0);
  fill_storage(u_wind_x, 10.0);
  fill_storage(u_wind_y, -5.0);
  fill_storage(v_wind_x, 7.0);
  fill_storage(v_wind_y, -11.0);
  const auto u_before = u_target;
  const auto v_before = v_target;

  auto report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y, v_wind_x,
            v_wind_y),
      {.dt_seconds = 4.0, .dx_m = 2.0, .dy_m = 3.0});
  expect(report.status == tywrf::dynamics::WindTendencyStatus::ok,
         "zero-gradient status");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], u_before.data()[idx]),
           "zero-gradient leaves U unchanged");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], v_before.data()[idx]),
           "zero-gradient leaves V unchanged");
  }

  fill_x_gradient(u_source, 0.0, 2.0);
  fill_y_gradient(v_source, 0.0, -3.0);
  fill_storage(u_wind_x, 0.0);
  fill_storage(u_wind_y, 0.0);
  fill_storage(v_wind_x, 0.0);
  fill_storage(v_wind_y, 0.0);
  report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y, v_wind_x,
            v_wind_y),
      {.dt_seconds = 4.0, .dx_m = 2.0, .dy_m = 3.0});
  expect(report.status == tywrf::dynamics::WindTendencyStatus::ok, "zero-wind status");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], u_before.data()[idx]),
           "zero-wind leaves U unchanged");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], v_before.data()[idx]),
           "zero-wind leaves V unchanged");
  }
}

void test_simple_gradient_sign_and_magnitude() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> u_source(grid.u_layout());
  tywrf::FieldStorage3D<double> v_source(grid.v_layout());
  tywrf::FieldStorage3D<double> u_wind_x(grid.u_layout());
  tywrf::FieldStorage3D<double> u_wind_y(grid.u_layout());
  tywrf::FieldStorage3D<double> v_wind_x(grid.v_layout());
  tywrf::FieldStorage3D<double> v_wind_y(grid.v_layout());
  fill_storage(u_target, 50.0);
  fill_storage(v_target, 70.0);
  fill_x_gradient(u_source, 5.0, 4.0);
  fill_y_gradient(v_source, 9.0, -6.0);
  fill_storage(u_wind_x, 1.5);
  fill_storage(u_wind_y, 0.0);
  fill_storage(v_wind_x, 0.0);
  fill_storage(v_wind_y, -2.0);

  const auto report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y, v_wind_x,
            v_wind_y),
      {.dt_seconds = 5.0, .dx_m = 2.0, .dy_m = 3.0});

  expect(report.status == tywrf::dynamics::WindTendencyStatus::ok,
         "simple-gradient status");
  expect(nearly_equal(
             u_target.view()(
                 grid.u_layout().i_begin(),
                 grid.u_layout().j_begin(),
                 grid.u_layout().k_begin()),
             35.0),
         "positive x wind against positive U gradient has negative tendency");
  expect(nearly_equal(
             v_target.view()(
                 grid.v_layout().i_begin(),
                 grid.v_layout().j_begin(),
                 grid.v_layout().k_begin()),
             50.0),
         "negative y wind with negative V gradient has negative tendency");
}

void test_centered_default_equivalence() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> default_u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> default_v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> u_source(grid.u_layout());
  tywrf::FieldStorage3D<double> v_source(grid.v_layout());
  tywrf::FieldStorage3D<double> u_wind_x(grid.u_layout());
  tywrf::FieldStorage3D<double> u_wind_y(grid.u_layout());
  tywrf::FieldStorage3D<double> v_wind_x(grid.v_layout());
  tywrf::FieldStorage3D<double> v_wind_y(grid.v_layout());
  fill_storage(default_u_target, 50.0);
  fill_storage(default_v_target, 70.0);
  fill_quadratic_xy(u_source);
  fill_quadratic_xy(v_source);
  fill_storage(u_wind_x, -1.25);
  fill_storage(u_wind_y, 0.75);
  fill_storage(v_wind_x, 1.5);
  fill_storage(v_wind_y, -0.5);

  auto centered_u_target = default_u_target;
  auto centered_v_target = default_v_target;

  const auto default_report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(default_u_target, default_v_target, u_source, v_source, u_wind_x,
            u_wind_y, v_wind_x, v_wind_y),
      {.dt_seconds = 2.0, .dx_m = 2.0, .dy_m = 4.0});
  const auto centered_report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(centered_u_target, centered_v_target, u_source, v_source, u_wind_x,
            u_wind_y, v_wind_x, v_wind_y),
      {.dt_seconds = 2.0,
       .dx_m = 2.0,
       .dy_m = 4.0,
       .advection_form = tywrf::dynamics::WindTendencyAdvectionForm::centered});

  expect(default_report.status == tywrf::dynamics::WindTendencyStatus::ok,
         "default centered status");
  expect(centered_report.status == tywrf::dynamics::WindTendencyStatus::ok,
         "explicit centered status");
  for (std::size_t idx = 0; idx < default_u_target.size(); ++idx) {
    expect(nearly_equal(default_u_target.data()[idx], centered_u_target.data()[idx]),
           "default matches explicit centered U");
  }
  for (std::size_t idx = 0; idx < default_v_target.size(); ++idx) {
    expect(nearly_equal(default_v_target.data()[idx], centered_v_target.data()[idx]),
           "default matches explicit centered V");
  }
}

void test_upwind_sign_dependent_behavior() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> upwind_u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> upwind_v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> centered_u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> centered_v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> u_source(grid.u_layout());
  tywrf::FieldStorage3D<double> v_source(grid.v_layout());
  tywrf::FieldStorage3D<double> u_wind_x(grid.u_layout());
  tywrf::FieldStorage3D<double> u_wind_y(grid.u_layout());
  tywrf::FieldStorage3D<double> v_wind_x(grid.v_layout());
  tywrf::FieldStorage3D<double> v_wind_y(grid.v_layout());
  fill_storage(upwind_u_target, 100.0);
  fill_storage(upwind_v_target, 200.0);
  fill_storage(centered_u_target, 100.0);
  fill_storage(centered_v_target, 200.0);
  fill_quadratic_xy(u_source);
  fill_quadratic_xy(v_source);
  fill_storage(u_wind_x, 0.0);
  fill_storage(u_wind_y, 0.0);
  fill_storage(v_wind_x, 0.0);
  fill_storage(v_wind_y, 0.0);

  const auto ui = grid.u_layout().i_begin() + 2;
  const auto uj = grid.u_layout().j_begin() + 1;
  const auto uk = grid.u_layout().k_begin();
  u_wind_x.view()(ui, uj, uk) = 2.0;
  u_wind_x.view()(ui - 1, uj, uk) = -3.0;

  const auto vi = grid.v_layout().i_begin() + 1;
  const auto vj = grid.v_layout().j_begin() + 2;
  const auto vk = grid.v_layout().k_begin();
  v_wind_y.view()(vi, vj, vk) = -4.0;
  v_wind_y.view()(vi, vj - 1, vk) = 1.5;

  const auto centered_report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(centered_u_target, centered_v_target, u_source, v_source, u_wind_x,
            u_wind_y, v_wind_x, v_wind_y),
      {.dt_seconds = 1.0,
       .dx_m = 1.0,
       .dy_m = 1.0,
       .advection_form = tywrf::dynamics::WindTendencyAdvectionForm::centered});
  const auto upwind_report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(upwind_u_target, upwind_v_target, u_source, v_source, u_wind_x,
            u_wind_y, v_wind_x, v_wind_y),
      {.dt_seconds = 1.0,
       .dx_m = 1.0,
       .dy_m = 1.0,
       .advection_form = tywrf::dynamics::WindTendencyAdvectionForm::upwind});

  expect(centered_report.status == tywrf::dynamics::WindTendencyStatus::ok,
         "centered comparison status");
  expect(upwind_report.status == tywrf::dynamics::WindTendencyStatus::ok,
         "upwind status");
  const auto u_q = u_source.view()(ui, uj, uk);
  const auto expected_u = 100.0 - (2.0 * u_q - (-3.0 * u_q));
  expect(nearly_equal(upwind_u_target.view()(ui, uj, uk), expected_u),
         "upwind U uses positive right and negative left donors");
  const auto v_q_up = v_source.view()(vi, vj + 1, vk);
  const auto v_q_down = v_source.view()(vi, vj - 1, vk);
  const auto expected_v = 200.0 - ((-4.0 * v_q_up) - (1.5 * v_q_down));
  expect(nearly_equal(upwind_v_target.view()(vi, vj, vk), expected_v),
         "upwind V uses negative up and positive down donors");
  expect(
      !nearly_equal(
          upwind_u_target.view()(ui, uj, uk),
          centered_u_target.view()(ui, uj, uk)),
      "upwind U sensitivity differs from centered");
  expect(
      !nearly_equal(
          upwind_v_target.view()(vi, vj, vk),
          centered_v_target.view()(vi, vj, vk)),
      "upwind V sensitivity differs from centered");
}

void test_invalid_contracts_fail_closed() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> u_target(grid.u_layout());
  tywrf::FieldStorage3D<double> v_target(grid.v_layout());
  tywrf::FieldStorage3D<double> u_source(grid.u_layout());
  tywrf::FieldStorage3D<double> v_source(grid.v_layout());
  tywrf::FieldStorage3D<double> bad_u_source(grid.mass_layout());
  tywrf::FieldStorage3D<double> u_wind_x(grid.u_layout());
  tywrf::FieldStorage3D<double> u_wind_y(grid.u_layout());
  tywrf::FieldStorage3D<double> v_wind_x(grid.v_layout());
  tywrf::FieldStorage3D<double> v_wind_y(grid.v_layout());
  fill_storage(u_target, 11.0);
  fill_storage(v_target, 12.0);
  fill_x_gradient(u_source, 0.0, 2.0);
  fill_y_gradient(v_source, 0.0, 3.0);
  fill_storage(bad_u_source, 1.0);
  fill_storage(u_wind_x, 3.0);
  fill_storage(u_wind_y, 0.0);
  fill_storage(v_wind_x, 0.0);
  fill_storage(v_wind_y, 2.0);

  auto bad_views =
      views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y, v_wind_x,
            v_wind_y);
  const auto& const_bad_u_source = bad_u_source;
  bad_views.u.source = const_bad_u_source.view();
  auto report = tywrf::dynamics::apply_horizontal_wind_tendency(
      bad_views,
      {.dt_seconds = 4.0, .dx_m = 2.0, .dy_m = 3.0});
  expect(report.status == tywrf::dynamics::WindTendencyStatus::mismatched_source_layout,
         "mismatched source status");
  expect(report.active_points == 0, "mismatch active count zero");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], 11.0), "mismatch leaves U unchanged");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], 12.0), "mismatch leaves V unchanged");
  }

  bad_views = views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y,
                    v_wind_x, v_wind_y);
  bad_views.v.target.data = nullptr;
  report = tywrf::dynamics::apply_horizontal_wind_tendency(
      bad_views,
      {.dt_seconds = 4.0, .dx_m = 2.0, .dy_m = 3.0});
  expect(report.status == tywrf::dynamics::WindTendencyStatus::null_target,
         "null target status");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], 11.0), "null target leaves U unchanged");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], 12.0), "null target leaves V unchanged");
  }

  report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y, v_wind_x,
            v_wind_y),
      {.dt_seconds = 4.0, .dx_m = 0.0, .dy_m = 3.0});
  expect(report.status == tywrf::dynamics::WindTendencyStatus::invalid_config,
         "invalid config status");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], 11.0),
           "invalid config leaves U unchanged");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], 12.0),
           "invalid config leaves V unchanged");
  }

  report = tywrf::dynamics::apply_horizontal_wind_tendency(
      views(u_target, v_target, u_source, v_source, u_wind_x, u_wind_y, v_wind_x,
            v_wind_y),
      {.dt_seconds = 4.0,
       .dx_m = 2.0,
       .dy_m = 3.0,
       .advection_form =
           static_cast<tywrf::dynamics::WindTendencyAdvectionForm>(99)});
  expect(report.status == tywrf::dynamics::WindTendencyStatus::invalid_config,
         "invalid advection form status");
  for (std::size_t idx = 0; idx < u_target.size(); ++idx) {
    expect(nearly_equal(u_target.data()[idx], 11.0),
           "invalid advection form leaves U unchanged");
  }
  for (std::size_t idx = 0; idx < v_target.size(); ++idx) {
    expect(nearly_equal(v_target.data()[idx], 12.0),
           "invalid advection form leaves V unchanged");
  }
}

}  // namespace

int main() {
  test_reports_pod_flags_and_active_counts();
  test_active_only_writes_and_halo_preservation();
  test_i_contiguous_expectation();
  test_zero_gradient_and_zero_wind_noop();
  test_simple_gradient_sign_and_magnitude();
  test_centered_default_equivalence();
  test_upwind_sign_dependent_behavior();
  test_invalid_contracts_fail_closed();

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated CUDA-ready U/V wind tendency skeleton\n";
  return 0;
}
