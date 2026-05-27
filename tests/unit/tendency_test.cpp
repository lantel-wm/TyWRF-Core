#include "tywrf/dynamics/tendency.hpp"

#include "tywrf/grid.hpp"
#include "tywrf/state.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string_view>

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
      .halo = {1, 2, 1, 1, 1, 0},
  });
}

template <typename Storage>
void fill_storage(Storage& storage, const double value) {
  for (std::size_t idx = 0; idx < storage.size(); ++idx) {
    storage.data()[idx] = value;
  }
}

void test_apply_3d_updates_only_active_region() {
  const auto layout = make_grid().mass_layout();
  tywrf::FieldStorage3D<double> field(layout);
  tywrf::FieldStorage3D<double> tendency(layout);
  fill_storage(field, -999.0);
  fill_storage(tendency, 2.0);

  auto field_view = field.view();
  const auto& const_tendency = tendency;
  const auto report = tywrf::dynamics::apply_tendency(
      field_view,
      const_tendency.view(),
      3.0);

  expect(report.status == tywrf::dynamics::TendencyApplyStatus::ok, "3D apply status");
  expect(report.active_points == 24, "3D active count");

  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t k = 0; k < layout.nz; ++k) {
      for (std::int32_t i = 0; i < layout.nx; ++i) {
        const bool active = i >= layout.i_begin() && i < layout.i_end() &&
                            j >= layout.j_begin() && j < layout.j_end() &&
                            k >= layout.k_begin() && k < layout.k_end();
        const auto value = field_view(i, j, k);
        expect(nearly_equal(value, active ? -993.0 : -999.0),
               active ? "3D active value changed" : "3D halo value preserved");
      }
    }
  }

  expect(field_view.index(layout.i_begin() + 1, layout.j_begin(), layout.k_begin()) ==
             field_view.index(layout.i_begin(), layout.j_begin(), layout.k_begin()) + 1,
         "3D i dimension is contiguous");
}

void test_apply_2d_updates_only_active_region() {
  const auto layout = make_grid().surface_layout();
  tywrf::FieldStorage2D<double> field(layout);
  tywrf::FieldStorage2D<double> tendency(layout);
  fill_storage(field, 10.0);
  fill_storage(tendency, -0.5);

  auto field_view = field.view();
  const auto& const_tendency = tendency;
  const auto report = tywrf::dynamics::apply_tendency(
      field_view,
      const_tendency.view(),
      4.0);

  expect(report.status == tywrf::dynamics::TendencyApplyStatus::ok, "2D apply status");
  expect(report.active_points == 12, "2D active count");

  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t i = 0; i < layout.nx; ++i) {
      const bool active = i >= layout.i_begin() && i < layout.i_end() &&
                          j >= layout.j_begin() && j < layout.j_end();
      const auto value = field_view(i, j);
      expect(nearly_equal(value, active ? 8.0 : 10.0),
             active ? "2D active value changed" : "2D halo value preserved");
    }
  }

  expect(field_view.index(layout.i_begin() + 1, layout.j_begin()) ==
             field_view.index(layout.i_begin(), layout.j_begin()) + 1,
         "2D i dimension is contiguous");
}

void test_zero_tendency_is_noop_but_counts_active_region() {
  const auto layout = make_grid().mass_layout();
  tywrf::FieldStorage3D<double> field(layout);
  for (std::size_t idx = 0; idx < field.size(); ++idx) {
    field.data()[idx] = static_cast<double>(idx);
  }

  const auto report = tywrf::dynamics::apply_zero_tendency(field.view());
  expect(report.status == tywrf::dynamics::TendencyApplyStatus::ok, "zero apply status");
  expect(report.active_points == 24, "zero active count");

  for (std::size_t idx = 0; idx < field.size(); ++idx) {
    expect(nearly_equal(field.data()[idx], static_cast<double>(idx)),
           "zero tendency leaves storage unchanged");
  }
}

void test_state_apply_uses_existing_state_views() {
  const auto grid = make_grid();
  tywrf::State<double> state(grid);
  tywrf::State<double> tendency(grid);
  fill_storage(state.t, 300.0);
  fill_storage(tendency.t, 0.25);
  fill_storage(state.mu, 1.0);
  fill_storage(tendency.mu, 2.0);

  const auto report = tywrf::dynamics::apply_state_tendencies(
      state.view(),
      static_cast<const tywrf::State<double>&>(tendency).view(),
      8.0);

  const auto expected_active_points =
      5 * 3 * 2 + 4 * 4 * 2 + 3 * 4 * 3 * 3 + 11 * 4 * 3 * 2 + 9 * 4 * 3;
  expect(report.status == tywrf::dynamics::TendencyApplyStatus::ok, "state apply status");
  expect(report.active_points == expected_active_points, "state active count");

  const auto mass = grid.mass_layout();
  const auto surface = grid.surface_layout();
  expect(nearly_equal(state.t.view()(mass.i_begin(), mass.j_begin(), mass.k_begin()), 302.0),
         "state 3D field updated");
  expect(nearly_equal(state.mu.view()(surface.i_begin(), surface.j_begin()), 17.0),
         "state 2D field updated");
  expect(nearly_equal(state.t.view()(0, 0, 0), 300.0), "state 3D halo preserved");
  expect(nearly_equal(state.mu.view()(0, 0), 1.0), "state 2D halo preserved");
}

void test_mismatched_layout_fails_without_writing() {
  const auto grid = make_grid();
  tywrf::FieldStorage3D<double> field(grid.mass_layout());
  tywrf::FieldStorage3D<double> tendency(grid.u_layout());
  fill_storage(field, 5.0);
  fill_storage(tendency, 1.0);

  auto field_view = field.view();
  const auto& const_tendency = tendency;
  const auto report = tywrf::dynamics::apply_tendency(
      field_view,
      const_tendency.view(),
      10.0);

  expect(report.status == tywrf::dynamics::TendencyApplyStatus::mismatched_layout,
         "mismatched layout status");
  expect(report.active_points == 0, "mismatched layout count");
  for (std::size_t idx = 0; idx < field.size(); ++idx) {
    expect(nearly_equal(field.data()[idx], 5.0), "mismatched layout does not write");
  }
}

void test_state_apply_validates_before_writing() {
  const auto grid = make_grid();
  tywrf::State<double> state(grid);
  tywrf::State<double> tendency(grid);
  fill_storage(state.u, 4.0);
  fill_storage(state.t, 300.0);
  fill_storage(tendency.u, 1.0);
  fill_storage(tendency.t, 1.0);

  auto bad_tendency = static_cast<const tywrf::State<double>&>(tendency).view();
  bad_tendency.qnrain = bad_tendency.u;
  const auto report = tywrf::dynamics::apply_state_tendencies(
      state.view(),
      bad_tendency,
      5.0);

  const auto u_layout = grid.u_layout();
  const auto mass = grid.mass_layout();
  expect(report.status == tywrf::dynamics::TendencyApplyStatus::mismatched_layout,
         "bad state tendency status");
  expect(report.active_points == 0, "bad state tendency count");
  expect(nearly_equal(state.u.view()(u_layout.i_begin(), u_layout.j_begin(), u_layout.k_begin()),
                      4.0),
         "bad state tendency leaves early 3D field unchanged");
  expect(nearly_equal(state.t.view()(mass.i_begin(), mass.j_begin(), mass.k_begin()), 300.0),
         "bad state tendency leaves later 3D field unchanged");
}

}  // namespace

int main() {
  test_apply_3d_updates_only_active_region();
  test_apply_2d_updates_only_active_region();
  test_zero_tendency_is_noop_but_counts_active_region();
  test_state_apply_uses_existing_state_views();
  test_mismatched_layout_fails_without_writing();
  test_state_apply_validates_before_writing();

  if (failures != 0) {
    return 1;
  }
  std::cout << "Validated CUDA-ready tendency apply skeleton\n";
  return 0;
}
