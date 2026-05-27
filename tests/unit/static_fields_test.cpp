#include "tywrf/nest/nest_interface.hpp"
#include "tywrf/nest/static_fields.hpp"
#include "tywrf/state.hpp"

#include <cassert>
#include <cmath>
#include <string_view>

namespace {

constexpr float kTolerance = 1.0e-4F;

[[nodiscard]] tywrf::nest::ParentChildDescriptor make_descriptor() {
  return {
      {1, 10'000, 8, 8, 9, 9},
      {2, 2'000, 10, 10, 11, 11},
      5,
      5,
  };
}

[[nodiscard]] tywrf::FieldStorage2D<float> make_field(const int nx, const int ny) {
  return tywrf::FieldStorage2D<float>(
      tywrf::ActiveShape2D{nx, ny}, tywrf::uniform_halo_2d(0));
}

void fill_linear(
    tywrf::FieldStorage2D<float>& field,
    const float base,
    const float ax,
    const float ay) {
  auto view = field.view();
  const auto layout = field.layout();
  for (int j = 0; j < layout.active_ny(); ++j) {
    for (int i = 0; i < layout.active_nx(); ++i) {
      view(i, j) = base + ax * static_cast<float>(i) + ay * static_cast<float>(j);
    }
  }
}

[[nodiscard]] float linear(
    const float base,
    const float ax,
    const float ay,
    const double x,
    const double y) {
  return base + ax * static_cast<float>(x) + ay * static_cast<float>(y);
}

void expect_close(const float actual, const float expected, const std::string_view label) {
  if (std::fabs(actual - expected) > kTolerance) {
    (void)label;
    assert(false);
  }
}

void run_forward_refresh_case() {
  const auto descriptor = make_descriptor();
  const auto from = tywrf::nest::make_nest_pose(
      descriptor, {2, 2, tywrf::nest::IndexBase::one_based});
  const auto to = tywrf::nest::make_nest_pose(
      descriptor, {3, 2, tywrf::nest::IndexBase::one_based});
  const auto plan = tywrf::nest::build_remap_plan(from, to);
  assert(plan.ok());

  auto xlat_start = make_field(10, 10);
  auto xlong_start = make_field(10, 10);
  auto hgt_start = make_field(10, 10);
  auto parent_hgt = make_field(8, 8);
  auto xlat_out = make_field(10, 10);
  auto xlong_out = make_field(10, 10);
  auto hgt_out = make_field(10, 10);
  fill_linear(xlat_start, 100.0F, 1.0F, 10.0F);
  fill_linear(xlong_start, 200.0F, 2.0F, 5.0F);
  fill_linear(hgt_start, 1'000.0F, 3.0F, 7.0F);
  fill_linear(parent_hgt, 500.0F, 11.0F, 13.0F);

  const auto& const_xlat_start = xlat_start;
  const auto& const_xlong_start = xlong_start;
  const auto& const_hgt_start = hgt_start;
  const auto& const_parent_hgt = parent_hgt;
  const auto report = tywrf::nest::refresh_moving_nest_static_fields(
      descriptor,
      {3, 2, tywrf::nest::IndexBase::one_based},
      plan,
      const_xlat_start.view(),
      const_xlong_start.view(),
      const_hgt_start.view(),
      const_parent_hgt.view(),
      xlat_out.view(),
      xlong_out.view(),
      hgt_out.view());
  assert(report.ok());
  assert(!report.uses_reference_end);
  assert(report.overlap_cell_count == 50);
  assert(report.exposed_cell_count == 50);
  assert(report.coordinate_extrapolated_cell_count == 50);
  assert(report.parent_hgt_interpolated_cell_count == 50);

  expect_close(xlat_out.view()(0, 2), linear(100.0F, 1.0F, 10.0F, 5.0, 2.0), "XLAT overlap");
  expect_close(xlong_out.view()(0, 2), linear(200.0F, 2.0F, 5.0F, 5.0, 2.0), "XLONG overlap");
  expect_close(hgt_out.view()(0, 2), linear(1'000.0F, 3.0F, 7.0F, 5.0, 2.0), "HGT overlap");

  const auto exposed_xlat = xlat_out.view()(9, 0);
  const auto stale_xlat = xlat_start.view()(9, 0);
  assert(std::isfinite(exposed_xlat));
  assert(std::fabs(exposed_xlat - stale_xlat) > kTolerance);
  expect_close(exposed_xlat, linear(100.0F, 1.0F, 10.0F, 14.0, 0.0), "XLAT exposed");
  expect_close(
      xlong_out.view()(9, 0),
      linear(200.0F, 2.0F, 5.0F, 14.0, 0.0),
      "XLONG exposed");
  expect_close(
      hgt_out.view()(9, 0),
      linear(500.0F, 11.0F, 13.0F, 3.8, 1.0),
      "HGT parent exposed");
}

void run_reverse_refresh_case() {
  const auto descriptor = make_descriptor();
  const auto from = tywrf::nest::make_nest_pose(
      descriptor, {3, 3, tywrf::nest::IndexBase::one_based});
  const auto to = tywrf::nest::make_nest_pose(
      descriptor, {2, 2, tywrf::nest::IndexBase::one_based});
  const auto plan = tywrf::nest::build_remap_plan(from, to);
  assert(plan.ok());

  auto xlat_start = make_field(10, 10);
  auto xlong_start = make_field(10, 10);
  auto hgt_start = make_field(10, 10);
  auto parent_hgt = make_field(8, 8);
  auto xlat_out = make_field(10, 10);
  auto xlong_out = make_field(10, 10);
  auto hgt_out = make_field(10, 10);
  fill_linear(xlat_start, 100.0F, 1.0F, 10.0F);
  fill_linear(xlong_start, 200.0F, 2.0F, 5.0F);
  fill_linear(hgt_start, 1'000.0F, 3.0F, 7.0F);
  fill_linear(parent_hgt, 500.0F, 11.0F, 13.0F);

  const auto& const_xlat_start = xlat_start;
  const auto& const_xlong_start = xlong_start;
  const auto& const_hgt_start = hgt_start;
  const auto& const_parent_hgt = parent_hgt;
  const auto report = tywrf::nest::refresh_moving_nest_static_fields(
      descriptor,
      {2, 2, tywrf::nest::IndexBase::one_based},
      plan,
      const_xlat_start.view(),
      const_xlong_start.view(),
      const_hgt_start.view(),
      const_parent_hgt.view(),
      xlat_out.view(),
      xlong_out.view(),
      hgt_out.view());
  assert(report.ok());
  assert(report.overlap_cell_count == 25);
  assert(report.exposed_cell_count == 75);

  expect_close(xlat_out.view()(5, 5), linear(100.0F, 1.0F, 10.0F, 0.0, 0.0), "reverse overlap");
  assert(std::isfinite(xlat_out.view()(0, 0)));
  assert(std::isfinite(xlong_out.view()(0, 0)));
  expect_close(xlat_out.view()(0, 0), linear(100.0F, 1.0F, 10.0F, -5.0, -5.0), "reverse xlat");
  expect_close(hgt_out.view()(0, 0), linear(500.0F, 11.0F, 13.0F, 1.0, 1.0), "reverse hgt");
}

}  // namespace

int main() {
  run_forward_refresh_case();
  run_reverse_refresh_case();
  return 0;
}
