#include "tywrf/nest/state_remap.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <type_traits>

namespace {

constexpr float kSentinel = -9'999.0F;
constexpr float kParentFill = 123'456.0F;

tywrf::Grid make_child_grid() {
  return tywrf::Grid({
      .mass_nx = 4,
      .mass_ny = 3,
      .mass_nz = 2,
      .full_nz = 3,
      .halo = tywrf::Halo3D{1, 1, 1, 1, 1, 1},
  });
}

tywrf::nest::ParentChildDescriptor make_small_descriptor() {
  return {
      tywrf::nest::HorizontalDomainDescriptor{1, 2'000, 8, 8, 9, 9},
      tywrf::nest::HorizontalDomainDescriptor{2, 2'000, 4, 3, 5, 4},
      1,
      1,
  };
}

template <typename Storage>
void fill_storage(Storage& storage, const float value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

void fill_state(tywrf::State<float>& state, const float value) {
  fill_storage(state.u, value);
  fill_storage(state.v, value);
  fill_storage(state.w, value);
  fill_storage(state.ph, value);
  fill_storage(state.phb, value);
  fill_storage(state.t, value);
  fill_storage(state.p, value);
  fill_storage(state.pb, value);
  fill_storage(state.qvapor, value);
  fill_storage(state.qcloud, value);
  fill_storage(state.qrain, value);
  fill_storage(state.qice, value);
  fill_storage(state.qsnow, value);
  fill_storage(state.qgraup, value);
  fill_storage(state.qnice, value);
  fill_storage(state.qnrain, value);
  fill_storage(state.mu, value);
  fill_storage(state.mub, value);
  fill_storage(state.psfc, value);
  fill_storage(state.u10, value);
  fill_storage(state.v10, value);
  fill_storage(state.t2, value);
  fill_storage(state.q2, value);
  fill_storage(state.rainc, value);
  fill_storage(state.rainnc, value);
}

void fill_indexed_2d(tywrf::FieldStorage2D<float>& field, const float base) {
  const auto layout = field.layout();
  auto view = field.view();
  for (std::int32_t j = layout.j_begin(); j < layout.j_end(); ++j) {
    for (std::int32_t i = layout.i_begin(); i < layout.i_end(); ++i) {
      view(i, j) = base + static_cast<float>(100 * i + j);
    }
  }
}

void fill_indexed_3d(tywrf::FieldStorage3D<float>& field, const float base) {
  const auto layout = field.layout();
  auto view = field.view();
  for (std::int32_t j = layout.j_begin(); j < layout.j_end(); ++j) {
    for (std::int32_t k = layout.k_begin(); k < layout.k_end(); ++k) {
      for (std::int32_t i = layout.i_begin(); i < layout.i_end(); ++i) {
        view(i, j, k) = base + static_cast<float>(10'000 * i + 100 * j + k);
      }
    }
  }
}

void fill_indexed_state(tywrf::State<float>& state) {
  fill_indexed_3d(state.u, 1'000.0F);
  fill_indexed_3d(state.v, 2'000.0F);
  fill_indexed_3d(state.w, 3'000.0F);
  fill_indexed_3d(state.ph, 4'000.0F);
  fill_indexed_3d(state.phb, 5'000.0F);
  fill_indexed_3d(state.t, 6'000.0F);
  fill_indexed_3d(state.p, 7'000.0F);
  fill_indexed_3d(state.pb, 8'000.0F);
  fill_indexed_3d(state.qvapor, 9'000.0F);
  fill_indexed_3d(state.qcloud, 10'000.0F);
  fill_indexed_3d(state.qrain, 11'000.0F);
  fill_indexed_3d(state.qice, 12'000.0F);
  fill_indexed_3d(state.qsnow, 13'000.0F);
  fill_indexed_3d(state.qgraup, 14'000.0F);
  fill_indexed_3d(state.qnice, 15'000.0F);
  fill_indexed_3d(state.qnrain, 16'000.0F);
  fill_indexed_2d(state.mu, 17'000.0F);
  fill_indexed_2d(state.mub, 18'000.0F);
  fill_indexed_2d(state.psfc, 19'000.0F);
  fill_indexed_2d(state.u10, 20'000.0F);
  fill_indexed_2d(state.v10, 21'000.0F);
  fill_indexed_2d(state.t2, 22'000.0F);
  fill_indexed_2d(state.q2, 23'000.0F);
  fill_indexed_2d(state.rainc, 24'000.0F);
  fill_indexed_2d(state.rainnc, 25'000.0F);
}

[[nodiscard]] bool in_window(
    const tywrf::nest::RemapWindow& window,
    const std::int32_t i,
    const std::int32_t j) {
  return i >= window.new_i_begin && i < window.new_i_begin + window.extent_i &&
         j >= window.new_j_begin && j < window.new_j_begin + window.extent_j;
}

void expect_value(
    const float actual,
    const float expected,
    const std::string_view label) {
  if (actual != expected) {
    std::cerr << label << " value mismatch: got " << actual << ", expected "
              << expected << '\n';
    assert(false);
  }
}

void expect_overlap_2d(
    const tywrf::FieldStorage2D<float>& old_field,
    const tywrf::FieldStorage2D<float>& new_field,
    const tywrf::nest::RemapWindow& window,
    const std::string_view label) {
  const auto old_layout = old_field.layout();
  const auto new_layout = new_field.layout();
  const auto old_view = old_field.view();
  const auto new_view = new_field.view();

  for (std::int32_t j = 0; j < new_layout.active_ny(); ++j) {
    for (std::int32_t i = 0; i < new_layout.active_nx(); ++i) {
      const auto new_i = new_layout.i_begin() + i;
      const auto new_j = new_layout.j_begin() + j;
      if (in_window(window, i, j)) {
        const auto old_i = old_layout.i_begin() + i + window.old_i_offset_from_new;
        const auto old_j = old_layout.j_begin() + j + window.old_j_offset_from_new;
        expect_value(new_view(new_i, new_j), old_view(old_i, old_j), label);
      } else {
        expect_value(new_view(new_i, new_j), kSentinel, label);
      }
    }
  }
}

void expect_parent_fill_spots_2d(
    const tywrf::FieldStorage2D<float>& old_field,
    const tywrf::FieldStorage2D<float>& parent_field,
    const tywrf::FieldStorage2D<float>& new_field,
    const tywrf::nest::RemapWindow& window,
    const std::string_view label) {
  const auto old_layout = old_field.layout();
  const auto parent_layout = parent_field.layout();
  const auto new_layout = new_field.layout();
  const auto old_view = old_field.view();
  const auto parent_view = parent_field.view();
  const auto new_view = new_field.view();

  const auto overlap_new_i = new_layout.i_begin() + window.new_i_begin;
  const auto overlap_new_j = new_layout.j_begin() + window.new_j_begin;
  const auto overlap_old_i = old_layout.i_begin() + window.old_i_begin;
  const auto overlap_old_j = old_layout.j_begin() + window.old_j_begin;
  expect_value(
      new_view(overlap_new_i, overlap_new_j),
      old_view(overlap_old_i, overlap_old_j),
      label);

  const auto exposed_i = new_layout.active_nx() - 1;
  const auto exposed_j = new_layout.active_ny() - 1;
  assert(!in_window(window, exposed_i, exposed_j));
  expect_value(
      new_view(new_layout.i_begin() + exposed_i, new_layout.j_begin() + exposed_j),
      parent_view(parent_layout.i_begin() + exposed_i, parent_layout.j_begin() + exposed_j),
      label);
}

void expect_overlap_3d(
    const tywrf::FieldStorage3D<float>& old_field,
    const tywrf::FieldStorage3D<float>& new_field,
    const tywrf::nest::RemapWindow& window,
    const std::string_view label) {
  const auto old_layout = old_field.layout();
  const auto new_layout = new_field.layout();
  const auto old_view = old_field.view();
  const auto new_view = new_field.view();

  for (std::int32_t j = 0; j < new_layout.active_ny(); ++j) {
    for (std::int32_t k = 0; k < new_layout.active_nz(); ++k) {
      for (std::int32_t i = 0; i < new_layout.active_nx(); ++i) {
        const auto new_i = new_layout.i_begin() + i;
        const auto new_j = new_layout.j_begin() + j;
        const auto new_k = new_layout.k_begin() + k;
        if (in_window(window, i, j)) {
          const auto old_i = old_layout.i_begin() + i + window.old_i_offset_from_new;
          const auto old_j = old_layout.j_begin() + j + window.old_j_offset_from_new;
          const auto old_k = old_layout.k_begin() + k;
          expect_value(new_view(new_i, new_j, new_k), old_view(old_i, old_j, old_k), label);
        } else {
          expect_value(new_view(new_i, new_j, new_k), kSentinel, label);
        }
      }
    }
  }
}

void expect_parent_fill_spots_3d(
    const tywrf::FieldStorage3D<float>& old_field,
    const tywrf::FieldStorage3D<float>& parent_field,
    const tywrf::FieldStorage3D<float>& new_field,
    const tywrf::nest::RemapWindow& window,
    const std::string_view label) {
  const auto old_layout = old_field.layout();
  const auto parent_layout = parent_field.layout();
  const auto new_layout = new_field.layout();
  const auto old_view = old_field.view();
  const auto parent_view = parent_field.view();
  const auto new_view = new_field.view();

  const auto overlap_new_i = new_layout.i_begin() + window.new_i_begin;
  const auto overlap_new_j = new_layout.j_begin() + window.new_j_begin;
  const auto overlap_new_k = new_layout.k_begin();
  const auto overlap_old_i = old_layout.i_begin() + window.old_i_begin;
  const auto overlap_old_j = old_layout.j_begin() + window.old_j_begin;
  const auto overlap_old_k = old_layout.k_begin();
  expect_value(
      new_view(overlap_new_i, overlap_new_j, overlap_new_k),
      old_view(overlap_old_i, overlap_old_j, overlap_old_k),
      label);

  const auto exposed_i = new_layout.active_nx() - 1;
  const auto exposed_j = new_layout.active_ny() - 1;
  assert(!in_window(window, exposed_i, exposed_j));
  expect_value(
      new_view(
          new_layout.i_begin() + exposed_i,
          new_layout.j_begin() + exposed_j,
          new_layout.k_begin()),
      parent_view(
          parent_layout.i_begin() + exposed_i,
          parent_layout.j_begin() + exposed_j,
          parent_layout.k_begin()),
      label);
}

[[nodiscard]] std::uint64_t points_2d(
    const tywrf::nest::RemapWindow& window) {
  return static_cast<std::uint64_t>(window.extent_i) *
         static_cast<std::uint64_t>(window.extent_j);
}

[[nodiscard]] std::uint64_t points_3d(
    const tywrf::nest::RemapWindow& window,
    const std::int32_t nz) {
  return points_2d(window) * static_cast<std::uint64_t>(nz);
}

[[nodiscard]] std::uint64_t expected_copied_points(
    const tywrf::nest::RemapPlan& plan,
    const tywrf::Grid& grid) {
  const auto mass_nz = grid.config().mass_nz;
  const auto full_nz = grid.config().full_nz;
  return points_3d(plan.u, mass_nz) + points_3d(plan.v, mass_nz) +
         3U * points_3d(plan.w_full, full_nz) +
         11U * points_3d(plan.mass, mass_nz) + 9U * points_2d(plan.surface);
}

[[nodiscard]] std::uint64_t exposed_horizontal_cells(
    const tywrf::nest::RemapWindow& window,
    const std::int32_t nx,
    const std::int32_t ny) {
  return static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny) -
         points_2d(window);
}

void expect_all_supported_fields(
    const tywrf::State<float>& old_state,
    const tywrf::State<float>& new_state,
    const tywrf::nest::RemapPlan& plan) {
  expect_overlap_3d(old_state.u, new_state.u, plan.u, "U");
  expect_overlap_3d(old_state.v, new_state.v, plan.v, "V");
  expect_overlap_3d(old_state.w, new_state.w, plan.w_full, "W");
  expect_overlap_3d(old_state.ph, new_state.ph, plan.w_full, "PH");
  expect_overlap_3d(old_state.phb, new_state.phb, plan.w_full, "PHB");
  expect_overlap_3d(old_state.t, new_state.t, plan.mass, "T");
  expect_overlap_3d(old_state.p, new_state.p, plan.mass, "P");
  expect_overlap_3d(old_state.pb, new_state.pb, plan.mass, "PB");
  expect_overlap_3d(old_state.qvapor, new_state.qvapor, plan.mass, "QVAPOR");
  expect_overlap_3d(old_state.qcloud, new_state.qcloud, plan.mass, "QCLOUD");
  expect_overlap_3d(old_state.qrain, new_state.qrain, plan.mass, "QRAIN");
  expect_overlap_3d(old_state.qice, new_state.qice, plan.mass, "QICE");
  expect_overlap_3d(old_state.qsnow, new_state.qsnow, plan.mass, "QSNOW");
  expect_overlap_3d(old_state.qgraup, new_state.qgraup, plan.mass, "QGRAUP");
  expect_overlap_3d(old_state.qnice, new_state.qnice, plan.mass, "QNICE");
  expect_overlap_3d(old_state.qnrain, new_state.qnrain, plan.mass, "QNRAIN");
  expect_overlap_2d(old_state.mu, new_state.mu, plan.surface, "MU");
  expect_overlap_2d(old_state.mub, new_state.mub, plan.surface, "MUB");
  expect_overlap_2d(old_state.psfc, new_state.psfc, plan.surface, "PSFC");
  expect_overlap_2d(old_state.u10, new_state.u10, plan.surface, "U10");
  expect_overlap_2d(old_state.v10, new_state.v10, plan.surface, "V10");
  expect_overlap_2d(old_state.t2, new_state.t2, plan.surface, "T2");
  expect_overlap_2d(old_state.q2, new_state.q2, plan.surface, "Q2");
  expect_overlap_2d(old_state.rainc, new_state.rainc, plan.surface, "RAINC");
  expect_overlap_2d(old_state.rainnc, new_state.rainnc, plan.surface, "RAINNC");
}

void expect_representative_parent_fill_fields(
    const tywrf::State<float>& old_state,
    const tywrf::State<float>& parent_state,
    const tywrf::State<float>& new_state,
    const tywrf::nest::RemapPlan& plan) {
  expect_parent_fill_spots_3d(old_state.u, parent_state.u, new_state.u, plan.u, "U");
  expect_parent_fill_spots_3d(old_state.v, parent_state.v, new_state.v, plan.v, "V");
  expect_parent_fill_spots_3d(old_state.w, parent_state.w, new_state.w, plan.w_full, "W");
  expect_parent_fill_spots_3d(old_state.t, parent_state.t, new_state.t, plan.mass, "T");
  expect_parent_fill_spots_2d(old_state.mu, parent_state.mu, new_state.mu, plan.surface, "MU");
}

void expect_representative_halos(const tywrf::State<float>& state) {
  const auto u_layout = state.u.layout();
  const auto u = state.u.view();
  const auto t_layout = state.t.layout();
  const auto t = state.t.view();
  const auto mu_layout = state.mu.layout();
  const auto mu = state.mu.view();
  expect_value(u(0, u_layout.j_begin(), u_layout.k_begin()), kSentinel, "U halo");
  expect_value(t(0, t_layout.j_begin(), t_layout.k_begin()), kSentinel, "T halo");
  expect_value(mu(0, mu_layout.j_begin()), kSentinel, "MU halo");
}

}  // namespace

int main() {
  static_assert(std::is_standard_layout_v<tywrf::nest::ChildStateRemapReport>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::ChildStateRemapReport>);

  const auto grid = make_child_grid();
  tywrf::State<float> old_state(grid);
  tywrf::State<float> new_state(grid);
  fill_state(old_state, kSentinel);
  fill_state(new_state, kSentinel);
  fill_indexed_state(old_state);

  tywrf::State<float> parent_filled_state(grid);
  fill_state(parent_filled_state, kParentFill);

  const auto descriptor = make_small_descriptor();
  const auto from_pose = tywrf::nest::make_domain_pose(
      descriptor, {1, 1, tywrf::nest::IndexBase::one_based});
  const auto to_pose = tywrf::nest::make_domain_pose(
      descriptor, {2, 2, tywrf::nest::IndexBase::one_based});
  assert(from_pose.result.ok());
  assert(to_pose.result.ok());

  const auto plan = tywrf::nest::build_remap_plan(from_pose, to_pose);
  assert(plan.ok());
  assert(plan.delta.child_di == 1);
  assert(plan.delta.child_dj == 1);
  assert(plan.mass.extent_i == 3);
  assert(plan.mass.extent_j == 2);
  assert(plan.u.extent_i == 4);
  assert(plan.u.extent_j == 2);
  assert(plan.v.extent_i == 3);
  assert(plan.v.extent_j == 3);
  const auto mass_exposed = exposed_horizontal_cells(
      plan.mass, grid.mass_shape().nx, grid.mass_shape().ny);
  const auto u_exposed = exposed_horizontal_cells(
      plan.u, grid.u_shape().nx, grid.u_shape().ny);
  const auto v_exposed = exposed_horizontal_cells(
      plan.v, grid.v_shape().nx, grid.v_shape().ny);
  const auto w_exposed = exposed_horizontal_cells(
      plan.w_full, grid.w_shape().nx, grid.w_shape().ny);
  const auto surface_exposed = exposed_horizontal_cells(
      plan.surface, grid.surface_shape().nx, grid.surface_shape().ny);
  assert(mass_exposed == 6);
  assert(u_exposed == 7);
  assert(v_exposed == 7);
  assert(w_exposed == 6);
  assert(surface_exposed == 6);
  const auto expected_parent_fill_points =
      u_exposed * static_cast<std::uint64_t>(grid.u_shape().nz) +
      v_exposed * static_cast<std::uint64_t>(grid.v_shape().nz) +
      3U * w_exposed * static_cast<std::uint64_t>(grid.w_shape().nz) +
      11U * mass_exposed * static_cast<std::uint64_t>(grid.mass_shape().nz) +
      9U * surface_exposed;

  const auto report =
      tywrf::nest::remap_child_state_overlap_only(plan, old_state, new_state);
  assert(report.ok());
  assert(report.copied_field_count == 25);
  assert(report.copied_point_count == expected_copied_points(plan, grid));
  assert(report.parent_fill_field_count == 0);
  assert(report.parent_fill_point_count == 0);
  assert(report.needs_parent_fill);
  assert(!report.filled_exposed_cells);
  assert(!report.copied_halo_cells);
  expect_all_supported_fields(old_state, new_state, plan);
  expect_representative_halos(new_state);

  tywrf::State<float> filled_new_state(grid);
  fill_state(filled_new_state, kSentinel);
  const auto fill_report = tywrf::nest::remap_child_state_overlap_with_parent_fill(
      plan, old_state, parent_filled_state, filled_new_state);
  assert(fill_report.ok());
  assert(fill_report.copied_field_count == 25);
  assert(fill_report.copied_point_count == expected_copied_points(plan, grid));
  assert(fill_report.parent_fill_field_count == 25);
  assert(fill_report.parent_fill_point_count == expected_parent_fill_points);
  assert(!fill_report.needs_parent_fill);
  assert(fill_report.filled_exposed_cells);
  assert(!fill_report.copied_halo_cells);
  expect_representative_parent_fill_fields(
      old_state, parent_filled_state, filled_new_state, plan);
  expect_representative_halos(filled_new_state);

  tywrf::State<float> full_new_state(grid);
  fill_state(full_new_state, kSentinel);
  const auto full_plan = tywrf::nest::build_remap_plan(from_pose, from_pose);
  assert(full_plan.ok());
  const auto full_report =
      tywrf::nest::remap_child_state_overlap_only(full_plan, old_state, full_new_state);
  assert(full_report.ok());
  assert(full_report.copied_field_count == 25);
  assert(full_report.copied_point_count == expected_copied_points(full_plan, grid));
  assert(full_report.parent_fill_field_count == 0);
  assert(full_report.parent_fill_point_count == 0);
  assert(!full_report.needs_parent_fill);
  assert(!full_report.filled_exposed_cells);
  assert(!full_report.copied_halo_cells);
  expect_all_supported_fields(old_state, full_new_state, full_plan);
  expect_representative_halos(full_new_state);

  tywrf::State<float> full_filled_new_state(grid);
  fill_state(full_filled_new_state, kSentinel);
  const auto full_fill_report =
      tywrf::nest::remap_child_state_overlap_with_parent_fill(
          full_plan, old_state, parent_filled_state, full_filled_new_state);
  assert(full_fill_report.ok());
  assert(full_fill_report.copied_field_count == 25);
  assert(full_fill_report.copied_point_count ==
         expected_copied_points(full_plan, grid));
  assert(full_fill_report.parent_fill_field_count == 0);
  assert(full_fill_report.parent_fill_point_count == 0);
  assert(!full_fill_report.needs_parent_fill);
  assert(!full_fill_report.filled_exposed_cells);
  assert(!full_fill_report.copied_halo_cells);
  expect_all_supported_fields(old_state, full_filled_new_state, full_plan);
  expect_representative_halos(full_filled_new_state);

  return 0;
}
