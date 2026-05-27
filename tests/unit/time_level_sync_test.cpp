#include "tywrf/nest/time_level_sync.hpp"
#include "tywrf/state.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <type_traits>

namespace {

constexpr float kSentinel = -9'999.0F;

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

struct SyncFields {
  explicit SyncFields(const tywrf::Grid& grid)
      : u(grid.u_layout()),
        v(grid.v_layout()),
        t(grid.mass_layout()),
        w(grid.w_layout()),
        ph(grid.w_layout()),
        mu(grid.surface_layout()),
        tke(grid.mass_layout()) {}

  tywrf::FieldStorage3D<float> u;
  tywrf::FieldStorage3D<float> v;
  tywrf::FieldStorage3D<float> t;
  tywrf::FieldStorage3D<float> w;
  tywrf::FieldStorage3D<float> ph;
  tywrf::FieldStorage2D<float> mu;
  tywrf::FieldStorage3D<float> tke;
};

template <typename Storage>
void fill_storage(Storage& storage, const float value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

void fill_all(SyncFields& fields, const float value) {
  fill_storage(fields.u, value);
  fill_storage(fields.v, value);
  fill_storage(fields.t, value);
  fill_storage(fields.w, value);
  fill_storage(fields.ph, value);
  fill_storage(fields.mu, value);
  fill_storage(fields.tke, value);
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

void fill_indexed_all(SyncFields& fields) {
  fill_indexed_3d(fields.u, 1'000.0F);
  fill_indexed_3d(fields.v, 2'000.0F);
  fill_indexed_3d(fields.t, 3'000.0F);
  fill_indexed_3d(fields.w, 4'000.0F);
  fill_indexed_3d(fields.ph, 5'000.0F);
  fill_indexed_2d(fields.mu, 6'000.0F);
  fill_indexed_3d(fields.tke, 7'000.0F);
}

[[nodiscard]] tywrf::nest::PostStartDomainLevel2Views level2_views(
    const SyncFields& fields) noexcept {
  return {
      fields.u.view(),
      fields.v.view(),
      fields.t.view(),
      fields.w.view(),
      fields.ph.view(),
      fields.mu.view()};
}

[[nodiscard]] tywrf::nest::PostStartDomainLevel1Views level1_views(
    SyncFields& fields) noexcept {
  return {
      fields.u.view(),
      fields.v.view(),
      fields.t.view(),
      fields.w.view(),
      fields.ph.view(),
      fields.mu.view()};
}

[[nodiscard]] tywrf::nest::OptionalTimeLevelSyncTke tke_views(
    const SyncFields& level2,
    SyncFields& level1) noexcept {
  return {level2.tke.view(), level1.tke.view(), true};
}

[[nodiscard]] bool in_window(
    const tywrf::nest::RemapWindow& window,
    const std::int32_t i,
    const std::int32_t j) noexcept {
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

void expect_exposed_copy_2d(
    const tywrf::FieldStorage2D<float>& level2,
    const tywrf::FieldStorage2D<float>& level1,
    const tywrf::nest::RemapWindow& window,
    const std::string_view label) {
  const auto src_layout = level2.layout();
  const auto dst_layout = level1.layout();
  const auto src = level2.view();
  const auto dst = level1.view();

  for (std::int32_t j = 0; j < dst_layout.active_ny(); ++j) {
    for (std::int32_t i = 0; i < dst_layout.active_nx(); ++i) {
      const auto dst_i = dst_layout.i_begin() + i;
      const auto dst_j = dst_layout.j_begin() + j;
      if (in_window(window, i, j)) {
        expect_value(dst(dst_i, dst_j), kSentinel, label);
      } else {
        expect_value(
            dst(dst_i, dst_j),
            src(src_layout.i_begin() + i, src_layout.j_begin() + j),
            label);
      }
    }
  }
}

void expect_exposed_copy_3d(
    const tywrf::FieldStorage3D<float>& level2,
    const tywrf::FieldStorage3D<float>& level1,
    const tywrf::nest::RemapWindow& window,
    const std::string_view label) {
  const auto src_layout = level2.layout();
  const auto dst_layout = level1.layout();
  const auto src = level2.view();
  const auto dst = level1.view();

  for (std::int32_t j = 0; j < dst_layout.active_ny(); ++j) {
    for (std::int32_t k = 0; k < dst_layout.active_nz(); ++k) {
      for (std::int32_t i = 0; i < dst_layout.active_nx(); ++i) {
        const auto dst_i = dst_layout.i_begin() + i;
        const auto dst_j = dst_layout.j_begin() + j;
        const auto dst_k = dst_layout.k_begin() + k;
        if (in_window(window, i, j)) {
          expect_value(dst(dst_i, dst_j, dst_k), kSentinel, label);
        } else {
          expect_value(
              dst(dst_i, dst_j, dst_k),
              src(
                  src_layout.i_begin() + i,
                  src_layout.j_begin() + j,
                  src_layout.k_begin() + k),
              label);
        }
      }
    }
  }
}

void expect_all_sentinel_3d(
    const tywrf::FieldStorage3D<float>& field,
    const std::string_view label) {
  const auto layout = field.layout();
  const auto view = field.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t k = 0; k < layout.nz; ++k) {
      for (std::int32_t i = 0; i < layout.nx; ++i) {
        expect_value(view(i, j, k), kSentinel, label);
      }
    }
  }
}

void expect_halo_2d(
    const tywrf::FieldStorage2D<float>& field,
    const std::string_view label) {
  const auto layout = field.layout();
  const auto view = field.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t i = 0; i < layout.nx; ++i) {
      const bool halo =
          i < layout.i_begin() || i >= layout.i_end() ||
          j < layout.j_begin() || j >= layout.j_end();
      if (halo) {
        expect_value(view(i, j), kSentinel, label);
      }
    }
  }
}

void expect_halo_3d(
    const tywrf::FieldStorage3D<float>& field,
    const std::string_view label) {
  const auto layout = field.layout();
  const auto view = field.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t k = 0; k < layout.nz; ++k) {
      for (std::int32_t i = 0; i < layout.nx; ++i) {
        const bool halo =
            i < layout.i_begin() || i >= layout.i_end() ||
            j < layout.j_begin() || j >= layout.j_end() ||
            k < layout.k_begin() || k >= layout.k_end();
        if (halo) {
          expect_value(view(i, j, k), kSentinel, label);
        }
      }
    }
  }
}

void expect_all_halos(const SyncFields& fields) {
  expect_halo_3d(fields.u, "U halo");
  expect_halo_3d(fields.v, "V halo");
  expect_halo_3d(fields.t, "T halo");
  expect_halo_3d(fields.w, "W halo");
  expect_halo_3d(fields.ph, "PH halo");
  expect_halo_2d(fields.mu, "MU halo");
  expect_halo_3d(fields.tke, "TKE halo");
}

[[nodiscard]] std::uint64_t exposed_2d_points(
    const tywrf::nest::RemapWindow& window,
    const tywrf::FieldLayout2D& layout) noexcept {
  return static_cast<std::uint64_t>(layout.active_nx()) *
             static_cast<std::uint64_t>(layout.active_ny()) -
         static_cast<std::uint64_t>(window.extent_i) *
             static_cast<std::uint64_t>(window.extent_j);
}

[[nodiscard]] std::uint64_t exposed_3d_points(
    const tywrf::nest::RemapWindow& window,
    const tywrf::FieldLayout3D& layout) noexcept {
  return (static_cast<std::uint64_t>(layout.active_nx()) *
              static_cast<std::uint64_t>(layout.active_ny()) -
          static_cast<std::uint64_t>(window.extent_i) *
              static_cast<std::uint64_t>(window.extent_j)) *
         static_cast<std::uint64_t>(layout.active_nz());
}

}  // namespace

int main() {
  static_assert(std::is_standard_layout_v<tywrf::nest::PostStartDomainLevel2Views>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::PostStartDomainLevel2Views>);
  static_assert(std::is_standard_layout_v<tywrf::nest::PostStartDomainLevel1Views>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::PostStartDomainLevel1Views>);
  static_assert(std::is_standard_layout_v<tywrf::nest::OptionalTimeLevelSyncTke>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::OptionalTimeLevelSyncTke>);
  static_assert(std::is_standard_layout_v<tywrf::nest::TimeLevelSyncReport>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::TimeLevelSyncReport>);

  const auto grid = make_child_grid();
  const auto descriptor = make_small_descriptor();
  const auto from_pose = tywrf::nest::make_domain_pose(
      descriptor, {1, 1, tywrf::nest::IndexBase::one_based});
  const auto to_pose = tywrf::nest::make_domain_pose(
      descriptor, {2, 2, tywrf::nest::IndexBase::one_based});
  assert(from_pose.result.ok());
  assert(to_pose.result.ok());

  const auto plan = tywrf::nest::build_remap_plan(from_pose, to_pose);
  assert(plan.ok());
  assert(plan.mass.extent_i == 3);
  assert(plan.mass.extent_j == 2);
  assert(plan.u.extent_i == 4);
  assert(plan.u.extent_j == 2);
  assert(plan.v.extent_i == 3);
  assert(plan.v.extent_j == 3);

  SyncFields level2(grid);
  SyncFields level1(grid);
  fill_all(level2, kSentinel);
  fill_all(level1, kSentinel);
  fill_indexed_all(level2);

  const auto report =
      tywrf::nest::copy_post_start_domain_level2_to_level1_exposed_cells(
          plan, level2_views(level2), level1_views(level1));
  assert(report.ok());
  assert(report.copied_field_count == 6);
  assert(report.u_point_count == 14);
  assert(report.v_point_count == 14);
  assert(report.t_point_count == 12);
  assert(report.w_point_count == 18);
  assert(report.ph_point_count == 18);
  assert(report.mu_point_count == 6);
  assert(report.tke_point_count == 0);
  assert(report.copied_point_count == 82);
  assert(!report.copied_optional_tke);
  assert(!report.copied_halo_cells);

  expect_exposed_copy_3d(level2.u, level1.u, plan.u, "U exposed");
  expect_exposed_copy_3d(level2.v, level1.v, plan.v, "V exposed");
  expect_exposed_copy_3d(level2.t, level1.t, plan.mass, "T exposed");
  expect_exposed_copy_3d(level2.w, level1.w, plan.w_full, "W exposed");
  expect_exposed_copy_3d(level2.ph, level1.ph, plan.w_full, "PH exposed");
  expect_exposed_copy_2d(level2.mu, level1.mu, plan.surface, "MU exposed");
  expect_all_sentinel_3d(level1.tke, "absent TKE");
  expect_all_halos(level1);

  SyncFields level1_with_tke(grid);
  fill_all(level1_with_tke, kSentinel);
  const auto tke_report =
      tywrf::nest::copy_post_start_domain_level2_to_level1_exposed_cells(
          plan,
          level2_views(level2),
          level1_views(level1_with_tke),
          tke_views(level2, level1_with_tke));
  assert(tke_report.ok());
  assert(tke_report.copied_field_count == 7);
  assert(tke_report.tke_point_count == 12);
  assert(tke_report.copied_point_count == 94);
  assert(tke_report.copied_optional_tke);
  expect_exposed_copy_3d(level2.tke, level1_with_tke.tke, plan.mass, "TKE exposed");
  expect_all_halos(level1_with_tke);

  SyncFields full_overlap_level1(grid);
  fill_all(full_overlap_level1, kSentinel);
  const auto full_plan = tywrf::nest::build_remap_plan(from_pose, from_pose);
  assert(full_plan.ok());
  const auto full_report =
      tywrf::nest::copy_post_start_domain_level2_to_level1_exposed_cells(
          full_plan,
          level2_views(level2),
          level1_views(full_overlap_level1),
          tke_views(level2, full_overlap_level1));
  assert(full_report.ok());
  assert(full_report.copied_field_count == 0);
  assert(full_report.copied_point_count == 0);
  assert(full_report.u_point_count == 0);
  assert(full_report.v_point_count == 0);
  assert(full_report.t_point_count == 0);
  assert(full_report.w_point_count == 0);
  assert(full_report.ph_point_count == 0);
  assert(full_report.mu_point_count == 0);
  assert(full_report.tke_point_count == 0);
  assert(!full_report.copied_optional_tke);
  expect_exposed_copy_3d(level2.u, full_overlap_level1.u, full_plan.u, "full U");
  expect_exposed_copy_3d(level2.tke, full_overlap_level1.tke, full_plan.mass, "full TKE");
  expect_all_halos(full_overlap_level1);

  const auto expected_u = exposed_3d_points(plan.u, level1.u.layout());
  const auto expected_v = exposed_3d_points(plan.v, level1.v.layout());
  const auto expected_t = exposed_3d_points(plan.mass, level1.t.layout());
  const auto expected_w = exposed_3d_points(plan.w_full, level1.w.layout());
  const auto expected_mu = exposed_2d_points(plan.surface, level1.mu.layout());
  assert(report.u_point_count == expected_u);
  assert(report.v_point_count == expected_v);
  assert(report.t_point_count == expected_t);
  assert(report.w_point_count == expected_w);
  assert(report.ph_point_count == expected_w);
  assert(report.mu_point_count == expected_mu);

  auto invalid_level2 = level2_views(level2);
  invalid_level2.u.data = nullptr;
  const auto null_report =
      tywrf::nest::copy_post_start_domain_level2_to_level1_exposed_cells(
          plan, invalid_level2, level1_views(level1));
  assert(!null_report.ok());
  assert(null_report.result.status == tywrf::nest::NestStatus::invalid_contract);

  auto invalid_level1 = level1_views(level1);
  ++invalid_level1.t.nx;
  const auto shape_report =
      tywrf::nest::copy_post_start_domain_level2_to_level1_exposed_cells(
          plan, level2_views(level2), invalid_level1);
  assert(!shape_report.ok());
  assert(shape_report.result.status == tywrf::nest::NestStatus::invalid_contract);

  auto invalid_tke = tke_views(level2, level1);
  ++invalid_tke.level1.nz;
  const auto tke_shape_report =
      tywrf::nest::copy_post_start_domain_level2_to_level1_exposed_cells(
          plan, level2_views(level2), level1_views(level1), invalid_tke);
  assert(!tke_shape_report.ok());
  assert(tke_shape_report.result.status == tywrf::nest::NestStatus::invalid_contract);

  return 0;
}
