#include "tywrf/nest/parent_child_interpolation.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <type_traits>

namespace {

constexpr float kSentinel = -9'999.0F;
constexpr float kTolerance = 1.0e-3F;

tywrf::Grid make_parent_grid() {
  return tywrf::Grid({
      .mass_nx = 8,
      .mass_ny = 8,
      .mass_nz = 3,
      .full_nz = 4,
      .halo = tywrf::Halo3D{1, 1, 1, 1, 1, 1},
  });
}

tywrf::Grid make_child_grid() {
  return tywrf::Grid({
      .mass_nx = 10,
      .mass_ny = 10,
      .mass_nz = 3,
      .full_nz = 4,
      .halo = tywrf::Halo3D{1, 1, 1, 1, 1, 1},
  });
}

tywrf::nest::ParentChildDescriptor make_ratio5_descriptor() {
  return {
      tywrf::nest::HorizontalDomainDescriptor{1, 10'000, 8, 8, 9, 9},
      tywrf::nest::HorizontalDomainDescriptor{2, 2'000, 10, 10, 11, 11},
      5,
      5,
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

float linear_value(
    const float base,
    const float ax,
    const float ay,
    const float ak,
    const double x,
    const double y,
    const std::int32_t k) {
  return base + ax * static_cast<float>(x) + ay * static_cast<float>(y) +
         ak * static_cast<float>(k);
}

void fill_gradient_3d(
    tywrf::FieldStorage3D<float>& field,
    const float base,
    const float ax,
    const float ay,
    const float ak) {
  fill_storage(field, kSentinel);
  const auto layout = field.layout();
  auto view = field.view();
  for (std::int32_t j = layout.j_begin(); j < layout.j_end(); ++j) {
    const auto active_j = j - layout.j_begin();
    for (std::int32_t k = layout.k_begin(); k < layout.k_end(); ++k) {
      const auto active_k = k - layout.k_begin();
      for (std::int32_t i = layout.i_begin(); i < layout.i_end(); ++i) {
        const auto active_i = i - layout.i_begin();
        view(i, j, k) = linear_value(
            base, ax, ay, ak, active_i, active_j, active_k);
      }
    }
  }
}

void fill_gradient_2d(
    tywrf::FieldStorage2D<float>& field,
    const float base,
    const float ax,
    const float ay) {
  fill_storage(field, kSentinel);
  const auto layout = field.layout();
  auto view = field.view();
  for (std::int32_t j = layout.j_begin(); j < layout.j_end(); ++j) {
    const auto active_j = j - layout.j_begin();
    for (std::int32_t i = layout.i_begin(); i < layout.i_end(); ++i) {
      const auto active_i = i - layout.i_begin();
      view(i, j) = linear_value(base, ax, ay, 0.0F, active_i, active_j, 0);
    }
  }
}

void fill_supported_parent_gradients(tywrf::State<float>& parent) {
  fill_state(parent, kSentinel);
  fill_gradient_3d(parent.u, 1'000.0F, 7.0F, 13.0F, 101.0F);
  fill_gradient_3d(parent.v, 2'000.0F, 11.0F, 17.0F, 103.0F);
  fill_gradient_2d(parent.mu, 3'000.0F, 19.0F, 23.0F);
  fill_gradient_3d(parent.qvapor, 4'000.0F, 29.0F, 31.0F, 107.0F);
  fill_gradient_3d(parent.t, 5'000.0F, 5.0F, 3.0F, 10.0F);
  fill_gradient_3d(parent.ph, 6'000.0F, 2.0F, 4.0F, 50.0F);
}

tywrf::nest::StateExchangePlan make_exchange_plan(
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const tywrf::State<float>& child) {
  const auto from_pose = tywrf::nest::make_nest_pose(
      descriptor, {2, 2, tywrf::nest::IndexBase::one_based});
  const auto to_pose = tywrf::nest::make_nest_pose(
      descriptor, {3, 2, tywrf::nest::IndexBase::one_based});
  assert(from_pose.child.result.ok());
  assert(to_pose.child.result.ok());

  const auto remap = tywrf::nest::build_remap_plan(from_pose, to_pose);
  assert(remap.ok());
  assert(remap.delta.child_di == 5);
  assert(remap.delta.child_dj == 0);
  return tywrf::nest::build_exposed_child_state_exchange_plan(
      remap, child.view());
}

[[nodiscard]] bool in_exposed_region(
    const tywrf::nest::FieldStateExchangePlan& field,
    const std::int32_t i,
    const std::int32_t j) {
  for (std::uint8_t region_index = 0; region_index < field.exposed_region_count;
       ++region_index) {
    const auto& region = field.exposed_regions[region_index];
    const bool inside_i = i >= region.child_i_begin &&
                          i < region.child_i_begin + region.extent_i;
    const bool inside_j = j >= region.child_j_begin &&
                          j < region.child_j_begin + region.extent_j;
    if (inside_i && inside_j) {
      return true;
    }
  }
  return false;
}

void expect_close(
    const float actual,
    const float expected,
    const std::string_view label) {
  if (std::fabs(actual - expected) > kTolerance) {
    std::cerr << label << " mismatch: got " << actual << ", expected "
              << expected << '\n';
    assert(false);
  }
}

void expect_value(
    const float actual,
    const float expected,
    const std::string_view label) {
  if (actual != expected) {
    std::cerr << label << " mismatch: got " << actual << ", expected "
              << expected << '\n';
    assert(false);
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
        const bool is_halo =
            i < layout.i_begin() || i >= layout.i_end() ||
            j < layout.j_begin() || j >= layout.j_end() ||
            k < layout.k_begin() || k >= layout.k_end();
        if (is_halo) {
          expect_value(view(i, j, k), kSentinel, label);
        }
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
      const bool is_halo =
          i < layout.i_begin() || i >= layout.i_end() ||
          j < layout.j_begin() || j >= layout.j_end();
      if (is_halo) {
        expect_value(view(i, j), kSentinel, label);
      }
    }
  }
}

void expect_interpolated_3d(
    const tywrf::FieldStorage3D<float>& child_field,
    const tywrf::nest::FieldStateExchangePlan& field_plan,
    const float base,
    const float ax,
    const float ay,
    const float ak,
    const std::string_view label) {
  const auto layout = child_field.layout();
  const auto view = child_field.view();
  constexpr auto ratio = 5;
  constexpr auto parent_i_start_zero = 2;
  constexpr auto parent_j_start_zero = 1;

  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    const auto y = static_cast<double>(parent_j_start_zero) +
                   static_cast<double>(j) / static_cast<double>(ratio);
    for (std::int32_t k = 0; k < layout.active_nz(); ++k) {
      for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
        const auto storage_i = layout.i_begin() + i;
        const auto storage_j = layout.j_begin() + j;
        const auto storage_k = layout.k_begin() + k;
        if (in_exposed_region(field_plan, i, j)) {
          const auto x = static_cast<double>(parent_i_start_zero) +
                         static_cast<double>(i) / static_cast<double>(ratio);
          expect_close(
              view(storage_i, storage_j, storage_k),
              linear_value(base, ax, ay, ak, x, y, k),
              label);
        } else {
          expect_value(view(storage_i, storage_j, storage_k), kSentinel, label);
        }
      }
    }
  }
  expect_halo_3d(child_field, label);
}

void expect_interpolated_2d(
    const tywrf::FieldStorage2D<float>& child_field,
    const tywrf::nest::FieldStateExchangePlan& field_plan,
    const float base,
    const float ax,
    const float ay,
    const std::string_view label) {
  const auto layout = child_field.layout();
  const auto view = child_field.view();
  constexpr auto ratio = 5;
  constexpr auto parent_i_start_zero = 2;
  constexpr auto parent_j_start_zero = 1;

  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    const auto y = static_cast<double>(parent_j_start_zero) +
                   static_cast<double>(j) / static_cast<double>(ratio);
    for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
      const auto storage_i = layout.i_begin() + i;
      const auto storage_j = layout.j_begin() + j;
      if (in_exposed_region(field_plan, i, j)) {
        const auto x = static_cast<double>(parent_i_start_zero) +
                       static_cast<double>(i) / static_cast<double>(ratio);
        expect_close(
            view(storage_i, storage_j),
            linear_value(base, ax, ay, 0.0F, x, y, 0),
            label);
      } else {
        expect_value(view(storage_i, storage_j), kSentinel, label);
      }
    }
  }
  expect_halo_2d(child_field, label);
}

void test_bilinear_exposed_fill_for_selected_staggers() {
  const auto descriptor = make_ratio5_descriptor();
  tywrf::State<float> parent(make_parent_grid());
  tywrf::State<float> child(make_child_grid());
  fill_supported_parent_gradients(parent);
  fill_state(child, kSentinel);

  const auto exchange = make_exchange_plan(descriptor, child);
  assert(exchange.ok());
  assert(exchange.field_count == 6);
  assert(exchange.report.exposed_region_count == 6);
  assert(exchange.fields[0].active_nx == 11);
  assert(exchange.fields[0].active_ny == 10);
  assert(exchange.fields[1].active_nx == 10);
  assert(exchange.fields[1].active_ny == 11);
  assert(exchange.fields[2].active_nx == 10);
  assert(exchange.fields[2].active_ny == 10);
  assert(exchange.fields[3].active_nx == 10);
  assert(exchange.fields[3].active_ny == 10);
  assert(exchange.fields[4].active_nx == 10);
  assert(exchange.fields[4].active_ny == 10);
  assert(exchange.fields[4].active_k_count == 3);
  assert(exchange.fields[5].active_nx == 10);
  assert(exchange.fields[5].active_ny == 10);
  assert(exchange.fields[5].active_k_count == 4);

  const auto report = tywrf::nest::interpolate_parent_to_exposed_child(
      descriptor,
      {3, 2, tywrf::nest::IndexBase::one_based},
      exchange,
      static_cast<const tywrf::State<float>&>(parent).view(),
      child.view());
  assert(report.ok());
  assert(report.method == tywrf::nest::ParentChildInterpolationMethod::bilinear);
  assert(report.requested_field_count == 6);
  assert(report.interpolated_field_count == 6);
  assert(report.interpolated_region_count == 6);
  assert(report.interpolated_horizontal_cell_count == 305);
  assert(report.interpolated_point_count == 865);
  assert(!report.wrote_overlap);
  assert(!report.wrote_halo);

  expect_interpolated_3d(child.u, exchange.fields[0], 1'000.0F, 7.0F, 13.0F, 101.0F, "U");
  expect_interpolated_3d(child.v, exchange.fields[1], 2'000.0F, 11.0F, 17.0F, 103.0F, "V");
  expect_interpolated_2d(child.mu, exchange.fields[2], 3'000.0F, 19.0F, 23.0F, "MU");
  expect_interpolated_3d(
      child.qvapor, exchange.fields[3], 4'000.0F, 29.0F, 31.0F, 107.0F, "QVAPOR");
  expect_interpolated_3d(child.t, exchange.fields[4], 5'000.0F, 5.0F, 3.0F, 10.0F, "T");
  expect_interpolated_3d(child.ph, exchange.fields[5], 6'000.0F, 2.0F, 4.0F, 50.0F, "PH");
}

void test_plan_selects_only_requested_field() {
  const auto descriptor = make_ratio5_descriptor();
  tywrf::State<float> parent(make_parent_grid());
  tywrf::State<float> child(make_child_grid());
  fill_supported_parent_gradients(parent);
  fill_state(child, kSentinel);

  const auto full_exchange = make_exchange_plan(descriptor, child);
  tywrf::nest::StateExchangePlan qvapor_only{};
  qvapor_only.result = {tywrf::nest::NestStatus::ok, "ok"};
  qvapor_only.operation =
      tywrf::nest::ExchangeOperation::parent_to_child_interpolation;
  qvapor_only.field_count = 1;
  qvapor_only.fields[0] = full_exchange.fields[3];

  const auto report = tywrf::nest::interpolate_parent_to_exposed_child(
      descriptor,
      {3, 2, tywrf::nest::IndexBase::one_based},
      qvapor_only,
      static_cast<const tywrf::State<float>&>(parent).view(),
      child.view());
  assert(report.ok());
  assert(report.requested_field_count == 1);
  assert(report.interpolated_field_count == 1);
  assert(report.interpolated_point_count == 150);

  expect_halo_3d(child.u, "U");
  expect_halo_3d(child.v, "V");
  expect_halo_2d(child.mu, "MU");
  expect_halo_3d(child.t, "T");
  expect_halo_3d(child.ph, "PH");
  expect_interpolated_3d(
      child.qvapor,
      full_exchange.fields[3],
      4'000.0F,
      29.0F,
      31.0F,
      107.0F,
      "QVAPOR selected");

  const auto u_layout = child.u.layout();
  const auto u = child.u.view();
  expect_value(
      u(u_layout.i_begin() + u_layout.active_nx() - 1,
        u_layout.j_begin(),
        u_layout.k_begin()),
      kSentinel,
      "unselected U active");
  const auto t_layout = child.t.layout();
  const auto t = child.t.view();
  expect_value(
      t(t_layout.i_begin() + t_layout.active_nx() - 1,
        t_layout.j_begin(),
        t_layout.k_begin()),
      kSentinel,
      "unselected T active");
  const auto ph_layout = child.ph.layout();
  const auto ph = child.ph.view();
  expect_value(
      ph(ph_layout.i_begin() + ph_layout.active_nx() - 1,
         ph_layout.j_begin(),
         ph_layout.k_begin()),
      kSentinel,
      "unselected PH active");
}

void test_rejects_non_krosa_ratio_or_resolution() {
  auto descriptor = make_ratio5_descriptor();
  descriptor.parent_grid_ratio = 4;

  tywrf::State<float> parent(make_parent_grid());
  tywrf::State<float> child(make_child_grid());
  fill_supported_parent_gradients(parent);
  fill_state(child, kSentinel);

  const auto exchange = make_exchange_plan(make_ratio5_descriptor(), child);
  const auto ratio_report = tywrf::nest::interpolate_parent_to_exposed_child(
      descriptor,
      {3, 2, tywrf::nest::IndexBase::one_based},
      exchange,
      static_cast<const tywrf::State<float>&>(parent).view(),
      child.view());
  assert(!ratio_report.ok());
  assert(ratio_report.result.status == tywrf::nest::NestStatus::unsupported_resolution);

  descriptor = make_ratio5_descriptor();
  descriptor.child.grid_spacing_m = 4'000;
  const auto resolution_report = tywrf::nest::interpolate_parent_to_exposed_child(
      descriptor,
      {3, 2, tywrf::nest::IndexBase::one_based},
      exchange,
      static_cast<const tywrf::State<float>&>(parent).view(),
      child.view());
  assert(!resolution_report.ok());
  assert(resolution_report.result.status == tywrf::nest::NestStatus::unsupported_resolution);
}

void test_rejects_unowned_overlap_or_halo_plan() {
  const auto descriptor = make_ratio5_descriptor();
  tywrf::State<float> parent(make_parent_grid());
  tywrf::State<float> child(make_child_grid());
  fill_supported_parent_gradients(parent);
  fill_state(child, kSentinel);

  auto exchange = make_exchange_plan(descriptor, child);
  exchange.fields[3].exposed_regions[0].owns_overlap = true;

  const auto report = tywrf::nest::interpolate_parent_to_exposed_child(
      descriptor,
      {3, 2, tywrf::nest::IndexBase::one_based},
      exchange,
      static_cast<const tywrf::State<float>&>(parent).view(),
      child.view());
  assert(!report.ok());
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(!report.wrote_overlap);
  assert(!report.wrote_halo);
}

}  // namespace

int main() {
  static_assert(
      std::is_standard_layout_v<tywrf::nest::ParentChildInterpolationConfig>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::nest::ParentChildInterpolationConfig>);
  static_assert(
      std::is_standard_layout_v<tywrf::nest::ParentChildInterpolationReport>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::nest::ParentChildInterpolationReport>);

  test_bilinear_exposed_fill_for_selected_staggers();
  test_plan_selects_only_requested_field();
  test_rejects_non_krosa_ratio_or_resolution();
  test_rejects_unowned_overlap_or_halo_plan();
  return 0;
}
