#include "tywrf/nest/base_state_exchange.hpp"

#include <array>
#include <cstdint>

namespace tywrf::nest {
namespace {

struct ExposedRegionSet {
  std::array<RemapWindow, 4> regions{};
  std::uint8_t count = 0;
  std::uint64_t horizontal_cell_count = 0;
};

[[nodiscard]] constexpr NestResult ok_result() noexcept {
  return {NestStatus::ok, "ok"};
}

[[nodiscard]] constexpr NestResult result(
    const NestStatus status,
    const char* message) noexcept {
  return {status, message};
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(
    const FieldView2D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(
    const FieldView2D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(
    const FieldView3D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(
    const FieldView3D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nz(
    const FieldView3D<Real>& field) noexcept {
  return field.nz - field.halo.k_lower - field.halo.k_upper;
}

template <typename Real>
[[nodiscard]] constexpr bool valid_canonical_view(
    const FieldView2D<Real>& field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 &&
         field.stride_i == 1 && field.stride_j == field.nx &&
         field.halo.i_lower >= 0 && field.halo.i_upper >= 0 &&
         field.halo.j_lower >= 0 && field.halo.j_upper >= 0 &&
         active_nx(field) > 0 && active_ny(field) > 0;
}

template <typename Real>
[[nodiscard]] constexpr bool valid_canonical_view(
    const FieldView3D<Real>& field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 &&
         field.nz > 0 && field.stride_i == 1 && field.stride_k == field.nx &&
         field.stride_j == field.nx * field.nz && field.halo.i_lower >= 0 &&
         field.halo.i_upper >= 0 && field.halo.j_lower >= 0 &&
         field.halo.j_upper >= 0 && field.halo.k_lower >= 0 &&
         field.halo.k_upper >= 0 && active_nx(field) > 0 &&
         active_ny(field) > 0 && active_nz(field) > 0;
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_active_shape(
    const FieldView2D<LhsReal>& lhs,
    const FieldView2D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_active_shape(
    const FieldView3D<LhsReal>& lhs,
    const FieldView3D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs) &&
         active_nz(lhs) == active_nz(rhs);
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_horizontal_active_shape(
    const FieldView2D<LhsReal>& lhs,
    const FieldView3D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_horizontal_active_shape(
    const FieldView3D<LhsReal>& lhs,
    const FieldView2D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_horizontal_active_shape(
    const FieldView3D<LhsReal>& lhs,
    const FieldView3D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

[[nodiscard]] constexpr bool window_fits_active_shape(
    const RemapWindow& window,
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0 &&
         window.new_i_begin + window.extent_i <= nx &&
         window.new_j_begin + window.extent_j <= ny;
}

void append_region(
    ExposedRegionSet& set,
    const HorizontalStagger stagger,
    const std::int32_t i_begin,
    const std::int32_t j_begin,
    const std::int32_t extent_i,
    const std::int32_t extent_j) noexcept {
  if (extent_i <= 0 || extent_j <= 0 || set.count >= set.regions.size()) {
    return;
  }

  auto& region = set.regions[set.count];
  region.stagger = stagger;
  region.new_i_begin = i_begin;
  region.new_j_begin = j_begin;
  region.extent_i = extent_i;
  region.extent_j = extent_j;
  region.old_i_begin = 0;
  region.old_j_begin = 0;
  region.old_i_offset_from_new = 0;
  region.old_j_offset_from_new = 0;
  ++set.count;
  set.horizontal_cell_count +=
      static_cast<std::uint64_t>(extent_i) * static_cast<std::uint64_t>(extent_j);
}

[[nodiscard]] ExposedRegionSet exposed_regions_from_overlap(
    const RemapWindow& overlap,
    const std::int32_t active_nx_value,
    const std::int32_t active_ny_value) noexcept {
  ExposedRegionSet set{};
  const auto overlap_i0 = overlap.new_i_begin;
  const auto overlap_j0 = overlap.new_j_begin;
  const auto overlap_i1 = overlap.new_i_begin + overlap.extent_i;
  const auto overlap_j1 = overlap.new_j_begin + overlap.extent_j;

  append_region(set, overlap.stagger, 0, 0, active_nx_value, overlap_j0);
  append_region(
      set,
      overlap.stagger,
      0,
      overlap_j1,
      active_nx_value,
      active_ny_value - overlap_j1);
  append_region(set, overlap.stagger, 0, overlap_j0, overlap_i0, overlap.extent_j);
  append_region(
      set,
      overlap.stagger,
      overlap_i1,
      overlap_j0,
      active_nx_value - overlap_i1,
      overlap.extent_j);
  return set;
}

template <typename Real>
[[nodiscard]] NestResult validate_window(
    const RemapWindow& window,
    const FieldView2D<Real>& field,
    const char* message) noexcept {
  if (!window_fits_active_shape(window, active_nx(field), active_ny(field))) {
    return result(NestStatus::invalid_contract, message);
  }
  return ok_result();
}

template <typename Real>
[[nodiscard]] NestResult validate_window(
    const RemapWindow& window,
    const FieldView3D<Real>& field,
    const char* message) noexcept {
  if (!window_fits_active_shape(window, active_nx(field), active_ny(field))) {
    return result(NestStatus::invalid_contract, message);
  }
  return ok_result();
}

[[nodiscard]] NestResult validate_views(
    const RemapPlan& remap_plan,
    const ExposedBaseStateViews<const float>& source,
    const ExposedBaseStateViews<float>& child) noexcept {
  if (!valid_canonical_view(source.phb) || !valid_canonical_view(child.phb) ||
      !valid_canonical_view(source.mub) || !valid_canonical_view(child.mub) ||
      !valid_canonical_view(source.ht) || !valid_canonical_view(child.ht) ||
      !valid_canonical_view(child.pb) || !valid_canonical_view(child.t_init) ||
      !valid_canonical_view(child.alb)) {
    return result(
        NestStatus::invalid_contract,
        "exposed base-state exchange requires non-null canonical views");
  }

  if (!same_active_shape(source.phb, child.phb) ||
      !same_active_shape(source.mub, child.mub) ||
      !same_active_shape(source.ht, child.ht)) {
    return result(
        NestStatus::invalid_contract,
        "exposed base-state source and child views must share active shapes");
  }

  if (!same_active_shape(child.pb, child.t_init) ||
      !same_active_shape(child.pb, child.alb) ||
      !same_horizontal_active_shape(child.pb, child.mub) ||
      !same_horizontal_active_shape(child.pb, child.ht) ||
      !same_horizontal_active_shape(child.pb, child.phb)) {
    return result(
        NestStatus::invalid_contract,
        "exposed base-state child views must share mass-grid horizontal shape");
  }

  if (const auto validation = validate_window(
          remap_plan.w_full,
          child.phb,
          "exposed base-state PHB remap window is outside child active extents");
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_window(
          remap_plan.surface,
          child.mub,
          "exposed base-state MUB remap window is outside child active extents");
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_window(
          remap_plan.surface,
          child.ht,
          "exposed base-state HT remap window is outside child active extents");
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_window(
          remap_plan.mass,
          child.pb,
          "exposed base-state mass remap window is outside child active extents");
      !validation.ok()) {
    return validation;
  }

  return ok_result();
}

void copy_exposed_3d(
    const ExposedRegionSet& regions,
    const FieldView3D<const float>& source,
    const FieldView3D<float>& child) noexcept {
  for (std::uint8_t region_index = 0; region_index < regions.count; ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t k = 0; k < active_nz(child); ++k) {
        for (std::int32_t i = region.new_i_begin;
             i < region.new_i_begin + region.extent_i;
             ++i) {
          child(
              child.halo.i_lower + i,
              child.halo.j_lower + j,
              child.halo.k_lower + k) =
              source(
                  source.halo.i_lower + i,
                  source.halo.j_lower + j,
                  source.halo.k_lower + k);
        }
      }
    }
  }
}

void copy_exposed_2d(
    const ExposedRegionSet& regions,
    const FieldView2D<const float>& source,
    const FieldView2D<float>& child) noexcept {
  for (std::uint8_t region_index = 0; region_index < regions.count; ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t i = region.new_i_begin;
           i < region.new_i_begin + region.extent_i;
           ++i) {
        child(child.halo.i_lower + i, child.halo.j_lower + j) =
            source(source.halo.i_lower + i, source.halo.j_lower + j);
      }
    }
  }
}

}  // namespace

ExposedBaseStateExchangeReport apply_exposed_base_state_exchange(
    const RemapPlan& remap_plan,
    const ExposedBaseStateViews<const float>& source_interpolated_base_state,
    const ExposedBaseStateViews<float>& child_base_state,
    const ExposedBaseStateExchangeOptions options) noexcept {
  ExposedBaseStateExchangeReport report{};
  report.rebalance_requested = options.rebalance;

  if (!remap_plan.ok()) {
    report.result = remap_plan.result;
    return report;
  }

  if (options.rebalance) {
    report.result = result(
        NestStatus::not_implemented,
        "diagnostic exposed base-state exchange does not implement rebalancing");
    return report;
  }

  if (const auto validation =
          validate_views(
              remap_plan, source_interpolated_base_state, child_base_state);
      !validation.ok()) {
    report.result = validation;
    return report;
  }

  const auto mass_regions = exposed_regions_from_overlap(
      remap_plan.mass, active_nx(child_base_state.pb), active_ny(child_base_state.pb));
  const auto surface_regions = exposed_regions_from_overlap(
      remap_plan.surface,
      active_nx(child_base_state.mub),
      active_ny(child_base_state.mub));
  const auto w_full_regions = exposed_regions_from_overlap(
      remap_plan.w_full,
      active_nx(child_base_state.phb),
      active_ny(child_base_state.phb));

  report.result = ok_result();
  report.exposed_region_count = static_cast<std::uint32_t>(
      mass_regions.count + surface_regions.count + w_full_regions.count);
  report.exposed_mass_cell_count = mass_regions.horizontal_cell_count;
  report.exposed_surface_cell_count = surface_regions.horizontal_cell_count;
  report.exposed_w_full_column_count = w_full_regions.horizontal_cell_count;

  if (w_full_regions.horizontal_cell_count != 0) {
    copy_exposed_3d(
        w_full_regions, source_interpolated_base_state.phb, child_base_state.phb);
    report.wrote_phb = true;
    report.phb_written_point_count =
        w_full_regions.horizontal_cell_count *
        static_cast<std::uint64_t>(active_nz(child_base_state.phb));
  }

  if (surface_regions.horizontal_cell_count != 0) {
    copy_exposed_2d(
        surface_regions, source_interpolated_base_state.mub, child_base_state.mub);
    copy_exposed_2d(
        surface_regions, source_interpolated_base_state.ht, child_base_state.ht);
    report.wrote_mub = true;
    report.wrote_ht = true;
    report.mub_written_cell_count = surface_regions.horizontal_cell_count;
    report.ht_written_cell_count = surface_regions.horizontal_cell_count;
  }

  const auto mass_point_count =
      mass_regions.horizontal_cell_count *
      static_cast<std::uint64_t>(active_nz(child_base_state.pb));
  report.pb_recompute_mark_count = mass_point_count;
  report.t_init_recompute_mark_count = mass_point_count;
  report.alb_recompute_mark_count = mass_point_count;
  report.direct_write_point_count =
      report.phb_written_point_count + report.mub_written_cell_count +
      report.ht_written_cell_count;
  report.recompute_mark_count =
      report.pb_recompute_mark_count + report.t_init_recompute_mark_count +
      report.alb_recompute_mark_count;

  return report;
}

}  // namespace tywrf::nest
