#include "tywrf/nest/base_state_source_staging.hpp"

#include <algorithm>
#include <array>
#include <cmath>
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

[[nodiscard]] constexpr std::uint64_t point_count_2d_u64(
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny);
}

[[nodiscard]] constexpr bool valid_grid_config(const GridConfig& config) noexcept {
  return config.mass_nx > 0 && config.mass_ny > 0 && config.mass_nz > 0 &&
         config.full_nz == config.mass_nz + 1 && config.halo.i_lower >= 0 &&
         config.halo.i_upper >= 0 && config.halo.j_lower >= 0 &&
         config.halo.j_upper >= 0 && config.halo.k_lower >= 0 &&
         config.halo.k_upper >= 0;
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

template <typename Real>
[[nodiscard]] constexpr bool matches_active_shape(
    const FieldView2D<Real>& field,
    const ActiveShape2D expected) noexcept {
  return active_nx(field) == expected.nx && active_ny(field) == expected.ny;
}

template <typename Real>
[[nodiscard]] constexpr bool matches_active_shape(
    const FieldView3D<Real>& field,
    const ActiveShape3D expected) noexcept {
  return active_nx(field) == expected.nx && active_ny(field) == expected.ny &&
         active_nz(field) == expected.nz;
}

[[nodiscard]] constexpr bool window_fits_active_shape(
    const RemapWindow& window,
    const HorizontalStagger expected_stagger,
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return window.stagger == expected_stagger && window.new_i_begin >= 0 &&
         window.new_j_begin >= 0 && window.extent_i > 0 && window.extent_j > 0 &&
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
  ++set.count;
  set.horizontal_cell_count += point_count_2d_u64(extent_i, extent_j);
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

template <typename Storage>
void fill_storage(Storage& storage, const float value) noexcept {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

void copy_exposed_3d(
    const ExposedRegionSet& regions,
    const FieldView3D<const float>& source,
    const FieldView3D<float>& destination) noexcept {
  for (std::uint8_t region_index = 0; region_index < regions.count; ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t k = 0; k < active_nz(destination); ++k) {
        for (std::int32_t i = region.new_i_begin;
             i < region.new_i_begin + region.extent_i;
             ++i) {
          destination(
              destination.halo.i_lower + i,
              destination.halo.j_lower + j,
              destination.halo.k_lower + k) =
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
    const FieldView2D<float>& destination) noexcept {
  for (std::uint8_t region_index = 0; region_index < regions.count; ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t i = region.new_i_begin;
           i < region.new_i_begin + region.extent_i;
           ++i) {
        destination(destination.halo.i_lower + i, destination.halo.j_lower + j) =
            source(source.halo.i_lower + i, source.halo.j_lower + j);
      }
    }
  }
}

[[nodiscard]] std::uint64_t count_nonfinite_exposed_3d(
    const ExposedRegionSet& regions,
    const FieldView3D<const float>& source) noexcept {
  std::uint64_t count = 0;
  for (std::uint8_t region_index = 0; region_index < regions.count; ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t k = 0; k < active_nz(source); ++k) {
        for (std::int32_t i = region.new_i_begin;
             i < region.new_i_begin + region.extent_i;
             ++i) {
          if (!std::isfinite(
                  source(
                      source.halo.i_lower + i,
                      source.halo.j_lower + j,
                      source.halo.k_lower + k))) {
            ++count;
          }
        }
      }
    }
  }
  return count;
}

[[nodiscard]] std::uint64_t count_nonfinite_exposed_2d(
    const ExposedRegionSet& regions,
    const FieldView2D<const float>& source) noexcept {
  std::uint64_t count = 0;
  for (std::uint8_t region_index = 0; region_index < regions.count; ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t i = region.new_i_begin;
           i < region.new_i_begin + region.extent_i;
           ++i) {
        if (!std::isfinite(
                source(source.halo.i_lower + i, source.halo.j_lower + j))) {
          ++count;
        }
      }
    }
  }
  return count;
}

[[nodiscard]] BaseStateSourceStagingReport base_report(
    const Grid& child_grid) noexcept {
  const auto& config = child_grid.config();
  BaseStateSourceStagingReport report{};
  report.active_nx = config.mass_nx;
  report.active_ny = config.mass_ny;
  report.mass_nz = config.mass_nz;
  report.full_nz = config.full_nz;
  report.expected_mass_level_count = config.mass_nz;
  report.expected_full_level_count = config.mass_nz + 1;

  if (config.mass_nx > 0 && config.mass_ny > 0 && config.mass_nz > 0 &&
      config.full_nz > 0) {
    report.active_mass_cell_count =
        point_count_2d_u64(config.mass_nx, config.mass_ny);
    report.active_mass_point_count =
        report.active_mass_cell_count *
        static_cast<std::uint64_t>(config.mass_nz);
    report.active_surface_cell_count = report.active_mass_cell_count;
    report.active_w_full_point_count =
        report.active_mass_cell_count *
        static_cast<std::uint64_t>(config.full_nz);
  }
  return report;
}

void populate_region_counts(
    BaseStateSourceStagingReport& report,
    const Grid& child_grid,
    const ExposedRegionSet& mass_regions,
    const ExposedRegionSet& surface_regions,
    const ExposedRegionSet& w_full_regions) noexcept {
  report.exposed_region_count = static_cast<std::uint32_t>(
      mass_regions.count + surface_regions.count + w_full_regions.count);
  report.exposed_mass_cell_count = mass_regions.horizontal_cell_count;
  report.exposed_mass_point_count =
      mass_regions.horizontal_cell_count *
      static_cast<std::uint64_t>(child_grid.config().mass_nz);
  report.exposed_surface_cell_count = surface_regions.horizontal_cell_count;
  report.exposed_w_full_column_count = w_full_regions.horizontal_cell_count;
  report.exposed_w_full_point_count =
      w_full_regions.horizontal_cell_count *
      static_cast<std::uint64_t>(child_grid.config().full_nz);

  report.masked_mass_cell_count =
      report.active_mass_cell_count - report.exposed_mass_cell_count;
  report.masked_mass_point_count =
      report.active_mass_point_count - report.exposed_mass_point_count;
  report.masked_surface_cell_count =
      report.active_surface_cell_count - report.exposed_surface_cell_count;
  report.masked_w_full_column_count =
      report.active_mass_cell_count - report.exposed_w_full_column_count;
  report.masked_w_full_point_count =
      report.active_w_full_point_count - report.exposed_w_full_point_count;

  report.staged_phb_point_count = report.exposed_w_full_point_count;
  report.staged_mub_cell_count = report.exposed_surface_cell_count;
  report.staged_ht_cell_count = report.exposed_surface_cell_count;
  report.staged_pb_point_count = report.exposed_mass_point_count;
  report.staged_t_init_point_count = report.exposed_mass_point_count;
  report.staged_alb_point_count = report.exposed_mass_point_count;
  report.staged_value_count =
      report.staged_phb_point_count + report.staged_mub_cell_count +
      report.staged_ht_cell_count + report.staged_pb_point_count +
      report.staged_t_init_point_count + report.staged_alb_point_count;
  report.active_masked_value_count =
      report.masked_w_full_point_count +
      (2U * report.masked_surface_cell_count) +
      (3U * report.masked_mass_point_count);

  const auto mass_allocation = child_grid.mass_layout().allocation_size();
  const auto surface_allocation = child_grid.surface_layout().allocation_size();
  const auto w_allocation = child_grid.w_layout().allocation_size();
  const auto active_staged_layout_value_count =
      report.active_w_full_point_count +
      (2U * report.active_surface_cell_count) +
      (3U * report.active_mass_point_count);
  const auto allocated_staging_value_count =
      static_cast<std::uint64_t>(w_allocation) +
      (2U * static_cast<std::uint64_t>(surface_allocation)) +
      (3U * static_cast<std::uint64_t>(mass_allocation));
  report.halo_masked_value_count =
      allocated_staging_value_count - active_staged_layout_value_count;
}

void populate_invalid_exposed_counts(
    BaseStateSourceStagingReport& report,
    const ExposedBaseStateViews<const float>& source,
    const ExposedRegionSet& mass_regions,
    const ExposedRegionSet& surface_regions,
    const ExposedRegionSet& w_full_regions) noexcept {
  report.invalid_exposed_phb_point_count =
      count_nonfinite_exposed_3d(w_full_regions, source.phb);
  report.invalid_exposed_mub_cell_count =
      count_nonfinite_exposed_2d(surface_regions, source.mub);
  report.invalid_exposed_ht_cell_count =
      count_nonfinite_exposed_2d(surface_regions, source.ht);
  report.invalid_exposed_pb_point_count =
      count_nonfinite_exposed_3d(mass_regions, source.pb);
  report.invalid_exposed_t_init_point_count =
      count_nonfinite_exposed_3d(mass_regions, source.t_init);
  report.invalid_exposed_alb_point_count =
      count_nonfinite_exposed_3d(mass_regions, source.alb);
  report.invalid_exposed_value_count =
      report.invalid_exposed_phb_point_count +
      report.invalid_exposed_mub_cell_count +
      report.invalid_exposed_ht_cell_count +
      report.invalid_exposed_pb_point_count +
      report.invalid_exposed_t_init_point_count +
      report.invalid_exposed_alb_point_count;
}

[[nodiscard]] NestResult validate_source_views(
    const Grid& child_grid,
    const ExposedBaseStateViews<const float>& source) noexcept {
  if (!valid_canonical_view(source.phb) || !valid_canonical_view(source.mub) ||
      !valid_canonical_view(source.pb) || !valid_canonical_view(source.t_init) ||
      !valid_canonical_view(source.alb) || !valid_canonical_view(source.ht)) {
    return result(
        NestStatus::invalid_contract,
        "base-state source staging requires non-null canonical explicit source views");
  }

  if (!matches_active_shape(source.phb, child_grid.w_shape()) ||
      !matches_active_shape(source.pb, child_grid.mass_shape()) ||
      !matches_active_shape(source.t_init, child_grid.mass_shape()) ||
      !matches_active_shape(source.alb, child_grid.mass_shape()) ||
      !matches_active_shape(source.mub, child_grid.surface_shape()) ||
      !matches_active_shape(source.ht, child_grid.surface_shape())) {
    return result(
        NestStatus::invalid_contract,
        "base-state source staging explicit source active shapes must match child shapes");
  }

  return ok_result();
}

[[nodiscard]] NestResult validate_remap_windows(
    const RemapPlan& remap_plan,
    const Grid& child_grid) noexcept {
  const auto& config = child_grid.config();
  if (!window_fits_active_shape(
          remap_plan.mass,
          HorizontalStagger::mass,
          config.mass_nx,
          config.mass_ny) ||
      !window_fits_active_shape(
          remap_plan.surface,
          HorizontalStagger::surface,
          config.mass_nx,
          config.mass_ny) ||
      !window_fits_active_shape(
          remap_plan.w_full,
          HorizontalStagger::w_full,
          config.mass_nx,
          config.mass_ny)) {
    return result(
        NestStatus::invalid_contract,
        "base-state source staging remap windows must fit child active extents");
  }
  return ok_result();
}

}  // namespace

BaseStateSourceStagingReport BaseStateSourceStagingProvider::stage(
    const Grid& child_grid,
    const RemapPlan& remap_plan,
    const ExposedBaseStateViews<const float>& explicit_source_views,
    const BaseStateSourceStagingOptions options) {
  clear_views();
  auto report = base_report(child_grid);

  if (!valid_grid_config(child_grid.config())) {
    report.result = result(
        NestStatus::invalid_configuration,
        "base-state source staging requires positive child shape and full_nz = mass_nz + 1");
    return report;
  }

  if (!remap_plan.ok()) {
    report.result = remap_plan.result;
    return report;
  }

  if (const auto validation = validate_remap_windows(remap_plan, child_grid);
      !validation.ok()) {
    report.result = validation;
    return report;
  }

  if (const auto validation = validate_source_views(child_grid, explicit_source_views);
      !validation.ok()) {
    report.result = validation;
    return report;
  }

  const auto mass_regions = exposed_regions_from_overlap(
      remap_plan.mass, child_grid.config().mass_nx, child_grid.config().mass_ny);
  const auto surface_regions = exposed_regions_from_overlap(
      remap_plan.surface, child_grid.config().mass_nx, child_grid.config().mass_ny);
  const auto w_full_regions = exposed_regions_from_overlap(
      remap_plan.w_full, child_grid.config().mass_nx, child_grid.config().mass_ny);

  populate_region_counts(
      report, child_grid, mass_regions, surface_regions, w_full_regions);
  populate_invalid_exposed_counts(
      report,
      explicit_source_views,
      mass_regions,
      surface_regions,
      w_full_regions);
  if (report.invalid_exposed_value_count != 0) {
    report.result = result(
        NestStatus::invalid_contract,
        "base-state source staging exposed source values must be finite");
    return report;
  }

  phb_ = FieldStorage3D<float>(child_grid.w_layout());
  mub_ = FieldStorage2D<float>(child_grid.surface_layout());
  pb_ = FieldStorage3D<float>(child_grid.mass_layout());
  t_init_ = FieldStorage3D<float>(child_grid.mass_layout());
  alb_ = FieldStorage3D<float>(child_grid.mass_layout());
  ht_ = FieldStorage2D<float>(child_grid.surface_layout());

  fill_storage(phb_, options.mask_value);
  fill_storage(mub_, options.mask_value);
  fill_storage(pb_, options.mask_value);
  fill_storage(t_init_, options.mask_value);
  fill_storage(alb_, options.mask_value);
  fill_storage(ht_, options.mask_value);

  copy_exposed_3d(w_full_regions, explicit_source_views.phb, phb_.view());
  copy_exposed_2d(surface_regions, explicit_source_views.mub, mub_.view());
  copy_exposed_2d(surface_regions, explicit_source_views.ht, ht_.view());
  copy_exposed_3d(mass_regions, explicit_source_views.pb, pb_.view());
  copy_exposed_3d(mass_regions, explicit_source_views.t_init, t_init_.view());
  copy_exposed_3d(mass_regions, explicit_source_views.alb, alb_.view());

  report.result = ok_result();
  report.owns_staging_buffers = true;
  report.allocated_buffers = true;
  return report;
}

void BaseStateSourceStagingProvider::clear_views() noexcept {
  phb_ = FieldStorage3D<float>{};
  mub_ = FieldStorage2D<float>{};
  pb_ = FieldStorage3D<float>{};
  t_init_ = FieldStorage3D<float>{};
  alb_ = FieldStorage3D<float>{};
  ht_ = FieldStorage2D<float>{};
}

ExposedBaseStateViews<const float> BaseStateSourceStagingProvider::views()
    const noexcept {
  return {phb_.view(), mub_.view(), pb_.view(), t_init_.view(), alb_.view(), ht_.view()};
}

}  // namespace tywrf::nest
