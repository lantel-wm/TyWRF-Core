#include "tywrf/nest/time_level_sync.hpp"

#include <cstdint>

namespace tywrf::nest {
namespace {

[[nodiscard]] constexpr NestResult ok_result() noexcept {
  return {NestStatus::ok, "ok"};
}

[[nodiscard]] constexpr NestResult result(
    const NestStatus status,
    const char* message) noexcept {
  return {status, message};
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(const FieldView2D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(const FieldView2D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(const FieldView3D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(const FieldView3D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nz(const FieldView3D<Real>& field) noexcept {
  return field.nz - field.halo.k_lower - field.halo.k_upper;
}

template <typename Real>
[[nodiscard]] constexpr bool valid_view(const FieldView2D<Real>& field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 &&
         field.stride_i == 1 && field.stride_j > 0 &&
         field.halo.i_lower >= 0 && field.halo.i_upper >= 0 &&
         field.halo.j_lower >= 0 && field.halo.j_upper >= 0 &&
         active_nx(field) > 0 && active_ny(field) > 0;
}

template <typename Real>
[[nodiscard]] constexpr bool valid_view(const FieldView3D<Real>& field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 &&
         field.nz > 0 && field.stride_i == 1 && field.stride_k > 0 &&
         field.stride_j > 0 && field.halo.i_lower >= 0 &&
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

template <typename SrcReal, typename DstReal>
[[nodiscard]] constexpr bool window_fits_pair(
    const RemapWindow& window,
    const HorizontalStagger expected_stagger,
    const FieldView2D<SrcReal>& level2,
    const FieldView2D<DstReal>& level1) noexcept {
  return window.stagger == expected_stagger && valid_view(level2) &&
         valid_view(level1) && same_active_shape(level2, level1) &&
         window.old_i_begin >= 0 && window.old_j_begin >= 0 &&
         window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0 &&
         window.old_i_begin + window.extent_i <= active_nx(level2) &&
         window.old_j_begin + window.extent_j <= active_ny(level2) &&
         window.new_i_begin + window.extent_i <= active_nx(level1) &&
         window.new_j_begin + window.extent_j <= active_ny(level1);
}

template <typename SrcReal, typename DstReal>
[[nodiscard]] constexpr bool window_fits_pair(
    const RemapWindow& window,
    const HorizontalStagger expected_stagger,
    const FieldView3D<SrcReal>& level2,
    const FieldView3D<DstReal>& level1) noexcept {
  return window.stagger == expected_stagger && valid_view(level2) &&
         valid_view(level1) && same_active_shape(level2, level1) &&
         window.old_i_begin >= 0 && window.old_j_begin >= 0 &&
         window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0 &&
         window.old_i_begin + window.extent_i <= active_nx(level2) &&
         window.old_j_begin + window.extent_j <= active_ny(level2) &&
         window.new_i_begin + window.extent_i <= active_nx(level1) &&
         window.new_j_begin + window.extent_j <= active_ny(level1);
}

[[nodiscard]] constexpr bool in_new_overlap_window(
    const RemapWindow& window,
    const std::int32_t i,
    const std::int32_t j) noexcept {
  return i >= window.new_i_begin && i < window.new_i_begin + window.extent_i &&
         j >= window.new_j_begin && j < window.new_j_begin + window.extent_j;
}

template <typename Real>
[[nodiscard]] std::uint64_t exposed_horizontal_cell_count(
    const RemapWindow& window,
    const FieldView2D<Real>& field) noexcept {
  return static_cast<std::uint64_t>(active_nx(field)) *
             static_cast<std::uint64_t>(active_ny(field)) -
         static_cast<std::uint64_t>(window.extent_i) *
             static_cast<std::uint64_t>(window.extent_j);
}

template <typename Real>
[[nodiscard]] std::uint64_t exposed_horizontal_cell_count(
    const RemapWindow& window,
    const FieldView3D<Real>& field) noexcept {
  return static_cast<std::uint64_t>(active_nx(field)) *
             static_cast<std::uint64_t>(active_ny(field)) -
         static_cast<std::uint64_t>(window.extent_i) *
             static_cast<std::uint64_t>(window.extent_j);
}

std::uint64_t copy_exposed_2d(
    const RemapWindow& window,
    const FieldView2D<const float>& level2,
    const FieldView2D<float>& level1) noexcept {
  const auto exposed_cells = exposed_horizontal_cell_count(window, level1);
  if (exposed_cells == 0) {
    return 0;
  }

  const auto src_i0 = level2.halo.i_lower;
  const auto src_j0 = level2.halo.j_lower;
  const auto dst_i0 = level1.halo.i_lower;
  const auto dst_j0 = level1.halo.j_lower;
  const auto nx = active_nx(level1);
  const auto ny = active_ny(level1);
  std::uint64_t copied_points = 0;

  for (std::int32_t j = 0; j < ny; ++j) {
    for (std::int32_t i = 0; i < nx; ++i) {
      if (in_new_overlap_window(window, i, j)) {
        continue;
      }
      level1(dst_i0 + i, dst_j0 + j) = level2(src_i0 + i, src_j0 + j);
      ++copied_points;
    }
  }

  return copied_points;
}

std::uint64_t copy_exposed_3d(
    const RemapWindow& window,
    const FieldView3D<const float>& level2,
    const FieldView3D<float>& level1) noexcept {
  const auto exposed_cells = exposed_horizontal_cell_count(window, level1);
  if (exposed_cells == 0) {
    return 0;
  }

  const auto src_i0 = level2.halo.i_lower;
  const auto src_j0 = level2.halo.j_lower;
  const auto src_k0 = level2.halo.k_lower;
  const auto dst_i0 = level1.halo.i_lower;
  const auto dst_j0 = level1.halo.j_lower;
  const auto dst_k0 = level1.halo.k_lower;
  const auto nx = active_nx(level1);
  const auto ny = active_ny(level1);
  const auto nz = active_nz(level1);
  std::uint64_t copied_points = 0;

  for (std::int32_t j = 0; j < ny; ++j) {
    for (std::int32_t k = 0; k < nz; ++k) {
      for (std::int32_t i = 0; i < nx; ++i) {
        if (in_new_overlap_window(window, i, j)) {
          continue;
        }
        level1(dst_i0 + i, dst_j0 + j, dst_k0 + k) =
            level2(src_i0 + i, src_j0 + j, src_k0 + k);
        ++copied_points;
      }
    }
  }

  return copied_points;
}

void add_field_points(
    TimeLevelSyncReport& report,
    std::uint64_t& field_point_count,
    const std::uint64_t copied_points) noexcept {
  field_point_count = copied_points;
  if (copied_points == 0) {
    return;
  }

  ++report.copied_field_count;
  report.copied_point_count += copied_points;
}

[[nodiscard]] NestResult validate_contract(
    const RemapPlan& plan,
    const PostStartDomainLevel2Views& level2,
    const PostStartDomainLevel1Views& level1,
    const OptionalTimeLevelSyncTke& optional_tke) noexcept {
  if (!plan.ok()) {
    return plan.result;
  }

  if (!window_fits_pair(plan.u, HorizontalStagger::u, level2.u, level1.u) ||
      !window_fits_pair(plan.v, HorizontalStagger::v, level2.v, level1.v) ||
      !window_fits_pair(plan.mass, HorizontalStagger::mass, level2.t, level1.t) ||
      !window_fits_pair(plan.w_full, HorizontalStagger::w_full, level2.w, level1.w) ||
      !window_fits_pair(plan.w_full, HorizontalStagger::w_full, level2.ph, level1.ph) ||
      !window_fits_pair(
          plan.surface, HorizontalStagger::surface, level2.mu, level1.mu)) {
    return result(
        NestStatus::invalid_contract,
        "time-level sync views must match remap windows and active shapes");
  }

  if (optional_tke.present &&
      (!window_fits_pair(
           plan.mass,
           HorizontalStagger::mass,
           optional_tke.level2,
           optional_tke.level1) ||
       !same_active_shape(level2.t, optional_tke.level2) ||
       !same_active_shape(level1.t, optional_tke.level1))) {
    return result(
        NestStatus::invalid_contract,
        "optional TKE time-level sync view must be mass shaped");
  }

  return ok_result();
}

}  // namespace

TimeLevelSyncReport copy_post_start_domain_level2_to_level1_exposed_cells(
    const RemapPlan& plan,
    const PostStartDomainLevel2Views& level2,
    const PostStartDomainLevel1Views& level1,
    const OptionalTimeLevelSyncTke optional_tke) noexcept {
  TimeLevelSyncReport report{};
  report.result = ok_result();

  if (const auto validation = validate_contract(plan, level2, level1, optional_tke);
      !validation.ok()) {
    report.result = validation;
    return report;
  }

  add_field_points(report, report.u_point_count, copy_exposed_3d(plan.u, level2.u, level1.u));
  add_field_points(report, report.v_point_count, copy_exposed_3d(plan.v, level2.v, level1.v));
  add_field_points(
      report, report.t_point_count, copy_exposed_3d(plan.mass, level2.t, level1.t));
  add_field_points(
      report, report.w_point_count, copy_exposed_3d(plan.w_full, level2.w, level1.w));
  add_field_points(
      report, report.ph_point_count, copy_exposed_3d(plan.w_full, level2.ph, level1.ph));
  add_field_points(
      report, report.mu_point_count, copy_exposed_2d(plan.surface, level2.mu, level1.mu));

  if (optional_tke.present) {
    add_field_points(
        report,
        report.tke_point_count,
        copy_exposed_3d(plan.mass, optional_tke.level2, optional_tke.level1));
    report.copied_optional_tke = report.tke_point_count > 0;
  }

  return report;
}

}  // namespace tywrf::nest
