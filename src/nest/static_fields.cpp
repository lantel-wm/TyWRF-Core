#include "tywrf/nest/static_fields.hpp"

#include <algorithm>
#include <cmath>
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
[[nodiscard]] constexpr bool positive_shape(const FieldView2D<Real>& field) noexcept {
  return active_nx(field) > 0 && active_ny(field) > 0;
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_active_shape(
    const FieldView2D<LhsReal>& lhs,
    const FieldView2D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

template <typename Real>
[[nodiscard]] constexpr bool matches_domain(
    const FieldView2D<Real>& field,
    const HorizontalDomainDescriptor& domain) noexcept {
  return active_nx(field) == domain.mass_nx && active_ny(field) == domain.mass_ny;
}

[[nodiscard]] constexpr bool surface_window_fits(
    const RemapWindow& window,
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return window.stagger == HorizontalStagger::surface && window.old_i_begin >= 0 &&
         window.old_j_begin >= 0 && window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0 &&
         window.old_i_begin + window.extent_i <= nx &&
         window.old_j_begin + window.extent_j <= ny &&
         window.new_i_begin + window.extent_i <= nx &&
         window.new_j_begin + window.extent_j <= ny;
}

[[nodiscard]] constexpr bool in_new_overlap_window(
    const RemapWindow& window,
    const std::int32_t i,
    const std::int32_t j) noexcept {
  return i >= window.new_i_begin && i < window.new_i_begin + window.extent_i &&
         j >= window.new_j_begin && j < window.new_j_begin + window.extent_j;
}

[[nodiscard]] constexpr std::int32_t zero_based_parent_start(
    const std::int32_t parent_start,
    const IndexBase index_base) noexcept {
  return index_base == IndexBase::one_based ? parent_start - 1 : parent_start;
}

template <typename Real>
[[nodiscard]] float value_at(
    const FieldView2D<Real>& field,
    const std::int32_t active_i,
    const std::int32_t active_j) noexcept {
  return static_cast<float>(
      field(field.halo.i_lower + active_i, field.halo.j_lower + active_j));
}

template <typename Real>
[[nodiscard]] float local_linear_extrapolated_2d(
    const FieldView2D<Real>& field,
    const double x,
    const double y) noexcept {
  const auto nx = active_nx(field);
  const auto ny = active_ny(field);
  const auto ix = std::clamp(static_cast<std::int32_t>(std::floor(x)), 0, nx - 1);
  const auto iy = std::clamp(static_cast<std::int32_t>(std::floor(y)), 0, ny - 1);
  const auto base = static_cast<double>(value_at(field, ix, iy));

  double grad_x = 0.0;
  if (nx > 1) {
    if (ix == 0) {
      grad_x = static_cast<double>(value_at(field, 1, iy)) -
               static_cast<double>(value_at(field, 0, iy));
    } else if (ix == nx - 1) {
      grad_x = static_cast<double>(value_at(field, nx - 1, iy)) -
               static_cast<double>(value_at(field, nx - 2, iy));
    } else {
      grad_x = 0.5 *
               (static_cast<double>(value_at(field, ix + 1, iy)) -
                static_cast<double>(value_at(field, ix - 1, iy)));
    }
  }

  double grad_y = 0.0;
  if (ny > 1) {
    if (iy == 0) {
      grad_y = static_cast<double>(value_at(field, ix, 1)) -
               static_cast<double>(value_at(field, ix, 0));
    } else if (iy == ny - 1) {
      grad_y = static_cast<double>(value_at(field, ix, ny - 1)) -
               static_cast<double>(value_at(field, ix, ny - 2));
    } else {
      grad_y = 0.5 *
               (static_cast<double>(value_at(field, ix, iy + 1)) -
                static_cast<double>(value_at(field, ix, iy - 1)));
    }
  }

  return static_cast<float>(
      base + (x - static_cast<double>(ix)) * grad_x +
      (y - static_cast<double>(iy)) * grad_y);
}

template <typename Real>
[[nodiscard]] float bilinear_clamped_2d(
    const FieldView2D<Real>& field,
    const double x,
    const double y) noexcept {
  const auto nx = active_nx(field);
  const auto ny = active_ny(field);
  const auto clamped_x = std::clamp(x, 0.0, static_cast<double>(nx - 1));
  const auto clamped_y = std::clamp(y, 0.0, static_cast<double>(ny - 1));
  const auto x0 = static_cast<std::int32_t>(std::floor(clamped_x));
  const auto y0 = static_cast<std::int32_t>(std::floor(clamped_y));
  const auto x1 = std::min(x0 + 1, nx - 1);
  const auto y1 = std::min(y0 + 1, ny - 1);
  const auto ax = x1 == x0 ? 0.0 : clamped_x - static_cast<double>(x0);
  const auto ay = y1 == y0 ? 0.0 : clamped_y - static_cast<double>(y0);
  const auto f00 = static_cast<double>(value_at(field, x0, y0));
  const auto f10 = static_cast<double>(value_at(field, x1, y0));
  const auto f01 = static_cast<double>(value_at(field, x0, y1));
  const auto f11 = static_cast<double>(value_at(field, x1, y1));
  const auto fx0 = f00 * (1.0 - ax) + f10 * ax;
  const auto fx1 = f01 * (1.0 - ax) + f11 * ax;
  return static_cast<float>(fx0 * (1.0 - ay) + fx1 * ay);
}

[[nodiscard]] NestResult validate_inputs(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& target_position,
    const RemapPlan& remap_plan,
    const FieldView2D<const float>& d02_start_xlat,
    const FieldView2D<const float>& d02_start_xlong,
    const FieldView2D<const float>& d02_start_hgt,
    const FieldView2D<const float>& d01_start_hgt,
    const FieldView2D<float>& output_xlat,
    const FieldView2D<float>& output_xlong,
    const FieldView2D<float>& output_hgt) noexcept {
  if (const auto descriptor_status = validate_parent_child_descriptor(descriptor);
      !descriptor_status.ok()) {
    return descriptor_status;
  }
  if (const auto position_status =
          validate_parent_child_position(descriptor, target_position);
      !position_status.ok()) {
    return position_status;
  }
  if (!remap_plan.ok()) {
    return remap_plan.result;
  }
  if (!positive_shape(d02_start_xlat) || !positive_shape(d02_start_xlong) ||
      !positive_shape(d02_start_hgt) || !positive_shape(d01_start_hgt) ||
      !positive_shape(output_xlat) || !positive_shape(output_xlong) ||
      !positive_shape(output_hgt)) {
    return result(NestStatus::invalid_contract, "static refresh fields must be non-empty");
  }
  if (!matches_domain(d01_start_hgt, descriptor.parent) ||
      !matches_domain(d02_start_xlat, descriptor.child) ||
      !matches_domain(d02_start_xlong, descriptor.child) ||
      !matches_domain(d02_start_hgt, descriptor.child)) {
    return result(
        NestStatus::invalid_contract,
        "static refresh input field shape does not match parent-child descriptor");
  }
  if (!same_active_shape(d02_start_xlat, d02_start_xlong) ||
      !same_active_shape(d02_start_xlat, d02_start_hgt) ||
      !same_active_shape(d02_start_xlat, output_xlat) ||
      !same_active_shape(d02_start_xlat, output_xlong) ||
      !same_active_shape(d02_start_xlat, output_hgt)) {
    return result(
        NestStatus::invalid_contract,
        "static refresh d02 input and output fields must share active shape");
  }
  if (!surface_window_fits(
          remap_plan.surface, active_nx(d02_start_xlat), active_ny(d02_start_xlat))) {
    return result(
        NestStatus::invalid_contract,
        "static refresh surface remap window is outside d02 active extents");
  }
  return ok_result();
}

}  // namespace

MovingNestStaticRefreshReport refresh_moving_nest_static_fields(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& target_position,
    const RemapPlan& remap_plan,
    const FieldView2D<const float>& d02_start_xlat,
    const FieldView2D<const float>& d02_start_xlong,
    const FieldView2D<const float>& d02_start_hgt,
    const FieldView2D<const float>& d01_start_hgt,
    const FieldView2D<float>& output_xlat,
    const FieldView2D<float>& output_xlong,
    const FieldView2D<float>& output_hgt) noexcept {
  MovingNestStaticRefreshReport report{};
  report.result = validate_inputs(
      descriptor,
      target_position,
      remap_plan,
      d02_start_xlat,
      d02_start_xlong,
      d02_start_hgt,
      d01_start_hgt,
      output_xlat,
      output_xlong,
      output_hgt);
  if (!report.ok()) {
    return report;
  }

  const auto& window = remap_plan.surface;
  const auto nx = active_nx(output_xlat);
  const auto ny = active_ny(output_xlat);
  const auto parent_i_start_zero =
      zero_based_parent_start(target_position.i_parent_start, target_position.index_base);
  const auto parent_j_start_zero =
      zero_based_parent_start(target_position.j_parent_start, target_position.index_base);

  for (std::int32_t j = 0; j < ny; ++j) {
    for (std::int32_t i = 0; i < nx; ++i) {
      const auto out_i = output_xlat.halo.i_lower + i;
      const auto out_j = output_xlat.halo.j_lower + j;
      const auto old_i = i + window.old_i_offset_from_new;
      const auto old_j = j + window.old_j_offset_from_new;
      if (in_new_overlap_window(window, i, j)) {
        output_xlat(out_i, out_j) = value_at(d02_start_xlat, old_i, old_j);
        output_xlong(
            output_xlong.halo.i_lower + i,
            output_xlong.halo.j_lower + j) = value_at(d02_start_xlong, old_i, old_j);
        output_hgt(output_hgt.halo.i_lower + i, output_hgt.halo.j_lower + j) =
            value_at(d02_start_hgt, old_i, old_j);
        ++report.overlap_cell_count;
        continue;
      }

      output_xlat(out_i, out_j) =
          local_linear_extrapolated_2d(d02_start_xlat, old_i, old_j);
      output_xlong(output_xlong.halo.i_lower + i, output_xlong.halo.j_lower + j) =
          local_linear_extrapolated_2d(d02_start_xlong, old_i, old_j);

      const auto parent_x = static_cast<double>(parent_i_start_zero) +
                            static_cast<double>(i) /
                                static_cast<double>(descriptor.parent_grid_ratio);
      const auto parent_y = static_cast<double>(parent_j_start_zero) +
                            static_cast<double>(j) /
                                static_cast<double>(descriptor.parent_grid_ratio);
      output_hgt(output_hgt.halo.i_lower + i, output_hgt.halo.j_lower + j) =
          bilinear_clamped_2d(d01_start_hgt, parent_x, parent_y);
      ++report.exposed_cell_count;
    }
  }

  report.coordinate_extrapolated_cell_count = report.exposed_cell_count;
  report.parent_hgt_interpolated_cell_count = report.exposed_cell_count;
  report.copied_overlap = report.overlap_cell_count > 0;
  report.filled_exposed_coordinates = report.coordinate_extrapolated_cell_count > 0;
  report.filled_exposed_hgt = report.parent_hgt_interpolated_cell_count > 0;
  report.uses_reference_end = false;
  report.result = ok_result();
  return report;
}

}  // namespace tywrf::nest
