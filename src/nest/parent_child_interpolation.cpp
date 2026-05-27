#include "tywrf/nest/parent_child_interpolation.hpp"

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

[[nodiscard]] constexpr std::int32_t zero_based_parent_start(
    const std::int32_t parent_start,
    const IndexBase index_base) noexcept {
  return index_base == IndexBase::one_based ? parent_start - 1 : parent_start;
}

[[nodiscard]] constexpr bool positive_shape(
    const std::int32_t nx,
    const std::int32_t ny,
    const std::int32_t nz) noexcept {
  return nx > 0 && ny > 0 && nz > 0;
}

[[nodiscard]] constexpr bool region_fits_active_shape(
    const ExposedChildRegion& region,
    const std::int32_t nx,
    const std::int32_t ny,
    const std::int32_t nz) noexcept {
  return region.child_i_begin >= 0 && region.child_j_begin >= 0 &&
         region.extent_i >= 0 && region.extent_j >= 0 &&
         region.active_k_count >= 0 &&
         region.child_i_begin + region.extent_i <= nx &&
         region.child_j_begin + region.extent_j <= ny &&
         region.active_k_count <= nz;
}

[[nodiscard]] constexpr bool shape_matches_descriptor(
    const HorizontalDomainDescriptor& domain,
    const HorizontalStagger stagger,
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  switch (stagger) {
    case HorizontalStagger::mass:
    case HorizontalStagger::surface:
    case HorizontalStagger::w_full:
      return nx == domain.mass_nx && ny == domain.mass_ny;
    case HorizontalStagger::u:
      return nx == domain.mass_nx + 1 && ny == domain.mass_ny;
    case HorizontalStagger::v:
      return nx == domain.mass_nx && ny == domain.mass_ny + 1;
  }
  return false;
}

[[nodiscard]] constexpr bool field_plan_matches(
    const FieldStateExchangePlan& field,
    const StateExchangeField expected_field,
    const HorizontalStagger expected_stagger,
    const std::int32_t child_nx,
    const std::int32_t child_ny,
    const std::int32_t child_k_count,
    const bool three_dimensional) noexcept {
  return field.field == expected_field && field.stagger == expected_stagger &&
         field.active_nx == child_nx && field.active_ny == child_ny &&
         field.active_k_count == child_k_count &&
         field.three_dimensional == three_dimensional;
}

[[nodiscard]] constexpr bool supported_region_ownership(
    const FieldStateExchangePlan& field,
    const ParentChildInterpolationConfig& config) noexcept {
  if (field.owns_overlap && !config.allow_owned_overlap_regions) {
    return false;
  }
  if (field.owns_halo && !config.allow_owned_halo_regions) {
    return false;
  }
  for (std::uint8_t region_index = 0;
       region_index < field.exposed_region_count && region_index < field.exposed_regions.size();
       ++region_index) {
    const auto& region = field.exposed_regions[region_index];
    if (region.owns_overlap && !config.allow_owned_overlap_regions) {
      return false;
    }
    if (region.owns_halo && !config.allow_owned_halo_regions) {
      return false;
    }
  }
  return true;
}

void mark_written_ownership(
    const FieldStateExchangePlan& field,
    ParentChildInterpolationReport& report) noexcept {
  report.wrote_overlap = report.wrote_overlap || field.owns_overlap;
  report.wrote_halo = report.wrote_halo || field.owns_halo;
  for (std::uint8_t region_index = 0;
       region_index < field.exposed_region_count &&
       region_index < field.exposed_regions.size();
       ++region_index) {
    report.wrote_overlap =
        report.wrote_overlap || field.exposed_regions[region_index].owns_overlap;
    report.wrote_halo =
        report.wrote_halo || field.exposed_regions[region_index].owns_halo;
  }
}

template <typename ParentReal, typename ChildReal>
[[nodiscard]] constexpr NestResult validate_3d_field(
    const HorizontalDomainDescriptor& parent_domain,
    const HorizontalDomainDescriptor& child_domain,
    const HorizontalStagger stagger,
    const FieldView3D<ParentReal>& parent,
    const FieldView3D<ChildReal>& child) noexcept {
  const auto parent_nx = active_nx(parent);
  const auto parent_ny = active_ny(parent);
  const auto parent_nz = active_nz(parent);
  const auto child_nx = active_nx(child);
  const auto child_ny = active_ny(child);
  const auto child_nz = active_nz(child);

  if (!positive_shape(parent_nx, parent_ny, parent_nz) ||
      !positive_shape(child_nx, child_ny, child_nz)) {
    return result(NestStatus::invalid_contract, "parent-child interpolation field is empty");
  }
  if (parent_nz != child_nz) {
    return result(
        NestStatus::invalid_contract,
        "parent-child interpolation requires matching vertical active levels");
  }
  if (!shape_matches_descriptor(parent_domain, stagger, parent_nx, parent_ny) ||
      !shape_matches_descriptor(child_domain, stagger, child_nx, child_ny)) {
    return result(
        NestStatus::invalid_contract,
        "parent-child interpolation field shape does not match descriptor");
  }
  return ok_result();
}

template <typename ParentReal, typename ChildReal>
[[nodiscard]] constexpr NestResult validate_2d_field(
    const HorizontalDomainDescriptor& parent_domain,
    const HorizontalDomainDescriptor& child_domain,
    const HorizontalStagger stagger,
    const FieldView2D<ParentReal>& parent,
    const FieldView2D<ChildReal>& child) noexcept {
  const auto parent_nx = active_nx(parent);
  const auto parent_ny = active_ny(parent);
  const auto child_nx = active_nx(child);
  const auto child_ny = active_ny(child);

  if (!positive_shape(parent_nx, parent_ny, 1) ||
      !positive_shape(child_nx, child_ny, 1)) {
    return result(NestStatus::invalid_contract, "parent-child interpolation field is empty");
  }
  if (!shape_matches_descriptor(parent_domain, stagger, parent_nx, parent_ny) ||
      !shape_matches_descriptor(child_domain, stagger, child_nx, child_ny)) {
    return result(
        NestStatus::invalid_contract,
        "parent-child interpolation field shape does not match descriptor");
  }
  return ok_result();
}

[[nodiscard]] constexpr double parent_coordinate(
    const std::int32_t parent_start_zero_based,
    const std::int32_t child_active_index,
    const std::int32_t ratio) noexcept {
  return static_cast<double>(parent_start_zero_based) +
         static_cast<double>(child_active_index) / static_cast<double>(ratio);
}

[[nodiscard]] constexpr std::int32_t nearest_index(
    const double coordinate,
    const std::int32_t active_count) noexcept {
  const auto rounded = static_cast<std::int32_t>(coordinate + 0.5);
  return std::clamp(rounded, 0, active_count - 1);
}

template <typename Real>
[[nodiscard]] float sample_nearest_2d(
    const FieldView2D<Real>& field,
    const double x,
    const double y) noexcept {
  const auto i = field.halo.i_lower + nearest_index(x, active_nx(field));
  const auto j = field.halo.j_lower + nearest_index(y, active_ny(field));
  return static_cast<float>(field(i, j));
}

template <typename Real>
[[nodiscard]] float sample_bilinear_2d(
    const FieldView2D<Real>& field,
    const double x,
    const double y) noexcept {
  const auto nx = active_nx(field);
  const auto ny = active_ny(field);
  const auto x0 = std::clamp(static_cast<std::int32_t>(std::floor(x)), 0, nx - 1);
  const auto y0 = std::clamp(static_cast<std::int32_t>(std::floor(y)), 0, ny - 1);
  const auto x1 = std::min(x0 + 1, nx - 1);
  const auto y1 = std::min(y0 + 1, ny - 1);
  const auto ax = x1 == x0 ? 0.0 : x - static_cast<double>(x0);
  const auto ay = y1 == y0 ? 0.0 : y - static_cast<double>(y0);
  const auto i0 = field.halo.i_lower + x0;
  const auto i1 = field.halo.i_lower + x1;
  const auto j0 = field.halo.j_lower + y0;
  const auto j1 = field.halo.j_lower + y1;
  const auto f00 = static_cast<double>(field(i0, j0));
  const auto f10 = static_cast<double>(field(i1, j0));
  const auto f01 = static_cast<double>(field(i0, j1));
  const auto f11 = static_cast<double>(field(i1, j1));
  const auto fx0 = f00 * (1.0 - ax) + f10 * ax;
  const auto fx1 = f01 * (1.0 - ax) + f11 * ax;
  return static_cast<float>(fx0 * (1.0 - ay) + fx1 * ay);
}

template <typename Real>
[[nodiscard]] float sample_2d(
    const FieldView2D<Real>& field,
    const double x,
    const double y,
    const ParentChildInterpolationMethod method) noexcept {
  if (method == ParentChildInterpolationMethod::nearest_neighbor) {
    return sample_nearest_2d(field, x, y);
  }
  return sample_bilinear_2d(field, x, y);
}

template <typename Real>
[[nodiscard]] float sample_nearest_3d(
    const FieldView3D<Real>& field,
    const double x,
    const double y,
    const std::int32_t active_k) noexcept {
  const auto i = field.halo.i_lower + nearest_index(x, active_nx(field));
  const auto j = field.halo.j_lower + nearest_index(y, active_ny(field));
  const auto k = field.halo.k_lower + active_k;
  return static_cast<float>(field(i, j, k));
}

template <typename Real>
[[nodiscard]] float sample_bilinear_3d(
    const FieldView3D<Real>& field,
    const double x,
    const double y,
    const std::int32_t active_k) noexcept {
  const auto nx = active_nx(field);
  const auto ny = active_ny(field);
  const auto x0 = std::clamp(static_cast<std::int32_t>(std::floor(x)), 0, nx - 1);
  const auto y0 = std::clamp(static_cast<std::int32_t>(std::floor(y)), 0, ny - 1);
  const auto x1 = std::min(x0 + 1, nx - 1);
  const auto y1 = std::min(y0 + 1, ny - 1);
  const auto ax = x1 == x0 ? 0.0 : x - static_cast<double>(x0);
  const auto ay = y1 == y0 ? 0.0 : y - static_cast<double>(y0);
  const auto i0 = field.halo.i_lower + x0;
  const auto i1 = field.halo.i_lower + x1;
  const auto j0 = field.halo.j_lower + y0;
  const auto j1 = field.halo.j_lower + y1;
  const auto k = field.halo.k_lower + active_k;
  const auto f00 = static_cast<double>(field(i0, j0, k));
  const auto f10 = static_cast<double>(field(i1, j0, k));
  const auto f01 = static_cast<double>(field(i0, j1, k));
  const auto f11 = static_cast<double>(field(i1, j1, k));
  const auto fx0 = f00 * (1.0 - ax) + f10 * ax;
  const auto fx1 = f01 * (1.0 - ax) + f11 * ax;
  return static_cast<float>(fx0 * (1.0 - ay) + fx1 * ay);
}

template <typename Real>
[[nodiscard]] float sample_3d(
    const FieldView3D<Real>& field,
    const double x,
    const double y,
    const std::int32_t active_k,
    const ParentChildInterpolationMethod method) noexcept {
  if (method == ParentChildInterpolationMethod::nearest_neighbor) {
    return sample_nearest_3d(field, x, y, active_k);
  }
  return sample_bilinear_3d(field, x, y, active_k);
}

template <typename ParentReal, typename ChildReal>
[[nodiscard]] NestResult interpolate_regions_2d(
    const FieldStateExchangePlan& field_plan,
    const StateExchangeField expected_field,
    const HorizontalStagger expected_stagger,
    const ParentChildInterpolationConfig& config,
    const std::int32_t ratio,
    const std::int32_t parent_i_start_zero,
    const std::int32_t parent_j_start_zero,
    const FieldView2D<ParentReal>& parent,
    const FieldView2D<ChildReal>& child,
    ParentChildInterpolationReport& report) noexcept {
  if (!field_plan_matches(
          field_plan,
          expected_field,
          expected_stagger,
          active_nx(child),
          active_ny(child),
          1,
          false)) {
    return result(
        NestStatus::invalid_contract,
        "parent-child interpolation field plan does not match child field");
  }
  if (!supported_region_ownership(field_plan, config)) {
    return result(
        NestStatus::invalid_contract,
        "parent-child interpolation plan owns overlap or halo without permission");
  }
  mark_written_ownership(field_plan, report);

  std::uint64_t point_count = 0;
  for (std::uint8_t region_index = 0;
       region_index < field_plan.exposed_region_count &&
       region_index < field_plan.exposed_regions.size();
       ++region_index) {
    const auto& region = field_plan.exposed_regions[region_index];
    if (region.stagger != expected_stagger ||
        region.three_dimensional || !region_fits_active_shape(region, active_nx(child), active_ny(child), 1)) {
      return result(
          NestStatus::invalid_contract,
          "parent-child interpolation region is outside child active field");
    }
    const auto child_i0 = child.halo.i_lower + region.child_i_begin;
    const auto child_j0 = child.halo.j_lower + region.child_j_begin;
    for (std::int32_t j = 0; j < region.extent_j; ++j) {
      const auto child_active_j = region.child_j_begin + j;
      const auto y = parent_coordinate(parent_j_start_zero, child_active_j, ratio);
      for (std::int32_t i = 0; i < region.extent_i; ++i) {
        const auto child_active_i = region.child_i_begin + i;
        const auto x = parent_coordinate(parent_i_start_zero, child_active_i, ratio);
        child(child_i0 + i, child_j0 + j) = sample_2d(parent, x, y, config.method);
        ++point_count;
      }
    }
    ++report.interpolated_region_count;
    report.interpolated_horizontal_cell_count += region.horizontal_cell_count();
  }

  if (point_count > 0) {
    ++report.interpolated_field_count;
    report.interpolated_point_count += point_count;
  }
  return ok_result();
}

template <typename ParentReal, typename ChildReal>
[[nodiscard]] NestResult interpolate_regions_3d(
    const FieldStateExchangePlan& field_plan,
    const StateExchangeField expected_field,
    const HorizontalStagger expected_stagger,
    const ParentChildInterpolationConfig& config,
    const std::int32_t ratio,
    const std::int32_t parent_i_start_zero,
    const std::int32_t parent_j_start_zero,
    const FieldView3D<ParentReal>& parent,
    const FieldView3D<ChildReal>& child,
    ParentChildInterpolationReport& report) noexcept {
  if (!field_plan_matches(
          field_plan,
          expected_field,
          expected_stagger,
          active_nx(child),
          active_ny(child),
          active_nz(child),
          true)) {
    return result(
        NestStatus::invalid_contract,
        "parent-child interpolation field plan does not match child field");
  }
  if (!supported_region_ownership(field_plan, config)) {
    return result(
        NestStatus::invalid_contract,
        "parent-child interpolation plan owns overlap or halo without permission");
  }
  mark_written_ownership(field_plan, report);

  std::uint64_t point_count = 0;
  for (std::uint8_t region_index = 0;
       region_index < field_plan.exposed_region_count &&
       region_index < field_plan.exposed_regions.size();
       ++region_index) {
    const auto& region = field_plan.exposed_regions[region_index];
    if (region.stagger != expected_stagger || !region.three_dimensional ||
        !region_fits_active_shape(
            region, active_nx(child), active_ny(child), active_nz(child))) {
      return result(
          NestStatus::invalid_contract,
          "parent-child interpolation region is outside child active field");
    }
    const auto child_i0 = child.halo.i_lower + region.child_i_begin;
    const auto child_j0 = child.halo.j_lower + region.child_j_begin;
    const auto child_k0 = child.halo.k_lower;
    for (std::int32_t j = 0; j < region.extent_j; ++j) {
      const auto child_active_j = region.child_j_begin + j;
      const auto y = parent_coordinate(parent_j_start_zero, child_active_j, ratio);
      for (std::int32_t k = 0; k < region.active_k_count; ++k) {
        for (std::int32_t i = 0; i < region.extent_i; ++i) {
          const auto child_active_i = region.child_i_begin + i;
          const auto x = parent_coordinate(parent_i_start_zero, child_active_i, ratio);
          child(child_i0 + i, child_j0 + j, child_k0 + k) =
              sample_3d(parent, x, y, k, config.method);
          ++point_count;
        }
      }
    }
    ++report.interpolated_region_count;
    report.interpolated_horizontal_cell_count += region.horizontal_cell_count();
  }

  if (point_count > 0) {
    ++report.interpolated_field_count;
    report.interpolated_point_count += point_count;
  }
  return ok_result();
}

[[nodiscard]] NestResult validate_state_shapes(
    const ParentChildDescriptor& descriptor,
    const StateView<const float>& parent,
    const StateView<float>& child) noexcept {
  if (const auto validation = validate_3d_field(
          descriptor.parent, descriptor.child, HorizontalStagger::u, parent.u, child.u);
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_3d_field(
          descriptor.parent, descriptor.child, HorizontalStagger::v, parent.v, child.v);
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_2d_field(
          descriptor.parent,
          descriptor.child,
          HorizontalStagger::surface,
          parent.mu,
          child.mu);
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_3d_field(
          descriptor.parent,
          descriptor.child,
          HorizontalStagger::mass,
          parent.qvapor,
          child.qvapor);
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_3d_field(
          descriptor.parent,
          descriptor.child,
          HorizontalStagger::mass,
          parent.t,
          child.t);
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_3d_field(
          descriptor.parent,
          descriptor.child,
          HorizontalStagger::w_full,
          parent.ph,
          child.ph);
      !validation.ok()) {
    return validation;
  }
  return ok_result();
}

}  // namespace

ParentChildInterpolationReport interpolate_parent_to_exposed_child(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& target_position,
    const StateExchangePlan& exchange_plan,
    const StateView<const float>& parent_state,
    const StateView<float>& child_state,
    const ParentChildInterpolationConfig config) noexcept {
  ParentChildInterpolationReport report{};
  report.result = ok_result();
  report.method = config.method;

  if (const auto descriptor_status = validate_parent_child_descriptor(descriptor);
      !descriptor_status.ok()) {
    report.result = descriptor_status;
    return report;
  }
  if (const auto position_status =
          validate_parent_child_position(descriptor, target_position);
      !position_status.ok()) {
    report.result = position_status;
    return report;
  }
  if (!exchange_plan.ok()) {
    report.result = exchange_plan.result;
    return report;
  }
  if (exchange_plan.operation != ExchangeOperation::parent_to_child_interpolation) {
    report.result = result(
        NestStatus::invalid_contract,
        "state exchange plan is not parent-to-child interpolation");
    return report;
  }
  if (const auto shape_status = validate_state_shapes(descriptor, parent_state, child_state);
      !shape_status.ok()) {
    report.result = shape_status;
    return report;
  }

  const auto field_count = std::min<std::uint8_t>(
      exchange_plan.field_count, static_cast<std::uint8_t>(exchange_plan.fields.size()));
  report.requested_field_count = field_count;
  const auto parent_i_start_zero =
      zero_based_parent_start(target_position.i_parent_start, target_position.index_base);
  const auto parent_j_start_zero =
      zero_based_parent_start(target_position.j_parent_start, target_position.index_base);

  for (std::uint8_t field_index = 0; field_index < field_count; ++field_index) {
    const auto& field = exchange_plan.fields[field_index];
    NestResult field_result = ok_result();
    switch (field.field) {
      case StateExchangeField::u:
        field_result = interpolate_regions_3d(
            field,
            StateExchangeField::u,
            HorizontalStagger::u,
            config,
            descriptor.parent_grid_ratio,
            parent_i_start_zero,
            parent_j_start_zero,
            parent_state.u,
            child_state.u,
            report);
        break;
      case StateExchangeField::v:
        field_result = interpolate_regions_3d(
            field,
            StateExchangeField::v,
            HorizontalStagger::v,
            config,
            descriptor.parent_grid_ratio,
            parent_i_start_zero,
            parent_j_start_zero,
            parent_state.v,
            child_state.v,
            report);
        break;
      case StateExchangeField::mu:
        field_result = interpolate_regions_2d(
            field,
            StateExchangeField::mu,
            HorizontalStagger::surface,
            config,
            descriptor.parent_grid_ratio,
            parent_i_start_zero,
            parent_j_start_zero,
            parent_state.mu,
            child_state.mu,
            report);
        break;
      case StateExchangeField::qvapor:
        field_result = interpolate_regions_3d(
            field,
            StateExchangeField::qvapor,
            HorizontalStagger::mass,
            config,
            descriptor.parent_grid_ratio,
            parent_i_start_zero,
            parent_j_start_zero,
            parent_state.qvapor,
            child_state.qvapor,
            report);
        break;
      case StateExchangeField::t:
        field_result = interpolate_regions_3d(
            field,
            StateExchangeField::t,
            HorizontalStagger::mass,
            config,
            descriptor.parent_grid_ratio,
            parent_i_start_zero,
            parent_j_start_zero,
            parent_state.t,
            child_state.t,
            report);
        break;
      case StateExchangeField::ph:
        field_result = interpolate_regions_3d(
            field,
            StateExchangeField::ph,
            HorizontalStagger::w_full,
            config,
            descriptor.parent_grid_ratio,
            parent_i_start_zero,
            parent_j_start_zero,
            parent_state.ph,
            child_state.ph,
            report);
        break;
    }

    if (!field_result.ok()) {
      report.result = field_result;
      return report;
    }
  }

  report.result = ok_result();
  return report;
}

}  // namespace tywrf::nest
