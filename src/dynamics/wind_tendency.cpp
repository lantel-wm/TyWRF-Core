#include "tywrf/dynamics/wind_tendency.hpp"

#include "tywrf/dynamics/tendency.hpp"

#include <cstddef>
#include <cstdint>

namespace tywrf::dynamics {
namespace {

template <typename Real>
[[nodiscard]] bool positive(const Real value) noexcept {
  return value > static_cast<Real>(0);
}

template <typename Real>
[[nodiscard]] WindTendencyReport make_report(
    const WindTendencyConfig<Real> config,
    const WindTendencyStatus status) noexcept {
  WindTendencyReport report{};
  report.status = status;
  report.diagnostic_only = config.diagnostic_only;
  report.gate_candidate = config.gate_candidate;
  report.validation_gate_evidence = false;
  report.horizontal_advection_enabled = config.enable_horizontal_advection;
  return report;
}

template <typename Real>
[[nodiscard]] bool same_layout(
    const FieldView3D<Real> target,
    const FieldView3D<const Real> source) noexcept {
  return target.nx == source.nx && target.ny == source.ny &&
         target.nz == source.nz && target.stride_i == source.stride_i &&
         target.stride_k == source.stride_k &&
         target.stride_j == source.stride_j &&
         target.halo.i_lower == source.halo.i_lower &&
         target.halo.i_upper == source.halo.i_upper &&
         target.halo.j_lower == source.halo.j_lower &&
         target.halo.j_upper == source.halo.j_upper &&
         target.halo.k_lower == source.halo.k_lower &&
         target.halo.k_upper == source.halo.k_upper;
}

template <typename Real>
[[nodiscard]] bool has_one_cell_horizontal_halo(
    const FieldView3D<Real> view) noexcept {
  return view.halo.i_lower >= 1 && view.halo.i_upper >= 1 &&
         view.halo.j_lower >= 1 && view.halo.j_upper >= 1;
}

template <typename Real>
[[nodiscard]] bool has_one_cell_horizontal_halo(
    const FieldView3D<const Real> view) noexcept {
  return view.halo.i_lower >= 1 && view.halo.i_upper >= 1 &&
         view.halo.j_lower >= 1 && view.halo.j_upper >= 1;
}

template <typename Real>
[[nodiscard]] WindTendencyStatus validate_field(
    const WindTendencyFieldViews<Real> field) noexcept {
  if (field.target.data == nullptr) {
    return WindTendencyStatus::null_target;
  }
  if (field.source.data == nullptr || field.advect_x.data == nullptr ||
      field.advect_y.data == nullptr) {
    return WindTendencyStatus::null_source;
  }
  if (!valid_tendency_view(field.target) || !valid_tendency_view(field.source) ||
      !valid_tendency_view(field.advect_x) ||
      !valid_tendency_view(field.advect_y)) {
    return WindTendencyStatus::invalid_layout;
  }
  if (!has_one_cell_horizontal_halo(field.target) ||
      !has_one_cell_horizontal_halo(field.source)) {
    return WindTendencyStatus::insufficient_halo;
  }
  if (!same_layout(field.target, field.source)) {
    return WindTendencyStatus::mismatched_source_layout;
  }
  if (!same_layout(field.target, field.advect_x) ||
      !same_layout(field.target, field.advect_y)) {
    return WindTendencyStatus::mismatched_wind_layout;
  }
  return WindTendencyStatus::ok;
}

template <typename Real>
[[nodiscard]] WindTendencyStatus validate_contract(
    const WindTendencyViews<Real> views,
    const WindTendencyConfig<Real> config) noexcept {
  if (!positive(config.dx_m) || !positive(config.dy_m) ||
      config.dt_seconds < static_cast<Real>(0)) {
    return WindTendencyStatus::invalid_config;
  }

  const auto u_status = validate_field(views.u);
  if (u_status != WindTendencyStatus::ok) {
    return u_status;
  }
  const auto v_status = validate_field(views.v);
  if (v_status != WindTendencyStatus::ok) {
    return v_status;
  }
  return WindTendencyStatus::ok;
}

template <typename Real>
void apply_first_order_horizontal_advection(
    const WindTendencyFieldViews<Real> field,
    const WindTendencyConfig<Real> config) noexcept {
  const auto i_begin = field.target.halo.i_lower;
  const auto i_end = field.target.nx - field.target.halo.i_upper;
  const auto j_begin = field.target.halo.j_lower;
  const auto j_end = field.target.ny - field.target.halo.j_upper;
  const auto k_begin = field.target.halo.k_lower;
  const auto k_end = field.target.nz - field.target.halo.k_upper;
  const auto inv_dx = static_cast<Real>(1) / config.dx_m;
  const auto inv_dy = static_cast<Real>(1) / config.dy_m;

  for (std::int32_t j = j_begin; j < j_end; ++j) {
    for (std::int32_t k = k_begin; k < k_end; ++k) {
      const auto plane =
          static_cast<std::size_t>(j) *
              static_cast<std::size_t>(field.target.stride_j) +
          static_cast<std::size_t>(k) *
              static_cast<std::size_t>(field.target.stride_k);
      for (std::int32_t i = i_begin; i < i_end; ++i) {
        const auto idx = plane + static_cast<std::size_t>(i);
        const auto x_wind = field.advect_x.data[idx];
        const auto y_wind = field.advect_y.data[idx];
        const auto d_dx =
            x_wind >= static_cast<Real>(0)
                ? (field.source.data[idx] -
                   field.source.data[field.source.index(i - 1, j, k)]) *
                      inv_dx
                : (field.source.data[field.source.index(i + 1, j, k)] -
                   field.source.data[idx]) *
                      inv_dx;
        const auto d_dy =
            y_wind >= static_cast<Real>(0)
                ? (field.source.data[idx] -
                   field.source.data[field.source.index(i, j - 1, k)]) *
                      inv_dy
                : (field.source.data[field.source.index(i, j + 1, k)] -
                   field.source.data[idx]) *
                      inv_dy;
        field.target.data[idx] += config.dt_seconds * (-x_wind * d_dx - y_wind * d_dy);
      }
    }
  }
}

}  // namespace

template <typename Real>
WindTendencyReport apply_horizontal_wind_tendency(
    const WindTendencyViews<Real> views,
    const WindTendencyConfig<Real> config) noexcept {
  auto report = make_report(config, validate_contract(views, config));
  if (report.status != WindTendencyStatus::ok) {
    return report;
  }

  report.active_u_points = active_point_count(views.u.target);
  report.active_v_points = active_point_count(views.v.target);
  report.active_points = report.active_u_points + report.active_v_points;

  if (!config.enable_horizontal_advection || config.dt_seconds == static_cast<Real>(0)) {
    return report;
  }

  apply_first_order_horizontal_advection(views.u, config);
  apply_first_order_horizontal_advection(views.v, config);
  report.updated_u_points = report.active_u_points;
  report.updated_v_points = report.active_v_points;
  return report;
}

template WindTendencyReport apply_horizontal_wind_tendency<float>(
    WindTendencyViews<float>,
    WindTendencyConfig<float>) noexcept;
template WindTendencyReport apply_horizontal_wind_tendency<double>(
    WindTendencyViews<double>,
    WindTendencyConfig<double>) noexcept;

}  // namespace tywrf::dynamics
