#include "tywrf/dynamics/pressure_gradient_wind_tendency.hpp"

#include "tywrf/dynamics/tendency.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace tywrf::dynamics {
namespace {

template <typename Real>
[[nodiscard]] bool finite(const Real value) noexcept {
  return std::isfinite(value);
}

template <typename Real>
[[nodiscard]] bool finite_positive(const Real value) noexcept {
  return finite(value) && value > static_cast<Real>(0);
}

template <typename Real>
[[nodiscard]] std::int32_t active_nx(const FieldView3D<Real> view) noexcept {
  return view.nx - view.halo.i_lower - view.halo.i_upper;
}

template <typename Real>
[[nodiscard]] std::int32_t active_ny(const FieldView3D<Real> view) noexcept {
  return view.ny - view.halo.j_lower - view.halo.j_upper;
}

template <typename Real>
[[nodiscard]] std::int32_t active_nz(const FieldView3D<Real> view) noexcept {
  return view.nz - view.halo.k_lower - view.halo.k_upper;
}

template <typename Real>
[[nodiscard]] bool positive_active_shape(const FieldView3D<Real> view) noexcept {
  return active_nx(view) > 0 && active_ny(view) > 0 && active_nz(view) > 0;
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] bool same_layout(
    const FieldView3D<LhsReal> lhs,
    const FieldView3D<RhsReal> rhs) noexcept {
  return lhs.nx == rhs.nx && lhs.ny == rhs.ny && lhs.nz == rhs.nz &&
         lhs.stride_i == rhs.stride_i && lhs.stride_k == rhs.stride_k &&
         lhs.stride_j == rhs.stride_j &&
         lhs.halo.i_lower == rhs.halo.i_lower &&
         lhs.halo.i_upper == rhs.halo.i_upper &&
         lhs.halo.j_lower == rhs.halo.j_lower &&
         lhs.halo.j_upper == rhs.halo.j_upper &&
         lhs.halo.k_lower == rhs.halo.k_lower &&
         lhs.halo.k_upper == rhs.halo.k_upper;
}

template <typename Real>
[[nodiscard]] bool staggered_targets_match_mass_pressure(
    const PressureGradientWindTendencyViews<Real> views) noexcept {
  const auto mass_nx = active_nx(views.perturbation_pressure);
  const auto mass_ny = active_ny(views.perturbation_pressure);
  const auto mass_nz = active_nz(views.perturbation_pressure);

  return active_nx(views.u_target) == mass_nx + 1 &&
         active_ny(views.u_target) == mass_ny &&
         active_nz(views.u_target) == mass_nz &&
         active_nx(views.v_target) == mass_nx &&
         active_ny(views.v_target) == mass_ny + 1 &&
         active_nz(views.v_target) == mass_nz;
}

[[nodiscard]] bool valid_form(
    const PressureGradientWindTendencyForm form) noexcept {
  switch (form) {
    case PressureGradientWindTendencyForm::first_order_constant_specific_volume:
    case PressureGradientWindTendencyForm::wrf_exact:
      return true;
  }
  return false;
}

[[nodiscard]] bool valid_pressure_units(
    const PressureGradientPressureUnits units) noexcept {
  switch (units) {
    case PressureGradientPressureUnits::unspecified:
    case PressureGradientPressureUnits::pascal:
      return true;
  }
  return false;
}

[[nodiscard]] bool valid_specific_volume_units(
    const PressureGradientSpecificVolumeUnits units) noexcept {
  switch (units) {
    case PressureGradientSpecificVolumeUnits::unspecified:
    case PressureGradientSpecificVolumeUnits::cubic_meter_per_kilogram:
      return true;
  }
  return false;
}

template <typename Real>
[[nodiscard]] PressureGradientWindTendencyReport make_report(
    const PressureGradientWindTendencyConfig<Real> config,
    const PressureGradientWindTendencyStatus status) noexcept {
  PressureGradientWindTendencyReport report{};
  report.status = status;
  report.pressure_gradient_enabled = config.enable_pressure_gradient;
  report.diagnostic_only = config.diagnostic_only;
  report.gate_candidate = config.gate_candidate;
  report.validation_gate_evidence = false;
  return report;
}

template <typename Real>
[[nodiscard]] PressureGradientWindTendencyStatus validate_config(
    const PressureGradientWindTendencyConfig<Real> config) noexcept {
  if (!finite(config.dt_seconds) || config.dt_seconds < static_cast<Real>(0) ||
      !finite_positive(config.dx_m) || !finite_positive(config.dy_m)) {
    return PressureGradientWindTendencyStatus::invalid_config;
  }
  if (!valid_form(config.form) || !valid_pressure_units(config.pressure_units) ||
      !valid_specific_volume_units(config.specific_volume_units)) {
    return PressureGradientWindTendencyStatus::unsupported_config;
  }
  if (config.form == PressureGradientWindTendencyForm::wrf_exact) {
    return PressureGradientWindTendencyStatus::not_implemented;
  }
  if (!config.enable_pressure_gradient) {
    return PressureGradientWindTendencyStatus::ok;
  }
  if (config.pressure_units != PressureGradientPressureUnits::pascal ||
      config.specific_volume_units !=
          PressureGradientSpecificVolumeUnits::cubic_meter_per_kilogram) {
    return PressureGradientWindTendencyStatus::unsupported_config;
  }
  if (!finite_positive(config.constant_specific_volume_m3_per_kg)) {
    return PressureGradientWindTendencyStatus::invalid_config;
  }
  return PressureGradientWindTendencyStatus::ok;
}

template <typename Real>
[[nodiscard]] PressureGradientWindTendencyStatus validate_views(
    const PressureGradientWindTendencyViews<Real> views) noexcept {
  if (views.u_target.data == nullptr || views.v_target.data == nullptr) {
    return PressureGradientWindTendencyStatus::null_target;
  }
  if (views.perturbation_pressure.data == nullptr ||
      views.base_pressure.data == nullptr) {
    return PressureGradientWindTendencyStatus::null_pressure;
  }
  if (!valid_tendency_view(views.u_target) ||
      !valid_tendency_view(views.v_target) ||
      !valid_tendency_view(views.perturbation_pressure) ||
      !valid_tendency_view(views.base_pressure) ||
      !positive_active_shape(views.u_target) ||
      !positive_active_shape(views.v_target) ||
      !positive_active_shape(views.perturbation_pressure) ||
      !positive_active_shape(views.base_pressure)) {
    return PressureGradientWindTendencyStatus::invalid_layout;
  }
  if (!same_layout(views.perturbation_pressure, views.base_pressure)) {
    return PressureGradientWindTendencyStatus::mismatched_pressure_layout;
  }
  if (!staggered_targets_match_mass_pressure(views)) {
    return PressureGradientWindTendencyStatus::mismatched_target_layout;
  }
  return PressureGradientWindTendencyStatus::ok;
}

template <typename Real>
[[nodiscard]] PressureGradientWindTendencyStatus validate_contract(
    const PressureGradientWindTendencyViews<Real> views,
    const PressureGradientWindTendencyConfig<Real> config) noexcept {
  if (const auto config_status = validate_config(config);
      config_status != PressureGradientWindTendencyStatus::ok) {
    return config_status;
  }
  return validate_views(views);
}

template <typename Real>
[[nodiscard]] Real total_pressure_at(
    const FieldView3D<const Real> perturbation_pressure,
    const FieldView3D<const Real> base_pressure,
    const std::size_t idx) noexcept {
  return perturbation_pressure.data[idx] + base_pressure.data[idx];
}

template <typename Real>
[[nodiscard]] bool active_mass_pressure_is_valid(
    const PressureGradientWindTendencyViews<Real> views) noexcept {
  const auto mass_nx = active_nx(views.perturbation_pressure);
  const auto mass_ny = active_ny(views.perturbation_pressure);
  const auto mass_nz = active_nz(views.perturbation_pressure);

  for (std::int32_t active_j = 0; active_j < mass_ny; ++active_j) {
    const auto pj = views.perturbation_pressure.halo.j_lower + active_j;
    for (std::int32_t active_k = 0; active_k < mass_nz; ++active_k) {
      const auto pk = views.perturbation_pressure.halo.k_lower + active_k;
      const auto p_plane =
          static_cast<std::size_t>(pj) *
              static_cast<std::size_t>(views.perturbation_pressure.stride_j) +
          static_cast<std::size_t>(pk) *
              static_cast<std::size_t>(views.perturbation_pressure.stride_k);
      for (std::int32_t active_i = 0; active_i < mass_nx; ++active_i) {
        const auto pi = views.perturbation_pressure.halo.i_lower + active_i;
        const auto idx = p_plane + static_cast<std::size_t>(pi);
        const auto total_pressure = total_pressure_at(
            views.perturbation_pressure, views.base_pressure, idx);
        if (!finite(total_pressure) || total_pressure <= static_cast<Real>(0)) {
          return false;
        }
      }
    }
  }

  return true;
}

template <typename Real>
[[nodiscard]] std::int64_t interior_u_point_count(
    const PressureGradientWindTendencyViews<Real> views) noexcept {
  const auto mass_nx = active_nx(views.perturbation_pressure);
  const auto mass_ny = active_ny(views.perturbation_pressure);
  const auto mass_nz = active_nz(views.perturbation_pressure);
  return static_cast<std::int64_t>(mass_nx - 1) *
         static_cast<std::int64_t>(mass_ny) *
         static_cast<std::int64_t>(mass_nz);
}

template <typename Real>
[[nodiscard]] std::int64_t interior_v_point_count(
    const PressureGradientWindTendencyViews<Real> views) noexcept {
  const auto mass_nx = active_nx(views.perturbation_pressure);
  const auto mass_ny = active_ny(views.perturbation_pressure);
  const auto mass_nz = active_nz(views.perturbation_pressure);
  return static_cast<std::int64_t>(mass_nx) *
         static_cast<std::int64_t>(mass_ny - 1) *
         static_cast<std::int64_t>(mass_nz);
}

template <typename Real>
void apply_first_order_u_pressure_gradient(
    const PressureGradientWindTendencyViews<Real> views,
    const PressureGradientWindTendencyConfig<Real> config) noexcept {
  const auto mass_nx = active_nx(views.perturbation_pressure);
  const auto mass_ny = active_ny(views.perturbation_pressure);
  const auto mass_nz = active_nz(views.perturbation_pressure);
  const auto inv_dx = static_cast<Real>(1) / config.dx_m;
  const auto factor =
      -config.dt_seconds * config.constant_specific_volume_m3_per_kg * inv_dx;

  for (std::int32_t active_j = 0; active_j < mass_ny; ++active_j) {
    const auto uj = views.u_target.halo.j_lower + active_j;
    const auto pj = views.perturbation_pressure.halo.j_lower + active_j;
    for (std::int32_t active_k = 0; active_k < mass_nz; ++active_k) {
      const auto uk = views.u_target.halo.k_lower + active_k;
      const auto pk = views.perturbation_pressure.halo.k_lower + active_k;
      const auto u_plane =
          static_cast<std::size_t>(uj) *
              static_cast<std::size_t>(views.u_target.stride_j) +
          static_cast<std::size_t>(uk) *
              static_cast<std::size_t>(views.u_target.stride_k);
      const auto p_plane =
          static_cast<std::size_t>(pj) *
              static_cast<std::size_t>(views.perturbation_pressure.stride_j) +
          static_cast<std::size_t>(pk) *
              static_cast<std::size_t>(views.perturbation_pressure.stride_k);
      for (std::int32_t active_i = 1; active_i < mass_nx; ++active_i) {
        const auto ui = views.u_target.halo.i_lower + active_i;
        const auto right_i =
            views.perturbation_pressure.halo.i_lower + active_i;
        const auto left_i = right_i - 1;
        const auto target_idx = u_plane + static_cast<std::size_t>(ui);
        const auto right_idx = p_plane + static_cast<std::size_t>(right_i);
        const auto left_idx = p_plane + static_cast<std::size_t>(left_i);
        const auto gradient =
            total_pressure_at(
                views.perturbation_pressure, views.base_pressure, right_idx) -
            total_pressure_at(
                views.perturbation_pressure, views.base_pressure, left_idx);
        views.u_target.data[target_idx] += factor * gradient;
      }
    }
  }
}

template <typename Real>
void apply_first_order_v_pressure_gradient(
    const PressureGradientWindTendencyViews<Real> views,
    const PressureGradientWindTendencyConfig<Real> config) noexcept {
  const auto mass_nx = active_nx(views.perturbation_pressure);
  const auto mass_ny = active_ny(views.perturbation_pressure);
  const auto mass_nz = active_nz(views.perturbation_pressure);
  const auto inv_dy = static_cast<Real>(1) / config.dy_m;
  const auto factor =
      -config.dt_seconds * config.constant_specific_volume_m3_per_kg * inv_dy;

  for (std::int32_t active_j = 1; active_j < mass_ny; ++active_j) {
    const auto vj = views.v_target.halo.j_lower + active_j;
    const auto up_j = views.perturbation_pressure.halo.j_lower + active_j;
    const auto down_j = up_j - 1;
    for (std::int32_t active_k = 0; active_k < mass_nz; ++active_k) {
      const auto vk = views.v_target.halo.k_lower + active_k;
      const auto pk = views.perturbation_pressure.halo.k_lower + active_k;
      const auto v_plane =
          static_cast<std::size_t>(vj) *
              static_cast<std::size_t>(views.v_target.stride_j) +
          static_cast<std::size_t>(vk) *
              static_cast<std::size_t>(views.v_target.stride_k);
      const auto up_plane =
          static_cast<std::size_t>(up_j) *
              static_cast<std::size_t>(views.perturbation_pressure.stride_j) +
          static_cast<std::size_t>(pk) *
              static_cast<std::size_t>(views.perturbation_pressure.stride_k);
      const auto down_plane =
          static_cast<std::size_t>(down_j) *
              static_cast<std::size_t>(views.perturbation_pressure.stride_j) +
          static_cast<std::size_t>(pk) *
              static_cast<std::size_t>(views.perturbation_pressure.stride_k);
      for (std::int32_t active_i = 0; active_i < mass_nx; ++active_i) {
        const auto vi = views.v_target.halo.i_lower + active_i;
        const auto pi = views.perturbation_pressure.halo.i_lower + active_i;
        const auto target_idx = v_plane + static_cast<std::size_t>(vi);
        const auto up_idx = up_plane + static_cast<std::size_t>(pi);
        const auto down_idx = down_plane + static_cast<std::size_t>(pi);
        const auto gradient =
            total_pressure_at(
                views.perturbation_pressure, views.base_pressure, up_idx) -
            total_pressure_at(
                views.perturbation_pressure, views.base_pressure, down_idx);
        views.v_target.data[target_idx] += factor * gradient;
      }
    }
  }
}

template <typename Real>
void apply_first_order_pressure_gradient(
    const PressureGradientWindTendencyViews<Real> views,
    const PressureGradientWindTendencyConfig<Real> config) noexcept {
  apply_first_order_u_pressure_gradient(views, config);
  apply_first_order_v_pressure_gradient(views, config);
}

}  // namespace

template <typename Real>
PressureGradientWindTendencyReport
apply_horizontal_pressure_gradient_wind_tendency(
    const PressureGradientWindTendencyViews<Real> views,
    const PressureGradientWindTendencyConfig<Real> config) noexcept {
  auto report = make_report(config, validate_contract(views, config));
  if (report.status != PressureGradientWindTendencyStatus::ok) {
    return report;
  }

  report.active_u_points = active_point_count(views.u_target);
  report.active_v_points = active_point_count(views.v_target);
  report.active_points = report.active_u_points + report.active_v_points;
  report.used_base_pressure = true;

  if (!config.enable_pressure_gradient ||
      config.dt_seconds == static_cast<Real>(0)) {
    return report;
  }

  if (!active_mass_pressure_is_valid(views)) {
    report.status = PressureGradientWindTendencyStatus::invalid_pressure_value;
    return report;
  }

  apply_first_order_pressure_gradient(views, config);
  report.updated_u_points = interior_u_point_count(views);
  report.updated_v_points = interior_v_point_count(views);
  report.skipped_u_points = report.active_u_points - report.updated_u_points;
  report.skipped_v_points = report.active_v_points - report.updated_v_points;
  return report;
}

template PressureGradientWindTendencyReport
apply_horizontal_pressure_gradient_wind_tendency<float>(
    PressureGradientWindTendencyViews<float>,
    PressureGradientWindTendencyConfig<float>) noexcept;
template PressureGradientWindTendencyReport
apply_horizontal_pressure_gradient_wind_tendency<double>(
    PressureGradientWindTendencyViews<double>,
    PressureGradientWindTendencyConfig<double>) noexcept;

}  // namespace tywrf::dynamics
