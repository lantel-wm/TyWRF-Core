#include "tywrf/dynamics/mass_continuity_tendency.hpp"

#include "tywrf/dynamics/tendency.hpp"

#include <algorithm>
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
[[nodiscard]] std::int32_t active_nx(const FieldView2D<Real> view) noexcept {
  return view.nx - view.halo.i_lower - view.halo.i_upper;
}

template <typename Real>
[[nodiscard]] std::int32_t active_ny(const FieldView2D<Real> view) noexcept {
  return view.ny - view.halo.j_lower - view.halo.j_upper;
}

template <typename Real>
[[nodiscard]] bool positive_active_shape(const FieldView2D<Real> view) noexcept {
  return active_nx(view) > 0 && active_ny(view) > 0;
}

template <typename Real>
[[nodiscard]] bool finite_active_values(const FieldView2D<Real> view) noexcept {
  const auto i_begin = view.halo.i_lower;
  const auto i_end = view.nx - view.halo.i_upper;
  const auto j_begin = view.halo.j_lower;
  const auto j_end = view.ny - view.halo.j_upper;
  for (std::int32_t j = j_begin; j < j_end; ++j) {
    const auto row =
        static_cast<std::size_t>(j) * static_cast<std::size_t>(view.stride_j);
    for (std::int32_t i = i_begin; i < i_end; ++i) {
      if (!finite(view.data[row + static_cast<std::size_t>(i)])) {
        return false;
      }
    }
  }
  return true;
}

template <typename Real>
[[nodiscard]] bool valid_u_face_flux_shape(
    const FieldView2D<Real> mu_target,
    const FieldView2D<const Real> u_face_mass_flux) noexcept {
  const auto mu_nx = active_nx(mu_target);
  const auto mu_ny = active_ny(mu_target);
  const auto u_nx = active_nx(u_face_mass_flux);
  const auto u_ny = active_ny(u_face_mass_flux);
  return u_ny == mu_ny && (u_nx == mu_nx || u_nx == mu_nx + 1);
}

template <typename Real>
[[nodiscard]] bool valid_v_face_flux_shape(
    const FieldView2D<Real> mu_target,
    const FieldView2D<const Real> v_face_mass_flux) noexcept {
  const auto mu_nx = active_nx(mu_target);
  const auto mu_ny = active_ny(mu_target);
  const auto v_nx = active_nx(v_face_mass_flux);
  const auto v_ny = active_ny(v_face_mass_flux);
  return v_nx == mu_nx && (v_ny == mu_ny || v_ny == mu_ny + 1);
}

template <typename Real>
[[nodiscard]] std::int32_t updatable_nx(
    const MassContinuityTendencyViews<Real> views) noexcept {
  return std::max<std::int32_t>(
      0,
      std::min(active_nx(views.mu_target),
               std::min(active_nx(views.u_face_mass_flux) - 1,
                        active_nx(views.v_face_mass_flux))));
}

template <typename Real>
[[nodiscard]] std::int32_t updatable_ny(
    const MassContinuityTendencyViews<Real> views) noexcept {
  return std::max<std::int32_t>(
      0,
      std::min(active_ny(views.mu_target),
               std::min(active_ny(views.u_face_mass_flux),
                        active_ny(views.v_face_mass_flux) - 1)));
}

template <typename Real>
[[nodiscard]] std::int64_t updatable_point_count(
    const MassContinuityTendencyViews<Real> views) noexcept {
  return static_cast<std::int64_t>(updatable_nx(views)) *
         static_cast<std::int64_t>(updatable_ny(views));
}

template <typename Real>
[[nodiscard]] MassContinuityTendencyReport make_report(
    const MassContinuityTendencyConfig<Real> config,
    const MassContinuityTendencyStatus status) noexcept {
  MassContinuityTendencyReport report{};
  report.status = status;
  report.mass_continuity_enabled = config.enable_mass_continuity;
  report.diagnostic_only = config.diagnostic_only;
  report.gate_candidate = config.gate_candidate;
  report.validation_gate_evidence = false;
  return report;
}

template <typename Real>
[[nodiscard]] MassContinuityTendencyStatus validate_config(
    const MassContinuityTendencyConfig<Real> config) noexcept {
  if (!finite(config.dt_seconds) || config.dt_seconds < static_cast<Real>(0) ||
      !finite_positive(config.dx_m) || !finite_positive(config.dy_m)) {
    return MassContinuityTendencyStatus::invalid_config;
  }
  if (config.diagnostic_only && config.gate_candidate) {
    return MassContinuityTendencyStatus::invalid_config;
  }
  return MassContinuityTendencyStatus::ok;
}

template <typename Real>
[[nodiscard]] MassContinuityTendencyStatus validate_views(
    const MassContinuityTendencyViews<Real> views) noexcept {
  if (views.mu_target.data == nullptr) {
    return MassContinuityTendencyStatus::null_target;
  }
  if (views.u_face_mass_flux.data == nullptr ||
      views.v_face_mass_flux.data == nullptr) {
    return MassContinuityTendencyStatus::null_flux;
  }
  if (!valid_tendency_view(views.mu_target) ||
      !valid_tendency_view(views.u_face_mass_flux) ||
      !valid_tendency_view(views.v_face_mass_flux) ||
      !positive_active_shape(views.mu_target) ||
      !positive_active_shape(views.u_face_mass_flux) ||
      !positive_active_shape(views.v_face_mass_flux)) {
    return MassContinuityTendencyStatus::invalid_layout;
  }
  if (!valid_u_face_flux_shape(views.mu_target, views.u_face_mass_flux)) {
    return MassContinuityTendencyStatus::incompatible_u_flux_shape;
  }
  if (!valid_v_face_flux_shape(views.mu_target, views.v_face_mass_flux)) {
    return MassContinuityTendencyStatus::incompatible_v_flux_shape;
  }
  return MassContinuityTendencyStatus::ok;
}

template <typename Real>
[[nodiscard]] MassContinuityTendencyStatus validate_contract(
    const MassContinuityTendencyViews<Real> views,
    const MassContinuityTendencyConfig<Real> config) noexcept {
  if (const auto config_status = validate_config(config);
      config_status != MassContinuityTendencyStatus::ok) {
    return config_status;
  }
  return validate_views(views);
}

template <typename Real>
[[nodiscard]] MassContinuityTendencyStatus validate_finite_state(
    const MassContinuityTendencyViews<Real> views) noexcept {
  if (!finite_active_values(views.mu_target)) {
    return MassContinuityTendencyStatus::invalid_mu_value;
  }
  if (!finite_active_values(views.u_face_mass_flux) ||
      !finite_active_values(views.v_face_mass_flux)) {
    return MassContinuityTendencyStatus::invalid_flux_value;
  }
  return MassContinuityTendencyStatus::ok;
}

template <typename Real>
void apply_first_order_horizontal_divergence(
    const MassContinuityTendencyViews<Real> views,
    const MassContinuityTendencyConfig<Real> config) noexcept {
  const auto nx_update = updatable_nx(views);
  const auto ny_update = updatable_ny(views);
  const auto inv_dx = static_cast<Real>(1) / config.dx_m;
  const auto inv_dy = static_cast<Real>(1) / config.dy_m;

  for (std::int32_t active_j = 0; active_j < ny_update; ++active_j) {
    const auto mu_j = views.mu_target.halo.j_lower + active_j;
    const auto u_j = views.u_face_mass_flux.halo.j_lower + active_j;
    const auto v_south_j = views.v_face_mass_flux.halo.j_lower + active_j;
    const auto v_north_j = v_south_j + 1;
    const auto mu_row =
        static_cast<std::size_t>(mu_j) *
        static_cast<std::size_t>(views.mu_target.stride_j);
    const auto u_row =
        static_cast<std::size_t>(u_j) *
        static_cast<std::size_t>(views.u_face_mass_flux.stride_j);
    const auto v_south_row =
        static_cast<std::size_t>(v_south_j) *
        static_cast<std::size_t>(views.v_face_mass_flux.stride_j);
    const auto v_north_row =
        static_cast<std::size_t>(v_north_j) *
        static_cast<std::size_t>(views.v_face_mass_flux.stride_j);

    for (std::int32_t active_i = 0; active_i < nx_update; ++active_i) {
      const auto mu_i = views.mu_target.halo.i_lower + active_i;
      const auto u_west_i = views.u_face_mass_flux.halo.i_lower + active_i;
      const auto u_east_i = u_west_i + 1;
      const auto v_i = views.v_face_mass_flux.halo.i_lower + active_i;
      const auto mu_idx = mu_row + static_cast<std::size_t>(mu_i);
      const auto u_west_idx = u_row + static_cast<std::size_t>(u_west_i);
      const auto u_east_idx = u_row + static_cast<std::size_t>(u_east_i);
      const auto v_south_idx = v_south_row + static_cast<std::size_t>(v_i);
      const auto v_north_idx = v_north_row + static_cast<std::size_t>(v_i);
      const auto divergence =
          (views.u_face_mass_flux.data[u_east_idx] -
           views.u_face_mass_flux.data[u_west_idx]) *
              inv_dx +
          (views.v_face_mass_flux.data[v_north_idx] -
           views.v_face_mass_flux.data[v_south_idx]) *
              inv_dy;
      views.mu_target.data[mu_idx] -= config.dt_seconds * divergence;
    }
  }
}

}  // namespace

template <typename Real>
MassContinuityTendencyReport apply_dry_air_mass_continuity_mu_tendency(
    const MassContinuityTendencyViews<Real> views,
    const MassContinuityTendencyConfig<Real> config) noexcept {
  auto report = make_report(config, validate_contract(views, config));
  if (report.status != MassContinuityTendencyStatus::ok) {
    return report;
  }

  report.active_mu_points = active_point_count(views.mu_target);
  report.active_u_flux_points = active_point_count(views.u_face_mass_flux);
  report.active_v_flux_points = active_point_count(views.v_face_mass_flux);
  report.active_points = report.active_mu_points;

  if (const auto finite_status = validate_finite_state(views);
      finite_status != MassContinuityTendencyStatus::ok) {
    report.status = finite_status;
    return report;
  }

  if (!config.enable_mass_continuity ||
      config.dt_seconds == static_cast<Real>(0)) {
    report.skipped_mu_points = report.active_mu_points;
    if (!config.enable_mass_continuity) {
      report.skipped_disabled_points = report.active_mu_points;
    }
    return report;
  }

  apply_first_order_horizontal_divergence(views, config);
  report.updated_mu_points = updatable_point_count(views);
  report.skipped_boundary_points =
      report.active_mu_points - report.updated_mu_points;
  report.skipped_mu_points = report.skipped_boundary_points;
  return report;
}

template MassContinuityTendencyReport
apply_dry_air_mass_continuity_mu_tendency<float>(
    MassContinuityTendencyViews<float>,
    MassContinuityTendencyConfig<float>) noexcept;
template MassContinuityTendencyReport
apply_dry_air_mass_continuity_mu_tendency<double>(
    MassContinuityTendencyViews<double>,
    MassContinuityTendencyConfig<double>) noexcept;

}  // namespace tywrf::dynamics
