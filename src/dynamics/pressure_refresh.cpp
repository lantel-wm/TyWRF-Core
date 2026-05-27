#include "tywrf/dynamics/pressure_refresh.hpp"

#include <cmath>
#include <cstdint>

namespace tywrf::dynamics {
namespace {

[[nodiscard]] constexpr nest::NestResult ok_result() noexcept {
  return {nest::NestStatus::ok, "ok"};
}

[[nodiscard]] constexpr nest::NestResult result(
    const nest::NestStatus status,
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
    const FieldView3D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

[[nodiscard]] constexpr bool window_fits_mass_field(
    const nest::RemapWindow& window,
    const FieldView3D<float>& p) noexcept {
  return window.stagger == nest::HorizontalStagger::mass &&
         window.old_i_begin >= 0 && window.old_j_begin >= 0 &&
         window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0 &&
         window.old_i_begin + window.extent_i <= active_nx(p) &&
         window.old_j_begin + window.extent_j <= active_ny(p) &&
         window.new_i_begin + window.extent_i <= active_nx(p) &&
         window.new_j_begin + window.extent_j <= active_ny(p);
}

[[nodiscard]] bool finite_positive(const float value) noexcept {
  return std::isfinite(value) && value > 0.0F;
}

[[nodiscard]] constexpr bool valid_coefficients(
    const VerticalCoefficientView coefficients,
    const std::int32_t required_count) noexcept {
  return coefficients.values != nullptr && coefficients.count >= required_count;
}

[[nodiscard]] nest::NestResult validate_options(
    const KrosaPressureRefreshOptions& options) noexcept {
  switch (options.formula) {
    case PressureRefreshFormula::krosa_hypsometric_opt2_use_theta_m:
      break;
    case PressureRefreshFormula::wrf_start_em_general:
      return result(
          nest::NestStatus::not_implemented,
          "general WRF start_em pressure refresh is not implemented");
    default:
      return result(
          nest::NestStatus::not_implemented,
          "unsupported pressure refresh formula");
  }

  switch (options.region) {
    case PressureRefreshRegion::exposed_mass_cells:
    case PressureRefreshRegion::full_active_columns:
      break;
    default:
      return result(
          nest::NestStatus::not_implemented,
          "unsupported pressure refresh region");
  }

  switch (options.thermodynamic_mode) {
    case PressureRefreshThermodynamicMode::use_theta_m_1:
      break;
    case PressureRefreshThermodynamicMode::use_theta_m_0_qvapor:
      return result(
          nest::NestStatus::not_implemented,
          "USE_THETA_M=0 pressure refresh requires a separate QVAPOR path");
    default:
      return result(
          nest::NestStatus::not_implemented,
          "unsupported pressure refresh thermodynamic mode");
  }

  if (!finite_positive(options.dry_air_gas_constant) ||
      !finite_positive(options.reference_pressure_pa) ||
      !finite_positive(options.base_potential_temperature_k) ||
      !finite_positive(options.cp_over_cv) ||
      !finite_positive(options.min_total_pressure_pa) ||
      !finite_positive(options.min_log_pressure_ratio)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh constants must be finite and positive");
  }

  return ok_result();
}

[[nodiscard]] nest::NestResult validate_contract(
    const nest::RemapWindow& mass_window,
    const PressureRefreshInputs& inputs,
    const KrosaPressureRefreshOptions& options) noexcept {
  if (const auto option_status = validate_options(options); !option_status.ok()) {
    return option_status;
  }

  if (!valid_canonical_view(inputs.p) || !valid_canonical_view(inputs.pb) ||
      !valid_canonical_view(inputs.t) || !valid_canonical_view(inputs.alb) ||
      !valid_canonical_view(inputs.ph) || !valid_canonical_view(inputs.phb) ||
      !valid_canonical_view(inputs.mu) || !valid_canonical_view(inputs.mub)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh inputs must be non-null canonical field views");
  }

  if (!same_active_shape(inputs.p, inputs.pb) ||
      !same_active_shape(inputs.p, inputs.t) ||
      !same_active_shape(inputs.p, inputs.alb)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh mass views must have the same active shape");
  }

  if (!same_horizontal_active_shape(inputs.p, inputs.ph) ||
      !same_horizontal_active_shape(inputs.p, inputs.phb) ||
      active_nz(inputs.ph) != active_nz(inputs.p) + 1 ||
      active_nz(inputs.phb) != active_nz(inputs.p) + 1) {
    return result(
        nest::NestStatus::invalid_contract,
        "PH/PHB active shape must be mass shape with one extra vertical level");
  }

  if (!same_horizontal_active_shape(inputs.mu, inputs.p) ||
      !same_horizontal_active_shape(inputs.mub, inputs.p)) {
    return result(
        nest::NestStatus::invalid_contract,
        "MU/MUB active shape must match mass horizontal shape");
  }

  const auto mass_nz = active_nz(inputs.p);
  const auto full_nz = mass_nz + 1;
  if (!valid_coefficients(inputs.c3f, full_nz) ||
      !valid_coefficients(inputs.c4f, full_nz) ||
      !valid_coefficients(inputs.c3h, mass_nz) ||
      !valid_coefficients(inputs.c4h, mass_nz)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh vertical coefficients are missing or too short");
  }

  if (!std::isfinite(inputs.p_top_pa) || inputs.p_top_pa < 0.0F) {
    return result(
        nest::NestStatus::invalid_contract,
        "P_TOP must be finite and non-negative");
  }

  if (!window_fits_mass_field(mass_window, inputs.p)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh mass window is outside active field extents");
  }

  return ok_result();
}

[[nodiscard]] constexpr bool in_new_overlap_window(
    const nest::RemapWindow& window,
    const std::int32_t i,
    const std::int32_t j) noexcept {
  return i >= window.new_i_begin && i < window.new_i_begin + window.extent_i &&
         j >= window.new_j_begin && j < window.new_j_begin + window.extent_j;
}

template <typename Real>
[[nodiscard]] float mass_value(
    const FieldView3D<Real>& field,
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k) noexcept {
  return field(
      field.halo.i_lower + i,
      field.halo.j_lower + j,
      field.halo.k_lower + k);
}

template <typename Real>
[[nodiscard]] float full_value(
    const FieldView3D<Real>& field,
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k) noexcept {
  return field(
      field.halo.i_lower + i,
      field.halo.j_lower + j,
      field.halo.k_lower + k);
}

template <typename Real>
[[nodiscard]] float surface_value(
    const FieldView2D<Real>& field,
    const std::int32_t i,
    const std::int32_t j) noexcept {
  return field(field.halo.i_lower + i, field.halo.j_lower + j);
}

[[nodiscard]] bool compute_krosa_pressure(
    const PressureRefreshInputs& inputs,
    const KrosaPressureRefreshOptions& options,
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k,
    float& perturbation_pressure_pa) noexcept {
  const auto mu_total =
      static_cast<double>(surface_value(inputs.mu, i, j)) +
      static_cast<double>(surface_value(inputs.mub, i, j));
  const auto pfu =
      static_cast<double>(inputs.c3f.values[k + 1]) * mu_total +
      static_cast<double>(inputs.c4f.values[k + 1]) +
      static_cast<double>(inputs.p_top_pa);
  const auto pfd =
      static_cast<double>(inputs.c3f.values[k]) * mu_total +
      static_cast<double>(inputs.c4f.values[k]) +
      static_cast<double>(inputs.p_top_pa);
  const auto phm =
      static_cast<double>(inputs.c3h.values[k]) * mu_total +
      static_cast<double>(inputs.c4h.values[k]) +
      static_cast<double>(inputs.p_top_pa);

  if (!std::isfinite(mu_total) || mu_total <= 0.0 || !std::isfinite(pfu) ||
      !std::isfinite(pfd) || !std::isfinite(phm) ||
      pfu <= static_cast<double>(options.min_total_pressure_pa) ||
      pfd <= static_cast<double>(options.min_total_pressure_pa) ||
      phm <= static_cast<double>(options.min_total_pressure_pa) || pfd <= pfu) {
    return false;
  }

  const auto log_ratio = std::log(pfd / pfu);
  const auto phi_lower =
      static_cast<double>(full_value(inputs.ph, i, j, k)) +
      static_cast<double>(full_value(inputs.phb, i, j, k));
  const auto phi_upper =
      static_cast<double>(full_value(inputs.ph, i, j, k + 1)) +
      static_cast<double>(full_value(inputs.phb, i, j, k + 1));
  const auto delta_phi = phi_upper - phi_lower;
  const auto alb = static_cast<double>(mass_value(inputs.alb, i, j, k));
  const auto pb = static_cast<double>(mass_value(inputs.pb, i, j, k));
  const auto theta =
      static_cast<double>(options.base_potential_temperature_k) +
      static_cast<double>(mass_value(inputs.t, i, j, k));

  if (!std::isfinite(log_ratio) ||
      log_ratio <= static_cast<double>(options.min_log_pressure_ratio) ||
      !std::isfinite(delta_phi) || delta_phi <= 0.0 ||
      !std::isfinite(alb) || !std::isfinite(pb) ||
      pb <= static_cast<double>(options.min_total_pressure_pa) ||
      !std::isfinite(theta) || theta <= 0.0) {
    return false;
  }

  const auto alpha_total = delta_phi / (phm * log_ratio);
  const auto alpha_perturbation = alpha_total - alb;
  const auto alpha_from_wrf_branch = alpha_perturbation + alb;
  if (!std::isfinite(alpha_from_wrf_branch) || alpha_from_wrf_branch <= 0.0) {
    return false;
  }

  const auto pressure_base =
      static_cast<double>(options.dry_air_gas_constant) * theta /
      (static_cast<double>(options.reference_pressure_pa) * alpha_from_wrf_branch);
  if (!std::isfinite(pressure_base) || pressure_base <= 0.0) {
    return false;
  }

  const auto total_pressure =
      static_cast<double>(options.reference_pressure_pa) *
      std::pow(pressure_base, static_cast<double>(options.cp_over_cv));
  const auto perturbation = total_pressure - pb;
  if (!std::isfinite(total_pressure) ||
      total_pressure <= static_cast<double>(options.min_total_pressure_pa) ||
      !std::isfinite(perturbation)) {
    return false;
  }

  perturbation_pressure_pa = static_cast<float>(perturbation);
  return std::isfinite(perturbation_pressure_pa);
}

}  // namespace

PressureRefreshReport refresh_krosa_moving_nest_pressure(
    const nest::RemapWindow& mass_window,
    const PressureRefreshInputs& inputs,
    const KrosaPressureRefreshOptions options) noexcept {
  PressureRefreshReport report{};
  report.result = ok_result();
  report.full_column_mode =
      options.region == PressureRefreshRegion::full_active_columns;

  if (const auto validation = validate_contract(mass_window, inputs, options);
      !validation.ok()) {
    report.result = validation;
    return report;
  }

  report.used_krosa_hypsometric_opt2 = true;
  report.used_use_theta_m = true;
  report.used_qvapor = false;
  report.used_wrf_vertical_coefficients = true;

  const auto nx = active_nx(inputs.p);
  const auto ny = active_ny(inputs.p);
  const auto nz = active_nz(inputs.p);
  const auto p_i0 = inputs.p.halo.i_lower;
  const auto p_j0 = inputs.p.halo.j_lower;
  const auto p_k0 = inputs.p.halo.k_lower;

  for (std::int32_t j = 0; j < ny; ++j) {
    for (std::int32_t i = 0; i < nx; ++i) {
      const bool in_overlap = in_new_overlap_window(mass_window, i, j);
      if (options.region == PressureRefreshRegion::exposed_mass_cells &&
          in_overlap) {
        continue;
      }

      ++report.target_column_count;
      bool refreshed_column = false;
      for (std::int32_t k = 0; k < nz; ++k) {
        float refreshed_p = 0.0F;
        if (!compute_krosa_pressure(inputs, options, i, j, k, refreshed_p)) {
          ++report.invalid_point_count;
          ++report.skipped_point_count;
          continue;
        }

        inputs.p(p_i0 + i, p_j0 + j, p_k0 + k) = refreshed_p;
        ++report.refreshed_point_count;
        refreshed_column = true;
      }

      if (refreshed_column) {
        ++report.refreshed_column_count;
      }
      report.touched_overlap_cells = report.touched_overlap_cells ||
                                     (report.full_column_mode && in_overlap);
    }
  }

  return report;
}

PressureRefreshReport refresh_krosa_moving_nest_pressure(
    const nest::RemapPlan& plan,
    const PressureRefreshInputs& inputs,
    const KrosaPressureRefreshOptions options) noexcept {
  if (!plan.ok()) {
    PressureRefreshReport report{};
    report.result = plan.result;
    report.full_column_mode =
        options.region == PressureRefreshRegion::full_active_columns;
    return report;
  }

  return refresh_krosa_moving_nest_pressure(plan.mass, inputs, options);
}

}  // namespace tywrf::dynamics
