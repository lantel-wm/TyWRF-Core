#include "tywrf/dynamics/pressure_refresh.hpp"

#include <cmath>
#include <cstdint>
#include <limits>

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

[[nodiscard]] constexpr bool observation_requested(
    const PressureRefreshObservationView observation) noexcept {
  return observation.request_count > 0;
}

[[nodiscard]] constexpr std::int32_t observation_record_limit(
    const PressureRefreshObservationView observation) noexcept {
  if (!observation_requested(observation) || observation.records == nullptr ||
      observation.record_capacity <= 0) {
    return 0;
  }
  return observation.request_count < observation.record_capacity
             ? observation.request_count
             : observation.record_capacity;
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

  const auto observation = options.observation;
  if (observation.request_count < 0 || observation.record_capacity < 0) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh observation counts must be non-negative");
  }

  if (observation_requested(observation) &&
      (observation.requests == nullptr || observation.records == nullptr ||
       observation.record_capacity <= 0)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh observation requires request and record buffers");
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

[[nodiscard]] constexpr bool in_refresh_target_region(
    const nest::RemapWindow& window,
    const PressureRefreshRegion region,
    const std::int32_t i,
    const std::int32_t j) noexcept {
  return region == PressureRefreshRegion::full_active_columns ||
         !in_new_overlap_window(window, i, j);
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

[[nodiscard]] PressureRefreshFormulaObservation empty_observation(
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k) noexcept {
  const auto missing = std::numeric_limits<double>::quiet_NaN();
  PressureRefreshFormulaObservation observation{};
  observation.i = i;
  observation.j = j;
  observation.k = k;
  observation.status = PressureRefreshObservationStatus::not_recorded;
  observation.valid = 0;
  observation.mu_total = missing;
  observation.pfu = missing;
  observation.pfd = missing;
  observation.phm = missing;
  observation.log_ratio = missing;
  observation.phi_lower = missing;
  observation.phi_upper = missing;
  observation.delta_phi = missing;
  observation.alb = missing;
  observation.pb = missing;
  observation.theta = missing;
  observation.alpha_total = missing;
  observation.alpha_perturbation = missing;
  observation.alpha_from_wrf_branch = missing;
  observation.pressure_base = missing;
  observation.total_pressure = missing;
  observation.perturbation_pressure_pa = missing;
  return observation;
}

void finish_observation(
    PressureRefreshFormulaObservation* observation,
    const PressureRefreshObservationStatus status) noexcept {
  if (observation == nullptr) {
    return;
  }
  observation->status = status;
  observation->valid =
      status == PressureRefreshObservationStatus::recorded ? 1U : 0U;
}

[[nodiscard]] PressureRefreshObservationStatus compute_krosa_pressure(
    const PressureRefreshInputs& inputs,
    const KrosaPressureRefreshOptions& options,
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k,
    float& perturbation_pressure_pa,
    PressureRefreshFormulaObservation* observation = nullptr) noexcept {
  PressureRefreshFormulaObservation terms{};
  if (observation != nullptr) {
    terms = empty_observation(i, j, k);
  }

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

  if (observation != nullptr) {
    terms.mu_total = mu_total;
    terms.pfu = pfu;
    terms.pfd = pfd;
    terms.phm = phm;
  }

  if (!std::isfinite(mu_total) || mu_total <= 0.0) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_mu_total);
    return PressureRefreshObservationStatus::invalid_mu_total;
  }

  if (!std::isfinite(mu_total) || mu_total <= 0.0 || !std::isfinite(pfu) ||
      !std::isfinite(pfd) || !std::isfinite(phm) ||
      pfu <= static_cast<double>(options.min_total_pressure_pa) ||
      pfd <= static_cast<double>(options.min_total_pressure_pa) ||
      phm <= static_cast<double>(options.min_total_pressure_pa) || pfd <= pfu) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_pressure_levels);
    return PressureRefreshObservationStatus::invalid_pressure_levels;
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

  if (observation != nullptr) {
    terms.log_ratio = log_ratio;
    terms.phi_lower = phi_lower;
    terms.phi_upper = phi_upper;
    terms.delta_phi = delta_phi;
    terms.alb = alb;
    terms.pb = pb;
    terms.theta = theta;
  }

  if (!std::isfinite(log_ratio) ||
      log_ratio <= static_cast<double>(options.min_log_pressure_ratio)) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_log_ratio);
    return PressureRefreshObservationStatus::invalid_log_ratio;
  }

  if (!std::isfinite(delta_phi) || delta_phi <= 0.0) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_delta_phi);
    return PressureRefreshObservationStatus::invalid_delta_phi;
  }

  if (!std::isfinite(alb)) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_alb);
    return PressureRefreshObservationStatus::invalid_alb;
  }

  if (!std::isfinite(pb) ||
      pb <= static_cast<double>(options.min_total_pressure_pa)) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_pb);
    return PressureRefreshObservationStatus::invalid_pb;
  }

  if (!std::isfinite(theta) || theta <= 0.0) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_theta);
    return PressureRefreshObservationStatus::invalid_theta;
  }

  const auto alpha_total = delta_phi / (phm * log_ratio);
  const auto alpha_perturbation = alpha_total - alb;
  const auto alpha_from_wrf_branch = alpha_perturbation + alb;
  if (observation != nullptr) {
    terms.alpha_total = alpha_total;
    terms.alpha_perturbation = alpha_perturbation;
    terms.alpha_from_wrf_branch = alpha_from_wrf_branch;
  }
  if (!std::isfinite(alpha_from_wrf_branch) || alpha_from_wrf_branch <= 0.0) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_alpha);
    return PressureRefreshObservationStatus::invalid_alpha;
  }

  const auto pressure_base =
      static_cast<double>(options.dry_air_gas_constant) * theta /
      (static_cast<double>(options.reference_pressure_pa) * alpha_from_wrf_branch);
  if (observation != nullptr) {
    terms.pressure_base = pressure_base;
  }
  if (!std::isfinite(pressure_base) || pressure_base <= 0.0) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_pressure_base);
    return PressureRefreshObservationStatus::invalid_pressure_base;
  }

  const auto total_pressure =
      static_cast<double>(options.reference_pressure_pa) *
      std::pow(pressure_base, static_cast<double>(options.cp_over_cv));
  const auto perturbation = total_pressure - pb;
  if (observation != nullptr) {
    terms.total_pressure = total_pressure;
    terms.perturbation_pressure_pa = perturbation;
  }
  if (!std::isfinite(total_pressure) ||
      total_pressure <= static_cast<double>(options.min_total_pressure_pa) ||
      !std::isfinite(perturbation)) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_total_pressure);
    return PressureRefreshObservationStatus::invalid_total_pressure;
  }

  perturbation_pressure_pa = static_cast<float>(perturbation);
  if (!std::isfinite(perturbation_pressure_pa)) {
    if (observation != nullptr) {
      *observation = terms;
    }
    finish_observation(
        observation, PressureRefreshObservationStatus::invalid_total_pressure);
    return PressureRefreshObservationStatus::invalid_total_pressure;
  }

  if (observation != nullptr) {
    *observation = terms;
  }
  finish_observation(observation, PressureRefreshObservationStatus::recorded);
  return PressureRefreshObservationStatus::recorded;
}

[[nodiscard]] bool observation_is_pending_for_point(
    const PressureRefreshObservationView observation,
    const std::int32_t record_limit,
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k) noexcept {
  for (std::int32_t n = 0; n < record_limit; ++n) {
    const auto& record = observation.records[n];
    if (record.status == PressureRefreshObservationStatus::not_recorded &&
        record.i == i && record.j == j && record.k == k) {
      return true;
    }
  }
  return false;
}

void publish_observation_for_point(
    const PressureRefreshObservationView observation,
    const std::int32_t record_limit,
    const PressureRefreshFormulaObservation& computed,
    PressureRefreshReport& report) noexcept {
  for (std::int32_t n = 0; n < record_limit; ++n) {
    auto& record = observation.records[n];
    if (record.status == PressureRefreshObservationStatus::not_recorded &&
        record.i == computed.i && record.j == computed.j &&
        record.k == computed.k) {
      record = computed;
      if (record.valid != 0U) {
        ++report.observation_valid_count;
      } else {
        ++report.observation_invalid_count;
      }
    }
  }
}

void initialize_observations(
    const nest::RemapWindow& mass_window,
    const PressureRefreshInputs& inputs,
    const KrosaPressureRefreshOptions& options,
    PressureRefreshReport& report) noexcept {
  const auto observation = options.observation;
  if (!observation_requested(observation)) {
    return;
  }

  report.observation_request_count =
      static_cast<std::uint64_t>(observation.request_count);
  const auto record_limit = observation_record_limit(observation);
  report.observation_record_count = static_cast<std::uint64_t>(record_limit);
  report.observation_dropped_count = static_cast<std::uint64_t>(
      observation.request_count - record_limit);

  const auto nx = active_nx(inputs.p);
  const auto ny = active_ny(inputs.p);
  const auto nz = active_nz(inputs.p);

  for (std::int32_t n = 0; n < record_limit; ++n) {
    const auto request = observation.requests[n];
    auto record = empty_observation(request.i, request.j, request.k);
    if (request.i < 0 || request.i >= nx || request.j < 0 || request.j >= ny ||
        request.k < 0 || request.k >= nz) {
      record.status = PressureRefreshObservationStatus::request_out_of_bounds;
      ++report.observation_out_of_bounds_count;
    } else if (!in_refresh_target_region(
                   mass_window, options.region, request.i, request.j)) {
      record.status =
          PressureRefreshObservationStatus::request_outside_target_region;
      ++report.observation_outside_target_region_count;
    }
    observation.records[n] = record;
  }
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
  const auto observation = options.observation;
  const auto observation_limit = observation_record_limit(observation);
  initialize_observations(mass_window, inputs, options, report);

  for (std::int32_t j = 0; j < ny; ++j) {
    for (std::int32_t i = 0; i < nx; ++i) {
      const bool in_overlap = in_new_overlap_window(mass_window, i, j);
      if (!in_refresh_target_region(mass_window, options.region, i, j)) {
        continue;
      }

      ++report.target_column_count;
      bool refreshed_column = false;
      for (std::int32_t k = 0; k < nz; ++k) {
        float refreshed_p = 0.0F;
        PressureRefreshFormulaObservation observed{};
        const bool observe_point = observation_is_pending_for_point(
            observation, observation_limit, i, j, k);
        const auto status = compute_krosa_pressure(
            inputs, options, i, j, k, refreshed_p,
            observe_point ? &observed : nullptr);
        if (observe_point) {
          publish_observation_for_point(
              observation, observation_limit, observed, report);
        }
        if (status != PressureRefreshObservationStatus::recorded) {
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
