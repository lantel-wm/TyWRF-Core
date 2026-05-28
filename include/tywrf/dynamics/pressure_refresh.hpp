#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/nest/nest_interface.hpp"

#include <cstdint>

namespace tywrf::dynamics {

enum class PressureRefreshFormula : std::uint8_t {
  // Narrow KROSA branch matching WRF start_em HYPSOMETRIC_OPT=2 with
  // USE_THETA_M=1. QVAPOR is intentionally not an input for this branch.
  krosa_hypsometric_opt2_use_theta_m,

  // Explicit placeholders for unsupported future branches.
  wrf_start_em_general,
};

enum class PressureRefreshRegion : std::uint8_t {
  exposed_mass_cells,
  full_active_columns,
};

enum class PressureRefreshThermodynamicMode : std::uint8_t {
  use_theta_m_1,
  use_theta_m_0_qvapor,
};

struct VerticalCoefficientView {
  const float* values = nullptr;
  std::int32_t count = 0;
};

enum class PressureRefreshObservationStatus : std::uint8_t {
  not_recorded,
  recorded,
  request_out_of_bounds,
  request_outside_target_region,
  invalid_mu_total,
  invalid_pressure_levels,
  invalid_log_ratio,
  invalid_delta_phi,
  invalid_alb,
  invalid_pb,
  invalid_theta,
  invalid_alpha,
  invalid_pressure_base,
  invalid_total_pressure,
};

struct PressureRefreshObservationRequest {
  std::int32_t i = 0;
  std::int32_t j = 0;
  std::int32_t k = 0;
};

struct PressureRefreshFormulaObservation {
  std::int32_t i = 0;
  std::int32_t j = 0;
  std::int32_t k = 0;
  PressureRefreshObservationStatus status =
      PressureRefreshObservationStatus::not_recorded;
  std::uint8_t valid = 0;
  std::uint8_t reserved0 = 0;
  std::uint8_t reserved1 = 0;
  std::uint8_t reserved2 = 0;
  double mu_total = 0.0;
  double pfu = 0.0;
  double pfd = 0.0;
  double phm = 0.0;
  double log_ratio = 0.0;
  double phi_lower = 0.0;
  double phi_upper = 0.0;
  double delta_phi = 0.0;
  double alb = 0.0;
  double pb = 0.0;
  double theta = 0.0;
  double alpha_total = 0.0;
  double alpha_perturbation = 0.0;
  double alpha_from_wrf_branch = 0.0;
  double pressure_base = 0.0;
  double total_pressure = 0.0;
  double perturbation_pressure_pa = 0.0;
};

struct PressureRefreshObservationView {
  const PressureRefreshObservationRequest* requests = nullptr;
  std::int32_t request_count = 0;
  PressureRefreshFormulaObservation* records = nullptr;
  std::int32_t record_capacity = 0;
};

struct KrosaPressureRefreshOptions {
  PressureRefreshFormula formula =
      PressureRefreshFormula::krosa_hypsometric_opt2_use_theta_m;
  PressureRefreshRegion region = PressureRefreshRegion::exposed_mass_cells;
  PressureRefreshThermodynamicMode thermodynamic_mode =
      PressureRefreshThermodynamicMode::use_theta_m_1;
  float dry_air_gas_constant = 287.0F;
  float reference_pressure_pa = 100'000.0F;
  float base_potential_temperature_k = 300.0F;
  float cp_over_cv = 1004.5F / (1004.5F - 287.0F);
  float min_total_pressure_pa = 100.0F;
  float min_log_pressure_ratio = 1.0e-6F;
  PressureRefreshObservationView observation{};
};

struct PressureRefreshInputs {
  FieldView3D<float> p;
  FieldView3D<const float> pb;
  FieldView3D<const float> t;
  FieldView3D<const float> alb;
  FieldView3D<const float> ph;
  FieldView3D<const float> phb;
  FieldView2D<const float> mu;
  FieldView2D<const float> mub;
  VerticalCoefficientView c3f;
  VerticalCoefficientView c4f;
  VerticalCoefficientView c3h;
  VerticalCoefficientView c4h;
  float p_top_pa = 0.0F;
};

struct PressureRefreshReport {
  nest::NestResult result{nest::NestStatus::ok, "ok"};
  std::uint64_t target_column_count = 0;
  std::uint64_t refreshed_column_count = 0;
  std::uint64_t refreshed_point_count = 0;
  std::uint64_t skipped_point_count = 0;
  std::uint64_t invalid_point_count = 0;
  bool full_column_mode = false;
  bool touched_overlap_cells = false;
  bool touched_halo_cells = false;
  bool modified_pb = false;
  bool copied_parent_pressure = false;
  bool used_krosa_hypsometric_opt2 = false;
  bool used_use_theta_m = false;
  bool used_qvapor = false;
  bool used_wrf_vertical_coefficients = false;
  std::uint64_t observation_request_count = 0;
  std::uint64_t observation_record_count = 0;
  std::uint64_t observation_dropped_count = 0;
  std::uint64_t observation_valid_count = 0;
  std::uint64_t observation_invalid_count = 0;
  std::uint64_t observation_out_of_bounds_count = 0;
  std::uint64_t observation_outside_target_region_count = 0;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// Refreshes perturbation pressure P for KROSA moving-nest mass cells using the
// WRF start_em HYPSOMETRIC_OPT=2 / USE_THETA_M=1 branch:
//
//   pfu = C3F[k+1] * (MUB + MU) + C4F[k+1] + P_TOP
//   pfd = C3F[k]   * (MUB + MU) + C4F[k]   + P_TOP
//   phm = C3H[k]   * (MUB + MU) + C4H[k]   + P_TOP
//   alpha = (d(PH + PHB) / (phm * log(pfd / pfu)) - ALB) + ALB
//   P = p0 * (Rd * (t0 + T) / (p0 * alpha)) ** (cp / cv) - PB
//
// Overlap cells and halos are left unchanged in the default region mode. This
// is not a general WRF start_em replacement.
[[nodiscard]] PressureRefreshReport refresh_krosa_moving_nest_pressure(
    const nest::RemapWindow& mass_window,
    const PressureRefreshInputs& inputs,
    KrosaPressureRefreshOptions options = {}) noexcept;

[[nodiscard]] PressureRefreshReport refresh_krosa_moving_nest_pressure(
    const nest::RemapPlan& plan,
    const PressureRefreshInputs& inputs,
    KrosaPressureRefreshOptions options = {}) noexcept;

}  // namespace tywrf::dynamics
