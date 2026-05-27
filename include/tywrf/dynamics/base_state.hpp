#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/nest/nest_interface.hpp"

#include <cstdint>

namespace tywrf::dynamics {

struct BaseStateReconstructionInputs {
  FieldView3D<const float> pb;
  FieldView3D<const float> t_init;
  FieldView3D<float> alb;
};

struct BaseStateReconstructionOptions {
  float dry_air_gas_constant = 287.0F;
  float reference_pressure_pa = 100'000.0F;
  float base_potential_temperature_k = 300.0F;
  float specific_heat_cp = 1004.5F;
  // WRF module_model_constants defines cvpm = -cv/cp, where cv = cp - R_d.
  // This is not -R_d/Cp.
  float cvpm = -(1004.5F - 287.0F) / 1004.5F;
};

struct BaseStateReconstructionReport {
  nest::NestResult result{nest::NestStatus::ok, "ok"};
  std::int32_t active_nx = 0;
  std::int32_t active_ny = 0;
  std::int32_t active_nz = 0;
  std::uint64_t reconstructed_point_count = 0;
  std::uint64_t invalid_point_count = 0;
  bool requires_external_t_init_and_pb = true;
  bool used_wrf_inverse_base_density_alb = false;
  bool wrote_alb = false;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

struct BaseStateVerticalCoefficientView {
  const float* values = nullptr;
  std::int32_t count = 0;
};

struct KrosaMassBaseStateReconstructionInputs {
  FieldView2D<const float> terrain_height_m;
  BaseStateVerticalCoefficientView c3h;
  BaseStateVerticalCoefficientView c4h;
  BaseStateVerticalCoefficientView c3f;
  BaseStateVerticalCoefficientView c4f;
  float p_top_pa = 0.0F;
};

struct KrosaMassBaseStateReconstructionOutputs {
  FieldView3D<float> pb;
  FieldView3D<float> t_init;
  FieldView2D<float> mub;
  // Optional. Leave data null to reconstruct only PB/T_INIT/MUB.
  FieldView3D<float> alb;
  // Optional full-level PHB. Requires ALB output plus C3F/C4F full-level
  // coefficients; active_nz must be mass active_nz + 1.
  FieldView3D<float> phb;
};

struct KrosaMassBaseStateReconstructionOptions {
  float dry_air_gas_constant = 287.0F;
  float gravity = 9.81F;
  float sea_level_base_pressure_pa = 100'000.0F;
  float sea_level_base_temperature_k = 290.0F;
  float base_lapse_k = 50.0F;
  float isothermal_temperature_k = 200.0F;
  float stratosphere_base_pressure_pa = 0.0F;
  float stratosphere_lapse_k = -11.0F;
  float reference_pressure_pa = 100'000.0F;
  float base_potential_temperature_k = 300.0F;
  float specific_heat_cp = 1004.5F;
  // WRF module_model_constants defines cvpm = -cv/cp, where cv = cp - R_d.
  float cvpm = -(1004.5F - 287.0F) / 1004.5F;
};

struct KrosaMassBaseStateReconstructionReport {
  nest::NestResult result{nest::NestStatus::ok, "ok"};
  std::int32_t active_nx = 0;
  std::int32_t active_ny = 0;
  std::int32_t active_nz = 0;
  std::uint64_t reconstructed_column_count = 0;
  std::uint64_t reconstructed_point_count = 0;
  std::uint64_t invalid_column_count = 0;
  std::uint64_t invalid_point_count = 0;
  std::uint64_t unsupported_stratosphere_point_count = 0;
  bool used_wrf_surface_pressure_formula = false;
  bool used_wrf_hybrid_mass_coefficients = false;
  bool used_wrf_t_init_formula = false;
  bool reused_alb_helper = false;
  bool wrote_pb = false;
  bool wrote_t_init = false;
  bool wrote_mub = false;
  bool wrote_alb = false;
  bool wrote_phb = false;
  bool phb_full_level_reconstruction_implemented = false;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// Reconstructs WRF base-state inverse density ALB from externally supplied PB
// and T_INIT mass-level fields:
//
//   ALB = (Rd / p1000mb) * (T_INIT + T0) * (PB / p1000mb) ** cvpm
//
// This is intentionally narrow: it does not derive PB or T_INIT, does not read
// NetCDF, does not touch State, and must not be confused with surface ALBEDO.
[[nodiscard]] BaseStateReconstructionReport reconstruct_alb_from_pb_t_init(
    const BaseStateReconstructionInputs& inputs,
    BaseStateReconstructionOptions options = {}) noexcept;

// Reconstructs KROSA-scoped WRF mass-level base-state fields from terrain,
// WRF base constants, and mass-level hybrid coefficients:
//
//   p_surf = p00 * exp(-t00/a + sqrt((t00/a)^2 - 2*g*HT/(a*Rd)))
//   MUB    = p_surf - P_TOP
//   PB     = C3H[k] * MUB + C4H[k] + P_TOP
//   T_INIT = T_base(PB) * (p00/PB) ** (Rd/Cp) - T0
//
// If an ALB output view is supplied, ALB is derived with
// reconstruct_alb_from_pb_t_init after PB/T_INIT are available. If a full-level
// PHB output view is supplied, this uses the KROSA HYPSOMETRIC_OPT=2 log-linear
// branch only:
//
//   PHB[k=0] = HGT * g
//   pfu      = C3F[k+1] * MUB + C4F[k+1] + P_TOP
//   pfd      = C3F[k]   * MUB + C4F[k]   + P_TOP
//   phm      = C3H[k]   * MUB + C4H[k]   + P_TOP
//   PHB[k+1] = PHB[k] + ALB[k] * phm * log(pfd / pfu)
//
// PHB reconstruction requires an ALB output view; this helper does not read
// later restart ALB or any NetCDF source.
[[nodiscard]] KrosaMassBaseStateReconstructionReport
reconstruct_krosa_mass_base_state(
    const KrosaMassBaseStateReconstructionInputs& inputs,
    const KrosaMassBaseStateReconstructionOutputs& outputs,
    KrosaMassBaseStateReconstructionOptions options = {}) noexcept;

}  // namespace tywrf::dynamics
