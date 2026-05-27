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

}  // namespace tywrf::dynamics
