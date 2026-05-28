#pragma once

#include "tywrf/field_view.hpp"

#include <cstdint>

namespace tywrf::dynamics {

enum class PressureGradientWindTendencyStatus : std::uint8_t {
  ok,
  null_target,
  null_pressure,
  invalid_pressure_value,
  invalid_config,
  invalid_layout,
  insufficient_halo,
  mismatched_pressure_layout,
  mismatched_target_layout,
  unsupported_config,
  not_implemented,
};

enum class PressureGradientWindTendencyForm : std::uint8_t {
  first_order_constant_specific_volume,
  wrf_exact,
};

enum class PressureGradientPressureUnits : std::uint8_t {
  unspecified,
  pascal,
};

enum class PressureGradientSpecificVolumeUnits : std::uint8_t {
  unspecified,
  cubic_meter_per_kilogram,
};

template <typename Real>
struct PressureGradientWindTendencyViews {
  FieldView3D<Real> u_target;
  FieldView3D<Real> v_target;
  FieldView3D<const Real> perturbation_pressure;
  FieldView3D<const Real> base_pressure;
};

template <typename Real>
struct PressureGradientWindTendencyConfig {
  Real dt_seconds = 0;
  Real dx_m = 1;
  Real dy_m = 1;
  Real constant_specific_volume_m3_per_kg = 0;
  PressureGradientWindTendencyForm form =
      PressureGradientWindTendencyForm::first_order_constant_specific_volume;
  PressureGradientPressureUnits pressure_units =
      PressureGradientPressureUnits::unspecified;
  PressureGradientSpecificVolumeUnits specific_volume_units =
      PressureGradientSpecificVolumeUnits::unspecified;
  bool enable_pressure_gradient = false;
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool validation_gate_evidence = false;
};

struct PressureGradientWindTendencyReport {
  PressureGradientWindTendencyStatus status =
      PressureGradientWindTendencyStatus::ok;
  std::int64_t active_points = 0;
  std::int64_t active_u_points = 0;
  std::int64_t active_v_points = 0;
  std::int64_t updated_u_points = 0;
  std::int64_t updated_v_points = 0;
  std::int64_t skipped_u_points = 0;
  std::int64_t skipped_v_points = 0;
  bool halo_writes = false;
  bool pressure_gradient_enabled = false;
  bool used_base_pressure = false;
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool validation_gate_evidence = false;
};

template <typename Real>
[[nodiscard]] PressureGradientWindTendencyReport
apply_horizontal_pressure_gradient_wind_tendency(
    PressureGradientWindTendencyViews<Real> views,
    PressureGradientWindTendencyConfig<Real> config) noexcept;

}  // namespace tywrf::dynamics
