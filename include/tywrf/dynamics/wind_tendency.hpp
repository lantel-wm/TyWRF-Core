#pragma once

#include "tywrf/field_view.hpp"

#include <cstdint>

namespace tywrf::dynamics {

enum class WindTendencyStatus : std::uint8_t {
  ok,
  null_target,
  null_source,
  invalid_config,
  invalid_layout,
  insufficient_halo,
  mismatched_source_layout,
  mismatched_wind_layout,
};

enum class WindTendencyAdvectionForm : std::uint8_t {
  centered,
  upwind,
};

template <typename Real>
struct WindTendencyFieldViews {
  FieldView3D<Real> target;
  FieldView3D<const Real> source;
  FieldView3D<const Real> advect_x;
  FieldView3D<const Real> advect_y;
};

template <typename Real>
struct WindTendencyViews {
  WindTendencyFieldViews<Real> u;
  WindTendencyFieldViews<Real> v;
};

template <typename Real>
struct WindTendencyConfig {
  Real dt_seconds = 0;
  Real dx_m = 1;
  Real dy_m = 1;
  WindTendencyAdvectionForm advection_form = WindTendencyAdvectionForm::centered;
  bool enable_horizontal_advection = true;
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool validation_gate_evidence = false;
};

struct WindTendencyReport {
  WindTendencyStatus status = WindTendencyStatus::ok;
  std::int64_t active_points = 0;
  std::int64_t active_u_points = 0;
  std::int64_t active_v_points = 0;
  std::int64_t updated_u_points = 0;
  std::int64_t updated_v_points = 0;
  bool halo_writes = false;
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool validation_gate_evidence = false;
  bool horizontal_advection_enabled = false;
};

template <typename Real>
[[nodiscard]] WindTendencyReport apply_horizontal_wind_tendency(
    WindTendencyViews<Real> views,
    WindTendencyConfig<Real> config) noexcept;

}  // namespace tywrf::dynamics
