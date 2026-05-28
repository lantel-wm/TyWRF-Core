#pragma once

#include "tywrf/field_view.hpp"

#include <cstdint>

namespace tywrf::dynamics {

enum class MassContinuityTendencyStatus : std::uint8_t {
  ok,
  null_target,
  null_flux,
  invalid_mu_value,
  invalid_flux_value,
  invalid_config,
  invalid_layout,
  incompatible_u_flux_shape,
  incompatible_v_flux_shape,
};

template <typename Real>
struct MassContinuityTendencyViews {
  FieldView2D<Real> mu_target;
  FieldView2D<const Real> u_face_mass_flux;
  FieldView2D<const Real> v_face_mass_flux;
};

template <typename Real>
struct MassContinuityTendencyConfig {
  Real dt_seconds = 0;
  Real dx_m = 1;
  Real dy_m = 1;
  bool enable_mass_continuity = false;
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool validation_gate_evidence = false;
};

struct MassContinuityTendencyReport {
  MassContinuityTendencyStatus status = MassContinuityTendencyStatus::ok;
  std::int64_t active_points = 0;
  std::int64_t active_mu_points = 0;
  std::int64_t active_u_flux_points = 0;
  std::int64_t active_v_flux_points = 0;
  std::int64_t updated_mu_points = 0;
  std::int64_t skipped_mu_points = 0;
  std::int64_t skipped_boundary_points = 0;
  std::int64_t skipped_disabled_points = 0;
  bool halo_writes = false;
  bool mass_continuity_enabled = false;
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool validation_gate_evidence = false;
};

template <typename Real>
[[nodiscard]] MassContinuityTendencyReport
apply_dry_air_mass_continuity_mu_tendency(
    MassContinuityTendencyViews<Real> views,
    MassContinuityTendencyConfig<Real> config) noexcept;

}  // namespace tywrf::dynamics
