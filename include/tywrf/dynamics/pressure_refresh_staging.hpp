#pragma once

#include "tywrf/dynamics/pressure_refresh.hpp"
#include "tywrf/io/pressure_refresh_io.hpp"
#include "tywrf/nest/nest_interface.hpp"
#include "tywrf/state.hpp"

#include <cstdint>

namespace tywrf::dynamics {

// Stages the narrow KROSA pressure-refresh kernel inputs from an existing
// State plus externally supplied WRF ALB metadata. In this WRF branch ALB is
// inverse base density, not surface albedo. The d02 ALB buffer must come from a
// clean source, such as a future base-state reconstruction provider; this helper
// never reconstructs ALB and never treats later restart/history ALB as d02
// start-time truth.
struct PressureRefreshStagingReport {
  nest::NestResult result{nest::NestStatus::ok, "ok"};
  std::int32_t expected_mass_level_count = 0;
  std::int32_t expected_full_level_count = 0;
  std::int32_t c3f_count = 0;
  std::int32_t c4f_count = 0;
  std::int32_t c3h_count = 0;
  std::int32_t c4h_count = 0;
  bool alb_is_external = false;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

struct PressureRefreshStagingResult {
  PressureRefreshInputs inputs{};
  PressureRefreshStagingReport report{};

  [[nodiscard]] constexpr bool ok() const noexcept {
    return report.ok();
  }
};

[[nodiscard]] PressureRefreshStagingReport validate_krosa_pressure_refresh_staging(
    const StateView<float>& state,
    FieldView3D<const float> external_alb,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept;

[[nodiscard]] PressureRefreshStagingReport validate_krosa_pressure_refresh_staging(
    State<float>& state,
    const FieldStorage3D<float>& external_alb,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept;

[[nodiscard]] PressureRefreshStagingResult make_krosa_pressure_refresh_inputs(
    StateView<float> state,
    FieldView3D<const float> external_alb,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept;

[[nodiscard]] PressureRefreshStagingResult make_krosa_pressure_refresh_inputs(
    State<float>& state,
    const FieldStorage3D<float>& external_alb,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept;

}  // namespace tywrf::dynamics
