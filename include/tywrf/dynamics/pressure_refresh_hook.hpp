#pragma once

#include "tywrf/dynamics/base_state_provider.hpp"
#include "tywrf/dynamics/pressure_refresh.hpp"
#include "tywrf/dynamics/pressure_refresh_staging.hpp"
#include "tywrf/io/pressure_refresh_io.hpp"
#include "tywrf/nest/nest_interface.hpp"
#include "tywrf/state.hpp"

#include <cstdint>

namespace tywrf::dynamics {

struct KrosaPressureRefreshHookOptions {
  KrosaMassBaseStateReconstructionOptions base_state{};
  KrosaPressureRefreshOptions pressure_refresh{};
  const KrosaBaseStateProviderTerrainOverride* terrain_override = nullptr;
};

struct KrosaPressureRefreshHookReport {
  nest::NestResult result{nest::NestStatus::ok, "ok"};
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool integrator_output = false;
  bool provider_ok = false;
  bool staging_ok = false;
  bool pressure_refresh_applied = false;
  bool calls_pressure_refresh_compute = false;
  std::uint64_t synced_pb_point_count = 0;
  std::uint64_t synced_mub_point_count = 0;
  std::uint64_t synced_phb_point_count = 0;
  bool touched_overlap_cells = false;
  bool touched_halo_cells = false;
  KrosaBaseStateProviderReport provider_report{};
  PressureRefreshStagingReport staging_report{};
  PressureRefreshReport compute_report{};

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

[[nodiscard]] KrosaPressureRefreshHookReport
apply_krosa_moving_nest_pressure_refresh_hook(
    const nest::RemapPlan& plan,
    State<float>& new_child,
    const io::KrosaPressureRefreshMetadata& metadata,
    KrosaPressureRefreshHookOptions options = {});

}  // namespace tywrf::dynamics
