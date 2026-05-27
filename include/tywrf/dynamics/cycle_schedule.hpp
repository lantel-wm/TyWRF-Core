#pragma once

#include "tywrf/dynamics/dynamics_loop.hpp"

#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace tywrf::dynamics {

enum class CycleScheduleCallKind : std::uint8_t {
  boundary_input_refresh,
  spectral_nudging_input_refresh,
  d01_boundary_update,
  d01_spectral_nudging,
  moving_nest_move_check,
  vortex_center_recompute,
  parent_child_interpolation,
  two_way_feedback,
  history_output,
  moving_nest_position_update,
};

enum class SegmentEndpointPolicy : std::uint8_t {
  bracket_start_and_end,
};

enum class MovingNestTimingPolicy : std::uint8_t {
  snap_to_next_parent_step,
};

struct CycleScheduleConfig {
  std::int32_t parent_grid_spacing_m = 10'000;
  std::int32_t child_grid_spacing_m = 2'000;
  std::int32_t parent_time_step_seconds = 40;
  std::int32_t child_time_step_seconds = 8;
  std::int32_t parent_time_step_ratio = 5;
  std::int32_t segment_seconds = 21'600;
  std::int32_t boundary_refresh_interval_seconds = 21'600;
  std::int32_t spectral_nudging_input_interval_seconds = 21'600;
  std::int32_t moving_nest_interval_seconds = 900;
  std::int32_t history_interval_seconds = 21'600;
  SegmentEndpointPolicy endpoint_policy = SegmentEndpointPolicy::bracket_start_and_end;
  MovingNestTimingPolicy moving_nest_policy =
      MovingNestTimingPolicy::snap_to_next_parent_step;
};

struct CycleScheduleCall {
  CycleScheduleCallKind kind = CycleScheduleCallKind::d01_boundary_update;
  DomainId domain = DomainId::d01;
  std::int64_t sequence_index = 0;
  std::int64_t parent_step_index = 0;
  std::int32_t child_substep_index = -1;
  std::int64_t start_seconds = 0;
  std::int64_t end_seconds = 0;
  std::int64_t nominal_seconds = 0;

  [[nodiscard]] constexpr bool is_instant() const noexcept {
    return start_seconds == end_seconds;
  }
};

struct CycleScheduleSummary {
  std::int64_t parent_steps = 0;
  std::int64_t child_substeps = 0;
  std::int64_t boundary_input_refreshes = 0;
  std::int64_t spectral_nudging_input_refreshes = 0;
  std::int64_t d01_boundary_updates = 0;
  std::int64_t d01_spectral_nudging_calls = 0;
  std::int64_t moving_nest_move_checks = 0;
  std::int64_t vortex_center_recomputes = 0;
  std::int64_t moving_nest_position_updates = 0;
  std::int64_t parent_child_interpolations = 0;
  std::int64_t two_way_feedbacks = 0;
  std::int64_t history_outputs = 0;
};

class CycleSchedule {
 public:
  [[nodiscard]] static CycleSchedule build(CycleScheduleConfig config);

  [[nodiscard]] const CycleScheduleConfig& config() const noexcept;
  [[nodiscard]] const CycleScheduleSummary& summary() const noexcept;
  [[nodiscard]] std::span<const CycleScheduleCall> calls() const noexcept;

 private:
  CycleSchedule(CycleScheduleConfig config, CycleScheduleSummary summary,
                std::vector<CycleScheduleCall> calls);

  CycleScheduleConfig config_;
  CycleScheduleSummary summary_;
  std::vector<CycleScheduleCall> calls_;
};

[[nodiscard]] CycleScheduleConfig make_krosa_6h_cycle_schedule_config() noexcept;
[[nodiscard]] CycleSchedule build_krosa_6h_cycle_schedule();

[[nodiscard]] std::string_view cycle_schedule_call_name(
    CycleScheduleCallKind kind) noexcept;
[[nodiscard]] std::string_view segment_endpoint_policy_name(
    SegmentEndpointPolicy policy) noexcept;
[[nodiscard]] std::string_view moving_nest_timing_policy_name(
    MovingNestTimingPolicy policy) noexcept;

}  // namespace tywrf::dynamics
