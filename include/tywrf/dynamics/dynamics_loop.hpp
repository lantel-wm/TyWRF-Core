#pragma once

#include <cstdint>
#include <string_view>

namespace tywrf::dynamics {

enum class DomainId : std::uint8_t {
  d01 = 1,
  d02 = 2,
};

enum class LoopEventKind : std::uint8_t {
  boundary_input_refresh,
  spectral_nudging_input_refresh,
  boundary_update,
  spectral_nudging,
  moving_nest_position_update,
  zero_dynamics_tendency,
  physics,
  nest_interpolation,
  nest_feedback,
  history_output,
};

struct DomainDescriptor {
  DomainId id = DomainId::d01;
  std::int32_t grid_spacing_m = 0;
  std::int32_t time_step_seconds = 0;
  bool has_lateral_boundary = false;
  bool has_spectral_nudging = false;
  bool is_moving_nest = false;
};

struct TimeStepDescriptor {
  std::int32_t parent_time_step_seconds = 0;
  std::int32_t child_time_step_seconds = 0;
  std::int32_t parent_time_step_ratio = 0;
  std::int32_t segment_seconds = 0;
  std::int32_t history_interval_seconds = 0;
  std::int32_t boundary_refresh_interval_seconds = 0;
  std::int32_t spectral_nudging_input_interval_seconds = 0;
};

struct DynamicsLoopConfig {
  DomainDescriptor parent;
  DomainDescriptor child;
  TimeStepDescriptor timing;
};

struct LoopEvent {
  LoopEventKind kind = LoopEventKind::zero_dynamics_tendency;
  DomainId domain = DomainId::d01;
  std::int64_t schedule_sequence_index = -1;
  std::int64_t parent_step_index = 0;
  std::int32_t child_substep_index = -1;
  std::int64_t start_seconds = 0;
  std::int64_t end_seconds = 0;
  std::int64_t nominal_seconds = 0;
};

using LoopEventCallback = void (*)(void* user_data, const LoopEvent& event);

struct LoopEventSink {
  void* user_data = nullptr;
  LoopEventCallback callback = nullptr;
};

struct LoopSummary {
  std::int64_t parent_steps = 0;
  std::int64_t child_steps = 0;
  std::int64_t schedule_calls = 0;
  std::int64_t boundary_input_refreshes = 0;
  std::int64_t spectral_nudging_input_refreshes = 0;
  std::int64_t boundary_updates = 0;
  std::int64_t spectral_nudging_calls = 0;
  std::int64_t moving_nest_position_updates = 0;
  std::int64_t dynamics_tendency_calls = 0;
  std::int64_t physics_calls = 0;
  std::int64_t nest_interpolations = 0;
  std::int64_t nest_feedbacks = 0;
  std::int64_t history_outputs = 0;
};

class DynamicsLoopRunner {
 public:
  explicit DynamicsLoopRunner(DynamicsLoopConfig config);

  [[nodiscard]] const DynamicsLoopConfig& config() const noexcept;
  [[nodiscard]] LoopSummary run(const LoopEventSink& sink = {}) const;

 private:
  DynamicsLoopConfig config_;
};

[[nodiscard]] DynamicsLoopConfig make_krosa_phase4_loop_config() noexcept;

[[nodiscard]] std::string_view domain_name(DomainId domain) noexcept;
[[nodiscard]] std::string_view loop_event_name(LoopEventKind kind) noexcept;

}  // namespace tywrf::dynamics
