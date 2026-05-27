#include "tywrf/dynamics/cycle_schedule.hpp"

#include <stdexcept>
#include <utility>

namespace tywrf::dynamics {
namespace {

void require(const bool condition, const char* message) {
  if (!condition) {
    throw std::invalid_argument(message);
  }
}

[[nodiscard]] std::int64_t snap_to_next_parent_step_end(
    const std::int64_t nominal_seconds,
    const std::int64_t parent_time_step_seconds) noexcept {
  const auto snapped =
      ((nominal_seconds + parent_time_step_seconds - 1) /
       parent_time_step_seconds) *
      parent_time_step_seconds;
  if (snapped == 0) {
    return parent_time_step_seconds;
  }
  return snapped;
}

[[nodiscard]] std::int64_t instant_parent_step_index(
    const std::int64_t time_seconds,
    const std::int32_t parent_time_step_seconds) noexcept {
  if (time_seconds == 0) {
    return 0;
  }
  return (time_seconds / parent_time_step_seconds) - 1;
}

void validate_config(const CycleScheduleConfig& config) {
  require(config.parent_grid_spacing_m == 10'000,
          "KROSA schedule expects d01 horizontal spacing to be 10 km");
  require(config.child_grid_spacing_m == 2'000,
          "KROSA schedule requires d02 horizontal spacing to remain 2 km");
  require(config.parent_time_step_seconds > 0, "parent time step must be positive");
  require(config.child_time_step_seconds > 0, "child time step must be positive");
  require(config.parent_time_step_ratio > 0, "parent_time_step_ratio must be positive");
  require(config.parent_time_step_seconds ==
              config.child_time_step_seconds * config.parent_time_step_ratio,
          "parent time step must equal child time step multiplied by ratio");
  require(config.segment_seconds > 0, "segment length must be positive");
  require(config.segment_seconds % config.parent_time_step_seconds == 0,
          "segment length must be divisible by parent time step");
  require(config.boundary_refresh_interval_seconds == config.segment_seconds,
          "single-cycle schedule brackets boundary input at the 6 h segment endpoints");
  require(config.spectral_nudging_input_interval_seconds == config.segment_seconds,
          "single-cycle schedule brackets nudging input at the 6 h segment endpoints");
  require(config.history_interval_seconds == config.segment_seconds,
          "KROSA 6 h cycle schedule expects one history output at segment end");
  require(config.moving_nest_interval_seconds > 0,
          "moving nest interval must be positive");
  require(config.segment_seconds % config.moving_nest_interval_seconds == 0,
          "segment length must be divisible by moving nest interval");
}

void count_call(const CycleScheduleCallKind kind, CycleScheduleSummary& summary) noexcept {
  switch (kind) {
    case CycleScheduleCallKind::boundary_input_refresh:
      ++summary.boundary_input_refreshes;
      break;
    case CycleScheduleCallKind::spectral_nudging_input_refresh:
      ++summary.spectral_nudging_input_refreshes;
      break;
    case CycleScheduleCallKind::d01_boundary_update:
      ++summary.d01_boundary_updates;
      break;
    case CycleScheduleCallKind::d01_spectral_nudging:
      ++summary.d01_spectral_nudging_calls;
      break;
    case CycleScheduleCallKind::moving_nest_move_check:
      ++summary.moving_nest_move_checks;
      ++summary.moving_nest_position_updates;
      break;
    case CycleScheduleCallKind::vortex_center_recompute:
      ++summary.vortex_center_recomputes;
      break;
    case CycleScheduleCallKind::parent_child_interpolation:
      ++summary.parent_child_interpolations;
      break;
    case CycleScheduleCallKind::two_way_feedback:
      ++summary.two_way_feedbacks;
      break;
    case CycleScheduleCallKind::history_output:
      ++summary.history_outputs;
      break;
    case CycleScheduleCallKind::moving_nest_position_update:
      ++summary.moving_nest_move_checks;
      ++summary.moving_nest_position_updates;
      break;
  }
}

void emit_call(
    std::vector<CycleScheduleCall>& calls,
    CycleScheduleSummary& summary,
    const CycleScheduleCallKind kind,
    const DomainId domain,
    const std::int64_t parent_step_index,
    const std::int32_t child_substep_index,
    const std::int64_t start_seconds,
    const std::int64_t end_seconds,
    const std::int64_t nominal_seconds) {
  count_call(kind, summary);
  calls.push_back(
      CycleScheduleCall{kind, domain, static_cast<std::int64_t>(calls.size()),
                        parent_step_index, child_substep_index, start_seconds,
                        end_seconds, nominal_seconds});
}

}  // namespace

CycleSchedule::CycleSchedule(
    CycleScheduleConfig config,
    CycleScheduleSummary summary,
    std::vector<CycleScheduleCall> calls)
    : config_(config), summary_(summary), calls_(std::move(calls)) {}

CycleSchedule CycleSchedule::build(CycleScheduleConfig config) {
  validate_config(config);

  CycleScheduleSummary summary;
  std::vector<CycleScheduleCall> calls;

  const auto parent_steps =
      static_cast<std::int64_t>(config.segment_seconds / config.parent_time_step_seconds);
  const auto child_substeps =
      parent_steps * static_cast<std::int64_t>(config.parent_time_step_ratio);
  const auto vortex_center_recomputes =
      static_cast<std::int64_t>(config.segment_seconds / config.moving_nest_interval_seconds) + 1;
  const auto expected_calls =
      2 + parent_steps * 4 + child_substeps + vortex_center_recomputes + 2 + 2;
  calls.reserve(static_cast<std::size_t>(expected_calls));
  summary.parent_steps = parent_steps;
  summary.child_substeps = child_substeps;

  emit_call(calls, summary, CycleScheduleCallKind::boundary_input_refresh, DomainId::d01,
            0, -1, 0, 0, 0);
  emit_call(
      calls, summary, CycleScheduleCallKind::spectral_nudging_input_refresh,
      DomainId::d01, 0, -1, 0, 0, 0);
  for (std::int64_t parent_step = 0; parent_step < parent_steps; ++parent_step) {
    const auto parent_start =
        parent_step * static_cast<std::int64_t>(config.parent_time_step_seconds);
    const auto parent_end = parent_start + config.parent_time_step_seconds;

    emit_call(
        calls, summary, CycleScheduleCallKind::d01_boundary_update, DomainId::d01,
        parent_step, -1, parent_start, parent_end, parent_start);
    emit_call(
        calls, summary, CycleScheduleCallKind::d01_spectral_nudging, DomainId::d01,
        parent_step, -1, parent_start, parent_end, parent_start);

    for (std::int32_t child_substep = 0; child_substep < config.parent_time_step_ratio;
         ++child_substep) {
      const auto child_start =
          parent_start + child_substep *
                             static_cast<std::int64_t>(config.child_time_step_seconds);
      const auto child_end = child_start + config.child_time_step_seconds;
      emit_call(
          calls, summary, CycleScheduleCallKind::parent_child_interpolation,
          DomainId::d02, parent_step, child_substep, child_start, child_end,
          child_start);
    }

    emit_call(
        calls, summary, CycleScheduleCallKind::two_way_feedback, DomainId::d01,
        parent_step, -1, parent_start, parent_end, parent_end);
    emit_call(
        calls, summary, CycleScheduleCallKind::moving_nest_move_check,
        DomainId::d02, parent_step, -1, parent_end, parent_end, parent_end);

    for (std::int64_t nominal = 0;
         nominal <= config.segment_seconds;
         nominal += config.moving_nest_interval_seconds) {
      const auto scheduled =
          snap_to_next_parent_step_end(nominal, config.parent_time_step_seconds);
      if (scheduled == parent_end) {
        emit_call(
            calls, summary, CycleScheduleCallKind::vortex_center_recompute,
            DomainId::d02, parent_step, -1, scheduled, scheduled, nominal);
      }
    }

    if (parent_end % config.history_interval_seconds == 0) {
      emit_call(
          calls, summary, CycleScheduleCallKind::history_output, DomainId::d01,
          parent_step, -1, parent_end, parent_end, parent_end);
      emit_call(
          calls, summary, CycleScheduleCallKind::history_output, DomainId::d02,
          parent_step, -1, parent_end, parent_end, parent_end);
    }
  }

  emit_call(
      calls, summary, CycleScheduleCallKind::boundary_input_refresh, DomainId::d01,
      instant_parent_step_index(config.segment_seconds, config.parent_time_step_seconds), -1,
      config.segment_seconds, config.segment_seconds, config.segment_seconds);
  emit_call(
      calls, summary, CycleScheduleCallKind::spectral_nudging_input_refresh,
      DomainId::d01,
      instant_parent_step_index(config.segment_seconds, config.parent_time_step_seconds), -1,
      config.segment_seconds, config.segment_seconds, config.segment_seconds);

  return CycleSchedule(std::move(config), summary, std::move(calls));
}

const CycleScheduleConfig& CycleSchedule::config() const noexcept {
  return config_;
}

const CycleScheduleSummary& CycleSchedule::summary() const noexcept {
  return summary_;
}

std::span<const CycleScheduleCall> CycleSchedule::calls() const noexcept {
  return calls_;
}

CycleScheduleConfig make_krosa_6h_cycle_schedule_config() noexcept {
  return {};
}

CycleSchedule build_krosa_6h_cycle_schedule() {
  return CycleSchedule::build(make_krosa_6h_cycle_schedule_config());
}

std::string_view cycle_schedule_call_name(const CycleScheduleCallKind kind) noexcept {
  switch (kind) {
    case CycleScheduleCallKind::boundary_input_refresh:
      return "boundary_input_refresh";
    case CycleScheduleCallKind::spectral_nudging_input_refresh:
      return "spectral_nudging_input_refresh";
    case CycleScheduleCallKind::d01_boundary_update:
      return "d01_boundary_update";
    case CycleScheduleCallKind::d01_spectral_nudging:
      return "d01_spectral_nudging";
    case CycleScheduleCallKind::moving_nest_move_check:
      return "moving_nest_move_check";
    case CycleScheduleCallKind::vortex_center_recompute:
      return "vortex_center_recompute";
    case CycleScheduleCallKind::moving_nest_position_update:
      return "moving_nest_position_update";
    case CycleScheduleCallKind::parent_child_interpolation:
      return "parent_child_interpolation";
    case CycleScheduleCallKind::two_way_feedback:
      return "two_way_feedback";
    case CycleScheduleCallKind::history_output:
      return "history_output";
  }
  return "unknown";
}

std::string_view segment_endpoint_policy_name(
    const SegmentEndpointPolicy policy) noexcept {
  switch (policy) {
    case SegmentEndpointPolicy::bracket_start_and_end:
      return "bracket_start_and_end";
  }
  return "unknown";
}

std::string_view moving_nest_timing_policy_name(
    const MovingNestTimingPolicy policy) noexcept {
  switch (policy) {
    case MovingNestTimingPolicy::snap_to_next_parent_step:
      return "snap_to_next_parent_step";
  }
  return "unknown";
}

}  // namespace tywrf::dynamics
