#include "tywrf/dynamics/dynamics_loop.hpp"

#include "tywrf/dynamics/cycle_schedule.hpp"

#include <stdexcept>

namespace tywrf::dynamics {
namespace {

void require(const bool condition, const char* message) {
  if (!condition) {
    throw std::invalid_argument(message);
  }
}

void validate_config(const DynamicsLoopConfig& config) {
  const auto& parent = config.parent;
  const auto& child = config.child;
  const auto& timing = config.timing;

  require(parent.id == DomainId::d01, "Phase 4 dynamics skeleton expects d01 as parent");
  require(child.id == DomainId::d02, "Phase 4 dynamics skeleton expects d02 as child");
  require(parent.grid_spacing_m == 10'000,
          "KROSA dynamics schedule expects d01 horizontal spacing to be 10 km");
  require(child.grid_spacing_m == 2'000,
          "KROSA dynamics schedule requires d02 horizontal spacing to remain 2 km");
  require(parent.grid_spacing_m > 0, "parent grid spacing must be positive");
  require(child.grid_spacing_m > 0, "child grid spacing must be positive");
  require(parent.time_step_seconds > 0, "parent time step must be positive");
  require(child.time_step_seconds > 0, "child time step must be positive");
  require(parent.has_lateral_boundary,
          "KROSA dynamics schedule expects d01 lateral boundary updates");
  require(parent.has_spectral_nudging,
          "KROSA dynamics schedule expects d01 spectral nudging");
  require(child.is_moving_nest,
          "KROSA dynamics schedule expects d02 to remain a moving nest");
  require(timing.parent_time_step_seconds == parent.time_step_seconds,
          "parent descriptor and timing descriptor disagree");
  require(timing.child_time_step_seconds == child.time_step_seconds,
          "child descriptor and timing descriptor disagree");
  require(timing.parent_time_step_ratio > 0, "parent_time_step_ratio must be positive");
  require(parent.time_step_seconds ==
              child.time_step_seconds * timing.parent_time_step_ratio,
          "parent time step must equal child time step multiplied by parent_time_step_ratio");
  require(timing.segment_seconds > 0, "segment length must be positive");
  require(timing.segment_seconds % parent.time_step_seconds == 0,
          "segment length must be divisible by parent time step");
  require(timing.history_interval_seconds > 0, "history interval must be positive");
  require(timing.history_interval_seconds % parent.time_step_seconds == 0,
          "history interval must be divisible by parent time step");
  require(timing.boundary_refresh_interval_seconds > 0,
          "boundary refresh interval must be positive");
  require(timing.spectral_nudging_input_interval_seconds > 0,
          "spectral nudging input interval must be positive");
  require(timing.moving_nest_interval_seconds > 0,
          "moving nest interval must be positive");
}

[[nodiscard]] CycleScheduleConfig make_cycle_schedule_config(
    const DynamicsLoopConfig& config) {
  const auto& timing = config.timing;
  return CycleScheduleConfig{
      config.parent.grid_spacing_m,
      config.child.grid_spacing_m,
      timing.parent_time_step_seconds,
      timing.child_time_step_seconds,
      timing.parent_time_step_ratio,
      timing.segment_seconds,
      timing.boundary_refresh_interval_seconds,
      timing.spectral_nudging_input_interval_seconds,
      timing.moving_nest_interval_seconds,
      timing.history_interval_seconds,
      SegmentEndpointPolicy::bracket_start_and_end,
      MovingNestTimingPolicy::snap_to_next_parent_step,
  };
}

void count_event(const LoopEventKind kind, LoopSummary& summary) noexcept {
  switch (kind) {
    case LoopEventKind::boundary_input_refresh:
      ++summary.boundary_input_refreshes;
      break;
    case LoopEventKind::spectral_nudging_input_refresh:
      ++summary.spectral_nudging_input_refreshes;
      break;
    case LoopEventKind::boundary_update:
      ++summary.boundary_updates;
      break;
    case LoopEventKind::spectral_nudging:
      ++summary.spectral_nudging_calls;
      break;
    case LoopEventKind::moving_nest_move_check:
      ++summary.moving_nest_move_checks;
      ++summary.moving_nest_position_updates;
      break;
    case LoopEventKind::vortex_center_recompute:
      ++summary.vortex_center_recomputes;
      break;
    case LoopEventKind::moving_nest_position_update:
      ++summary.moving_nest_move_checks;
      ++summary.moving_nest_position_updates;
      break;
    case LoopEventKind::zero_dynamics_tendency:
      ++summary.dynamics_tendency_calls;
      break;
    case LoopEventKind::physics:
      ++summary.physics_calls;
      break;
    case LoopEventKind::nest_interpolation:
      ++summary.nest_interpolations;
      break;
    case LoopEventKind::nest_feedback:
      ++summary.nest_feedbacks;
      break;
    case LoopEventKind::history_output:
      ++summary.history_outputs;
      break;
  }
}

[[nodiscard]] LoopEventKind map_schedule_call_kind(
    const CycleScheduleCallKind kind) noexcept {
  switch (kind) {
    case CycleScheduleCallKind::boundary_input_refresh:
      return LoopEventKind::boundary_input_refresh;
    case CycleScheduleCallKind::spectral_nudging_input_refresh:
      return LoopEventKind::spectral_nudging_input_refresh;
    case CycleScheduleCallKind::d01_boundary_update:
      return LoopEventKind::boundary_update;
    case CycleScheduleCallKind::d01_spectral_nudging:
      return LoopEventKind::spectral_nudging;
    case CycleScheduleCallKind::moving_nest_move_check:
      return LoopEventKind::moving_nest_move_check;
    case CycleScheduleCallKind::vortex_center_recompute:
      return LoopEventKind::vortex_center_recompute;
    case CycleScheduleCallKind::moving_nest_position_update:
      return LoopEventKind::moving_nest_move_check;
    case CycleScheduleCallKind::parent_child_interpolation:
      return LoopEventKind::nest_interpolation;
    case CycleScheduleCallKind::two_way_feedback:
      return LoopEventKind::nest_feedback;
    case CycleScheduleCallKind::history_output:
      return LoopEventKind::history_output;
  }
  return LoopEventKind::zero_dynamics_tendency;
}

void emit_event(const LoopEventSink& sink, const LoopEvent& event, LoopSummary& summary) {
  count_event(event.kind, summary);
  if (sink.callback != nullptr) {
    sink.callback(sink.user_data, event);
  }
}

void emit_schedule_call(
    const LoopEventSink& sink,
    const CycleScheduleCall& call,
    LoopSummary& summary) {
  ++summary.schedule_calls;
  emit_event(
      sink,
      LoopEvent{
          map_schedule_call_kind(call.kind),
          call.domain,
          call.sequence_index,
          call.parent_step_index,
          call.child_substep_index,
          call.start_seconds,
          call.end_seconds,
          call.nominal_seconds,
      },
      summary);
}

void emit_skeleton_dynamics_and_physics(
    const LoopEventSink& sink,
    const CycleScheduleCall& call,
    LoopSummary& summary) {
  emit_event(
      sink,
      LoopEvent{
          LoopEventKind::zero_dynamics_tendency,
          call.domain,
          call.sequence_index,
          call.parent_step_index,
          call.child_substep_index,
          call.start_seconds,
          call.end_seconds,
          call.nominal_seconds,
      },
      summary);
  emit_event(
      sink,
      LoopEvent{
          LoopEventKind::physics,
          call.domain,
          call.sequence_index,
          call.parent_step_index,
          call.child_substep_index,
          call.start_seconds,
          call.end_seconds,
          call.nominal_seconds,
      },
      summary);
}

}  // namespace

DynamicsLoopRunner::DynamicsLoopRunner(DynamicsLoopConfig config) : config_(config) {
  validate_config(config_);
}

const DynamicsLoopConfig& DynamicsLoopRunner::config() const noexcept {
  return config_;
}

LoopSummary DynamicsLoopRunner::run(const LoopEventSink& sink) const {
  LoopSummary summary;
  const auto schedule = CycleSchedule::build(make_cycle_schedule_config(config_));
  summary.parent_steps = schedule.summary().parent_steps;
  summary.child_steps = schedule.summary().child_substeps;

  for (const auto& call : schedule.calls()) {
    emit_schedule_call(sink, call, summary);

    if (call.kind == CycleScheduleCallKind::d01_spectral_nudging ||
        call.kind == CycleScheduleCallKind::parent_child_interpolation) {
      emit_skeleton_dynamics_and_physics(sink, call, summary);
    }
  }

  return summary;
}

DynamicsLoopConfig make_krosa_phase4_loop_config() noexcept {
  return {
      DomainDescriptor{DomainId::d01, 10'000, 40, true, true, false},
      DomainDescriptor{DomainId::d02, 2'000, 8, false, false, true},
      TimeStepDescriptor{40, 8, 5, 21'600, 21'600, 21'600, 21'600, 900},
  };
}

DynamicsLoopConfig make_krosa_10min_validation_loop_config() noexcept {
  return {
      DomainDescriptor{DomainId::d01, 10'000, 40, true, true, false},
      DomainDescriptor{DomainId::d02, 2'000, 8, false, false, true},
      TimeStepDescriptor{40, 8, 5, 600, 600, 21'600, 21'600, 900},
  };
}

std::string_view domain_name(const DomainId domain) noexcept {
  switch (domain) {
    case DomainId::d01:
      return "d01";
    case DomainId::d02:
      return "d02";
  }
  return "unknown";
}

std::string_view loop_event_name(const LoopEventKind kind) noexcept {
  switch (kind) {
    case LoopEventKind::boundary_input_refresh:
      return "boundary_input_refresh";
    case LoopEventKind::spectral_nudging_input_refresh:
      return "spectral_nudging_input_refresh";
    case LoopEventKind::boundary_update:
      return "boundary_update";
    case LoopEventKind::spectral_nudging:
      return "spectral_nudging";
    case LoopEventKind::moving_nest_move_check:
      return "moving_nest_move_check";
    case LoopEventKind::vortex_center_recompute:
      return "vortex_center_recompute";
    case LoopEventKind::moving_nest_position_update:
      return "moving_nest_position_update";
    case LoopEventKind::zero_dynamics_tendency:
      return "zero_dynamics_tendency";
    case LoopEventKind::physics:
      return "physics";
    case LoopEventKind::nest_interpolation:
      return "nest_interpolation";
    case LoopEventKind::nest_feedback:
      return "nest_feedback";
    case LoopEventKind::history_output:
      return "history_output";
  }
  return "unknown";
}

}  // namespace tywrf::dynamics
