#include "tywrf/dynamics/dynamics_loop.hpp"

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
  require(parent.grid_spacing_m > 0, "parent grid spacing must be positive");
  require(child.grid_spacing_m > 0, "child grid spacing must be positive");
  require(parent.time_step_seconds > 0, "parent time step must be positive");
  require(child.time_step_seconds > 0, "child time step must be positive");
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
}

void count_event(const LoopEventKind kind, LoopSummary& summary) noexcept {
  switch (kind) {
    case LoopEventKind::boundary_update:
      ++summary.boundary_updates;
      break;
    case LoopEventKind::spectral_nudging:
      ++summary.spectral_nudging_calls;
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

void emit_event(const LoopEventSink& sink, const LoopEvent& event, LoopSummary& summary) {
  count_event(event.kind, summary);
  if (sink.callback != nullptr) {
    sink.callback(sink.user_data, event);
  }
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
  const auto& timing = config_.timing;

  const auto parent_steps =
      static_cast<std::int64_t>(timing.segment_seconds / timing.parent_time_step_seconds);

  for (std::int64_t parent_step = 0; parent_step < parent_steps; ++parent_step) {
    const auto parent_start =
        parent_step * static_cast<std::int64_t>(timing.parent_time_step_seconds);
    const auto parent_end = parent_start + timing.parent_time_step_seconds;

    if (config_.parent.has_lateral_boundary) {
      emit_event(
          sink,
          {LoopEventKind::boundary_update, DomainId::d01, parent_step, -1, parent_start,
           parent_end},
          summary);
    }

    if (config_.parent.has_spectral_nudging) {
      emit_event(
          sink,
          {LoopEventKind::spectral_nudging, DomainId::d01, parent_step, -1, parent_start,
           parent_end},
          summary);
    }

    emit_event(
        sink,
        {LoopEventKind::zero_dynamics_tendency, DomainId::d01, parent_step, -1,
         parent_start, parent_end},
        summary);
    emit_event(
        sink,
        {LoopEventKind::physics, DomainId::d01, parent_step, -1, parent_start, parent_end},
        summary);

    for (std::int32_t child_substep = 0; child_substep < timing.parent_time_step_ratio;
         ++child_substep) {
      const auto child_start =
          parent_start + child_substep * static_cast<std::int64_t>(timing.child_time_step_seconds);
      const auto child_end = child_start + timing.child_time_step_seconds;

      emit_event(
          sink,
          {LoopEventKind::nest_interpolation, DomainId::d02, parent_step, child_substep,
           child_start, child_end},
          summary);
      emit_event(
          sink,
          {LoopEventKind::zero_dynamics_tendency, DomainId::d02, parent_step, child_substep,
           child_start, child_end},
          summary);
      emit_event(
          sink,
          {LoopEventKind::physics, DomainId::d02, parent_step, child_substep, child_start,
           child_end},
          summary);

      ++summary.child_steps;
    }

    emit_event(
        sink,
        {LoopEventKind::nest_feedback, DomainId::d01, parent_step, -1, parent_start,
         parent_end},
        summary);

    if (parent_end % timing.history_interval_seconds == 0) {
      emit_event(
          sink,
          {LoopEventKind::history_output, DomainId::d01, parent_step, -1, parent_end,
           parent_end},
          summary);
      emit_event(
          sink,
          {LoopEventKind::history_output, DomainId::d02, parent_step, -1, parent_end,
           parent_end},
          summary);
    }

    ++summary.parent_steps;
  }

  return summary;
}

DynamicsLoopConfig make_krosa_phase4_loop_config() noexcept {
  return {
      DomainDescriptor{DomainId::d01, 10'000, 40, true, true, false},
      DomainDescriptor{DomainId::d02, 2'000, 8, false, false, true},
      TimeStepDescriptor{40, 8, 5, 21'600, 21'600, 21'600, 21'600},
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
    case LoopEventKind::boundary_update:
      return "boundary_update";
    case LoopEventKind::spectral_nudging:
      return "spectral_nudging";
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
