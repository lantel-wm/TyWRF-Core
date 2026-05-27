#include "tywrf/dynamics/dynamics_loop.hpp"

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

using tywrf::dynamics::DomainId;
using tywrf::dynamics::LoopEvent;
using tywrf::dynamics::LoopEventKind;

int failures = 0;

void expect(const bool condition, const std::string_view message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

void record_event(void* user_data, const LoopEvent& event) {
  auto* events = static_cast<std::vector<LoopEvent>*>(user_data);
  events->push_back(event);
}

void expect_event(
    const LoopEvent& event,
    const LoopEventKind kind,
    const DomainId domain,
    const std::int64_t parent_step,
    const std::int32_t child_substep,
    const std::int64_t start_seconds,
    const std::int64_t end_seconds,
    const std::string_view label) {
  expect(event.kind == kind, std::string(label) + " kind");
  expect(event.domain == domain, std::string(label) + " domain");
  expect(event.parent_step_index == parent_step, std::string(label) + " parent step");
  expect(event.child_substep_index == child_substep, std::string(label) + " child substep");
  expect(event.start_seconds == start_seconds, std::string(label) + " start time");
  expect(event.end_seconds == end_seconds, std::string(label) + " end time");
}

}  // namespace

int main() {
  const auto config = tywrf::dynamics::make_krosa_phase4_loop_config();
  expect(config.parent.grid_spacing_m == 10'000, "d01 grid spacing remains 10 km");
  expect(config.child.grid_spacing_m == 2'000, "d02 grid spacing remains 2 km");
  expect(config.timing.parent_time_step_seconds == 40, "d01 dt is 40 s");
  expect(config.timing.child_time_step_seconds == 8, "d02 dt is 8 s");
  expect(config.timing.parent_time_step_ratio == 5, "d02 subcycles five times per d01 step");
  expect(config.timing.segment_seconds == 21'600, "segment length is 6 h");

  std::vector<LoopEvent> events;
  const tywrf::dynamics::DynamicsLoopRunner runner(config);
  const auto summary = runner.run({&events, record_event});

  constexpr std::int64_t expected_parent_steps = 540;
  constexpr std::int64_t expected_child_steps = 2'700;
  constexpr std::int64_t expected_events = 10'802;

  expect(summary.parent_steps == expected_parent_steps, "6 h segment has 540 d01 steps");
  expect(summary.child_steps == expected_child_steps, "6 h segment has 2700 d02 steps");
  expect(summary.boundary_updates == expected_parent_steps, "d01 boundary update count");
  expect(summary.spectral_nudging_calls == expected_parent_steps, "d01 spectral nudging count");
  expect(summary.dynamics_tendency_calls == expected_parent_steps + expected_child_steps,
         "zero dynamics tendency count");
  expect(summary.physics_calls == expected_parent_steps + expected_child_steps, "physics count");
  expect(summary.nest_interpolations == expected_child_steps, "nest interpolation count");
  expect(summary.nest_feedbacks == expected_parent_steps, "nest feedback count");
  expect(summary.history_outputs == 2, "one 6 h history output per domain");
  expect(static_cast<std::int64_t>(events.size()) == expected_events, "recorded event count");

  if (events.size() >= 8) {
    expect_event(events[0], LoopEventKind::boundary_update, DomainId::d01, 0, -1, 0, 40,
                 "first event");
    expect_event(events[1], LoopEventKind::spectral_nudging, DomainId::d01, 0, -1, 0, 40,
                 "second event");
    expect_event(events[2], LoopEventKind::zero_dynamics_tendency, DomainId::d01, 0, -1, 0,
                 40, "third event");
    expect_event(events[3], LoopEventKind::physics, DomainId::d01, 0, -1, 0, 40,
                 "fourth event");
    expect_event(events[4], LoopEventKind::nest_interpolation, DomainId::d02, 0, 0, 0, 8,
                 "first d02 substep interpolation");
    expect_event(events[5], LoopEventKind::zero_dynamics_tendency, DomainId::d02, 0, 0, 0,
                 8, "first d02 substep tendency");
    expect_event(events[6], LoopEventKind::physics, DomainId::d02, 0, 0, 0, 8,
                 "first d02 substep physics");
    expect_event(events[7], LoopEventKind::nest_interpolation, DomainId::d02, 0, 1, 8, 16,
                 "second d02 substep interpolation");
  }

  if (events.size() >= 3) {
    const auto feedback_index = events.size() - 3;
    const auto history_d01_index = events.size() - 2;
    const auto history_d02_index = events.size() - 1;

    expect_event(events[feedback_index], LoopEventKind::nest_feedback, DomainId::d01, 539, -1,
                 21'560, 21'600, "final feedback before output");
    expect_event(events[history_d01_index], LoopEventKind::history_output, DomainId::d01, 539,
                 -1, 21'600, 21'600, "d01 6 h history output");
    expect_event(events[history_d02_index], LoopEventKind::history_output, DomainId::d02, 539,
                 -1, 21'600, 21'600, "d02 6 h history output");
  }

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated Phase 4 dynamics loop ordering for a 6 h KROSA segment\n";
  return 0;
}
