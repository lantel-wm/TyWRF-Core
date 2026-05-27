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
    const std::string_view label,
    const std::int64_t schedule_sequence = -1,
    const std::int64_t nominal_seconds = -1) {
  expect(event.kind == kind, std::string(label) + " kind");
  expect(event.domain == domain, std::string(label) + " domain");
  expect(event.parent_step_index == parent_step, std::string(label) + " parent step");
  expect(event.child_substep_index == child_substep, std::string(label) + " child substep");
  expect(event.start_seconds == start_seconds, std::string(label) + " start time");
  expect(event.end_seconds == end_seconds, std::string(label) + " end time");
  if (schedule_sequence >= 0) {
    expect(event.schedule_sequence_index == schedule_sequence,
           std::string(label) + " schedule sequence");
  }
  if (nominal_seconds >= 0) {
    expect(event.nominal_seconds == nominal_seconds, std::string(label) + " nominal time");
  }
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
  constexpr std::int64_t expected_moving_nest_updates = 25;
  constexpr std::int64_t expected_schedule_calls =
      2 + 2 + expected_parent_steps * 2 + expected_child_steps +
      expected_parent_steps + expected_moving_nest_updates + 2;
  constexpr std::int64_t expected_events =
      expected_schedule_calls + 2 * (expected_parent_steps + expected_child_steps);

  expect(summary.parent_steps == expected_parent_steps, "6 h segment has 540 d01 steps");
  expect(summary.child_steps == expected_child_steps, "6 h segment has 2700 d02 steps");
  expect(summary.schedule_calls == expected_schedule_calls,
         "runner traverses every KROSA 6 h schedule call");
  expect(summary.boundary_input_refreshes == 2,
         "boundary input refreshes bracket the segment");
  expect(summary.spectral_nudging_input_refreshes == 2,
         "spectral nudging input refreshes bracket the segment");
  expect(summary.boundary_updates == expected_parent_steps, "d01 boundary update count");
  expect(summary.spectral_nudging_calls == expected_parent_steps, "d01 spectral nudging count");
  expect(summary.moving_nest_position_updates == expected_moving_nest_updates,
         "moving nest schedule count");
  expect(summary.dynamics_tendency_calls == expected_parent_steps + expected_child_steps,
         "zero dynamics tendency count");
  expect(summary.physics_calls == expected_parent_steps + expected_child_steps, "physics count");
  expect(summary.nest_interpolations == expected_child_steps, "nest interpolation count");
  expect(summary.nest_feedbacks == expected_parent_steps, "nest feedback count");
  expect(summary.history_outputs == 2, "one 6 h history output per domain");
  expect(static_cast<std::int64_t>(events.size()) == expected_events, "recorded event count");

  if (events.size() >= 11) {
    expect_event(events[0], LoopEventKind::boundary_input_refresh, DomainId::d01, 0, -1, 0,
                 0, "start boundary input refresh", 0, 0);
    expect_event(events[1], LoopEventKind::spectral_nudging_input_refresh, DomainId::d01, 0,
                 -1, 0, 0, "start spectral nudging input refresh", 1, 0);
    expect_event(events[2], LoopEventKind::moving_nest_position_update, DomainId::d02, 0, -1,
                 0, 0, "initial moving nest position update", 2, 0);
    expect_event(events[3], LoopEventKind::boundary_update, DomainId::d01, 0, -1, 0, 40,
                 "first scheduled d01 boundary update", 3, 0);
    expect_event(events[4], LoopEventKind::spectral_nudging, DomainId::d01, 0, -1, 0, 40,
                 "first scheduled d01 spectral nudging", 4, 0);
    expect_event(events[5], LoopEventKind::zero_dynamics_tendency, DomainId::d01, 0, -1, 0,
                 40, "first d01 zero tendency placeholder", 4, 0);
    expect_event(events[6], LoopEventKind::physics, DomainId::d01, 0, -1, 0, 40,
                 "first d01 physics placeholder", 4, 0);
    expect_event(events[7], LoopEventKind::nest_interpolation, DomainId::d02, 0, 0, 0, 8,
                 "first scheduled d02 interpolation", 5, 0);
    expect_event(events[8], LoopEventKind::zero_dynamics_tendency, DomainId::d02, 0, 0, 0,
                 8, "first d02 zero tendency placeholder", 5, 0);
    expect_event(events[9], LoopEventKind::physics, DomainId::d02, 0, 0, 0, 8,
                 "first d02 physics placeholder", 5, 0);
    expect_event(events[10], LoopEventKind::nest_interpolation, DomainId::d02, 0, 1, 8, 16,
                 "second scheduled d02 interpolation", 6, 8);
  }

  if (events.size() >= 5) {
    const auto moving_nest_index = events.size() - 5;
    const auto history_d01_index = events.size() - 4;
    const auto history_d02_index = events.size() - 3;
    const auto boundary_refresh_index = events.size() - 2;
    const auto nudging_refresh_index = events.size() - 1;

    expect_event(events[moving_nest_index], LoopEventKind::moving_nest_position_update,
                 DomainId::d02, 539, -1, 21'600, 21'600,
                 "final moving nest update before output", expected_schedule_calls - 5,
                 21'600);
    expect_event(events[history_d01_index], LoopEventKind::history_output, DomainId::d01, 539,
                 -1, 21'600, 21'600, "d01 6 h history output",
                 expected_schedule_calls - 4, 21'600);
    expect_event(events[history_d02_index], LoopEventKind::history_output, DomainId::d02, 539,
                 -1, 21'600, 21'600, "d02 6 h history output",
                 expected_schedule_calls - 3, 21'600);
    expect_event(events[boundary_refresh_index], LoopEventKind::boundary_input_refresh,
                 DomainId::d01, 539, -1, 21'600, 21'600,
                 "end boundary input refresh", expected_schedule_calls - 2, 21'600);
    expect_event(events[nudging_refresh_index],
                 LoopEventKind::spectral_nudging_input_refresh, DomainId::d01, 539,
                 -1, 21'600, 21'600, "end spectral nudging input refresh",
                 expected_schedule_calls - 1, 21'600);
  }

  expect(tywrf::dynamics::loop_event_name(LoopEventKind::moving_nest_position_update) ==
             "moving_nest_position_update",
         "moving nest loop event name");

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated schedule-driven dynamics skeleton for a 6 h KROSA segment\n";
  return 0;
}
