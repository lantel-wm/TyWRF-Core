#include "tywrf/dynamics/cycle_schedule.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using tywrf::dynamics::CycleScheduleCall;
using tywrf::dynamics::CycleScheduleCallKind;
using tywrf::dynamics::DomainId;

int failures = 0;

void expect(const bool condition, const std::string_view message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

void expect_call(
    const CycleScheduleCall& call,
    const CycleScheduleCallKind kind,
    const DomainId domain,
    const std::int64_t parent_step,
    const std::int32_t child_substep,
    const std::int64_t start_seconds,
    const std::int64_t end_seconds,
    const std::int64_t nominal_seconds,
    const std::string_view label) {
  expect(call.kind == kind, std::string(label) + " kind");
  expect(call.domain == domain, std::string(label) + " domain");
  expect(call.parent_step_index == parent_step, std::string(label) + " parent step");
  expect(call.child_substep_index == child_substep,
         std::string(label) + " child substep");
  expect(call.start_seconds == start_seconds, std::string(label) + " start");
  expect(call.end_seconds == end_seconds, std::string(label) + " end");
  expect(call.nominal_seconds == nominal_seconds, std::string(label) + " nominal");
}

const CycleScheduleCall* find_call(
    const std::vector<CycleScheduleCall>& calls,
    const CycleScheduleCallKind kind,
    const std::int64_t nominal_seconds) {
  for (const auto& call : calls) {
    if (call.kind == kind && call.nominal_seconds == nominal_seconds) {
      return &call;
    }
  }
  return nullptr;
}

template <typename Fn>
void expect_invalid(Fn&& fn, const std::string_view label) {
  try {
    fn();
  } catch (const std::invalid_argument&) {
    return;
  }
  expect(false, std::string(label) + " should reject invalid config");
}

}  // namespace

int main() {
  const auto config = tywrf::dynamics::make_krosa_6h_cycle_schedule_config();
  expect(config.parent_grid_spacing_m == 10'000, "d01 grid spacing is 10 km");
  expect(config.child_grid_spacing_m == 2'000, "d02 grid spacing remains 2 km");
  expect(config.parent_time_step_seconds == 40, "d01 step is 40 s");
  expect(config.child_time_step_seconds == 8, "d02 step is 8 s");
  expect(config.parent_time_step_ratio == 5, "parent_time_step_ratio is 5");
  expect(config.segment_seconds == 21'600, "segment is 6 h");
  expect(config.boundary_refresh_interval_seconds == 21'600,
         "boundary refresh interval is 6 h");
  expect(config.spectral_nudging_input_interval_seconds == 21'600,
         "nudging input interval is 6 h");
  expect(config.moving_nest_interval_seconds == 900,
         "moving nest interval is 15 minutes");
  expect(config.history_interval_seconds == 21'600, "history interval is 6 h");

  const auto schedule = tywrf::dynamics::build_krosa_6h_cycle_schedule();
  const auto summary = schedule.summary();
  const std::vector<CycleScheduleCall> calls{
      schedule.calls().begin(), schedule.calls().end()};

  constexpr std::int64_t expected_parent_steps = 540;
  constexpr std::int64_t expected_child_substeps = 2'700;
  constexpr std::int64_t expected_moving_updates = 25;
  constexpr std::int64_t expected_calls =
      2 + 2 + expected_parent_steps * 2 + expected_child_substeps +
      expected_parent_steps + expected_moving_updates + 2;

  expect(summary.parent_steps == expected_parent_steps, "6 h has 540 d01 steps");
  expect(summary.child_substeps == expected_child_substeps,
         "6 h has 2700 d02 substeps");
  expect(summary.boundary_input_refreshes == 2,
         "boundary source is refreshed at segment start and end");
  expect(summary.spectral_nudging_input_refreshes == 2,
         "nudging source is refreshed at segment start and end");
  expect(summary.d01_boundary_updates == expected_parent_steps,
         "d01 boundary update is called each parent step");
  expect(summary.d01_spectral_nudging_calls == expected_parent_steps,
         "d01 spectral nudging is called each parent step");
  expect(summary.moving_nest_position_updates == expected_moving_updates,
         "moving nest has inclusive 15 minute endpoint schedule");
  expect(summary.parent_child_interpolations == expected_child_substeps,
         "parent-child interpolation is called once per d02 substep");
  expect(summary.two_way_feedbacks == expected_parent_steps,
         "two-way feedback is called once per d01 step");
  expect(summary.history_outputs == 2, "one history output call per domain");
  expect(static_cast<std::int64_t>(calls.size()) == expected_calls,
         "cycle schedule call count");

  if (calls.size() >= 9) {
    expect_call(
        calls[0], CycleScheduleCallKind::boundary_input_refresh, DomainId::d01,
        0, -1, 0, 0, 0, "start boundary refresh");
    expect_call(
        calls[1], CycleScheduleCallKind::spectral_nudging_input_refresh,
        DomainId::d01, 0, -1, 0, 0, 0, "start nudging refresh");
    expect_call(
        calls[2], CycleScheduleCallKind::moving_nest_position_update,
        DomainId::d02, 0, -1, 0, 0, 0, "initial moving nest position");
    expect_call(
        calls[3], CycleScheduleCallKind::d01_boundary_update, DomainId::d01,
        0, -1, 0, 40, 0, "first d01 boundary update");
    expect_call(
        calls[4], CycleScheduleCallKind::d01_spectral_nudging, DomainId::d01,
        0, -1, 0, 40, 0, "first d01 spectral nudging");
    expect_call(
        calls[5], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 0, 0, 8, 0, "first d02 interpolation");
    expect_call(
        calls[6], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 1, 8, 16, 8, "second d02 interpolation");
    expect_call(
        calls[7], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 2, 16, 24, 16, "third d02 interpolation");
    expect_call(
        calls[8], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 3, 24, 32, 24, "fourth d02 interpolation");
  }

  const auto* first_nominal_move =
      find_call(calls, CycleScheduleCallKind::moving_nest_position_update, 900);
  expect(first_nominal_move != nullptr, "first 15 minute moving nest call exists");
  if (first_nominal_move != nullptr) {
    expect_call(
        *first_nominal_move, CycleScheduleCallKind::moving_nest_position_update,
        DomainId::d02, 22, -1, 920, 920, 900,
        "900 s moving nest call snaps to next 40 s boundary");
  }

  const auto* final_nominal_move =
      find_call(calls, CycleScheduleCallKind::moving_nest_position_update, 21'600);
  expect(final_nominal_move != nullptr, "final moving nest endpoint call exists");
  if (final_nominal_move != nullptr) {
    expect_call(
        *final_nominal_move, CycleScheduleCallKind::moving_nest_position_update,
        DomainId::d02, 539, -1, 21'600, 21'600, 21'600,
        "final moving nest call");
  }

  if (calls.size() >= 5) {
    const auto first_final_event = calls.size() - 5;
    expect_call(
        calls[first_final_event],
        CycleScheduleCallKind::moving_nest_position_update, DomainId::d02,
        539, -1, 21'600, 21'600, 21'600, "final movement before output");
    expect_call(
        calls[first_final_event + 1], CycleScheduleCallKind::history_output,
        DomainId::d01, 539, -1, 21'600, 21'600, 21'600, "d01 history");
    expect_call(
        calls[first_final_event + 2], CycleScheduleCallKind::history_output,
        DomainId::d02, 539, -1, 21'600, 21'600, 21'600, "d02 history");
    expect_call(
        calls[first_final_event + 3],
        CycleScheduleCallKind::boundary_input_refresh, DomainId::d01, 539, -1,
        21'600, 21'600, 21'600, "end boundary refresh");
  }

  if (calls.size() >= 2) {
    expect_call(
        calls[calls.size() - 1],
        CycleScheduleCallKind::spectral_nudging_input_refresh, DomainId::d01,
        539, -1, 21'600, 21'600, 21'600, "end nudging refresh");
  }

  expect(
      tywrf::dynamics::cycle_schedule_call_name(
          CycleScheduleCallKind::two_way_feedback) == "two_way_feedback",
      "call kind name");
  expect(
      tywrf::dynamics::segment_endpoint_policy_name(
          tywrf::dynamics::SegmentEndpointPolicy::bracket_start_and_end) ==
          "bracket_start_and_end",
      "endpoint policy name");
  expect(
      tywrf::dynamics::moving_nest_timing_policy_name(
          tywrf::dynamics::MovingNestTimingPolicy::snap_to_next_parent_step) ==
          "snap_to_next_parent_step",
      "moving nest policy name");

  expect_invalid(
      [] {
        auto invalid = tywrf::dynamics::make_krosa_6h_cycle_schedule_config();
        invalid.child_grid_spacing_m = 3'000;
        (void)tywrf::dynamics::CycleSchedule::build(invalid);
      },
      "d02 resolution guard");
  expect_invalid(
      [] {
        auto invalid = tywrf::dynamics::make_krosa_6h_cycle_schedule_config();
        invalid.child_time_step_seconds = 10;
        (void)tywrf::dynamics::CycleSchedule::build(invalid);
      },
      "parent child ratio guard");
  expect_invalid(
      [] {
        auto invalid = tywrf::dynamics::make_krosa_6h_cycle_schedule_config();
        invalid.boundary_refresh_interval_seconds = 10'800;
        (void)tywrf::dynamics::CycleSchedule::build(invalid);
      },
      "endpoint bracketing guard");

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated KROSA 6 h boundary/nudging/nesting schedule contract\n";
  return 0;
}
