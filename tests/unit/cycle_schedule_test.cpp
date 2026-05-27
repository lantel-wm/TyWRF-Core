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

std::vector<const CycleScheduleCall*> find_calls_by_kind(
    const std::vector<CycleScheduleCall>& calls,
    const CycleScheduleCallKind kind) {
  std::vector<const CycleScheduleCall*> matches;
  for (const auto& call : calls) {
    if (call.kind == kind) {
      matches.push_back(&call);
    }
  }
  return matches;
}

const CycleScheduleCall* find_parent_step_call(
    const std::vector<CycleScheduleCall>& calls,
    const CycleScheduleCallKind kind,
    const std::int64_t parent_step) {
  for (const auto& call : calls) {
    if (call.kind == kind && call.parent_step_index == parent_step) {
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
  constexpr std::int64_t expected_move_checks = expected_parent_steps;
  constexpr std::int64_t expected_vortex_recomputes = 25;
  constexpr std::int64_t expected_calls =
      2 + expected_parent_steps * 4 + expected_child_substeps +
      expected_vortex_recomputes + 2 + 2;

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
  expect(summary.moving_nest_move_checks == expected_move_checks,
         "moving nest move check is called once per parent step");
  expect(summary.vortex_center_recomputes == expected_vortex_recomputes,
         "vortex center recompute has inclusive 15 minute endpoint schedule");
  expect(summary.parent_child_interpolations == expected_child_substeps,
         "parent-child interpolation is called once per d02 substep");
  expect(summary.two_way_feedbacks == expected_parent_steps,
         "two-way feedback is called once per d01 step");
  expect(summary.history_outputs == 2, "one history output call per domain");
  expect(static_cast<std::int64_t>(calls.size()) == expected_calls,
         "cycle schedule call count");

  if (calls.size() >= 11) {
    expect_call(
        calls[0], CycleScheduleCallKind::boundary_input_refresh, DomainId::d01,
        0, -1, 0, 0, 0, "start boundary refresh");
    expect_call(
        calls[1], CycleScheduleCallKind::spectral_nudging_input_refresh,
        DomainId::d01, 0, -1, 0, 0, 0, "start nudging refresh");
    expect_call(
        calls[2], CycleScheduleCallKind::d01_boundary_update, DomainId::d01,
        0, -1, 0, 40, 0, "first d01 boundary update");
    expect_call(
        calls[3], CycleScheduleCallKind::d01_spectral_nudging, DomainId::d01,
        0, -1, 0, 40, 0, "first d01 spectral nudging");
    expect_call(
        calls[4], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 0, 0, 8, 0, "first d02 interpolation");
    expect_call(
        calls[5], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 1, 8, 16, 8, "second d02 interpolation");
    expect_call(
        calls[6], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 2, 16, 24, 16, "third d02 interpolation");
    expect_call(
        calls[7], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 3, 24, 32, 24, "fourth d02 interpolation");
    expect_call(
        calls[8], CycleScheduleCallKind::parent_child_interpolation,
        DomainId::d02, 0, 4, 32, 40, 32, "fifth d02 interpolation");
    expect_call(
        calls[9], CycleScheduleCallKind::two_way_feedback, DomainId::d01,
        0, -1, 0, 40, 40, "first two-way feedback");
    expect_call(
        calls[10], CycleScheduleCallKind::moving_nest_move_check,
        DomainId::d02, 0, -1, 40, 40, 40,
        "first moving nest move check");
  }

  const auto move_checks =
      find_calls_by_kind(calls, CycleScheduleCallKind::moving_nest_move_check);
  const auto legacy_position_updates =
      find_calls_by_kind(calls, CycleScheduleCallKind::moving_nest_position_update);
  const auto vortex_recomputes =
      find_calls_by_kind(calls, CycleScheduleCallKind::vortex_center_recompute);
  expect(static_cast<std::int64_t>(move_checks.size()) == expected_move_checks,
         "moving nest move check count matches parent steps");
  expect(static_cast<std::int64_t>(legacy_position_updates.size()) == 0,
         "legacy position-update alias is absent from the default remap schedule");
  expect(static_cast<std::int64_t>(vortex_recomputes.size()) ==
             expected_vortex_recomputes,
         "15 minute events are vortex recomputes, not legacy position updates");
  expect(summary.parent_child_interpolations == expected_child_substeps,
         "d02 remap/interpolation sequencing does not require legacy position updates");
  if (!move_checks.empty()) {
    expect_call(
        *move_checks.front(), CycleScheduleCallKind::moving_nest_move_check,
        DomainId::d02, 0, -1, 40, 40, 40,
        "first actual moving nest check is at the first parent step end");
  }
  for (const auto* move_check : move_checks) {
    const auto* feedback = find_parent_step_call(
        calls, CycleScheduleCallKind::two_way_feedback,
        move_check->parent_step_index);
    expect(feedback != nullptr, "feedback exists before moving nest check");
    if (feedback != nullptr) {
      expect(feedback->sequence_index < move_check->sequence_index,
             "feedback precedes moving nest check within parent step");
      expect(feedback->end_seconds == move_check->end_seconds,
             "feedback and moving nest check share the same parent step end");
    }
    expect(move_check->end_seconds > 0, "moving nest check is not emitted at t=0");
  }

  const auto* initial_recompute =
      find_call(calls, CycleScheduleCallKind::vortex_center_recompute, 0);
  expect(initial_recompute != nullptr, "initial vortex recompute call exists");
  if (initial_recompute != nullptr) {
    expect_call(
        *initial_recompute, CycleScheduleCallKind::vortex_center_recompute,
        DomainId::d02, 0, -1, 40, 40, 0,
        "nominal 0 s vortex recompute snaps to first parent step end");
  }

  const auto* first_nominal_recompute =
      find_call(calls, CycleScheduleCallKind::vortex_center_recompute, 900);
  expect(first_nominal_recompute != nullptr,
         "first 15 minute vortex recompute call exists");
  if (first_nominal_recompute != nullptr) {
    expect_call(
        *first_nominal_recompute, CycleScheduleCallKind::vortex_center_recompute,
        DomainId::d02, 22, -1, 920, 920, 900,
        "900 s vortex recompute snaps to next 40 s boundary");
  }

  const auto* second_nominal_recompute =
      find_call(calls, CycleScheduleCallKind::vortex_center_recompute, 1'800);
  expect(second_nominal_recompute != nullptr,
         "30 minute vortex recompute call exists");
  if (second_nominal_recompute != nullptr) {
    expect_call(
        *second_nominal_recompute,
        CycleScheduleCallKind::vortex_center_recompute, DomainId::d02, 44, -1,
        1'800, 1'800, 1'800,
        "1800 s vortex recompute lands on exact parent step end");
  }

  const auto* final_nominal_recompute =
      find_call(calls, CycleScheduleCallKind::vortex_center_recompute, 21'600);
  expect(final_nominal_recompute != nullptr,
         "final vortex recompute endpoint call exists");
  if (final_nominal_recompute != nullptr) {
    expect_call(
        *final_nominal_recompute,
        CycleScheduleCallKind::vortex_center_recompute, DomainId::d02, 539, -1,
        21'600, 21'600, 21'600, "final vortex recompute call");
  }

  if (calls.size() >= 7) {
    const auto first_final_event = calls.size() - 7;
    expect_call(
        calls[first_final_event], CycleScheduleCallKind::two_way_feedback,
        DomainId::d01, 539, -1, 21'560, 21'600, 21'600,
        "final feedback before moving nest check");
    expect_call(
        calls[first_final_event + 1], CycleScheduleCallKind::moving_nest_move_check,
        DomainId::d02, 539, -1, 21'600, 21'600, 21'600,
        "final moving nest check before output");
    expect_call(
        calls[first_final_event + 2],
        CycleScheduleCallKind::vortex_center_recompute, DomainId::d02, 539, -1,
        21'600, 21'600, 21'600, "final vortex recompute before output");
    expect_call(
        calls[first_final_event + 3], CycleScheduleCallKind::history_output,
        DomainId::d01, 539, -1, 21'600, 21'600, 21'600, "d01 history");
    expect_call(
        calls[first_final_event + 4], CycleScheduleCallKind::history_output,
        DomainId::d02, 539, -1, 21'600, 21'600, 21'600, "d02 history");
    expect_call(
        calls[first_final_event + 5],
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
          CycleScheduleCallKind::moving_nest_move_check) ==
          "moving_nest_move_check",
      "call kind name");
  expect(
      tywrf::dynamics::cycle_schedule_call_name(
          CycleScheduleCallKind::vortex_center_recompute) ==
          "vortex_center_recompute",
      "vortex recompute call kind name");
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

  auto validation_config = tywrf::dynamics::make_krosa_6h_cycle_schedule_config();
  validation_config.segment_seconds = 600;
  validation_config.history_interval_seconds = 600;
  const auto validation_schedule =
      tywrf::dynamics::CycleSchedule::build(validation_config);
  const auto validation_summary = validation_schedule.summary();
  const std::vector<CycleScheduleCall> validation_calls{
      validation_schedule.calls().begin(), validation_schedule.calls().end()};

  expect(validation_summary.parent_steps == 15,
         "10 min validation segment has 15 d01 steps");
  expect(validation_summary.child_substeps == 75,
         "10 min validation segment has 75 d02 substeps");
  expect(validation_summary.boundary_input_refreshes == 1,
         "10 min validation segment does not refresh boundary input at 600 s");
  expect(validation_summary.spectral_nudging_input_refreshes == 1,
         "10 min validation segment does not refresh FDDA input at 600 s");
  expect(validation_summary.history_outputs == 2,
         "10 min validation segment writes history output at 600 s for both domains");
  expect(validation_summary.moving_nest_move_checks == 15,
         "10 min validation segment has one moving nest check per parent step");
  expect(validation_summary.vortex_center_recomputes == 1,
         "10 min validation segment has only the snapped initial vortex recompute");

  const auto validation_boundary_refreshes = find_calls_by_kind(
      validation_calls, CycleScheduleCallKind::boundary_input_refresh);
  const auto validation_nudging_refreshes = find_calls_by_kind(
      validation_calls, CycleScheduleCallKind::spectral_nudging_input_refresh);
  const auto validation_history_outputs =
      find_calls_by_kind(validation_calls, CycleScheduleCallKind::history_output);
  const auto validation_vortex_recomputes = find_calls_by_kind(
      validation_calls, CycleScheduleCallKind::vortex_center_recompute);
  expect(validation_boundary_refreshes.size() == 1,
         "10 min validation segment emits only the initial boundary refresh");
  expect(validation_nudging_refreshes.size() == 1,
         "10 min validation segment emits only the initial FDDA refresh");
  expect(validation_history_outputs.size() == 2,
         "10 min validation segment emits d01 and d02 history outputs");
  expect(validation_vortex_recomputes.size() == 1,
         "10 min validation segment emits one vortex recompute call");
  if (!validation_boundary_refreshes.empty()) {
    expect_call(
        *validation_boundary_refreshes.front(),
        CycleScheduleCallKind::boundary_input_refresh, DomainId::d01, 0, -1,
        0, 0, 0, "10 min initial boundary refresh");
  }
  if (!validation_nudging_refreshes.empty()) {
    expect_call(
        *validation_nudging_refreshes.front(),
        CycleScheduleCallKind::spectral_nudging_input_refresh, DomainId::d01,
        0, -1, 0, 0, 0, "10 min initial nudging refresh");
  }
  for (const auto* history_output : validation_history_outputs) {
    expect(history_output->end_seconds == 600,
           "10 min history output is emitted at 600 s");
  }
  if (!validation_vortex_recomputes.empty()) {
    expect_call(
        *validation_vortex_recomputes.front(),
        CycleScheduleCallKind::vortex_center_recompute, DomainId::d02, 0, -1,
        40, 40, 0, "10 min nominal 0 s vortex recompute");
  }
  expect(find_call(validation_calls, CycleScheduleCallKind::vortex_center_recompute,
                   900) == nullptr,
         "10 min validation segment does not emit the nominal 900 s recompute");

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
        invalid.boundary_refresh_interval_seconds = 0;
        (void)tywrf::dynamics::CycleSchedule::build(invalid);
      },
      "boundary refresh interval guard");

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated KROSA 6 h boundary/nudging/nesting schedule contract\n";
  return 0;
}
