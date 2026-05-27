#include "tywrf/dynamics/state_stepper.hpp"

#include "tywrf/grid.hpp"
#include "tywrf/state.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

int failures = 0;

void expect(const bool condition, const std::string_view message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

bool nearly_equal(const float lhs, const float rhs) {
  return std::abs(lhs - rhs) < 1.0e-5F;
}

tywrf::Grid make_grid(const std::int32_t nx, const std::int32_t ny) {
  return tywrf::Grid({
      .mass_nx = nx,
      .mass_ny = ny,
      .mass_nz = 2,
      .full_nz = 3,
      .halo = {1, 1, 1, 1, 0, 0},
  });
}

template <typename Storage>
void fill_storage(Storage& storage, float& next) {
  for (std::size_t idx = 0; idx < storage.size(); ++idx) {
    storage.data()[idx] = next;
    next += 0.25F;
  }
}

void fill_state(tywrf::State<float>& state, const float start) {
  float next = start;
  fill_storage(state.u, next);
  fill_storage(state.v, next);
  fill_storage(state.w, next);
  fill_storage(state.ph, next);
  fill_storage(state.phb, next);
  fill_storage(state.t, next);
  fill_storage(state.p, next);
  fill_storage(state.pb, next);
  fill_storage(state.qvapor, next);
  fill_storage(state.qcloud, next);
  fill_storage(state.qrain, next);
  fill_storage(state.qice, next);
  fill_storage(state.qsnow, next);
  fill_storage(state.qgraup, next);
  fill_storage(state.qnice, next);
  fill_storage(state.qnrain, next);
  fill_storage(state.mu, next);
  fill_storage(state.mub, next);
  fill_storage(state.psfc, next);
  fill_storage(state.u10, next);
  fill_storage(state.v10, next);
  fill_storage(state.t2, next);
  fill_storage(state.q2, next);
  fill_storage(state.rainc, next);
  fill_storage(state.rainnc, next);
}

template <typename Storage>
void fill_storage_constant(Storage& storage, const float value) {
  for (std::size_t idx = 0; idx < storage.size(); ++idx) {
    storage.data()[idx] = value;
  }
}

template <typename Storage>
void append_storage(const Storage& storage, std::vector<float>& values) {
  values.insert(values.end(), storage.data(), storage.data() + storage.size());
}

std::vector<float> snapshot_state(const tywrf::State<float>& state) {
  std::vector<float> values;
  append_storage(state.u, values);
  append_storage(state.v, values);
  append_storage(state.w, values);
  append_storage(state.ph, values);
  append_storage(state.phb, values);
  append_storage(state.t, values);
  append_storage(state.p, values);
  append_storage(state.pb, values);
  append_storage(state.qvapor, values);
  append_storage(state.qcloud, values);
  append_storage(state.qrain, values);
  append_storage(state.qice, values);
  append_storage(state.qsnow, values);
  append_storage(state.qgraup, values);
  append_storage(state.qnice, values);
  append_storage(state.qnrain, values);
  append_storage(state.mu, values);
  append_storage(state.mub, values);
  append_storage(state.psfc, values);
  append_storage(state.u10, values);
  append_storage(state.v10, values);
  append_storage(state.t2, values);
  append_storage(state.q2, values);
  append_storage(state.rainc, values);
  append_storage(state.rainnc, values);
  return values;
}

[[nodiscard]] std::int64_t active_points(const tywrf::FieldLayout2D layout) {
  return static_cast<std::int64_t>(layout.active_nx()) *
         static_cast<std::int64_t>(layout.active_ny());
}

[[nodiscard]] std::int64_t active_points(const tywrf::FieldLayout3D layout) {
  return static_cast<std::int64_t>(layout.active_nx()) *
         static_cast<std::int64_t>(layout.active_ny()) *
         static_cast<std::int64_t>(layout.active_nz());
}

[[nodiscard]] std::int64_t state_active_points(const tywrf::Grid& grid) {
  return active_points(grid.u_layout()) + active_points(grid.v_layout()) +
         3 * active_points(grid.w_layout()) +
         11 * active_points(grid.mass_layout()) +
         9 * active_points(grid.surface_layout());
}

void test_zero_tendency_stepper_preserves_state_and_counts_events() {
  tywrf::State<float> d01(make_grid(4, 3));
  tywrf::State<float> d02(make_grid(5, 4));
  fill_state(d01, 1.0F);
  fill_state(d02, 10'000.0F);
  const auto d01_before = snapshot_state(d01);
  const auto d02_before = snapshot_state(d02);

  const tywrf::dynamics::SkeletonStateStepper stepper(
      tywrf::dynamics::make_krosa_phase4_loop_config());
  const auto report = stepper.run(d01, d02);

  expect(report.status == tywrf::dynamics::StateStepStatus::ok, "stepper status ok");
  expect(!report.layout_or_status_failure, "no layout/status failure");
  expect(report.loop.parent_steps == 540, "d01 parent step count");
  expect(report.loop.child_steps == 2'700, "d02 child substep count");
  expect(report.d01.zero_tendency_apply_count == report.loop.parent_steps,
         "d01 zero-tendency applies once per parent step");
  expect(report.d02.zero_tendency_apply_count == report.loop.child_steps,
         "d02 zero-tendency applies once per child substep");
  expect(report.d01.explicit_tendency_apply_count == 0,
         "default d01 has no explicit tendency applies");
  expect(report.d02.explicit_tendency_apply_count == 0,
         "default d02 has no explicit tendency applies");
  expect(report.d01.tendency_dt_seconds == 21'600,
         "d01 zero tendency reports integrated dt");
  expect(report.d02.tendency_dt_seconds == 21'600,
         "d02 zero tendency reports integrated dt");
  expect(report.d01.zero_tendency_apply_count + report.d02.zero_tendency_apply_count ==
             report.loop.dynamics_tendency_calls,
         "domain apply counts match loop dynamics event count");
  expect(report.loop.physics_calls == report.loop.dynamics_tendency_calls,
         "executor observes existing physics placeholders but does not run physics");

  const auto d01_expected_points =
      state_active_points(d01.grid) * report.d01.zero_tendency_apply_count;
  const auto d02_expected_points =
      state_active_points(d02.grid) * report.d02.zero_tendency_apply_count;
  expect(d01_expected_points > 0, "d01 active points are positive");
  expect(d02_expected_points > 0, "d02 active points are positive");
  expect(report.d01.active_points == d01_expected_points,
         "d01 active points accumulate by apply count");
  expect(report.d02.active_points == d02_expected_points,
         "d02 active points accumulate by apply count");
  expect(report.total_active_points == d01_expected_points + d02_expected_points,
         "total active points accumulate both domains");

  expect(snapshot_state(d01) == d01_before, "d01 state remains unchanged");
  expect(snapshot_state(d02) == d02_before, "d02 state remains unchanged");
}

void test_zero_tendency_reports_600_second_validation_segment_dt() {
  tywrf::State<float> d01(make_grid(4, 3));
  tywrf::State<float> d02(make_grid(5, 4));
  fill_state(d01, 1.0F);
  fill_state(d02, 10'000.0F);
  const auto d01_before = snapshot_state(d01);
  const auto d02_before = snapshot_state(d02);

  const tywrf::dynamics::SkeletonStateStepper stepper(
      tywrf::dynamics::make_krosa_10min_validation_loop_config());
  const auto report = stepper.run(d01, d02);

  expect(report.status == tywrf::dynamics::StateStepStatus::ok,
         "600 s zero tendency status ok");
  expect(report.loop.parent_steps == 15, "600 s d01 parent step count");
  expect(report.loop.child_steps == 75, "600 s d02 child substep count");
  expect(report.d01.zero_tendency_apply_count == 15,
         "600 s d01 zero-tendency apply count");
  expect(report.d02.zero_tendency_apply_count == 75,
         "600 s d02 zero-tendency apply count");
  expect(report.d01.tendency_dt_seconds == 600,
         "600 s d01 zero tendency accumulates event dt");
  expect(report.d02.tendency_dt_seconds == 600,
         "600 s d02 zero tendency accumulates event dt");
  expect(snapshot_state(d01) == d01_before, "600 s d01 zero tendency is no-op");
  expect(snapshot_state(d02) == d02_before, "600 s d02 zero tendency is no-op");
}

void test_explicit_tendency_hook_applies_dt_scaled_state_tendency() {
  const auto d01_grid = make_grid(4, 3);
  tywrf::State<float> d01(d01_grid);
  tywrf::State<float> d02(make_grid(5, 4));
  tywrf::State<float> d01_tendency(d01_grid);
  fill_state(d01, 1.0F);
  fill_state(d02, 10'000.0F);
  fill_storage_constant(d01_tendency.t, 0.25F);
  const auto d02_before = snapshot_state(d02);

  const auto mass = d01_grid.mass_layout();
  const auto active_i = mass.i_begin();
  const auto active_j = mass.j_begin();
  const auto active_k = mass.k_begin();
  const auto active_before = d01.t.view()(active_i, active_j, active_k);
  const auto halo_before = d01.t.view()(0, 0, 0);
  const auto u_before = d01.u.view()(
      d01_grid.u_layout().i_begin(),
      d01_grid.u_layout().j_begin(),
      d01_grid.u_layout().k_begin());

  const tywrf::dynamics::SkeletonStateStepper stepper(
      tywrf::dynamics::make_krosa_10min_validation_loop_config());
  const auto report = stepper.run_with_explicit_tendencies(
      d01,
      d02,
      tywrf::dynamics::ExplicitStateTendencySet{&d01_tendency, nullptr});

  expect(report.status == tywrf::dynamics::StateStepStatus::ok,
         "explicit tendency status ok");
  expect(report.d01.explicit_tendency_apply_count == report.loop.parent_steps,
         "d01 explicit tendency applies once per parent step");
  expect(report.d01.zero_tendency_apply_count == 0,
         "d01 explicit tendency replaces zero tendency");
  expect(report.d01.tendency_dt_seconds == 600,
         "d01 explicit tendency accumulates 600 s dt");
  expect(report.d02.zero_tendency_apply_count == report.loop.child_steps,
         "d02 remains on default zero tendency");
  expect(report.d02.explicit_tendency_apply_count == 0,
         "d02 has no explicit tendency applies");
  expect(report.d02.tendency_dt_seconds == 600,
         "d02 zero tendency still accumulates 600 s dt");

  expect(nearly_equal(d01.t.view()(active_i, active_j, active_k),
                      active_before + 150.0F),
         "d01 active T receives dt-scaled explicit tendency");
  expect(nearly_equal(d01.t.view()(0, 0, 0), halo_before),
         "d01 T halo remains unchanged");
  expect(nearly_equal(d01.u.view()(
                          d01_grid.u_layout().i_begin(),
                          d01_grid.u_layout().j_begin(),
                          d01_grid.u_layout().k_begin()),
                      u_before),
         "d01 unrelated U field remains unchanged");
  expect(snapshot_state(d02) == d02_before, "d02 state remains zero-tendency no-op");
}

void test_missing_d02_state_reports_failure() {
  tywrf::State<float> d01(make_grid(4, 3));
  fill_state(d01, 2.0F);

  const tywrf::dynamics::SkeletonStateStepper stepper(
      tywrf::dynamics::make_krosa_phase4_loop_config());
  const auto report = stepper.run(&d01, nullptr);

  expect(report.layout_or_status_failure, "missing d02 reports a status failure");
  expect(report.status == tywrf::dynamics::StateStepStatus::missing_d02_state,
         "missing d02 status");
  expect(report.failure_domain == tywrf::dynamics::DomainId::d02,
         "missing d02 failure domain");
  expect(report.d02.missing_state, "d02 domain report marks missing state");
  expect(report.d02.failed, "d02 domain report marks failure");
  expect(report.d02.tendency_status == tywrf::dynamics::TendencyApplyStatus::null_field,
         "missing d02 domain tendency status");
  expect(report.d02.tendency_dt_seconds == 0,
         "missing d02 has no applied tendency dt");
  expect(report.d02.zero_tendency_apply_count == 0, "missing d02 is not applied");
  expect(report.d01.zero_tendency_apply_count == report.loop.parent_steps,
         "d01 continues to apply zero tendency");
  expect(report.loop.dynamics_tendency_calls ==
             report.loop.parent_steps + report.loop.child_steps,
         "loop summary still records the schedule");
}

}  // namespace

int main() {
  test_zero_tendency_stepper_preserves_state_and_counts_events();
  test_zero_tendency_reports_600_second_validation_segment_dt();
  test_explicit_tendency_hook_applies_dt_scaled_state_tendency();
  test_missing_d02_state_reports_failure();

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated skeleton state stepper zero-tendency executor\n";
  return 0;
}
