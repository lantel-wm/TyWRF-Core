#include "tywrf/dynamics/pressure_refresh_staging.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <type_traits>
#include <vector>

namespace {

tywrf::Grid make_grid() {
  return tywrf::Grid({4, 3, 2, 3, tywrf::uniform_halo_3d(1)});
}

tywrf::io::KrosaPressureRefreshMetadata make_metadata() {
  tywrf::io::KrosaPressureRefreshMetadata metadata{};
  metadata.p_top_pa = 5000.0F;
  metadata.p_top_source = tywrf::io::PressureRefreshPTopSource::time_variable;
  metadata.c3f = {1.0F, 0.5F, 0.0F};
  metadata.c4f = {0.0F, 0.0F, 0.0F};
  metadata.c3h = {0.75F, 0.25F};
  metadata.c4h = {0.0F, 0.0F};
  return metadata;
}

void test_success_stages_existing_buffers() {
  auto grid = make_grid();
  tywrf::State<float> state(grid);
  tywrf::FieldStorage3D<float> external_alb(grid.mass_layout());
  auto metadata = make_metadata();

  std::fill(state.p.data(), state.p.data() + state.p.size(), 1234.0F);
  std::fill(external_alb.data(), external_alb.data() + external_alb.size(), 0.75F);
  const std::vector<float> p_before(state.p.data(), state.p.data() + state.p.size());

  const auto staging = tywrf::dynamics::make_krosa_pressure_refresh_inputs(
      state,
      external_alb,
      metadata);

  assert(staging.ok());
  assert(staging.report.expected_mass_level_count == 2);
  assert(staging.report.expected_full_level_count == 3);
  assert(staging.report.c3f_count == 3);
  assert(staging.report.c3h_count == 2);
  assert(staging.report.alb_is_external);
  assert(staging.inputs.p.data == state.p.data());
  assert(staging.inputs.pb.data == state.pb.data());
  assert(staging.inputs.t.data == state.t.data());
  assert(staging.inputs.alb.data == external_alb.data());
  assert(staging.inputs.alb.data != state.pb.data());
  assert(staging.inputs.c3f.values == metadata.c3f.data());
  assert(staging.inputs.c4f.values == metadata.c4f.data());
  assert(staging.inputs.c3h.values == metadata.c3h.data());
  assert(staging.inputs.c4h.values == metadata.c4h.data());
  assert(staging.inputs.c3f.count == 3);
  assert(staging.inputs.c3h.count == 2);
  assert(staging.inputs.p_top_pa == metadata.p_top_pa);
  assert(std::equal(state.p.data(), state.p.data() + state.p.size(), p_before.begin()));
}

void test_missing_coefficients_fail_without_inputs() {
  auto grid = make_grid();
  tywrf::State<float> state(grid);
  tywrf::FieldStorage3D<float> external_alb(grid.mass_layout());
  auto metadata = make_metadata();
  metadata.c4h.clear();

  const auto staging = tywrf::dynamics::make_krosa_pressure_refresh_inputs(
      state,
      external_alb,
      metadata);

  assert(!staging.ok());
  assert(
      staging.report.result.status ==
      tywrf::nest::NestStatus::invalid_contract);
  assert(staging.report.c4h_count == 0);
  assert(staging.inputs.p.data == nullptr);
  assert(staging.inputs.alb.data == nullptr);
}

void test_bad_external_alb_shape_fails() {
  auto grid = make_grid();
  tywrf::State<float> state(grid);
  tywrf::FieldStorage3D<float> bad_alb(
      tywrf::ActiveShape3D{5, 3, 2},
      tywrf::uniform_halo_3d(1));
  const auto metadata = make_metadata();

  const auto report = tywrf::dynamics::validate_krosa_pressure_refresh_staging(
      state,
      bad_alb,
      metadata);

  assert(!report.ok());
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
}

void test_alb_must_be_external_and_p_top_present() {
  auto grid = make_grid();
  tywrf::State<float> state(grid);
  tywrf::FieldStorage3D<float> external_alb(grid.mass_layout());
  auto metadata = make_metadata();
  metadata.p_top_source = tywrf::io::PressureRefreshPTopSource::missing;

  const auto staging = tywrf::dynamics::make_krosa_pressure_refresh_inputs(
      state,
      external_alb,
      metadata);

  assert(!staging.ok());
  assert(staging.report.alb_is_external);
  assert(staging.report.result.status == tywrf::nest::NestStatus::invalid_contract);
}

void test_staging_does_not_compute_or_modify_pressure() {
  auto grid = make_grid();
  tywrf::State<float> state(grid);
  tywrf::FieldStorage3D<float> external_alb(grid.mass_layout());
  auto metadata = make_metadata();

  std::fill(state.p.data(), state.p.data() + state.p.size(), -9999.0F);
  const auto staging = tywrf::dynamics::make_krosa_pressure_refresh_inputs(
      state,
      external_alb,
      metadata);

  assert(staging.ok());
  for (std::size_t index = 0; index < state.p.size(); ++index) {
    assert(state.p.data()[index] == -9999.0F);
  }
}

}  // namespace

int main() {
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::PressureRefreshStagingReport>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::PressureRefreshStagingReport>);
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::PressureRefreshStagingResult>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::PressureRefreshStagingResult>);

  test_success_stages_existing_buffers();
  test_missing_coefficients_fail_without_inputs();
  test_bad_external_alb_shape_fails();
  test_alb_must_be_external_and_p_top_present();
  test_staging_does_not_compute_or_modify_pressure();

  std::cout << "Validated KROSA pressure refresh staging contract\n";
  return 0;
}
