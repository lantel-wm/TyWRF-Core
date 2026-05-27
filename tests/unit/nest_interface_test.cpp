#include "tywrf/nest/nest_interface.hpp"

#include <cassert>
#include <iostream>
#include <string_view>
#include <type_traits>

namespace {

using tywrf::nest::ExchangeOperation;
using tywrf::nest::HorizontalStagger;
using tywrf::nest::NestStatus;

void expect_status(
    const tywrf::nest::NestResult result,
    const NestStatus expected,
    const std::string_view label) {
  if (result.status != expected) {
    std::cerr << label << " returned "
              << tywrf::nest::nest_status_name(result.status) << ", expected "
              << tywrf::nest::nest_status_name(expected) << '\n';
    assert(false);
  }
}

void expect_exchange_status(
    const tywrf::nest::ExchangeResult result,
    const NestStatus expected,
    const ExchangeOperation operation,
    const std::string_view label) {
  if (result.status != expected || result.operation != operation) {
    std::cerr << label << " returned "
              << tywrf::nest::nest_status_name(result.status) << " for "
              << tywrf::nest::exchange_operation_name(result.operation) << '\n';
    assert(false);
  }
}

void expect_window(
    const tywrf::nest::RemapWindow& window,
    const HorizontalStagger stagger,
    const std::int32_t extent_i,
    const std::int32_t extent_j,
    const std::string_view label) {
  if (window.stagger != stagger || window.extent_i != extent_i ||
      window.extent_j != extent_j || window.new_i_begin != 0 ||
      window.new_j_begin != 0 || window.old_i_begin != 55 ||
      window.old_j_begin != 45 || window.old_i_offset_from_new != 55 ||
      window.old_j_offset_from_new != 45) {
    std::cerr << label << " remap window mismatch for "
              << tywrf::nest::horizontal_stagger_name(window.stagger) << '\n';
    assert(false);
  }
}

}  // namespace

int main() {
  static_assert(std::is_standard_layout_v<tywrf::nest::NestResult>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::NestResult>);
  static_assert(std::is_standard_layout_v<tywrf::nest::ExchangeResult>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::ExchangeResult>);
  static_assert(std::is_standard_layout_v<tywrf::nest::DomainPose>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::DomainPose>);
  static_assert(std::is_standard_layout_v<tywrf::nest::NestPose>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::NestPose>);
  static_assert(std::is_standard_layout_v<tywrf::nest::PoseDelta>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::PoseDelta>);
  static_assert(std::is_standard_layout_v<tywrf::nest::RemapWindow>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::RemapWindow>);
  static_assert(std::is_standard_layout_v<tywrf::nest::RemapPlan>);
  static_assert(std::is_trivially_copyable_v<tywrf::nest::RemapPlan>);

  const auto descriptor = tywrf::nest::make_krosa_parent_child_descriptor();
  expect_status(
      tywrf::nest::validate_parent_child_descriptor(descriptor), NestStatus::ok,
      "KROSA parent-child descriptor");

  assert(descriptor.parent.domain_id == 1);
  assert(descriptor.child.domain_id == 2);
  assert(descriptor.parent.grid_spacing_m == 10'000);
  assert(descriptor.child.grid_spacing_m == 2'000);
  assert(descriptor.parent_grid_ratio == 5);
  assert(descriptor.parent_time_step_ratio == 5);
  assert(descriptor.child.mass_nx == 210);
  assert(descriptor.child.mass_ny == 210);
  assert(descriptor.child.namelist_e_we == 211);
  assert(descriptor.child.namelist_e_sn == 211);

  auto reduced_resolution = descriptor;
  reduced_resolution.child.grid_spacing_m = 3'000;
  expect_status(
      tywrf::nest::validate_parent_child_descriptor(reduced_resolution),
      NestStatus::unsupported_resolution, "d02 resolution invariant");
  expect_status(
      tywrf::nest::make_domain_pose(
          reduced_resolution, {114, 96, tywrf::nest::IndexBase::one_based})
          .result,
      NestStatus::unsupported_resolution, "d02 pose resolution invariant");

  auto non_divisible_child = descriptor;
  non_divisible_child.child.mass_nx = 209;
  non_divisible_child.child.namelist_e_we = 210;
  expect_status(
      tywrf::nest::validate_parent_child_descriptor(non_divisible_child),
      NestStatus::invalid_configuration, "child ratio divisibility");

  const auto initial_position = tywrf::nest::make_krosa_initial_d02_position();
  assert(initial_position.i_parent_start == 114);
  assert(initial_position.j_parent_start == 96);

  const auto footprint = tywrf::nest::parent_child_footprint(descriptor, initial_position);
  assert(footprint.span_i_parent_cells == 42);
  assert(footprint.span_j_parent_cells == 42);
  assert(footprint.i_parent_end == 156);
  assert(footprint.j_parent_end == 138);

  const auto from_pose = tywrf::nest::make_domain_pose(descriptor, initial_position);
  expect_status(from_pose.result, NestStatus::ok, "initial d02 domain pose");
  assert(from_pose.parent.domain_id == 1);
  assert(from_pose.domain.domain_id == 2);
  assert(from_pose.domain.grid_spacing_m == 2'000);
  assert(from_pose.parent_grid_ratio == 5);
  assert(from_pose.parent_footprint.i_parent_start == 114);
  assert(from_pose.parent_footprint.j_parent_start == 96);

  const tywrf::nest::ParentChildPosition moved_position{
      125, 105, tywrf::nest::IndexBase::one_based};
  const auto to_pose = tywrf::nest::make_domain_pose(descriptor, moved_position);
  expect_status(to_pose.result, NestStatus::ok, "moved d02 domain pose");

  const auto delta = tywrf::nest::pose_delta(
      from_pose, to_pose, descriptor.parent_grid_ratio);
  expect_status(delta.result, NestStatus::ok, "KROSA pose delta");
  assert(delta.parent_di == 11);
  assert(delta.parent_dj == 9);
  assert(delta.child_di == 55);
  assert(delta.child_dj == 45);
  assert(delta.parent_to_child_ratio == 5);

  const auto from_nest_pose = tywrf::nest::make_nest_pose(descriptor, initial_position);
  const auto to_nest_pose = tywrf::nest::make_nest_pose(descriptor, moved_position);
  const auto nest_delta = tywrf::nest::pose_delta(from_nest_pose, to_nest_pose);
  expect_status(nest_delta.result, NestStatus::ok, "KROSA nest pose delta");
  assert(nest_delta.child_di == 55);
  assert(nest_delta.child_dj == 45);

  const auto remap_plan = tywrf::nest::build_remap_plan(from_pose, to_pose);
  expect_status(remap_plan.result, NestStatus::ok, "KROSA d02 remap plan");
  assert(remap_plan.delta.child_di == 55);
  assert(remap_plan.delta.child_dj == 45);
  expect_window(remap_plan.mass, HorizontalStagger::mass, 155, 165, "mass");
  expect_window(remap_plan.surface, HorizontalStagger::surface, 155, 165, "surface");
  expect_window(remap_plan.w_full, HorizontalStagger::w_full, 155, 165, "w_full");
  expect_window(remap_plan.u, HorizontalStagger::u, 156, 165, "U");
  expect_window(remap_plan.v, HorizontalStagger::v, 155, 166, "V");

  const auto nest_remap_plan =
      tywrf::nest::build_remap_plan(from_nest_pose, to_nest_pose);
  expect_status(nest_remap_plan.result, NestStatus::ok, "KROSA nest remap plan");
  assert(nest_remap_plan.mass.extent_i == 155);
  assert(nest_remap_plan.mass.extent_j == 165);

  const auto no_overlap_pose = tywrf::nest::make_domain_pose(
      descriptor, {156, 96, tywrf::nest::IndexBase::one_based});
  expect_status(no_overlap_pose.result, NestStatus::ok, "non-overlapping d02 pose");
  expect_status(
      tywrf::nest::build_remap_plan(from_pose, no_overlap_pose).result,
      NestStatus::out_of_bounds, "non-overlapping remap plan");

  auto invalid_footprint_pose = from_pose;
  invalid_footprint_pose.parent_footprint.span_i_parent_cells = 41;
  expect_status(
      tywrf::nest::build_remap_plan(invalid_footprint_pose, to_pose).result,
      NestStatus::invalid_contract, "invalid remap footprint");

  const auto moving_config = tywrf::nest::make_krosa_moving_nest_config();
  assert(moving_config.parent.vortex_interval_minutes == 15);
  assert(moving_config.child.vortex_interval_minutes == 15);
  assert(moving_config.parent.corral_dist_parent_cells == 10);
  assert(moving_config.child.corral_dist_parent_cells == 30);
  assert(moving_config.child.max_vortex_speed_mps == 30.0);
  assert(moving_config.child.track_level_pa == 70'000.0);
  assert(moving_config.child.time_to_move_seconds == 0);

  expect_status(
      tywrf::nest::validate_parent_child_position(
          descriptor, initial_position, moving_config.child.corral_dist_parent_cells),
      NestStatus::ok, "initial d02 corral position");
  expect_status(
      tywrf::nest::validate_parent_child_position(
          descriptor, {30, 96, tywrf::nest::IndexBase::one_based},
          moving_config.child.corral_dist_parent_cells),
      NestStatus::corral_violation, "west corral violation");
  expect_status(
      tywrf::nest::validate_parent_child_position(
          descriptor, {250, 96, tywrf::nest::IndexBase::one_based},
          moving_config.child.corral_dist_parent_cells),
      NestStatus::out_of_bounds, "parent bounds violation");

  expect_status(
      tywrf::nest::validate_movement_proposal(
          descriptor, moving_config,
          {initial_position, {116, 96, tywrf::nest::IndexBase::one_based}, 900}),
      NestStatus::ok, "valid 15 minute d02 move");
  expect_status(
      tywrf::nest::validate_movement_proposal(
          descriptor, moving_config,
          {initial_position, {117, 96, tywrf::nest::IndexBase::one_based}, 900}),
      NestStatus::movement_too_fast, "max_vortex_speed guard");
  expect_status(
      tywrf::nest::validate_movement_proposal(
          descriptor, moving_config,
          {initial_position, {116, 96, tywrf::nest::IndexBase::one_based}, 600}),
      NestStatus::invalid_contract, "vortex_interval guard");

  const auto spectral = tywrf::nest::make_krosa_spectral_nudging_config();
  expect_status(
      tywrf::nest::validate_spectral_nudging_config(spectral), NestStatus::ok,
      "KROSA spectral nudging config");
  assert(spectral.grid_fdda == 2);
  assert(spectral.input_template == "wrffdda_d<domain>");
  assert(spectral.input_interval_seconds == 21'600);
  assert(spectral.guv == 0.0003);
  assert(spectral.xwavenum == 2);
  assert(spectral.ywavenum == 4);
  assert(spectral.fields[0].old_variable == "U_NDG_OLD");
  assert(spectral.fields[0].new_variable == "U_NDG_NEW");
  assert(spectral.fields[3].field == tywrf::nest::NudgingField::q);
  assert(spectral.fields[3].old_variable == "Q_NDG_OLD");
  assert(spectral.fields[3].new_variable == "Q_NDG_NEW");
  assert(spectral.fields[5].field == tywrf::nest::NudgingField::mu);
  assert(spectral.fields[5].old_variable == "MU_NDG_OLD");
  assert(spectral.fields[5].new_variable == "MU_NDG_NEW");
  assert(!spectral.fields[5].three_dimensional);

  auto invalid_spectral = spectral;
  invalid_spectral.grid_fdda = 1;
  expect_status(
      tywrf::nest::validate_spectral_nudging_config(invalid_spectral),
      NestStatus::invalid_configuration, "unsupported grid_fdda");

  const tywrf::nest::ExchangeContract interpolation_contract{
      descriptor, initial_position, 7, 3, 280, 288};
  expect_exchange_status(
      tywrf::nest::interpolate_parent_to_child(interpolation_contract),
      NestStatus::not_implemented, ExchangeOperation::parent_to_child_interpolation,
      "parent-child interpolation stub");

  auto invalid_interpolation = interpolation_contract;
  invalid_interpolation.child_substep_index = descriptor.parent_time_step_ratio;
  expect_exchange_status(
      tywrf::nest::interpolate_parent_to_child(invalid_interpolation),
      NestStatus::invalid_contract, ExchangeOperation::parent_to_child_interpolation,
      "interpolation child substep guard");

  const tywrf::nest::ExchangeContract feedback_contract{
      descriptor, initial_position, 7, -1, 280, 320};
  expect_exchange_status(
      tywrf::nest::apply_child_feedback(feedback_contract), NestStatus::not_implemented,
      ExchangeOperation::child_to_parent_feedback, "child feedback stub");

  auto invalid_feedback = feedback_contract;
  invalid_feedback.child_substep_index = 0;
  expect_exchange_status(
      tywrf::nest::apply_child_feedback(invalid_feedback), NestStatus::invalid_contract,
      ExchangeOperation::child_to_parent_feedback, "feedback substep guard");

  std::cout << "Validated KROSA nest and spectral nudging interface baseline\n";
  return 0;
}
