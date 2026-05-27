#include "tywrf/nest/nest_interface.hpp"

#include <cmath>

namespace tywrf::nest {
namespace {

[[nodiscard]] constexpr NestResult ok_result() noexcept {
  return {NestStatus::ok, "ok"};
}

[[nodiscard]] constexpr NestResult result(
    const NestStatus status,
    const char* message) noexcept {
  return {status, message};
}

[[nodiscard]] constexpr ExchangeResult exchange_result(
    const NestStatus status,
    const ExchangeOperation operation,
    const char* message) noexcept {
  return {status, operation, message};
}

[[nodiscard]] constexpr bool positive_domain(
    const HorizontalDomainDescriptor& domain) noexcept {
  return domain.domain_id > 0 && domain.grid_spacing_m > 0 && domain.mass_nx > 0 &&
         domain.mass_ny > 0 && domain.namelist_e_we > 0 && domain.namelist_e_sn > 0;
}

[[nodiscard]] constexpr bool same_horizontal_domain(
    const HorizontalDomainDescriptor& lhs,
    const HorizontalDomainDescriptor& rhs) noexcept {
  return lhs.domain_id == rhs.domain_id && lhs.grid_spacing_m == rhs.grid_spacing_m &&
         lhs.mass_nx == rhs.mass_nx && lhs.mass_ny == rhs.mass_ny &&
         lhs.namelist_e_we == rhs.namelist_e_we &&
         lhs.namelist_e_sn == rhs.namelist_e_sn;
}

[[nodiscard]] constexpr bool same_parent_child_descriptor(
    const ParentChildDescriptor& lhs,
    const ParentChildDescriptor& rhs) noexcept {
  return same_horizontal_domain(lhs.parent, rhs.parent) &&
         same_horizontal_domain(lhs.child, rhs.child) &&
         lhs.parent_grid_ratio == rhs.parent_grid_ratio &&
         lhs.parent_time_step_ratio == rhs.parent_time_step_ratio;
}

[[nodiscard]] constexpr std::int32_t lower_bound(
    const IndexBase index_base,
    const std::int32_t corral_dist_parent_cells) noexcept {
  return index_base == IndexBase::one_based ? 1 + corral_dist_parent_cells
                                            : corral_dist_parent_cells;
}

[[nodiscard]] constexpr std::int32_t upper_bound(
    const std::int32_t namelist_extent,
    const IndexBase index_base,
    const std::int32_t corral_dist_parent_cells) noexcept {
  return index_base == IndexBase::one_based ? namelist_extent - corral_dist_parent_cells
                                            : namelist_extent - 1 - corral_dist_parent_cells;
}

[[nodiscard]] NestResult validate_domain_pose(const DomainPose& pose) noexcept {
  if (!pose.result.ok()) {
    return pose.result;
  }
  if (!positive_domain(pose.parent) || !positive_domain(pose.domain)) {
    return result(NestStatus::invalid_configuration, "domain pose descriptor must be positive");
  }
  if (pose.parent_grid_ratio <= 0) {
    return result(NestStatus::invalid_configuration, "domain pose ratio must be positive");
  }
  if (pose.parent.grid_spacing_m != pose.domain.grid_spacing_m * pose.parent_grid_ratio) {
    return result(
        NestStatus::unsupported_resolution,
        "domain pose spacing must match parent_grid_ratio");
  }
  if (pose.domain.grid_spacing_m != 2'000) {
    return result(NestStatus::unsupported_resolution, "TyWRF-Core v1 requires d02 at 2 km");
  }
  if (pose.domain.mass_nx % pose.parent_grid_ratio != 0 ||
      pose.domain.mass_ny % pose.parent_grid_ratio != 0) {
    return result(
        NestStatus::invalid_configuration,
        "domain pose mass dimensions must be divisible by parent_grid_ratio");
  }
  if (pose.parent_position.index_base != pose.parent_footprint.index_base) {
    return result(NestStatus::invalid_contract, "domain pose index bases must match");
  }

  const auto expected_span_i = pose.domain.mass_nx / pose.parent_grid_ratio;
  const auto expected_span_j = pose.domain.mass_ny / pose.parent_grid_ratio;
  if (pose.parent_footprint.span_i_parent_cells != expected_span_i ||
      pose.parent_footprint.span_j_parent_cells != expected_span_j ||
      pose.parent_footprint.i_parent_start != pose.parent_position.i_parent_start ||
      pose.parent_footprint.j_parent_start != pose.parent_position.j_parent_start ||
      pose.parent_footprint.i_parent_end !=
          pose.parent_footprint.i_parent_start + expected_span_i ||
      pose.parent_footprint.j_parent_end !=
          pose.parent_footprint.j_parent_start + expected_span_j) {
    return result(NestStatus::invalid_contract, "domain pose footprint is inconsistent");
  }

  const auto i_lower = lower_bound(pose.parent_position.index_base, 0);
  const auto j_lower = lower_bound(pose.parent_position.index_base, 0);
  const auto i_upper = upper_bound(pose.parent.namelist_e_we, pose.parent_position.index_base, 0);
  const auto j_upper = upper_bound(pose.parent.namelist_e_sn, pose.parent_position.index_base, 0);
  if (pose.parent_footprint.i_parent_start < i_lower ||
      pose.parent_footprint.j_parent_start < j_lower ||
      pose.parent_footprint.i_parent_end > i_upper ||
      pose.parent_footprint.j_parent_end > j_upper) {
    return result(NestStatus::out_of_bounds, "domain pose footprint is outside parent domain");
  }

  return ok_result();
}

[[nodiscard]] NestResult validate_remap_pair(
    const DomainPose& from,
    const DomainPose& to,
    const std::int32_t parent_to_child_ratio) noexcept {
  if (const auto from_status = validate_domain_pose(from); !from_status.ok()) {
    return from_status;
  }
  if (const auto to_status = validate_domain_pose(to); !to_status.ok()) {
    return to_status;
  }
  if (!same_horizontal_domain(from.parent, to.parent) ||
      !same_horizontal_domain(from.domain, to.domain) ||
      from.parent_grid_ratio != to.parent_grid_ratio) {
    return result(NestStatus::invalid_contract, "remap poses must describe the same child grid");
  }
  if (from.parent_position.index_base != to.parent_position.index_base) {
    return result(NestStatus::invalid_contract, "remap poses must use the same index base");
  }
  if (parent_to_child_ratio <= 0) {
    return result(NestStatus::invalid_configuration, "pose delta ratio must be positive");
  }
  if (parent_to_child_ratio != from.parent_grid_ratio) {
    return result(NestStatus::invalid_contract, "pose delta ratio must match parent_grid_ratio");
  }
  return ok_result();
}

[[nodiscard]] RemapWindow make_window(
    const HorizontalStagger stagger,
    const std::int32_t mass_nx,
    const std::int32_t mass_ny,
    const PoseDelta& delta) noexcept {
  auto nx = mass_nx;
  auto ny = mass_ny;
  if (stagger == HorizontalStagger::u) {
    ++nx;
  } else if (stagger == HorizontalStagger::v) {
    ++ny;
  }

  RemapWindow window{};
  window.stagger = stagger;
  if (delta.child_di >= 0) {
    window.old_i_begin = delta.child_di;
    window.new_i_begin = 0;
    window.extent_i = nx - delta.child_di;
  } else {
    window.old_i_begin = 0;
    window.new_i_begin = -delta.child_di;
    window.extent_i = nx + delta.child_di;
  }

  if (delta.child_dj >= 0) {
    window.old_j_begin = delta.child_dj;
    window.new_j_begin = 0;
    window.extent_j = ny - delta.child_dj;
  } else {
    window.old_j_begin = 0;
    window.new_j_begin = -delta.child_dj;
    window.extent_j = ny + delta.child_dj;
  }

  window.old_i_offset_from_new = window.old_i_begin - window.new_i_begin;
  window.old_j_offset_from_new = window.old_j_begin - window.new_j_begin;
  return window;
}

[[nodiscard]] NestResult validate_exchange_common(
    const ExchangeContract& contract,
    const bool require_child_substep) noexcept {
  if (const auto descriptor_status =
          validate_parent_child_descriptor(contract.relationship);
      !descriptor_status.ok()) {
    return descriptor_status;
  }

  if (const auto position_status =
          validate_parent_child_position(contract.relationship, contract.position);
      !position_status.ok()) {
    return position_status;
  }

  if (contract.parent_step_index < 0) {
    return result(NestStatus::invalid_contract, "parent step index must be non-negative");
  }
  if (contract.end_seconds <= contract.start_seconds) {
    return result(NestStatus::invalid_contract, "exchange end time must be after start time");
  }

  if (require_child_substep) {
    if (contract.child_substep_index < 0 ||
        contract.child_substep_index >= contract.relationship.parent_time_step_ratio) {
      return result(NestStatus::invalid_contract, "child substep index is outside ratio");
    }
  } else if (contract.child_substep_index != -1) {
    return result(NestStatus::invalid_contract, "feedback contract must not name a child substep");
  }

  return ok_result();
}

}  // namespace

ParentChildDescriptor make_krosa_parent_child_descriptor() noexcept {
  return {
      HorizontalDomainDescriptor{1, 10'000, 265, 429, 266, 430},
      HorizontalDomainDescriptor{2, 2'000, 210, 210, 211, 211},
      5,
      5,
  };
}

ParentChildPosition make_krosa_initial_d02_position() noexcept {
  return {114, 96, IndexBase::one_based};
}

MovingNestConfig make_krosa_moving_nest_config() noexcept {
  return {
      DomainMovingNestConfig{1, 15, 30.0, 10, 70'000.0, 0},
      DomainMovingNestConfig{2, 15, 30.0, 30, 70'000.0, 0},
  };
}

SpectralNudgingConfig make_krosa_spectral_nudging_config() noexcept {
  return {
      2,
      "wrffdda_d<domain>",
      21'600,
      336 * 3'600,
      0.0003,
      2,
      4,
      std::array<NudgingFieldBinding, 6>{
          NudgingFieldBinding{NudgingField::u, "U_NDG_OLD", "U_NDG_NEW", true},
          NudgingFieldBinding{NudgingField::v, "V_NDG_OLD", "V_NDG_NEW", true},
          NudgingFieldBinding{NudgingField::t, "T_NDG_OLD", "T_NDG_NEW", true},
          NudgingFieldBinding{NudgingField::q, "Q_NDG_OLD", "Q_NDG_NEW", true},
          NudgingFieldBinding{NudgingField::ph, "PH_NDG_OLD", "PH_NDG_NEW", true},
          NudgingFieldBinding{NudgingField::mu, "MU_NDG_OLD", "MU_NDG_NEW", false},
      },
  };
}

DomainPose make_domain_pose(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& position) noexcept {
  DomainPose pose{};
  pose.parent = descriptor.parent;
  pose.domain = descriptor.child;
  pose.parent_position = position;
  pose.parent_grid_ratio = descriptor.parent_grid_ratio;

  if (const auto descriptor_status = validate_parent_child_descriptor(descriptor);
      !descriptor_status.ok()) {
    pose.result = descriptor_status;
    return pose;
  }

  pose.parent_footprint = parent_child_footprint(descriptor, position);
  if (const auto position_status = validate_parent_child_position(descriptor, position);
      !position_status.ok()) {
    pose.result = position_status;
    return pose;
  }

  pose.result = ok_result();
  return pose;
}

NestPose make_nest_pose(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& position) noexcept {
  return {descriptor, make_domain_pose(descriptor, position)};
}

NestResult validate_parent_child_descriptor(
    const ParentChildDescriptor& descriptor) noexcept {
  if (!positive_domain(descriptor.parent) || !positive_domain(descriptor.child)) {
    return result(NestStatus::invalid_configuration, "domain descriptors must be positive");
  }
  if (descriptor.parent.domain_id == descriptor.child.domain_id) {
    return result(NestStatus::invalid_configuration, "parent and child domain ids must differ");
  }
  if (descriptor.parent_grid_ratio <= 0 || descriptor.parent_time_step_ratio <= 0) {
    return result(NestStatus::invalid_configuration, "nest ratios must be positive");
  }
  if (descriptor.parent.grid_spacing_m !=
      descriptor.child.grid_spacing_m * descriptor.parent_grid_ratio) {
    return result(
        NestStatus::unsupported_resolution,
        "parent and child grid spacing must match parent_grid_ratio");
  }
  if (descriptor.child.grid_spacing_m != 2'000) {
    return result(NestStatus::unsupported_resolution, "TyWRF-Core v1 requires d02 at 2 km");
  }
  if (descriptor.child.mass_nx % descriptor.parent_grid_ratio != 0 ||
      descriptor.child.mass_ny % descriptor.parent_grid_ratio != 0) {
    return result(
        NestStatus::invalid_configuration,
        "child mass dimensions must be divisible by parent_grid_ratio");
  }
  if (descriptor.parent.namelist_e_we != descriptor.parent.mass_nx + 1 ||
      descriptor.parent.namelist_e_sn != descriptor.parent.mass_ny + 1 ||
      descriptor.child.namelist_e_we != descriptor.child.mass_nx + 1 ||
      descriptor.child.namelist_e_sn != descriptor.child.mass_ny + 1) {
    return result(
        NestStatus::invalid_configuration,
        "namelist horizontal extents must equal mass dimensions plus one");
  }
  return ok_result();
}

ParentChildFootprint parent_child_footprint(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& position) noexcept {
  const auto span_i = descriptor.child.mass_nx / descriptor.parent_grid_ratio;
  const auto span_j = descriptor.child.mass_ny / descriptor.parent_grid_ratio;
  return {
      position.i_parent_start,
      position.j_parent_start,
      position.i_parent_start + span_i,
      position.j_parent_start + span_j,
      span_i,
      span_j,
      position.index_base,
  };
}

NestResult validate_parent_child_position(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& position,
    const std::int32_t corral_dist_parent_cells) noexcept {
  if (const auto descriptor_status = validate_parent_child_descriptor(descriptor);
      !descriptor_status.ok()) {
    return descriptor_status;
  }
  if (corral_dist_parent_cells < 0) {
    return result(NestStatus::invalid_configuration, "corral distance must be non-negative");
  }

  const auto footprint = parent_child_footprint(descriptor, position);
  const auto i_lower = lower_bound(position.index_base, 0);
  const auto j_lower = lower_bound(position.index_base, 0);
  const auto i_upper = upper_bound(descriptor.parent.namelist_e_we, position.index_base, 0);
  const auto j_upper = upper_bound(descriptor.parent.namelist_e_sn, position.index_base, 0);

  if (footprint.i_parent_start < i_lower || footprint.j_parent_start < j_lower ||
      footprint.i_parent_end > i_upper || footprint.j_parent_end > j_upper) {
    return result(NestStatus::out_of_bounds, "child footprint is outside parent domain");
  }

  if (corral_dist_parent_cells == 0) {
    return ok_result();
  }

  const auto i_corral_lower = lower_bound(position.index_base, corral_dist_parent_cells);
  const auto j_corral_lower = lower_bound(position.index_base, corral_dist_parent_cells);
  const auto i_corral_upper =
      upper_bound(descriptor.parent.namelist_e_we, position.index_base, corral_dist_parent_cells);
  const auto j_corral_upper =
      upper_bound(descriptor.parent.namelist_e_sn, position.index_base, corral_dist_parent_cells);

  if (footprint.i_parent_start < i_corral_lower ||
      footprint.j_parent_start < j_corral_lower ||
      footprint.i_parent_end > i_corral_upper ||
      footprint.j_parent_end > j_corral_upper) {
    return result(NestStatus::corral_violation, "child footprint violates corral distance");
  }

  return ok_result();
}

NestResult validate_movement_proposal(
    const ParentChildDescriptor& descriptor,
    const MovingNestConfig& moving_config,
    const MovementProposal& proposal) noexcept {
  if (proposal.from.index_base != proposal.to.index_base) {
    return result(NestStatus::invalid_contract, "movement endpoints must use the same index base");
  }
  if (moving_config.parent.domain_id != descriptor.parent.domain_id ||
      moving_config.child.domain_id != descriptor.child.domain_id) {
    return result(NestStatus::invalid_configuration, "moving config domains do not match nest");
  }
  if (moving_config.child.vortex_interval_minutes <= 0 ||
      moving_config.child.max_vortex_speed_mps <= 0.0 ||
      moving_config.child.track_level_pa <= 0.0 ||
      moving_config.child.time_to_move_seconds < 0) {
    return result(NestStatus::invalid_configuration, "moving nest config contains invalid values");
  }
  if (proposal.elapsed_seconds <= 0) {
    return result(NestStatus::invalid_contract, "movement elapsed time must be positive");
  }
  if (proposal.elapsed_seconds != moving_config.child.vortex_interval_minutes * 60) {
    return result(
        NestStatus::invalid_contract, "movement elapsed time must match vortex_interval");
  }

  if (const auto from_status = validate_parent_child_position(
          descriptor, proposal.from, moving_config.child.corral_dist_parent_cells);
      !from_status.ok()) {
    return from_status;
  }
  if (const auto to_status = validate_parent_child_position(
          descriptor, proposal.to, moving_config.child.corral_dist_parent_cells);
      !to_status.ok()) {
    return to_status;
  }

  const auto di =
      static_cast<double>(proposal.to.i_parent_start - proposal.from.i_parent_start);
  const auto dj =
      static_cast<double>(proposal.to.j_parent_start - proposal.from.j_parent_start);
  const auto distance_cells = std::sqrt(di * di + dj * dj);
  const auto speed_mps =
      distance_cells * static_cast<double>(descriptor.parent.grid_spacing_m) /
      static_cast<double>(proposal.elapsed_seconds);

  if (speed_mps > moving_config.child.max_vortex_speed_mps) {
    return result(NestStatus::movement_too_fast, "moving nest proposal exceeds max_vortex_speed");
  }

  return ok_result();
}

NestResult validate_spectral_nudging_config(
    const SpectralNudgingConfig& config) noexcept {
  if (config.grid_fdda != 2) {
    return result(NestStatus::invalid_configuration, "KROSA v1 expects grid_fdda = 2");
  }
  if (config.input_template.empty() || config.input_interval_seconds <= 0 ||
      config.end_seconds < 0) {
    return result(NestStatus::invalid_configuration, "spectral nudging input timing is invalid");
  }
  if (config.guv <= 0.0 || config.xwavenum <= 0 || config.ywavenum <= 0) {
    return result(NestStatus::invalid_configuration, "spectral nudging constants are invalid");
  }
  for (const auto& field : config.fields) {
    if (field.old_variable.empty() || field.new_variable.empty()) {
      return result(NestStatus::invalid_configuration, "nudging field binding is missing names");
    }
  }
  return ok_result();
}

PoseDelta pose_delta(
    const DomainPose& from,
    const DomainPose& to,
    const std::int32_t parent_to_child_ratio) noexcept {
  PoseDelta delta{};
  delta.parent_to_child_ratio = parent_to_child_ratio;

  if (const auto pair_status = validate_remap_pair(from, to, parent_to_child_ratio);
      !pair_status.ok()) {
    delta.result = pair_status;
    return delta;
  }

  delta.parent_di = to.parent_position.i_parent_start - from.parent_position.i_parent_start;
  delta.parent_dj = to.parent_position.j_parent_start - from.parent_position.j_parent_start;
  delta.child_di = delta.parent_di * parent_to_child_ratio;
  delta.child_dj = delta.parent_dj * parent_to_child_ratio;
  delta.result = ok_result();
  return delta;
}

PoseDelta pose_delta(const NestPose& from, const NestPose& to) noexcept {
  if (!same_parent_child_descriptor(from.relationship, to.relationship)) {
    PoseDelta delta{};
    delta.result =
        result(NestStatus::invalid_contract, "nest poses must use the same relationship");
    return delta;
  }
  return pose_delta(from.child, to.child, from.relationship.parent_grid_ratio);
}

RemapPlan build_remap_plan(const DomainPose& from, const DomainPose& to) noexcept {
  RemapPlan plan{};
  plan.delta = pose_delta(from, to, from.parent_grid_ratio);
  if (!plan.delta.ok()) {
    plan.result = plan.delta.result;
    return plan;
  }

  plan.mass = make_window(
      HorizontalStagger::mass, from.domain.mass_nx, from.domain.mass_ny, plan.delta);
  plan.u = make_window(
      HorizontalStagger::u, from.domain.mass_nx, from.domain.mass_ny, plan.delta);
  plan.v = make_window(
      HorizontalStagger::v, from.domain.mass_nx, from.domain.mass_ny, plan.delta);
  plan.w_full = make_window(
      HorizontalStagger::w_full, from.domain.mass_nx, from.domain.mass_ny, plan.delta);
  plan.surface = make_window(
      HorizontalStagger::surface, from.domain.mass_nx, from.domain.mass_ny, plan.delta);

  if (plan.mass.empty() || plan.u.empty() || plan.v.empty() || plan.w_full.empty() ||
      plan.surface.empty()) {
    plan.result = result(NestStatus::out_of_bounds, "remap poses do not overlap");
    return plan;
  }

  plan.result = ok_result();
  return plan;
}

RemapPlan build_remap_plan(const NestPose& from, const NestPose& to) noexcept {
  if (!same_parent_child_descriptor(from.relationship, to.relationship)) {
    RemapPlan plan{};
    plan.result =
        result(NestStatus::invalid_contract, "nest poses must use the same relationship");
    return plan;
  }
  return build_remap_plan(from.child, to.child);
}

ExchangeResult interpolate_parent_to_child(const ExchangeContract& contract) noexcept {
  if (const auto common = validate_exchange_common(contract, true); !common.ok()) {
    return exchange_result(
        common.status, ExchangeOperation::parent_to_child_interpolation, common.message);
  }
  return exchange_result(
      NestStatus::not_implemented,
      ExchangeOperation::parent_to_child_interpolation,
      "parent-child numerical interpolation is not implemented yet");
}

ExchangeResult apply_child_feedback(const ExchangeContract& contract) noexcept {
  if (const auto common = validate_exchange_common(contract, false); !common.ok()) {
    return exchange_result(
        common.status, ExchangeOperation::child_to_parent_feedback, common.message);
  }
  return exchange_result(
      NestStatus::not_implemented,
      ExchangeOperation::child_to_parent_feedback,
      "two-way nesting feedback is not implemented yet");
}

std::string_view nest_status_name(const NestStatus status) noexcept {
  switch (status) {
    case NestStatus::ok:
      return "ok";
    case NestStatus::invalid_configuration:
      return "invalid_configuration";
    case NestStatus::invalid_contract:
      return "invalid_contract";
    case NestStatus::unsupported_resolution:
      return "unsupported_resolution";
    case NestStatus::out_of_bounds:
      return "out_of_bounds";
    case NestStatus::corral_violation:
      return "corral_violation";
    case NestStatus::movement_too_fast:
      return "movement_too_fast";
    case NestStatus::not_implemented:
      return "not_implemented";
  }
  return "unknown";
}

std::string_view exchange_operation_name(const ExchangeOperation operation) noexcept {
  switch (operation) {
    case ExchangeOperation::parent_to_child_interpolation:
      return "parent_to_child_interpolation";
    case ExchangeOperation::child_to_parent_feedback:
      return "child_to_parent_feedback";
  }
  return "unknown";
}

std::string_view nudging_field_name(const NudgingField field) noexcept {
  switch (field) {
    case NudgingField::u:
      return "U";
    case NudgingField::v:
      return "V";
    case NudgingField::t:
      return "T";
    case NudgingField::q:
      return "Q";
    case NudgingField::ph:
      return "PH";
    case NudgingField::mu:
      return "MU";
  }
  return "unknown";
}

std::string_view horizontal_stagger_name(const HorizontalStagger stagger) noexcept {
  switch (stagger) {
    case HorizontalStagger::mass:
      return "mass";
    case HorizontalStagger::u:
      return "U";
    case HorizontalStagger::v:
      return "V";
    case HorizontalStagger::w_full:
      return "W_full";
    case HorizontalStagger::surface:
      return "surface";
  }
  return "unknown";
}

}  // namespace tywrf::nest
