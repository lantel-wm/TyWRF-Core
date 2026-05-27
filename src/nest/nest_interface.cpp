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

}  // namespace tywrf::nest
