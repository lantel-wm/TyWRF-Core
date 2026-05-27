#pragma once

#include <array>
#include <cstdint>
#include <string_view>

namespace tywrf::nest {

enum class NestStatus : std::uint8_t {
  ok,
  invalid_configuration,
  invalid_contract,
  unsupported_resolution,
  out_of_bounds,
  corral_violation,
  movement_too_fast,
  not_implemented,
};

enum class ExchangeOperation : std::uint8_t {
  parent_to_child_interpolation,
  child_to_parent_feedback,
};

enum class IndexBase : std::uint8_t {
  zero_based,
  one_based,
};

struct NestResult {
  NestStatus status = NestStatus::ok;
  const char* message = "";

  [[nodiscard]] constexpr bool ok() const noexcept {
    return status == NestStatus::ok;
  }
};

struct HorizontalDomainDescriptor {
  std::int32_t domain_id = 0;
  std::int32_t grid_spacing_m = 0;
  std::int32_t mass_nx = 0;
  std::int32_t mass_ny = 0;
  std::int32_t namelist_e_we = 0;
  std::int32_t namelist_e_sn = 0;
};

struct ParentChildDescriptor {
  HorizontalDomainDescriptor parent;
  HorizontalDomainDescriptor child;
  std::int32_t parent_grid_ratio = 0;
  std::int32_t parent_time_step_ratio = 0;
};

struct ParentChildPosition {
  std::int32_t i_parent_start = 0;
  std::int32_t j_parent_start = 0;
  IndexBase index_base = IndexBase::one_based;
};

struct ParentChildFootprint {
  std::int32_t i_parent_start = 0;
  std::int32_t j_parent_start = 0;
  std::int32_t i_parent_end = 0;
  std::int32_t j_parent_end = 0;
  std::int32_t span_i_parent_cells = 0;
  std::int32_t span_j_parent_cells = 0;
  IndexBase index_base = IndexBase::one_based;
};

struct DomainMovingNestConfig {
  std::int32_t domain_id = 0;
  std::int32_t vortex_interval_minutes = 0;
  double max_vortex_speed_mps = 0.0;
  std::int32_t corral_dist_parent_cells = 0;
  double track_level_pa = 0.0;
  std::int32_t time_to_move_seconds = 0;
};

struct MovingNestConfig {
  DomainMovingNestConfig parent;
  DomainMovingNestConfig child;
};

struct MovementProposal {
  ParentChildPosition from;
  ParentChildPosition to;
  std::int32_t elapsed_seconds = 0;
};

struct ExchangeContract {
  ParentChildDescriptor relationship;
  ParentChildPosition position;
  std::int64_t parent_step_index = 0;
  std::int32_t child_substep_index = -1;
  std::int64_t start_seconds = 0;
  std::int64_t end_seconds = 0;
};

struct ExchangeResult {
  NestStatus status = NestStatus::ok;
  ExchangeOperation operation = ExchangeOperation::parent_to_child_interpolation;
  const char* message = "";

  [[nodiscard]] constexpr bool ok() const noexcept {
    return status == NestStatus::ok;
  }
};

enum class NudgingField : std::uint8_t {
  u,
  v,
  t,
  q,
  ph,
  mu,
};

struct NudgingFieldBinding {
  NudgingField field = NudgingField::u;
  std::string_view old_variable;
  std::string_view new_variable;
  bool three_dimensional = true;
};

struct SpectralNudgingConfig {
  std::int32_t grid_fdda = 0;
  std::string_view input_template;
  std::int32_t input_interval_seconds = 0;
  std::int32_t end_seconds = 0;
  double guv = 0.0;
  std::int32_t xwavenum = 0;
  std::int32_t ywavenum = 0;
  std::array<NudgingFieldBinding, 6> fields{};
};

[[nodiscard]] ParentChildDescriptor make_krosa_parent_child_descriptor() noexcept;
[[nodiscard]] ParentChildPosition make_krosa_initial_d02_position() noexcept;
[[nodiscard]] MovingNestConfig make_krosa_moving_nest_config() noexcept;
[[nodiscard]] SpectralNudgingConfig make_krosa_spectral_nudging_config() noexcept;

[[nodiscard]] NestResult validate_parent_child_descriptor(
    const ParentChildDescriptor& descriptor) noexcept;

[[nodiscard]] ParentChildFootprint parent_child_footprint(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& position) noexcept;

[[nodiscard]] NestResult validate_parent_child_position(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& position,
    std::int32_t corral_dist_parent_cells = 0) noexcept;

[[nodiscard]] NestResult validate_movement_proposal(
    const ParentChildDescriptor& descriptor,
    const MovingNestConfig& moving_config,
    const MovementProposal& proposal) noexcept;

[[nodiscard]] NestResult validate_spectral_nudging_config(
    const SpectralNudgingConfig& config) noexcept;

[[nodiscard]] ExchangeResult interpolate_parent_to_child(
    const ExchangeContract& contract) noexcept;

[[nodiscard]] ExchangeResult apply_child_feedback(const ExchangeContract& contract) noexcept;

[[nodiscard]] std::string_view nest_status_name(NestStatus status) noexcept;
[[nodiscard]] std::string_view exchange_operation_name(ExchangeOperation operation) noexcept;
[[nodiscard]] std::string_view nudging_field_name(NudgingField field) noexcept;

}  // namespace tywrf::nest
