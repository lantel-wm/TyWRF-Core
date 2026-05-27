#pragma once

#include <array>
#include <cstddef>
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

enum class HorizontalStagger : std::uint8_t {
  mass,
  u,
  v,
  w_full,
  surface,
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

struct MovingNestPoseEvent {
  std::int64_t model_seconds = 0;
  std::int64_t parent_step_index = 0;
  ParentChildPosition position;
  bool moved = false;
};

struct MovingNestPoseTimeline {
  const MovingNestPoseEvent* events = nullptr;
  std::size_t event_count = 0;
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

struct DomainPose {
  NestResult result{NestStatus::ok, "ok"};
  HorizontalDomainDescriptor parent;
  HorizontalDomainDescriptor domain;
  ParentChildPosition parent_position;
  ParentChildFootprint parent_footprint;
  std::int32_t parent_grid_ratio = 0;
};

struct NestPose {
  ParentChildDescriptor relationship;
  DomainPose child;
};

struct PoseDelta {
  NestResult result{NestStatus::ok, "ok"};
  std::int32_t parent_di = 0;
  std::int32_t parent_dj = 0;
  std::int32_t child_di = 0;
  std::int32_t child_dj = 0;
  std::int32_t parent_to_child_ratio = 0;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

struct RemapWindow {
  HorizontalStagger stagger = HorizontalStagger::mass;
  std::int32_t old_i_begin = 0;
  std::int32_t old_j_begin = 0;
  std::int32_t new_i_begin = 0;
  std::int32_t new_j_begin = 0;
  std::int32_t extent_i = 0;
  std::int32_t extent_j = 0;
  std::int32_t old_i_offset_from_new = 0;
  std::int32_t old_j_offset_from_new = 0;

  [[nodiscard]] constexpr bool empty() const noexcept {
    return extent_i <= 0 || extent_j <= 0;
  }
};

struct RemapPlan {
  NestResult result{NestStatus::ok, "ok"};
  PoseDelta delta;
  RemapWindow mass;
  RemapWindow u;
  RemapWindow v;
  RemapWindow w_full;
  RemapWindow surface;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
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
[[nodiscard]] MovingNestPoseTimeline make_krosa_first_10min_d02_pose_timeline() noexcept;
[[nodiscard]] SpectralNudgingConfig make_krosa_spectral_nudging_config() noexcept;
[[nodiscard]] DomainPose make_domain_pose(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& position) noexcept;
[[nodiscard]] NestPose make_nest_pose(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& position) noexcept;

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

[[nodiscard]] NestResult validate_moving_nest_pose_timeline(
    const ParentChildDescriptor& descriptor,
    const MovingNestPoseTimeline& timeline) noexcept;

[[nodiscard]] const MovingNestPoseEvent* moving_nest_pose_at_model_second(
    const MovingNestPoseTimeline& timeline,
    std::int64_t model_seconds) noexcept;
[[nodiscard]] const MovingNestPoseEvent* moving_nest_final_pose(
    const MovingNestPoseTimeline& timeline) noexcept;
[[nodiscard]] std::size_t moving_nest_movement_count(
    const MovingNestPoseTimeline& timeline) noexcept;

[[nodiscard]] NestResult validate_spectral_nudging_config(
    const SpectralNudgingConfig& config) noexcept;

[[nodiscard]] PoseDelta pose_delta(
    const DomainPose& from,
    const DomainPose& to,
    std::int32_t parent_to_child_ratio) noexcept;
[[nodiscard]] PoseDelta pose_delta(
    const NestPose& from,
    const NestPose& to) noexcept;

[[nodiscard]] RemapPlan build_remap_plan(
    const DomainPose& from,
    const DomainPose& to) noexcept;
[[nodiscard]] RemapPlan build_remap_plan(
    const NestPose& from,
    const NestPose& to) noexcept;

[[nodiscard]] ExchangeResult interpolate_parent_to_child(
    const ExchangeContract& contract) noexcept;

[[nodiscard]] ExchangeResult apply_child_feedback(const ExchangeContract& contract) noexcept;

[[nodiscard]] std::string_view nest_status_name(NestStatus status) noexcept;
[[nodiscard]] std::string_view exchange_operation_name(ExchangeOperation operation) noexcept;
[[nodiscard]] std::string_view nudging_field_name(NudgingField field) noexcept;
[[nodiscard]] std::string_view horizontal_stagger_name(HorizontalStagger stagger) noexcept;

}  // namespace tywrf::nest
