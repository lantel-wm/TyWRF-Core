#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/io/forcing_frames.hpp"
#include "tywrf/state.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace tywrf::io {

enum class BoundarySide {
  i_lower,
  i_upper,
  j_lower,
  j_upper,
};

enum class ForcingApplyStatus {
  ok,
  shape_mismatch,
  invalid_range,
  unsupported_field,
};

enum class ForcingApplyOperation {
  validation_only,
  synthetic_nudging_delta,
  direct_boundary_copy_skeleton,
};

class ForcingApplyError : public std::runtime_error {
 public:
  ForcingApplyError(ForcingApplyStatus status, const std::string& message)
      : std::runtime_error(message), status_(status) {}

  [[nodiscard]] ForcingApplyStatus status() const noexcept {
    return status_;
  }

 private:
  ForcingApplyStatus status_ = ForcingApplyStatus::ok;
};

struct ForcingApplyReport {
  std::size_t field_count = 0;
  std::size_t side_count = 0;
  std::size_t point_count = 0;
  ForcingApplyStatus status = ForcingApplyStatus::ok;
  ForcingApplyOperation operation = ForcingApplyOperation::validation_only;
  bool would_modify_state = false;
  bool synthetic = false;
  std::string detail;
};

struct SyntheticNudgingDeltaConfig {
  std::int32_t active_i_begin = 0;
  std::int32_t active_j_begin = 0;
  std::int32_t active_k_begin = 0;
};

[[nodiscard]] std::string_view forcing_apply_status_name(
    ForcingApplyStatus status) noexcept;

[[nodiscard]] std::string_view forcing_apply_operation_name(
    ForcingApplyOperation operation) noexcept;

[[nodiscard]] ForcingApplyReport validate_boundary_pack_for_view(
    FieldView3D<const float> target,
    BoundarySide side,
    const PackedForcingField& packed);

[[nodiscard]] ForcingApplyReport validate_boundary_pack_for_view(
    FieldView2D<const float> target,
    BoundarySide side,
    const PackedForcingField& packed);

[[nodiscard]] ForcingApplyReport validate_boundary_pack_for_state(
    const State<float>& state,
    std::string_view field_name,
    BoundarySide side,
    const PackedForcingField& packed);

ForcingApplyReport apply_boundary_copy_skeleton_to_view(
    FieldView3D<float> target,
    BoundarySide side,
    const PackedForcingField& packed);

ForcingApplyReport apply_boundary_copy_skeleton_to_view(
    FieldView2D<float> target,
    BoundarySide side,
    const PackedForcingField& packed);

ForcingApplyReport apply_boundary_copy_skeleton_to_state(
    State<float>& state,
    std::string_view field_name,
    BoundarySide side,
    const PackedForcingField& packed);

ForcingApplyReport apply_krosa_boundary_copy_skeleton_to_state(
    State<float>& state,
    const KrosaForcingReader& reader,
    std::string_view field_name,
    std::size_t record_index);

ForcingApplyReport apply_synthetic_nudging_delta(
    State<float>& state,
    std::string_view field_name,
    const PackedForcingField& synthetic_delta,
    SyntheticNudgingDeltaConfig config = {});

}  // namespace tywrf::io
