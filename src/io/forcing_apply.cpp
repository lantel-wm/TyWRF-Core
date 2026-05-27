#include "tywrf/io/forcing_apply.hpp"

#include <array>
#include <limits>
#include <string>
#include <string_view>
#include <utility>

namespace tywrf::io {
namespace {

struct KrosaBoundarySideMapping {
  BoundarySide side = BoundarySide::i_lower;
  std::string_view suffix;
};

struct ActiveTargetShape {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
};

struct KrosaPackedBoundarySide {
  BoundarySide side = BoundarySide::i_lower;
  std::string variable_name;
  PackedForcingField packed;
};

constexpr std::array<KrosaBoundarySideMapping, 4> kKrosaBoundarySideMappings{{
    {BoundarySide::i_lower, "_BXS"},
    {BoundarySide::i_upper, "_BXE"},
    {BoundarySide::j_lower, "_BYS"},
    {BoundarySide::j_upper, "_BYE"},
}};

[[noreturn]] void fail(
    const ForcingApplyStatus status,
    const std::string& message) {
  throw ForcingApplyError(status, message);
}

void require(const bool condition, const ForcingApplyStatus status, const char* message) {
  if (!condition) {
    fail(status, message);
  }
}

[[nodiscard]] std::size_t checked_subtract(
    const std::int32_t total,
    const std::int32_t lower,
    const std::int32_t upper,
    const char* message) {
  require(total >= 0 && lower >= 0 && upper >= 0,
          ForcingApplyStatus::shape_mismatch,
          message);
  require(total >= lower + upper, ForcingApplyStatus::shape_mismatch, message);
  return static_cast<std::size_t>(total - lower - upper);
}

[[nodiscard]] ActiveTargetShape active_shape(const FieldView3D<const float> target) {
  return {
      checked_subtract(target.nx, target.halo.i_lower, target.halo.i_upper,
                       "3D target has invalid i extent or halo"),
      checked_subtract(target.ny, target.halo.j_lower, target.halo.j_upper,
                       "3D target has invalid j extent or halo"),
      checked_subtract(target.nz, target.halo.k_lower, target.halo.k_upper,
                       "3D target has invalid k extent or halo"),
  };
}

[[nodiscard]] ActiveTargetShape active_shape(const FieldView2D<const float> target) {
  return {
      checked_subtract(target.nx, target.halo.i_lower, target.halo.i_upper,
                       "2D target has invalid i extent or halo"),
      checked_subtract(target.ny, target.halo.j_lower, target.halo.j_upper,
                       "2D target has invalid j extent or halo"),
      1,
  };
}

void validate_packed_field(const PackedForcingField& packed) {
  require(packed.layout.nx > 0 && packed.layout.ny > 0 && packed.layout.nz > 0,
          ForcingApplyStatus::shape_mismatch,
          "packed forcing field requires non-empty canonical extents");
  require(packed.layout.value_count > 0,
          ForcingApplyStatus::shape_mismatch,
          "packed forcing field requires a positive value count");
  require(packed.values.size() == packed.layout.value_count,
          ForcingApplyStatus::shape_mismatch,
          "packed forcing value count does not match canonical layout");
}

[[nodiscard]] ForcingApplyReport make_ok_report(
    const std::size_t point_count,
    const ForcingApplyOperation operation,
    const bool would_modify_state,
    const bool synthetic,
    std::string detail,
    const std::size_t side_count = 1) {
  return {
      .field_count = 1,
      .side_count = side_count,
      .point_count = point_count,
      .status = ForcingApplyStatus::ok,
      .operation = operation,
      .would_modify_state = would_modify_state,
      .synthetic = synthetic,
      .detail = std::move(detail),
  };
}

void validate_boundary_shape(
    const ActiveTargetShape target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  validate_packed_field(packed);

  if (side == BoundarySide::i_lower || side == BoundarySide::i_upper) {
    require(packed.layout.nx <= target.nx && packed.layout.ny == target.ny &&
                packed.layout.nz == target.nz,
            ForcingApplyStatus::shape_mismatch,
            "i-side boundary pack must be [bdy_width, target_ny, target_nz]");
    return;
  }

  require(packed.layout.nx == target.nx && packed.layout.ny <= target.ny &&
              packed.layout.nz == target.nz,
          ForcingApplyStatus::shape_mismatch,
          "j-side boundary pack must be [target_nx, bdy_width, target_nz]");
}

[[nodiscard]] std::size_t boundary_active_i_begin(
    const ActiveTargetShape target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  if (side == BoundarySide::i_upper) {
    return target.nx - packed.layout.nx;
  }
  return 0;
}

[[nodiscard]] std::size_t boundary_active_j_begin(
    const ActiveTargetShape target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  if (side == BoundarySide::j_upper) {
    return target.ny - packed.layout.ny;
  }
  return 0;
}

[[nodiscard]] std::size_t checked_begin(
    const std::int32_t value,
    const char* message) {
  require(value >= 0, ForcingApplyStatus::invalid_range, message);
  return static_cast<std::size_t>(value);
}

void validate_delta_range(
    const ActiveTargetShape target,
    const PackedForcingField& delta,
    const SyntheticNudgingDeltaConfig config) {
  validate_packed_field(delta);
  const auto i_begin =
      checked_begin(config.active_i_begin, "synthetic delta i begin must be non-negative");
  const auto j_begin =
      checked_begin(config.active_j_begin, "synthetic delta j begin must be non-negative");
  const auto k_begin =
      checked_begin(config.active_k_begin, "synthetic delta k begin must be non-negative");

  require(i_begin <= target.nx && delta.layout.nx <= target.nx - i_begin,
          ForcingApplyStatus::invalid_range,
          "synthetic delta i range exceeds target active extent");
  require(j_begin <= target.ny && delta.layout.ny <= target.ny - j_begin,
          ForcingApplyStatus::invalid_range,
          "synthetic delta j range exceeds target active extent");
  require(k_begin <= target.nz && delta.layout.nz <= target.nz - k_begin,
          ForcingApplyStatus::invalid_range,
          "synthetic delta k range exceeds target active extent");
}

[[nodiscard]] std::int32_t checked_storage_index(
    const std::int32_t halo_begin,
    const std::int32_t active_begin,
    const std::size_t offset,
    const char* message) {
  const auto max_offset =
      static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max()) -
      static_cast<std::int64_t>(halo_begin) -
      static_cast<std::int64_t>(active_begin);
  require(max_offset >= 0, ForcingApplyStatus::invalid_range, message);
  require(offset <= static_cast<std::size_t>(max_offset),
          ForcingApplyStatus::invalid_range,
          message);
  return halo_begin + active_begin + static_cast<std::int32_t>(offset);
}

ForcingApplyReport apply_delta_to_view(
    FieldView3D<float> target,
    const PackedForcingField& delta,
    const SyntheticNudgingDeltaConfig config) {
  const auto target_shape = active_shape(FieldView3D<const float>{
      target.data,
      target.nx,
      target.ny,
      target.nz,
      target.stride_i,
      target.stride_k,
      target.stride_j,
      target.halo});
  validate_delta_range(target_shape, delta, config);

  for (std::size_t j = 0; j < delta.layout.ny; ++j) {
    const auto target_j = checked_storage_index(
        target.halo.j_lower, config.active_j_begin, j,
        "synthetic delta j index exceeds int32 storage index");
    for (std::size_t k = 0; k < delta.layout.nz; ++k) {
      const auto target_k = checked_storage_index(
          target.halo.k_lower, config.active_k_begin, k,
          "synthetic delta k index exceeds int32 storage index");
      for (std::size_t i = 0; i < delta.layout.nx; ++i) {
        const auto target_i = checked_storage_index(
            target.halo.i_lower, config.active_i_begin, i,
            "synthetic delta i index exceeds int32 storage index");
        const auto source = ((j * delta.layout.nz) + k) * delta.layout.nx + i;
        target(target_i, target_j, target_k) += delta.values[source];
      }
    }
  }

  return make_ok_report(
      delta.layout.value_count,
      ForcingApplyOperation::synthetic_nudging_delta,
      true,
      true,
      "applied explicit synthetic nudging delta; no WRF spectral nudging formula");
}

ForcingApplyReport apply_delta_to_view(
    FieldView2D<float> target,
    const PackedForcingField& delta,
    const SyntheticNudgingDeltaConfig config) {
  require(delta.layout.nz == 1,
          ForcingApplyStatus::shape_mismatch,
          "2D synthetic delta must use canonical nz=1");
  const auto target_shape = active_shape(FieldView2D<const float>{
      target.data,
      target.nx,
      target.ny,
      target.stride_i,
      target.stride_j,
      target.halo});
  validate_delta_range(target_shape, delta, config);

  for (std::size_t j = 0; j < delta.layout.ny; ++j) {
    const auto target_j = checked_storage_index(
        target.halo.j_lower, config.active_j_begin, j,
        "synthetic delta j index exceeds int32 storage index");
    for (std::size_t i = 0; i < delta.layout.nx; ++i) {
      const auto target_i = checked_storage_index(
          target.halo.i_lower, config.active_i_begin, i,
          "synthetic delta i index exceeds int32 storage index");
      const auto source = (j * delta.layout.nz) * delta.layout.nx + i;
      target(target_i, target_j) += delta.values[source];
    }
  }

  return make_ok_report(
      delta.layout.value_count,
      ForcingApplyOperation::synthetic_nudging_delta,
      true,
      true,
      "applied explicit synthetic 2D nudging delta; no WRF spectral nudging formula");
}

ForcingApplyReport copy_boundary_to_view(
    FieldView3D<float> target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  const auto target_shape = active_shape(FieldView3D<const float>{
      target.data,
      target.nx,
      target.ny,
      target.nz,
      target.stride_i,
      target.stride_k,
      target.stride_j,
      target.halo});
  validate_boundary_shape(target_shape, side, packed);

  const auto active_i_begin = boundary_active_i_begin(target_shape, side, packed);
  const auto active_j_begin = boundary_active_j_begin(target_shape, side, packed);
  for (std::size_t j = 0; j < packed.layout.ny; ++j) {
    const auto target_j = checked_storage_index(
        target.halo.j_lower, static_cast<std::int32_t>(active_j_begin), j,
        "boundary copy j index exceeds int32 storage index");
    for (std::size_t k = 0; k < packed.layout.nz; ++k) {
      const auto target_k = checked_storage_index(
          target.halo.k_lower, 0, k,
          "boundary copy k index exceeds int32 storage index");
      for (std::size_t i = 0; i < packed.layout.nx; ++i) {
        const auto target_i = checked_storage_index(
            target.halo.i_lower, static_cast<std::int32_t>(active_i_begin), i,
            "boundary copy i index exceeds int32 storage index");
        const auto source = ((j * packed.layout.nz) + k) * packed.layout.nx + i;
        target(target_i, target_j, target_k) = packed.values[source];
      }
    }
  }

  return make_ok_report(
      packed.layout.value_count,
      ForcingApplyOperation::direct_boundary_copy_skeleton,
      true,
      false,
      "copied packed boundary values directly into active edge window; skeleton only, no WRF lateral relaxation formula");
}

ForcingApplyReport copy_boundary_to_view(
    FieldView2D<float> target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  require(packed.layout.nz == 1,
          ForcingApplyStatus::shape_mismatch,
          "2D boundary copy pack must use canonical nz=1");
  const auto target_shape = active_shape(FieldView2D<const float>{
      target.data,
      target.nx,
      target.ny,
      target.stride_i,
      target.stride_j,
      target.halo});
  validate_boundary_shape(target_shape, side, packed);

  const auto active_i_begin = boundary_active_i_begin(target_shape, side, packed);
  const auto active_j_begin = boundary_active_j_begin(target_shape, side, packed);
  for (std::size_t j = 0; j < packed.layout.ny; ++j) {
    const auto target_j = checked_storage_index(
        target.halo.j_lower, static_cast<std::int32_t>(active_j_begin), j,
        "boundary copy j index exceeds int32 storage index");
    for (std::size_t i = 0; i < packed.layout.nx; ++i) {
      const auto target_i = checked_storage_index(
          target.halo.i_lower, static_cast<std::int32_t>(active_i_begin), i,
          "boundary copy i index exceeds int32 storage index");
      const auto source = (j * packed.layout.nz) * packed.layout.nx + i;
      target(target_i, target_j) = packed.values[source];
    }
  }

  return make_ok_report(
      packed.layout.value_count,
      ForcingApplyOperation::direct_boundary_copy_skeleton,
      true,
      false,
      "copied packed 2D boundary values directly into active edge window; skeleton only, no WRF lateral relaxation formula");
}

[[nodiscard]] bool same_field(
    const std::string_view lhs,
    const std::string_view rhs) noexcept {
  return lhs == rhs;
}

[[nodiscard]] bool is_x_boundary_side(const BoundarySide side) noexcept {
  return side == BoundarySide::i_lower || side == BoundarySide::i_upper;
}

[[nodiscard]] std::string krosa_boundary_variable_name(
    const std::string_view field_name,
    const std::string_view suffix) {
  std::string name(field_name);
  name += suffix;
  return name;
}

[[nodiscard]] PackedForcingField pack_krosa_boundary_slice(
    const ForcingTimeSlice& slice,
    const BoundarySide side) {
  if (is_x_boundary_side(side)) {
    return pack_boundary_x_side_raw_to_canonical(
        slice.values,
        slice.metadata.slice_shape);
  }

  return pack_boundary_y_side_raw_to_canonical(
      slice.values,
      slice.metadata.slice_shape);
}

}  // namespace

std::string_view forcing_apply_status_name(const ForcingApplyStatus status) noexcept {
  switch (status) {
    case ForcingApplyStatus::ok:
      return "ok";
    case ForcingApplyStatus::shape_mismatch:
      return "shape_mismatch";
    case ForcingApplyStatus::invalid_range:
      return "invalid_range";
    case ForcingApplyStatus::unsupported_field:
      return "unsupported_field";
  }
  return "unknown";
}

std::string_view forcing_apply_operation_name(
    const ForcingApplyOperation operation) noexcept {
  switch (operation) {
    case ForcingApplyOperation::validation_only:
      return "validation_only";
    case ForcingApplyOperation::synthetic_nudging_delta:
      return "synthetic_nudging_delta";
    case ForcingApplyOperation::direct_boundary_copy_skeleton:
      return "direct_boundary_copy_skeleton";
  }
  return "unknown";
}

ForcingApplyReport validate_boundary_pack_for_view(
    const FieldView3D<const float> target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  validate_boundary_shape(active_shape(target), side, packed);
  return make_ok_report(
      packed.layout.value_count,
      ForcingApplyOperation::validation_only,
      false,
      false,
      "validated boundary pack shape against 3D target; no state write");
}

ForcingApplyReport validate_boundary_pack_for_view(
    const FieldView2D<const float> target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  validate_boundary_shape(active_shape(target), side, packed);
  return make_ok_report(
      packed.layout.value_count,
      ForcingApplyOperation::validation_only,
      false,
      false,
      "validated boundary pack shape against 2D target; no state write");
}

ForcingApplyReport validate_boundary_pack_for_state(
    const State<float>& state,
    const std::string_view field_name,
    const BoundarySide side,
    const PackedForcingField& packed) {
  const auto view = state.view();
  if (same_field(field_name, "U")) {
    return validate_boundary_pack_for_view(view.u, side, packed);
  }
  if (same_field(field_name, "V")) {
    return validate_boundary_pack_for_view(view.v, side, packed);
  }
  if (same_field(field_name, "W")) {
    return validate_boundary_pack_for_view(view.w, side, packed);
  }
  if (same_field(field_name, "PH")) {
    return validate_boundary_pack_for_view(view.ph, side, packed);
  }
  if (same_field(field_name, "PHB")) {
    return validate_boundary_pack_for_view(view.phb, side, packed);
  }
  if (same_field(field_name, "T")) {
    return validate_boundary_pack_for_view(view.t, side, packed);
  }
  if (same_field(field_name, "P")) {
    return validate_boundary_pack_for_view(view.p, side, packed);
  }
  if (same_field(field_name, "PB")) {
    return validate_boundary_pack_for_view(view.pb, side, packed);
  }
  if (same_field(field_name, "QVAPOR")) {
    return validate_boundary_pack_for_view(view.qvapor, side, packed);
  }
  if (same_field(field_name, "QCLOUD")) {
    return validate_boundary_pack_for_view(view.qcloud, side, packed);
  }
  if (same_field(field_name, "QRAIN")) {
    return validate_boundary_pack_for_view(view.qrain, side, packed);
  }
  if (same_field(field_name, "QICE")) {
    return validate_boundary_pack_for_view(view.qice, side, packed);
  }
  if (same_field(field_name, "QSNOW")) {
    return validate_boundary_pack_for_view(view.qsnow, side, packed);
  }
  if (same_field(field_name, "QGRAUP")) {
    return validate_boundary_pack_for_view(view.qgraup, side, packed);
  }
  if (same_field(field_name, "QNICE")) {
    return validate_boundary_pack_for_view(view.qnice, side, packed);
  }
  if (same_field(field_name, "QNRAIN")) {
    return validate_boundary_pack_for_view(view.qnrain, side, packed);
  }
  if (same_field(field_name, "MU")) {
    return validate_boundary_pack_for_view(view.mu, side, packed);
  }
  if (same_field(field_name, "MUB")) {
    return validate_boundary_pack_for_view(view.mub, side, packed);
  }
  if (same_field(field_name, "PSFC")) {
    return validate_boundary_pack_for_view(view.psfc, side, packed);
  }
  if (same_field(field_name, "U10")) {
    return validate_boundary_pack_for_view(view.u10, side, packed);
  }
  if (same_field(field_name, "V10")) {
    return validate_boundary_pack_for_view(view.v10, side, packed);
  }
  if (same_field(field_name, "T2")) {
    return validate_boundary_pack_for_view(view.t2, side, packed);
  }
  if (same_field(field_name, "Q2")) {
    return validate_boundary_pack_for_view(view.q2, side, packed);
  }
  if (same_field(field_name, "RAINC")) {
    return validate_boundary_pack_for_view(view.rainc, side, packed);
  }
  if (same_field(field_name, "RAINNC")) {
    return validate_boundary_pack_for_view(view.rainnc, side, packed);
  }

  fail(ForcingApplyStatus::unsupported_field,
       "forcing apply does not know the requested state field");
}

ForcingApplyReport apply_boundary_copy_skeleton_to_view(
    const FieldView3D<float> target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  return copy_boundary_to_view(target, side, packed);
}

ForcingApplyReport apply_boundary_copy_skeleton_to_view(
    const FieldView2D<float> target,
    const BoundarySide side,
    const PackedForcingField& packed) {
  return copy_boundary_to_view(target, side, packed);
}

ForcingApplyReport apply_boundary_copy_skeleton_to_state(
    State<float>& state,
    const std::string_view field_name,
    const BoundarySide side,
    const PackedForcingField& packed) {
  auto view = state.view();
  if (same_field(field_name, "U")) {
    return apply_boundary_copy_skeleton_to_view(view.u, side, packed);
  }
  if (same_field(field_name, "V")) {
    return apply_boundary_copy_skeleton_to_view(view.v, side, packed);
  }
  if (same_field(field_name, "W")) {
    return apply_boundary_copy_skeleton_to_view(view.w, side, packed);
  }
  if (same_field(field_name, "PH")) {
    return apply_boundary_copy_skeleton_to_view(view.ph, side, packed);
  }
  if (same_field(field_name, "PHB")) {
    return apply_boundary_copy_skeleton_to_view(view.phb, side, packed);
  }
  if (same_field(field_name, "T")) {
    return apply_boundary_copy_skeleton_to_view(view.t, side, packed);
  }
  if (same_field(field_name, "P")) {
    return apply_boundary_copy_skeleton_to_view(view.p, side, packed);
  }
  if (same_field(field_name, "PB")) {
    return apply_boundary_copy_skeleton_to_view(view.pb, side, packed);
  }
  if (same_field(field_name, "QVAPOR")) {
    return apply_boundary_copy_skeleton_to_view(view.qvapor, side, packed);
  }
  if (same_field(field_name, "QCLOUD")) {
    return apply_boundary_copy_skeleton_to_view(view.qcloud, side, packed);
  }
  if (same_field(field_name, "QRAIN")) {
    return apply_boundary_copy_skeleton_to_view(view.qrain, side, packed);
  }
  if (same_field(field_name, "QICE")) {
    return apply_boundary_copy_skeleton_to_view(view.qice, side, packed);
  }
  if (same_field(field_name, "QSNOW")) {
    return apply_boundary_copy_skeleton_to_view(view.qsnow, side, packed);
  }
  if (same_field(field_name, "QGRAUP")) {
    return apply_boundary_copy_skeleton_to_view(view.qgraup, side, packed);
  }
  if (same_field(field_name, "QNICE")) {
    return apply_boundary_copy_skeleton_to_view(view.qnice, side, packed);
  }
  if (same_field(field_name, "QNRAIN")) {
    return apply_boundary_copy_skeleton_to_view(view.qnrain, side, packed);
  }
  if (same_field(field_name, "MU")) {
    return apply_boundary_copy_skeleton_to_view(view.mu, side, packed);
  }
  if (same_field(field_name, "MUB")) {
    return apply_boundary_copy_skeleton_to_view(view.mub, side, packed);
  }
  if (same_field(field_name, "PSFC")) {
    return apply_boundary_copy_skeleton_to_view(view.psfc, side, packed);
  }
  if (same_field(field_name, "U10")) {
    return apply_boundary_copy_skeleton_to_view(view.u10, side, packed);
  }
  if (same_field(field_name, "V10")) {
    return apply_boundary_copy_skeleton_to_view(view.v10, side, packed);
  }
  if (same_field(field_name, "T2")) {
    return apply_boundary_copy_skeleton_to_view(view.t2, side, packed);
  }
  if (same_field(field_name, "Q2")) {
    return apply_boundary_copy_skeleton_to_view(view.q2, side, packed);
  }
  if (same_field(field_name, "RAINC")) {
    return apply_boundary_copy_skeleton_to_view(view.rainc, side, packed);
  }
  if (same_field(field_name, "RAINNC")) {
    return apply_boundary_copy_skeleton_to_view(view.rainnc, side, packed);
  }

  fail(ForcingApplyStatus::unsupported_field,
       "forcing apply does not know the requested state field");
}

ForcingApplyReport apply_krosa_boundary_copy_skeleton_to_state(
    State<float>& state,
    const KrosaForcingReader& reader,
    const std::string_view field_name,
    const std::size_t record_index) {
  if (reader.kind() != KrosaForcingKind::boundary) {
    fail(ForcingApplyStatus::unsupported_field,
         "KROSA boundary copy skeleton requires a wrfbdy boundary reader");
  }
  if (!same_field(field_name, "U")) {
    fail(ForcingApplyStatus::unsupported_field,
         "KROSA boundary copy skeleton currently maps U_BXS/U_BXE/U_BYS/U_BYE only");
  }

  std::array<KrosaPackedBoundarySide, kKrosaBoundarySideMappings.size()> sides{};
  for (std::size_t index = 0; index < kKrosaBoundarySideMappings.size(); ++index) {
    const auto mapping = kKrosaBoundarySideMappings.at(index);
    auto variable_name = krosa_boundary_variable_name(field_name, mapping.suffix);
    auto slice = reader.read_float_time_slice(variable_name, record_index);
    sides.at(index) = {
        .side = mapping.side,
        .variable_name = std::move(variable_name),
        .packed = pack_krosa_boundary_slice(slice, mapping.side),
    };
  }

  for (const auto& side : sides) {
    (void)validate_boundary_pack_for_state(
        state,
        field_name,
        side.side,
        side.packed);
  }

  std::size_t point_count = 0;
  for (const auto& side : sides) {
    const auto report = apply_boundary_copy_skeleton_to_state(
        state,
        field_name,
        side.side,
        side.packed);
    point_count += report.point_count;
  }

  std::string detail = "applied KROSA ";
  detail += field_name;
  detail += " wrfbdy sides as direct_boundary_copy_skeleton: ";
  for (std::size_t index = 0; index < sides.size(); ++index) {
    if (index != 0) {
      detail += ", ";
    }
    detail += sides.at(index).variable_name;
  }
  detail +=
      "; skeleton direct copy only, not synthetic nudging and not WRF lateral relaxation";

  return {
      .field_count = 1,
      .side_count = sides.size(),
      .point_count = point_count,
      .status = ForcingApplyStatus::ok,
      .operation = ForcingApplyOperation::direct_boundary_copy_skeleton,
      .would_modify_state = true,
      .synthetic = false,
      .detail = std::move(detail),
  };
}

ForcingApplyReport apply_synthetic_nudging_delta(
    State<float>& state,
    const std::string_view field_name,
    const PackedForcingField& synthetic_delta,
    const SyntheticNudgingDeltaConfig config) {
  auto view = state.view();
  if (same_field(field_name, "U")) {
    return apply_delta_to_view(view.u, synthetic_delta, config);
  }
  if (same_field(field_name, "V")) {
    return apply_delta_to_view(view.v, synthetic_delta, config);
  }
  if (same_field(field_name, "W")) {
    return apply_delta_to_view(view.w, synthetic_delta, config);
  }
  if (same_field(field_name, "PH")) {
    return apply_delta_to_view(view.ph, synthetic_delta, config);
  }
  if (same_field(field_name, "PHB")) {
    return apply_delta_to_view(view.phb, synthetic_delta, config);
  }
  if (same_field(field_name, "T")) {
    return apply_delta_to_view(view.t, synthetic_delta, config);
  }
  if (same_field(field_name, "P")) {
    return apply_delta_to_view(view.p, synthetic_delta, config);
  }
  if (same_field(field_name, "PB")) {
    return apply_delta_to_view(view.pb, synthetic_delta, config);
  }
  if (same_field(field_name, "QVAPOR")) {
    return apply_delta_to_view(view.qvapor, synthetic_delta, config);
  }
  if (same_field(field_name, "QCLOUD")) {
    return apply_delta_to_view(view.qcloud, synthetic_delta, config);
  }
  if (same_field(field_name, "QRAIN")) {
    return apply_delta_to_view(view.qrain, synthetic_delta, config);
  }
  if (same_field(field_name, "QICE")) {
    return apply_delta_to_view(view.qice, synthetic_delta, config);
  }
  if (same_field(field_name, "QSNOW")) {
    return apply_delta_to_view(view.qsnow, synthetic_delta, config);
  }
  if (same_field(field_name, "QGRAUP")) {
    return apply_delta_to_view(view.qgraup, synthetic_delta, config);
  }
  if (same_field(field_name, "QNICE")) {
    return apply_delta_to_view(view.qnice, synthetic_delta, config);
  }
  if (same_field(field_name, "QNRAIN")) {
    return apply_delta_to_view(view.qnrain, synthetic_delta, config);
  }
  if (same_field(field_name, "MU")) {
    return apply_delta_to_view(view.mu, synthetic_delta, config);
  }
  if (same_field(field_name, "MUB")) {
    return apply_delta_to_view(view.mub, synthetic_delta, config);
  }
  if (same_field(field_name, "PSFC")) {
    return apply_delta_to_view(view.psfc, synthetic_delta, config);
  }
  if (same_field(field_name, "U10")) {
    return apply_delta_to_view(view.u10, synthetic_delta, config);
  }
  if (same_field(field_name, "V10")) {
    return apply_delta_to_view(view.v10, synthetic_delta, config);
  }
  if (same_field(field_name, "T2")) {
    return apply_delta_to_view(view.t2, synthetic_delta, config);
  }
  if (same_field(field_name, "Q2")) {
    return apply_delta_to_view(view.q2, synthetic_delta, config);
  }
  if (same_field(field_name, "RAINC")) {
    return apply_delta_to_view(view.rainc, synthetic_delta, config);
  }
  if (same_field(field_name, "RAINNC")) {
    return apply_delta_to_view(view.rainnc, synthetic_delta, config);
  }

  fail(ForcingApplyStatus::unsupported_field,
       "forcing apply does not know the requested state field");
}

}  // namespace tywrf::io
