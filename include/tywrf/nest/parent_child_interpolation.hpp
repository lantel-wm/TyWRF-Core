#pragma once

#include "tywrf/nest/state_exchange.hpp"
#include "tywrf/state.hpp"

#include <cstdint>

namespace tywrf::nest {

enum class ParentChildInterpolationMethod : std::uint8_t {
  nearest_neighbor,
  bilinear,
};

struct ParentChildInterpolationConfig {
  ParentChildInterpolationMethod method = ParentChildInterpolationMethod::bilinear;
  bool allow_owned_overlap_regions = false;
  bool allow_owned_halo_regions = false;
};

struct ParentChildInterpolationReport {
  NestResult result{NestStatus::ok, "ok"};
  ParentChildInterpolationMethod method = ParentChildInterpolationMethod::bilinear;
  std::uint32_t requested_field_count = 0;
  std::uint32_t interpolated_field_count = 0;
  std::uint32_t interpolated_region_count = 0;
  std::uint64_t interpolated_horizontal_cell_count = 0;
  std::uint64_t interpolated_point_count = 0;
  bool wrote_overlap = false;
  bool wrote_halo = false;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// Host-side first-slice parent-to-child fill for moving-nest exposed cells.
// The exchange plan selects target child active regions. Coordinates are mapped
// from child active indices to parent active indices using the target nest
// position and parent_grid_ratio, then sampled by nearest-neighbor or bilinear
// horizontal interpolation. Vertical levels are copied by matching active k.
// Supported fields are U, V, MU, and QVAPOR. No NetCDF, logging, allocation, or
// overlap/halo writes are performed unless the plan explicitly owns those
// regions and config allows them.
[[nodiscard]] ParentChildInterpolationReport interpolate_parent_to_exposed_child(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& target_position,
    const StateExchangePlan& exchange_plan,
    const StateView<const float>& parent_state,
    const StateView<float>& child_state,
    ParentChildInterpolationConfig config = {}) noexcept;

}  // namespace tywrf::nest
