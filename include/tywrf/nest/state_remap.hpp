#pragma once

#include "tywrf/nest/nest_interface.hpp"
#include "tywrf/state.hpp"

#include <cstdint>

namespace tywrf::nest {

struct ChildStateRemapReport {
  NestResult result{NestStatus::ok, "ok"};
  std::uint32_t copied_field_count = 0;
  std::uint64_t copied_point_count = 0;
  bool needs_parent_fill = false;
  bool copied_halo_cells = false;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// Copies only the overlapping active-array footprint described by plan.
// Halo cells and newly exposed child cells are left untouched; callers must
// fill exposed cells from the parent before treating the new state as physical.
[[nodiscard]] ChildStateRemapReport remap_child_state_overlap_only(
    const RemapPlan& plan,
    const State<float>& old_child,
    State<float>& new_child) noexcept;

}  // namespace tywrf::nest
