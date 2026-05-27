#pragma once

#include "tywrf/nest/nest_interface.hpp"
#include "tywrf/state.hpp"

#include <cstdint>

namespace tywrf::nest {

struct ChildStateRemapReport {
  NestResult result{NestStatus::ok, "ok"};
  std::uint32_t copied_field_count = 0;
  std::uint64_t copied_point_count = 0;
  std::uint32_t parent_fill_field_count = 0;
  std::uint64_t parent_fill_point_count = 0;
  bool needs_parent_fill = false;
  bool filled_exposed_cells = false;
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

// Staging helper for moving nests: copies the old-child overlap, then fills
// newly exposed active cells from an already interpolated child-shaped parent
// state. Halo cells are left untouched; no parent-child interpolation is done
// here.
[[nodiscard]] ChildStateRemapReport remap_child_state_overlap_with_parent_fill(
    const RemapPlan& plan,
    const State<float>& old_child,
    const State<float>& parent_filled_child,
    State<float>& new_child) noexcept;

}  // namespace tywrf::nest
