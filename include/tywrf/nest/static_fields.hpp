#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/nest/nest_interface.hpp"

#include <cstdint>

namespace tywrf::nest {

struct MovingNestStaticRefreshReport {
  NestResult result{NestStatus::ok, "ok"};
  std::uint64_t overlap_cell_count = 0;
  std::uint64_t exposed_cell_count = 0;
  std::uint64_t coordinate_extrapolated_cell_count = 0;
  std::uint64_t parent_hgt_interpolated_cell_count = 0;
  bool copied_overlap = false;
  bool filled_exposed_coordinates = false;
  bool filled_exposed_hgt = false;
  bool uses_reference_end = false;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// Refreshes moving-nest mass-grid static fields without using end-state truth.
//
// The overlap footprint follows remap_plan.surface exactly: XLAT, XLONG, and
// HGT are shifted from the d02 start-state. Newly exposed XLAT/XLONG cells are
// filled by local linear extrapolation of the d02 start grid at the shifted
// source coordinate. Newly exposed HGT cells are filled from d01 start-state
// HGT by parent-grid bilinear interpolation at the target pose.
[[nodiscard]] MovingNestStaticRefreshReport refresh_moving_nest_static_fields(
    const ParentChildDescriptor& descriptor,
    const ParentChildPosition& target_position,
    const RemapPlan& remap_plan,
    const FieldView2D<const float>& d02_start_xlat,
    const FieldView2D<const float>& d02_start_xlong,
    const FieldView2D<const float>& d02_start_hgt,
    const FieldView2D<const float>& d01_start_hgt,
    const FieldView2D<float>& output_xlat,
    const FieldView2D<float>& output_xlong,
    const FieldView2D<float>& output_hgt) noexcept;

}  // namespace tywrf::nest
