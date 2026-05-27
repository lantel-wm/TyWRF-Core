#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/nest/nest_interface.hpp"

#include <cstdint>

namespace tywrf::nest {

struct PostStartDomainLevel2Views {
  FieldView3D<const float> u;
  FieldView3D<const float> v;
  FieldView3D<const float> t;
  FieldView3D<const float> w;
  FieldView3D<const float> ph;
  FieldView2D<const float> mu;
};

struct PostStartDomainLevel1Views {
  FieldView3D<float> u;
  FieldView3D<float> v;
  FieldView3D<float> t;
  FieldView3D<float> w;
  FieldView3D<float> ph;
  FieldView2D<float> mu;
};

struct OptionalTimeLevelSyncTke {
  FieldView3D<const float> level2;
  FieldView3D<float> level1;
  bool present = false;
};

struct TimeLevelSyncReport {
  NestResult result{NestStatus::ok, "ok"};
  std::uint32_t copied_field_count = 0;
  std::uint64_t copied_point_count = 0;
  std::uint64_t u_point_count = 0;
  std::uint64_t v_point_count = 0;
  std::uint64_t t_point_count = 0;
  std::uint64_t w_point_count = 0;
  std::uint64_t ph_point_count = 0;
  std::uint64_t mu_point_count = 0;
  std::uint64_t tke_point_count = 0;
  bool copied_optional_tke = false;
  bool copied_halo_cells = false;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// WRF start_domain initializes level 2 for newly exposed moving-nest cells.
// This helper mirrors WRF's post-start-domain level2 -> level1 copy only for
// active cells outside the new overlap window. Overlap cells and halos are left
// untouched.
[[nodiscard]] TimeLevelSyncReport copy_post_start_domain_level2_to_level1_exposed_cells(
    const RemapPlan& plan,
    const PostStartDomainLevel2Views& level2,
    const PostStartDomainLevel1Views& level1,
    OptionalTimeLevelSyncTke optional_tke = {}) noexcept;

}  // namespace tywrf::nest
