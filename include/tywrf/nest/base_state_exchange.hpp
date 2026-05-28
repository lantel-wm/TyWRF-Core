#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/nest/nest_interface.hpp"

#include <cstdint>
#include <string_view>

namespace tywrf::nest {

template <typename Real>
struct ExposedBaseStateViews {
  FieldView3D<Real> phb;
  FieldView2D<Real> mub;
  FieldView3D<Real> pb;
  FieldView3D<Real> t_init;
  FieldView3D<Real> alb;
  FieldView2D<Real> ht;
};

struct ExposedBaseStateExchangeOptions {
  bool rebalance = false;
};

struct ExposedBaseStateExchangeReport {
  NestResult result{NestStatus::ok, "ok"};
  std::string_view source = "diagnostic_exposed_base_state_exchange";
  std::string_view disposition =
      "diagnostic_only_no_gate_no_selected_field_numerics";
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool integrator_output = false;
  bool selected_field_numerics_enabled = false;
  bool enables_selected_field_numerics = false;
  bool rebalance_requested = false;
  bool rebalance_applied = false;

  std::uint32_t exposed_region_count = 0;
  std::uint64_t exposed_mass_cell_count = 0;
  std::uint64_t exposed_surface_cell_count = 0;
  std::uint64_t exposed_w_full_column_count = 0;

  bool wrote_overlap = false;
  bool wrote_halo = false;
  bool wrote_phb = false;
  bool wrote_mub = false;
  bool wrote_ht = false;
  bool wrote_pb = false;
  bool wrote_t_init = false;
  bool wrote_alb = false;

  std::uint64_t phb_written_point_count = 0;
  std::uint64_t mub_written_cell_count = 0;
  std::uint64_t ht_written_cell_count = 0;
  std::uint64_t direct_write_point_count = 0;

  std::uint64_t pb_recompute_mark_count = 0;
  std::uint64_t t_init_recompute_mark_count = 0;
  std::uint64_t alb_recompute_mark_count = 0;
  std::uint64_t recompute_mark_count = 0;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// Diagnostic-only helper for the WRF moving-nest exposed-cell base-state slice.
// It writes only child-shaped exposed PHB, MUB, and HT from already interpolated
// source views. PB, T_INIT, and ALB are not copied; their exposed mass points are
// reported as recompute marks for a future provider-backed path. This helper
// performs no NetCDF I/O and does not enable selected-field numerics or produce
// a gate candidate.
[[nodiscard]] ExposedBaseStateExchangeReport apply_exposed_base_state_exchange(
    const RemapPlan& remap_plan,
    const ExposedBaseStateViews<const float>& source_interpolated_base_state,
    const ExposedBaseStateViews<float>& child_base_state,
    ExposedBaseStateExchangeOptions options = {}) noexcept;

}  // namespace tywrf::nest
