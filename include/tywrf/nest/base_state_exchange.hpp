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

template <typename Real>
struct NormalCandidateBaseStateSourceViews {
  FieldView3D<Real> phb;
  FieldView2D<Real> mub;
  FieldView3D<Real> pb;
  FieldView2D<Real> ht;
};

template <typename Real>
struct NormalCandidateBaseStateTargetViews {
  FieldView3D<Real> phb;
  FieldView2D<Real> mub;
  FieldView3D<Real> pb;
  FieldView2D<Real> ht;
};

struct NormalCandidateBaseStateExchangeReport {
  NestResult result{NestStatus::ok, "ok"};
  std::string_view source = "normal_selected_field_base_state_producer";
  std::string_view disposition =
      "non_oracle_candidate_base_state_fields_only_no_direct_p_no_gate_claim";
  std::string_view source_origin =
      "krosa_base_state_provider+moved_candidate_HGT+same_cycle_vertical_metadata";
  std::string_view written_fields = "PHB,MUB,PB,HGT";

  bool diagnostic_only = false;
  bool normal_candidate_producer = true;
  bool gate_candidate = false;
  bool integrator_output = false;
  bool writes_candidate = true;
  bool writes_netcdf = false;
  bool selected_field_numerics_enabled = true;
  bool enables_selected_field_numerics = false;
  bool uses_reference_end_truth = false;
  bool uses_direct_p_shortcut = false;
  bool reads_direct_p = false;
  bool writes_p = false;
  bool no_gate_pass_claim = true;
  bool no_00_20_progression = true;

  std::uint32_t exposed_region_count = 0;
  std::uint64_t exposed_mass_cell_count = 0;
  std::uint64_t exposed_surface_cell_count = 0;
  std::uint64_t exposed_w_full_column_count = 0;

  bool wrote_overlap = false;
  bool wrote_halo = false;
  bool wrote_phb = false;
  bool wrote_mub = false;
  bool wrote_pb = false;
  bool wrote_ht = false;
  bool wrote_t_init = false;
  bool wrote_alb = false;

  std::uint64_t phb_written_point_count = 0;
  std::uint64_t mub_written_cell_count = 0;
  std::uint64_t pb_written_point_count = 0;
  std::uint64_t ht_written_cell_count = 0;
  std::uint64_t direct_write_point_count = 0;

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

// Normal selected-field candidate base-state producer for moving-nest exposed
// cells. Sources must already be child-shaped and non-oracle, for example from
// KROSA base-state reconstruction using moved candidate HGT plus same-cycle
// vertical metadata. This writes only normal candidate PHB/MUB/PB and the
// normal static HT/HGT view for exposed cells. It has no perturbation-P view and
// cannot directly patch State::P.
[[nodiscard]] NormalCandidateBaseStateExchangeReport
apply_normal_candidate_base_state_exchange(
    const RemapPlan& remap_plan,
    const NormalCandidateBaseStateSourceViews<const float>& source_base_state,
    const NormalCandidateBaseStateTargetViews<float>& candidate_base_state) noexcept;

}  // namespace tywrf::nest
