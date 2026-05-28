#pragma once

#include "tywrf/dynamics/base_state.hpp"
#include "tywrf/nest/base_state_exchange.hpp"

#include <cstdint>
#include <string_view>

namespace tywrf::nest {

struct ExposedBaseStateExchangeAdapterInputs {
  dynamics::BaseStateVerticalCoefficientView c3h;
  dynamics::BaseStateVerticalCoefficientView c4h;
  float p_top_pa = 0.0F;
};

struct ExposedBaseStateExchangeAdapterOptions {
  ExposedBaseStateExchangeOptions exchange{};
  dynamics::KrosaMassBaseStateReconstructionOptions recompute{};
};

struct ExposedBaseStateExchangeAdapterReport {
  NestResult result{NestStatus::ok, "ok"};
  NestResult exchange_result{NestStatus::ok, "ok"};
  NestResult recompute_result{NestStatus::ok, "ok"};
  std::string_view source = "diagnostic_exposed_base_state_exchange_adapter";
  std::string_view disposition =
      "diagnostic_only_staging_no_gate_no_integrator_no_selected_field_numerics";
  std::string_view exchange_source = "diagnostic_exposed_base_state_exchange";
  std::string_view recompute_source = "exposed_mub_base_state_recompute_provider";

  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool integrator_output = false;
  bool selected_field_numerics_enabled = false;
  bool enables_selected_field_numerics = false;
  bool writes_netcdf = false;
  bool writes_candidate = false;

  bool called_d68_exchange = false;
  bool called_d69_recompute = false;
  bool rebalance_requested = false;
  bool rebalance_applied = false;

  std::string_view ht_source_name = "HT";
  std::string_view ht_diagnostic_label = "HGT";
  bool ht_is_hgt_alias = true;
  bool creates_terrain_owner = false;
  bool terrain_owner_created = false;

  std::int32_t active_nx = 0;
  std::int32_t active_ny = 0;
  std::int32_t active_nz = 0;
  std::int32_t c3h_count = 0;
  std::int32_t c4h_count = 0;

  std::uint32_t exchange_exposed_region_count = 0;
  std::uint32_t recompute_exposed_region_count = 0;
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

  std::uint64_t exchange_recompute_mark_count = 0;
  std::uint64_t recomputed_point_count = 0;
  std::uint64_t pb_recomputed_point_count = 0;
  std::uint64_t t_init_recomputed_point_count = 0;
  std::uint64_t alb_recomputed_point_count = 0;
  std::uint64_t invalid_column_count = 0;
  std::uint64_t invalid_point_count = 0;
  std::uint64_t unsupported_stratosphere_point_count = 0;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// Diagnostic-only adapter that stages the D68 exposed base-state exchange
// (PHB/MUB/HT) and then uses the D69 exposed-MUB provider recompute for
// provider-owned PB/T_INIT/ALB. HT is reported only as the WRF HGT diagnostic
// alias; this adapter does not create a second terrain owner, write NetCDF,
// produce a gate candidate, or enable selected-field numerics.
[[nodiscard]] ExposedBaseStateExchangeAdapterReport
apply_exposed_base_state_exchange_adapter(
    const RemapPlan& remap_plan,
    const ExposedBaseStateViews<const float>& source_interpolated_base_state,
    const ExposedBaseStateViews<float>& child_base_state,
    const ExposedBaseStateExchangeAdapterInputs& inputs,
    ExposedBaseStateExchangeAdapterOptions options = {}) noexcept;

}  // namespace tywrf::nest
