#pragma once

#include "tywrf/grid.hpp"
#include "tywrf/nest/base_state_exchange.hpp"
#include "tywrf/nest/nest_interface.hpp"
#include "tywrf/state.hpp"

#include <cstdint>
#include <limits>
#include <string_view>

namespace tywrf::nest {

struct BaseStateSourceStagingOptions {
  float mask_value = std::numeric_limits<float>::quiet_NaN();
};

struct BaseStateSourceStagingReport {
  NestResult result{NestStatus::ok, "ok"};
  std::string_view source = "explicit_base_state_source_staging_provider";
  std::string_view disposition =
      "diagnostic_only_source_staging_no_gate_no_integrator_no_candidate";
  std::string_view source_shape = "child_shaped_explicit_source_views";

  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool integrator_output = false;
  bool selected_field_numerics_enabled = false;
  bool enables_selected_field_numerics = false;
  bool writes_netcdf = false;
  bool writes_candidate = false;
  bool candidate_buffers_preserved = true;
  bool owns_staging_buffers = false;
  bool allocated_buffers = false;
  bool uses_reference_end_truth = false;
  bool uses_direct_p_shortcut = false;
  bool reads_direct_p = false;

  std::int32_t active_nx = 0;
  std::int32_t active_ny = 0;
  std::int32_t mass_nz = 0;
  std::int32_t full_nz = 0;
  std::int32_t expected_mass_level_count = 0;
  std::int32_t expected_full_level_count = 0;

  std::uint32_t exposed_region_count = 0;
  std::uint64_t active_mass_cell_count = 0;
  std::uint64_t active_mass_point_count = 0;
  std::uint64_t active_surface_cell_count = 0;
  std::uint64_t active_w_full_point_count = 0;
  std::uint64_t exposed_mass_cell_count = 0;
  std::uint64_t exposed_mass_point_count = 0;
  std::uint64_t exposed_surface_cell_count = 0;
  std::uint64_t exposed_w_full_column_count = 0;
  std::uint64_t exposed_w_full_point_count = 0;
  std::uint64_t masked_mass_cell_count = 0;
  std::uint64_t masked_mass_point_count = 0;
  std::uint64_t masked_surface_cell_count = 0;
  std::uint64_t masked_w_full_column_count = 0;
  std::uint64_t masked_w_full_point_count = 0;

  std::uint64_t staged_phb_point_count = 0;
  std::uint64_t staged_mub_cell_count = 0;
  std::uint64_t staged_ht_cell_count = 0;
  std::uint64_t staged_pb_point_count = 0;
  std::uint64_t staged_t_init_point_count = 0;
  std::uint64_t staged_alb_point_count = 0;
  std::uint64_t staged_value_count = 0;
  std::uint64_t active_masked_value_count = 0;
  std::uint64_t halo_masked_value_count = 0;

  std::uint64_t invalid_exposed_value_count = 0;
  std::uint64_t invalid_exposed_phb_point_count = 0;
  std::uint64_t invalid_exposed_mub_cell_count = 0;
  std::uint64_t invalid_exposed_ht_cell_count = 0;
  std::uint64_t invalid_exposed_pb_point_count = 0;
  std::uint64_t invalid_exposed_t_init_point_count = 0;
  std::uint64_t invalid_exposed_alb_point_count = 0;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

class BaseStateSourceStagingProvider {
 public:
  BaseStateSourceStagingProvider() = default;

  [[nodiscard]] BaseStateSourceStagingReport stage(
      const Grid& child_grid,
      const RemapPlan& remap_plan,
      const ExposedBaseStateViews<const float>& explicit_source_views,
      BaseStateSourceStagingOptions options = {});

  [[nodiscard]] ExposedBaseStateViews<const float> views() const noexcept;

 private:
  void clear_views() noexcept;

  FieldStorage3D<float> phb_;
  FieldStorage2D<float> mub_;
  FieldStorage3D<float> pb_;
  FieldStorage3D<float> t_init_;
  FieldStorage3D<float> alb_;
  FieldStorage2D<float> ht_;
};

}  // namespace tywrf::nest
