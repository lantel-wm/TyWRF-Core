#pragma once

#include "tywrf/dynamics/base_state.hpp"
#include "tywrf/grid.hpp"
#include "tywrf/io/pressure_refresh_io.hpp"
#include "tywrf/nest/nest_interface.hpp"
#include "tywrf/state.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace tywrf::dynamics {

struct KrosaBaseStateProviderViews {
  FieldView3D<const float> pb;
  FieldView3D<const float> t_init;
  FieldView2D<const float> mub;
  FieldView3D<const float> alb;
  FieldView3D<const float> phb;
};

struct KrosaBaseStateProviderTerrainOverride {
  FieldView2D<const float> terrain_height_m{};
  std::string_view source_name = "moved_candidate_HGT";
  std::string_view provenance = "override:moved_candidate_HGT";
};

struct KrosaBaseStateProviderReport {
  nest::NestResult result{nest::NestStatus::ok, "ok"};
  std::string_view source = "base_state_reconstruction_provider";
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool integrator_output = false;
  bool calls_pressure_refresh_compute = false;

  std::int32_t active_nx = 0;
  std::int32_t active_ny = 0;
  std::int32_t mass_nz = 0;
  std::int32_t full_nz = 0;
  std::int32_t expected_mass_level_count = 0;
  std::int32_t expected_full_level_count = 0;
  std::int32_t c3f_count = 0;
  std::int32_t c4f_count = 0;
  std::int32_t c3h_count = 0;
  std::int32_t c4h_count = 0;
  std::int32_t terrain_nx = 0;
  std::int32_t terrain_ny = 0;
  std::size_t expected_terrain_point_count = 0;
  std::size_t terrain_point_count = 0;
  bool terrain_override_used = false;
  std::string terrain_source_name;
  std::string terrain_provenance;
  bool p_top_present = false;
  bool allocated_buffers = false;
  bool wrote_pb = false;
  bool wrote_t_init = false;
  bool wrote_mub = false;
  bool wrote_alb = false;
  bool wrote_phb = false;
  KrosaMassBaseStateReconstructionReport reconstruction{};

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

struct KrosaExposedBaseStateRecomputeInputs {
  FieldView2D<const float> mub;
  BaseStateVerticalCoefficientView c3h;
  BaseStateVerticalCoefficientView c4h;
  float p_top_pa = 0.0F;
};

struct KrosaExposedBaseStateRecomputeOutputs {
  FieldView3D<float> pb;
  FieldView3D<float> t_init;
  FieldView3D<float> alb;
};

struct KrosaExposedBaseStateRecomputeReport {
  nest::NestResult result{nest::NestStatus::ok, "ok"};
  const char* source = "exposed_mub_base_state_recompute_provider";
  const char* disposition =
      "provider_buffer_only_no_gate_no_selected_field_numerics";
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool integrator_output = false;
  bool selected_field_numerics_enabled = false;
  bool calls_pressure_refresh_compute = false;
  bool read_interpolated_mub = false;
  bool reads_hgt = false;
  bool touched_phb = false;
  bool wrote_pb = false;
  bool wrote_t_init = false;
  bool wrote_alb = false;
  bool wrote_state_t = false;
  bool wrote_state_alb = false;

  std::int32_t active_nx = 0;
  std::int32_t active_ny = 0;
  std::int32_t active_nz = 0;
  std::int32_t c3h_count = 0;
  std::int32_t c4h_count = 0;
  std::uint32_t exposed_region_count = 0;
  std::uint64_t exposed_mass_cell_count = 0;
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

class KrosaBaseStateProvider {
 public:
  KrosaBaseStateProvider() = default;

  [[nodiscard]] KrosaBaseStateProviderReport reconstruct(
      const Grid& grid,
      const io::KrosaPressureRefreshMetadata& metadata,
      KrosaMassBaseStateReconstructionOptions options = {});

  [[nodiscard]] KrosaBaseStateProviderReport reconstruct(
      const Grid& grid,
      const io::KrosaPressureRefreshMetadata& metadata,
      const KrosaBaseStateProviderTerrainOverride& terrain_override,
      KrosaMassBaseStateReconstructionOptions options = {});

  [[nodiscard]] KrosaBaseStateProviderViews views() const noexcept;

 private:
  [[nodiscard]] KrosaBaseStateProviderReport reconstruct_with_terrain(
      const Grid& grid,
      const io::KrosaPressureRefreshMetadata& metadata,
      FieldView2D<const float> terrain_height_m,
      bool terrain_override_used,
      std::string_view terrain_source_name,
      std::string_view terrain_provenance,
      KrosaMassBaseStateReconstructionOptions options);

  void clear_views() noexcept;

  FieldStorage3D<float> pb_;
  FieldStorage3D<float> t_init_;
  FieldStorage2D<float> mub_;
  FieldStorage3D<float> alb_;
  FieldStorage3D<float> phb_;
};

// Recomputes provider-owned PB/T_INIT/ALB only for child mass cells exposed by
// a moving-nest remap overlap window. The input MUB is assumed to be the
// already exposed-interpolated child-shaped MUB. This helper does not derive
// MUB from HGT, does not rebuild or synchronize PHB, does not touch State::t,
// does not write State ALB, and is not a gate/selected-field production path.
[[nodiscard]] KrosaExposedBaseStateRecomputeReport
recompute_exposed_base_state_from_mub(
    const nest::RemapWindow& overlap_window,
    const KrosaExposedBaseStateRecomputeInputs& inputs,
    const KrosaExposedBaseStateRecomputeOutputs& outputs,
    KrosaMassBaseStateReconstructionOptions options = {}) noexcept;

}  // namespace tywrf::dynamics
