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

}  // namespace tywrf::dynamics
