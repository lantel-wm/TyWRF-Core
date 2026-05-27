#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/grid.hpp"
#include "tywrf/state.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace tywrf::io {

class PressureRefreshIoError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

enum class PressureRefreshPTopSource : std::uint8_t {
  missing,
  global_attribute,
  scalar_variable,
  time_variable,
};

struct KrosaPressureRefreshReadOptions {
  std::size_t time_index = 0;
};

struct KrosaPressureRefreshInputReport {
  std::filesystem::path source_path;
  std::size_t time_index = 0;
  std::vector<std::string> missing_names;
  std::vector<std::string> missing_base_state_reconstruction_names;

  std::size_t expected_mass_level_count = 0;
  std::size_t expected_full_level_count = 0;
  std::size_t c3f_count = 0;
  std::size_t c4f_count = 0;
  std::size_t c3h_count = 0;
  std::size_t c4h_count = 0;

  std::int32_t expected_alb_nx = 0;
  std::int32_t expected_alb_ny = 0;
  std::int32_t expected_alb_nz = 0;
  std::int32_t alb_nx = 0;
  std::int32_t alb_ny = 0;
  std::int32_t alb_nz = 0;
  std::size_t expected_alb_point_count = 0;
  std::size_t alb_point_count = 0;

  std::int32_t expected_terrain_nx = 0;
  std::int32_t expected_terrain_ny = 0;
  std::int32_t terrain_nx = 0;
  std::int32_t terrain_ny = 0;
  std::size_t expected_terrain_point_count = 0;
  std::size_t terrain_point_count = 0;
  std::string terrain_source_name;

  PressureRefreshPTopSource p_top_source = PressureRefreshPTopSource::missing;
  bool alb_loaded = false;
  bool terrain_loaded = false;

  [[nodiscard]] bool p_top_present() const noexcept {
    return p_top_source != PressureRefreshPTopSource::missing;
  }

  [[nodiscard]] bool ok() const noexcept {
    return missing_names.empty() && p_top_present() &&
           c3f_count == expected_full_level_count &&
           c4f_count == expected_full_level_count &&
           c3h_count == expected_mass_level_count &&
           c4h_count == expected_mass_level_count &&
           alb_point_count == expected_alb_point_count;
  }

  [[nodiscard]] bool base_state_reconstruction_inputs_ready() const noexcept {
    return missing_base_state_reconstruction_names.empty() && p_top_present() &&
           c3f_count == expected_full_level_count &&
           c4f_count == expected_full_level_count &&
           c3h_count == expected_mass_level_count &&
           c4h_count == expected_mass_level_count &&
           terrain_point_count == expected_terrain_point_count;
  }
};

struct KrosaPressureRefreshMetadata {
  std::filesystem::path source_path;
  std::size_t time_index = 0;
  float p_top_pa = 0.0F;
  PressureRefreshPTopSource p_top_source = PressureRefreshPTopSource::missing;
  std::vector<float> c3f;
  std::vector<float> c4f;
  std::vector<float> c3h;
  std::vector<float> c4h;
  std::string terrain_source_name;
  std::int32_t terrain_nx = 0;
  std::int32_t terrain_ny = 0;
  std::vector<float> terrain_height_m;
};

struct KrosaPressureRefreshReadResult {
  KrosaPressureRefreshMetadata metadata;
  KrosaPressureRefreshInputReport report;

  [[nodiscard]] bool ok() const noexcept {
    return report.ok();
  }

  [[nodiscard]] bool base_state_reconstruction_inputs_ready() const noexcept {
    return report.base_state_reconstruction_inputs_ready();
  }
};

[[nodiscard]] std::vector<std::string> krosa_pressure_refresh_required_names();

[[nodiscard]] KrosaPressureRefreshInputReport inspect_krosa_pressure_refresh_inputs(
    const std::filesystem::path& path,
    const Grid& grid,
    const KrosaPressureRefreshReadOptions& options = {});

[[nodiscard]] KrosaPressureRefreshReadResult read_krosa_pressure_refresh_inputs(
    const std::filesystem::path& path,
    const Grid& grid,
    FieldView3D<float> alb,
    const KrosaPressureRefreshReadOptions& options = {});

[[nodiscard]] KrosaPressureRefreshReadResult read_krosa_pressure_refresh_inputs(
    const std::filesystem::path& path,
    const Grid& grid,
    FieldStorage3D<float>& alb,
    const KrosaPressureRefreshReadOptions& options = {});

}  // namespace tywrf::io
