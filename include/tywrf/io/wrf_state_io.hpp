#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/grid.hpp"
#include "tywrf/io/netcdf_schema.hpp"
#include "tywrf/state.hpp"

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace tywrf::io {

class WrfStateIoError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

struct WrfStateReadOptions {
  std::size_t time_index = 0;
  std::vector<std::string> variables;
};

struct WrfStateWriteOptions {
  std::size_t time_index = 0;
  std::vector<std::string> variables;
};

[[nodiscard]] std::vector<std::string> wrf_state_core_field_names();

[[nodiscard]] Grid derive_grid_from_wrf_schema(
    const DatasetSchema& schema,
    Halo3D halo = {});

[[nodiscard]] Grid derive_grid_from_wrf_file(
    const std::filesystem::path& path,
    Halo3D halo = {});

void load_wrf_state(
    const std::filesystem::path& path,
    State<float>& state,
    const WrfStateReadOptions& options = {});

void write_wrf_state(
    const std::filesystem::path& path,
    const State<float>& state,
    const WrfStateWriteOptions& options = {});

}  // namespace tywrf::io
