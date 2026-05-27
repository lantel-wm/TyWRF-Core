#pragma once

#include "tywrf/io/netcdf_schema.hpp"

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace tywrf::io {

class ForcingIoError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

enum class KrosaForcingKind {
  boundary,
  spectral_nudging,
};

struct ForcingVariableMetadata {
  std::string name;
  std::string type_name;
  std::vector<std::string> dimensions;
  std::vector<std::size_t> shape;
  std::vector<std::string> slice_dimensions;
  std::vector<std::size_t> slice_shape;
  std::size_t time_count = 0;
  std::size_t values_per_time_slice = 0;
};

struct ForcingTimeSlice {
  ForcingVariableMetadata metadata;
  std::size_t time_index = 0;
  std::vector<float> values;
};

[[nodiscard]] std::vector<std::string> krosa_wrfbdy_required_variable_names();
[[nodiscard]] std::vector<std::string> krosa_wrffdda_required_variable_names();

class KrosaForcingReader {
 public:
  KrosaForcingReader(std::filesystem::path path, KrosaForcingKind kind);

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] KrosaForcingKind kind() const noexcept;
  [[nodiscard]] const DatasetSchema& schema() const noexcept;
  [[nodiscard]] std::size_t time_count() const noexcept;
  [[nodiscard]] std::size_t date_string_length() const noexcept;

  [[nodiscard]] bool has_variable(std::string_view name) const noexcept;
  [[nodiscard]] std::vector<std::string> missing_required_variables() const;
  [[nodiscard]] ForcingVariableMetadata variable_metadata(std::string_view name) const;
  [[nodiscard]] std::string read_char_time_slice(
      std::string_view name,
      std::size_t time_index) const;
  [[nodiscard]] std::string read_time_string(std::size_t time_index) const;
  [[nodiscard]] std::vector<std::string> read_time_strings() const;
  [[nodiscard]] ForcingTimeSlice read_float_time_slice(
      std::string_view name,
      std::size_t time_index) const;

 private:
  std::filesystem::path path_;
  KrosaForcingKind kind_ = KrosaForcingKind::boundary;
  DatasetSchema schema_;
  std::size_t time_count_ = 0;
  std::size_t date_string_length_ = 0;
};

}  // namespace tywrf::io
