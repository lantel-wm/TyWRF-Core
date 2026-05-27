#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace tywrf::io {

struct DimensionInfo {
  std::string name;
  std::size_t size = 0;
  bool unlimited = false;
};

struct VariableInfo {
  std::string name;
  std::string type_name;
  std::vector<std::string> dimensions;
  std::vector<std::size_t> shape;
};

struct DatasetSchema {
  std::filesystem::path path;
  std::vector<DimensionInfo> dimensions;
  std::vector<VariableInfo> variables;
};

enum class SchemaSeverity {
  error,
  warning,
};

struct SchemaIssue {
  SchemaSeverity severity = SchemaSeverity::error;
  std::string message;
};

[[nodiscard]] DatasetSchema read_netcdf_schema(
    const std::filesystem::path& path,
    const std::vector<std::string>& selected_variables = {});

[[nodiscard]] const DimensionInfo* find_dimension(
    const DatasetSchema& schema,
    std::string_view name) noexcept;

[[nodiscard]] const VariableInfo* find_variable(
    const DatasetSchema& schema,
    std::string_view name) noexcept;

[[nodiscard]] std::vector<SchemaIssue> validate_wrfinput_d01_schema(
    const DatasetSchema& schema);

[[nodiscard]] std::vector<SchemaIssue> validate_wrfbdy_d01_schema(
    const DatasetSchema& schema);

[[nodiscard]] std::vector<SchemaIssue> validate_wrffdda_d01_schema(
    const DatasetSchema& schema);

[[nodiscard]] bool has_schema_errors(const std::vector<SchemaIssue>& issues) noexcept;

}  // namespace tywrf::io
