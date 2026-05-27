#include "tywrf/io/netcdf_schema.hpp"

#include <netcdf.h>

#include <array>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

std::filesystem::path reference_dir() {
  if (const char* env = std::getenv("TYWRF_REFERENCE_DIR")) {
    return env;
  }
  return "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF";
}

void print_issues(
    const std::string_view label,
    const std::vector<tywrf::io::SchemaIssue>& issues) {
  for (const auto& issue : issues) {
    std::cerr << label << ": " << issue.message << '\n';
  }
}

void check_nc(const int status, const std::string_view operation) {
  if (status == NC_NOERR) {
    return;
  }
  throw std::runtime_error(std::string(operation) + ": " + nc_strerror(status));
}

void create_synthetic_schema_file(const std::filesystem::path& path) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create synthetic netcdf");

  int time_dim = -1;
  int date_dim = -1;
  int west_east_dim = -1;
  check_nc(nc_def_dim(file_id, "Time", NC_UNLIMITED, &time_dim), "define Time");
  check_nc(nc_def_dim(file_id, "DateStrLen", 19, &date_dim), "define DateStrLen");
  check_nc(nc_def_dim(file_id, "west_east", 4, &west_east_dim), "define west_east");

  int times_var = -1;
  std::array<int, 2> times_dims = {time_dim, date_dim};
  check_nc(
      nc_def_var(file_id, "Times", NC_CHAR, 2, times_dims.data(), &times_var),
      "define Times");

  int sample_var = -1;
  std::array<int, 2> sample_dims = {time_dim, west_east_dim};
  check_nc(
      nc_def_var(file_id, "SAMPLE", NC_FLOAT, 2, sample_dims.data(), &sample_var),
      "define SAMPLE");

  check_nc(nc_close(file_id), "close synthetic netcdf");
}

bool contains_issue(
    const std::vector<tywrf::io::SchemaIssue>& issues,
    const std::string_view needle) {
  for (const auto& issue : issues) {
    if (issue.message.find(needle) != std::string::npos) {
      return true;
    }
  }
  return false;
}

int run_local_schema_reader_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_schema_reader_test.nc";
  std::filesystem::remove(path);
  create_synthetic_schema_file(path);

  const auto schema = tywrf::io::read_netcdf_schema(path, {"Times"});
  std::filesystem::remove(path);

  if (schema.variables.size() != 1) {
    std::cerr << "selected synthetic schema returned " << schema.variables.size()
              << " variables, expected 1\n";
    return 1;
  }
  const auto* time = tywrf::io::find_dimension(schema, "Time");
  if (time == nullptr || !time->unlimited || time->size != 0) {
    std::cerr << "synthetic Time dimension was not read as an empty unlimited dimension\n";
    return 1;
  }
  const auto* times = tywrf::io::find_variable(schema, "Times");
  if (times == nullptr || times->type_name != "char" ||
      times->dimensions != std::vector<std::string>{"Time", "DateStrLen"} ||
      times->shape != std::vector<std::size_t>{0, 19}) {
    std::cerr << "synthetic Times variable schema did not match expected metadata\n";
    return 1;
  }
  return 0;
}

int run_local_validator_negative_test() {
  tywrf::io::DatasetSchema schema;
  schema.dimensions = {
      {"Time", 1, false},
      {"DateStrLen", 19, false},
      {"west_east", 265, false},
      {"south_north", 429, false},
      {"bottom_top", 59, false},
      {"bottom_top_stag", 60, false},
      {"west_east_stag", 266, false},
      {"south_north_stag", 430, false},
  };
  schema.variables = {
      {"Times", "char", {"Time", "DateStrLen"}, {1, 19}},
      {"U", "double", {"Time", "bottom_top", "south_north", "west_east"}, {1, 59, 429, 265}},
  };

  const auto issues = tywrf::io::validate_wrfinput_d01_schema(schema);
  if (!tywrf::io::has_schema_errors(issues)) {
    std::cerr << "broken wrfinput_d01 schema did not produce errors\n";
    return 1;
  }
  if (!contains_issue(issues, "dimension Time is fixed, expected unlimited")) {
    std::cerr << "broken wrfinput_d01 schema did not report fixed Time dimension\n";
    return 1;
  }
  if (!contains_issue(issues, "variable U has type double, expected float")) {
    std::cerr << "broken wrfinput_d01 schema did not report U type mismatch\n";
    return 1;
  }
  if (!contains_issue(issues, "variable U has dimensions Time,bottom_top,south_north,west_east")) {
    std::cerr << "broken wrfinput_d01 schema did not report U dimension mismatch\n";
    return 1;
  }
  return 0;
}

}  // namespace

int main() {
  try {
    if (const int status = run_local_schema_reader_test(); status != 0) {
      return status;
    }
  } catch (const std::exception& error) {
    std::cerr << "local NetCDF schema reader test failed: " << error.what() << '\n';
    return 1;
  }

  if (const int status = run_local_validator_negative_test(); status != 0) {
    return status;
  }

  const auto root = reference_dir();
  const auto wrfinput = root / "wrfinput_d01";
  const auto wrfbdy = root / "wrfbdy_d01";
  const auto wrffdda = root / "wrffdda_d01";

  if (!std::filesystem::exists(wrfinput) || !std::filesystem::exists(wrfbdy) ||
      !std::filesystem::exists(wrffdda)) {
    std::cout << "Skipping NetCDF schema test; reference inputs are not present under "
              << root << '\n';
    return 0;
  }

  const auto input_schema = tywrf::io::read_netcdf_schema(wrfinput);
  const auto boundary_schema = tywrf::io::read_netcdf_schema(wrfbdy);
  const auto fdda_schema = tywrf::io::read_netcdf_schema(wrffdda);

  const auto input_issues = tywrf::io::validate_wrfinput_d01_schema(input_schema);
  const auto boundary_issues = tywrf::io::validate_wrfbdy_d01_schema(boundary_schema);
  const auto fdda_issues = tywrf::io::validate_wrffdda_d01_schema(fdda_schema);

  print_issues("wrfinput_d01", input_issues);
  print_issues("wrfbdy_d01", boundary_issues);
  print_issues("wrffdda_d01", fdda_issues);

  if (tywrf::io::has_schema_errors(input_issues) ||
      tywrf::io::has_schema_errors(boundary_issues) ||
      tywrf::io::has_schema_errors(fdda_issues)) {
    return 1;
  }

  const auto selected = tywrf::io::read_netcdf_schema(
      wrffdda,
      {"Times", "U_NDG_OLD", "Q_NDG_OLD", "MU_NDG_NEW"});
  if (selected.variables.size() != 4) {
    std::cerr << "selected variable schema returned " << selected.variables.size()
              << " variables, expected 4\n";
    return 1;
  }

  const auto* q_nudging = tywrf::io::find_variable(selected, "Q_NDG_OLD");
  if (q_nudging == nullptr || q_nudging->shape != std::vector<std::size_t>{28, 59, 429, 265}) {
    std::cerr << "Q_NDG_OLD schema shape did not match the KROSA reference file\n";
    return 1;
  }

  std::cout << "Validated KROSA NetCDF input schemas under " << root << '\n';
  return 0;
}
