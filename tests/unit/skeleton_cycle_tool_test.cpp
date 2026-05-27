#include <netcdf.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

constexpr std::string_view kCycleStart = "2025-07-26_00:00:00";
constexpr std::string_view kTenMinuteEnd = "2025-07-26_00:10:00";
constexpr std::string_view kPressureRefreshFormulaStatus =
    "staged_helper_available_not_applied";
constexpr std::string_view kPressureRefreshFormulaStagingName =
    "krosa_hypsometric_opt2_use_theta_m";
constexpr std::string_view kPressureRefreshRegionStagingName = "exposed_mass_cells";
constexpr std::string_view kPressureRefreshThermodynamicMode = "use_theta_m_1";
constexpr std::string_view kPressureRefreshIntegrationStatus = "helper_available_not_invoked";
constexpr std::string_view kPressureRefreshOptInIntegrationStatus =
    "hook_compute_invoked_diagnostic_non_gate";
constexpr std::string_view kPressureRefreshRequiredInputsCsv =
    "ALB,C3F,C4F,C3H,C4H,P_TOP,MU,MUB,T,PB,PH,PHB";
constexpr std::string_view kPressureRefreshRequiredInputsJson =
    "[\"ALB\", \"C3F\", \"C4F\", \"C3H\", \"C4H\", \"P_TOP\", \"MU\", \"MUB\", "
    "\"T\", \"PB\", \"PH\", \"PHB\"]";
const std::filesystem::path kKrosaTenMinuteReferenceDir =
    "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF_1h_10min_20260527_172838";

void check_nc(const int status, const std::string_view operation) {
  if (status == NC_NOERR) {
    return;
  }
  throw std::runtime_error(std::string(operation) + ": " + nc_strerror(status));
}

float value_3d(
    const int time_index,
    const int k,
    const int j,
    const int i,
    const float offset) {
  return offset + 1000.0F * static_cast<float>(time_index) +
         100.0F * static_cast<float>(k) + 10.0F * static_cast<float>(j) +
         static_cast<float>(i);
}

float value_2d(const int time_index, const int j, const int i, const float offset) {
  return offset + 1000.0F * static_cast<float>(time_index) +
         10.0F * static_cast<float>(j) + static_cast<float>(i);
}

void put_times(
    const int file_id,
    const int variable_id,
    const std::string_view value) {
  constexpr std::size_t date_str_len = 19;
  std::vector<char> buffer(date_str_len, ' ');
  std::copy(value.begin(), value.end(), buffer.begin());
  const std::size_t start[2] = {0, 0};
  const std::size_t count[2] = {1, date_str_len};
  check_nc(nc_put_vara_text(file_id, variable_id, start, count, buffer.data()), "write Times");
}

void put_3d(
    const int file_id,
    const int variable_id,
    const int nz,
    const int ny,
    const int nx,
    const float offset) {
  std::vector<float> values(static_cast<std::size_t>(nz * ny * nx));
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        values[(static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
                static_cast<std::size_t>(j)) *
                   static_cast<std::size_t>(nx) +
               static_cast<std::size_t>(i)] = value_3d(0, k, j, i, offset);
      }
    }
  }
  const std::size_t start[4] = {0, 0, 0, 0};
  const std::size_t count[4] = {
      1,
      static_cast<std::size_t>(nz),
      static_cast<std::size_t>(ny),
      static_cast<std::size_t>(nx)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, count, values.data()), "write 3D var");
}

void put_2d(
    const int file_id,
    const int variable_id,
    const int ny,
    const int nx,
    const float offset) {
  std::vector<float> values(static_cast<std::size_t>(ny * nx));
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      values[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
             static_cast<std::size_t>(i)] = value_2d(0, j, i, offset);
    }
  }
  const std::size_t start[3] = {0, 0, 0};
  const std::size_t count[3] = {1, static_cast<std::size_t>(ny), static_cast<std::size_t>(nx)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, count, values.data()), "write 2D var");
}

void put_vector(
    const int file_id,
    const int variable_id,
    const std::vector<float>& values) {
  const std::size_t start[2] = {0, 0};
  const std::size_t count[2] = {1, values.size()};
  check_nc(
      nc_put_vara_float(file_id, variable_id, start, count, values.data()),
      "write vector var");
}

void put_alb(
    const int file_id,
    const int variable_id,
    const int nz,
    const int ny,
    const int nx) {
  std::vector<float> values(static_cast<std::size_t>(nz * ny * nx), 0.0F);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        values[(static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
                static_cast<std::size_t>(j)) *
                   static_cast<std::size_t>(nx) +
               static_cast<std::size_t>(i)] =
            0.8F + 0.01F * static_cast<float>(k + j + i);
      }
    }
  }
  const std::size_t start[4] = {0, 0, 0, 0};
  const std::size_t count[4] = {
      1,
      static_cast<std::size_t>(nz),
      static_cast<std::size_t>(ny),
      static_cast<std::size_t>(nx)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, count, values.data()), "write ALB");
}

void create_synthetic_wrfout(
    const std::filesystem::path& path,
    const double dx,
    const std::string_view time_value,
    const bool state_fields,
    const float xlat_offset = 40000.0F,
    const float xlong_offset = 50000.0F,
    const float hgt_offset = 60000.0F,
    const int mass_nx = 4,
    const int mass_ny = 3,
    const int mass_nz = 2,
    const int full_nz = 3,
    const float u_offset = 10000.0F,
    const float t_offset = 20000.0F,
    const float mu_offset = 30000.0F,
    const float p_offset = 40000.0F,
    const float pb_offset = 80000.0F,
    const float ph_offset = 1000.0F,
    const float phb_offset = 2000.0F,
    const float mub_offset = 90000.0F,
    const bool pressure_refresh_metadata = false,
    const bool include_direct_alb = false) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create synthetic wrfout");
  check_nc(nc_put_att_double(file_id, NC_GLOBAL, "DX", NC_DOUBLE, 1, &dx), "write DX");
  check_nc(nc_put_att_double(file_id, NC_GLOBAL, "DY", NC_DOUBLE, 1, &dx), "write DY");
  if (pressure_refresh_metadata) {
    const float p_top = 5000.0F;
    check_nc(nc_put_att_float(file_id, NC_GLOBAL, "P_TOP", NC_FLOAT, 1, &p_top), "write P_TOP");
  }

  int time_dim = -1;
  int date_str_len_dim = -1;
  int bottom_top_dim = -1;
  int bottom_top_stag_dim = -1;
  int south_north_dim = -1;
  int south_north_stag_dim = -1;
  int west_east_dim = -1;
  int west_east_stag_dim = -1;
  check_nc(nc_def_dim(file_id, "Time", NC_UNLIMITED, &time_dim), "define Time");
  check_nc(nc_def_dim(file_id, "DateStrLen", 19, &date_str_len_dim), "define DateStrLen");
  check_nc(nc_def_dim(file_id, "bottom_top", mass_nz, &bottom_top_dim), "define bottom_top");
  check_nc(
      nc_def_dim(file_id, "bottom_top_stag", full_nz, &bottom_top_stag_dim),
      "define bottom_top_stag");
  check_nc(nc_def_dim(file_id, "south_north", mass_ny, &south_north_dim), "define south_north");
  check_nc(
      nc_def_dim(file_id, "south_north_stag", mass_ny + 1, &south_north_stag_dim),
      "define south_north_stag");
  check_nc(nc_def_dim(file_id, "west_east", mass_nx, &west_east_dim), "define west_east");
  check_nc(
      nc_def_dim(file_id, "west_east_stag", mass_nx + 1, &west_east_stag_dim),
      "define west_east_stag");

  int times_var = -1;
  int xlat_var = -1;
  int xlong_var = -1;
  int hgt_var = -1;
  int u_var = -1;
  int t_var = -1;
  int mu_var = -1;
  int p_var = -1;
  int pb_var = -1;
  int ph_var = -1;
  int phb_var = -1;
  int mub_var = -1;
  int c3f_var = -1;
  int c4f_var = -1;
  int c3h_var = -1;
  int c4h_var = -1;
  int alb_var = -1;
  const int times_dims[2] = {time_dim, date_str_len_dim};
  const int static_dims[3] = {time_dim, south_north_dim, west_east_dim};
  const int u_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_stag_dim};
  const int w_dims[4] = {time_dim, bottom_top_stag_dim, south_north_dim, west_east_dim};
  const int t_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_dim};
  const int mu_dims[3] = {time_dim, south_north_dim, west_east_dim};
  const int full_coeff_dims[2] = {time_dim, bottom_top_stag_dim};
  const int mass_coeff_dims[2] = {time_dim, bottom_top_dim};
  check_nc(nc_def_var(file_id, "Times", NC_CHAR, 2, times_dims, &times_var), "define Times");
  check_nc(nc_def_var(file_id, "XLAT", NC_FLOAT, 3, static_dims, &xlat_var), "define XLAT");
  check_nc(nc_def_var(file_id, "XLONG", NC_FLOAT, 3, static_dims, &xlong_var), "define XLONG");
  check_nc(nc_def_var(file_id, "HGT", NC_FLOAT, 3, static_dims, &hgt_var), "define HGT");
  if (state_fields) {
    check_nc(nc_def_var(file_id, "U", NC_FLOAT, 4, u_dims, &u_var), "define U");
    check_nc(nc_def_var(file_id, "T", NC_FLOAT, 4, t_dims, &t_var), "define T");
    check_nc(nc_def_var(file_id, "MU", NC_FLOAT, 3, mu_dims, &mu_var), "define MU");
    check_nc(nc_def_var(file_id, "P", NC_FLOAT, 4, t_dims, &p_var), "define P");
    check_nc(nc_def_var(file_id, "PB", NC_FLOAT, 4, t_dims, &pb_var), "define PB");
    check_nc(nc_def_var(file_id, "PH", NC_FLOAT, 4, w_dims, &ph_var), "define PH");
    check_nc(nc_def_var(file_id, "PHB", NC_FLOAT, 4, w_dims, &phb_var), "define PHB");
    check_nc(nc_def_var(file_id, "MUB", NC_FLOAT, 3, mu_dims, &mub_var), "define MUB");
  }
  if (pressure_refresh_metadata) {
    check_nc(nc_def_var(file_id, "C3F", NC_FLOAT, 2, full_coeff_dims, &c3f_var), "define C3F");
    check_nc(nc_def_var(file_id, "C4F", NC_FLOAT, 2, full_coeff_dims, &c4f_var), "define C4F");
    check_nc(nc_def_var(file_id, "C3H", NC_FLOAT, 2, mass_coeff_dims, &c3h_var), "define C3H");
    check_nc(nc_def_var(file_id, "C4H", NC_FLOAT, 2, mass_coeff_dims, &c4h_var), "define C4H");
    if (include_direct_alb) {
      check_nc(nc_def_var(file_id, "ALB", NC_FLOAT, 4, t_dims, &alb_var), "define ALB");
    }
  }
  check_nc(nc_enddef(file_id), "end definitions");

  put_times(file_id, times_var, time_value);
  put_2d(file_id, xlat_var, mass_ny, mass_nx, xlat_offset);
  put_2d(file_id, xlong_var, mass_ny, mass_nx, xlong_offset);
  put_2d(file_id, hgt_var, mass_ny, mass_nx, hgt_offset);
  if (state_fields) {
    put_3d(file_id, u_var, mass_nz, mass_ny, mass_nx + 1, u_offset);
    put_3d(file_id, t_var, mass_nz, mass_ny, mass_nx, t_offset);
    put_2d(file_id, mu_var, mass_ny, mass_nx, mu_offset);
    put_3d(file_id, p_var, mass_nz, mass_ny, mass_nx, p_offset);
    put_3d(file_id, pb_var, mass_nz, mass_ny, mass_nx, pb_offset);
    put_3d(file_id, ph_var, full_nz, mass_ny, mass_nx, ph_offset);
    put_3d(file_id, phb_var, full_nz, mass_ny, mass_nx, phb_offset);
    put_2d(file_id, mub_var, mass_ny, mass_nx, mub_offset);
  }
  if (pressure_refresh_metadata) {
    put_vector(file_id, c3f_var, {1.0F, 0.5F, 0.0F});
    put_vector(file_id, c4f_var, {0.0F, 0.0F, 0.0F});
    put_vector(file_id, c3h_var, {0.75F, 0.25F});
    put_vector(file_id, c4h_var, {0.0F, 0.0F});
    if (include_direct_alb) {
      put_alb(file_id, alb_var, mass_nz, mass_ny, mass_nx);
    }
  }

  check_nc(nc_close(file_id), "close synthetic wrfout");
}

[[nodiscard]] std::string shell_quote(const std::filesystem::path& path) {
  return "'" + path.string() + "'";
}

[[nodiscard]] std::string shell_quote(const std::string& value) {
  return "'" + value + "'";
}

void run_command(const std::string& command) {
  const int status = std::system(command.c_str());
  if (status != 0) {
    throw std::runtime_error("command failed: " + command);
  }
}

void run_command_expect_failure(const std::string& command) {
  const int status = std::system(command.c_str());
  if (status == 0) {
    throw std::runtime_error("command unexpectedly passed: " + command);
  }
}

void run_command_capture_expect_failure(
    const std::string& command,
    const std::filesystem::path& output) {
  const int status = std::system((command + " > " + shell_quote(output)).c_str());
  if (status == 0) {
    throw std::runtime_error("command unexpectedly passed: " + command);
  }
}

[[nodiscard]] std::string read_text_attr(
    const int file_id,
    const std::string_view name) {
  std::size_t length = 0;
  check_nc(nc_inq_attlen(file_id, NC_GLOBAL, std::string(name).c_str(), &length), "read attr len");
  std::string value(length, '\0');
  check_nc(nc_get_att_text(file_id, NC_GLOBAL, std::string(name).c_str(), value.data()), "read attr");
  return value;
}

[[nodiscard]] double read_double_attr(
    const int file_id,
    const std::string_view name) {
  double value = 0.0;
  check_nc(nc_get_att_double(file_id, NC_GLOBAL, std::string(name).c_str(), &value), "read attr");
  return value;
}

[[nodiscard]] bool has_global_attr(
    const int file_id,
    const std::string_view name) {
  const int status = nc_inq_att(file_id, NC_GLOBAL, std::string(name).c_str(), nullptr, nullptr);
  if (status == NC_NOERR) {
    return true;
  }
  if (status == NC_ENOTATT) {
    return false;
  }
  check_nc(status, "inquire attr");
  return false;
}

[[nodiscard]] std::string read_times(const int file_id) {
  int variable_id = -1;
  check_nc(nc_inq_varid(file_id, "Times", &variable_id), "inquire Times");
  std::vector<char> buffer(19, '\0');
  const std::size_t start[2] = {0, 0};
  const std::size_t count[2] = {1, 19};
  check_nc(nc_get_vara_text(file_id, variable_id, start, count, buffer.data()), "read Times");
  return std::string(buffer.begin(), buffer.end());
}

[[nodiscard]] float read_3d_value(
    const int file_id,
    const std::string_view name,
    const int k,
    const int j,
    const int i) {
  int variable_id = -1;
  check_nc(nc_inq_varid(file_id, std::string(name).c_str(), &variable_id), "inquire 3D var");
  float value = 0.0F;
  const std::size_t start[4] = {0, static_cast<std::size_t>(k), static_cast<std::size_t>(j), static_cast<std::size_t>(i)};
  const std::size_t count[4] = {1, 1, 1, 1};
  check_nc(nc_get_vara_float(file_id, variable_id, start, count, &value), "read 3D var");
  return value;
}

[[nodiscard]] float read_2d_value(
    const int file_id,
    const std::string_view name,
    const int j,
    const int i) {
  int variable_id = -1;
  check_nc(nc_inq_varid(file_id, std::string(name).c_str(), &variable_id), "inquire 2D var");
  float value = 0.0F;
  const std::size_t start[3] = {0, static_cast<std::size_t>(j), static_cast<std::size_t>(i)};
  const std::size_t count[3] = {1, 1, 1};
  check_nc(nc_get_vara_float(file_id, variable_id, start, count, &value), "read 2D var");
  return value;
}

void assert_parent_fill_pressure_refresh_attrs(const int file_id) {
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_REQUIRED") == "true");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_APPLIED") == "false");
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_REQUIREMENT_STATUS") ==
      "required_for_exposed_parent_fill_cells");
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS") ==
      std::string(kPressureRefreshIntegrationStatus));
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_FORMULA_STATUS") ==
      std::string(kPressureRefreshFormulaStatus));
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_FORMULA_STAGING_NAME") ==
      std::string(kPressureRefreshFormulaStagingName));
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_REGION_STAGING_NAME") ==
      std::string(kPressureRefreshRegionStagingName));
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_THERMODYNAMIC_MODE") ==
      std::string(kPressureRefreshThermodynamicMode));
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_REQUIRED_INPUTS") ==
      std::string(kPressureRefreshRequiredInputsCsv));
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_OUTPUT_FIELD") == "P");
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_HELPER_NAME") ==
      "refresh_krosa_moving_nest_pressure");
}

void assert_parent_fill_pressure_refresh_report(const std::string& report) {
  assert(report.find("\"pressure_refresh_required\": true") != std::string::npos);
  assert(report.find("\"pressure_refresh_applied\": false") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_requirement_status\": "
          "\"required_for_exposed_parent_fill_cells\"") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_integration_status\": \"" +
          std::string(kPressureRefreshIntegrationStatus) + "\"") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_formula_status\": \"" +
          std::string(kPressureRefreshFormulaStatus) + "\"") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_formula_staging_name\": \"" +
          std::string(kPressureRefreshFormulaStagingName) + "\"") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_region_staging_name\": \"" +
          std::string(kPressureRefreshRegionStagingName) + "\"") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_thermodynamic_mode\": \"" +
          std::string(kPressureRefreshThermodynamicMode) + "\"") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_required_inputs\": " +
          std::string(kPressureRefreshRequiredInputsJson)) != std::string::npos);
  assert(report.find("\"pressure_refresh_output_field\": \"P\"") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_helper_name\": \"refresh_krosa_moving_nest_pressure\"") !=
      std::string::npos);
}

void assert_parent_fill_pressure_refresh_opt_in_attrs(
    const int file_id,
    const std::filesystem::path& metadata_source) {
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_REQUIRED") == "true");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_APPLIED") == "true");
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS") ==
      std::string(kPressureRefreshOptInIntegrationStatus));
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_ALB_SOURCE") ==
      "base_state_reconstruction_provider");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_PROVIDER_OK") == "true");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_STAGING_OK") == "true");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED") == "true");
  assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS") > 0.0);
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE") == metadata_source.string());
  assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_METADATA_TIME_INDEX") == 0.0);
}

void assert_parent_fill_pressure_refresh_opt_in_report(
    const std::string& report,
    const std::filesystem::path& metadata_source) {
  assert(report.find("\"pressure_refresh_required\": true") != std::string::npos);
  assert(report.find("\"pressure_refresh_applied\": true") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_integration_status\": \"" +
          std::string(kPressureRefreshOptInIntegrationStatus) + "\"") != std::string::npos);
  assert(
      report.find(
          "\"pressure_refresh_alb_source\": "
          "\"base_state_reconstruction_provider\"") != std::string::npos);
  assert(report.find("\"pressure_refresh_provider_ok\": true") != std::string::npos);
  assert(report.find("\"pressure_refresh_staging_ok\": true") != std::string::npos);
  assert(report.find("\"pressure_refresh_compute_called\": true") != std::string::npos);
  assert(report.find("\"pressure_refresh_synced_pb_points\":") != std::string::npos);
  assert(report.find("\"pressure_refresh_synced_mub_points\":") != std::string::npos);
  assert(report.find("\"pressure_refresh_synced_phb_points\":") != std::string::npos);
  assert(report.find("\"pressure_refresh_refreshed_p_points\":") != std::string::npos);
  assert(report.find("\"pressure_refresh_metadata_source\": \"" + metadata_source.string()) !=
         std::string::npos);
  assert(report.find("\"pressure_refresh_metadata_time_index\": 0") != std::string::npos);
}

void assert_output(const std::filesystem::path& output) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open output");
  assert(read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") == "cpp_skeleton_candidate");
  assert(read_text_attr(file_id, "TYWRF_SKELETON") == "true");
  assert(read_text_attr(file_id, "TYWRF_NOT_PHYSICAL") == "true");
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") == "false");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "true");
  assert(read_text_attr(file_id, "TYWRF_D02_RESOLUTION_CHECK") == "d02_2km");
  assert(!has_global_attr(file_id, "TYWRF_DIAGNOSTIC_REMAP_OVERLAP"));
  assert(read_text_attr(file_id, "TYWRF_STATIC_SOURCE").find("wrfout_d02_2025-07-26_00:00:00") !=
         std::string::npos);
  assert(read_text_attr(file_id, "TYWRF_TIMES_SOURCE") == "--times:2025-07-26_00:10:00");
  assert(read_text_attr(file_id, "TYWRF_STATIC_COORDS_MATCH_STATE_SOURCE") == "same_file");
  assert(
      read_text_attr(file_id, "TYWRF_STATE_STATIC_CONSISTENCY") ==
      "static_coords_same_file_as_state_source");
  assert(read_double_attr(file_id, "TYWRF_SUGGESTED_GATE_INTERVAL_MINUTES") == 10.0);
  assert(read_text_attr(file_id, "TYWRF_SUGGESTED_GATE_INTERVAL_OPTION") == "--interval-minutes 10");
  assert(read_double_attr(file_id, "DX") == 2000.0);
  assert(read_double_attr(file_id, "DY") == 2000.0);
  assert(read_times(file_id).substr(0, 19) == kTenMinuteEnd);
  assert(read_3d_value(file_id, "U", 1, 2, 4) == value_3d(0, 1, 2, 4, 10000.0F));
  assert(read_3d_value(file_id, "T", 1, 2, 3) == value_3d(0, 1, 2, 3, 20000.0F));
  assert(read_2d_value(file_id, "MU", 2, 3) == value_2d(0, 2, 3, 30000.0F));
  assert(read_2d_value(file_id, "XLAT", 2, 3) == value_2d(0, 2, 3, 40000.0F));
  assert(read_2d_value(file_id, "XLAT", 2, 3) != value_2d(0, 2, 3, 70000.0F));
  assert(read_2d_value(file_id, "XLONG", 1, 2) == value_2d(0, 1, 2, 50000.0F));
  assert(read_2d_value(file_id, "XLONG", 1, 2) != value_2d(0, 1, 2, 80000.0F));
  assert(read_2d_value(file_id, "HGT", 0, 1) == value_2d(0, 0, 1, 60000.0F));
  assert(read_2d_value(file_id, "HGT", 0, 1) != value_2d(0, 0, 1, 90000.0F));
  check_nc(nc_close(file_id), "close output");
}

void assert_remap_output(
    const std::filesystem::path& output,
    const std::filesystem::path& report_path) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open remap output");
  assert(read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") == "cpp_skeleton_remap_overlap_diagnostic");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_REMAP_OVERLAP") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ONLY") == "true");
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") == "false");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "false");
  assert(read_text_attr(file_id, "TYWRF_NOT_PHYSICAL") == "true");
  assert(read_text_attr(file_id, "TYWRF_NEEDS_PARENT_FILL") == "true");
  assert(read_text_attr(file_id, "TYWRF_REMAP_EXPOSED_CELLS_FILLED_BY_PARENT") == "false");
  assert(read_text_attr(file_id, "TYWRF_UNFILLED_EXPOSED_CELLS") == "nan_sentinel_pending_parent_fill");
  assert(read_text_attr(file_id, "TYWRF_REMAP_FROM_PARENT_START") == "114,96");
  assert(read_text_attr(file_id, "TYWRF_REMAP_TO_PARENT_START") == "115,96");
  assert(read_double_attr(file_id, "TYWRF_REMAP_COPIED_FIELD_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_REMAP_COPIED_POINT_COUNT") > 0.0);
  assert(read_3d_value(file_id, "U", 0, 0, 0) == value_3d(0, 0, 0, 5, 10000.0F));
  assert(std::isnan(read_3d_value(file_id, "U", 0, 0, 210)));
  check_nc(nc_close(file_id), "close remap output");

  std::ifstream report_file(report_path);
  std::ostringstream report_buffer;
  report_buffer << report_file.rdbuf();
  const auto report = report_buffer.str();
  assert(report.find("\"diagnostic_remap_overlap\": true") != std::string::npos);
  assert(report.find("\"needs_parent_fill\": true") != std::string::npos);
  assert(report.find("\"gate_candidate\": false") != std::string::npos);
  assert(report.find("\"remap_copied_point_count\":") != std::string::npos);
}

void assert_parent_fill_output(
    const std::filesystem::path& output,
    const std::filesystem::path& parent_fill_state,
    const std::filesystem::path& report_path) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open parent-fill output");
  assert(
      read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") ==
      "cpp_skeleton_remap_parent_fill_diagnostic");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_REMAP_OVERLAP") == "false");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_REMAP_PARENT_FILL") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ONLY") == "true");
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") == "false");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "false");
  assert(read_text_attr(file_id, "TYWRF_EXPECTED_TO_MEET_THRESHOLDS") == "false");
  assert(read_text_attr(file_id, "TYWRF_NOT_PHYSICAL") == "true");
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") == "false");
  assert(read_text_attr(file_id, "TYWRF_MINIMUM_STATIC_REFRESH_FIELDS") == "XLAT,XLONG,HGT");
  assert(
      read_text_attr(file_id, "TYWRF_STAGGERED_STATIC_COORDS_STATUS") ==
      "pending_unless_emitted_later");
  assert(
      read_text_attr(file_id, "TYWRF_P_DERIVED_REFRESH_STATUS") ==
      "pending_derive_or_recompute_after_parent_fill_not_direct_wrf_parent_fill");
  assert_parent_fill_pressure_refresh_attrs(file_id);
  assert(
      read_text_attr(file_id, "TYWRF_DIRECT_WRF_END_STATE_ORACLE_STATUS") ==
      "diagnostic_only_nonphysical_non_gate");
  assert(read_text_attr(file_id, "TYWRF_NEEDS_PARENT_FILL") == "false");
  assert(read_text_attr(file_id, "TYWRF_REMAP_EXPOSED_CELLS_FILLED_BY_PARENT") == "true");
  assert(read_text_attr(file_id, "TYWRF_UNFILLED_EXPOSED_CELLS") == "none_parent_fill_diagnostic");
  assert(read_text_attr(file_id, "TYWRF_REMAP_FROM_PARENT_START") == "114,96");
  assert(read_text_attr(file_id, "TYWRF_REMAP_TO_PARENT_START") == "115,96");
  assert(read_text_attr(file_id, "TYWRF_PARENT_FILL_SOURCE") == parent_fill_state.string());
  assert(read_double_attr(file_id, "TYWRF_PARENT_FILL_TIME_INDEX") == 0.0);
  assert(read_double_attr(file_id, "TYWRF_REMAP_COPIED_FIELD_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_REMAP_COPIED_POINT_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_REMAP_PARENT_FILL_FIELD_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_REMAP_PARENT_FILL_POINT_COUNT") > 0.0);
  assert(read_3d_value(file_id, "U", 0, 0, 0) == value_3d(0, 0, 0, 5, 10000.0F));
  assert(read_3d_value(file_id, "U", 0, 0, 210) == value_3d(0, 0, 0, 210, 110000.0F));
  assert(read_3d_value(file_id, "T", 0, 0, 0) == value_3d(0, 0, 0, 5, 20000.0F));
  assert(read_3d_value(file_id, "T", 0, 0, 209) == value_3d(0, 0, 0, 209, 120000.0F));
  assert(read_2d_value(file_id, "MU", 0, 0) == value_2d(0, 0, 5, 30000.0F));
  assert(read_2d_value(file_id, "MU", 0, 209) == value_2d(0, 0, 209, 130000.0F));
  assert(read_2d_value(file_id, "XLAT", 2, 3) == value_2d(0, 2, 3, 140000.0F));
  assert(read_2d_value(file_id, "XLONG", 1, 2) == value_2d(0, 1, 2, 150000.0F));
  assert(read_2d_value(file_id, "HGT", 0, 1) == value_2d(0, 0, 1, 160000.0F));
  check_nc(nc_close(file_id), "close parent-fill output");

  std::ifstream report_file(report_path);
  std::ostringstream report_buffer;
  report_buffer << report_file.rdbuf();
  const auto report = report_buffer.str();
  assert(report.find("\"diagnostic_remap_overlap\": false") != std::string::npos);
  assert(report.find("\"diagnostic_remap_parent_fill\": true") != std::string::npos);
  assert(report.find("\"needs_parent_fill\": false") != std::string::npos);
  assert(report.find("\"exposed_cells_filled_by_parent\": true") != std::string::npos);
  assert(report.find("\"gate_candidate\": false") != std::string::npos);
  assert(report.find("\"integrator_output\": false") != std::string::npos);
  assert(report.find("\"minimum_static_refresh_fields\": [\"XLAT\", \"XLONG\", \"HGT\"]") !=
         std::string::npos);
  assert(report.find("\"staggered_static_coords_status\": \"pending_unless_emitted_later\"") !=
         std::string::npos);
  assert(
      report.find(
          "\"p_derived_refresh_status\": "
          "\"pending_derive_or_recompute_after_parent_fill_not_direct_wrf_parent_fill\"") !=
      std::string::npos);
  assert_parent_fill_pressure_refresh_report(report);
  assert(
      report.find(
          "\"direct_wrf_end_state_oracle_status\": "
          "\"diagnostic_only_nonphysical_non_gate\"") != std::string::npos);
  assert(report.find("\"parent_fill_source\": \"" + parent_fill_state.string()) !=
         std::string::npos);
  assert(report.find("\"remap_parent_fill_point_count\":") != std::string::npos);
}

void assert_parent_fill_pressure_refresh_output(
    const std::filesystem::path& output,
    const std::filesystem::path& parent_fill_state,
    const std::filesystem::path& metadata_source,
    const std::filesystem::path& report_path) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open pressure-refresh output");
  assert(
      read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") ==
      "cpp_skeleton_remap_parent_fill_diagnostic");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_REMAP_PARENT_FILL") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ONLY") == "true");
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") == "false");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "false");
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") == "false");
  assert(read_text_attr(file_id, "TYWRF_PARENT_FILL_SOURCE") == parent_fill_state.string());
  assert_parent_fill_pressure_refresh_opt_in_attrs(file_id, metadata_source);

  assert(read_3d_value(file_id, "P", 0, 0, 0) == value_3d(0, 0, 0, 5, 40000.0F));
  const auto exposed_parent_fill_p = value_3d(0, 0, 0, 209, 240000.0F);
  const auto exposed_refreshed_p = read_3d_value(file_id, "P", 0, 0, 209);
  assert(std::isfinite(exposed_refreshed_p));
  assert(exposed_refreshed_p != exposed_parent_fill_p);
  assert(read_3d_value(file_id, "PB", 0, 0, 209) != value_3d(0, 0, 0, 209, 280000.0F));
  assert(read_2d_value(file_id, "MUB", 0, 209) != value_2d(0, 0, 209, 290000.0F));
  check_nc(nc_close(file_id), "close pressure-refresh output");

  std::ifstream report_file(report_path);
  std::ostringstream report_buffer;
  report_buffer << report_file.rdbuf();
  const auto report = report_buffer.str();
  assert(report.find("\"gate_candidate\": false") != std::string::npos);
  assert(report.find("\"integrator_output\": false") != std::string::npos);
  assert_parent_fill_pressure_refresh_opt_in_report(report, metadata_source);
}

std::string base_command(
    const std::filesystem::path& executable,
    const std::filesystem::path& state,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output) {
  return shell_quote(executable) + " --state " + shell_quote(state) + " --template " +
         shell_quote(template_path) + " --output " + shell_quote(output) +
         " --cycle-start " + shell_quote(std::string(kCycleStart)) +
         " --cycle-end " + shell_quote(std::string(kTenMinuteEnd)) +
         " --times " + shell_quote(std::string(kTenMinuteEnd)) +
         " --variables Times,XLAT,XLONG,HGT,U,T,MU --pretty";
}

void run_real_krosa_smoke_if_available(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto state = kKrosaTenMinuteReferenceDir / "wrfout_d02_2025-07-26_00:00:00";
  const auto reference_end = kKrosaTenMinuteReferenceDir / "wrfout_d02_2025-07-26_00:10:00";
  if (!std::filesystem::exists(state) || !std::filesystem::exists(reference_end)) {
    std::cerr << "skipping real KROSA 10 min smoke; reference files are not present\n";
    return;
  }

  const auto candidate_dir = root / "real_krosa_candidate";
  const auto output = candidate_dir / "wrfout_d02_2025-07-26_00:10:00";
  const auto skeleton_report = root / "real_krosa_skeleton_report.json";
  const auto gate_report = root / "real_krosa_strict_gate_report.json";
  std::filesystem::create_directories(candidate_dir);
  std::filesystem::remove(output);
  std::filesystem::remove(skeleton_report);
  std::filesystem::remove(gate_report);

  run_command(
      base_command(executable, state, state, output) + " > " + shell_quote(skeleton_report));

  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open real KROSA output");
  assert(read_text_attr(file_id, "TYWRF_NOT_PHYSICAL") == "true");
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") == "false");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "true");
  assert(read_double_attr(file_id, "DX") == 2000.0);
  assert(read_double_attr(file_id, "DY") == 2000.0);
  assert(read_double_attr(file_id, "TYWRF_SUGGESTED_GATE_INTERVAL_MINUTES") == 10.0);
  assert(read_text_attr(file_id, "TYWRF_SUGGESTED_GATE_INTERVAL_OPTION") == "--interval-minutes 10");
  assert(read_times(file_id).substr(0, 19) == kTenMinuteEnd);
  check_nc(nc_close(file_id), "close real KROSA output");

  const auto source_root = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
  const std::string gate_command =
      "cd " +
      shell_quote(source_root) +
      " && UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python "
      "tools/cycle_gate.py --reference-dir " +
      shell_quote(kKrosaTenMinuteReferenceDir) + " --candidate-dir " + shell_quote(candidate_dir) +
      " --start " + shell_quote(std::string(kCycleStart)) + " --end " +
      shell_quote(std::string(kTenMinuteEnd)) + " --domain d02 --interval-minutes 10 --pretty";
  run_command_capture_expect_failure(gate_command, gate_report);

  std::ifstream report_file(gate_report);
  std::ostringstream report_buffer;
  report_buffer << report_file.rdbuf();
  const auto report = report_buffer.str();
  assert(report.find("\"status\": \"failed\"") != std::string::npos);
  assert(report.find("\"interval_minutes\": 10") != std::string::npos);
  assert(report.find("\"first_failure\":") != std::string::npos);
  assert(report.find("\"end_time\": \"2025-07-26_00:10:00\"") != std::string::npos);
  assert(report.find("TYWRF_VALIDATION_GATE_ONLY=true") != std::string::npos);

  const auto parent_fill_candidate_dir = root / "real_krosa_parent_fill_candidate";
  const auto parent_fill_output =
      parent_fill_candidate_dir / "wrfout_d02_2025-07-26_00:10:00_parent_fill";
  const auto parent_fill_report = root / "real_krosa_parent_fill_report.json";
  std::filesystem::create_directories(parent_fill_candidate_dir);
  std::filesystem::remove(parent_fill_output);
  std::filesystem::remove(parent_fill_report);

  run_command(
      base_command(executable, state, reference_end, parent_fill_output) +
      " --diagnostic-remap-parent-fill --from-parent-start 114,96 --to-parent-start 126,103 "
      "--parent-fill-state " +
      shell_quote(reference_end) + " > " + shell_quote(parent_fill_report));

  check_nc(
      nc_open(parent_fill_output.string().c_str(), NC_NOWRITE, &file_id),
      "open real KROSA parent-fill output");
  assert(
      read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") ==
      "cpp_skeleton_remap_parent_fill_diagnostic");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_REMAP_PARENT_FILL") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_REMAP_OVERLAP") == "false");
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") == "false");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "false");
  assert(read_text_attr(file_id, "TYWRF_NOT_PHYSICAL") == "true");
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") == "false");
  assert(read_text_attr(file_id, "TYWRF_MINIMUM_STATIC_REFRESH_FIELDS") == "XLAT,XLONG,HGT");
  assert(
      read_text_attr(file_id, "TYWRF_STAGGERED_STATIC_COORDS_STATUS") ==
      "pending_unless_emitted_later");
  assert(
      read_text_attr(file_id, "TYWRF_P_DERIVED_REFRESH_STATUS") ==
      "pending_derive_or_recompute_after_parent_fill_not_direct_wrf_parent_fill");
  assert_parent_fill_pressure_refresh_attrs(file_id);
  assert(
      read_text_attr(file_id, "TYWRF_DIRECT_WRF_END_STATE_ORACLE_STATUS") ==
      "diagnostic_only_nonphysical_non_gate");
  assert(read_text_attr(file_id, "TYWRF_PARENT_FILL_SOURCE") == reference_end.string());
  assert(read_double_attr(file_id, "TYWRF_REMAP_PARENT_FILL_FIELD_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_REMAP_PARENT_FILL_POINT_COUNT") > 0.0);
  assert(read_double_attr(file_id, "DX") == 2000.0);
  assert(read_double_attr(file_id, "DY") == 2000.0);
  assert(read_times(file_id).substr(0, 19) == kTenMinuteEnd);
  check_nc(nc_close(file_id), "close real KROSA parent-fill output");

  std::ifstream parent_fill_report_file(parent_fill_report);
  std::ostringstream parent_fill_report_buffer;
  parent_fill_report_buffer << parent_fill_report_file.rdbuf();
  const auto parent_fill_report_text = parent_fill_report_buffer.str();
  assert(
      parent_fill_report_text.find(
          "\"minimum_static_refresh_fields\": [\"XLAT\", \"XLONG\", \"HGT\"]") !=
      std::string::npos);
  assert(
      parent_fill_report_text.find(
          "\"staggered_static_coords_status\": \"pending_unless_emitted_later\"") !=
      std::string::npos);
  assert(
      parent_fill_report_text.find(
          "\"p_derived_refresh_status\": "
          "\"pending_derive_or_recompute_after_parent_fill_not_direct_wrf_parent_fill\"") !=
      std::string::npos);
  assert_parent_fill_pressure_refresh_report(parent_fill_report_text);
  assert(
      parent_fill_report_text.find(
          "\"direct_wrf_end_state_oracle_status\": "
          "\"diagnostic_only_nonphysical_non_gate\"") != std::string::npos);

  const auto pressure_candidate_dir = root / "real_krosa_pressure_refresh_candidate";
  const auto pressure_output =
      pressure_candidate_dir / "wrfout_d02_2025-07-26_00:10:00_pressure_refresh";
  const auto pressure_report = root / "real_krosa_pressure_refresh_report.json";
  std::filesystem::create_directories(pressure_candidate_dir);
  std::filesystem::remove(pressure_output);
  std::filesystem::remove(pressure_report);

  run_command(
      base_command(executable, state, reference_end, pressure_output) +
      " --diagnostic-remap-parent-fill --diagnostic-remap-pressure-refresh "
      "--from-parent-start 114,96 --to-parent-start 126,103 --parent-fill-state " +
      shell_quote(state) + " > " + shell_quote(pressure_report));

  check_nc(
      nc_open(pressure_output.string().c_str(), NC_NOWRITE, &file_id),
      "open real KROSA pressure-refresh output");
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") == "false");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "false");
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") == "false");
  assert_parent_fill_pressure_refresh_opt_in_attrs(file_id, reference_end);
  check_nc(nc_close(file_id), "close real KROSA pressure-refresh output");

  std::ifstream pressure_report_file(pressure_report);
  std::ostringstream pressure_report_buffer;
  pressure_report_buffer << pressure_report_file.rdbuf();
  const auto pressure_report_text = pressure_report_buffer.str();
  assert(pressure_report_text.find("\"gate_candidate\": false") != std::string::npos);
  assert(pressure_report_text.find("\"integrator_output\": false") != std::string::npos);
  assert_parent_fill_pressure_refresh_opt_in_report(pressure_report_text, reference_end);
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: skeleton_cycle_tool_test <tywrf_skeleton_cycle path>\n";
      return 2;
    }
    const auto executable = std::filesystem::path(argv[1]);
    const auto root = std::filesystem::temp_directory_path() / "tywrf_skeleton_cycle_tool_test";
    std::filesystem::create_directories(root);

    const auto state = root / "wrfout_d02_2025-07-26_00:00:00";
    const auto end_template = root / "wrfout_d02_2025-07-26_00:10:00";
    const auto output = root / "tywrf_cpp_skeleton_wrfout_d02_2025-07-26_00:10:00";
    const auto remap_state = root / "wrfout_d02_krosa_remap_state";
    const auto remap_output = root / "tywrf_cpp_skeleton_remap_wrfout_d02";
    const auto remap_report = root / "tywrf_cpp_skeleton_remap_report.json";
    const auto parent_fill_state = root / "wrfout_d02_parent_fill_state";
    const auto parent_fill_template = root / "wrfout_d02_parent_fill_template";
    const auto parent_fill_output = root / "tywrf_cpp_skeleton_parent_fill_wrfout_d02";
    const auto parent_fill_report = root / "tywrf_cpp_skeleton_parent_fill_report.json";
    const auto pressure_template = root / "wrfout_d02_pressure_refresh_template";
    const auto pressure_output =
        root / "tywrf_cpp_skeleton_parent_fill_pressure_refresh_wrfout_d02";
    const auto pressure_report =
        root / "tywrf_cpp_skeleton_parent_fill_pressure_refresh_report.json";
    std::filesystem::remove(state);
    std::filesystem::remove(end_template);
    std::filesystem::remove(output);
    std::filesystem::remove(remap_state);
    std::filesystem::remove(remap_output);
    std::filesystem::remove(remap_report);
    std::filesystem::remove(parent_fill_state);
    std::filesystem::remove(parent_fill_template);
    std::filesystem::remove(parent_fill_output);
    std::filesystem::remove(parent_fill_report);
    std::filesystem::remove(pressure_template);
    std::filesystem::remove(pressure_output);
    std::filesystem::remove(pressure_report);
    create_synthetic_wrfout(state, 2000.0, "2025-07-26_00:00:00", true);
    create_synthetic_wrfout(
        end_template,
        2000.0,
        "2025-07-26_00:10:00",
        false,
        70000.0F,
        80000.0F,
        90000.0F);

    run_command(base_command(executable, state, state, output));
    assert_output(output);

    create_synthetic_wrfout(
        remap_state,
        2000.0,
        "2025-07-26_00:00:00",
        true,
        40000.0F,
        50000.0F,
        60000.0F,
        210,
        210,
        2,
        3);
    run_command(
        base_command(executable, remap_state, remap_state, remap_output) +
        " --diagnostic-remap-overlap --from-parent-start 114,96 --to-parent-start 115,96 > " +
        shell_quote(remap_report));
    assert_remap_output(remap_output, remap_report);

    create_synthetic_wrfout(
        parent_fill_state,
        2000.0,
        "2025-07-26_00:10:00",
        true,
        170000.0F,
        180000.0F,
        190000.0F,
        210,
        210,
        2,
        3,
        110000.0F,
        120000.0F,
        130000.0F,
        240000.0F,
        280000.0F,
        2000.0F,
        3000.0F,
        290000.0F);
    create_synthetic_wrfout(
        parent_fill_template,
        2000.0,
        "2025-07-26_00:10:00",
        false,
        140000.0F,
        150000.0F,
        160000.0F,
        210,
        210,
        2,
        3);
    run_command(
        base_command(executable, remap_state, parent_fill_template, parent_fill_output) +
        " --diagnostic-remap-parent-fill --from-parent-start 114,96 --to-parent-start 115,96 "
        "--parent-fill-state " +
        shell_quote(parent_fill_state) + " --parent-fill-time-index 0 > " +
        shell_quote(parent_fill_report));
    assert_parent_fill_output(parent_fill_output, parent_fill_state, parent_fill_report);

    create_synthetic_wrfout(
        pressure_template,
        2000.0,
        "2025-07-26_00:10:00",
        false,
        240000.0F,
        250000.0F,
        20.0F,
        210,
        210,
        2,
        3,
        10000.0F,
        20000.0F,
        30000.0F,
        40000.0F,
        80000.0F,
        1000.0F,
        2000.0F,
        90000.0F,
        true,
        false);
    run_command(
        base_command(executable, remap_state, pressure_template, pressure_output) +
        " --diagnostic-remap-parent-fill --diagnostic-remap-pressure-refresh "
        "--from-parent-start 114,96 --to-parent-start 115,96 --parent-fill-state " +
        shell_quote(parent_fill_state) + " --parent-fill-time-index 0 > " +
        shell_quote(pressure_report));
    assert_parent_fill_pressure_refresh_output(
        pressure_output,
        parent_fill_state,
        pressure_template,
        pressure_report);

    const auto bad_state = root / "wrfout_d02_bad_dx";
    const auto bad_output = root / "tywrf_bad_dx_output";
    std::filesystem::remove(bad_state);
    std::filesystem::remove(bad_output);
    create_synthetic_wrfout(bad_state, 1000.0, "2025-07-26_00:00:00", true);
    run_command_expect_failure(base_command(executable, bad_state, state, bad_output));

    run_real_krosa_smoke_if_available(executable, root);

    std::filesystem::remove_all(root);
  } catch (const std::exception& error) {
    std::cerr << "skeleton cycle tool test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
