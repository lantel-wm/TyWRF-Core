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
    const int full_nz = 3) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create synthetic wrfout");
  check_nc(nc_put_att_double(file_id, NC_GLOBAL, "DX", NC_DOUBLE, 1, &dx), "write DX");
  check_nc(nc_put_att_double(file_id, NC_GLOBAL, "DY", NC_DOUBLE, 1, &dx), "write DY");

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
  const int times_dims[2] = {time_dim, date_str_len_dim};
  const int static_dims[3] = {time_dim, south_north_dim, west_east_dim};
  const int u_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_stag_dim};
  const int t_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_dim};
  const int mu_dims[3] = {time_dim, south_north_dim, west_east_dim};
  check_nc(nc_def_var(file_id, "Times", NC_CHAR, 2, times_dims, &times_var), "define Times");
  check_nc(nc_def_var(file_id, "XLAT", NC_FLOAT, 3, static_dims, &xlat_var), "define XLAT");
  check_nc(nc_def_var(file_id, "XLONG", NC_FLOAT, 3, static_dims, &xlong_var), "define XLONG");
  check_nc(nc_def_var(file_id, "HGT", NC_FLOAT, 3, static_dims, &hgt_var), "define HGT");
  if (state_fields) {
    check_nc(nc_def_var(file_id, "U", NC_FLOAT, 4, u_dims, &u_var), "define U");
    check_nc(nc_def_var(file_id, "T", NC_FLOAT, 4, t_dims, &t_var), "define T");
    check_nc(nc_def_var(file_id, "MU", NC_FLOAT, 3, mu_dims, &mu_var), "define MU");
  }
  check_nc(nc_enddef(file_id), "end definitions");

  put_times(file_id, times_var, time_value);
  put_2d(file_id, xlat_var, mass_ny, mass_nx, xlat_offset);
  put_2d(file_id, xlong_var, mass_ny, mass_nx, xlong_offset);
  put_2d(file_id, hgt_var, mass_ny, mass_nx, hgt_offset);
  if (state_fields) {
    put_3d(file_id, u_var, mass_nz, mass_ny, mass_nx + 1, 10000.0F);
    put_3d(file_id, t_var, mass_nz, mass_ny, mass_nx, 20000.0F);
    put_2d(file_id, mu_var, mass_ny, mass_nx, 30000.0F);
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
    std::filesystem::remove(state);
    std::filesystem::remove(end_template);
    std::filesystem::remove(output);
    std::filesystem::remove(remap_state);
    std::filesystem::remove(remap_output);
    std::filesystem::remove(remap_report);
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
