#include <netcdf.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
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

struct StateOffsets {
  float u = 10000.0F;
  float v = 11000.0F;
  float w = 12000.0F;
  float t = 20000.0F;
  float qvapor = 21000.0F;
  float mu = 30000.0F;
};

void check_nc(const int status, const std::string_view operation) {
  if (status == NC_NOERR) {
    return;
  }
  throw std::runtime_error(std::string(operation) + ": " + nc_strerror(status));
}

float value_3d(const int k, const int j, const int i, const float offset) {
  return offset + 100.0F * static_cast<float>(k) + 10.0F * static_cast<float>(j) +
         static_cast<float>(i);
}

float value_2d(const int j, const int i, const float offset) {
  return offset + 10.0F * static_cast<float>(j) + static_cast<float>(i);
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
               static_cast<std::size_t>(i)] = value_3d(k, j, i, offset);
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
             static_cast<std::size_t>(i)] = value_2d(j, i, offset);
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
    const StateOffsets offsets = {},
    const float xlat_offset = 40000.0F,
    const float xlong_offset = 50000.0F,
    const float hgt_offset = 60000.0F) {
  constexpr int mass_nx = 4;
  constexpr int mass_ny = 3;
  constexpr int mass_nz = 2;
  constexpr int full_nz = 3;

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
  int v_var = -1;
  int w_var = -1;
  int t_var = -1;
  int qvapor_var = -1;
  int mu_var = -1;
  const int times_dims[2] = {time_dim, date_str_len_dim};
  const int static_dims[3] = {time_dim, south_north_dim, west_east_dim};
  const int u_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_stag_dim};
  const int v_dims[4] = {time_dim, bottom_top_dim, south_north_stag_dim, west_east_dim};
  const int w_dims[4] = {time_dim, bottom_top_stag_dim, south_north_dim, west_east_dim};
  const int mass_3d_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_dim};
  const int mass_2d_dims[3] = {time_dim, south_north_dim, west_east_dim};
  check_nc(nc_def_var(file_id, "Times", NC_CHAR, 2, times_dims, &times_var), "define Times");
  check_nc(nc_def_var(file_id, "XLAT", NC_FLOAT, 3, static_dims, &xlat_var), "define XLAT");
  check_nc(nc_def_var(file_id, "XLONG", NC_FLOAT, 3, static_dims, &xlong_var), "define XLONG");
  check_nc(nc_def_var(file_id, "HGT", NC_FLOAT, 3, static_dims, &hgt_var), "define HGT");
  if (state_fields) {
    check_nc(nc_def_var(file_id, "U", NC_FLOAT, 4, u_dims, &u_var), "define U");
    check_nc(nc_def_var(file_id, "V", NC_FLOAT, 4, v_dims, &v_var), "define V");
    check_nc(nc_def_var(file_id, "W", NC_FLOAT, 4, w_dims, &w_var), "define W");
    check_nc(nc_def_var(file_id, "T", NC_FLOAT, 4, mass_3d_dims, &t_var), "define T");
    check_nc(
        nc_def_var(file_id, "QVAPOR", NC_FLOAT, 4, mass_3d_dims, &qvapor_var),
        "define QVAPOR");
    check_nc(nc_def_var(file_id, "MU", NC_FLOAT, 3, mass_2d_dims, &mu_var), "define MU");
  }
  check_nc(nc_enddef(file_id), "end definitions");

  put_times(file_id, times_var, time_value);
  put_2d(file_id, xlat_var, mass_ny, mass_nx, xlat_offset);
  put_2d(file_id, xlong_var, mass_ny, mass_nx, xlong_offset);
  put_2d(file_id, hgt_var, mass_ny, mass_nx, hgt_offset);
  if (state_fields) {
    put_3d(file_id, u_var, mass_nz, mass_ny, mass_nx + 1, offsets.u);
    put_3d(file_id, v_var, mass_nz, mass_ny + 1, mass_nx, offsets.v);
    put_3d(file_id, w_var, full_nz, mass_ny, mass_nx, offsets.w);
    put_3d(file_id, t_var, mass_nz, mass_ny, mass_nx, offsets.t);
    put_3d(file_id, qvapor_var, mass_nz, mass_ny, mass_nx, offsets.qvapor);
    put_2d(file_id, mu_var, mass_ny, mass_nx, offsets.mu);
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
  const std::size_t start[4] = {
      0,
      static_cast<std::size_t>(k),
      static_cast<std::size_t>(j),
      static_cast<std::size_t>(i)};
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

bool nearly_equal(const float lhs, const float rhs) {
  return std::abs(lhs - rhs) < 1.0e-3F;
}

void assert_synthetic_output(const std::filesystem::path& output, const StateOffsets& end_offsets) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open output");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ONLY") == "true");
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") == "false");
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") == "false");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "true");
  assert(read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") == "reference_delta_oracle");
  assert(read_text_attr(file_id, "TYWRF_NOT_PHYSICAL") == "true");
  assert(read_text_attr(file_id, "TYWRF_CYCLE_START") == kCycleStart);
  assert(read_text_attr(file_id, "TYWRF_CYCLE_END") == kTenMinuteEnd);
  assert(read_double_attr(file_id, "TYWRF_DT_SECONDS") == 600.0);
  assert(read_double_attr(file_id, "TYWRF_D02_TENDENCY_DT_SECONDS") == 600.0);
  assert(read_double_attr(file_id, "TYWRF_D02_EXPLICIT_TENDENCY_APPLY_COUNT") == 75.0);
  assert(read_double_attr(file_id, "DX") == 2000.0);
  assert(read_double_attr(file_id, "DY") == 2000.0);
  assert(read_times(file_id).substr(0, 19) == kTenMinuteEnd);

  assert(nearly_equal(read_3d_value(file_id, "U", 1, 2, 4), value_3d(1, 2, 4, end_offsets.u)));
  assert(nearly_equal(read_3d_value(file_id, "V", 1, 3, 2), value_3d(1, 3, 2, end_offsets.v)));
  assert(nearly_equal(read_3d_value(file_id, "W", 2, 1, 3), value_3d(2, 1, 3, end_offsets.w)));
  assert(nearly_equal(read_3d_value(file_id, "T", 1, 2, 3), value_3d(1, 2, 3, end_offsets.t)));
  assert(nearly_equal(
      read_3d_value(file_id, "QVAPOR", 0, 1, 2),
      value_3d(0, 1, 2, end_offsets.qvapor)));
  assert(nearly_equal(read_2d_value(file_id, "MU", 2, 3), value_2d(2, 3, end_offsets.mu)));
  assert(read_2d_value(file_id, "XLAT", 2, 3) == value_2d(2, 3, 70000.0F));
  assert(read_2d_value(file_id, "XLONG", 1, 2) == value_2d(1, 2, 80000.0F));
  assert(read_2d_value(file_id, "HGT", 0, 1) == value_2d(0, 1, 90000.0F));
  check_nc(nc_close(file_id), "close output");
}

std::string synthetic_command(
    const std::filesystem::path& executable,
    const std::filesystem::path& start_state,
    const std::filesystem::path& end_state,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output) {
  return shell_quote(executable) + " --start-state " + shell_quote(start_state) +
         " --end-state " + shell_quote(end_state) + " --template " + shell_quote(template_path) +
         " --output " + shell_quote(output) + " --cycle-start " +
         shell_quote(std::string(kCycleStart)) + " --cycle-end " +
         shell_quote(std::string(kTenMinuteEnd)) + " --times " +
         shell_quote(std::string(kTenMinuteEnd)) +
         " --domain d02 --variables Times,XLAT,XLONG,HGT,U,V,W,T,QVAPOR,MU --pretty";
}

void run_real_krosa_smoke_if_available(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto start_state = kKrosaTenMinuteReferenceDir / "wrfout_d02_2025-07-26_00:00:00";
  const auto end_state = kKrosaTenMinuteReferenceDir / "wrfout_d02_2025-07-26_00:10:00";
  if (!std::filesystem::exists(start_state) || !std::filesystem::exists(end_state)) {
    std::cerr << "skipping real KROSA tendency-cycle smoke; reference files are not present\n";
    return;
  }

  const auto output = root / "real_krosa_oracle" / "wrfout_d02_2025-07-26_00:10:00";
  std::filesystem::create_directories(output.parent_path());
  std::filesystem::remove(output);
  run_command(
      shell_quote(executable) + " --start-state " + shell_quote(start_state) +
      " --end-state " + shell_quote(end_state) + " --template " + shell_quote(start_state) +
      " --output " + shell_quote(output) + " --cycle-start " +
      shell_quote(std::string(kCycleStart)) + " --cycle-end " +
      shell_quote(std::string(kTenMinuteEnd)) + " --times " +
      shell_quote(std::string(kTenMinuteEnd)) +
      " --domain d02 --variables Times,XLAT,XLONG,HGT,U,V,T,PH,MU,P,QVAPOR");

  int output_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &output_id), "open real output");
  assert(read_text_attr(output_id, "TYWRF_DIAGNOSTIC_ONLY") == "true");
  assert(read_text_attr(output_id, "TYWRF_GATE_CANDIDATE") == "false");
  assert(read_text_attr(output_id, "TYWRF_CANDIDATE_KIND") == "reference_delta_oracle");
  assert(read_double_attr(output_id, "TYWRF_DT_SECONDS") == 600.0);
  const float output_u = read_3d_value(output_id, "U", 0, 0, 0);
  check_nc(nc_close(output_id), "close real output");

  int end_id = -1;
  check_nc(nc_open(end_state.string().c_str(), NC_NOWRITE, &end_id), "open real end");
  const float end_u = read_3d_value(end_id, "U", 0, 0, 0);
  check_nc(nc_close(end_id), "close real end");
  assert(nearly_equal(output_u, end_u));
  std::cerr << "real KROSA tendency-cycle smoke wrote " << output << '\n';
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: tendency_cycle_tool_test <tywrf_tendency_cycle path>\n";
      return 2;
    }
    const auto executable = std::filesystem::path(argv[1]);
    const auto root = std::filesystem::temp_directory_path() / "tywrf_tendency_cycle_tool_test";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);

    const auto start_state = root / "wrfout_d02_2025-07-26_00:00:00";
    const auto end_state = root / "wrfout_d02_2025-07-26_00:10:00";
    const auto template_path = root / "wrfout_d02_template_static";
    const auto output = root / "tywrf_reference_delta_oracle_d02_2025-07-26_00:10:00";
    const StateOffsets start_offsets{};
    const StateOffsets end_offsets{
        .u = start_offsets.u + 600.0F,
        .v = start_offsets.v + 1200.0F,
        .w = start_offsets.w + 1800.0F,
        .t = start_offsets.t + 300.0F,
        .qvapor = start_offsets.qvapor + 75.0F,
        .mu = start_offsets.mu + 150.0F,
    };
    create_synthetic_wrfout(start_state, 2000.0, kCycleStart, true, start_offsets);
    create_synthetic_wrfout(end_state, 2000.0, kTenMinuteEnd, true, end_offsets);
    create_synthetic_wrfout(
        template_path,
        2000.0,
        kTenMinuteEnd,
        false,
        {},
        70000.0F,
        80000.0F,
        90000.0F);

    run_command(synthetic_command(executable, start_state, end_state, template_path, output));
    assert_synthetic_output(output, end_offsets);
    run_real_krosa_smoke_if_available(executable, root);

    std::filesystem::remove_all(root);
  } catch (const std::exception& error) {
    std::cerr << "tendency cycle tool test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
