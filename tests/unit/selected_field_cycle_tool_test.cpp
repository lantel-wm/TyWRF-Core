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
constexpr float kTolerance = 1.0e-3F;

void check_nc(const int status, const std::string_view operation) {
  if (status == NC_NOERR) {
    return;
  }
  throw std::runtime_error(std::string(operation) + ": " + nc_strerror(status));
}

float linear_3d(
    const int k,
    const int j,
    const int i,
    const float base,
    const float ax,
    const float ay,
    const float ak) {
  return base + ax * static_cast<float>(i) + ay * static_cast<float>(j) +
         ak * static_cast<float>(k);
}

float linear_2d(const int j, const int i, const float base, const float ax, const float ay) {
  return base + ax * static_cast<float>(i) + ay * static_cast<float>(j);
}

float parent_linear(
    const float base,
    const float ax,
    const float ay,
    const float ak,
    const double x,
    const double y,
    const int k) {
  return base + ax * static_cast<float>(x) + ay * static_cast<float>(y) +
         ak * static_cast<float>(k);
}

void put_times(const int file_id, const int variable_id, const std::string_view value) {
  constexpr std::size_t date_str_len = 19;
  std::vector<char> buffer(date_str_len, ' ');
  std::copy(value.begin(), value.end(), buffer.begin());
  const std::size_t start[2] = {0, 0};
  const std::size_t count[2] = {1, date_str_len};
  check_nc(nc_put_vara_text(file_id, variable_id, start, count, buffer.data()), "write Times");
}

void put_3d_linear(
    const int file_id,
    const int variable_id,
    const int nz,
    const int ny,
    const int nx,
    const float base,
    const float ax,
    const float ay,
    const float ak) {
  std::vector<float> values(static_cast<std::size_t>(nz * ny * nx));
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        values[(static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
                static_cast<std::size_t>(j)) *
                   static_cast<std::size_t>(nx) +
               static_cast<std::size_t>(i)] = linear_3d(k, j, i, base, ax, ay, ak);
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

void put_2d_linear(
    const int file_id,
    const int variable_id,
    const int ny,
    const int nx,
    const float base,
    const float ax,
    const float ay) {
  std::vector<float> values(static_cast<std::size_t>(ny * nx));
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      values[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
             static_cast<std::size_t>(i)] = linear_2d(j, i, base, ax, ay);
    }
  }
  const std::size_t start[3] = {0, 0, 0};
  const std::size_t count[3] = {1, static_cast<std::size_t>(ny), static_cast<std::size_t>(nx)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, count, values.data()), "write 2D var");
}

struct FixtureShape {
  int mass_nx = 0;
  int mass_ny = 0;
  int mass_nz = 2;
  int full_nz = 3;
};

struct DefinedFixtureVars {
  int times = -1;
  int xlat = -1;
  int xlong = -1;
  int hgt = -1;
  int u = -1;
  int v = -1;
  int t = -1;
  int ph = -1;
  int mu = -1;
  int p = -1;
  int qvapor = -1;
};

DefinedFixtureVars define_fixture_vars(
    const int file_id,
    const FixtureShape shape,
    const bool state_fields,
    const int time_dim,
    const int date_str_len_dim,
    const int bottom_top_dim,
    const int bottom_top_stag_dim,
    const int south_north_dim,
    const int south_north_stag_dim,
    const int west_east_dim,
    const int west_east_stag_dim) {
  DefinedFixtureVars vars;
  const int times_dims[2] = {time_dim, date_str_len_dim};
  const int static_dims[3] = {time_dim, south_north_dim, west_east_dim};
  const int u_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_stag_dim};
  const int v_dims[4] = {time_dim, bottom_top_dim, south_north_stag_dim, west_east_dim};
  const int w_dims[4] = {time_dim, bottom_top_stag_dim, south_north_dim, west_east_dim};
  const int mass_3d_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_dim};
  const int mass_2d_dims[3] = {time_dim, south_north_dim, west_east_dim};
  check_nc(nc_def_var(file_id, "Times", NC_CHAR, 2, times_dims, &vars.times), "define Times");
  check_nc(nc_def_var(file_id, "XLAT", NC_FLOAT, 3, static_dims, &vars.xlat), "define XLAT");
  check_nc(nc_def_var(file_id, "XLONG", NC_FLOAT, 3, static_dims, &vars.xlong), "define XLONG");
  check_nc(nc_def_var(file_id, "HGT", NC_FLOAT, 3, static_dims, &vars.hgt), "define HGT");
  if (state_fields) {
    check_nc(nc_def_var(file_id, "U", NC_FLOAT, 4, u_dims, &vars.u), "define U");
    check_nc(nc_def_var(file_id, "V", NC_FLOAT, 4, v_dims, &vars.v), "define V");
    check_nc(nc_def_var(file_id, "T", NC_FLOAT, 4, mass_3d_dims, &vars.t), "define T");
    check_nc(nc_def_var(file_id, "PH", NC_FLOAT, 4, w_dims, &vars.ph), "define PH");
    check_nc(nc_def_var(file_id, "MU", NC_FLOAT, 3, mass_2d_dims, &vars.mu), "define MU");
    check_nc(nc_def_var(file_id, "P", NC_FLOAT, 4, mass_3d_dims, &vars.p), "define P");
    check_nc(
        nc_def_var(file_id, "QVAPOR", NC_FLOAT, 4, mass_3d_dims, &vars.qvapor),
        "define QVAPOR");
  }
  (void)shape;
  return vars;
}

void create_wrf_fixture(
    const std::filesystem::path& path,
    const double dx,
    const FixtureShape shape,
    const bool state_fields,
    const bool parent_values) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create fixture");
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
  check_nc(nc_def_dim(file_id, "bottom_top", shape.mass_nz, &bottom_top_dim), "define bottom_top");
  check_nc(
      nc_def_dim(file_id, "bottom_top_stag", shape.full_nz, &bottom_top_stag_dim),
      "define bottom_top_stag");
  check_nc(
      nc_def_dim(file_id, "south_north", shape.mass_ny, &south_north_dim),
      "define south_north");
  check_nc(
      nc_def_dim(file_id, "south_north_stag", shape.mass_ny + 1, &south_north_stag_dim),
      "define south_north_stag");
  check_nc(nc_def_dim(file_id, "west_east", shape.mass_nx, &west_east_dim), "define west_east");
  check_nc(
      nc_def_dim(file_id, "west_east_stag", shape.mass_nx + 1, &west_east_stag_dim),
      "define west_east_stag");

  const auto vars = define_fixture_vars(
      file_id,
      shape,
      state_fields,
      time_dim,
      date_str_len_dim,
      bottom_top_dim,
      bottom_top_stag_dim,
      south_north_dim,
      south_north_stag_dim,
      west_east_dim,
      west_east_stag_dim);
  check_nc(nc_enddef(file_id), "end definitions");

  put_times(file_id, vars.times, state_fields ? kCycleStart : kTenMinuteEnd);
  put_2d_linear(file_id, vars.xlat, shape.mass_ny, shape.mass_nx, 40'000.0F, 1.0F, 10.0F);
  put_2d_linear(file_id, vars.xlong, shape.mass_ny, shape.mass_nx, 50'000.0F, 1.0F, 10.0F);
  put_2d_linear(file_id, vars.hgt, shape.mass_ny, shape.mass_nx, 60'000.0F, 1.0F, 10.0F);
  if (state_fields) {
    if (parent_values) {
      put_3d_linear(file_id, vars.u, shape.mass_nz, shape.mass_ny, shape.mass_nx + 1, 1'000.0F, 7.0F, 13.0F, 101.0F);
      put_3d_linear(file_id, vars.v, shape.mass_nz, shape.mass_ny + 1, shape.mass_nx, 2'000.0F, 11.0F, 17.0F, 103.0F);
      put_3d_linear(file_id, vars.t, shape.mass_nz, shape.mass_ny, shape.mass_nx, 3'000.0F, 5.0F, 3.0F, 10.0F);
      put_3d_linear(file_id, vars.ph, shape.full_nz, shape.mass_ny, shape.mass_nx, 4'000.0F, 2.0F, 4.0F, 50.0F);
      put_2d_linear(file_id, vars.mu, shape.mass_ny, shape.mass_nx, 5'000.0F, 19.0F, 23.0F);
      put_3d_linear(file_id, vars.p, shape.mass_nz, shape.mass_ny, shape.mass_nx, 6'000.0F, 3.0F, 2.0F, 7.0F);
      put_3d_linear(file_id, vars.qvapor, shape.mass_nz, shape.mass_ny, shape.mass_nx, 7'000.0F, 29.0F, 31.0F, 107.0F);
    } else {
      put_3d_linear(file_id, vars.u, shape.mass_nz, shape.mass_ny, shape.mass_nx + 1, 10'000.0F, 1.0F, 10.0F, 100.0F);
      put_3d_linear(file_id, vars.v, shape.mass_nz, shape.mass_ny + 1, shape.mass_nx, 20'000.0F, 1.0F, 10.0F, 100.0F);
      put_3d_linear(file_id, vars.t, shape.mass_nz, shape.mass_ny, shape.mass_nx, 30'000.0F, 1.0F, 10.0F, 100.0F);
      put_3d_linear(file_id, vars.ph, shape.full_nz, shape.mass_ny, shape.mass_nx, 40'000.0F, 1.0F, 10.0F, 100.0F);
      put_2d_linear(file_id, vars.mu, shape.mass_ny, shape.mass_nx, 50'000.0F, 1.0F, 10.0F);
      put_3d_linear(file_id, vars.p, shape.mass_nz, shape.mass_ny, shape.mass_nx, 60'000.0F, 1.0F, 10.0F, 100.0F);
      put_3d_linear(file_id, vars.qvapor, shape.mass_nz, shape.mass_ny, shape.mass_nx, 70'000.0F, 1.0F, 10.0F, 100.0F);
    }
  }
  check_nc(nc_close(file_id), "close fixture");
}

[[nodiscard]] std::string shell_quote(const std::filesystem::path& path) {
  return "'" + path.string() + "'";
}

[[nodiscard]] std::string shell_quote(const std::string& value) {
  return "'" + value + "'";
}

int run_command_status(const std::string& command) {
  return std::system(command.c_str());
}

void run_command(const std::string& command) {
  const int status = run_command_status(command);
  if (status != 0) {
    throw std::runtime_error("command failed: " + command);
  }
}

[[nodiscard]] std::string base_command(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output) {
  return shell_quote(executable) + " --d01-start-state " + shell_quote(d01_start) +
         " --d02-start-state " + shell_quote(d02_start) + " --template " +
         shell_quote(template_path) + " --output " + shell_quote(output) +
         " --cycle-start " + shell_quote(std::string(kCycleStart)) + " --cycle-end " +
         shell_quote(std::string(kTenMinuteEnd)) + " --from-parent-start 2,2"
         " --to-parent-start 3,2 --times " +
         shell_quote(std::string(kTenMinuteEnd)) + " --pretty";
}

[[nodiscard]] std::string read_text_attr(const int file_id, const std::string_view name) {
  std::size_t length = 0;
  check_nc(nc_inq_attlen(file_id, NC_GLOBAL, std::string(name).c_str(), &length), "read attr len");
  std::string value(length, '\0');
  check_nc(nc_get_att_text(file_id, NC_GLOBAL, std::string(name).c_str(), value.data()), "read attr");
  return value;
}

[[nodiscard]] double read_double_attr(const int file_id, const std::string_view name) {
  double value = 0.0;
  check_nc(nc_get_att_double(file_id, NC_GLOBAL, std::string(name).c_str(), &value), "read attr");
  return value;
}

[[nodiscard]] bool has_text_attr(const int file_id, const std::string_view name) {
  return nc_inq_att(file_id, NC_GLOBAL, std::string(name).c_str(), nullptr, nullptr) == NC_NOERR;
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

void expect_close(const float actual, const float expected, const std::string_view label) {
  if (std::fabs(actual - expected) > kTolerance) {
    std::cerr << label << " mismatch: got " << actual << ", expected " << expected << '\n';
    assert(false);
  }
}

void assert_successful_candidate(const std::filesystem::path& output) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open output");
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") == "true");
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") == "true");
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "false");
  assert(read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") == "selected_field_integrator_v0");
  assert(read_text_attr(file_id, "TYWRF_D02_RESOLUTION_CHECK") == "d02_2km");
  assert(!has_text_attr(file_id, "TYWRF_NOT_PHYSICAL"));
  assert(read_double_attr(file_id, "DX") == 2000.0);
  assert(read_double_attr(file_id, "DY") == 2000.0);
  assert(read_double_attr(file_id, "TYWRF_PARENT_GRID_RATIO") == 5.0);
  assert(read_double_attr(file_id, "TYWRF_SELECTED_FIELD_CHANGED_POINTS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_INTERPOLATED_POINTS") > 0.0);
  assert(read_times(file_id).substr(0, 19) == kTenMinuteEnd);

  expect_close(read_3d_value(file_id, "U", 1, 2, 0), linear_3d(1, 2, 5, 10'000.0F, 1.0F, 10.0F, 100.0F), "U overlap");
  expect_close(read_3d_value(file_id, "V", 1, 2, 0), linear_3d(1, 2, 5, 20'000.0F, 1.0F, 10.0F, 100.0F), "V overlap");
  expect_close(read_2d_value(file_id, "MU", 2, 0), linear_2d(2, 5, 50'000.0F, 1.0F, 10.0F), "MU overlap");
  expect_close(read_3d_value(file_id, "QVAPOR", 1, 2, 0), linear_3d(1, 2, 5, 70'000.0F, 1.0F, 10.0F, 100.0F), "QVAPOR overlap");

  expect_close(read_3d_value(file_id, "U", 1, 0, 10), parent_linear(1'000.0F, 7.0F, 13.0F, 101.0F, 4.0, 1.0, 1), "U exposed");
  expect_close(read_3d_value(file_id, "V", 1, 10, 9), parent_linear(2'000.0F, 11.0F, 17.0F, 103.0F, 3.8, 3.0, 1), "V exposed");
  expect_close(read_2d_value(file_id, "MU", 0, 9), parent_linear(5'000.0F, 19.0F, 23.0F, 0.0F, 3.8, 1.0, 0), "MU exposed");
  expect_close(read_3d_value(file_id, "QVAPOR", 1, 9, 9), parent_linear(7'000.0F, 29.0F, 31.0F, 107.0F, 3.8, 2.8, 1), "QVAPOR exposed");

  expect_close(read_3d_value(file_id, "T", 1, 9, 9), linear_3d(1, 9, 9, 30'000.0F, 1.0F, 10.0F, 100.0F), "T preserved exposed");
  expect_close(read_3d_value(file_id, "PH", 2, 9, 9), linear_3d(2, 9, 9, 40'000.0F, 1.0F, 10.0F, 100.0F), "PH preserved exposed");
  expect_close(read_3d_value(file_id, "P", 1, 9, 9), linear_3d(1, 9, 9, 60'000.0F, 1.0F, 10.0F, 100.0F), "P preserved exposed");
  check_nc(nc_close(file_id), "close output");
}

void run_rejection_tests(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& root) {
  const auto oracle_output = root / "oracle_rejected";
  const auto oracle_status = run_command_status(
      base_command(executable, d01_start, d02_start, template_path, oracle_output) +
      " --end-state " + shell_quote(d02_start) + " >/dev/null 2>&1");
  assert(oracle_status != 0);
  assert(!std::filesystem::exists(oracle_output));

  const auto bad_template = root / "bad_template";
  create_wrf_fixture(bad_template, 4000.0, FixtureShape{10, 10}, false, false);
  const auto bad_output = root / "bad_resolution_output";
  const auto bad_status = run_command_status(
      base_command(executable, d01_start, d02_start, bad_template, bad_output) +
      " >/dev/null 2>&1");
  assert(bad_status != 0);
  assert(!std::filesystem::exists(bad_output));

  const auto missing_static_output = root / "missing_static_output";
  const auto missing_static_status = run_command_status(
      base_command(executable, d01_start, d02_start, template_path, missing_static_output) +
      " --variables U,V,T,PH,MU,P,QVAPOR >/dev/null 2>&1");
  assert(missing_static_status != 0);
  assert(!std::filesystem::exists(missing_static_output));
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: selected_field_cycle_tool_test <tywrf_selected_field_cycle path>\n";
      return 2;
    }
    const auto executable = std::filesystem::path(argv[1]);
    const auto root = std::filesystem::temp_directory_path() / "tywrf_selected_field_cycle_tool_test";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);

    const auto d01_start = root / "wrfout_d01_2025-07-26_00:00:00";
    const auto d02_start = root / "wrfout_d02_2025-07-26_00:00:00";
    const auto template_path = root / "wrfout_d02_template";
    const auto output = root / "tywrf_selected_field_d02_2025-07-26_00:10:00";
    create_wrf_fixture(d01_start, 10000.0, FixtureShape{8, 8}, true, true);
    create_wrf_fixture(d02_start, 2000.0, FixtureShape{10, 10}, true, false);
    create_wrf_fixture(template_path, 2000.0, FixtureShape{10, 10}, false, false);

    run_command(base_command(executable, d01_start, d02_start, template_path, output));
    assert_successful_candidate(output);
    run_rejection_tests(executable, d01_start, d02_start, template_path, root);

    std::filesystem::remove_all(root);
  } catch (const std::exception& error) {
    std::cerr << "selected field cycle tool test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
