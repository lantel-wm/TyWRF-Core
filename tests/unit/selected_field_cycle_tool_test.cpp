#include <netcdf.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
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
constexpr float kTolerance = 1.0e-3F;
constexpr float kSyntheticHgtBase = 600.0F;

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

void put_time_vector(
    const int file_id,
    const int variable_id,
    const std::vector<float>& values) {
  const std::size_t start[2] = {0, 0};
  const std::size_t count[2] = {1, values.size()};
  check_nc(nc_put_vara_float(file_id, variable_id, start, count, values.data()), "write vector");
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
  int pb = -1;
  int phb = -1;
  int mub = -1;
  int psfc = -1;
  int u10 = -1;
  int v10 = -1;
  int t2 = -1;
  int q2 = -1;
  int rainc = -1;
  int rainnc = -1;
};

DefinedFixtureVars define_fixture_vars(
    const int file_id,
    const FixtureShape shape,
    const bool state_fields,
    const bool optional_state_fields,
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
    if (optional_state_fields) {
      check_nc(nc_def_var(file_id, "PB", NC_FLOAT, 4, mass_3d_dims, &vars.pb), "define PB");
      check_nc(nc_def_var(file_id, "PHB", NC_FLOAT, 4, w_dims, &vars.phb), "define PHB");
      check_nc(nc_def_var(file_id, "MUB", NC_FLOAT, 3, mass_2d_dims, &vars.mub), "define MUB");
      check_nc(nc_def_var(file_id, "PSFC", NC_FLOAT, 3, mass_2d_dims, &vars.psfc), "define PSFC");
      check_nc(nc_def_var(file_id, "U10", NC_FLOAT, 3, mass_2d_dims, &vars.u10), "define U10");
      check_nc(nc_def_var(file_id, "V10", NC_FLOAT, 3, mass_2d_dims, &vars.v10), "define V10");
      check_nc(nc_def_var(file_id, "T2", NC_FLOAT, 3, mass_2d_dims, &vars.t2), "define T2");
      check_nc(nc_def_var(file_id, "Q2", NC_FLOAT, 3, mass_2d_dims, &vars.q2), "define Q2");
      check_nc(nc_def_var(file_id, "RAINC", NC_FLOAT, 3, mass_2d_dims, &vars.rainc), "define RAINC");
      check_nc(nc_def_var(file_id, "RAINNC", NC_FLOAT, 3, mass_2d_dims, &vars.rainnc), "define RAINNC");
    }
  }
  (void)shape;
  return vars;
}

void create_wrf_fixture(
    const std::filesystem::path& path,
    const double dx,
    const FixtureShape shape,
    const bool state_fields,
    const bool parent_values,
    const bool optional_state_fields = false,
    const bool pressure_refresh_metadata = false,
    const bool invalid_parent_pressure_thermodynamics = false,
    const bool pressure_metadata_hgt_values = true) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create fixture");
  check_nc(nc_put_att_double(file_id, NC_GLOBAL, "DX", NC_DOUBLE, 1, &dx), "write DX");
  check_nc(nc_put_att_double(file_id, NC_GLOBAL, "DY", NC_DOUBLE, 1, &dx), "write DY");
  if (pressure_refresh_metadata) {
    const float p_top = 5'000.0F;
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
      optional_state_fields,
      time_dim,
      date_str_len_dim,
      bottom_top_dim,
      bottom_top_stag_dim,
      south_north_dim,
      south_north_stag_dim,
      west_east_dim,
      west_east_stag_dim);
  int c3f = -1;
  int c4f = -1;
  int c3h = -1;
  int c4h = -1;
  if (pressure_refresh_metadata) {
    const int full_dims[2] = {time_dim, bottom_top_stag_dim};
    const int mass_dims[2] = {time_dim, bottom_top_dim};
    check_nc(nc_def_var(file_id, "C3F", NC_FLOAT, 2, full_dims, &c3f), "define C3F");
    check_nc(nc_def_var(file_id, "C4F", NC_FLOAT, 2, full_dims, &c4f), "define C4F");
    check_nc(nc_def_var(file_id, "C3H", NC_FLOAT, 2, mass_dims, &c3h), "define C3H");
    check_nc(nc_def_var(file_id, "C4H", NC_FLOAT, 2, mass_dims, &c4h), "define C4H");
  }
  check_nc(nc_enddef(file_id), "end definitions");

  put_times(file_id, vars.times, state_fields ? kCycleStart : kTenMinuteEnd);
  put_2d_linear(file_id, vars.xlat, shape.mass_ny, shape.mass_nx, 40'000.0F, 1.0F, 10.0F);
  put_2d_linear(file_id, vars.xlong, shape.mass_ny, shape.mass_nx, 50'000.0F, 1.0F, 10.0F);
  if (pressure_refresh_metadata && pressure_metadata_hgt_values) {
    put_2d_linear(file_id, vars.hgt, shape.mass_ny, shape.mass_nx, 20.0F, 1.0F, 5.0F);
  } else {
    put_2d_linear(
        file_id, vars.hgt, shape.mass_ny, shape.mass_nx, kSyntheticHgtBase, 1.0F, 10.0F);
  }
  if (state_fields) {
    if (parent_values) {
      put_3d_linear(file_id, vars.u, shape.mass_nz, shape.mass_ny, shape.mass_nx + 1, 1'000.0F, 7.0F, 13.0F, 101.0F);
      put_3d_linear(file_id, vars.v, shape.mass_nz, shape.mass_ny + 1, shape.mass_nx, 2'000.0F, 11.0F, 17.0F, 103.0F);
      put_3d_linear(
          file_id,
          vars.t,
          shape.mass_nz,
          shape.mass_ny,
          shape.mass_nx,
          invalid_parent_pressure_thermodynamics ? -500.0F : 3'000.0F,
          5.0F,
          3.0F,
          10.0F);
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
    if (optional_state_fields) {
      put_3d_linear(file_id, vars.pb, shape.mass_nz, shape.mass_ny, shape.mass_nx, 80'000.0F, 1.0F, 10.0F, 100.0F);
      put_3d_linear(file_id, vars.phb, shape.full_nz, shape.mass_ny, shape.mass_nx, 90'000.0F, 1.0F, 10.0F, 100.0F);
      put_2d_linear(file_id, vars.mub, shape.mass_ny, shape.mass_nx, 100'000.0F, 1.0F, 10.0F);
      put_2d_linear(file_id, vars.psfc, shape.mass_ny, shape.mass_nx, 110'000.0F, 1.0F, 10.0F);
      put_2d_linear(file_id, vars.u10, shape.mass_ny, shape.mass_nx, 120'000.0F, 1.0F, 10.0F);
      put_2d_linear(file_id, vars.v10, shape.mass_ny, shape.mass_nx, 130'000.0F, 1.0F, 10.0F);
      put_2d_linear(file_id, vars.t2, shape.mass_ny, shape.mass_nx, 140'000.0F, 1.0F, 10.0F);
      put_2d_linear(file_id, vars.q2, shape.mass_ny, shape.mass_nx, 150'000.0F, 1.0F, 10.0F);
      put_2d_linear(file_id, vars.rainc, shape.mass_ny, shape.mass_nx, 160'000.0F, 1.0F, 10.0F);
      put_2d_linear(file_id, vars.rainnc, shape.mass_ny, shape.mass_nx, 170'000.0F, 1.0F, 10.0F);
    }
  }
  if (pressure_refresh_metadata) {
    put_time_vector(file_id, c3f, {1.0F, 0.5F, 0.0F});
    put_time_vector(file_id, c4f, {0.0F, 0.0F, 0.0F});
    put_time_vector(file_id, c3h, {0.75F, 0.25F});
    put_time_vector(file_id, c4h, {0.0F, 0.0F});
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

[[nodiscard]] std::string read_file(const std::filesystem::path& path) {
  std::ifstream stream(path);
  if (!stream) {
    throw std::runtime_error("failed to read file: " + path.string());
  }
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  return buffer.str();
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

[[nodiscard]] int read_int_attr(const int file_id, const std::string_view name) {
  int value = 0;
  check_nc(nc_get_att_int(file_id, NC_GLOBAL, std::string(name).c_str(), &value), "read attr");
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

[[nodiscard]] std::vector<float> read_float_variable(
    const int file_id,
    const std::string_view name) {
  int variable_id = -1;
  check_nc(nc_inq_varid(file_id, std::string(name).c_str(), &variable_id), "inquire float var");
  int ndim = 0;
  check_nc(nc_inq_varndims(file_id, variable_id, &ndim), "inquire float var ndims");
  std::vector<int> dim_ids(static_cast<std::size_t>(ndim));
  check_nc(nc_inq_vardimid(file_id, variable_id, dim_ids.data()), "inquire float var dims");

  std::size_t value_count = 1;
  for (const int dim_id : dim_ids) {
    std::size_t dim_len = 0;
    check_nc(nc_inq_dimlen(file_id, dim_id, &dim_len), "inquire float var dim len");
    value_count *= dim_len;
  }
  std::vector<float> values(value_count);
  check_nc(nc_get_var_float(file_id, variable_id, values.data()), "read float var");
  return values;
}

void expect_close(const float actual, const float expected, const std::string_view label) {
  if (std::fabs(actual - expected) > kTolerance) {
    std::cerr << label << " mismatch: got " << actual << ", expected " << expected << '\n';
    assert(false);
  }
}

void expect_float_variable_match(
    const int actual_file_id,
    const int expected_file_id,
    const std::string_view name,
    const std::string_view label) {
  const auto actual = read_float_variable(actual_file_id, name);
  const auto expected = read_float_variable(expected_file_id, name);
  assert(actual.size() == expected.size());
  for (std::size_t index = 0; index < actual.size(); ++index) {
    if (std::fabs(actual[index] - expected[index]) > kTolerance) {
      std::cerr << label << " " << name << " mismatch at flat index " << index << ": got "
                << actual[index] << ", expected " << expected[index] << '\n';
      assert(false);
    }
  }
}

[[nodiscard]] float max_abs_float_variable_diff(
    const int actual_file_id,
    const int expected_file_id,
    const std::string_view name) {
  const auto actual = read_float_variable(actual_file_id, name);
  const auto expected = read_float_variable(expected_file_id, name);
  assert(actual.size() == expected.size());
  float max_abs_diff = 0.0F;
  for (std::size_t index = 0; index < actual.size(); ++index) {
    max_abs_diff = std::max(max_abs_diff, std::fabs(actual[index] - expected[index]));
  }
  return max_abs_diff;
}

[[nodiscard]] bool contains_csv_value(const std::string& values, const std::string_view expected) {
  std::size_t begin = 0;
  while (begin <= values.size()) {
    const auto end = values.find(',', begin);
    const auto token = values.substr(begin, end == std::string::npos ? values.size() - begin : end - begin);
    if (token == expected) {
      return true;
    }
    if (end == std::string::npos) {
      break;
    }
    begin = end + 1;
  }
  return false;
}

void assert_contains_in_order(
    const std::string& value,
    const std::vector<std::string_view>& markers) {
  std::size_t search_begin = 0;
  for (const auto marker : markers) {
    const auto found = value.find(marker, search_begin);
    assert(found != std::string::npos);
    search_begin = found + marker.size();
  }
}

[[nodiscard]] std::string timeline_field_value(
    const std::string& events,
    const std::string_view event_name,
    const std::string_view field_name) {
  const std::string event_marker = ":" + std::string(event_name) + "(";
  const auto event_begin = events.find(event_marker);
  assert(event_begin != std::string::npos);
  const auto fields_begin = event_begin + event_marker.size();
  const auto event_end = events.find(')', fields_begin);
  assert(event_end != std::string::npos);
  const std::string field_marker = std::string(field_name) + "=";
  const auto field_begin = events.find(field_marker, fields_begin);
  assert(field_begin != std::string::npos);
  assert(field_begin < event_end);
  const auto value_begin = field_begin + field_marker.size();
  auto value_end = events.find(',', value_begin);
  if (value_end == std::string::npos || value_end > event_end) {
    value_end = event_end;
  }
  assert(value_end > value_begin);
  return events.substr(value_begin, value_end - value_begin);
}

[[nodiscard]] std::uint64_t timeline_u64_field(
    const std::string& events,
    const std::string_view event_name,
    const std::string_view field_name) {
  return static_cast<std::uint64_t>(
      std::stoull(timeline_field_value(events, event_name, field_name)));
}

void assert_selected_field_timeline_attrs(
    const int file_id,
    const bool pressure_refresh,
    const bool experimental_pressure_refresh,
    const bool pressure_column_probe = false) {
  const std::string expected_names =
      std::string(
          "cycle_start,move_from_to_parent_start,overlap_remap,exchange_plan_build,"
          "parent_interpolation,selected_field_change_summary,static_refresh,"
          "normal_base_state_producer,pressure_refresh_readiness,pressure_refresh_apply,") +
      (pressure_column_probe ? "pressure_column_probe," : "") +
      "cycle_end,output_write_preparation";
  assert(read_text_attr(file_id, "TYWRF_SELECTED_FIELD_TIMELINE_VERSION") == "runtime_v0");
  assert(read_text_attr(file_id, "TYWRF_SELECTED_FIELD_TIMELINE_EVIDENCE_ONLY") == "true");
  assert(
      read_int_attr(file_id, "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_COUNT") ==
      (pressure_column_probe ? 13 : 12));
  assert(read_text_attr(file_id, "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_NAMES") == expected_names);
  const auto events = read_text_attr(file_id, "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS");
  std::vector<std::string_view> expected_order = {
      ":cycle_start(",
      ":move_from_to_parent_start(",
      ":overlap_remap(",
      ":exchange_plan_build(",
      ":parent_interpolation(",
      ":selected_field_change_summary(",
      ":static_refresh(",
      ":normal_base_state_producer(",
      ":pressure_refresh_readiness(",
      ":pressure_refresh_apply(",
  };
  if (pressure_column_probe) {
    expected_order.push_back(":pressure_column_probe(");
  }
  expected_order.push_back(":cycle_end(");
  expected_order.push_back(":output_write_preparation(");
  assert_contains_in_order(events, expected_order);
  assert(timeline_field_value(events, "cycle_start", "cycle_start") == kCycleStart);
  assert(timeline_field_value(events, "cycle_start", "cycle_end") == kTenMinuteEnd);
  assert(timeline_field_value(events, "cycle_end", "cycle_start") == kCycleStart);
  assert(timeline_field_value(events, "cycle_end", "cycle_end") == kTenMinuteEnd);
  assert(timeline_field_value(events, "move_from_to_parent_start", "from_i") == "2");
  assert(timeline_field_value(events, "move_from_to_parent_start", "from_j") == "2");
  assert(timeline_field_value(events, "move_from_to_parent_start", "to_i") == "3");
  assert(timeline_field_value(events, "move_from_to_parent_start", "to_j") == "2");
  assert(timeline_field_value(events, "move_from_to_parent_start", "parent_grid_ratio") == "5");
  assert(timeline_field_value(events, "move_from_to_parent_start", "child_delta_i") == "5");
  assert(timeline_field_value(events, "move_from_to_parent_start", "child_delta_j") == "0");
  assert(
      timeline_u64_field(events, "selected_field_change_summary", "changed_points") ==
      static_cast<std::uint64_t>(
          read_double_attr(file_id, "TYWRF_SELECTED_FIELD_CHANGED_POINTS")));
  assert(
      timeline_u64_field(events, "exchange_plan_build", "exchange_points") ==
      static_cast<std::uint64_t>(read_double_attr(file_id, "TYWRF_EXPOSED_EXCHANGE_POINTS")));
  assert(
      timeline_u64_field(events, "parent_interpolation", "interpolated_points") ==
      static_cast<std::uint64_t>(read_double_attr(file_id, "TYWRF_INTERPOLATED_POINTS")));
  assert(
      timeline_u64_field(events, "static_refresh", "overlap_cells") ==
      static_cast<std::uint64_t>(read_double_attr(file_id, "TYWRF_STATIC_REFRESH_OVERLAP_CELLS")));
  assert(
      timeline_u64_field(events, "static_refresh", "exposed_cells") ==
      static_cast<std::uint64_t>(read_double_attr(file_id, "TYWRF_STATIC_REFRESH_EXPOSED_CELLS")));
  assert(
      timeline_u64_field(events, "static_refresh", "hgt_parent_interpolated_cells") ==
      static_cast<std::uint64_t>(
          read_double_attr(file_id, "TYWRF_STATIC_REFRESH_HGT_PARENT_INTERPOLATED_CELLS")));
  assert(timeline_field_value(events, "static_refresh", "uses_reference_end") == "false");
  assert(timeline_field_value(events, "normal_base_state_producer", "source") ==
         "normal_selected_field_base_state_producer");
  assert(
      timeline_field_value(events, "normal_base_state_producer", "source_origin") ==
      "krosa_base_state_provider+moved_candidate_HGT+same_cycle_vertical_metadata");
  assert(timeline_field_value(events, "normal_base_state_producer", "diagnostic_only") == "false");
  assert(
      timeline_field_value(events, "normal_base_state_producer", "normal_candidate_producer") ==
      "true");
  assert(timeline_field_value(events, "normal_base_state_producer", "writes_candidate") == "true");
  assert(
      timeline_field_value(events, "normal_base_state_producer", "uses_reference_end_truth") ==
      "false");
  assert(
      timeline_field_value(events, "normal_base_state_producer", "uses_direct_p_shortcut") ==
      "false");
  assert(timeline_field_value(events, "normal_base_state_producer", "reads_direct_p") == "false");
  assert(timeline_field_value(events, "normal_base_state_producer", "writes_p") == "false");
  assert(
      timeline_field_value(events, "normal_base_state_producer", "no_gate_pass_claim") ==
      "true");
  assert(
      timeline_field_value(events, "normal_base_state_producer", "no_00_20_progression") ==
      "true");
  assert(
      timeline_field_value(events, "normal_base_state_producer", "terrain_source") ==
      "moved_candidate_HGT");
  assert(
      timeline_field_value(events, "normal_base_state_producer", "terrain_provenance") ==
      "override:moved_candidate_HGT");
  assert(timeline_u64_field(events, "normal_base_state_producer", "changed_p_points") == 0);
  assert(timeline_u64_field(events, "normal_base_state_producer", "changed_pb_points") > 0);
  assert(timeline_u64_field(events, "normal_base_state_producer", "changed_mub_points") > 0);
  assert(timeline_u64_field(events, "normal_base_state_producer", "changed_phb_points") > 0);
  assert(timeline_field_value(events, "output_write_preparation", "times") == kTenMinuteEnd);
  assert(timeline_u64_field(events, "output_write_preparation", "variable_count") > 0);
  if (pressure_refresh) {
    if (!pressure_column_probe) {
      assert(events.find("helper") == std::string::npos);
      assert(events.find("dry_run") == std::string::npos);
      assert(events.find("staging") == std::string::npos);
      assert(events.find("experimental") == std::string::npos);
    }
    assert(timeline_field_value(events, "pressure_refresh_readiness", "opt_in") == "true");
    assert(timeline_field_value(events, "pressure_refresh_readiness", "ready") == "true");
    assert(
        timeline_field_value(
            events, "pressure_refresh_readiness", "base_state_source_sync_readiness_check") ==
        "true");
    assert(
        timeline_u64_field(
            events, "pressure_refresh_readiness", "source_sync_planned_pb_points") > 0);
    assert(
        timeline_u64_field(
            events, "pressure_refresh_readiness", "source_sync_planned_mub_points") > 0);
    assert(
        timeline_u64_field(
            events, "pressure_refresh_readiness", "source_sync_planned_phb_points") > 0);
    assert(
        timeline_u64_field(
            events, "pressure_refresh_readiness", "production_refresh_planned_p_points") > 0);
    assert(
        timeline_u64_field(events, "pressure_refresh_readiness", "readiness_invalid_p_points") ==
        0);
    assert(
        timeline_u64_field(events, "pressure_refresh_readiness", "readiness_skipped_p_points") ==
        0);
    assert(timeline_field_value(events, "pressure_refresh_apply", "opt_in") == "true");
    assert(timeline_field_value(events, "pressure_refresh_apply", "applied") == "true");
    assert(
        timeline_u64_field(events, "pressure_refresh_apply", "refreshed_points") ==
        static_cast<std::uint64_t>(
            read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT")));
    assert(
        timeline_u64_field(events, "pressure_refresh_apply", "synced_pb_points") ==
        static_cast<std::uint64_t>(
            read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS")));
    assert(
        timeline_u64_field(events, "pressure_refresh_apply", "synced_mub_points") ==
        static_cast<std::uint64_t>(
            read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS")));
    assert(
        timeline_u64_field(events, "pressure_refresh_apply", "synced_phb_points") ==
        static_cast<std::uint64_t>(
            read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS")));
    assert(
        timeline_u64_field(events, "pressure_refresh_apply", "changed_p_points") ==
        static_cast<std::uint64_t>(
            read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS")));
  } else {
    assert(timeline_field_value(events, "pressure_refresh_readiness", "opt_in") == "false");
    assert(timeline_field_value(events, "pressure_refresh_readiness", "status") == "skipped");
    assert(timeline_field_value(events, "pressure_refresh_apply", "opt_in") == "false");
    assert(timeline_field_value(events, "pressure_refresh_apply", "applied") == "false");
    assert(timeline_field_value(events, "pressure_refresh_apply", "status") == "skipped");
  }
  if (pressure_column_probe) {
    assert(timeline_field_value(events, "pressure_column_probe", "enabled") == "true");
    assert(timeline_field_value(events, "pressure_column_probe", "evidence_only") == "true");
    assert(timeline_u64_field(events, "pressure_column_probe", "record_count") > 0);
  }
  assert(
      timeline_field_value(events, "output_write_preparation", "metadata_write") == "pending");
  (void)experimental_pressure_refresh;
}

void assert_pressure_refresh_readiness_attrs(const int file_id) {
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_READINESS_READY") == "true");
  assert(
      read_text_attr(file_id, "TYWRF_THERMODYNAMIC_BASE_STATE_CONSISTENCY_READY") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT") ==
      "true");
  assert(read_text_attr(file_id, "TYWRF_PROVIDER_BASE_STATE_RECONSTRUCT_OK") == "true");
  assert(read_text_attr(file_id, "TYWRF_BASE_STATE_SYNC_CONTRACT_OK") == "true");
  assert(read_text_attr(file_id, "TYWRF_BASE_STATE_SOURCE_SYNC_READINESS_CHECK") == "true");
  assert(read_text_attr(file_id, "TYWRF_BASE_STATE_SOURCE_SYNC_APPLIED") == "false");
  assert(read_double_attr(file_id, "TYWRF_SOURCE_SYNC_PLANNED_PB_POINT_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_SOURCE_SYNC_PLANNED_MUB_POINT_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_SOURCE_SYNC_PLANNED_PHB_POINT_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_SOURCE_SYNC_OVERLAP_WRITE_COUNT") == 0.0);
  assert(read_double_attr(file_id, "TYWRF_SOURCE_SYNC_HALO_WRITE_COUNT") == 0.0);
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_READINESS_COMPUTE_CALLED") == "false");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK") == "true");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK_CALLED") == "true");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK_OK") == "true");
  const auto would_refresh_p_points =
      read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_PLANNED_P_POINT_COUNT");
  const auto dry_run_refreshed_points =
      read_double_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_REFRESHED_POINT_COUNT");
  assert(would_refresh_p_points > 0.0);
  assert(dry_run_refreshed_points == would_refresh_p_points);
  assert(read_double_attr(file_id, "TYWRF_PRESSURE_READINESS_INVALID_P_POINT_COUNT") == 0.0);
  assert(read_double_attr(file_id, "TYWRF_PRESSURE_READINESS_SKIPPED_P_POINT_COUNT") == 0.0);
  assert(
      read_double_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_TARGET_COLUMN_COUNT") >
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_INVALID_POINT_COUNT") ==
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_SKIPPED_POINT_COUNT") ==
      0.0);
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_TOUCHED_OVERLAP_CELLS") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_TOUCHED_HALO_CELLS") ==
      "false");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_READINESS_APPLIED") == "false");
  assert(!has_text_attr(file_id, "TYWRF_BASE_STATE_SYNC_DRY_RUN"));
  assert(!has_text_attr(file_id, "TYWRF_PRESSURE_COMPUTE_DRY_RUN"));
  assert(!has_text_attr(file_id, "TYWRF_DRY_RUN_INVALID_P_POINT_COUNT"));
  assert(!has_text_attr(file_id, "TYWRF_DRY_RUN_SKIPPED_P_POINT_COUNT"));
  assert(!has_text_attr(file_id, "TYWRF_PRESSURE_COMPUTE_DRY_RUN_REPORT_TARGET_COLUMN_COUNT"));
  assert(!has_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_DRY_RUN_APPLIED"));
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_SOURCE") ==
      "moved_candidate_HGT");
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_PROVENANCE") ==
      "override:moved_candidate_HGT");
}

void assert_normal_base_state_producer_attrs(const int file_id) {
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_VERSION") ==
      "a79_normal_v0");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_SOURCE") ==
      "normal_selected_field_base_state_producer");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_SOURCE_ORIGIN") ==
      "krosa_base_state_provider+moved_candidate_HGT+same_cycle_vertical_metadata");
  assert(read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_OK") == "true");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_DIAGNOSTIC_ONLY") ==
      "false");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_NORMAL_BASE_STATE_PRODUCER_NORMAL_CANDIDATE_PRODUCER") == "true");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_GATE_CANDIDATE") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_INTEGRATOR_OUTPUT") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_WRITES_CANDIDATE") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_WRITES_NETCDF") ==
      "false");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_NORMAL_BASE_STATE_PRODUCER_USES_REFERENCE_END_TRUTH") == "false");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_NORMAL_BASE_STATE_PRODUCER_NO_REFERENCE_END_TRUTH") == "true");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_NORMAL_BASE_STATE_PRODUCER_USES_DIRECT_P_SHORTCUT") == "false");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_NORMAL_BASE_STATE_PRODUCER_NO_DIRECT_P_SHORTCUT") == "true");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_READS_DIRECT_P") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_WRITES_P") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_NO_GATE_PASS_CLAIM") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_NO_00_20_PROGRESSION") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_WRITTEN_FIELDS") ==
      "PHB,MUB,PB,HGT");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_PROVIDER_SOURCE") ==
      "base_state_reconstruction_provider");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_TERRAIN_SOURCE") ==
      "moved_candidate_HGT");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_TERRAIN_PROVENANCE") ==
      "override:moved_candidate_HGT");
  assert(
      read_text_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_TERRAIN_OVERRIDE_USED") ==
      "true");
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_EXPOSED_MASS_CELL_COUNT") >
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_PB_WRITTEN_POINT_COUNT") >
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_MUB_WRITTEN_CELL_COUNT") >
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_PHB_WRITTEN_POINT_COUNT") >
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_HT_WRITTEN_CELL_COUNT") >
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_P_POINTS") ==
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_PB_POINTS") >
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_MUB_POINTS") >
      0.0);
  assert(
      read_double_attr(file_id, "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_PHB_POINTS") >
      0.0);
}

void assert_no_pressure_formula_observation_attrs(const int file_id) {
  assert(!has_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_ENABLED"));
  assert(!has_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALUES"));
}

void assert_no_wind_tendency_attrs(const int file_id) {
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_OPT_IN"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_APPLIED"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_SOURCE_KIND"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_GATE_EVIDENCE"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_USES_REFERENCE_END_TRUTH"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ZERO_OR_IDENTITY_ONLY"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_WRITTEN_FIELDS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_STATUS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_VELOCITY_MODE"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_SUBSTEP_COUNT"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_SUBSTEP_DT_SECONDS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_TOTAL_SECONDS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ACTIVE_U_POINTS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ACTIVE_V_POINTS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_UPDATED_U_POINTS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_UPDATED_V_POINTS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_CHANGED_U_POINTS"));
  assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_CHANGED_V_POINTS"));
}

void assert_successful_candidate(
    const std::filesystem::path& output,
    const bool pressure_refresh = false,
    const std::filesystem::path& pressure_refresh_metadata_source = {},
    const bool experimental_pressure_refresh = false,
    const bool pressure_column_probe = false) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open output");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ONLY") ==
         (experimental_pressure_refresh ? "true" : "false"));
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") ==
         (experimental_pressure_refresh ? "false" : "true"));
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") ==
         (experimental_pressure_refresh ? "false" : "true"));
  assert(read_text_attr(file_id, "TYWRF_VALIDATION_GATE_ONLY") == "false");
  assert(
      read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") ==
      (experimental_pressure_refresh ? "selected_field_pressure_refresh_experimental_apply_v0"
                                     : "selected_field_integrator_v0"));
  if (experimental_pressure_refresh) {
    assert(read_text_attr(file_id, "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY") == "true");
  } else {
    assert(!has_text_attr(file_id, "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY"));
  }
  assert(read_text_attr(file_id, "TYWRF_D02_RESOLUTION_CHECK") == "d02_2km");
  assert(!has_text_attr(file_id, "TYWRF_NOT_PHYSICAL"));
  assert(read_double_attr(file_id, "DX") == 2000.0);
  assert(read_double_attr(file_id, "DY") == 2000.0);
  assert(read_double_attr(file_id, "TYWRF_PARENT_GRID_RATIO") == 5.0);
  assert(read_double_attr(file_id, "TYWRF_SELECTED_FIELD_CHANGED_POINTS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_INTERPOLATED_POINTS") > 0.0);
  assert(read_int_attr(file_id, "I_PARENT_START") == 3);
  assert(read_int_attr(file_id, "J_PARENT_START") == 2);
  assert(std::isfinite(read_double_attr(file_id, "CEN_LAT")));
  assert(std::isfinite(read_double_attr(file_id, "CEN_LON")));
  assert(read_text_attr(file_id, "TYWRF_STATIC_REFRESH_APPLIED") == "true");
  assert(read_text_attr(file_id, "TYWRF_STATIC_REFRESH_USES_REFERENCE_END") == "false");
  assert(read_text_attr(file_id, "TYWRF_STATIC_REFRESH_D02_START_SOURCE").find("wrfout_d02") != std::string::npos);
  assert(read_text_attr(file_id, "TYWRF_STATIC_REFRESH_D01_HGT_SOURCE").find("wrfout_d01") != std::string::npos);
  assert(read_double_attr(file_id, "TYWRF_STATIC_REFRESH_OVERLAP_CELLS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_STATIC_REFRESH_EXPOSED_CELLS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_STATIC_REFRESH_COORD_EXTRAPOLATED_CELLS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_STATIC_REFRESH_HGT_PARENT_INTERPOLATED_CELLS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_STATIC_REFRESH_CHANGED_TEMPLATE_POINTS") > 0.0);
  assert_normal_base_state_producer_attrs(file_id);
  if (pressure_refresh) {
    assert_pressure_refresh_readiness_attrs(file_id);
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_OPT_IN") == "true");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_APPLIED") == "true");
    assert(
        read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS") ==
        (experimental_pressure_refresh ? "experimental_apply_test_only" : "applied_to_candidate"));
    if (experimental_pressure_refresh) {
      assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY") == "true");
    } else {
      assert(!has_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY"));
    }
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_PROVIDER_OK") == "true");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_SOURCE_SYNC_OK") == "true");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED") == "true");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_TERRAIN_OVERRIDE_USED") == "true");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE") ==
           "moved_candidate_HGT");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE") ==
           "override:moved_candidate_HGT");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE") == pressure_refresh_metadata_source.string());
    assert(
        read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_PRODUCTION_SOURCE") ==
        "krosa_moving_nest_pressure_refresh");
    assert(!has_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_STAGING_OK"));
    assert(!has_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_HELPER_NAME"));
    assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS") > 0.0);
    assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS") > 0.0);
    assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS") > 0.0);
    const auto target_columns =
        read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT");
    const auto refreshed_columns =
        read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT");
    const auto refreshed_points =
        read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT");
    const auto refreshed_p_points =
        read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS");
    const auto changed_p_points =
        read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS");
    assert(target_columns > 0.0);
    assert(refreshed_columns == target_columns);
    assert(refreshed_points > 0.0);
    assert(refreshed_p_points == refreshed_points);
    assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT") == 0.0);
    assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT") == 0.0);
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS") == "false");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_TOUCHED_HALO_CELLS") == "false");
    assert(changed_p_points == refreshed_points);
    assert(changed_p_points > 0.0);
    assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS") >= 0.0);
    assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS") >= 0.0);
    assert(read_double_attr(file_id, "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS") >= 0.0);
    assert(
        read_text_attr(
            file_id,
            "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT") == "true");
    assert(
        read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO") ==
        "true");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED") == "true");
  } else {
    assert(!has_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_OPT_IN"));
    assert(!has_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_APPLIED"));
  }
  if (pressure_column_probe) {
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_VERSION") == "runtime_v0");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_ENABLED") == "true");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_EVIDENCE_ONLY") == "true");
  } else {
    assert(!has_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_ENABLED"));
  }
  if (!(pressure_refresh && pressure_column_probe)) {
    assert_no_pressure_formula_observation_attrs(file_id);
  }
  assert_no_wind_tendency_attrs(file_id);
  assert_selected_field_timeline_attrs(
      file_id, pressure_refresh, experimental_pressure_refresh, pressure_column_probe);
  const auto state_variables = read_text_attr(file_id, "TYWRF_STATE_VARIABLES");
  assert(contains_csv_value(state_variables, "PB"));
  assert(contains_csv_value(state_variables, "PHB"));
  assert(contains_csv_value(state_variables, "MUB"));
  assert(contains_csv_value(state_variables, "U10"));
  assert(contains_csv_value(state_variables, "V10"));
  assert(contains_csv_value(state_variables, "RAINC"));
  assert(contains_csv_value(state_variables, "RAINNC"));
  assert(!contains_csv_value(state_variables, "SLP"));
  const auto interpolated_variables =
      read_text_attr(file_id, "TYWRF_PARENT_INTERPOLATED_STATE_VARIABLES");
  assert(contains_csv_value(interpolated_variables, "U"));
  assert(contains_csv_value(interpolated_variables, "V"));
  assert(contains_csv_value(interpolated_variables, "T"));
  assert(contains_csv_value(interpolated_variables, "PH"));
  assert(contains_csv_value(interpolated_variables, "MU"));
  assert(contains_csv_value(interpolated_variables, "QVAPOR"));
  assert(!contains_csv_value(interpolated_variables, "P"));
  const auto candidate_message = read_text_attr(file_id, "TYWRF_CANDIDATE_MESSAGE");
  if (experimental_pressure_refresh) {
    assert(candidate_message.find("Experimental selected-field pressure-refresh apply seam") !=
           std::string::npos);
    assert(candidate_message.find("not a validation gate pass or normal integrator output") !=
           std::string::npos);
    assert(candidate_message.find("moved_candidate_HGT terrain override") != std::string::npos);
    assert(candidate_message.find("refreshed exposed P/PB/PHB/MUB") != std::string::npos);
  } else if (pressure_refresh) {
    assert(candidate_message.find("pressure refresh applied to exposed P/PB/PHB/MUB") !=
           std::string::npos);
  } else {
    assert(candidate_message.find("U/V/T/PH/MU/QVAPOR exposed cells are parent interpolated") !=
           std::string::npos);
    assert(candidate_message.find("exposed PB/PHB/MUB are updated by the normal") !=
           std::string::npos);
    assert(candidate_message.find("P remains finite d02 start-state ownership") !=
           std::string::npos);
  }
  assert(read_times(file_id).substr(0, 19) == kTenMinuteEnd);

  expect_close(read_2d_value(file_id, "XLAT", 2, 0), linear_2d(2, 5, 40'000.0F, 1.0F, 10.0F), "XLAT overlap");
  expect_close(read_2d_value(file_id, "XLONG", 2, 0), linear_2d(2, 5, 50'000.0F, 1.0F, 10.0F), "XLONG overlap");
  expect_close(
      read_2d_value(file_id, "HGT", 2, 0),
      linear_2d(2, 5, kSyntheticHgtBase, 1.0F, 10.0F),
      "HGT overlap");
  const auto exposed_xlat = read_2d_value(file_id, "XLAT", 0, 9);
  const auto stale_xlat = linear_2d(0, 9, 40'000.0F, 1.0F, 10.0F);
  assert(std::isfinite(exposed_xlat));
  assert(std::fabs(exposed_xlat - stale_xlat) > kTolerance);
  expect_close(exposed_xlat, linear_2d(0, 14, 40'000.0F, 1.0F, 10.0F), "XLAT exposed");
  expect_close(read_2d_value(file_id, "XLONG", 0, 9), linear_2d(0, 14, 50'000.0F, 1.0F, 10.0F), "XLONG exposed");
  expect_close(read_2d_value(file_id, "HGT", 0, 9), 613.8F, "HGT parent exposed");
  if (experimental_pressure_refresh) {
    const auto moved_candidate_hgt = read_2d_value(file_id, "HGT", 0, 9);
    const auto metadata_hgt = linear_2d(0, 9, 20.0F, 1.0F, 5.0F);
    assert(std::fabs(moved_candidate_hgt - metadata_hgt) > kTolerance);
  }

  expect_close(read_3d_value(file_id, "U", 1, 2, 0), linear_3d(1, 2, 5, 10'000.0F, 1.0F, 10.0F, 100.0F), "U overlap");
  expect_close(read_3d_value(file_id, "V", 1, 2, 0), linear_3d(1, 2, 5, 20'000.0F, 1.0F, 10.0F, 100.0F), "V overlap");
  expect_close(read_2d_value(file_id, "MU", 2, 0), linear_2d(2, 5, 50'000.0F, 1.0F, 10.0F), "MU overlap");
  expect_close(read_3d_value(file_id, "QVAPOR", 1, 2, 0), linear_3d(1, 2, 5, 70'000.0F, 1.0F, 10.0F, 100.0F), "QVAPOR overlap");

  expect_close(read_3d_value(file_id, "U", 1, 0, 10), parent_linear(1'000.0F, 7.0F, 13.0F, 101.0F, 4.0, 0.6, 1), "U exposed");
  expect_close(read_3d_value(file_id, "V", 1, 10, 9), parent_linear(2'000.0F, 11.0F, 17.0F, 103.0F, 3.4, 3.0, 1), "V exposed");
  expect_close(read_2d_value(file_id, "MU", 0, 9), parent_linear(5'000.0F, 19.0F, 23.0F, 0.0F, 3.4, 0.6, 0), "MU exposed");
  expect_close(read_3d_value(file_id, "QVAPOR", 1, 9, 9), parent_linear(7'000.0F, 29.0F, 31.0F, 107.0F, 3.4, 2.4, 1), "QVAPOR exposed");

  expect_close(read_3d_value(file_id, "T", 1, 9, 9), parent_linear(3'000.0F, 5.0F, 3.0F, 10.0F, 3.4, 2.4, 1), "T exposed");
  expect_close(read_3d_value(file_id, "PH", 2, 9, 9), parent_linear(4'000.0F, 2.0F, 4.0F, 50.0F, 3.4, 2.4, 2), "PH exposed");
  const auto preserved_exposed_p = linear_3d(1, 9, 9, 60'000.0F, 1.0F, 10.0F, 100.0F);
  if (pressure_refresh) {
    const auto refreshed_p = read_3d_value(file_id, "P", 1, 9, 9);
    assert(std::isfinite(refreshed_p));
    assert(std::fabs(refreshed_p - preserved_exposed_p) > kTolerance);
  } else {
    expect_close(read_3d_value(file_id, "P", 1, 9, 9), preserved_exposed_p, "P preserved exposed");
  }
  const auto preserved_pb = linear_3d(1, 9, 9, 80'000.0F, 1.0F, 10.0F, 100.0F);
  const auto preserved_phb = linear_3d(2, 9, 9, 90'000.0F, 1.0F, 10.0F, 100.0F);
  const auto preserved_mub = linear_2d(9, 9, 100'000.0F, 1.0F, 10.0F);
  if (pressure_refresh) {
    assert(std::isfinite(read_3d_value(file_id, "PB", 1, 9, 9)));
    assert(std::fabs(read_3d_value(file_id, "PB", 1, 9, 9) - preserved_pb) > kTolerance);
    assert(std::isfinite(read_3d_value(file_id, "PHB", 2, 9, 9)));
    assert(std::fabs(read_3d_value(file_id, "PHB", 2, 9, 9) - preserved_phb) > kTolerance);
    assert(std::isfinite(read_2d_value(file_id, "MUB", 9, 9)));
    assert(std::fabs(read_2d_value(file_id, "MUB", 9, 9) - preserved_mub) > kTolerance);
  } else {
    assert(std::isfinite(read_3d_value(file_id, "PB", 1, 9, 9)));
    assert(std::fabs(read_3d_value(file_id, "PB", 1, 9, 9) - preserved_pb) > kTolerance);
    assert(std::isfinite(read_3d_value(file_id, "PHB", 2, 9, 9)));
    assert(std::fabs(read_3d_value(file_id, "PHB", 2, 9, 9) - preserved_phb) > kTolerance);
    assert(std::isfinite(read_2d_value(file_id, "MUB", 9, 9)));
    assert(std::fabs(read_2d_value(file_id, "MUB", 9, 9) - preserved_mub) > kTolerance);
  }
  expect_close(read_2d_value(file_id, "PSFC", 9, 9), linear_2d(9, 9, 110'000.0F, 1.0F, 10.0F), "PSFC preserved");
  expect_close(read_2d_value(file_id, "U10", 9, 9), linear_2d(9, 9, 120'000.0F, 1.0F, 10.0F), "U10 preserved");
  expect_close(read_2d_value(file_id, "V10", 9, 9), linear_2d(9, 9, 130'000.0F, 1.0F, 10.0F), "V10 preserved");
  expect_close(read_2d_value(file_id, "T2", 9, 9), linear_2d(9, 9, 140'000.0F, 1.0F, 10.0F), "T2 preserved");
  expect_close(read_2d_value(file_id, "Q2", 9, 9), linear_2d(9, 9, 150'000.0F, 1.0F, 10.0F), "Q2 preserved");
  expect_close(read_2d_value(file_id, "RAINC", 9, 9), linear_2d(9, 9, 160'000.0F, 1.0F, 10.0F), "RAINC preserved");
  expect_close(read_2d_value(file_id, "RAINNC", 9, 9), linear_2d(9, 9, 170'000.0F, 1.0F, 10.0F), "RAINNC preserved");
  check_nc(nc_close(file_id), "close output");
}

void assert_pressure_column_probe_output(
    const std::filesystem::path& output,
    const std::filesystem::path& log_path,
    const bool pressure_refresh) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open probe output");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_VERSION") == "runtime_v0");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_ENABLED") == "true");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_EVIDENCE_ONLY") == "true");
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_INDEX_BASE") ==
      "zero_based_mass_grid");
  assert(read_int_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_COLUMN_COUNT") == 1);
  assert(read_int_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_LEVEL_COUNT") == 2);
  assert(read_int_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_PHASE_COUNT") == 2);
  assert(read_int_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_RECORD_COUNT") == 4);
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_COLUMNS") == "9,9");
  assert(read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_LEVELS") == "0,1");
  assert(
      read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_PHASES") ==
      (pressure_refresh ? "post_static_refresh,post_pressure_refresh"
                        : "post_static_refresh,pressure_refresh_skipped"));
  const auto fields = read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_FIELDS");
  assert(contains_csv_value(fields, "P"));
  assert(contains_csv_value(fields, "PB"));
  assert(contains_csv_value(fields, "P+PB"));
  assert(contains_csv_value(fields, "MU"));
  assert(contains_csv_value(fields, "MUB"));
  assert(contains_csv_value(fields, "MU+MUB"));
  assert(contains_csv_value(fields, "PH"));
  assert(contains_csv_value(fields, "PHB"));
  assert(contains_csv_value(fields, "PH+PHB"));
  assert(contains_csv_value(fields, "T"));
  assert(contains_csv_value(fields, "QVAPOR"));
  assert(contains_csv_value(fields, "HGT"));
  const auto not_available =
      read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_NOT_AVAILABLE");
  assert(contains_csv_value(not_available, "ALB"));
  assert(contains_csv_value(not_available, "C3F"));
  assert(contains_csv_value(not_available, "C4F"));
  assert(contains_csv_value(not_available, "C3H"));
  assert(contains_csv_value(not_available, "C4H"));
  assert(contains_csv_value(not_available, "P_TOP"));
  assert(contains_csv_value(not_available, "theta_m"));
  const auto values = read_text_attr(file_id, "TYWRF_PRESSURE_COLUMN_PROBE_VALUES");
  assert(values.find("phase=post_static_refresh;i=9;j=9;k=0") != std::string::npos);
  assert(values.find("phase=post_static_refresh;i=9;j=9;k=1") != std::string::npos);
  assert(
      values.find(pressure_refresh ? "phase=post_pressure_refresh;i=9;j=9;k=1"
                                   : "phase=pressure_refresh_skipped;i=9;j=9;k=1") !=
      std::string::npos);
  assert(values.find(";P=") != std::string::npos);
  assert(values.find(";PB=") != std::string::npos);
  assert(values.find(";P_PLUS_PB=") != std::string::npos);
  assert(values.find(";MU=") != std::string::npos);
  assert(values.find(";MUB=") != std::string::npos);
  assert(values.find(";MU_PLUS_MUB=") != std::string::npos);
  assert(values.find(";PH=") != std::string::npos);
  assert(values.find(";PHB=") != std::string::npos);
  assert(values.find(";PH_PLUS_PHB=") != std::string::npos);
  assert(values.find(";T=") != std::string::npos);
  assert(values.find(";QVAPOR=") != std::string::npos);
  assert(values.find(";HGT=") != std::string::npos);
  if (pressure_refresh) {
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_VERSION") == "runtime_v0");
    assert(read_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_ENABLED") == "true");
    assert(
        read_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_EVIDENCE_ONLY") ==
        "true");
    assert(
        read_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_INDEX_BASE") ==
        "zero_based_mass_grid");
    assert(read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_REQUEST_COUNT") == 2);
    assert(read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_RECORD_COUNT") == 2);
    assert(read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALID_COUNT") == 2);
    assert(read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_INVALID_COUNT") == 0);
    assert(
        read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_OUT_OF_BOUNDS_COUNT") ==
        0);
    assert(
        read_int_attr(
            file_id,
            "TYWRF_PRESSURE_FORMULA_OBSERVATION_OUTSIDE_TARGET_REGION_COUNT") == 0);
    const auto formula_fields =
        read_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_FIELDS");
    assert(contains_csv_value(formula_fields, "status"));
    assert(contains_csv_value(formula_fields, "valid"));
    assert(contains_csv_value(formula_fields, "mu_total"));
    assert(contains_csv_value(formula_fields, "ALB"));
    assert(contains_csv_value(formula_fields, "PB"));
    assert(contains_csv_value(formula_fields, "theta"));
    assert(contains_csv_value(formula_fields, "alpha_from_wrf_branch"));
    assert(contains_csv_value(formula_fields, "perturbation_pressure_pa"));
    const auto formula_values =
        read_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALUES");
    assert(formula_values.find("status=recorded;valid=1;i=9;j=9;k=0") !=
           std::string::npos);
    assert(formula_values.find("status=recorded;valid=1;i=9;j=9;k=1") !=
           std::string::npos);
    assert(formula_values.find(";mu_total=") != std::string::npos);
    assert(formula_values.find(";pfu=") != std::string::npos);
    assert(formula_values.find(";pfd=") != std::string::npos);
    assert(formula_values.find(";phm=") != std::string::npos);
    assert(formula_values.find(";log_ratio=") != std::string::npos);
    assert(formula_values.find(";phi_lower=") != std::string::npos);
    assert(formula_values.find(";phi_upper=") != std::string::npos);
    assert(formula_values.find(";delta_phi=") != std::string::npos);
    assert(formula_values.find(";ALB=") != std::string::npos);
    assert(formula_values.find(";PB=") != std::string::npos);
    assert(formula_values.find(";theta=") != std::string::npos);
    assert(formula_values.find(";alpha_total=") != std::string::npos);
    assert(formula_values.find(";alpha_perturbation=") != std::string::npos);
    assert(formula_values.find(";alpha_from_wrf_branch=") != std::string::npos);
    assert(formula_values.find(";pressure_base=") != std::string::npos);
    assert(formula_values.find(";total_pressure=") != std::string::npos);
    assert(formula_values.find(";perturbation_pressure_pa=") != std::string::npos);
  } else {
    assert_no_pressure_formula_observation_attrs(file_id);
  }
  check_nc(nc_close(file_id), "close probe output");

  const auto log = read_file(log_path);
  assert(log.find("\"pressure_column_probe_enabled\": true") != std::string::npos);
  assert(log.find("\"pressure_column_probe_version\": \"runtime_v0\"") != std::string::npos);
  assert(
      log.find("\"pressure_column_probe_index_base\": \"zero_based_mass_grid\"") !=
      std::string::npos);
  assert(log.find("\"pressure_column_probe_column_count\": 1") != std::string::npos);
  assert(log.find("\"pressure_column_probe_level_count\": 2") != std::string::npos);
  assert(log.find("\"pressure_column_probe_record_count\": 4") != std::string::npos);
  assert(log.find("\"pressure_column_probe_columns\": \"9,9\"") != std::string::npos);
  assert(log.find("\"pressure_column_probe_levels\": \"0,1\"") != std::string::npos);
  assert(log.find("\"pressure_column_probe_not_available\":") != std::string::npos);
  assert(log.find("P_PLUS_PB=") != std::string::npos);
  if (pressure_refresh) {
    assert(
        log.find("\"pressure_formula_observation_enabled\": true") !=
        std::string::npos);
    assert(
        log.find("\"pressure_formula_observation_version\": \"runtime_v0\"") !=
        std::string::npos);
    assert(
        log.find("\"pressure_formula_observation_evidence_only\": true") !=
        std::string::npos);
    assert(
        log.find("\"pressure_formula_observation_index_base\": \"zero_based_mass_grid\"") !=
        std::string::npos);
    assert(
        log.find("\"pressure_formula_observation_request_count\": 2") !=
        std::string::npos);
    assert(
        log.find("\"pressure_formula_observation_record_count\": 2") !=
        std::string::npos);
    assert(
        log.find("\"pressure_formula_observation_valid_count\": 2") !=
        std::string::npos);
    assert(log.find("status=recorded;valid=1;i=9;j=9;k=0") != std::string::npos);
    assert(log.find("alpha_from_wrf_branch=") != std::string::npos);
    assert(log.find("perturbation_pressure_pa=") != std::string::npos);
  } else {
    assert(log.find("pressure_formula_observation_enabled") == std::string::npos);
  }
}

[[nodiscard]] std::uint64_t log_u64_field(
    const std::string& log,
    const std::string_view name) {
  const std::string marker = std::string(name) + "=";
  const auto value_begin = log.find(marker);
  assert(value_begin != std::string::npos);
  std::size_t digit_begin = value_begin + marker.size();
  std::size_t digit_end = digit_begin;
  while (digit_end < log.size() && std::isdigit(static_cast<unsigned char>(log[digit_end]))) {
    ++digit_end;
  }
  assert(digit_end > digit_begin);
  return static_cast<std::uint64_t>(std::stoull(log.substr(digit_begin, digit_end - digit_begin)));
}

[[nodiscard]] double json_number_field(
    const std::string& json,
    const std::string_view name) {
  const std::string marker = "\"" + std::string(name) + "\":";
  const auto value_begin = json.find(marker);
  assert(value_begin != std::string::npos);
  std::size_t digit_begin = value_begin + marker.size();
  while (digit_begin < json.size() &&
         std::isspace(static_cast<unsigned char>(json[digit_begin]))) {
    ++digit_begin;
  }
  std::size_t digit_end = digit_begin;
  while (digit_end < json.size() && json[digit_end] != ',' && json[digit_end] != '\n' &&
         json[digit_end] != '}') {
    ++digit_end;
  }
  assert(digit_end > digit_begin);
  const auto value = std::stod(json.substr(digit_begin, digit_end - digit_begin));
  assert(std::isfinite(value));
  return value;
}

[[nodiscard]] bool json_bool_field(
    const std::string& json,
    const std::string_view name) {
  const std::string marker = "\"" + std::string(name) + "\":";
  const auto value_begin = json.find(marker);
  assert(value_begin != std::string::npos);
  std::size_t bool_begin = value_begin + marker.size();
  while (bool_begin < json.size() &&
         std::isspace(static_cast<unsigned char>(json[bool_begin]))) {
    ++bool_begin;
  }
  if (json.compare(bool_begin, 4, "true") == 0) {
    return true;
  }
  assert(json.compare(bool_begin, 5, "false") == 0);
  return false;
}

struct DispositionExpectation {
  bool diagnostic_only = false;
  bool gate_candidate = false;
  bool integrator_output = false;
  std::string_view candidate_kind;
};

[[nodiscard]] std::string bool_string(const bool value) {
  return value ? "true" : "false";
}

void assert_disposition_report_matches_metadata(
    const std::filesystem::path& output,
    const std::filesystem::path& log_path,
    const DispositionExpectation expected) {
  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open disposition output");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ONLY") ==
         bool_string(expected.diagnostic_only));
  assert(read_text_attr(file_id, "TYWRF_GATE_CANDIDATE") ==
         bool_string(expected.gate_candidate));
  assert(read_text_attr(file_id, "TYWRF_INTEGRATOR_OUTPUT") ==
         bool_string(expected.integrator_output));
  assert(read_text_attr(file_id, "TYWRF_CANDIDATE_KIND") ==
         std::string(expected.candidate_kind));
  check_nc(nc_close(file_id), "close disposition output");

  const auto log = read_file(log_path);
  assert(json_bool_field(log, "diagnostic_only") == expected.diagnostic_only);
  assert(json_bool_field(log, "gate_candidate") == expected.gate_candidate);
  assert(json_bool_field(log, "integrator_output") == expected.integrator_output);
  assert(
      log.find("\"candidate_kind\": \"" + std::string(expected.candidate_kind) + "\"") !=
      std::string::npos);
}

void assert_candidate_disposition_self_test(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto log_path = root / "candidate_disposition_self_test.log";
  const auto status = run_command_status(
      shell_quote(executable) + " --candidate-disposition-self-test >" +
      shell_quote(log_path) + " 2>&1");
  assert(status == 0);
  assert(read_file(log_path).empty());
}

void assert_pressure_refresh_readiness_json(const std::string& log) {
  assert(json_bool_field(log, "pressure_refresh_readiness_ready"));
  assert(json_bool_field(log, "thermodynamic_base_state_consistency_ready"));
  assert(json_bool_field(log, "provider_terrain_uses_moved_candidate_hgt"));
  assert(json_bool_field(log, "provider_base_state_reconstruct_ok"));
  assert(json_bool_field(log, "base_state_sync_contract_ok"));
  assert(json_bool_field(log, "base_state_source_sync_readiness_check"));
  assert(!json_bool_field(log, "base_state_source_sync_applied"));
  assert(json_number_field(log, "source_sync_planned_pb_point_count") > 0.0);
  assert(json_number_field(log, "source_sync_planned_mub_point_count") > 0.0);
  assert(json_number_field(log, "source_sync_planned_phb_point_count") > 0.0);
  assert(json_number_field(log, "source_sync_overlap_write_count") == 0.0);
  assert(json_number_field(log, "source_sync_halo_write_count") == 0.0);
  assert(!json_bool_field(log, "pressure_refresh_readiness_compute_called"));
  assert(json_bool_field(log, "pressure_compute_readiness_check"));
  assert(json_bool_field(log, "pressure_compute_readiness_check_called"));
  assert(json_bool_field(log, "pressure_compute_readiness_check_ok"));
  const auto would_refresh_p_points =
      json_number_field(log, "pressure_refresh_planned_p_point_count");
  assert(would_refresh_p_points > 0.0);
  assert(json_number_field(log, "pressure_readiness_invalid_p_point_count") == 0.0);
  assert(json_number_field(log, "pressure_readiness_skipped_p_point_count") == 0.0);
  assert(json_number_field(log, "pressure_compute_readiness_report_target_column_count") > 0.0);
  assert(
      json_number_field(log, "pressure_compute_readiness_report_refreshed_point_count") ==
      would_refresh_p_points);
  assert(json_number_field(log, "pressure_compute_readiness_report_invalid_point_count") == 0.0);
  assert(json_number_field(log, "pressure_compute_readiness_report_skipped_point_count") == 0.0);
  assert(!json_bool_field(log, "pressure_compute_readiness_report_touched_overlap_cells"));
  assert(!json_bool_field(log, "pressure_compute_readiness_report_touched_halo_cells"));
  assert(!json_bool_field(log, "pressure_refresh_readiness_applied"));
  assert(
      log.find("\"pressure_refresh_readiness_provider_terrain_source\": "
               "\"moved_candidate_HGT\"") != std::string::npos);
  assert(
      log.find("\"pressure_refresh_readiness_provider_terrain_provenance\": "
               "\"override:moved_candidate_HGT\"") != std::string::npos);
}

void assert_normal_pressure_refresh_apply(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output,
    const std::filesystem::path& log_path) {
  std::filesystem::remove(output);
  assert(!std::filesystem::exists(output));
  run_command(
      base_command(executable, d01_start, d02_start, template_path, output) +
      " --pressure-refresh >" + shell_quote(log_path) + " 2>&1");
  assert(std::filesystem::exists(output));
  const auto log = read_file(log_path);
  assert(log.find("\"status\": \"selected_field_candidate_generated\"") !=
         std::string::npos);
  assert(log.find("\"candidate_kind\": \"selected_field_integrator_v0\"") !=
         std::string::npos);
  assert(log.find("\"gate_candidate\": true") != std::string::npos);
  assert(log.find("\"integrator_output\": true") != std::string::npos);
  assert(log.find("\"experimental_pressure_refresh_apply\": true") ==
         std::string::npos);
  assert(log.find("\"pressure_refresh_experimental_apply\"") == std::string::npos);
  assert(log.find("\"pressure_refresh_applied\": true") != std::string::npos);
  assert(log.find("\"pressure_refresh_source_sync_ok\": true") != std::string::npos);
  assert(log.find("\"pressure_refresh_integration_status\": \"applied_to_candidate\"") !=
         std::string::npos);
  assert_pressure_refresh_readiness_json(log);
  assert(log.find("\"pressure_refresh_terrain_override_used\": true") != std::string::npos);
  assert(log.find("\"pressure_refresh_terrain_source\": \"moved_candidate_HGT\"") !=
         std::string::npos);
  assert(log.find("\"pressure_refresh_terrain_provenance\": \"override:moved_candidate_HGT\"") !=
         std::string::npos);
  const auto target_columns = json_number_field(log, "pressure_refresh_target_column_count");
  const auto refreshed_columns =
      json_number_field(log, "pressure_refresh_refreshed_column_count");
  const auto refreshed_points =
      json_number_field(log, "pressure_refresh_refreshed_point_count");
  const auto refreshed_p_points =
      json_number_field(log, "pressure_refresh_refreshed_p_points");
  assert(target_columns > 0.0);
  assert(refreshed_columns == target_columns);
  assert(refreshed_points > 0.0);
  assert(refreshed_p_points == refreshed_points);
  assert(json_number_field(log, "pressure_refresh_skipped_point_count") == 0.0);
  assert(json_number_field(log, "pressure_refresh_invalid_point_count") == 0.0);
  assert(!json_bool_field(log, "pressure_refresh_touched_overlap_cells"));
  assert(!json_bool_field(log, "pressure_refresh_touched_halo_cells"));
  assert(json_bool_field(log, "pressure_refresh_changed_p_matches_refreshed_point_count"));
  assert(json_bool_field(log, "pressure_refresh_invalid_and_skipped_points_zero"));
  assert(json_bool_field(log, "pressure_refresh_overlap_halo_untouched"));
  assert(log.find("\"normal_base_state_producer_source\": "
                  "\"normal_selected_field_base_state_producer\"") != std::string::npos);
  assert(!json_bool_field(log, "normal_base_state_producer_diagnostic_only"));
  assert(json_bool_field(log, "normal_base_state_producer_normal_candidate_producer"));
  assert(json_bool_field(log, "normal_base_state_producer_writes_candidate"));
  assert(!json_bool_field(log, "normal_base_state_producer_uses_reference_end_truth"));
  assert(json_bool_field(log, "normal_base_state_producer_no_reference_end_truth"));
  assert(!json_bool_field(log, "normal_base_state_producer_uses_direct_p_shortcut"));
  assert(json_bool_field(log, "normal_base_state_producer_no_direct_p_shortcut"));
  assert(!json_bool_field(log, "normal_base_state_producer_reads_direct_p"));
  assert(!json_bool_field(log, "normal_base_state_producer_writes_p"));
  assert(json_bool_field(log, "normal_base_state_producer_no_gate_pass_claim"));
  assert(json_bool_field(log, "normal_base_state_producer_no_00_20_progression"));
  assert(json_number_field(log, "normal_base_state_producer_changed_p_points") == 0.0);
  assert(json_number_field(log, "normal_base_state_producer_changed_pb_points") > 0.0);
  assert(json_number_field(log, "normal_base_state_producer_changed_mub_points") > 0.0);
  assert(json_number_field(log, "normal_base_state_producer_changed_phb_points") > 0.0);
  assert_successful_candidate(output, true, template_path);
  assert_disposition_report_matches_metadata(
      output,
      log_path,
      {false, true, true, "selected_field_integrator_v0"});
}

void assert_pressure_column_probe_default_path(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output,
    const std::filesystem::path& log_path) {
  std::filesystem::remove(output);
  assert(!std::filesystem::exists(output));
  run_command(
      base_command(executable, d01_start, d02_start, template_path, output) +
      " --pressure-column-probe 9,9 --pressure-column-levels 0,1 >" +
      shell_quote(log_path) + " 2>&1");
  assert(std::filesystem::exists(output));
  assert_successful_candidate(output, false, std::filesystem::path{}, false, true);
  assert_pressure_column_probe_output(output, log_path, false);
}

void assert_pressure_column_probe_pressure_refresh_path(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output,
    const std::filesystem::path& log_path) {
  std::filesystem::remove(output);
  assert(!std::filesystem::exists(output));
  run_command(
      base_command(executable, d01_start, d02_start, template_path, output) +
      " --pressure-refresh --pressure-column-probe 9,9 --pressure-column-levels 0,1 >" +
      shell_quote(log_path) + " 2>&1");
  assert(std::filesystem::exists(output));
  assert_successful_candidate(output, true, template_path, false, true);
  assert_pressure_column_probe_output(output, log_path, true);
}

void assert_pressure_formula_observation_outside_target_status(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output,
    const std::filesystem::path& log_path) {
  std::filesystem::remove(output);
  assert(!std::filesystem::exists(output));
  run_command(
      base_command(executable, d01_start, d02_start, template_path, output) +
      " --pressure-refresh --pressure-column-probe " +
      shell_quote(std::string("9,9;0,0")) + " --pressure-column-levels 0 >" +
      shell_quote(log_path) + " 2>&1");
  assert(std::filesystem::exists(output));
  assert_successful_candidate(output, true, template_path, false, true);

  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open formula output");
  assert(read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_REQUEST_COUNT") == 2);
  assert(read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_RECORD_COUNT") == 2);
  assert(read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALID_COUNT") == 1);
  assert(read_int_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_INVALID_COUNT") == 0);
  assert(
      read_int_attr(
          file_id,
          "TYWRF_PRESSURE_FORMULA_OBSERVATION_OUTSIDE_TARGET_REGION_COUNT") == 1);
  const auto values =
      read_text_attr(file_id, "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALUES");
  assert(values.find("status=recorded;valid=1;i=9;j=9;k=0") != std::string::npos);
  assert(
      values.find("status=request_outside_target_region;valid=0;i=0;j=0;k=0") !=
      std::string::npos);
  check_nc(nc_close(file_id), "close formula output");

  const auto log = read_file(log_path);
  assert(
      log.find("\"pressure_formula_observation_valid_count\": 1") !=
      std::string::npos);
  assert(
      log.find("\"pressure_formula_observation_outside_target_region_count\": 1") !=
      std::string::npos);
  assert(
      log.find("status=request_outside_target_region;valid=0;i=0;j=0;k=0") !=
      std::string::npos);
}

void assert_pressure_refresh_invalid_dry_run_fail_closed(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output,
    const std::filesystem::path& log_path) {
  std::filesystem::remove(output);
  assert(!std::filesystem::exists(output));
  const auto status = run_command_status(
      base_command(executable, d01_start, d02_start, template_path, output) +
      " --pressure-refresh >" + shell_quote(log_path) + " 2>&1");
  assert(status != 0);
  assert(!std::filesystem::exists(output));
  const auto log = read_file(log_path);
  assert(log.find("pressure_refresh_dry_run_contract_failed") != std::string::npos);
  assert(log.find("provider_ok=true") != std::string::npos);
  assert(log.find("base_state_sync_contract_ok=true") != std::string::npos);
  assert(log.find("base_state_sync_dry_run=true") != std::string::npos);
  assert(log.find("base_state_sync_applied=false") != std::string::npos);
  assert(log.find("pressure_compute_dry_run=true") != std::string::npos);
  assert(log.find("pressure_compute_dry_run_called=true") != std::string::npos);
  assert(log.find("pressure_compute_dry_run_ok=false") != std::string::npos);
  assert(log_u64_field(log, "dry_run_invalid_p_point_count") > 0);
  assert(log_u64_field(log, "dry_run_skipped_p_point_count") > 0);
  assert(log_u64_field(log, "pressure_compute_dry_run_report_invalid_point_count") > 0);
  assert(log_u64_field(log, "pressure_compute_dry_run_report_skipped_point_count") > 0);
  assert(log.find("pressure_compute_dry_run_report_touched_overlap_cells=false") !=
         std::string::npos);
  assert(log.find("pressure_compute_dry_run_report_touched_halo_cells=false") !=
         std::string::npos);
  assert(log.find("pressure_refresh_applied=false") != std::string::npos);
}

void assert_hidden_apply_flag_absent_from_help(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto help_log = root / "help.log";
  const auto status = run_command_status(
      shell_quote(executable) + " --help >" + shell_quote(help_log) + " 2>&1");
  assert(status == 0);
  const auto help = read_file(help_log);
  assert(help.find("--wind-tendency-source") != std::string::npos);
  assert(help.find("--wind-tendency-substeps") != std::string::npos);
  assert(help.find("--wind-tendency-advecting-velocity") != std::string::npos);
  assert(help.find("--wind-tendency-advecting-components") != std::string::npos);
  assert(help.find("--pressure-refresh") != std::string::npos);
  assert(help.find("--pressure-column-probe") != std::string::npos);
  assert(help.find("--pressure-column-levels") != std::string::npos);
  assert(help.find("--experimental-pressure-refresh-apply") == std::string::npos);
  assert(help.find("--diagnostic-base-state-adapter-report") == std::string::npos);
  assert(help.find("--candidate-disposition-self-test") == std::string::npos);
}

void assert_experimental_pressure_refresh_apply(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output,
    const std::filesystem::path& log_path) {
  std::filesystem::remove(output);
  assert(!std::filesystem::exists(output));
  run_command(
      base_command(executable, d01_start, d02_start, template_path, output) +
      " --pressure-refresh --experimental-pressure-refresh-apply >" + shell_quote(log_path) +
      " 2>&1");
  assert(std::filesystem::exists(output));
  const auto log = read_file(log_path);
  assert(
      log.find("\"status\": \"selected_field_pressure_refresh_experimental_apply_generated\"") !=
      std::string::npos);
  assert(log.find("\"candidate_kind\": \"selected_field_pressure_refresh_experimental_apply_v0\"") !=
         std::string::npos);
  assert(log.find("\"gate_candidate\": false") != std::string::npos);
  assert(log.find("\"integrator_output\": false") != std::string::npos);
  assert(log.find("\"experimental_pressure_refresh_apply\": true") != std::string::npos);
  assert(log.find("\"pressure_refresh_applied\": true") != std::string::npos);
  assert(log.find("\"pressure_refresh_experimental_apply\": true") != std::string::npos);
  assert(log.find("\"pressure_refresh_integration_status\": \"experimental_apply_test_only\"") !=
         std::string::npos);
  assert_pressure_refresh_readiness_json(log);
  assert(log.find("\"pressure_refresh_terrain_override_used\": true") != std::string::npos);
  assert(log.find("\"pressure_refresh_terrain_source\": \"moved_candidate_HGT\"") !=
         std::string::npos);
  assert(log.find("\"pressure_refresh_terrain_source\": \"HGT\"") == std::string::npos);
  assert(log.find("\"pressure_refresh_terrain_provenance\": \"override:moved_candidate_HGT\"") !=
         std::string::npos);
  const auto target_columns = json_number_field(log, "pressure_refresh_target_column_count");
  const auto refreshed_columns =
      json_number_field(log, "pressure_refresh_refreshed_column_count");
  const auto refreshed_points =
      json_number_field(log, "pressure_refresh_refreshed_point_count");
  const auto refreshed_p_points =
      json_number_field(log, "pressure_refresh_refreshed_p_points");
  const auto skipped_points = json_number_field(log, "pressure_refresh_skipped_point_count");
  const auto invalid_points = json_number_field(log, "pressure_refresh_invalid_point_count");
  const auto changed_p_points = json_number_field(log, "pressure_refresh_changed_p_points");
  const auto changed_pb_points = json_number_field(log, "pressure_refresh_changed_pb_points");
  const auto changed_mub_points = json_number_field(log, "pressure_refresh_changed_mub_points");
  const auto changed_phb_points = json_number_field(log, "pressure_refresh_changed_phb_points");
  assert(target_columns > 0.0);
  assert(refreshed_columns == target_columns);
  assert(refreshed_points > 0.0);
  assert(refreshed_p_points == refreshed_points);
  assert(skipped_points == 0.0);
  assert(invalid_points == 0.0);
  assert(!json_bool_field(log, "pressure_refresh_touched_overlap_cells"));
  assert(!json_bool_field(log, "pressure_refresh_touched_halo_cells"));
  assert(changed_p_points == refreshed_points);
  assert(changed_pb_points >= 0.0);
  assert(changed_mub_points >= 0.0);
  assert(changed_phb_points >= 0.0);
  assert(json_bool_field(log, "pressure_refresh_changed_p_matches_refreshed_point_count"));
  assert(json_bool_field(log, "pressure_refresh_invalid_and_skipped_points_zero"));
  assert(json_bool_field(log, "pressure_refresh_overlap_halo_untouched"));
  assert_successful_candidate(output, true, template_path, true);
  assert_disposition_report_matches_metadata(
      output,
      log_path,
      {true, false, false, "selected_field_pressure_refresh_experimental_apply_v0"});
}

void assert_diagnostic_adapter_provider_source_log(const std::string& log) {
  assert(
      log.find("\"diagnostic_adapter_provider_source_version\": "
               "\"d77_provider_source_v0\"") != std::string::npos);
  assert(
      log.find("\"diagnostic_adapter_provider_source_origin\": "
               "\"base_state_reconstruction_provider+moved_candidate_HGT\"") !=
      std::string::npos);
  assert(
      log.find("\"diagnostic_adapter_provider_source_source_origin\": "
               "\"base_state_reconstruction_provider+moved_candidate_HGT\"") !=
      std::string::npos);
  assert(
      log.find("\"diagnostic_adapter_provider_source_provider_source\": "
               "\"base_state_reconstruction_provider\"") != std::string::npos);
  assert(
      log.find("\"diagnostic_adapter_provider_source_terrain_source\": "
               "\"moved_candidate_HGT\"") != std::string::npos);
  assert(
      log.find("\"diagnostic_adapter_provider_source_terrain_provenance\": "
               "\"override:moved_candidate_HGT\"") != std::string::npos);
  assert(
      log.find("\"diagnostic_adapter_provider_source_ht_source\": "
               "\"output_static.hgt\"") != std::string::npos);
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_terrain_override_used"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_provider_ok"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_diagnostic_only"));
  assert(!json_bool_field(log, "diagnostic_adapter_provider_source_gate_candidate"));
  assert(!json_bool_field(log, "diagnostic_adapter_provider_source_integrator_output"));
  assert(!json_bool_field(log, "diagnostic_adapter_provider_source_writes_candidate"));
  assert(!json_bool_field(log, "diagnostic_adapter_provider_source_writes_netcdf"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_no_candidate_write"));
  assert(!json_bool_field(log, "diagnostic_adapter_provider_source_uses_reference_end_truth"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_no_reference_end_truth"));
  assert(!json_bool_field(log, "diagnostic_adapter_provider_source_uses_direct_p_shortcut"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_no_direct_p_shortcut"));
  assert(!json_bool_field(log, "diagnostic_adapter_provider_source_reads_direct_p"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_wrote_pb"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_wrote_t_init"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_wrote_mub"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_wrote_alb"));
  assert(json_bool_field(log, "diagnostic_adapter_provider_source_wrote_phb"));
  assert(json_bool_field(
      log,
      "diagnostic_adapter_provider_source_provider_reconstructed_phb_not_wrf_rebalance_validated"));
}

void assert_diagnostic_adapter_source_child_delta_json_field(
    const std::string& log,
    const std::string_view field) {
  const std::string prefix =
      "diagnostic_adapter_source_child_delta_" + std::string(field);
  assert(json_number_field(log, prefix + "_compared_value_count") > 0.0);
  assert(json_number_field(log, prefix + "_differing_value_count") >= 0.0);
  assert(json_number_field(log, prefix + "_max_abs_diff") >= 0.0);
}

void assert_diagnostic_adapter_source_child_delta_log(const std::string& log) {
  assert(
      log.find("\"diagnostic_adapter_source_child_delta_version\": "
               "\"a76_source_child_delta_v0\"") != std::string::npos);
  assert(json_bool_field(log, "diagnostic_adapter_source_child_delta_diagnostic_only"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_child_delta_gate_candidate"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_child_delta_integrator_output"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_child_delta_writes_candidate"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_child_delta_writes_netcdf"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_child_delta_values_identical"));
  const auto compared =
      json_number_field(log, "diagnostic_adapter_source_child_delta_compared_value_count");
  assert(compared > 0.0);
  assert(
      compared ==
      json_number_field(log, "diagnostic_adapter_source_staging_staged_value_count"));
  assert(
      json_number_field(log, "diagnostic_adapter_source_child_delta_differing_value_count") >
      0.0);
  assert(json_number_field(log, "diagnostic_adapter_source_child_delta_max_abs_diff") > 0.0);
  assert_diagnostic_adapter_source_child_delta_json_field(log, "phb");
  assert_diagnostic_adapter_source_child_delta_json_field(log, "mub");
  assert_diagnostic_adapter_source_child_delta_json_field(log, "ht");
  assert_diagnostic_adapter_source_child_delta_json_field(log, "pb");
  assert_diagnostic_adapter_source_child_delta_json_field(log, "t_init");
  assert_diagnostic_adapter_source_child_delta_json_field(log, "alb");
}

void assert_diagnostic_adapter_source_child_delta_attr_field(
    const int file_id,
    const std::string_view field) {
  const std::string prefix =
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_" + std::string(field);
  assert(read_double_attr(file_id, prefix + "_COMPARED_VALUE_COUNT") > 0.0);
  assert(read_double_attr(file_id, prefix + "_DIFFERING_VALUE_COUNT") >= 0.0);
  assert(read_double_attr(file_id, prefix + "_MAX_ABS_DIFF") >= 0.0);
}

void assert_diagnostic_adapter_provider_source_attrs(const int file_id) {
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_VERSION") ==
      "d77_provider_source_v0");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_ORIGIN") ==
      "base_state_reconstruction_provider+moved_candidate_HGT");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_SOURCE_ORIGIN") ==
      "base_state_reconstruction_provider+moved_candidate_HGT");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_PROVIDER_SOURCE") ==
      "base_state_reconstruction_provider");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_TERRAIN_SOURCE") ==
      "moved_candidate_HGT");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_TERRAIN_PROVENANCE") ==
      "override:moved_candidate_HGT");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_TERRAIN_OVERRIDE_USED") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_HT_SOURCE") ==
      "output_static.hgt");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_PROVIDER_OK") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_DIAGNOSTIC_ONLY") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_GATE_CANDIDATE") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_INTEGRATOR_OUTPUT") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WRITES_CANDIDATE") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WRITES_NETCDF") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_NO_CANDIDATE_WRITE") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_USES_REFERENCE_END_TRUTH") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_NO_REFERENCE_END_TRUTH") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_USES_DIRECT_P_SHORTCUT") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_NO_DIRECT_P_SHORTCUT") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_READS_DIRECT_P") ==
      "false");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_PB") == "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_T_INIT") ==
      "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_MUB") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_ALB") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_PHB") == "true");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_PROVIDER_RECONSTRUCTED_PHB_NOT_WRF_REBALANCE_VALIDATED") ==
      "true");
}

void assert_diagnostic_adapter_source_child_delta_attrs(const int file_id) {
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_VERSION") ==
      "a76_source_child_delta_v0");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_SOURCE") ==
      "BaseStateSourceStagingProvider_vs_child_staging_pre_adapter");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_SCOPE") ==
      "exposed_base_state_values_only");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_FIELDS") ==
      "PHB,MUB,HT,PB,T_INIT,ALB");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_DIAGNOSTIC_ONLY") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_GATE_CANDIDATE") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_INTEGRATOR_OUTPUT") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_WRITES_CANDIDATE") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_WRITES_NETCDF") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_VALUES_IDENTICAL") ==
      "false");
  const auto compared = read_double_attr(
      file_id,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_COMPARED_VALUE_COUNT");
  assert(compared > 0.0);
  assert(
      compared ==
      read_double_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_VALUE_COUNT"));
  assert(
      read_double_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_DIFFERING_VALUE_COUNT") > 0.0);
  assert(
      read_double_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_MAX_ABS_DIFF") >
      0.0);
  assert_diagnostic_adapter_source_child_delta_attr_field(file_id, "PHB");
  assert_diagnostic_adapter_source_child_delta_attr_field(file_id, "MUB");
  assert_diagnostic_adapter_source_child_delta_attr_field(file_id, "HT");
  assert_diagnostic_adapter_source_child_delta_attr_field(file_id, "PB");
  assert_diagnostic_adapter_source_child_delta_attr_field(file_id, "T_INIT");
  assert_diagnostic_adapter_source_child_delta_attr_field(file_id, "ALB");
}

void assert_diagnostic_adapter_report_path(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& output,
    const std::filesystem::path& log_path) {
  std::filesystem::remove(output);
  assert(!std::filesystem::exists(output));
  run_command(
      base_command(executable, d01_start, d02_start, template_path, output) +
      " --diagnostic-base-state-adapter-report >" + shell_quote(log_path) + " 2>&1");
  assert(std::filesystem::exists(output));

  const auto log = read_file(log_path);
  assert(log.find("\"status\": \"selected_field_diagnostic_adapter_reported\"") !=
         std::string::npos);
  assert(log.find("\"candidate_kind\": \"selected_field_diagnostic_adapter_v0\"") !=
         std::string::npos);
  assert(log.find("\"diagnostic_only\": true") != std::string::npos);
  assert(log.find("\"gate_candidate\": false") != std::string::npos);
  assert(log.find("\"integrator_output\": false") != std::string::npos);
  assert(log.find("\"pressure_refresh_applied\": true") == std::string::npos);
  assert(log.find("\"pressure_refresh_integration_status\"") == std::string::npos);
  assert(log.find("\"diagnostic_adapter_opt_in\": true") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_ok\": true") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_diagnostic_only\": true") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_gate_candidate\": false") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_integrator_output\": false") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_writes_netcdf\": false") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_writes_candidate\": false") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_called_d68_exchange\": true") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_called_d69_recompute\": true") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_source_staging_version\": "
                  "\"d75_provider_report_v0\"") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_source_staging_provider_kind\": "
                  "\"BaseStateSourceStagingProvider\"") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_source_staging_source\": "
                  "\"explicit_base_state_source_staging_provider\"") != std::string::npos);
  assert(log.find("\"diagnostic_adapter_source_staging_source\": "
                  "\"selected_field_integrator_v0\"") == std::string::npos);
  assert_diagnostic_adapter_provider_source_log(log);
  assert(json_bool_field(log, "diagnostic_adapter_source_staging_ok"));
  assert(json_bool_field(log, "diagnostic_adapter_source_staging_diagnostic_only"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_staging_gate_candidate"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_staging_integrator_output"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_staging_writes_candidate"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_staging_writes_netcdf"));
  assert(json_bool_field(log, "diagnostic_adapter_source_staging_candidate_buffers_preserved"));
  assert(json_bool_field(log, "diagnostic_adapter_source_staging_owns_staging_buffers"));
  assert(json_bool_field(log, "diagnostic_adapter_source_staging_allocated_buffers"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_staging_uses_reference_end_truth"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_staging_uses_direct_p_shortcut"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_staging_reads_direct_p"));
  assert(!json_bool_field(log, "diagnostic_adapter_source_staging_aliases_child"));
  assert(log.find("\"diagnostic_adapter_integration_status\": "
                  "\"staging_report_only_no_gate_no_integrator\"") != std::string::npos);
  assert(json_number_field(log, "diagnostic_adapter_c3h_count") == 2.0);
  assert(json_number_field(log, "diagnostic_adapter_c4h_count") == 2.0);
  assert(json_number_field(log, "diagnostic_adapter_exposed_mass_cell_count") > 0.0);
  assert(json_number_field(log, "diagnostic_adapter_recomputed_point_count") > 0.0);
  assert(json_number_field(log, "diagnostic_adapter_invalid_point_count") == 0.0);
  assert(json_number_field(log, "diagnostic_adapter_source_staging_exposed_region_count") > 0.0);
  assert(json_number_field(log, "diagnostic_adapter_source_staging_exposed_mass_cell_count") > 0.0);
  assert(json_number_field(log, "diagnostic_adapter_source_staging_exposed_mass_point_count") > 0.0);
  assert(json_number_field(log, "diagnostic_adapter_source_staging_masked_mass_cell_count") > 0.0);
  assert(json_number_field(log, "diagnostic_adapter_source_staging_staged_value_count") > 0.0);
  assert(json_number_field(log, "diagnostic_adapter_source_staging_invalid_exposed_value_count") == 0.0);
  assert_diagnostic_adapter_source_child_delta_log(log);
  assert_disposition_report_matches_metadata(
      output,
      log_path,
      {true, false, false, "selected_field_diagnostic_adapter_v0"});

  int file_id = -1;
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open diagnostic adapter output");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_OPT_IN") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_OK") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_DIAGNOSTIC_ONLY") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_GATE_CANDIDATE") == "false");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_INTEGRATOR_OUTPUT") == "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SELECTED_FIELD_NUMERICS_ENABLED") ==
      "false");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_WRITES_NETCDF") == "false");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_WRITES_CANDIDATE") == "false");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_CALLED_D68_EXCHANGE") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_CALLED_D69_RECOMPUTE") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_HT_DIAGNOSTIC_LABEL") == "HGT");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_HT_IS_HGT_ALIAS") == "true");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_CREATES_TERRAIN_OWNER") == "false");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_TERRAIN_OWNER_CREATED") == "false");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_METADATA_SOURCE") ==
         template_path.string());
  assert(read_double_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_C3H_COUNT") == 2.0);
  assert(read_double_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_C4H_COUNT") == 2.0);
  assert(read_double_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_EXPOSED_MASS_CELL_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_RECOMPUTED_POINT_COUNT") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_INVALID_POINT_COUNT") == 0.0);
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_VERSION") ==
      "d75_provider_report_v0");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_PROVIDER_KIND") ==
      "BaseStateSourceStagingProvider");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_SOURCE") ==
      "explicit_base_state_source_staging_provider");
  assert(read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_OK") == "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_DIAGNOSTIC_ONLY") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_GATE_CANDIDATE") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INTEGRATOR_OUTPUT") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_WRITES_CANDIDATE") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_WRITES_NETCDF") ==
      "false");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_CANDIDATE_BUFFERS_PRESERVED") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_OWNS_STAGING_BUFFERS") ==
      "true");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_ALLOCATED_BUFFERS") ==
      "true");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_USES_REFERENCE_END_TRUTH") == "false");
  assert(
      read_text_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_USES_DIRECT_P_SHORTCUT") == "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_READS_DIRECT_P") ==
      "false");
  assert(
      read_text_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_ALIASES_CHILD") ==
      "false");
  assert(
      read_double_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_REGION_COUNT") > 0.0);
  assert(
      read_double_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_MASS_CELL_COUNT") > 0.0);
  assert(
      read_double_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_MASS_POINT_COUNT") > 0.0);
  assert(
      read_double_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_MASKED_MASS_CELL_COUNT") > 0.0);
  assert(
      read_double_attr(file_id, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_VALUE_COUNT") >
      0.0);
  assert(
      read_double_attr(
          file_id,
          "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INVALID_EXPOSED_VALUE_COUNT") == 0.0);
  assert_diagnostic_adapter_provider_source_attrs(file_id);
  assert_diagnostic_adapter_source_child_delta_attrs(file_id);
  assert(!has_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_OPT_IN"));
  assert(!has_text_attr(file_id, "TYWRF_PRESSURE_REFRESH_APPLIED"));
  const auto events = read_text_attr(file_id, "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS");
  assert(events.find(":diagnostic_adapter_provider_source(") != std::string::npos);
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "origin") ==
      "base_state_reconstruction_provider+moved_candidate_HGT");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "source_origin") ==
      "base_state_reconstruction_provider+moved_candidate_HGT");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "provider_ok") ==
      "true");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "gate_candidate") ==
      "false");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "integrator_output") ==
      "false");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "writes_candidate") ==
      "false");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "reads_direct_p") ==
      "false");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "terrain_source") ==
      "moved_candidate_HGT");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "terrain_provenance") ==
      "override:moved_candidate_HGT");
  assert(
      timeline_field_value(events, "diagnostic_adapter_provider_source", "wrote_phb") ==
      "true");
  assert(
      timeline_field_value(
          events,
          "diagnostic_adapter_provider_source",
          "provider_reconstructed_phb_not_wrf_rebalance_validated") == "true");
  assert(events.find(":diagnostic_adapter_source_staging(") != std::string::npos);
  assert(
      timeline_field_value(events, "diagnostic_adapter_source_staging", "provider") ==
      "BaseStateSourceStagingProvider");
  assert(timeline_field_value(events, "diagnostic_adapter_source_staging", "ok") == "true");
  assert(
      timeline_field_value(events, "diagnostic_adapter_source_staging", "gate_candidate") ==
      "false");
  assert(
      timeline_field_value(events, "diagnostic_adapter_source_staging", "writes_candidate") ==
      "false");
  assert(
      timeline_field_value(events, "diagnostic_adapter_source_staging", "aliases_child") ==
      "false");
  assert(events.find(":diagnostic_adapter_source_child_delta(") != std::string::npos);
  assert(
      timeline_field_value(events, "diagnostic_adapter_source_child_delta", "diagnostic_only") ==
      "true");
  assert(
      timeline_field_value(events, "diagnostic_adapter_source_child_delta", "gate_candidate") ==
      "false");
  assert(
      timeline_field_value(
          events,
          "diagnostic_adapter_source_child_delta",
          "integrator_output") == "false");
  assert(
      timeline_field_value(
          events,
          "diagnostic_adapter_source_child_delta",
          "writes_candidate") == "false");
  assert(
      timeline_field_value(events, "diagnostic_adapter_source_child_delta", "writes_netcdf") ==
      "false");
  assert(
      timeline_field_value(
          events,
          "diagnostic_adapter_source_child_delta",
          "values_identical") == "false");
  assert(
      timeline_u64_field(events, "diagnostic_adapter_source_child_delta", "compared_values") ==
      static_cast<std::uint64_t>(
          read_double_attr(
              file_id,
              "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_COMPARED_VALUE_COUNT")));
  assert(
      timeline_u64_field(events, "diagnostic_adapter_source_child_delta", "differing_values") ==
      static_cast<std::uint64_t>(
          read_double_attr(
              file_id,
              "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_DIFFERING_VALUE_COUNT")));
  assert(timeline_u64_field(events, "diagnostic_adapter_source_child_delta", "differing_values") >
         0U);
  assert(
      std::stod(
          timeline_field_value(
              events,
              "diagnostic_adapter_source_child_delta",
              "max_abs_diff")) > 0.0);
  assert(events.find(":diagnostic_adapter_report(") != std::string::npos);
  assert(timeline_field_value(events, "diagnostic_adapter_report", "opt_in") == "true");
  assert(timeline_field_value(events, "diagnostic_adapter_report", "diagnostic_only") == "true");
  assert(timeline_field_value(events, "diagnostic_adapter_report", "gate_candidate") == "false");
  assert(timeline_field_value(events, "diagnostic_adapter_report", "integrator_output") == "false");
  assert(timeline_field_value(events, "diagnostic_adapter_report", "writes_candidate") == "false");
  assert(timeline_field_value(events, "pressure_refresh_apply", "applied") == "false");
  const auto candidate_message = read_text_attr(file_id, "TYWRF_CANDIDATE_MESSAGE");
  assert(candidate_message.find("Diagnostic-only selected-field base-state adapter") !=
         std::string::npos);
  assert(candidate_message.find("not a validation gate pass or normal integrator output") !=
         std::string::npos);
  assert(candidate_message.find("staging buffers only") != std::string::npos);

  const auto preserved_exposed_p = linear_3d(1, 9, 9, 60'000.0F, 1.0F, 10.0F, 100.0F);
  const auto preserved_pb = linear_3d(1, 9, 9, 80'000.0F, 1.0F, 10.0F, 100.0F);
  const auto preserved_phb = linear_3d(2, 9, 9, 90'000.0F, 1.0F, 10.0F, 100.0F);
  const auto preserved_mub = linear_2d(9, 9, 100'000.0F, 1.0F, 10.0F);
  const auto preserved_t =
      parent_linear(3'000.0F, 5.0F, 3.0F, 10.0F, 3.4, 2.4, 1);
  expect_close(read_3d_value(file_id, "P", 1, 9, 9), preserved_exposed_p, "adapter P preserved");
  expect_close(read_3d_value(file_id, "PB", 1, 9, 9), preserved_pb, "adapter PB preserved");
  expect_close(read_3d_value(file_id, "PHB", 2, 9, 9), preserved_phb, "adapter PHB preserved");
  expect_close(read_2d_value(file_id, "MUB", 9, 9), preserved_mub, "adapter MUB preserved");
  expect_close(read_3d_value(file_id, "T", 1, 9, 9), preserved_t, "adapter T preserved");
  check_nc(nc_close(file_id), "close diagnostic adapter output");
}

void assert_nonwind_fields_match(
    const int expected_file_id,
    const int actual_file_id,
    const std::string_view label) {
  assert(read_times(actual_file_id) == read_times(expected_file_id));
  const std::vector<std::string_view> nonwind_variables = {
      "XLAT", "XLONG", "HGT",   "T",    "PH",   "MU", "P",
      "QVAPOR", "PB",   "PHB",  "MUB",  "PSFC", "U10", "V10",
      "T2",   "Q2",    "RAINC", "RAINNC"};
  for (const auto variable : nonwind_variables) {
    expect_float_variable_match(actual_file_id, expected_file_id, variable, label);
  }
}

void assert_wind_tendency_opt_in_output(
    const std::filesystem::path& baseline_output,
    const std::filesystem::path& output,
    const std::filesystem::path& log_path,
    const std::string_view source_kind,
    const bool expect_uv_changes,
    const bool expect_gate_evidence,
    const bool expect_zero_or_identity_only,
    const int expected_self_advection_substeps = 0,
    const std::string_view expected_advecting_velocity_mode = "",
    const std::string_view expected_advecting_components = "",
    const std::string_view expected_advecting_collocation = "") {
  int baseline_file_id = -1;
  int file_id = -1;
  check_nc(nc_open(baseline_output.string().c_str(), NC_NOWRITE, &baseline_file_id), "open wind baseline");
  check_nc(nc_open(output.string().c_str(), NC_NOWRITE, &file_id), "open wind output");
  const bool expect_self_advection = expected_self_advection_substeps > 0;
  const std::string expected_advecting_velocity_mode_text =
      expected_advecting_velocity_mode.empty() && expect_self_advection
          ? "refreshed"
          : std::string(expected_advecting_velocity_mode);
  const std::string expected_advecting_components_text =
      expected_advecting_components.empty() && expect_self_advection
          ? "same_component"
          : std::string(expected_advecting_components);
  const std::string expected_advecting_collocation_text =
      expected_advecting_collocation.empty() && expect_self_advection
          ? (expected_advecting_components_text == "cross_component" ? "average" : "same_grid")
          : std::string(expected_advecting_collocation);

  assert(read_text_attr(file_id, "TYWRF_WIND_TENDENCY_OPT_IN") == "true");
  assert(read_text_attr(file_id, "TYWRF_WIND_TENDENCY_APPLIED") == "true");
  assert(read_text_attr(file_id, "TYWRF_WIND_TENDENCY_SOURCE_KIND") ==
         std::string(source_kind));
  assert(
      read_text_attr(file_id, "TYWRF_WIND_TENDENCY_GATE_EVIDENCE") ==
      (expect_gate_evidence ? "true" : "false"));
  assert(
      read_text_attr(file_id, "TYWRF_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE") ==
      (expect_gate_evidence ? "true" : "false"));
  assert(read_text_attr(file_id, "TYWRF_WIND_TENDENCY_USES_REFERENCE_END_TRUTH") == "false");
  assert(
      read_text_attr(file_id, "TYWRF_WIND_TENDENCY_ZERO_OR_IDENTITY_ONLY") ==
      (expect_zero_or_identity_only ? "true" : "false"));
  assert(read_text_attr(file_id, "TYWRF_WIND_TENDENCY_WRITTEN_FIELDS") == "U,V");
  assert(read_text_attr(file_id, "TYWRF_WIND_TENDENCY_STATUS") == "ok");
  assert(read_double_attr(file_id, "TYWRF_WIND_TENDENCY_ACTIVE_U_POINTS") > 0.0);
  assert(read_double_attr(file_id, "TYWRF_WIND_TENDENCY_ACTIVE_V_POINTS") > 0.0);
  if (expect_self_advection) {
    assert(
        read_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_VELOCITY_MODE") ==
        expected_advecting_velocity_mode_text);
    assert(
        read_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS") ==
        expected_advecting_components_text);
    assert(
        read_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION") ==
        expected_advecting_collocation_text);
    assert(
        read_int_attr(file_id, "TYWRF_WIND_TENDENCY_SUBSTEP_COUNT") ==
        expected_self_advection_substeps);
    assert(read_double_attr(file_id, "TYWRF_WIND_TENDENCY_SUBSTEP_DT_SECONDS") == 8.0);
    assert(
        read_double_attr(file_id, "TYWRF_WIND_TENDENCY_TOTAL_SECONDS") ==
        8.0 * static_cast<double>(expected_self_advection_substeps));
    assert(
        read_double_attr(file_id, "TYWRF_WIND_TENDENCY_UPDATED_U_POINTS") ==
        read_double_attr(file_id, "TYWRF_WIND_TENDENCY_ACTIVE_U_POINTS") *
            static_cast<double>(expected_self_advection_substeps));
    assert(
        read_double_attr(file_id, "TYWRF_WIND_TENDENCY_UPDATED_V_POINTS") ==
        read_double_attr(file_id, "TYWRF_WIND_TENDENCY_ACTIVE_V_POINTS") *
            static_cast<double>(expected_self_advection_substeps));
  } else {
    assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_VELOCITY_MODE"));
    assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_COMPONENTS"));
    assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_ADVECTING_COLLOCATION"));
    assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_SUBSTEP_COUNT"));
    assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_SUBSTEP_DT_SECONDS"));
    assert(!has_text_attr(file_id, "TYWRF_WIND_TENDENCY_TOTAL_SECONDS"));
  }

  const auto changed_u = read_double_attr(file_id, "TYWRF_WIND_TENDENCY_CHANGED_U_POINTS");
  const auto changed_v = read_double_attr(file_id, "TYWRF_WIND_TENDENCY_CHANGED_V_POINTS");
  if (expect_uv_changes) {
    assert(changed_u > 0.0);
    assert(changed_v > 0.0);
    assert(
        std::fabs(
            read_3d_value(file_id, "U", 1, 2, 0) -
            read_3d_value(baseline_file_id, "U", 1, 2, 0)) > kTolerance);
    assert(
        std::fabs(
            read_3d_value(file_id, "V", 1, 10, 9) -
            read_3d_value(baseline_file_id, "V", 1, 10, 9)) > kTolerance);
  } else {
    assert(changed_u == 0.0);
    assert(changed_v == 0.0);
    expect_close(
        read_3d_value(file_id, "U", 1, 2, 0),
        read_3d_value(baseline_file_id, "U", 1, 2, 0),
        "zero wind U unchanged");
    expect_close(
        read_3d_value(file_id, "V", 1, 10, 9),
        read_3d_value(baseline_file_id, "V", 1, 10, 9),
        "zero wind V unchanged");
  }
  assert_nonwind_fields_match(baseline_file_id, file_id, source_kind);

  const auto events = read_text_attr(file_id, "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS");
  assert_contains_in_order(
      events,
      {
          ":parent_interpolation(",
          ":wind_tendency_apply(",
          ":selected_field_change_summary(",
      });
  assert(timeline_field_value(events, "wind_tendency_apply", "opt_in") == "true");
  assert(timeline_field_value(events, "wind_tendency_apply", "applied") == "true");
  assert(timeline_field_value(events, "wind_tendency_apply", "source_kind") ==
         std::string(source_kind));
  assert(timeline_field_value(events, "wind_tendency_apply", "fields") == "U_V");
  assert(
      timeline_field_value(events, "wind_tendency_apply", "zero_or_identity_only") ==
      (expect_zero_or_identity_only ? "true" : "false"));
  assert(
      timeline_field_value(events, "wind_tendency_apply", "gate_evidence") ==
      (expect_gate_evidence ? "true" : "false"));
  assert(
      timeline_field_value(events, "wind_tendency_apply", "validation_gate_evidence") ==
      (expect_gate_evidence ? "true" : "false"));
  assert(
      timeline_field_value(events, "wind_tendency_apply", "uses_reference_end_truth") ==
      "false");
  if (expect_self_advection) {
    assert(
        timeline_field_value(
            events, "wind_tendency_apply", "advecting_velocity_mode") ==
        expected_advecting_velocity_mode_text);
    assert(
        timeline_field_value(events, "wind_tendency_apply", "advecting_components") ==
        expected_advecting_components_text);
    assert(
        timeline_field_value(events, "wind_tendency_apply", "advecting_collocation") ==
        expected_advecting_collocation_text);
    assert(
        timeline_field_value(events, "wind_tendency_apply", "substep_count") ==
        std::to_string(expected_self_advection_substeps));
    assert(timeline_field_value(events, "wind_tendency_apply", "substep_dt_seconds") == "8");
    assert(
        timeline_field_value(events, "wind_tendency_apply", "total_seconds") ==
        std::to_string(8 * expected_self_advection_substeps));
  }
  assert(
      timeline_u64_field(events, "wind_tendency_apply", "changed_u_points") ==
      static_cast<std::uint64_t>(changed_u));
  assert(
      timeline_u64_field(events, "wind_tendency_apply", "changed_v_points") ==
      static_cast<std::uint64_t>(changed_v));

  check_nc(nc_close(file_id), "close wind output");
  check_nc(nc_close(baseline_file_id), "close wind baseline");

  const auto log = read_file(log_path);
  assert(log.find("\"wind_tendency_opt_in\": true") != std::string::npos);
  assert(log.find("\"wind_tendency_applied\": true") != std::string::npos);
  assert(
      log.find("\"wind_tendency_source_kind\": \"" + std::string(source_kind) + "\"") !=
      std::string::npos);
  assert(
      log.find(
          "\"wind_tendency_gate_evidence\": " +
          std::string(expect_gate_evidence ? "true" : "false")) != std::string::npos);
  assert(
      log.find(
          "\"wind_tendency_validation_gate_evidence\": " +
          std::string(expect_gate_evidence ? "true" : "false")) != std::string::npos);
  assert(
      log.find("\"wind_tendency_uses_reference_end_truth\": false") !=
      std::string::npos);
  assert(
      log.find(
          "\"wind_tendency_zero_or_identity_only\": " +
          std::string(expect_zero_or_identity_only ? "true" : "false")) !=
      std::string::npos);
  if (expect_self_advection) {
    assert(
        log.find(
            "\"wind_tendency_advecting_velocity_mode\": \"" +
            expected_advecting_velocity_mode_text + "\"") != std::string::npos);
    assert(
        log.find(
            "\"wind_tendency_advecting_components\": \"" +
            expected_advecting_components_text + "\"") != std::string::npos);
    assert(
        log.find(
            "\"wind_tendency_advecting_collocation\": \"" +
            expected_advecting_collocation_text + "\"") != std::string::npos);
    assert(
        log.find(
            "\"wind_tendency_substep_count\": " +
            std::to_string(expected_self_advection_substeps)) != std::string::npos);
    assert(log.find("\"wind_tendency_substep_dt_seconds\": 8") != std::string::npos);
    assert(
        log.find(
            "\"wind_tendency_total_seconds\": " +
            std::to_string(8 * expected_self_advection_substeps)) != std::string::npos);
  } else {
    assert(log.find("wind_tendency_substep_count") == std::string::npos);
    assert(log.find("wind_tendency_advecting_velocity_mode") == std::string::npos);
    assert(log.find("wind_tendency_advecting_components") == std::string::npos);
    assert(log.find("wind_tendency_advecting_collocation") == std::string::npos);
  }
}

void assert_wind_tendency_paths(
    const std::filesystem::path& executable,
    const std::filesystem::path& d01_start,
    const std::filesystem::path& d02_start,
    const std::filesystem::path& template_path,
    const std::filesystem::path& baseline_output,
    const std::filesystem::path& root) {
  const auto explicit_none_output = root / "tywrf_selected_field_wind_none";
  const auto zero_output = root / "tywrf_selected_field_wind_zero";
  const auto identity_output = root / "tywrf_selected_field_wind_identity";
  const auto self_advection_output = root / "tywrf_selected_field_wind_self_advection";
  const auto self_advection_explicit_same_output =
      root / "tywrf_selected_field_wind_self_advection_explicit_same";
  const auto self_advection_alias_output =
      root / "tywrf_selected_field_wind_self_advection_alias";
  const auto self_advection_substeps_output =
      root / "tywrf_selected_field_wind_self_advection_substeps";
  const auto self_advection_refreshed_substeps_output =
      root / "tywrf_selected_field_wind_self_advection_refreshed_substeps";
  const auto self_advection_cross_refreshed_substeps_output =
      root / "tywrf_selected_field_wind_self_advection_cross_refreshed_substeps";
  const auto self_advection_cross_frozen_substeps_output =
      root / "tywrf_selected_field_wind_self_advection_cross_frozen_substeps";
  const auto self_advection_frozen_substeps_output =
      root / "tywrf_selected_field_wind_self_advection_frozen_substeps";
  const auto explicit_none_log = root / "wind_none.log";
  const auto zero_log = root / "wind_zero.log";
  const auto identity_log = root / "wind_identity.log";
  const auto self_advection_log = root / "wind_self_advection.log";
  const auto self_advection_explicit_same_log = root / "wind_self_advection_explicit_same.log";
  const auto self_advection_alias_log = root / "wind_self_advection_alias.log";
  const auto self_advection_substeps_log = root / "wind_self_advection_substeps.log";
  const auto self_advection_refreshed_substeps_log =
      root / "wind_self_advection_refreshed_substeps.log";
  const auto self_advection_cross_refreshed_substeps_log =
      root / "wind_self_advection_cross_refreshed_substeps.log";
  const auto self_advection_cross_frozen_substeps_log =
      root / "wind_self_advection_cross_frozen_substeps.log";
  const auto self_advection_frozen_substeps_log =
      root / "wind_self_advection_frozen_substeps.log";

  run_command(
      base_command(executable, d01_start, d02_start, template_path, explicit_none_output) +
      " --wind-tendency-source none >" + shell_quote(explicit_none_log) + " 2>&1");
  int baseline_file_id = -1;
  int none_file_id = -1;
  check_nc(nc_open(baseline_output.string().c_str(), NC_NOWRITE, &baseline_file_id), "open baseline");
  check_nc(nc_open(explicit_none_output.string().c_str(), NC_NOWRITE, &none_file_id), "open explicit none");
  assert_no_wind_tendency_attrs(none_file_id);
  assert_nonwind_fields_match(baseline_file_id, none_file_id, "explicit none");
  expect_close(
      read_3d_value(none_file_id, "U", 1, 2, 0),
      read_3d_value(baseline_file_id, "U", 1, 2, 0),
      "explicit none U");
  expect_close(
      read_3d_value(none_file_id, "V", 1, 10, 9),
      read_3d_value(baseline_file_id, "V", 1, 10, 9),
      "explicit none V");
  check_nc(nc_close(none_file_id), "close explicit none");
  check_nc(nc_close(baseline_file_id), "close baseline");
  assert(read_file(explicit_none_log).find("wind_tendency_opt_in") == std::string::npos);

  run_command(
      base_command(executable, d01_start, d02_start, template_path, zero_output) +
      " --wind-tendency-source zero >" + shell_quote(zero_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output, zero_output, zero_log, "zero", false, false, true);

  run_command(
      base_command(executable, d01_start, d02_start, template_path, identity_output) +
      " --wind-tendency-source identity >" + shell_quote(identity_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output, identity_output, identity_log, "identity", true, false, true);

  run_command(
      base_command(executable, d01_start, d02_start, template_path, self_advection_output) +
      " --wind-tendency-source self-advection >" + shell_quote(self_advection_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output,
      self_advection_output,
      self_advection_log,
      "self_advection",
      true,
      true,
      false,
      1);

  run_command(
      base_command(executable, d01_start, d02_start, template_path, self_advection_explicit_same_output) +
      " --wind-tendency-source self_advection"
      " --wind-tendency-advecting-components same_component >" +
      shell_quote(self_advection_explicit_same_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output,
      self_advection_explicit_same_output,
      self_advection_explicit_same_log,
      "self_advection",
      true,
      true,
      false,
      1,
      "refreshed",
      "same_component",
      "same_grid");
  int default_same_file_id = -1;
  int explicit_same_file_id = -1;
  check_nc(
      nc_open(self_advection_output.string().c_str(), NC_NOWRITE, &default_same_file_id),
      "open default same-component self advection");
  check_nc(
      nc_open(
          self_advection_explicit_same_output.string().c_str(),
          NC_NOWRITE,
          &explicit_same_file_id),
      "open explicit same-component self advection");
  expect_float_variable_match(
      explicit_same_file_id, default_same_file_id, "U", "explicit same-component default U");
  expect_float_variable_match(
      explicit_same_file_id, default_same_file_id, "V", "explicit same-component default V");
  assert_nonwind_fields_match(
      default_same_file_id, explicit_same_file_id, "explicit same-component default");
  check_nc(nc_close(explicit_same_file_id), "close explicit same-component self advection");
  check_nc(nc_close(default_same_file_id), "close default same-component self advection");

  run_command(
      base_command(executable, d01_start, d02_start, template_path, self_advection_alias_output) +
      " --wind-tendency-source self_advection >" + shell_quote(self_advection_alias_log) +
      " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output,
      self_advection_alias_output,
      self_advection_alias_log,
      "self_advection",
      true,
      true,
      false,
      1);

  run_command(
      base_command(executable, d01_start, d02_start, template_path, self_advection_substeps_output) +
      " --wind-tendency-source self_advection --wind-tendency-substeps 3 >" +
      shell_quote(self_advection_substeps_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output,
      self_advection_substeps_output,
      self_advection_substeps_log,
      "self_advection",
      true,
      true,
      false,
      3);
  run_command(
      base_command(
          executable, d01_start, d02_start, template_path, self_advection_refreshed_substeps_output) +
      " --wind-tendency-source self_advection --wind-tendency-substeps 3"
      " --wind-tendency-advecting-velocity refreshed >" +
      shell_quote(self_advection_refreshed_substeps_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output,
      self_advection_refreshed_substeps_output,
      self_advection_refreshed_substeps_log,
      "self_advection",
      true,
      true,
      false,
      3,
      "refreshed");
  run_command(
      base_command(
          executable,
          d01_start,
          d02_start,
          template_path,
          self_advection_cross_refreshed_substeps_output) +
      " --wind-tendency-source self_advection --wind-tendency-substeps 3"
      " --wind-tendency-advecting-velocity refreshed"
      " --wind-tendency-advecting-components cross_component >" +
      shell_quote(self_advection_cross_refreshed_substeps_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output,
      self_advection_cross_refreshed_substeps_output,
      self_advection_cross_refreshed_substeps_log,
      "self_advection",
      true,
      true,
      false,
      3,
      "refreshed",
      "cross_component",
      "average");
  run_command(
      base_command(
          executable, d01_start, d02_start, template_path, self_advection_frozen_substeps_output) +
      " --wind-tendency-source self_advection --wind-tendency-substeps 3"
      " --wind-tendency-advecting-velocity frozen >" +
      shell_quote(self_advection_frozen_substeps_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output,
      self_advection_frozen_substeps_output,
      self_advection_frozen_substeps_log,
      "self_advection",
      true,
      true,
      false,
      3,
      "frozen");
  run_command(
      base_command(
          executable,
          d01_start,
          d02_start,
          template_path,
          self_advection_cross_frozen_substeps_output) +
      " --wind-tendency-source self_advection --wind-tendency-substeps 3"
      " --wind-tendency-advecting-velocity frozen"
      " --wind-tendency-advecting-components cross_component >" +
      shell_quote(self_advection_cross_frozen_substeps_log) + " 2>&1");
  assert_wind_tendency_opt_in_output(
      baseline_output,
      self_advection_cross_frozen_substeps_output,
      self_advection_cross_frozen_substeps_log,
      "self_advection",
      true,
      true,
      false,
      3,
      "frozen",
      "cross_component",
      "average");
  int one_step_file_id = -1;
  int substeps_file_id = -1;
  int refreshed_substeps_file_id = -1;
  int cross_refreshed_substeps_file_id = -1;
  int cross_frozen_substeps_file_id = -1;
  int frozen_substeps_file_id = -1;
  check_nc(nc_open(self_advection_output.string().c_str(), NC_NOWRITE, &one_step_file_id), "open one-step self advection");
  check_nc(nc_open(self_advection_substeps_output.string().c_str(), NC_NOWRITE, &substeps_file_id), "open substep self advection");
  check_nc(
      nc_open(
          self_advection_refreshed_substeps_output.string().c_str(),
          NC_NOWRITE,
          &refreshed_substeps_file_id),
      "open refreshed substep self advection");
  check_nc(
      nc_open(
          self_advection_frozen_substeps_output.string().c_str(),
          NC_NOWRITE,
          &frozen_substeps_file_id),
      "open frozen substep self advection");
  check_nc(
      nc_open(
          self_advection_cross_refreshed_substeps_output.string().c_str(),
          NC_NOWRITE,
          &cross_refreshed_substeps_file_id),
      "open cross-component refreshed substep self advection");
  check_nc(
      nc_open(
          self_advection_cross_frozen_substeps_output.string().c_str(),
          NC_NOWRITE,
          &cross_frozen_substeps_file_id),
      "open cross-component frozen substep self advection");
  assert(
      std::fabs(
          read_3d_value(substeps_file_id, "U", 1, 2, 0) -
          read_3d_value(one_step_file_id, "U", 1, 2, 0)) > kTolerance);
  assert(
      std::fabs(
          read_3d_value(substeps_file_id, "V", 1, 10, 9) -
          read_3d_value(one_step_file_id, "V", 1, 10, 9)) > kTolerance);
  assert_nonwind_fields_match(one_step_file_id, substeps_file_id, "self-advection substeps");
  expect_float_variable_match(
      refreshed_substeps_file_id, substeps_file_id, "U", "explicit refreshed substeps");
  expect_float_variable_match(
      refreshed_substeps_file_id, substeps_file_id, "V", "explicit refreshed substeps");
  assert_nonwind_fields_match(
      substeps_file_id, refreshed_substeps_file_id, "explicit refreshed substeps");
  assert(max_abs_float_variable_diff(frozen_substeps_file_id, substeps_file_id, "U") > kTolerance);
  assert(max_abs_float_variable_diff(frozen_substeps_file_id, substeps_file_id, "V") > kTolerance);
  assert_nonwind_fields_match(
      substeps_file_id, frozen_substeps_file_id, "frozen versus refreshed substeps");
  assert(
      max_abs_float_variable_diff(cross_refreshed_substeps_file_id, refreshed_substeps_file_id, "U") >
      kTolerance);
  assert(
      max_abs_float_variable_diff(cross_refreshed_substeps_file_id, refreshed_substeps_file_id, "V") >
      kTolerance);
  assert_nonwind_fields_match(
      refreshed_substeps_file_id,
      cross_refreshed_substeps_file_id,
      "cross-component refreshed versus same-component refreshed");
  assert_nonwind_fields_match(
      frozen_substeps_file_id,
      cross_frozen_substeps_file_id,
      "cross-component frozen versus same-component frozen");
  check_nc(
      nc_close(cross_frozen_substeps_file_id),
      "close cross-component frozen substep self advection");
  check_nc(
      nc_close(cross_refreshed_substeps_file_id),
      "close cross-component refreshed substep self advection");
  check_nc(nc_close(frozen_substeps_file_id), "close frozen substep self advection");
  check_nc(nc_close(refreshed_substeps_file_id), "close refreshed substep self advection");
  check_nc(nc_close(substeps_file_id), "close substep self advection");
  check_nc(nc_close(one_step_file_id), "close one-step self advection");
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

  const std::vector<std::string> oracle_options = {
      "--end-state=" + d02_start.string(),
      "--reference-end " + shell_quote(d02_start),
      "--reference-end=" + d02_start.string(),
      "--d01-end-state " + shell_quote(d01_start),
      "--d01-end-state=" + d01_start.string(),
      "--d02-end-state " + shell_quote(d02_start),
      "--d02-end-state=" + d02_start.string(),
  };
  for (std::size_t index = 0; index < oracle_options.size(); ++index) {
    const auto output = root / ("oracle_rejected_" + std::to_string(index));
    const auto status = run_command_status(
        base_command(executable, d01_start, d02_start, template_path, output) + " " +
        oracle_options[index] + " >/dev/null 2>&1");
    assert(status != 0);
    assert(!std::filesystem::exists(output));
  }

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

  const auto no_change_output = root / "no_change_output";
  const auto no_change_status = run_command_status(
      base_command(executable, d01_start, d02_start, template_path, no_change_output) +
      " --to-parent-start 2,2 >/dev/null 2>&1");
  assert(no_change_status != 0);
  assert(!std::filesystem::exists(no_change_output));

  const auto hidden_without_pressure_output = root / "hidden_without_pressure_output";
  const auto hidden_without_pressure_log = root / "hidden_without_pressure.log";
  const auto hidden_without_pressure_status = run_command_status(
      base_command(
          executable,
          d01_start,
          d02_start,
          template_path,
          hidden_without_pressure_output) +
      " --experimental-pressure-refresh-apply >" + shell_quote(hidden_without_pressure_log) +
      " 2>&1");
  assert(hidden_without_pressure_status != 0);
  assert(!std::filesystem::exists(hidden_without_pressure_output));
  const auto hidden_without_pressure = read_file(hidden_without_pressure_log);
  assert(hidden_without_pressure.find(
             "--experimental-pressure-refresh-apply requires --pressure-refresh") !=
         std::string::npos);

  const auto adapter_pressure_output = root / "adapter_pressure_rejected";
  const auto adapter_pressure_status = run_command_status(
      base_command(executable, d01_start, d02_start, template_path, adapter_pressure_output) +
      " --pressure-refresh --diagnostic-base-state-adapter-report >/dev/null 2>&1");
  assert(adapter_pressure_status != 0);
  assert(!std::filesystem::exists(adapter_pressure_output));

  const auto bad_wind_source_output = root / "bad_wind_source_output";
  const auto bad_wind_source_status = run_command_status(
      base_command(executable, d01_start, d02_start, template_path, bad_wind_source_output) +
      " --wind-tendency-source bogus >/dev/null 2>&1");
  assert(bad_wind_source_status != 0);
  assert(!std::filesystem::exists(bad_wind_source_output));

  const std::vector<std::string> bad_wind_substep_options = {
      "--wind-tendency-source self_advection --wind-tendency-substeps 0",
      "--wind-tendency-source self_advection --wind-tendency-substeps -1",
      "--wind-tendency-source self_advection --wind-tendency-substeps 601",
      "--wind-tendency-source none --wind-tendency-substeps 2",
      "--wind-tendency-source zero --wind-tendency-substeps 2",
      "--wind-tendency-source identity --wind-tendency-substeps 2",
  };
  for (std::size_t index = 0; index < bad_wind_substep_options.size(); ++index) {
    const auto output = root / ("bad_wind_substeps_" + std::to_string(index));
    const auto status = run_command_status(
        base_command(executable, d01_start, d02_start, template_path, output) + " " +
        bad_wind_substep_options[index] + " >/dev/null 2>&1");
    assert(status != 0);
    assert(!std::filesystem::exists(output));
  }

  const std::vector<std::string> bad_wind_advecting_velocity_options = {
      "--wind-tendency-advecting-velocity refreshed",
      "--wind-tendency-source none --wind-tendency-advecting-velocity refreshed",
      "--wind-tendency-source zero --wind-tendency-advecting-velocity refreshed",
      "--wind-tendency-source identity --wind-tendency-advecting-velocity frozen",
      "--wind-tendency-source self_advection --wind-tendency-advecting-velocity stale",
      "--wind-tendency-source self_advection --wind-tendency-advecting-velocity=",
  };
  for (std::size_t index = 0; index < bad_wind_advecting_velocity_options.size(); ++index) {
    const auto output = root / ("bad_wind_advecting_velocity_" + std::to_string(index));
    const auto status = run_command_status(
        base_command(executable, d01_start, d02_start, template_path, output) + " " +
        bad_wind_advecting_velocity_options[index] + " >/dev/null 2>&1");
    assert(status != 0);
    assert(!std::filesystem::exists(output));
  }

  const std::vector<std::string> bad_wind_advecting_components_options = {
      "--wind-tendency-advecting-components same_component",
      "--wind-tendency-source none --wind-tendency-advecting-components same_component",
      "--wind-tendency-source zero --wind-tendency-advecting-components same_component",
      "--wind-tendency-source identity --wind-tendency-advecting-components cross_component",
      "--wind-tendency-source self_advection --wind-tendency-advecting-components stale",
      "--wind-tendency-source self_advection --wind-tendency-advecting-components=",
  };
  for (std::size_t index = 0; index < bad_wind_advecting_components_options.size(); ++index) {
    const auto output = root / ("bad_wind_advecting_components_" + std::to_string(index));
    const auto status = run_command_status(
        base_command(executable, d01_start, d02_start, template_path, output) + " " +
        bad_wind_advecting_components_options[index] + " >/dev/null 2>&1");
    assert(status != 0);
    assert(!std::filesystem::exists(output));
  }

  const std::vector<std::string> bad_probe_options = {
      "--pressure-column-probe 9",
      "--pressure-column-probe bad",
      "--pressure-column-probe -1,0",
      "--pressure-column-probe 99,99",
      "--pressure-column-probe 9,9 --pressure-column-levels 99",
      "--pressure-column-levels 0,1",
  };
  for (std::size_t index = 0; index < bad_probe_options.size(); ++index) {
    const auto output = root / ("bad_pressure_column_probe_" + std::to_string(index));
    const auto status = run_command_status(
        base_command(executable, d01_start, d02_start, template_path, output) + " " +
        bad_probe_options[index] + " >/dev/null 2>&1");
    assert(status != 0);
    assert(!std::filesystem::exists(output));
  }
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
    const auto d01_invalid_pressure_start =
        root / "wrfout_d01_invalid_pressure_2025-07-26_00:00:00";
    const auto d02_start = root / "wrfout_d02_2025-07-26_00:00:00";
    const auto template_path = root / "wrfout_d02_template";
    const auto pressure_template_path = root / "wrfout_d02_pressure_template";
    const auto output = root / "tywrf_selected_field_d02_2025-07-26_00:10:00";
    const auto pressure_output =
        root / "tywrf_selected_field_pressure_d02_2025-07-26_00:10:00";
    const auto pressure_column_probe_output =
        root / "tywrf_selected_field_pressure_column_probe_d02_2025-07-26_00:10:00";
    const auto pressure_column_probe_refresh_output =
        root / "tywrf_selected_field_pressure_column_probe_refresh_d02_2025-07-26_00:10:00";
    const auto pressure_formula_outside_output =
        root / "tywrf_selected_field_pressure_formula_outside_d02_2025-07-26_00:10:00";
    const auto experimental_pressure_output =
        root / "tywrf_selected_field_experimental_pressure_d02_2025-07-26_00:10:00";
    const auto diagnostic_adapter_output =
        root / "tywrf_selected_field_diagnostic_adapter_d02_2025-07-26_00:10:00";
    const auto invalid_pressure_output =
        root / "tywrf_selected_field_invalid_pressure_d02_2025-07-26_00:10:00";
    const auto default_log = root / "selected_field_default.log";
    const auto pressure_log = root / "pressure_refresh_normal_apply.log";
    const auto pressure_column_probe_log = root / "pressure_column_probe.log";
    const auto pressure_column_probe_refresh_log = root / "pressure_column_probe_refresh.log";
    const auto pressure_formula_outside_log = root / "pressure_formula_outside.log";
    const auto experimental_pressure_log = root / "pressure_refresh_experimental_apply.log";
    const auto diagnostic_adapter_log = root / "diagnostic_adapter_report.log";
    const auto invalid_pressure_log = root / "pressure_refresh_invalid_dry_run.log";
    create_wrf_fixture(d01_start, 10000.0, FixtureShape{8, 8}, true, true);
    create_wrf_fixture(
        d01_invalid_pressure_start,
        10000.0,
        FixtureShape{8, 8},
        true,
        true,
        false,
        false,
        true);
    create_wrf_fixture(
        d02_start,
        2000.0,
        FixtureShape{10, 10},
        true,
        false,
        true,
        true,
        false,
        false);
    create_wrf_fixture(
        template_path,
        2000.0,
        FixtureShape{10, 10},
        false,
        false,
        false,
        true);
    create_wrf_fixture(
        pressure_template_path,
        2000.0,
        FixtureShape{10, 10},
        false,
        false,
        false,
        true);

    assert_hidden_apply_flag_absent_from_help(executable, root);
    assert_candidate_disposition_self_test(executable, root);

    run_command(
        base_command(executable, d01_start, d02_start, template_path, output) + " >" +
        shell_quote(default_log) + " 2>&1");
    assert_successful_candidate(output);
    assert_disposition_report_matches_metadata(
        output,
        default_log,
        {false, true, true, "selected_field_integrator_v0"});
    assert_wind_tendency_paths(
        executable,
        d01_start,
        d02_start,
        template_path,
        output,
        root);

    assert_normal_pressure_refresh_apply(
        executable,
        d01_start,
        d02_start,
        pressure_template_path,
        pressure_output,
        pressure_log);
    assert_pressure_column_probe_default_path(
        executable,
        d01_start,
        d02_start,
        template_path,
        pressure_column_probe_output,
        pressure_column_probe_log);
    assert_pressure_column_probe_pressure_refresh_path(
        executable,
        d01_start,
        d02_start,
        pressure_template_path,
        pressure_column_probe_refresh_output,
        pressure_column_probe_refresh_log);
    assert_pressure_formula_observation_outside_target_status(
        executable,
        d01_start,
        d02_start,
        pressure_template_path,
        pressure_formula_outside_output,
        pressure_formula_outside_log);
    assert_experimental_pressure_refresh_apply(
        executable,
        d01_start,
        d02_start,
        pressure_template_path,
        experimental_pressure_output,
        experimental_pressure_log);
    assert_diagnostic_adapter_report_path(
        executable,
        d01_start,
        d02_start,
        pressure_template_path,
        diagnostic_adapter_output,
        diagnostic_adapter_log);
    assert_pressure_refresh_invalid_dry_run_fail_closed(
        executable,
        d01_invalid_pressure_start,
        d02_start,
        pressure_template_path,
        invalid_pressure_output,
        invalid_pressure_log);
    run_rejection_tests(executable, d01_start, d02_start, template_path, root);

    std::filesystem::remove_all(root);
  } catch (const std::exception& error) {
    std::cerr << "selected field cycle tool test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
