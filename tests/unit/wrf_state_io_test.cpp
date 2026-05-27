#include "tywrf/io/wrf_state_io.hpp"

#include <netcdf.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void check_nc(const int status, const std::string_view operation) {
  if (status == NC_NOERR) {
    return;
  }
  throw std::runtime_error(std::string(operation) + ": " + nc_strerror(status));
}

std::filesystem::path reference_dir() {
  if (const char* env = std::getenv("TYWRF_REFERENCE_DIR")) {
    return env;
  }
  return "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF";
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

std::string synthetic_time_string(const int time_index) {
  return time_index == 0 ? "2025-07-26_00:00:00" : "2025-07-26_06:00:00";
}

void put_times(
    const int file_id,
    const int var_id,
    const int time_index,
    const std::string_view value) {
  constexpr std::size_t date_str_len = 19;
  assert(value.size() <= date_str_len);
  std::vector<char> buffer(date_str_len, ' ');
  std::copy(value.begin(), value.end(), buffer.begin());
  const std::size_t start[2] = {static_cast<std::size_t>(time_index), 0};
  const std::size_t count[2] = {1, date_str_len};
  check_nc(nc_put_vara_text(file_id, var_id, start, count, buffer.data()), "write Times");
}

void put_3d(
    const int file_id,
    const int var_id,
    const int time_index,
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
               static_cast<std::size_t>(i)] = value_3d(time_index, k, j, i, offset);
      }
    }
  }
  const std::size_t start[4] = {static_cast<std::size_t>(time_index), 0, 0, 0};
  const std::size_t count[4] = {
      1,
      static_cast<std::size_t>(nz),
      static_cast<std::size_t>(ny),
      static_cast<std::size_t>(nx)};
  check_nc(nc_put_vara_float(file_id, var_id, start, count, values.data()), "write 3D var");
}

void put_2d(
    const int file_id,
    const int var_id,
    const int time_index,
    const int ny,
    const int nx,
    const float offset) {
  std::vector<float> values(static_cast<std::size_t>(ny * nx));
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      values[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
             static_cast<std::size_t>(i)] = value_2d(time_index, j, i, offset);
    }
  }
  const std::size_t start[3] = {static_cast<std::size_t>(time_index), 0, 0};
  const std::size_t count[3] = {1, static_cast<std::size_t>(ny), static_cast<std::size_t>(nx)};
  check_nc(nc_put_vara_float(file_id, var_id, start, count, values.data()), "write 2D var");
}

void fill_3d(
    tywrf::FieldStorage3D<float>& field,
    const int time_index,
    const float offset) {
  const auto layout = field.layout();
  auto view = field.view();
  for (int j = 0; j < layout.active_ny(); ++j) {
    for (int k = 0; k < layout.active_nz(); ++k) {
      for (int i = 0; i < layout.active_nx(); ++i) {
        view(i + layout.i_begin(), j + layout.j_begin(), k + layout.k_begin()) =
            value_3d(time_index, k, j, i, offset);
      }
    }
  }
}

void fill_2d(
    tywrf::FieldStorage2D<float>& field,
    const int time_index,
    const float offset) {
  const auto layout = field.layout();
  auto view = field.view();
  for (int j = 0; j < layout.active_ny(); ++j) {
    for (int i = 0; i < layout.active_nx(); ++i) {
      view(i + layout.i_begin(), j + layout.j_begin()) =
          value_2d(time_index, j, i, offset);
    }
  }
}

void create_synthetic_wrfout(const std::filesystem::path& path) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create synthetic wrfout");

  int time_dim = -1;
  int bottom_top_dim = -1;
  int bottom_top_stag_dim = -1;
  int south_north_dim = -1;
  int south_north_stag_dim = -1;
  int west_east_dim = -1;
  int west_east_stag_dim = -1;
  int date_str_len_dim = -1;
  check_nc(nc_def_dim(file_id, "Time", NC_UNLIMITED, &time_dim), "define Time");
  check_nc(nc_def_dim(file_id, "DateStrLen", 19, &date_str_len_dim), "define DateStrLen");
  check_nc(nc_def_dim(file_id, "bottom_top", 2, &bottom_top_dim), "define bottom_top");
  check_nc(nc_def_dim(file_id, "bottom_top_stag", 3, &bottom_top_stag_dim), "define bottom_top_stag");
  check_nc(nc_def_dim(file_id, "south_north", 3, &south_north_dim), "define south_north");
  check_nc(
      nc_def_dim(file_id, "south_north_stag", 4, &south_north_stag_dim),
      "define south_north_stag");
  check_nc(nc_def_dim(file_id, "west_east", 4, &west_east_dim), "define west_east");
  check_nc(
      nc_def_dim(file_id, "west_east_stag", 5, &west_east_stag_dim),
      "define west_east_stag");

  int times_var = -1;
  int u_var = -1;
  int t_var = -1;
  int mu_var = -1;
  int xlat_var = -1;
  int xlong_var = -1;
  int hgt_var = -1;
  const int times_dims[2] = {time_dim, date_str_len_dim};
  const int u_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_stag_dim};
  const int t_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_dim};
  const int mu_dims[3] = {time_dim, south_north_dim, west_east_dim};
  const int static_dims[3] = {time_dim, south_north_dim, west_east_dim};
  check_nc(nc_def_var(file_id, "Times", NC_CHAR, 2, times_dims, &times_var), "define Times");
  check_nc(nc_def_var(file_id, "U", NC_FLOAT, 4, u_dims, &u_var), "define U");
  check_nc(nc_def_var(file_id, "T", NC_FLOAT, 4, t_dims, &t_var), "define T");
  check_nc(nc_def_var(file_id, "MU", NC_FLOAT, 3, mu_dims, &mu_var), "define MU");
  check_nc(nc_def_var(file_id, "XLAT", NC_FLOAT, 3, static_dims, &xlat_var), "define XLAT");
  check_nc(nc_def_var(file_id, "XLONG", NC_FLOAT, 3, static_dims, &xlong_var), "define XLONG");
  check_nc(nc_def_var(file_id, "HGT", NC_FLOAT, 3, static_dims, &hgt_var), "define HGT");
  check_nc(nc_enddef(file_id), "end definitions");

  for (int time_index = 0; time_index < 2; ++time_index) {
    put_times(file_id, times_var, time_index, synthetic_time_string(time_index));
    put_3d(file_id, u_var, time_index, 2, 3, 5, 10000.0F);
    put_3d(file_id, t_var, time_index, 2, 3, 4, 20000.0F);
    put_2d(file_id, mu_var, time_index, 3, 4, 30000.0F);
    put_2d(file_id, xlat_var, time_index, 3, 4, 40000.0F);
    put_2d(file_id, xlong_var, time_index, 3, 4, 50000.0F);
    put_2d(file_id, hgt_var, time_index, 3, 4, 60000.0F);
  }

  check_nc(nc_close(file_id), "close synthetic wrfout");
}

std::string read_times(const std::filesystem::path& path, const int time_index) {
  int file_id = -1;
  check_nc(nc_open(path.string().c_str(), NC_NOWRITE, &file_id), "open Times reader");

  int times_var = -1;
  int date_str_len_dim = -1;
  check_nc(nc_inq_varid(file_id, "Times", &times_var), "inquire Times");
  check_nc(nc_inq_dimid(file_id, "DateStrLen", &date_str_len_dim), "inquire DateStrLen");
  std::size_t date_str_len = 0;
  check_nc(nc_inq_dimlen(file_id, date_str_len_dim, &date_str_len), "read DateStrLen");

  std::vector<char> buffer(date_str_len, '\0');
  const std::size_t start[2] = {static_cast<std::size_t>(time_index), 0};
  const std::size_t count[2] = {1, date_str_len};
  check_nc(nc_get_vara_text(file_id, times_var, start, count, buffer.data()), "read Times");
  check_nc(nc_close(file_id), "close Times reader");
  return std::string(buffer.begin(), buffer.end());
}

float read_2d_value(
    const std::filesystem::path& path,
    const std::string_view name,
    const int time_index,
    const int j,
    const int i) {
  int file_id = -1;
  check_nc(nc_open(path.string().c_str(), NC_NOWRITE, &file_id), "open static reader");

  int var_id = -1;
  check_nc(nc_inq_varid(file_id, std::string(name).c_str(), &var_id), "inquire static var");
  float value = 0.0F;
  const std::size_t start[3] = {
      static_cast<std::size_t>(time_index), static_cast<std::size_t>(j), static_cast<std::size_t>(i)};
  const std::size_t count[3] = {1, 1, 1};
  check_nc(nc_get_vara_float(file_id, var_id, start, count, &value), "read static var");
  check_nc(nc_close(file_id), "close static reader");
  return value;
}

void create_broken_wrfout(const std::filesystem::path& path) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create broken wrfout");

  int time_dim = -1;
  int bottom_top_dim = -1;
  int bottom_top_stag_dim = -1;
  int south_north_dim = -1;
  int south_north_stag_dim = -1;
  int west_east_dim = -1;
  int west_east_stag_dim = -1;
  check_nc(nc_def_dim(file_id, "Time", NC_UNLIMITED, &time_dim), "define Time");
  check_nc(nc_def_dim(file_id, "bottom_top", 2, &bottom_top_dim), "define bottom_top");
  check_nc(nc_def_dim(file_id, "bottom_top_stag", 3, &bottom_top_stag_dim), "define bottom_top_stag");
  check_nc(nc_def_dim(file_id, "south_north", 3, &south_north_dim), "define south_north");
  check_nc(
      nc_def_dim(file_id, "south_north_stag", 4, &south_north_stag_dim),
      "define south_north_stag");
  check_nc(nc_def_dim(file_id, "west_east", 4, &west_east_dim), "define west_east");
  check_nc(
      nc_def_dim(file_id, "west_east_stag", 5, &west_east_stag_dim),
      "define west_east_stag");

  int t_var = -1;
  const int wrong_t_dims[4] = {time_dim, bottom_top_dim, west_east_dim, south_north_dim};
  check_nc(nc_def_var(file_id, "T", NC_FLOAT, 4, wrong_t_dims, &t_var), "define broken T");
  check_nc(nc_close(file_id), "close broken wrfout");
}

bool message_contains(const tywrf::io::WrfStateIoError& error, const std::string_view needle) {
  return std::string_view(error.what()).find(needle) != std::string_view::npos;
}

int run_synthetic_loader_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_wrf_state_io_test.nc";
  std::filesystem::remove(path);
  create_synthetic_wrfout(path);

  const auto grid = tywrf::io::derive_grid_from_wrf_file(path, tywrf::uniform_halo_3d(1));
  assert(grid.config().mass_nx == 4);
  assert(grid.config().mass_ny == 3);
  assert(grid.config().mass_nz == 2);
  assert(grid.config().full_nz == 3);

  tywrf::State<float> state(grid);
  tywrf::io::load_wrf_state(path, state, {.time_index = 1, .variables = {"U", "T", "MU"}});

  const auto u = state.u.view();
  const auto t = state.t.view();
  const auto mu = state.mu.view();
  const auto u_layout = state.u.layout();
  const auto t_layout = state.t.layout();
  const auto mu_layout = state.mu.layout();

  assert(u(
             u_layout.i_begin() + 4,
             u_layout.j_begin() + 2,
             u_layout.k_begin() + 1) == value_3d(1, 1, 2, 4, 10000.0F));
  assert(t(
             t_layout.i_begin() + 3,
             t_layout.j_begin() + 2,
             t_layout.k_begin() + 1) == value_3d(1, 1, 2, 3, 20000.0F));
  assert(mu(mu_layout.i_begin() + 3, mu_layout.j_begin() + 2) ==
         value_2d(1, 2, 3, 30000.0F));
  assert(state.t.data()[t_layout.index(t_layout.i_begin(), t_layout.j_begin(), t_layout.k_begin())] ==
         value_3d(1, 0, 0, 0, 20000.0F));
  assert(state.t.data()[t_layout.index(t_layout.i_begin() - 1, t_layout.j_begin(), t_layout.k_begin())] ==
         0.0F);

  std::filesystem::remove(path);
  return 0;
}

int run_synthetic_writer_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_wrf_state_io_writer_test.nc";
  const auto all_path =
      std::filesystem::temp_directory_path() / "tywrf_wrf_state_io_writer_all_test.nc";
  std::filesystem::remove(path);
  std::filesystem::remove(all_path);

  const int write_time_index = 1;
  const tywrf::Grid grid({4, 3, 2, 3, tywrf::uniform_halo_3d(1)});
  tywrf::State<float> state(grid);
  fill_3d(state.u, write_time_index, 11000.0F);
  fill_3d(state.t, write_time_index, 22000.0F);
  fill_2d(state.mu, write_time_index, 33000.0F);

  tywrf::io::write_wrf_state(
      path,
      state,
      {.time_index = static_cast<std::size_t>(write_time_index), .variables = {"U", "T", "MU"}});

  const auto writable = tywrf::io::wrf_state_writable_field_names();
  assert(std::find(writable.begin(), writable.end(), "U") != writable.end());
  assert(std::find(writable.begin(), writable.end(), "XLAT") == writable.end());
  const auto missing = tywrf::io::wrf_state_writer_missing_field_names();
  assert(std::find(missing.begin(), missing.end(), "XLAT") != missing.end());

  const auto written_grid = tywrf::io::derive_grid_from_wrf_file(path);
  assert(written_grid.config().mass_nx == 4);
  assert(written_grid.config().mass_ny == 3);
  assert(written_grid.config().mass_nz == 2);
  assert(written_grid.config().full_nz == 3);

  tywrf::State<float> roundtrip(written_grid);
  tywrf::io::load_wrf_state(
      path,
      roundtrip,
      {.time_index = static_cast<std::size_t>(write_time_index), .variables = {"U", "T", "MU"}});

  const auto u = roundtrip.u.view();
  const auto t = roundtrip.t.view();
  const auto mu = roundtrip.mu.view();
  assert(u(4, 2, 1) == value_3d(write_time_index, 1, 2, 4, 11000.0F));
  assert(t(3, 2, 1) == value_3d(write_time_index, 1, 2, 3, 22000.0F));
  assert(mu(3, 2) == value_2d(write_time_index, 2, 3, 33000.0F));

  tywrf::io::write_wrf_state(all_path, state);
  const auto all_grid = tywrf::io::derive_grid_from_wrf_file(all_path);
  tywrf::State<float> all_roundtrip(all_grid);
  tywrf::io::load_wrf_state(
      all_path,
      all_roundtrip,
      {.time_index = 0, .variables = tywrf::io::wrf_state_writable_field_names()});

  std::filesystem::remove(path);
  std::filesystem::remove(all_path);
  return 0;
}

int run_template_writer_test() {
  const auto template_path =
      std::filesystem::temp_directory_path() / "tywrf_wrf_state_io_template_source_test.nc";
  const auto output_path =
      std::filesystem::temp_directory_path() / "tywrf_wrf_state_io_template_writer_test.nc";
  std::filesystem::remove(template_path);
  std::filesystem::remove(output_path);
  create_synthetic_wrfout(template_path);

  const int output_time_index = 0;
  const int template_time_index = 1;
  const tywrf::Grid grid({4, 3, 2, 3, tywrf::uniform_halo_3d(1)});
  tywrf::State<float> state(grid);
  fill_3d(state.u, output_time_index, 71000.0F);
  fill_3d(state.t, output_time_index, 72000.0F);
  fill_2d(state.mu, output_time_index, 73000.0F);

  const std::string cycle_end = "2025-07-26_12:00:00";
  tywrf::io::write_wrf_state(
      output_path,
      state,
      {
          .time_index = static_cast<std::size_t>(output_time_index),
          .variables = {"Times", "XLAT", "XLONG", "HGT", "U", "T", "MU"},
          .template_path = template_path,
          .template_time_index = static_cast<std::size_t>(template_time_index),
          .times_value = cycle_end,
      });

  assert(read_times(output_path, output_time_index) == cycle_end);
  assert(
      read_2d_value(output_path, "XLAT", output_time_index, 2, 3) ==
      value_2d(template_time_index, 2, 3, 40000.0F));
  assert(
      read_2d_value(output_path, "XLONG", output_time_index, 1, 2) ==
      value_2d(template_time_index, 1, 2, 50000.0F));
  assert(
      read_2d_value(output_path, "HGT", output_time_index, 0, 1) ==
      value_2d(template_time_index, 0, 1, 60000.0F));

  const auto written_grid = tywrf::io::derive_grid_from_wrf_file(output_path);
  tywrf::State<float> roundtrip(written_grid);
  tywrf::io::load_wrf_state(
      output_path,
      roundtrip,
      {.time_index = static_cast<std::size_t>(output_time_index), .variables = {"U", "T", "MU"}});
  assert(roundtrip.u.view()(4, 2, 1) == value_3d(output_time_index, 1, 2, 4, 71000.0F));
  assert(roundtrip.t.view()(3, 2, 1) == value_3d(output_time_index, 1, 2, 3, 72000.0F));
  assert(roundtrip.mu.view()(3, 2) == value_2d(output_time_index, 2, 3, 73000.0F));

  std::filesystem::remove(template_path);
  std::filesystem::remove(output_path);
  return 0;
}

int run_error_contract_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_wrf_state_io_broken_test.nc";
  std::filesystem::remove(path);
  create_broken_wrfout(path);

  const auto grid = tywrf::io::derive_grid_from_wrf_file(path);
  tywrf::State<float> state(grid);

  try {
    tywrf::io::load_wrf_state(path, state, {.time_index = 0, .variables = {"T"}});
    std::cerr << "broken T layout did not throw\n";
    return 1;
  } catch (const tywrf::io::WrfStateIoError& error) {
    if (!message_contains(error, "has dimensions")) {
      std::cerr << "broken T layout error did not describe dimension mismatch: " << error.what()
                << '\n';
      return 1;
    }
  }

  try {
    tywrf::io::write_wrf_state(path, state, {.time_index = 0, .variables = {"XLAT"}});
    std::cerr << "unsupported writer variable did not throw\n";
    return 1;
  } catch (const tywrf::io::WrfStateIoError& error) {
    if (!message_contains(error, "unsupported WRF state writer variable selection") ||
        !message_contains(error, "not yet written")) {
      std::cerr << "unsupported writer variable error was not explicit: " << error.what()
                << '\n';
      return 1;
    }
  }

  std::filesystem::remove(path);
  return 0;
}

int run_reference_smoke_test() {
  const auto root = reference_dir();
  const auto wrfinput_d01 = root / "wrfinput_d01";
  const auto d01 = root / "wrfout_d01_2025-07-26_00:00:00";
  const auto d02 = root / "wrfout_d02_2025-07-26_00:00:00";
  if (!std::filesystem::exists(d01) || !std::filesystem::exists(d02)) {
    std::cout << "Skipping WRF state I/O reference smoke; wrfout d01/d02 not present under "
              << root << '\n';
    return 0;
  }

  const auto d01_grid = tywrf::io::derive_grid_from_wrf_file(d01);
  const auto d02_grid = tywrf::io::derive_grid_from_wrf_file(d02);
  assert(d01_grid.config().mass_nx == 265);
  assert(d01_grid.config().mass_ny == 429);
  assert(d01_grid.config().mass_nz == 59);
  assert(d01_grid.config().full_nz == 60);
  assert(d02_grid.config().mass_nx == 210);
  assert(d02_grid.config().mass_ny == 210);
  assert(d02_grid.config().mass_nz == 59);
  assert(d02_grid.config().full_nz == 60);

  if (std::filesystem::exists(wrfinput_d01)) {
    const auto input_grid = tywrf::io::derive_grid_from_wrf_file(wrfinput_d01);
    assert(input_grid.config().mass_nx == 265);
    assert(input_grid.config().mass_ny == 429);
    assert(input_grid.config().mass_nz == 59);
    assert(input_grid.config().full_nz == 60);
  }

  tywrf::State<float> d02_state(d02_grid);
  tywrf::io::load_wrf_state(d02, d02_state, {.time_index = 0, .variables = {"MU"}});
  return 0;
}

}  // namespace

int main() {
  try {
    if (const int status = run_synthetic_loader_test(); status != 0) {
      return status;
    }
    if (const int status = run_synthetic_writer_test(); status != 0) {
      return status;
    }
    if (const int status = run_template_writer_test(); status != 0) {
      return status;
    }
    if (const int status = run_error_contract_test(); status != 0) {
      return status;
    }
    if (const int status = run_reference_smoke_test(); status != 0) {
      return status;
    }
  } catch (const std::exception& error) {
    std::cerr << "WRF state I/O test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
