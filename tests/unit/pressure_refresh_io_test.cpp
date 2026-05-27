#include "tywrf/io/pressure_refresh_io.hpp"

#include <netcdf.h>

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

enum class PTopMode {
  global_attribute,
  scalar_variable,
  time_variable,
  missing,
};

struct SyntheticOptions {
  PTopMode p_top_mode = PTopMode::time_variable;
  bool include_coefficients = true;
  bool include_alb = true;
  bool bad_alb_shape = false;
};

void check_nc(const int status, const std::string_view operation) {
  if (status == NC_NOERR) {
    return;
  }
  throw std::runtime_error(std::string(operation) + ": " + nc_strerror(status));
}

tywrf::Grid make_grid() {
  return tywrf::Grid({4, 3, 2, 3, tywrf::uniform_halo_3d(1)});
}

float alb_value(const int k, const int j, const int i) {
  return 100.0F * static_cast<float>(k) + 10.0F * static_cast<float>(j) +
         static_cast<float>(i);
}

void write_vector(
    const int file_id,
    const int variable_id,
    const int count,
    const float offset) {
  std::vector<float> values(static_cast<std::size_t>(count), 0.0F);
  for (int index = 0; index < count; ++index) {
    values[static_cast<std::size_t>(index)] = offset + static_cast<float>(index);
  }
  const std::size_t start[2] = {0, 0};
  const std::size_t counts[2] = {1, static_cast<std::size_t>(count)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, counts, values.data()), "write vector");
}

void write_alb(
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
               static_cast<std::size_t>(i)] = alb_value(k, j, i);
      }
    }
  }
  const std::size_t start[4] = {0, 0, 0, 0};
  const std::size_t counts[4] = {
      1,
      static_cast<std::size_t>(nz),
      static_cast<std::size_t>(ny),
      static_cast<std::size_t>(nx)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, counts, values.data()), "write ALB");
}

void create_synthetic_file(const std::filesystem::path& path, const SyntheticOptions& options) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create synthetic file");

  int time_dim = -1;
  int bottom_top_dim = -1;
  int bottom_top_stag_dim = -1;
  int south_north_dim = -1;
  int west_east_dim = -1;
  int west_east_bad_dim = -1;
  check_nc(nc_def_dim(file_id, "Time", NC_UNLIMITED, &time_dim), "define Time");
  check_nc(nc_def_dim(file_id, "bottom_top", 2, &bottom_top_dim), "define bottom_top");
  check_nc(
      nc_def_dim(file_id, "bottom_top_stag", 3, &bottom_top_stag_dim),
      "define bottom_top_stag");
  check_nc(nc_def_dim(file_id, "south_north", 3, &south_north_dim), "define south_north");
  check_nc(nc_def_dim(file_id, "west_east", 4, &west_east_dim), "define west_east");
  check_nc(nc_def_dim(file_id, "west_east_bad", 5, &west_east_bad_dim), "define bad west_east");

  if (options.p_top_mode == PTopMode::global_attribute) {
    const float p_top = 5000.0F;
    check_nc(nc_put_att_float(file_id, NC_GLOBAL, "P_TOP", NC_FLOAT, 1, &p_top), "P_TOP attr");
  }

  int p_top_var = -1;
  if (options.p_top_mode == PTopMode::scalar_variable) {
    check_nc(nc_def_var(file_id, "P_TOP", NC_FLOAT, 0, nullptr, &p_top_var), "define P_TOP");
  } else if (options.p_top_mode == PTopMode::time_variable) {
    const int dims[1] = {time_dim};
    check_nc(nc_def_var(file_id, "P_TOP", NC_FLOAT, 1, dims, &p_top_var), "define P_TOP");
  }

  int c3f_var = -1;
  int c4f_var = -1;
  int c3h_var = -1;
  int c4h_var = -1;
  if (options.include_coefficients) {
    const int full_dims[2] = {time_dim, bottom_top_stag_dim};
    const int mass_dims[2] = {time_dim, bottom_top_dim};
    check_nc(nc_def_var(file_id, "C3F", NC_FLOAT, 2, full_dims, &c3f_var), "define C3F");
    check_nc(nc_def_var(file_id, "C4F", NC_FLOAT, 2, full_dims, &c4f_var), "define C4F");
    check_nc(nc_def_var(file_id, "C3H", NC_FLOAT, 2, mass_dims, &c3h_var), "define C3H");
    check_nc(nc_def_var(file_id, "C4H", NC_FLOAT, 2, mass_dims, &c4h_var), "define C4H");
  }

  int alb_var = -1;
  if (options.include_alb) {
    const int alb_dims[4] = {
        time_dim,
        bottom_top_dim,
        south_north_dim,
        options.bad_alb_shape ? west_east_bad_dim : west_east_dim};
    check_nc(nc_def_var(file_id, "ALB", NC_FLOAT, 4, alb_dims, &alb_var), "define ALB");
  }

  check_nc(nc_enddef(file_id), "end definitions");

  if (options.p_top_mode == PTopMode::scalar_variable) {
    const float p_top = 5100.0F;
    check_nc(nc_put_var_float(file_id, p_top_var, &p_top), "write scalar P_TOP");
  } else if (options.p_top_mode == PTopMode::time_variable) {
    const float p_top = 5200.0F;
    const std::size_t start[1] = {0};
    const std::size_t count[1] = {1};
    check_nc(nc_put_vara_float(file_id, p_top_var, start, count, &p_top), "write P_TOP");
  }

  if (options.include_coefficients) {
    write_vector(file_id, c3f_var, 3, 10.0F);
    write_vector(file_id, c4f_var, 3, 20.0F);
    write_vector(file_id, c3h_var, 2, 30.0F);
    write_vector(file_id, c4h_var, 2, 40.0F);
  }
  if (options.include_alb) {
    write_alb(file_id, alb_var, 2, 3, options.bad_alb_shape ? 5 : 4);
  }

  check_nc(nc_close(file_id), "close synthetic file");
}

bool contains_name(const std::vector<std::string>& names, const std::string_view name) {
  return std::find(names.begin(), names.end(), name) != names.end();
}

bool message_contains(
    const tywrf::io::PressureRefreshIoError& error,
    const std::string_view needle) {
  return std::string_view(error.what()).find(needle) != std::string_view::npos;
}

int run_success_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_pressure_refresh_io_success.nc";
  std::filesystem::remove(path);
  create_synthetic_file(path, {.p_top_mode = PTopMode::time_variable});

  const auto grid = make_grid();
  tywrf::FieldStorage3D<float> alb(grid.mass_layout());
  const auto result = tywrf::io::read_krosa_pressure_refresh_inputs(path, grid, alb);

  assert(result.ok());
  assert(result.metadata.p_top_pa == 5200.0F);
  assert(result.metadata.p_top_source == tywrf::io::PressureRefreshPTopSource::time_variable);
  assert(result.metadata.c3f.size() == 3);
  assert(result.metadata.c4f.size() == 3);
  assert(result.metadata.c3h.size() == 2);
  assert(result.metadata.c4h.size() == 2);
  assert(result.metadata.c3f[2] == 12.0F);
  assert(result.metadata.c4h[1] == 41.0F);
  assert(result.report.alb_loaded);
  assert(result.report.source_path == path);
  assert(result.report.alb_point_count == 24);

  const auto view = alb.view();
  const auto layout = alb.layout();
  assert(
      view(layout.i_begin() + 3, layout.j_begin() + 2, layout.k_begin() + 1) ==
      alb_value(1, 2, 3));
  assert(view(layout.i_begin() - 1, layout.j_begin(), layout.k_begin()) == 0.0F);

  std::filesystem::remove(path);
  return 0;
}

int run_missing_metadata_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_pressure_refresh_io_missing.nc";
  std::filesystem::remove(path);
  create_synthetic_file(
      path,
      {.p_top_mode = PTopMode::missing, .include_coefficients = false, .include_alb = false});

  const auto report = tywrf::io::inspect_krosa_pressure_refresh_inputs(path, make_grid());
  assert(!report.ok());
  assert(contains_name(report.missing_names, "P_TOP"));
  assert(contains_name(report.missing_names, "C3F"));
  assert(contains_name(report.missing_names, "C4F"));
  assert(contains_name(report.missing_names, "C3H"));
  assert(contains_name(report.missing_names, "C4H"));
  assert(contains_name(report.missing_names, "ALB"));
  assert(report.expected_alb_point_count == 24);
  assert(report.alb_point_count == 0);

  std::filesystem::remove(path);
  return 0;
}

int run_history_without_alb_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_pressure_refresh_io_no_alb.nc";
  std::filesystem::remove(path);
  create_synthetic_file(path, {.p_top_mode = PTopMode::time_variable, .include_alb = false});

  tywrf::FieldStorage3D<float> alb(make_grid().mass_layout());
  const auto result = tywrf::io::read_krosa_pressure_refresh_inputs(path, make_grid(), alb);
  assert(!result.ok());
  assert(contains_name(result.report.missing_names, "ALB"));
  assert(result.metadata.p_top_pa == 5200.0F);
  assert(result.metadata.c3f.size() == 3);
  assert(!result.report.alb_loaded);

  std::filesystem::remove(path);
  return 0;
}

int run_bad_alb_shape_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_pressure_refresh_io_bad_alb.nc";
  std::filesystem::remove(path);
  create_synthetic_file(path, {.bad_alb_shape = true});

  tywrf::FieldStorage3D<float> alb(make_grid().mass_layout());
  try {
    (void)tywrf::io::read_krosa_pressure_refresh_inputs(path, make_grid(), alb);
    std::cerr << "bad ALB shape did not throw\n";
    return 1;
  } catch (const tywrf::io::PressureRefreshIoError& error) {
    if (!message_contains(error, "ALB") || !message_contains(error, "west_east=4")) {
      std::cerr << "bad ALB shape error was not explicit: " << error.what() << '\n';
      return 1;
    }
  }

  std::filesystem::remove(path);
  return 0;
}

int run_p_top_attr_and_scalar_variable_test() {
  const auto attr_path =
      std::filesystem::temp_directory_path() / "tywrf_pressure_refresh_io_p_top_attr.nc";
  const auto scalar_path =
      std::filesystem::temp_directory_path() / "tywrf_pressure_refresh_io_p_top_scalar.nc";
  std::filesystem::remove(attr_path);
  std::filesystem::remove(scalar_path);
  create_synthetic_file(attr_path, {.p_top_mode = PTopMode::global_attribute});
  create_synthetic_file(scalar_path, {.p_top_mode = PTopMode::scalar_variable});

  tywrf::FieldStorage3D<float> attr_alb(make_grid().mass_layout());
  tywrf::FieldStorage3D<float> scalar_alb(make_grid().mass_layout());
  const auto attr_result =
      tywrf::io::read_krosa_pressure_refresh_inputs(attr_path, make_grid(), attr_alb);
  const auto scalar_result =
      tywrf::io::read_krosa_pressure_refresh_inputs(scalar_path, make_grid(), scalar_alb);

  assert(attr_result.ok());
  assert(attr_result.metadata.p_top_pa == 5000.0F);
  assert(
      attr_result.metadata.p_top_source ==
      tywrf::io::PressureRefreshPTopSource::global_attribute);
  assert(scalar_result.ok());
  assert(scalar_result.metadata.p_top_pa == 5100.0F);
  assert(
      scalar_result.metadata.p_top_source ==
      tywrf::io::PressureRefreshPTopSource::scalar_variable);

  std::filesystem::remove(attr_path);
  std::filesystem::remove(scalar_path);
  return 0;
}

}  // namespace

int main() {
  try {
    if (run_success_test() != 0 || run_missing_metadata_test() != 0 ||
        run_history_without_alb_test() != 0 || run_bad_alb_shape_test() != 0 ||
        run_p_top_attr_and_scalar_variable_test() != 0) {
      return 1;
    }
  } catch (const std::exception& error) {
    std::cerr << "pressure-refresh I/O test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
