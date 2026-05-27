#include "tywrf/io/forcing_io.hpp"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
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

std::string synthetic_time_string(const int time_index) {
  return time_index == 0 ? "2025-07-26_00:00:00" : "2025-07-26_06:00:00";
}

std::size_t product(const std::vector<std::size_t>& shape) {
  std::size_t count = 1;
  for (const auto extent : shape) {
    count *= extent;
  }
  return count;
}

float sample_value(const int time_index, const std::size_t linear_index, const float base) {
  return base + 1000.0F * static_cast<float>(time_index) +
         static_cast<float>(linear_index);
}

int define_dim(
    const int file_id,
    const std::string_view name,
    const std::size_t length) {
  int dim_id = -1;
  check_nc(
      nc_def_dim(
          file_id,
          std::string(name).c_str(),
          length == 0 ? NC_UNLIMITED : length,
          &dim_id),
      "define dimension");
  return dim_id;
}

int define_float_var(
    const int file_id,
    const std::string_view name,
    const std::vector<int>& dims) {
  int var_id = -1;
  check_nc(
      nc_def_var(
          file_id,
          std::string(name).c_str(),
          NC_FLOAT,
          static_cast<int>(dims.size()),
          dims.data(),
          &var_id),
      "define float variable");
  return var_id;
}

void write_times(
    const int file_id,
    const int var_id,
    const int time_index,
    const std::string_view value) {
  std::vector<char> buffer(19, ' ');
  std::copy(value.begin(), value.end(), buffer.begin());
  const std::size_t start[2] = {static_cast<std::size_t>(time_index), 0};
  const std::size_t count[2] = {1, buffer.size()};
  check_nc(nc_put_vara_text(file_id, var_id, start, count, buffer.data()), "write Times");
}

void write_float_time_slice(
    const int file_id,
    const int var_id,
    const int time_index,
    const std::vector<std::size_t>& slice_shape,
    const float base) {
  std::vector<float> values(product(slice_shape), 0.0F);
  for (std::size_t index = 0; index < values.size(); ++index) {
    values[index] = sample_value(time_index, index, base);
  }

  std::vector<std::size_t> start(slice_shape.size() + 1, 0);
  std::vector<std::size_t> count(slice_shape.size() + 1, 1);
  start.front() = static_cast<std::size_t>(time_index);
  for (std::size_t dim = 0; dim < slice_shape.size(); ++dim) {
    count.at(dim + 1) = slice_shape.at(dim);
  }
  check_nc(
      nc_put_vara_float(file_id, var_id, start.data(), count.data(), values.data()),
      "write float slice");
}

void create_synthetic_wrfbdy(const std::filesystem::path& path) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create wrfbdy");

  const int time_dim = define_dim(file_id, "Time", 0);
  const int date_dim = define_dim(file_id, "DateStrLen", 19);
  const int bdy_width_dim = define_dim(file_id, "bdy_width", 2);
  const int bottom_top_dim = define_dim(file_id, "bottom_top", 2);
  const int south_north_dim = define_dim(file_id, "south_north", 3);
  const int west_east_dim = define_dim(file_id, "west_east", 5);

  int times_var = -1;
  const int times_dims[2] = {time_dim, date_dim};
  check_nc(nc_def_var(file_id, "Times", NC_CHAR, 2, times_dims, &times_var), "define Times");
  const int u_bxs_var =
      define_float_var(file_id, "U_BXS", {time_dim, bdy_width_dim, bottom_top_dim, south_north_dim});
  const int u_btxs_var =
      define_float_var(file_id, "U_BTXS", {time_dim, bdy_width_dim, bottom_top_dim, south_north_dim});
  const int mu_bys_var =
      define_float_var(file_id, "MU_BYS", {time_dim, bdy_width_dim, west_east_dim});
  const int qvapor_btye_var =
      define_float_var(file_id, "QVAPOR_BTYE", {time_dim, bdy_width_dim, bottom_top_dim, west_east_dim});

  check_nc(nc_enddef(file_id), "end wrfbdy definitions");
  for (int time_index = 0; time_index < 2; ++time_index) {
    write_times(file_id, times_var, time_index, synthetic_time_string(time_index));
    write_float_time_slice(file_id, u_bxs_var, time_index, {2, 2, 3}, 100.0F);
    write_float_time_slice(file_id, u_btxs_var, time_index, {2, 2, 3}, 200.0F);
    write_float_time_slice(file_id, mu_bys_var, time_index, {2, 5}, 300.0F);
    write_float_time_slice(file_id, qvapor_btye_var, time_index, {2, 2, 5}, 400.0F);
  }
  check_nc(nc_close(file_id), "close wrfbdy");
}

void create_synthetic_wrffdda(const std::filesystem::path& path) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create wrffdda");

  const int time_dim = define_dim(file_id, "Time", 0);
  const int date_dim = define_dim(file_id, "DateStrLen", 19);
  const int bottom_top_dim = define_dim(file_id, "bottom_top", 2);
  const int one_stag_dim = define_dim(file_id, "one_stag", 1);
  const int south_north_dim = define_dim(file_id, "south_north", 3);
  const int west_east_dim = define_dim(file_id, "west_east", 5);

  int times_var = -1;
  const int times_dims[2] = {time_dim, date_dim};
  check_nc(nc_def_var(file_id, "Times", NC_CHAR, 2, times_dims, &times_var), "define Times");
  const int u_old_var =
      define_float_var(file_id, "U_NDG_OLD", {time_dim, bottom_top_dim, south_north_dim, west_east_dim});
  const int u_new_var =
      define_float_var(file_id, "U_NDG_NEW", {time_dim, bottom_top_dim, south_north_dim, west_east_dim});
  const int q_new_var =
      define_float_var(file_id, "Q_NDG_NEW", {time_dim, bottom_top_dim, south_north_dim, west_east_dim});
  const int mu_old_var =
      define_float_var(file_id, "MU_NDG_OLD", {time_dim, one_stag_dim, south_north_dim, west_east_dim});

  check_nc(nc_enddef(file_id), "end wrffdda definitions");
  for (int time_index = 0; time_index < 2; ++time_index) {
    write_times(file_id, times_var, time_index, synthetic_time_string(time_index));
    write_float_time_slice(file_id, u_old_var, time_index, {2, 3, 5}, 500.0F);
    write_float_time_slice(file_id, u_new_var, time_index, {2, 3, 5}, 600.0F);
    write_float_time_slice(file_id, q_new_var, time_index, {2, 3, 5}, 700.0F);
    write_float_time_slice(file_id, mu_old_var, time_index, {1, 3, 5}, 800.0F);
  }
  check_nc(nc_close(file_id), "close wrffdda");
}

bool contains(const std::vector<std::string>& values, const std::string_view needle) {
  return std::any_of(values.begin(), values.end(), [&](const std::string& value) {
    return value == needle;
  });
}

template <typename Fn>
bool throws_forcing_error(Fn&& fn) {
  try {
    fn();
  } catch (const tywrf::io::ForcingIoError&) {
    return true;
  }
  return false;
}

int run_boundary_reader_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_forcing_wrfbdy_test.nc";
  std::filesystem::remove(path);
  create_synthetic_wrfbdy(path);

  const tywrf::io::KrosaForcingReader reader(path, tywrf::io::KrosaForcingKind::boundary);
  if (reader.time_count() != 2 || reader.date_string_length() != 19) {
    std::cerr << "synthetic wrfbdy time metadata mismatch\n";
    return 1;
  }
  if (reader.read_time_string(1) != "2025-07-26_06:00:00") {
    std::cerr << "synthetic wrfbdy Times value mismatch\n";
    return 1;
  }
  if (reader.read_char_time_slice("Times", 1) != "2025-07-26_06:00:00") {
    std::cerr << "synthetic wrfbdy generic char time-slice mismatch\n";
    return 1;
  }
  const auto time_strings = reader.read_time_strings();
  if (time_strings.size() != 2 || time_strings.front() != "2025-07-26_00:00:00" ||
      time_strings.back() != "2025-07-26_06:00:00") {
    std::cerr << "synthetic wrfbdy Times vector mismatch\n";
    return 1;
  }
  if (!reader.has_variable("U_BXS") || reader.has_variable("V_BXS")) {
    std::cerr << "synthetic wrfbdy variable existence query mismatch\n";
    return 1;
  }

  const auto metadata = reader.variable_metadata("U_BXS");
  if (metadata.dimensions != std::vector<std::string>{"Time", "bdy_width", "bottom_top", "south_north"} ||
      metadata.shape != std::vector<std::size_t>{2, 2, 2, 3} ||
      metadata.slice_shape != std::vector<std::size_t>{2, 2, 3} ||
      metadata.values_per_time_slice != 12) {
    std::cerr << "synthetic U_BXS metadata mismatch\n";
    return 1;
  }

  const auto slice = reader.read_float_time_slice("U_BXS", 1);
  if (slice.values.size() != 12 || slice.values.front() != sample_value(1, 0, 100.0F) ||
      slice.values.back() != sample_value(1, 11, 100.0F)) {
    std::cerr << "synthetic U_BXS slice values are not contiguous as expected\n";
    return 1;
  }

  const auto missing = reader.missing_required_variables();
  if (!contains(missing, "V_BXS") || contains(missing, "U_BXS")) {
    std::cerr << "synthetic wrfbdy missing-required list mismatch\n";
    return 1;
  }

  std::filesystem::remove(path);
  return 0;
}

int run_nudging_reader_test() {
  const auto path = std::filesystem::temp_directory_path() / "tywrf_forcing_wrffdda_test.nc";
  std::filesystem::remove(path);
  create_synthetic_wrffdda(path);

  const auto required = tywrf::io::krosa_wrffdda_required_variable_names();
  if (required.size() != 13 || !contains(required, "U_NDG_OLD") ||
      !contains(required, "MU_NDG_NEW")) {
    std::cerr << "KROSA wrffdda required variable list is incomplete\n";
    return 1;
  }

  const tywrf::io::KrosaForcingReader reader(
      path,
      tywrf::io::KrosaForcingKind::spectral_nudging);
  const auto metadata = reader.variable_metadata("MU_NDG_OLD");
  if (metadata.dimensions != std::vector<std::string>{"Time", "one_stag", "south_north", "west_east"} ||
      metadata.shape != std::vector<std::size_t>{2, 1, 3, 5} ||
      metadata.values_per_time_slice != 15) {
    std::cerr << "synthetic MU_NDG_OLD metadata mismatch\n";
    return 1;
  }

  const auto slice = reader.read_float_time_slice("U_NDG_NEW", 1);
  if (slice.metadata.slice_shape != std::vector<std::size_t>{2, 3, 5} ||
      slice.values.size() != 30 || slice.values.front() != sample_value(1, 0, 600.0F) ||
      slice.values.back() != sample_value(1, 29, 600.0F)) {
    std::cerr << "synthetic U_NDG_NEW slice values mismatch\n";
    return 1;
  }

  if (!throws_forcing_error([&] {
        const auto ignored = reader.variable_metadata("V_NDG_OLD");
        (void)ignored;
      }) ||
      !throws_forcing_error([&] {
        const auto ignored = reader.read_float_time_slice("Times", 0);
        (void)ignored;
      }) ||
      !throws_forcing_error([&] {
        const auto ignored = reader.read_char_time_slice("U_NDG_OLD", 0);
        (void)ignored;
      }) ||
      !throws_forcing_error([&] {
        const auto ignored = reader.read_float_time_slice("U_NDG_OLD", 2);
        (void)ignored;
      })) {
    std::cerr << "synthetic wrffdda error handling did not throw ForcingIoError\n";
    return 1;
  }

  std::filesystem::remove(path);
  return 0;
}

int run_reference_smoke_test() {
  const auto root = reference_dir();
  const auto wrfbdy = root / "wrfbdy_d01";
  const auto wrffdda = root / "wrffdda_d01";
  if (!std::filesystem::exists(wrfbdy) || !std::filesystem::exists(wrffdda)) {
    std::cout << "Skipping forcing I/O reference smoke; KROSA inputs are not present under "
              << root << '\n';
    return 0;
  }

  const tywrf::io::KrosaForcingReader boundary(
      wrfbdy,
      tywrf::io::KrosaForcingKind::boundary);
  const auto missing_boundary = boundary.missing_required_variables();
  if (!missing_boundary.empty()) {
    std::cerr << "KROSA wrfbdy missing required variable " << missing_boundary.front() << '\n';
    return 1;
  }
  const auto u_bxs = boundary.variable_metadata("U_BXS");
  if (boundary.time_count() != 28 ||
      u_bxs.shape != std::vector<std::size_t>{28, 5, 59, 429}) {
    std::cerr << "KROSA wrfbdy U_BXS shape/time metadata mismatch\n";
    return 1;
  }
  const auto mu_bxs_slice = boundary.read_float_time_slice("MU_BXS", 0);
  if (mu_bxs_slice.values.size() != 5 * 429 ||
      !std::isfinite(mu_bxs_slice.values.front()) ||
      !std::isfinite(mu_bxs_slice.values.back())) {
    std::cerr << "KROSA wrfbdy MU_BXS slice read failed\n";
    return 1;
  }

  const tywrf::io::KrosaForcingReader nudging(
      wrffdda,
      tywrf::io::KrosaForcingKind::spectral_nudging);
  const auto missing_nudging = nudging.missing_required_variables();
  if (!missing_nudging.empty()) {
    std::cerr << "KROSA wrffdda missing required variable " << missing_nudging.front() << '\n';
    return 1;
  }
  const auto q_new = nudging.variable_metadata("Q_NDG_NEW");
  if (nudging.time_count() != 28 ||
      q_new.shape != std::vector<std::size_t>{28, 59, 429, 265}) {
    std::cerr << "KROSA wrffdda Q_NDG_NEW shape/time metadata mismatch\n";
    return 1;
  }
  const auto mu_new_slice = nudging.read_float_time_slice("MU_NDG_NEW", 1);
  if (mu_new_slice.values.size() != 429 * 265 ||
      !std::isfinite(mu_new_slice.values.front()) ||
      !std::isfinite(mu_new_slice.values.back())) {
    std::cerr << "KROSA wrffdda MU_NDG_NEW slice read failed\n";
    return 1;
  }

  std::cout << "Validated KROSA forcing slice reads under " << root << '\n';
  return 0;
}

}  // namespace

int main() {
  try {
    if (const int status = run_boundary_reader_test(); status != 0) {
      return status;
    }
    if (const int status = run_nudging_reader_test(); status != 0) {
      return status;
    }
    if (const int status = run_reference_smoke_test(); status != 0) {
      return status;
    }
  } catch (const std::exception& error) {
    std::cerr << "forcing I/O test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
