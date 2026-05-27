#include <netcdf.h>
#include <sys/wait.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
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

enum class PTopMode {
  global_attribute,
  time_variable,
};

struct SyntheticOptions {
  PTopMode p_top_mode = PTopMode::time_variable;
  bool include_alb = true;
  bool bad_alb_shape = false;
};

struct CommandResult {
  int exit_code = 0;
  std::string output;
};

void check_nc(const int status, const std::string_view operation) {
  if (status == NC_NOERR) {
    return;
  }
  throw std::runtime_error(std::string(operation) + ": " + nc_strerror(status));
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
               static_cast<std::size_t>(i)] =
            100.0F * static_cast<float>(k) + 10.0F * static_cast<float>(j) +
            static_cast<float>(i);
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

  int p_top_var = -1;
  if (options.p_top_mode == PTopMode::global_attribute) {
    const float p_top = 5000.0F;
    check_nc(nc_put_att_float(file_id, NC_GLOBAL, "P_TOP", NC_FLOAT, 1, &p_top), "P_TOP attr");
  } else {
    const int dims[1] = {time_dim};
    check_nc(nc_def_var(file_id, "P_TOP", NC_FLOAT, 1, dims, &p_top_var), "define P_TOP");
  }

  int c3f_var = -1;
  int c4f_var = -1;
  int c3h_var = -1;
  int c4h_var = -1;
  const int full_dims[2] = {time_dim, bottom_top_stag_dim};
  const int mass_dims[2] = {time_dim, bottom_top_dim};
  check_nc(nc_def_var(file_id, "C3F", NC_FLOAT, 2, full_dims, &c3f_var), "define C3F");
  check_nc(nc_def_var(file_id, "C4F", NC_FLOAT, 2, full_dims, &c4f_var), "define C4F");
  check_nc(nc_def_var(file_id, "C3H", NC_FLOAT, 2, mass_dims, &c3h_var), "define C3H");
  check_nc(nc_def_var(file_id, "C4H", NC_FLOAT, 2, mass_dims, &c4h_var), "define C4H");

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

  if (options.p_top_mode == PTopMode::time_variable) {
    const float p_top = 5200.0F;
    const std::size_t start[1] = {0};
    const std::size_t count[1] = {1};
    check_nc(nc_put_vara_float(file_id, p_top_var, start, count, &p_top), "write P_TOP");
  }

  write_vector(file_id, c3f_var, 3, 10.0F);
  write_vector(file_id, c4f_var, 3, 20.0F);
  write_vector(file_id, c3h_var, 2, 30.0F);
  write_vector(file_id, c4h_var, 2, 40.0F);
  if (options.include_alb) {
    write_alb(file_id, alb_var, 2, 3, options.bad_alb_shape ? 5 : 4);
  }

  check_nc(nc_close(file_id), "close synthetic file");
}

[[nodiscard]] std::string shell_quote(const std::filesystem::path& path) {
  return "'" + path.string() + "'";
}

[[nodiscard]] CommandResult run_probe(
    const std::filesystem::path& executable,
    const std::filesystem::path& input,
    const std::string_view extra_args = {}) {
  const std::string command =
      shell_quote(executable) + " --input " + shell_quote(input) +
      " --mass-nx 4 --mass-ny 3 --mass-nz 2 --full-nz 3 --halo 1 --pretty " +
      std::string(extra_args) + " 2>&1";

  FILE* pipe = popen(command.c_str(), "r");
  if (pipe == nullptr) {
    throw std::runtime_error("failed to run pressure-refresh probe");
  }

  std::string output;
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
    output += buffer;
  }
  const int status = pclose(pipe);
  int exit_code = status;
  if (WIFEXITED(status)) {
    exit_code = WEXITSTATUS(status);
  }
  return {.exit_code = exit_code, .output = output};
}

void require_contains(const std::string& haystack, const std::string_view needle) {
  if (haystack.find(needle) == std::string::npos) {
    throw std::runtime_error("expected output to contain: " + std::string(needle) +
                             "\noutput was:\n" + haystack);
  }
}

void require_not_contains(const std::string& haystack, const std::string_view needle) {
  if (haystack.find(needle) != std::string::npos) {
    throw std::runtime_error("expected output not to contain: " + std::string(needle) +
                             "\noutput was:\n" + haystack);
  }
}

void run_complete_file_test(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto path = root / "complete_pressure_refresh_inputs.nc";
  create_synthetic_file(path, {});

  const auto result = run_probe(executable, path);
  assert(result.exit_code == 0);
  require_contains(result.output, "\"status\": \"ok\"");
  require_contains(result.output, "\"input_file\": \"" + path.string() + "\"");
  require_contains(result.output, "\"time_index\": 0");
  require_contains(result.output, "\"mass_levels\": 2");
  require_contains(result.output, "\"full_levels\": 3");
  require_contains(result.output, "\"alb_points\": 24");
  require_contains(result.output, "\"alb\": [4,3,2]");
  require_contains(result.output, "\"missing_names\": []");
  require_contains(result.output, "\"p_top_present\": true");
  require_contains(result.output, "\"p_top_source\": \"time_variable\"");
  require_contains(result.output, "\"alb_available\": true");
  require_contains(result.output, "\"alb_loaded\": true");
  require_contains(result.output, "\"diagnostic_only\": true");
  require_contains(result.output, "\"calls_pressure_refresh_compute\": false");
  require_contains(result.output, "\"reads_or_writes_p\": false");
}

void run_history_like_missing_alb_test(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto path = root / "history_like_without_alb.nc";
  create_synthetic_file(path, {.p_top_mode = PTopMode::global_attribute, .include_alb = false});

  const auto result = run_probe(executable, path);
  assert(result.exit_code == 1);
  require_contains(result.output, "\"status\": \"missing\"");
  require_contains(result.output, "\"missing_names\": [");
  require_contains(result.output, "\"ALB\"");
  require_contains(result.output, "\"p_top_present\": true");
  require_contains(result.output, "\"p_top_source\": \"global_attribute\"");
  require_contains(result.output, "\"alb_available\": false");
  require_contains(result.output, "\"alb_loaded\": false");
  require_contains(result.output, "\"diagnostic_only\": true");
  require_contains(result.output, "\"calls_pressure_refresh_compute\": false");
  require_contains(result.output, "\"reads_or_writes_p\": false");
  require_not_contains(result.output, "segmentation");
}

void run_bad_shape_test(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto path = root / "bad_alb_shape.nc";
  create_synthetic_file(path, {.bad_alb_shape = true});

  const auto result = run_probe(executable, path);
  assert(result.exit_code == 1);
  require_contains(result.output, "\"status\": \"error\"");
  require_contains(result.output, "ALB");
  require_contains(result.output, "west_east=4");
  require_contains(result.output, "\"diagnostic_only\": true");
  require_contains(result.output, "\"calls_pressure_refresh_compute\": false");
  require_contains(result.output, "\"reads_or_writes_p\": false");
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: pressure_refresh_probe_tool_test <tywrf_pressure_refresh_probe path>\n";
      return 2;
    }

    const auto executable = std::filesystem::path(argv[1]);
    const auto root = std::filesystem::temp_directory_path() / "tywrf_pressure_refresh_probe_test";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);

    run_complete_file_test(executable, root);
    run_history_like_missing_alb_test(executable, root);
    run_bad_shape_test(executable, root);

    std::filesystem::remove_all(root);
  } catch (const std::exception& error) {
    std::cerr << "pressure-refresh probe tool test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
