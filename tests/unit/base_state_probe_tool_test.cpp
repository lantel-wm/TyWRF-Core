#include <netcdf.h>
#include <sys/wait.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

constexpr int kMassNx = 4;
constexpr int kMassNy = 3;
constexpr int kMassNz = 2;
constexpr int kFullNz = 3;
constexpr float kPTop = 5000.0F;
constexpr double kRd = 287.0;
constexpr double kCp = 1004.5;
constexpr double kGravity = 9.81;
constexpr double kP00 = 100000.0;
constexpr double kT00 = 290.0;
constexpr double kBaseLapse = 50.0;
constexpr double kTiso = 200.0;
constexpr double kThetaBase = 300.0;
constexpr std::array<float, kFullNz> kC3F = {1.0F, 0.5F, 0.0F};
constexpr std::array<float, kFullNz> kC4F = {0.0F, 0.0F, 0.0F};
constexpr std::array<float, kMassNz> kC3H = {0.75F, 0.25F};
constexpr std::array<float, kMassNz> kC4H = {0.0F, 0.0F};
const std::filesystem::path kKrosaWrfinputD01 =
    "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/model/WRFV4.6.1/run/wrfinput_d01";
const std::filesystem::path kKrosaWrfoutD02 =
    "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF_1h_10min_20260527_172838/wrfout_d02_2025-07-26_00:00:00";

struct SyntheticOptions {
  bool include_pb = true;
  bool include_t_init = true;
  bool include_alb = true;
  bool include_mub = true;
  bool include_phb = true;
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

[[nodiscard]] std::size_t index2(const int j, const int i) {
  return static_cast<std::size_t>(j) * static_cast<std::size_t>(kMassNx) +
         static_cast<std::size_t>(i);
}

[[nodiscard]] std::size_t index3(const int k, const int j, const int i) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(kMassNy) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(kMassNx) +
         static_cast<std::size_t>(i);
}

[[nodiscard]] float terrain(const int j, const int i) {
  return 10.0F * static_cast<float>(j) + static_cast<float>(i);
}

[[nodiscard]] double surface_pressure_pa(const double hgt) {
  const auto t00_over_lapse = kT00 / kBaseLapse;
  const auto discriminant =
      t00_over_lapse * t00_over_lapse -
      (2.0 * kGravity * hgt) / (kBaseLapse * kRd);
  return kP00 * std::exp(-t00_over_lapse + std::sqrt(discriminant));
}

[[nodiscard]] double base_temperature_k(const double pb) {
  return std::max(kT00 + kBaseLapse * std::log(pb / kP00), kTiso);
}

struct SyntheticFields {
  std::vector<float> hgt;
  std::vector<float> pb;
  std::vector<float> t_init;
  std::vector<float> alb;
  std::vector<float> mub;
  std::vector<float> phb;
};

[[nodiscard]] SyntheticFields make_synthetic_fields() {
  SyntheticFields fields;
  fields.hgt.resize(static_cast<std::size_t>(kMassNy * kMassNx));
  fields.pb.resize(static_cast<std::size_t>(kMassNz * kMassNy * kMassNx));
  fields.t_init.resize(static_cast<std::size_t>(kMassNz * kMassNy * kMassNx));
  fields.alb.resize(static_cast<std::size_t>(kMassNz * kMassNy * kMassNx));
  fields.mub.resize(static_cast<std::size_t>(kMassNy * kMassNx));
  fields.phb.resize(static_cast<std::size_t>(kFullNz * kMassNy * kMassNx));

  const auto rd_over_cp = kRd / kCp;
  const auto cvpm = -(kCp - kRd) / kCp;
  for (int j = 0; j < kMassNy; ++j) {
    for (int i = 0; i < kMassNx; ++i) {
      const auto hgt = static_cast<double>(terrain(j, i));
      fields.hgt[index2(j, i)] = static_cast<float>(hgt);
      const auto mub = surface_pressure_pa(hgt) - static_cast<double>(kPTop);
      fields.mub[index2(j, i)] = static_cast<float>(mub);

      for (int k = 0; k < kMassNz; ++k) {
        const auto pb = static_cast<double>(kC3H[static_cast<std::size_t>(k)]) * mub +
                        static_cast<double>(kC4H[static_cast<std::size_t>(k)]) +
                        static_cast<double>(kPTop);
        const auto t_init = base_temperature_k(pb) * std::pow(kP00 / pb, rd_over_cp) -
                            kThetaBase;
        const auto alb =
            (kRd / kP00) * (t_init + kThetaBase) * std::pow(pb / kP00, cvpm);
        fields.pb[index3(k, j, i)] = static_cast<float>(pb);
        fields.t_init[index3(k, j, i)] = static_cast<float>(t_init);
        fields.alb[index3(k, j, i)] = static_cast<float>(alb);
      }

      fields.phb[index3(0, j, i)] = static_cast<float>(hgt * kGravity);
      for (int k = 0; k < kMassNz; ++k) {
        const auto pfu =
            static_cast<double>(kC3F[static_cast<std::size_t>(k + 1)]) * mub +
            static_cast<double>(kC4F[static_cast<std::size_t>(k + 1)]) +
            static_cast<double>(kPTop);
        const auto pfd =
            static_cast<double>(kC3F[static_cast<std::size_t>(k)]) * mub +
            static_cast<double>(kC4F[static_cast<std::size_t>(k)]) +
            static_cast<double>(kPTop);
        const auto phm =
            static_cast<double>(kC3H[static_cast<std::size_t>(k)]) * mub +
            static_cast<double>(kC4H[static_cast<std::size_t>(k)]) +
            static_cast<double>(kPTop);
        const auto previous = static_cast<double>(fields.phb[index3(k, j, i)]);
        const auto next =
            previous + static_cast<double>(fields.alb[index3(k, j, i)]) * phm *
                           std::log(pfd / pfu);
        fields.phb[index3(k + 1, j, i)] = static_cast<float>(next);
      }
    }
  }
  return fields;
}

void write_time_vector(
    const int file_id,
    const int variable_id,
    const float* values,
    const int count) {
  const std::size_t start[2] = {0, 0};
  const std::size_t counts[2] = {1, static_cast<std::size_t>(count)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, counts, values), "write vector");
}

void write_time_2d(
    const int file_id,
    const int variable_id,
    const std::vector<float>& values) {
  const std::size_t start[3] = {0, 0, 0};
  const std::size_t counts[3] = {
      1,
      static_cast<std::size_t>(kMassNy),
      static_cast<std::size_t>(kMassNx)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, counts, values.data()), "write 2D");
}

void write_time_3d(
    const int file_id,
    const int variable_id,
    const std::vector<float>& values,
    const int nz) {
  const std::size_t start[4] = {0, 0, 0, 0};
  const std::size_t counts[4] = {
      1,
      static_cast<std::size_t>(nz),
      static_cast<std::size_t>(kMassNy),
      static_cast<std::size_t>(kMassNx)};
  check_nc(nc_put_vara_float(file_id, variable_id, start, counts, values.data()), "write 3D");
}

void create_synthetic_file(
    const std::filesystem::path& path,
    const SyntheticOptions& options) {
  const auto fields = make_synthetic_fields();
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id), "create synthetic file");

  int time_dim = -1;
  int bottom_top_dim = -1;
  int bottom_top_stag_dim = -1;
  int south_north_dim = -1;
  int west_east_dim = -1;
  check_nc(nc_def_dim(file_id, "Time", NC_UNLIMITED, &time_dim), "define Time");
  check_nc(nc_def_dim(file_id, "bottom_top", kMassNz, &bottom_top_dim), "define bottom_top");
  check_nc(
      nc_def_dim(file_id, "bottom_top_stag", kFullNz, &bottom_top_stag_dim),
      "define bottom_top_stag");
  check_nc(nc_def_dim(file_id, "south_north", kMassNy, &south_north_dim), "define south_north");
  check_nc(nc_def_dim(file_id, "west_east", kMassNx, &west_east_dim), "define west_east");

  int p_top_var = -1;
  const int p_top_dims[1] = {time_dim};
  check_nc(nc_def_var(file_id, "P_TOP", NC_FLOAT, 1, p_top_dims, &p_top_var), "define P_TOP");

  int c3f_var = -1;
  int c4f_var = -1;
  int c3h_var = -1;
  int c4h_var = -1;
  const int full_dims[2] = {time_dim, bottom_top_stag_dim};
  const int mass_coeff_dims[2] = {time_dim, bottom_top_dim};
  check_nc(nc_def_var(file_id, "C3F", NC_FLOAT, 2, full_dims, &c3f_var), "define C3F");
  check_nc(nc_def_var(file_id, "C4F", NC_FLOAT, 2, full_dims, &c4f_var), "define C4F");
  check_nc(nc_def_var(file_id, "C3H", NC_FLOAT, 2, mass_coeff_dims, &c3h_var), "define C3H");
  check_nc(nc_def_var(file_id, "C4H", NC_FLOAT, 2, mass_coeff_dims, &c4h_var), "define C4H");

  int hgt_var = -1;
  const int surface_dims[3] = {time_dim, south_north_dim, west_east_dim};
  check_nc(nc_def_var(file_id, "HGT", NC_FLOAT, 3, surface_dims, &hgt_var), "define HGT");

  int pb_var = -1;
  int t_init_var = -1;
  int alb_var = -1;
  int mub_var = -1;
  int phb_var = -1;
  const int mass_dims[4] = {time_dim, bottom_top_dim, south_north_dim, west_east_dim};
  const int full_field_dims[4] = {
      time_dim,
      bottom_top_stag_dim,
      south_north_dim,
      west_east_dim};
  if (options.include_pb) {
    check_nc(nc_def_var(file_id, "PB", NC_FLOAT, 4, mass_dims, &pb_var), "define PB");
  }
  if (options.include_t_init) {
    check_nc(nc_def_var(file_id, "T_INIT", NC_FLOAT, 4, mass_dims, &t_init_var), "define T_INIT");
  }
  if (options.include_alb) {
    check_nc(nc_def_var(file_id, "ALB", NC_FLOAT, 4, mass_dims, &alb_var), "define ALB");
  }
  if (options.include_mub) {
    check_nc(nc_def_var(file_id, "MUB", NC_FLOAT, 3, surface_dims, &mub_var), "define MUB");
  }
  if (options.include_phb) {
    check_nc(nc_def_var(file_id, "PHB", NC_FLOAT, 4, full_field_dims, &phb_var), "define PHB");
  }

  check_nc(nc_enddef(file_id), "end definitions");

  const std::size_t p_top_start[1] = {0};
  const std::size_t p_top_count[1] = {1};
  check_nc(nc_put_vara_float(file_id, p_top_var, p_top_start, p_top_count, &kPTop), "write P_TOP");
  write_time_vector(file_id, c3f_var, kC3F.data(), kFullNz);
  write_time_vector(file_id, c4f_var, kC4F.data(), kFullNz);
  write_time_vector(file_id, c3h_var, kC3H.data(), kMassNz);
  write_time_vector(file_id, c4h_var, kC4H.data(), kMassNz);
  write_time_2d(file_id, hgt_var, fields.hgt);
  if (options.include_pb) {
    write_time_3d(file_id, pb_var, fields.pb, kMassNz);
  }
  if (options.include_t_init) {
    write_time_3d(file_id, t_init_var, fields.t_init, kMassNz);
  }
  if (options.include_alb) {
    write_time_3d(file_id, alb_var, fields.alb, kMassNz);
  }
  if (options.include_mub) {
    write_time_2d(file_id, mub_var, fields.mub);
  }
  if (options.include_phb) {
    write_time_3d(file_id, phb_var, fields.phb, kFullNz);
  }

  check_nc(nc_close(file_id), "close synthetic file");
}

[[nodiscard]] std::string shell_quote(const std::filesystem::path& path) {
  return "'" + path.string() + "'";
}

[[nodiscard]] CommandResult run_probe(
    const std::filesystem::path& executable,
    const std::filesystem::path& input,
    const std::string_view domain,
    const int mass_nx,
    const int mass_ny,
    const int mass_nz,
    const int full_nz) {
  const std::string command =
      shell_quote(executable) + " --input " + shell_quote(input) + " --domain " +
      std::string(domain) + " --mass-nx " + std::to_string(mass_nx) +
      " --mass-ny " + std::to_string(mass_ny) + " --mass-nz " +
      std::to_string(mass_nz) + " --full-nz " + std::to_string(full_nz) +
      " --time-index 0 --pretty 2>&1";

  FILE* pipe = popen(command.c_str(), "r");
  if (pipe == nullptr) {
    throw std::runtime_error("failed to run base-state probe");
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

void run_complete_synthetic_test(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto path = root / "complete_base_state_inputs.nc";
  create_synthetic_file(path, {});

  const auto result =
      run_probe(executable, path, "d01", kMassNx, kMassNy, kMassNz, kFullNz);
  assert(result.exit_code == 0);
  require_contains(result.output, "\"status\": \"ok\"");
  require_contains(result.output, "\"domain\": \"d01\"");
  require_contains(result.output, "\"diagnostic_only\": true");
  require_contains(result.output, "\"gate_candidate\": false");
  require_contains(result.output, "\"integrator_output\": false");
  require_contains(result.output, "\"writes_netcdf\": false");
  require_contains(result.output, "\"modifies_state\": false");
  require_contains(result.output, "\"calls_pressure_refresh_compute\": false");
  require_contains(result.output, "\"uses_later_restart_truth\": false");
  require_contains(result.output, "\"wrote_phb\": true");
  require_contains(result.output, "\"PB\": {");
  require_contains(result.output, "\"T_INIT\": {");
  require_contains(result.output, "\"ALB\": {");
  require_contains(result.output, "\"MUB\": {");
  require_contains(result.output, "\"PHB\": {");
  require_contains(result.output, "\"status\": \"compared\"");
}

void run_missing_optional_direct_test(
    const std::filesystem::path& executable,
    const std::filesystem::path& root) {
  const auto path = root / "d02_like_missing_t_init_alb.nc";
  create_synthetic_file(
      path,
      {.include_pb = true,
       .include_t_init = false,
       .include_alb = false,
       .include_mub = true,
       .include_phb = true});

  const auto result =
      run_probe(executable, path, "d02", kMassNx, kMassNy, kMassNz, kFullNz);
  assert(result.exit_code == 0);
  require_contains(result.output, "\"status\": \"ok\"");
  require_contains(result.output, "\"domain\": \"d02\"");
  require_contains(result.output, "\"PB\": {");
  require_contains(result.output, "\"MUB\": {");
  require_contains(result.output, "\"PHB\": {");
  require_contains(result.output, "\"detail\": \"missing\"");
  require_contains(result.output, "\"status\": \"not_available\"");
}

void run_real_krosa_smoke_if_available(const std::filesystem::path& executable) {
  if (std::filesystem::exists(kKrosaWrfinputD01)) {
    const auto result = run_probe(executable, kKrosaWrfinputD01, "d01", 265, 429, 59, 60);
    assert(result.exit_code == 0);
    require_contains(result.output, "\"status\": \"ok\"");
    require_contains(result.output, "\"domain\": \"d01\"");
    require_contains(result.output, "\"PB\": {");
    require_contains(result.output, "\"T_INIT\": {");
    require_contains(result.output, "\"ALB\": {");
    require_contains(result.output, "\"MUB\": {");
    require_contains(result.output, "\"PHB\": {");
    std::cerr << "real KROSA d01 wrfinput base-state probe smoke read "
              << kKrosaWrfinputD01 << '\n';
  } else {
    std::cerr << "skipping real KROSA d01 wrfinput smoke; file is not present\n";
  }

  if (std::filesystem::exists(kKrosaWrfoutD02)) {
    const auto result = run_probe(executable, kKrosaWrfoutD02, "d02", 210, 210, 59, 60);
    assert(result.exit_code == 0);
    require_contains(result.output, "\"status\": \"ok\"");
    require_contains(result.output, "\"domain\": \"d02\"");
    require_contains(result.output, "\"PB\": {");
    require_contains(result.output, "\"MUB\": {");
    require_contains(result.output, "\"PHB\": {");
    require_contains(result.output, "\"T_INIT\": {");
    require_contains(result.output, "\"ALB\": {");
    require_contains(result.output, "\"status\": \"not_available\"");
    std::cerr << "real KROSA d02 start wrfout base-state probe smoke read "
              << kKrosaWrfoutD02 << '\n';
  } else {
    std::cerr << "skipping real KROSA d02 wrfout smoke; file is not present\n";
  }
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: base_state_probe_tool_test <tywrf_base_state_probe path>\n";
      return 2;
    }

    const auto executable = std::filesystem::path(argv[1]);
    const auto root = std::filesystem::temp_directory_path() / "tywrf_base_state_probe_test";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);

    run_complete_synthetic_test(executable, root);
    run_missing_optional_direct_test(executable, root);
    run_real_krosa_smoke_if_available(executable);

    std::filesystem::remove_all(root);
  } catch (const std::exception& error) {
    std::cerr << "base-state probe tool test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
