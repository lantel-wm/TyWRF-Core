#include "tywrf/dynamics/base_state_provider.hpp"
#include "tywrf/dynamics/pressure_refresh_hook.hpp"
#include "tywrf/io/pressure_refresh_io.hpp"
#include "tywrf/io/wrf_state_io.hpp"
#include "tywrf/nest/parent_child_interpolation.hpp"
#include "tywrf/nest/static_fields.hpp"
#include "tywrf/nest/state_exchange.hpp"
#include "tywrf/nest/state_remap.hpp"
#include "tywrf/state.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

constexpr double kD02TargetDxMeters = 2000.0;
constexpr double kD02ResolutionToleranceMeters = 0.5;

const std::vector<std::string>& template_variable_names() {
  static const std::vector<std::string> names = {"Times", "XLAT", "XLONG", "HGT"};
  return names;
}

const std::vector<std::string>& strict_state_variable_names() {
  static const std::vector<std::string> names = {"U", "V", "T", "PH", "MU", "P", "QVAPOR"};
  return names;
}

const std::vector<std::string>& optional_preserved_state_variable_names() {
  static const std::vector<std::string> names = {
      "PB", "PHB", "MUB", "PSFC", "U10", "V10", "T2", "Q2", "RAINC", "RAINNC"};
  return names;
}

const std::vector<std::string>& parent_interpolation_variable_names() {
  static const std::vector<std::string> names = {"U", "V", "T", "PH", "MU", "QVAPOR"};
  return names;
}

struct Options {
  std::filesystem::path d01_start_state_path;
  std::filesystem::path d02_start_state_path;
  std::filesystem::path template_path;
  std::filesystem::path output_path;
  std::string cycle_start;
  std::string cycle_end;
  std::string times_value;
  std::size_t d01_time_index = 0;
  std::size_t d02_time_index = 0;
  std::size_t template_time_index = 0;
  std::size_t output_time_index = 0;
  tywrf::nest::ParentChildPosition from_parent_start;
  tywrf::nest::ParentChildPosition to_parent_start;
  bool has_from_parent_start = false;
  bool has_to_parent_start = false;
  bool pressure_refresh = false;
  bool experimental_pressure_refresh_apply = false;
  bool pretty = false;
  std::vector<std::string> variables;
};

struct Resolution {
  double dx = 0.0;
  double dy = 0.0;
};

struct CandidateReport {
  std::uint64_t changed_selected_points = 0;
  std::uint64_t changed_static_template_points = 0;
  tywrf::nest::ChildStateRemapReport remap;
  tywrf::nest::RemapPlan remap_plan;
  tywrf::nest::StateExchangePlan exchange;
  tywrf::nest::ParentChildInterpolationReport interpolation;
  tywrf::nest::MovingNestStaticRefreshReport static_refresh;
  std::optional<tywrf::dynamics::KrosaBaseStateProviderReport>
      pressure_refresh_provider_probe;
  std::optional<tywrf::dynamics::KrosaPressureRefreshHookReport>
      pressure_refresh_dry_run_contract;
  std::optional<tywrf::dynamics::KrosaPressureRefreshHookReport> pressure_refresh;
  std::filesystem::path pressure_refresh_metadata_source;
  std::size_t pressure_refresh_metadata_time_index = 0;
  std::uint64_t pressure_refresh_changed_p_points = 0;
  std::uint64_t pressure_refresh_changed_pb_points = 0;
  std::uint64_t pressure_refresh_changed_mub_points = 0;
  std::uint64_t pressure_refresh_changed_phb_points = 0;
  double cen_lat = 0.0;
  double cen_lon = 0.0;
};

struct PressureRefreshReportParity {
  bool changed_p_matches_refreshed_point_count = false;
  bool invalid_and_skipped_points_zero = false;
  bool overlap_halo_untouched = false;
};

struct StaticFieldSet {
  explicit StaticFieldSet(const tywrf::Grid& grid)
      : xlat(grid.surface_layout()), xlong(grid.surface_layout()), hgt(grid.surface_layout()) {}

  tywrf::FieldStorage2D<float> xlat;
  tywrf::FieldStorage2D<float> xlong;
  tywrf::FieldStorage2D<float> hgt;
};

[[nodiscard]] const char* bool_text(const bool value) noexcept {
  return value ? "true" : "false";
}

[[nodiscard]] tywrf::dynamics::KrosaBaseStateProviderTerrainOverride
make_moved_candidate_hgt_terrain_override(const StaticFieldSet& output_static) {
  return {
      .terrain_height_m = output_static.hgt.view(),
      .source_name = "moved_candidate_HGT",
      .provenance = "override:moved_candidate_HGT"};
}

[[nodiscard]] PressureRefreshReportParity pressure_refresh_report_parity(
    const CandidateReport& candidate_report,
    const tywrf::dynamics::KrosaPressureRefreshHookReport& pressure_report) {
  return {
      .changed_p_matches_refreshed_point_count =
          candidate_report.pressure_refresh_changed_p_points ==
          pressure_report.compute_report.refreshed_point_count,
      .invalid_and_skipped_points_zero =
          pressure_report.compute_report.invalid_point_count == 0 &&
          pressure_report.compute_report.skipped_point_count == 0,
      .overlap_halo_untouched =
          !pressure_report.compute_report.touched_overlap_cells &&
          !pressure_report.compute_report.touched_halo_cells};
}

class NetcdfHandle {
 public:
  enum class Mode {
    read,
    write,
  };

  NetcdfHandle(const std::filesystem::path& path, const Mode mode) : path_(path) {
    const int status = mode == Mode::read
                           ? nc_open(path.string().c_str(), NC_NOWRITE, &id_)
                           : nc_open(path.string().c_str(), NC_WRITE, &id_);
    check(status, "open");
  }

  NetcdfHandle(const NetcdfHandle&) = delete;
  NetcdfHandle& operator=(const NetcdfHandle&) = delete;

  ~NetcdfHandle() {
    if (id_ >= 0) {
      nc_close(id_);
    }
  }

  [[nodiscard]] int id() const noexcept {
    return id_;
  }

  void check(const int status, const std::string_view operation) const {
    if (status == NC_NOERR) {
      return;
    }
    std::ostringstream message;
    message << "NetCDF " << operation << " failed for " << path_ << ": "
            << nc_strerror(status);
    throw std::runtime_error(message.str());
  }

 private:
  std::filesystem::path path_;
  int id_ = -1;
};

[[nodiscard]] std::string usage() {
  return R"(Usage:
  tywrf_selected_field_cycle --d01-start-state <wrfout_d01 start> --d02-start-state <wrfout_d02 start> --template <d02 template> --output <candidate wrfout> --cycle-start TEXT --cycle-end TEXT --from-parent-start I,J --to-parent-start I,J --times TEXT [options]

Options:
  --d01-time-index N             Time index read from --d01-start-state; default 0.
  --d02-time-index N             Time index read from --d02-start-state; default 0.
  --template-time-index N        Time index copied from --template for static fields; default 0.
  --output-time-index N          Time index written in --output; default 0.
  --variables A,B,C              Output variables; default strict fields, Times/XLAT/XLONG/HGT,
                                  plus available d02 PB/PHB/MUB/PSFC/U10/V10/T2/Q2/RAINC/RAINNC.
  --pressure-refresh             Opt in to provider-backed KROSA pressure refresh readiness check.
                                  Current selected-field state aborts before output because
                                  PB/PHB/MUB/P ownership is not yet thermodynamically consistent
                                  with the moved-pose T/PH and terrain producer.
  --pretty                       Pretty-print JSON report.
  --help                         Show this help.

Oracle inputs are not accepted: --end-state, --reference-end, --d01-end-state,
and --d02-end-state are rejected. This tool remaps the d02 start state from the
old moving-nest pose to the new pose, fills newly exposed U/V/MU/QVAPOR cells
from d01 start-state parent interpolation, fills newly exposed T/PH from the
same d01 start-state parent interpolation, and preserves P plus available d02
start-state diagnostic/base-state fields without using reference-end truth. It
is a selected-field integrator candidate, not a WRF-exact physics result.
)";
}

[[nodiscard]] std::size_t parse_size(const std::string& value, const std::string_view option) {
  std::size_t parsed_chars = 0;
  const auto result = std::stoull(value, &parsed_chars);
  if (parsed_chars != value.size()) {
    throw std::invalid_argument(std::string(option) + " expects an unsigned integer");
  }
  return static_cast<std::size_t>(result);
}

[[nodiscard]] std::vector<std::string> split_variables(const std::string& value) {
  std::vector<std::string> variables;
  std::size_t begin = 0;
  while (begin <= value.size()) {
    const auto end = value.find(',', begin);
    auto token = value.substr(begin, end == std::string::npos ? std::string::npos : end - begin);
    token.erase(token.begin(), std::find_if(token.begin(), token.end(), [](const unsigned char c) {
                  return !std::isspace(c);
                }));
    token.erase(
        std::find_if(token.rbegin(), token.rend(), [](const unsigned char c) {
          return !std::isspace(c);
        }).base(),
        token.end());
    if (!token.empty()) {
      variables.push_back(token);
    }
    if (end == std::string::npos) {
      break;
    }
    begin = end + 1;
  }
  return variables;
}

[[nodiscard]] tywrf::nest::ParentChildPosition parse_parent_start(
    const std::string& value,
    const std::string_view option) {
  const auto separator = value.find(',');
  if (separator == std::string::npos || value.find(',', separator + 1) != std::string::npos) {
    throw std::invalid_argument(std::string(option) + " must be I,J using one-based WRF indices");
  }
  const auto parse_component = [&](const std::string& component) {
    std::size_t parsed_chars = 0;
    const auto parsed = std::stol(component, &parsed_chars);
    if (parsed_chars != component.size() || parsed < 1 ||
        parsed > static_cast<long>(std::numeric_limits<std::int32_t>::max())) {
      throw std::invalid_argument(
          std::string(option) + " must contain positive one-based WRF indices");
    }
    return static_cast<std::int32_t>(parsed);
  };
  return {
      parse_component(value.substr(0, separator)),
      parse_component(value.substr(separator + 1)),
      tywrf::nest::IndexBase::one_based};
}

[[nodiscard]] bool is_oracle_input_option(const std::string& arg) {
  return arg == "--end-state" || arg == "--reference-end" || arg == "--d01-end-state" ||
         arg == "--d02-end-state" || arg.rfind("--end-state=", 0) == 0 ||
         arg.rfind("--reference-end=", 0) == 0 || arg.rfind("--d01-end-state=", 0) == 0 ||
         arg.rfind("--d02-end-state=", 0) == 0;
}

[[nodiscard]] Options parse_options(const int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string arg = argv[index];
    if (is_oracle_input_option(arg)) {
      throw std::invalid_argument("selected-field cycle rejects end-state/reference-end oracle input");
    }

    const auto require_value = [&](const std::string_view name) -> std::string {
      if (index + 1 >= argc) {
        throw std::invalid_argument(std::string(name) + " requires a value");
      }
      ++index;
      return argv[index];
    };

    if (arg == "--help" || arg == "-h") {
      std::cout << usage();
      std::exit(0);
    } else if (arg == "--d01-start-state") {
      options.d01_start_state_path = require_value(arg);
    } else if (arg == "--d02-start-state") {
      options.d02_start_state_path = require_value(arg);
    } else if (arg == "--template") {
      options.template_path = require_value(arg);
    } else if (arg == "--output") {
      options.output_path = require_value(arg);
    } else if (arg == "--cycle-start") {
      options.cycle_start = require_value(arg);
    } else if (arg == "--cycle-end") {
      options.cycle_end = require_value(arg);
    } else if (arg == "--from-parent-start") {
      options.from_parent_start = parse_parent_start(require_value(arg), arg);
      options.has_from_parent_start = true;
    } else if (arg == "--to-parent-start") {
      options.to_parent_start = parse_parent_start(require_value(arg), arg);
      options.has_to_parent_start = true;
    } else if (arg == "--times") {
      options.times_value = require_value(arg);
    } else if (arg == "--d01-time-index") {
      options.d01_time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--d02-time-index") {
      options.d02_time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--template-time-index") {
      options.template_time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--output-time-index") {
      options.output_time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--variables") {
      options.variables = split_variables(require_value(arg));
    } else if (arg == "--pressure-refresh") {
      options.pressure_refresh = true;
    } else if (arg == "--experimental-pressure-refresh-apply") {
      options.experimental_pressure_refresh_apply = true;
    } else if (arg == "--pretty") {
      options.pretty = true;
    } else {
      throw std::invalid_argument("unknown option: " + arg);
    }
  }

  if (options.d01_start_state_path.empty()) {
    throw std::invalid_argument("--d01-start-state is required");
  }
  if (options.d02_start_state_path.empty()) {
    throw std::invalid_argument("--d02-start-state is required");
  }
  if (options.template_path.empty()) {
    throw std::invalid_argument("--template is required");
  }
  if (options.output_path.empty()) {
    throw std::invalid_argument("--output is required");
  }
  if (options.cycle_start.empty()) {
    throw std::invalid_argument("--cycle-start is required");
  }
  if (options.cycle_end.empty()) {
    throw std::invalid_argument("--cycle-end is required");
  }
  if (options.times_value.empty()) {
    throw std::invalid_argument("--times is required");
  }
  if (!options.has_from_parent_start) {
    throw std::invalid_argument("--from-parent-start is required");
  }
  if (!options.has_to_parent_start) {
    throw std::invalid_argument("--to-parent-start is required");
  }
  if (options.experimental_pressure_refresh_apply && !options.pressure_refresh) {
    throw std::invalid_argument(
        "--experimental-pressure-refresh-apply requires --pressure-refresh");
  }
  return options;
}

[[nodiscard]] bool contains(
    const std::vector<std::string>& values,
    const std::string_view value) {
  return std::find(values.begin(), values.end(), value) != values.end();
}

void require_output_variables(const std::vector<std::string>& variables) {
  for (const auto& name : template_variable_names()) {
    if (!contains(variables, name)) {
      throw std::invalid_argument(
          "selected-field candidate output must include time/static field " + name);
    }
  }
  for (const auto& name : strict_state_variable_names()) {
    if (!contains(variables, name)) {
      throw std::invalid_argument("selected-field candidate output must include strict field " + name);
    }
  }
}

[[nodiscard]] bool has_variable(const NetcdfHandle& file, const std::string_view name) {
  int variable_id = -1;
  const int status = nc_inq_varid(file.id(), std::string(name).c_str(), &variable_id);
  if (status == NC_NOERR) {
    return true;
  }
  if (status == NC_ENOTVAR) {
    return false;
  }
  file.check(status, "inquire variable");
  return false;
}

[[nodiscard]] std::vector<std::string> available_optional_preserved_variables(
    const std::filesystem::path& path) {
  const NetcdfHandle file(path, NetcdfHandle::Mode::read);
  std::vector<std::string> variables;
  for (const auto& name : optional_preserved_state_variable_names()) {
    if (has_variable(file, name)) {
      variables.push_back(name);
    }
  }
  return variables;
}

[[nodiscard]] std::vector<std::string> default_output_variables(
    const std::filesystem::path& d02_start_state_path) {
  auto variables = template_variable_names();
  const auto& strict = strict_state_variable_names();
  variables.insert(variables.end(), strict.begin(), strict.end());
  const auto optional = available_optional_preserved_variables(d02_start_state_path);
  variables.insert(variables.end(), optional.begin(), optional.end());
  return variables;
}

[[nodiscard]] std::vector<std::string> selected_d02_read_variables(
    const std::filesystem::path& d02_start_state_path) {
  auto variables = strict_state_variable_names();
  const auto optional = available_optional_preserved_variables(d02_start_state_path);
  variables.insert(variables.end(), optional.begin(), optional.end());
  return variables;
}

[[nodiscard]] std::optional<double> read_double_attr(
    const NetcdfHandle& file,
    const std::string_view name) {
  double value = 0.0;
  const int status = nc_get_att_double(file.id(), NC_GLOBAL, std::string(name).c_str(), &value);
  if (status == NC_ENOTATT) {
    return std::nullopt;
  }
  file.check(status, "read global attribute");
  return value;
}

[[nodiscard]] Resolution read_resolution(const std::filesystem::path& path) {
  const NetcdfHandle file(path, NetcdfHandle::Mode::read);
  const auto dx = read_double_attr(file, "DX");
  const auto dy = read_double_attr(file, "DY");
  if (!dx.has_value() || !dy.has_value()) {
    throw std::runtime_error("WRF input must carry DX and DY attributes: " + path.string());
  }
  return {*dx, *dy};
}

void require_d02_resolution(const std::filesystem::path& path, const Resolution resolution) {
  const auto near_target = [](const double value) {
    return std::abs(value - kD02TargetDxMeters) <= kD02ResolutionToleranceMeters;
  };
  if (near_target(resolution.dx) && near_target(resolution.dy)) {
    return;
  }
  std::ostringstream message;
  message << "d02 resolution must remain 2 km; " << path << " has DX=" << resolution.dx
          << " DY=" << resolution.dy;
  throw std::runtime_error(message.str());
}

[[nodiscard]] std::int32_t rounded_spacing_m(
    const Resolution resolution,
    const std::string_view label) {
  if (std::abs(resolution.dx - resolution.dy) > kD02ResolutionToleranceMeters) {
    throw std::runtime_error(std::string(label) + " DX/DY must match");
  }
  const auto rounded = std::lround(resolution.dx);
  if (rounded <= 0 || rounded > std::numeric_limits<std::int32_t>::max()) {
    throw std::runtime_error(std::string(label) + " has unsupported grid spacing");
  }
  return static_cast<std::int32_t>(rounded);
}

[[nodiscard]] tywrf::nest::HorizontalDomainDescriptor make_domain_descriptor(
    const int domain_id,
    const Resolution resolution,
    const tywrf::Grid& grid,
    const std::string_view label) {
  const auto& config = grid.config();
  return {
      domain_id,
      rounded_spacing_m(resolution, label),
      config.mass_nx,
      config.mass_ny,
      config.mass_nx + 1,
      config.mass_ny + 1};
}

[[nodiscard]] tywrf::nest::ParentChildDescriptor make_descriptor(
    const Resolution d01_resolution,
    const Resolution d02_resolution,
    const tywrf::Grid& d01_grid,
    const tywrf::Grid& d02_grid) {
  const auto parent_spacing = rounded_spacing_m(d01_resolution, "d01 start-state");
  const auto child_spacing = rounded_spacing_m(d02_resolution, "d02 start-state");
  if (parent_spacing % child_spacing != 0) {
    throw std::runtime_error("parent/child grid spacing must be an integer ratio");
  }
  const auto ratio = parent_spacing / child_spacing;
  tywrf::nest::ParentChildDescriptor descriptor{
      make_domain_descriptor(1, d01_resolution, d01_grid, "d01 start-state"),
      make_domain_descriptor(2, d02_resolution, d02_grid, "d02 start-state"),
      ratio,
      ratio,
  };
  if (const auto status = tywrf::nest::validate_parent_child_descriptor(descriptor);
      !status.ok()) {
    throw std::runtime_error(
        "invalid parent-child descriptor for selected-field cycle: " +
        std::string(status.message));
  }
  if (ratio != 5) {
    throw std::runtime_error("selected-field cycle currently supports KROSA parent_grid_ratio=5");
  }
  return descriptor;
}

void require_path_exists(const std::filesystem::path& path, const std::string_view label) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error(std::string(label) + " does not exist: " + path.string());
  }
}

void load_static_2d_variable(
    const std::filesystem::path& path,
    const std::string_view name,
    tywrf::FieldStorage2D<float>& field,
    const std::size_t time_index) {
  const NetcdfHandle file(path, NetcdfHandle::Mode::read);
  int variable_id = -1;
  file.check(nc_inq_varid(file.id(), std::string(name).c_str(), &variable_id), "inquire static variable");
  const auto layout = field.layout();
  const std::array<std::size_t, 3> start = {time_index, 0, 0};
  const std::array<std::size_t, 3> count = {
      1,
      static_cast<std::size_t>(layout.active_ny()),
      static_cast<std::size_t>(layout.active_nx())};
  std::vector<float> buffer(count[1] * count[2]);
  file.check(
      nc_get_vara_float(file.id(), variable_id, start.data(), count.data(), buffer.data()),
      "read static variable");

  auto view = field.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
      view(layout.i_begin() + i, layout.j_begin() + j) =
          buffer[static_cast<std::size_t>(j) * count[2] + static_cast<std::size_t>(i)];
    }
  }
}

[[nodiscard]] StaticFieldSet load_static_fields(
    const std::filesystem::path& path,
    const tywrf::Grid& grid,
    const std::size_t time_index) {
  StaticFieldSet fields(grid);
  load_static_2d_variable(path, "XLAT", fields.xlat, time_index);
  load_static_2d_variable(path, "XLONG", fields.xlong, time_index);
  load_static_2d_variable(path, "HGT", fields.hgt, time_index);
  return fields;
}

[[nodiscard]] tywrf::FieldStorage2D<float> load_hgt_field(
    const std::filesystem::path& path,
    const tywrf::Grid& grid,
    const std::size_t time_index) {
  tywrf::FieldStorage2D<float> hgt(grid.surface_layout());
  load_static_2d_variable(path, "HGT", hgt, time_index);
  return hgt;
}

void write_static_2d_variable(
    const NetcdfHandle& file,
    const std::string_view name,
    const tywrf::FieldStorage2D<float>& field,
    const std::size_t time_index) {
  int variable_id = -1;
  file.check(nc_inq_varid(file.id(), std::string(name).c_str(), &variable_id), "inquire output static variable");
  const auto layout = field.layout();
  const std::array<std::size_t, 3> start = {time_index, 0, 0};
  const std::array<std::size_t, 3> count = {
      1,
      static_cast<std::size_t>(layout.active_ny()),
      static_cast<std::size_t>(layout.active_nx())};
  std::vector<float> buffer(count[1] * count[2]);
  const auto view = field.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
      buffer[static_cast<std::size_t>(j) * count[2] + static_cast<std::size_t>(i)] =
          view(layout.i_begin() + i, layout.j_begin() + j);
    }
  }
  file.check(
      nc_put_vara_float(file.id(), variable_id, start.data(), count.data(), buffer.data()),
      "write output static variable");
}

void overwrite_output_static_fields(
    const std::filesystem::path& output_path,
    const StaticFieldSet& fields,
    const std::size_t time_index) {
  const NetcdfHandle file(output_path, NetcdfHandle::Mode::write);
  write_static_2d_variable(file, "XLAT", fields.xlat, time_index);
  write_static_2d_variable(file, "XLONG", fields.xlong, time_index);
  write_static_2d_variable(file, "HGT", fields.hgt, time_index);
}

template <typename BeforeStorage, typename AfterStorage>
[[nodiscard]] std::uint64_t changed_points(
    const BeforeStorage& before,
    const AfterStorage& after) {
  if (before.size() != after.size()) {
    throw std::runtime_error("storage size mismatch while checking selected-field changes");
  }
  std::uint64_t changed = 0;
  for (std::size_t index = 0; index < before.size(); ++index) {
    changed += before.data()[index] != after.data()[index] ? 1U : 0U;
  }
  return changed;
}

[[nodiscard]] std::uint64_t changed_selected_points(
    const tywrf::State<float>& before,
    const tywrf::State<float>& after) {
  return changed_points(before.u, after.u) + changed_points(before.v, after.v) +
         changed_points(before.t, after.t) + changed_points(before.ph, after.ph) +
         changed_points(before.mu, after.mu) + changed_points(before.qvapor, after.qvapor);
}

template <typename Storage>
void require_finite_storage(const std::string_view name, const Storage& storage) {
  for (std::size_t index = 0; index < storage.size(); ++index) {
    if (!std::isfinite(storage.data()[index])) {
      throw std::runtime_error("selected-field candidate contains non-finite " + std::string(name));
    }
  }
}

void require_finite_strict_fields(const tywrf::State<float>& state) {
  require_finite_storage("U", state.u);
  require_finite_storage("V", state.v);
  require_finite_storage("T", state.t);
  require_finite_storage("PH", state.ph);
  require_finite_storage("MU", state.mu);
  require_finite_storage("P", state.p);
  require_finite_storage("QVAPOR", state.qvapor);
}

void require_finite_static_fields(const StaticFieldSet& fields) {
  require_finite_storage("XLAT", fields.xlat);
  require_finite_storage("XLONG", fields.xlong);
  require_finite_storage("HGT", fields.hgt);
}

[[nodiscard]] std::uint64_t changed_static_points(
    const StaticFieldSet& before,
    const StaticFieldSet& after) {
  return changed_points(before.xlat, after.xlat) + changed_points(before.xlong, after.xlong) +
         changed_points(before.hgt, after.hgt);
}

[[nodiscard]] std::string join_variables(const std::vector<std::string>& variables);

struct PressureRefreshReadiness {
  bool static_refresh_applied = false;
  bool static_refresh_uses_reference_end = false;
  bool thermodynamic_base_state_consistency_ready = false;
  bool provider_terrain_uses_moved_candidate_hgt = false;
  bool provider_base_state_reconstruct_ok = false;
  bool base_state_sync_contract_ok = false;
  bool base_state_sync_dry_run = false;
  bool base_state_sync_applied = false;
  std::uint64_t would_sync_pb_point_count = 0;
  std::uint64_t would_sync_mub_point_count = 0;
  std::uint64_t would_sync_phb_point_count = 0;
  std::uint64_t sync_overlap_write_count = 0;
  std::uint64_t sync_halo_write_count = 0;
  bool pressure_refresh_compute_called = false;
  bool pressure_compute_dry_run = false;
  bool pressure_compute_dry_run_called = false;
  bool pressure_compute_dry_run_ok = false;
  std::uint64_t would_refresh_p_point_count = 0;
  std::uint64_t dry_run_invalid_p_point_count = 0;
  std::uint64_t dry_run_skipped_p_point_count = 0;
  std::uint64_t pressure_compute_dry_run_report_target_column_count = 0;
  std::uint64_t pressure_compute_dry_run_report_refreshed_point_count = 0;
  std::uint64_t pressure_compute_dry_run_report_invalid_point_count = 0;
  std::uint64_t pressure_compute_dry_run_report_skipped_point_count = 0;
  bool pressure_compute_dry_run_report_touched_overlap_cells = false;
  bool pressure_compute_dry_run_report_touched_halo_cells = false;
  bool pressure_refresh_applied = false;
  std::string provider_terrain_source_name;
  std::string provider_terrain_provenance;

  [[nodiscard]] bool ready() const noexcept {
    return static_refresh_applied && !static_refresh_uses_reference_end &&
           thermodynamic_base_state_consistency_ready &&
           provider_terrain_uses_moved_candidate_hgt &&
           provider_base_state_reconstruct_ok && base_state_sync_contract_ok &&
           base_state_sync_dry_run && !base_state_sync_applied &&
           would_sync_pb_point_count > 0 && would_sync_mub_point_count > 0 &&
           would_sync_phb_point_count > 0 && sync_overlap_write_count == 0 &&
           sync_halo_write_count == 0 && !pressure_refresh_compute_called &&
           pressure_compute_dry_run && pressure_compute_dry_run_called &&
           pressure_compute_dry_run_ok && would_refresh_p_point_count > 0 &&
           dry_run_invalid_p_point_count == 0 &&
           dry_run_skipped_p_point_count == 0 &&
           pressure_compute_dry_run_report_target_column_count > 0 &&
           pressure_compute_dry_run_report_refreshed_point_count ==
               would_refresh_p_point_count &&
           pressure_compute_dry_run_report_invalid_point_count == 0 &&
           pressure_compute_dry_run_report_skipped_point_count == 0 &&
           !pressure_compute_dry_run_report_touched_overlap_cells &&
           !pressure_compute_dry_run_report_touched_halo_cells &&
           !pressure_refresh_applied;
  }
};

[[nodiscard]] PressureRefreshReadiness evaluate_pressure_refresh_readiness(
    const CandidateReport& report) {
  PressureRefreshReadiness readiness;
  readiness.static_refresh_applied =
      report.static_refresh.ok() && report.changed_static_template_points > 0;
  readiness.static_refresh_uses_reference_end = false;
  if (report.pressure_refresh_provider_probe.has_value()) {
    const auto& provider = *report.pressure_refresh_provider_probe;
    readiness.provider_base_state_reconstruct_ok =
        provider.ok() && provider.allocated_buffers && provider.wrote_pb &&
        provider.wrote_t_init && provider.wrote_mub && provider.wrote_alb &&
        provider.wrote_phb;
    readiness.provider_terrain_uses_moved_candidate_hgt =
        provider.ok() && provider.terrain_override_used &&
        provider.terrain_source_name == "moved_candidate_HGT" &&
        provider.terrain_provenance == "override:moved_candidate_HGT";
    readiness.provider_terrain_source_name = provider.terrain_source_name;
    readiness.provider_terrain_provenance = provider.terrain_provenance;
  }
  if (report.pressure_refresh_dry_run_contract.has_value()) {
    const auto& dry_run = *report.pressure_refresh_dry_run_contract;
    readiness.base_state_sync_contract_ok = dry_run.provider_ok &&
                                            dry_run.base_state_sync_contract_ok;
    readiness.base_state_sync_dry_run = dry_run.base_state_sync_dry_run;
    readiness.base_state_sync_applied = dry_run.base_state_sync_applied;
    readiness.would_sync_pb_point_count = dry_run.would_sync_pb_point_count;
    readiness.would_sync_mub_point_count = dry_run.would_sync_mub_point_count;
    readiness.would_sync_phb_point_count = dry_run.would_sync_phb_point_count;
    readiness.sync_overlap_write_count = dry_run.sync_overlap_write_count;
    readiness.sync_halo_write_count = dry_run.sync_halo_write_count;
    readiness.pressure_refresh_compute_called = dry_run.calls_pressure_refresh_compute;
    readiness.pressure_compute_dry_run = dry_run.pressure_compute_dry_run;
    readiness.pressure_compute_dry_run_called =
        dry_run.pressure_compute_dry_run_called;
    readiness.pressure_compute_dry_run_ok = dry_run.pressure_compute_dry_run_ok;
    readiness.would_refresh_p_point_count = dry_run.would_refresh_p_point_count;
    readiness.dry_run_invalid_p_point_count =
        dry_run.dry_run_invalid_p_point_count;
    readiness.dry_run_skipped_p_point_count =
        dry_run.pressure_compute_dry_run_report.skipped_point_count;
    readiness.pressure_compute_dry_run_report_target_column_count =
        dry_run.pressure_compute_dry_run_report.target_column_count;
    readiness.pressure_compute_dry_run_report_refreshed_point_count =
        dry_run.pressure_compute_dry_run_report.refreshed_point_count;
    readiness.pressure_compute_dry_run_report_invalid_point_count =
        dry_run.pressure_compute_dry_run_report.invalid_point_count;
    readiness.pressure_compute_dry_run_report_skipped_point_count =
        dry_run.pressure_compute_dry_run_report.skipped_point_count;
    readiness.pressure_compute_dry_run_report_touched_overlap_cells =
        dry_run.pressure_compute_dry_run_report.touched_overlap_cells;
    readiness.pressure_compute_dry_run_report_touched_halo_cells =
        dry_run.pressure_compute_dry_run_report.touched_halo_cells;
    readiness.pressure_refresh_applied = dry_run.pressure_refresh_applied;
  }
  readiness.thermodynamic_base_state_consistency_ready =
      readiness.provider_terrain_uses_moved_candidate_hgt &&
      readiness.provider_base_state_reconstruct_ok &&
      readiness.base_state_sync_contract_ok &&
      readiness.base_state_sync_dry_run && !readiness.base_state_sync_applied &&
      readiness.would_sync_pb_point_count > 0 &&
      readiness.would_sync_mub_point_count > 0 &&
      readiness.would_sync_phb_point_count > 0 &&
      readiness.sync_overlap_write_count == 0 &&
      readiness.sync_halo_write_count == 0 &&
      !readiness.pressure_refresh_compute_called &&
      readiness.pressure_compute_dry_run &&
      readiness.pressure_compute_dry_run_called &&
      readiness.pressure_compute_dry_run_ok &&
      readiness.would_refresh_p_point_count > 0 &&
      readiness.dry_run_invalid_p_point_count == 0 &&
      readiness.dry_run_skipped_p_point_count == 0 &&
      readiness.pressure_compute_dry_run_report_target_column_count > 0 &&
      readiness.pressure_compute_dry_run_report_refreshed_point_count ==
          readiness.would_refresh_p_point_count &&
      readiness.pressure_compute_dry_run_report_invalid_point_count == 0 &&
      readiness.pressure_compute_dry_run_report_skipped_point_count == 0 &&
      !readiness.pressure_compute_dry_run_report_touched_overlap_cells &&
      !readiness.pressure_compute_dry_run_report_touched_halo_cells &&
      !readiness.pressure_refresh_applied;
  return readiness;
}

[[nodiscard]] std::string pressure_refresh_not_ready_message(
    const PressureRefreshReadiness readiness) {
  std::ostringstream message;
  message << "pressure_refresh_not_ready: static refreshed but thermodynamic/base-state "
             "consistency missing";
  message << "; static_refresh_applied="
          << (readiness.static_refresh_applied ? "true" : "false");
  message << "; static_refresh_uses_reference_end="
          << (readiness.static_refresh_uses_reference_end ? "true" : "false");
  message << "; thermodynamic_base_state_consistency_ready="
          << (readiness.thermodynamic_base_state_consistency_ready ? "true" : "false");
  message << "; exposed T/PH are parent interpolated from d01 start-state fields, but "
             "PB/PHB/MUB/P base-state ownership is still preserved d02 start-state data";
  message << "; provider_terrain_uses_moved_candidate_hgt="
          << (readiness.provider_terrain_uses_moved_candidate_hgt ? "true" : "false");
  message << "; provider_base_state_reconstruct_ok="
          << (readiness.provider_base_state_reconstruct_ok ? "true" : "false");
  if (!readiness.provider_terrain_source_name.empty()) {
    message << "; provider_terrain_source=" << readiness.provider_terrain_source_name;
  }
  if (!readiness.provider_terrain_provenance.empty()) {
    message << "; provider_terrain_provenance=" << readiness.provider_terrain_provenance;
  }
  message << "; base_state_sync_contract_ok="
          << (readiness.base_state_sync_contract_ok ? "true" : "false");
  message << "; base_state_sync_dry_run="
          << (readiness.base_state_sync_dry_run ? "true" : "false");
  message << "; base_state_sync_applied="
          << (readiness.base_state_sync_applied ? "true" : "false");
  message << "; would_sync_pb_point_count=" << readiness.would_sync_pb_point_count;
  message << "; would_sync_mub_point_count=" << readiness.would_sync_mub_point_count;
  message << "; would_sync_phb_point_count=" << readiness.would_sync_phb_point_count;
  message << "; sync_overlap_write_count=" << readiness.sync_overlap_write_count;
  message << "; sync_halo_write_count=" << readiness.sync_halo_write_count;
  message << "; pressure_refresh_compute_called="
          << (readiness.pressure_refresh_compute_called ? "true" : "false");
  message << "; pressure_compute_dry_run="
          << (readiness.pressure_compute_dry_run ? "true" : "false");
  message << "; pressure_compute_dry_run_called="
          << (readiness.pressure_compute_dry_run_called ? "true" : "false");
  message << "; pressure_compute_dry_run_ok="
          << (readiness.pressure_compute_dry_run_ok ? "true" : "false");
  message << "; would_refresh_p_point_count="
          << readiness.would_refresh_p_point_count;
  message << "; dry_run_invalid_p_point_count="
          << readiness.dry_run_invalid_p_point_count;
  message << "; dry_run_skipped_p_point_count="
          << readiness.dry_run_skipped_p_point_count;
  message << "; pressure_compute_dry_run_report_target_column_count="
          << readiness.pressure_compute_dry_run_report_target_column_count;
  message << "; pressure_compute_dry_run_report_refreshed_point_count="
          << readiness.pressure_compute_dry_run_report_refreshed_point_count;
  message << "; pressure_compute_dry_run_report_invalid_point_count="
          << readiness.pressure_compute_dry_run_report_invalid_point_count;
  message << "; pressure_compute_dry_run_report_skipped_point_count="
          << readiness.pressure_compute_dry_run_report_skipped_point_count;
  message << "; pressure_compute_dry_run_report_touched_overlap_cells="
          << (readiness.pressure_compute_dry_run_report_touched_overlap_cells ? "true"
                                                                              : "false");
  message << "; pressure_compute_dry_run_report_touched_halo_cells="
          << (readiness.pressure_compute_dry_run_report_touched_halo_cells ? "true" : "false");
  message << "; pressure_refresh_applied="
          << (readiness.pressure_refresh_applied ? "true" : "false");
  return message.str();
}

void require_pressure_refresh_ready_for_compute(const CandidateReport& report) {
  const auto readiness = evaluate_pressure_refresh_readiness(report);
  if (readiness.ready()) {
    return;
  }
  throw std::runtime_error(pressure_refresh_not_ready_message(readiness));
}

void require_pressure_refresh_inputs_ready(
    const tywrf::io::KrosaPressureRefreshReadResult& inputs) {
  if (inputs.ok() || inputs.base_state_reconstruction_inputs_ready()) {
    return;
  }

  std::ostringstream message;
  message << "pressure refresh metadata source is not ready for provider-backed selected-field "
             "candidate";
  if (!inputs.report.missing_names.empty()) {
    message << "; missing direct inputs=" << join_variables(inputs.report.missing_names);
  }
  if (!inputs.report.missing_base_state_reconstruction_names.empty()) {
    message << "; missing provider inputs="
            << join_variables(inputs.report.missing_base_state_reconstruction_names);
  }
  throw std::runtime_error(message.str());
}

void probe_pressure_refresh_provider_readiness(
    const Options& options,
    const tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    CandidateReport& report) {
  tywrf::FieldStorage3D<float> direct_alb(candidate.grid.mass_layout());
  const auto inputs = tywrf::io::read_krosa_pressure_refresh_inputs(
      options.template_path,
      candidate.grid,
      direct_alb,
      {.time_index = options.template_time_index});
  require_pressure_refresh_inputs_ready(inputs);

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto terrain_override = make_moved_candidate_hgt_terrain_override(output_static);
  auto provider_report =
      provider.reconstruct(candidate.grid, inputs.metadata, terrain_override);
  const bool uses_moved_candidate_hgt =
      provider_report.terrain_override_used &&
      provider_report.terrain_source_name == "moved_candidate_HGT";
  if (!provider_report.ok() || !provider_report.allocated_buffers ||
      !provider_report.wrote_pb || !provider_report.wrote_t_init ||
      !provider_report.wrote_mub || !provider_report.wrote_alb ||
      !provider_report.wrote_phb || !uses_moved_candidate_hgt) {
    std::ostringstream message;
    message << "pressure_refresh_provider_probe_failed";
    if (provider_report.result.message != nullptr &&
        provider_report.result.message[0] != '\0') {
      message << ": " << provider_report.result.message;
    }
    message << "; provider_terrain_uses_moved_candidate_hgt="
            << (uses_moved_candidate_hgt ? "true" : "false");
    message << "; provider_terrain_source=" << provider_report.terrain_source_name;
    message << "; provider_terrain_provenance=" << provider_report.terrain_provenance;
    throw std::runtime_error(message.str());
  }

  report.pressure_refresh_provider_probe = provider_report;
  report.pressure_refresh_metadata_source = options.template_path;
  report.pressure_refresh_metadata_time_index = options.template_time_index;
}

void probe_pressure_refresh_dry_run_contract(
    const Options& options,
    tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    CandidateReport& report) {
  tywrf::FieldStorage3D<float> direct_alb(candidate.grid.mass_layout());
  const auto inputs = tywrf::io::read_krosa_pressure_refresh_inputs(
      options.template_path,
      candidate.grid,
      direct_alb,
      {.time_index = options.template_time_index});
  require_pressure_refresh_inputs_ready(inputs);

  const auto terrain_override = make_moved_candidate_hgt_terrain_override(output_static);
  tywrf::dynamics::KrosaPressureRefreshHookOptions hook_options{};
  hook_options.terrain_override = &terrain_override;
  hook_options.base_state_sync_dry_run = true;
  hook_options.pressure_compute_dry_run = true;
  auto dry_run_report = tywrf::dynamics::apply_krosa_moving_nest_pressure_refresh_hook(
      report.remap_plan,
      candidate,
      inputs.metadata,
      hook_options);
  report.pressure_refresh_dry_run_contract = dry_run_report;
  report.pressure_refresh_metadata_source = options.template_path;
  report.pressure_refresh_metadata_time_index = options.template_time_index;

  const bool dry_run_contract_ok =
      dry_run_report.ok() && dry_run_report.provider_ok &&
      dry_run_report.base_state_sync_dry_run &&
      dry_run_report.base_state_sync_contract_ok &&
      !dry_run_report.base_state_sync_applied &&
      dry_run_report.would_sync_pb_point_count > 0 &&
      dry_run_report.would_sync_mub_point_count > 0 &&
      dry_run_report.would_sync_phb_point_count > 0 &&
      dry_run_report.sync_overlap_write_count == 0 &&
      dry_run_report.sync_halo_write_count == 0 &&
      !dry_run_report.calls_pressure_refresh_compute &&
      dry_run_report.pressure_compute_dry_run &&
      dry_run_report.pressure_compute_dry_run_called &&
      dry_run_report.pressure_compute_dry_run_ok &&
      dry_run_report.would_refresh_p_point_count > 0 &&
      dry_run_report.dry_run_invalid_p_point_count == 0 &&
      dry_run_report.pressure_compute_dry_run_report.target_column_count > 0 &&
      dry_run_report.pressure_compute_dry_run_report.refreshed_point_count ==
          dry_run_report.would_refresh_p_point_count &&
      dry_run_report.pressure_compute_dry_run_report.invalid_point_count == 0 &&
      dry_run_report.pressure_compute_dry_run_report.skipped_point_count == 0 &&
      !dry_run_report.pressure_compute_dry_run_report.touched_overlap_cells &&
      !dry_run_report.pressure_compute_dry_run_report.touched_halo_cells &&
      !dry_run_report.pressure_refresh_applied;
  if (dry_run_contract_ok) {
    return;
  }

  std::ostringstream message;
  message << "pressure_refresh_dry_run_contract_failed";
  if (dry_run_report.result.message != nullptr &&
      dry_run_report.result.message[0] != '\0') {
    message << ": " << dry_run_report.result.message;
  }
  message << "; provider_ok=" << (dry_run_report.provider_ok ? "true" : "false")
          << "; base_state_sync_contract_ok="
          << (dry_run_report.base_state_sync_contract_ok ? "true" : "false")
          << "; base_state_sync_dry_run="
          << (dry_run_report.base_state_sync_dry_run ? "true" : "false")
          << "; base_state_sync_applied="
          << (dry_run_report.base_state_sync_applied ? "true" : "false")
          << "; would_sync_pb_point_count=" << dry_run_report.would_sync_pb_point_count
          << "; would_sync_mub_point_count=" << dry_run_report.would_sync_mub_point_count
          << "; would_sync_phb_point_count=" << dry_run_report.would_sync_phb_point_count
          << "; sync_overlap_write_count=" << dry_run_report.sync_overlap_write_count
          << "; sync_halo_write_count=" << dry_run_report.sync_halo_write_count
          << "; pressure_refresh_compute_called="
          << (dry_run_report.calls_pressure_refresh_compute ? "true" : "false")
          << "; pressure_compute_dry_run="
          << (dry_run_report.pressure_compute_dry_run ? "true" : "false")
          << "; pressure_compute_dry_run_called="
          << (dry_run_report.pressure_compute_dry_run_called ? "true" : "false")
          << "; pressure_compute_dry_run_ok="
          << (dry_run_report.pressure_compute_dry_run_ok ? "true" : "false")
          << "; would_refresh_p_point_count=" << dry_run_report.would_refresh_p_point_count
          << "; dry_run_invalid_p_point_count="
          << dry_run_report.dry_run_invalid_p_point_count
          << "; dry_run_skipped_p_point_count="
          << dry_run_report.pressure_compute_dry_run_report.skipped_point_count
          << "; pressure_compute_dry_run_report_target_column_count="
          << dry_run_report.pressure_compute_dry_run_report.target_column_count
          << "; pressure_compute_dry_run_report_refreshed_point_count="
          << dry_run_report.pressure_compute_dry_run_report.refreshed_point_count
          << "; pressure_compute_dry_run_report_invalid_point_count="
          << dry_run_report.pressure_compute_dry_run_report.invalid_point_count
          << "; pressure_compute_dry_run_report_skipped_point_count="
          << dry_run_report.pressure_compute_dry_run_report.skipped_point_count
          << "; pressure_compute_dry_run_report_touched_overlap_cells="
          << (dry_run_report.pressure_compute_dry_run_report.touched_overlap_cells ? "true"
                                                                                   : "false")
          << "; pressure_compute_dry_run_report_touched_halo_cells="
          << (dry_run_report.pressure_compute_dry_run_report.touched_halo_cells ? "true"
                                                                                : "false")
          << "; pressure_refresh_applied="
          << (dry_run_report.pressure_refresh_applied ? "true" : "false");
  throw std::runtime_error(message.str());
}

void require_pressure_refresh_hook_success(
    const tywrf::dynamics::KrosaPressureRefreshHookReport& report) {
  const bool uses_moved_candidate_hgt =
      report.provider_report.terrain_override_used &&
      report.provider_report.terrain_source_name == "moved_candidate_HGT" &&
      report.provider_report.terrain_provenance == "override:moved_candidate_HGT";
  if (report.ok() && report.provider_ok && report.staging_ok &&
      report.calls_pressure_refresh_compute && report.pressure_refresh_applied &&
      report.base_state_sync_applied && report.synced_pb_point_count > 0 &&
      report.synced_mub_point_count > 0 && report.synced_phb_point_count > 0 &&
      report.compute_report.refreshed_point_count > 0 &&
      report.compute_report.invalid_point_count == 0 &&
      report.compute_report.skipped_point_count == 0 &&
      !report.touched_overlap_cells && !report.touched_halo_cells &&
      uses_moved_candidate_hgt) {
    return;
  }

  std::ostringstream message;
  message << "selected-field pressure refresh hook failed";
  if (report.result.message != nullptr && report.result.message[0] != '\0') {
    message << ": " << report.result.message;
  }
  message << " provider_ok=" << (report.provider_ok ? "true" : "false")
          << " staging_ok=" << (report.staging_ok ? "true" : "false")
          << " compute_called=" << (report.calls_pressure_refresh_compute ? "true" : "false")
          << " applied=" << (report.pressure_refresh_applied ? "true" : "false")
          << " base_state_sync_applied="
          << (report.base_state_sync_applied ? "true" : "false")
          << " synced_pb=" << report.synced_pb_point_count
          << " synced_mub=" << report.synced_mub_point_count
          << " synced_phb=" << report.synced_phb_point_count
          << " refreshed_points=" << report.compute_report.refreshed_point_count
          << " invalid_points=" << report.compute_report.invalid_point_count
          << " skipped_points=" << report.compute_report.skipped_point_count
          << " touched_overlap=" << (report.touched_overlap_cells ? "true" : "false")
          << " touched_halo=" << (report.touched_halo_cells ? "true" : "false")
          << " terrain_override_used="
          << (report.provider_report.terrain_override_used ? "true" : "false")
          << " terrain_source=" << report.provider_report.terrain_source_name
          << " terrain_provenance=" << report.provider_report.terrain_provenance;
  throw std::runtime_error(message.str());
}

void apply_pressure_refresh(
    const Options& options,
    tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    CandidateReport& report) {
  tywrf::FieldStorage3D<float> direct_alb(candidate.grid.mass_layout());
  const auto metadata = tywrf::io::read_krosa_pressure_refresh_inputs(
      options.template_path,
      candidate.grid,
      direct_alb,
      {.time_index = options.template_time_index});
  require_pressure_refresh_inputs_ready(metadata);

  const std::vector<float> p_before(
      candidate.p.data(), candidate.p.data() + candidate.p.size());
  const std::vector<float> pb_before(
      candidate.pb.data(), candidate.pb.data() + candidate.pb.size());
  const std::vector<float> mub_before(
      candidate.mub.data(), candidate.mub.data() + candidate.mub.size());
  const std::vector<float> phb_before(
      candidate.phb.data(), candidate.phb.data() + candidate.phb.size());
  const auto terrain_override = make_moved_candidate_hgt_terrain_override(output_static);
  tywrf::dynamics::KrosaPressureRefreshHookOptions hook_options{};
  hook_options.terrain_override = &terrain_override;
  auto hook_report = tywrf::dynamics::apply_krosa_moving_nest_pressure_refresh_hook(
      report.remap_plan,
      candidate,
      metadata.metadata,
      hook_options);
  require_pressure_refresh_hook_success(hook_report);
  require_finite_strict_fields(candidate);

  report.pressure_refresh_changed_p_points = changed_points(p_before, candidate.p);
  report.pressure_refresh_changed_pb_points = changed_points(pb_before, candidate.pb);
  report.pressure_refresh_changed_mub_points = changed_points(mub_before, candidate.mub);
  report.pressure_refresh_changed_phb_points = changed_points(phb_before, candidate.phb);
  if (report.pressure_refresh_changed_p_points == 0) {
    throw std::runtime_error("selected-field pressure refresh did not change any P point");
  }
  report.pressure_refresh_metadata_source = options.template_path;
  report.pressure_refresh_metadata_time_index = options.template_time_index;
  report.pressure_refresh = hook_report;
}

[[nodiscard]] CandidateReport build_candidate_state(
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const Options& options,
    const tywrf::State<float>& d01_start,
    const tywrf::State<float>& d02_start,
    tywrf::State<float>& candidate) {
  const auto from_pose = tywrf::nest::make_nest_pose(descriptor, options.from_parent_start);
  const auto to_pose = tywrf::nest::make_nest_pose(descriptor, options.to_parent_start);
  if (!from_pose.child.result.ok()) {
    throw std::runtime_error(
        "invalid --from-parent-start: " + std::string(from_pose.child.result.message));
  }
  if (!to_pose.child.result.ok()) {
    throw std::runtime_error(
        "invalid --to-parent-start: " + std::string(to_pose.child.result.message));
  }

  const auto remap_plan = tywrf::nest::build_remap_plan(from_pose, to_pose);
  if (!remap_plan.ok()) {
    throw std::runtime_error("failed to build remap plan: " + std::string(remap_plan.result.message));
  }

  candidate = d02_start;
  CandidateReport report;
  report.remap_plan = remap_plan;
  report.remap = tywrf::nest::remap_child_state_overlap_only(remap_plan, d02_start, candidate);
  if (!report.remap.ok()) {
    throw std::runtime_error("failed to remap d02 overlap: " + std::string(report.remap.result.message));
  }

  report.exchange =
      tywrf::nest::build_exposed_child_state_exchange_plan(
          remap_plan, static_cast<const tywrf::State<float>&>(candidate).view());
  if (!report.exchange.ok()) {
    throw std::runtime_error(
        "failed to build exposed-cell exchange plan: " +
        std::string(report.exchange.result.message));
  }
  if (!report.exchange.report.requires_parent_interpolation ||
      report.exchange.report.exchange_point_count == 0) {
    throw std::runtime_error("moving-nest pose change exposes no selected-field cells");
  }

  report.interpolation = tywrf::nest::interpolate_parent_to_exposed_child(
      descriptor,
      options.to_parent_start,
      report.exchange,
      d01_start.view(),
      candidate.view());
  if (!report.interpolation.ok()) {
    throw std::runtime_error(
        "failed parent-to-child selected-field interpolation: " +
        std::string(report.interpolation.result.message));
  }
  if (report.interpolation.wrote_overlap || report.interpolation.wrote_halo ||
      report.exchange.report.modifies_overlap || report.exchange.report.modifies_halo) {
    throw std::runtime_error("selected-field interpolation unexpectedly wrote overlap or halo cells");
  }

  report.changed_selected_points = changed_selected_points(d02_start, candidate);
  if (report.changed_selected_points == 0) {
    throw std::runtime_error("selected-field candidate did not change any selected-field point");
  }
  require_finite_strict_fields(candidate);
  return report;
}

void refresh_static_fields(
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const Options& options,
    const StaticFieldSet& d02_start_static,
    const StaticFieldSet& template_static,
    const tywrf::FieldStorage2D<float>& d01_start_hgt,
    StaticFieldSet& output_static,
    CandidateReport& report) {
  output_static = template_static;
  report.static_refresh = tywrf::nest::refresh_moving_nest_static_fields(
      descriptor,
      options.to_parent_start,
      report.remap_plan,
      d02_start_static.xlat.view(),
      d02_start_static.xlong.view(),
      d02_start_static.hgt.view(),
      d01_start_hgt.view(),
      output_static.xlat.view(),
      output_static.xlong.view(),
      output_static.hgt.view());
  if (!report.static_refresh.ok()) {
    throw std::runtime_error(
        "failed moving-nest static refresh: " +
        std::string(report.static_refresh.result.message));
  }
  require_finite_static_fields(output_static);
  report.changed_static_template_points = changed_static_points(template_static, output_static);

  const auto layout = output_static.xlat.layout();
  const auto center_i = layout.i_begin() + layout.active_nx() / 2;
  const auto center_j = layout.j_begin() + layout.active_ny() / 2;
  const auto xlat_view = output_static.xlat.view();
  const auto xlong_view = output_static.xlong.view();
  report.cen_lat = xlat_view(center_i, center_j);
  report.cen_lon = xlong_view(center_i, center_j);
}

void write_text_attr(
    const NetcdfHandle& file,
    const std::string_view name,
    const std::string_view value) {
  file.check(
      nc_put_att_text(
          file.id(),
          NC_GLOBAL,
          std::string(name).c_str(),
          value.size(),
          value.data()),
      "write global text attribute");
}

void write_double_attr(const NetcdfHandle& file, const std::string_view name, const double value) {
  file.check(
      nc_put_att_double(file.id(), NC_GLOBAL, std::string(name).c_str(), NC_DOUBLE, 1, &value),
      "write global numeric attribute");
}

void write_int_attr(
    const NetcdfHandle& file,
    const std::string_view name,
    const std::int32_t value) {
  file.check(
      nc_put_att_int(file.id(), NC_GLOBAL, std::string(name).c_str(), NC_INT, 1, &value),
      "write global integer attribute");
}

void write_pressure_refresh_readiness_attrs(
    const NetcdfHandle& file,
    const PressureRefreshReadiness& readiness) {
  write_text_attr(file, "TYWRF_PRESSURE_REFRESH_READINESS_READY", bool_text(readiness.ready()));
  write_text_attr(
      file,
      "TYWRF_THERMODYNAMIC_BASE_STATE_CONSISTENCY_READY",
      bool_text(readiness.thermodynamic_base_state_consistency_ready));
  write_text_attr(
      file,
      "TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT",
      bool_text(readiness.provider_terrain_uses_moved_candidate_hgt));
  write_text_attr(
      file,
      "TYWRF_PROVIDER_BASE_STATE_RECONSTRUCT_OK",
      bool_text(readiness.provider_base_state_reconstruct_ok));
  write_text_attr(
      file,
      "TYWRF_BASE_STATE_SYNC_CONTRACT_OK",
      bool_text(readiness.base_state_sync_contract_ok));
  write_text_attr(
      file, "TYWRF_BASE_STATE_SYNC_DRY_RUN", bool_text(readiness.base_state_sync_dry_run));
  write_text_attr(
      file, "TYWRF_BASE_STATE_SYNC_APPLIED", bool_text(readiness.base_state_sync_applied));
  write_double_attr(
      file, "TYWRF_WOULD_SYNC_PB_POINT_COUNT", static_cast<double>(readiness.would_sync_pb_point_count));
  write_double_attr(
      file,
      "TYWRF_WOULD_SYNC_MUB_POINT_COUNT",
      static_cast<double>(readiness.would_sync_mub_point_count));
  write_double_attr(
      file,
      "TYWRF_WOULD_SYNC_PHB_POINT_COUNT",
      static_cast<double>(readiness.would_sync_phb_point_count));
  write_double_attr(
      file,
      "TYWRF_SYNC_OVERLAP_WRITE_COUNT",
      static_cast<double>(readiness.sync_overlap_write_count));
  write_double_attr(
      file, "TYWRF_SYNC_HALO_WRITE_COUNT", static_cast<double>(readiness.sync_halo_write_count));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_REFRESH_DRY_RUN_COMPUTE_CALLED",
      bool_text(readiness.pressure_refresh_compute_called));
  write_text_attr(
      file, "TYWRF_PRESSURE_COMPUTE_DRY_RUN", bool_text(readiness.pressure_compute_dry_run));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_DRY_RUN_CALLED",
      bool_text(readiness.pressure_compute_dry_run_called));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_DRY_RUN_OK",
      bool_text(readiness.pressure_compute_dry_run_ok));
  write_double_attr(
      file,
      "TYWRF_WOULD_REFRESH_P_POINT_COUNT",
      static_cast<double>(readiness.would_refresh_p_point_count));
  write_double_attr(
      file,
      "TYWRF_DRY_RUN_INVALID_P_POINT_COUNT",
      static_cast<double>(readiness.dry_run_invalid_p_point_count));
  write_double_attr(
      file,
      "TYWRF_DRY_RUN_SKIPPED_P_POINT_COUNT",
      static_cast<double>(readiness.dry_run_skipped_p_point_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_DRY_RUN_REPORT_TARGET_COLUMN_COUNT",
      static_cast<double>(readiness.pressure_compute_dry_run_report_target_column_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_DRY_RUN_REPORT_REFRESHED_POINT_COUNT",
      static_cast<double>(readiness.pressure_compute_dry_run_report_refreshed_point_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_DRY_RUN_REPORT_INVALID_POINT_COUNT",
      static_cast<double>(readiness.pressure_compute_dry_run_report_invalid_point_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_DRY_RUN_REPORT_SKIPPED_POINT_COUNT",
      static_cast<double>(readiness.pressure_compute_dry_run_report_skipped_point_count));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_DRY_RUN_REPORT_TOUCHED_OVERLAP_CELLS",
      bool_text(readiness.pressure_compute_dry_run_report_touched_overlap_cells));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_DRY_RUN_REPORT_TOUCHED_HALO_CELLS",
      bool_text(readiness.pressure_compute_dry_run_report_touched_halo_cells));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_REFRESH_DRY_RUN_APPLIED",
      bool_text(readiness.pressure_refresh_applied));
  if (!readiness.provider_terrain_source_name.empty()) {
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_SOURCE",
        readiness.provider_terrain_source_name);
  }
  if (!readiness.provider_terrain_provenance.empty()) {
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_PROVENANCE",
        readiness.provider_terrain_provenance);
  }
}

[[nodiscard]] std::string join_variables(const std::vector<std::string>& variables) {
  std::ostringstream joined;
  for (std::size_t index = 0; index < variables.size(); ++index) {
    if (index != 0) {
      joined << ",";
    }
    joined << variables[index];
  }
  return joined.str();
}

void stamp_gate_metadata(
    const Options& options,
    const Resolution resolution,
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const CandidateReport& report) {
  const bool experimental_pressure_refresh_apply =
      options.experimental_pressure_refresh_apply && report.pressure_refresh.has_value();
  std::optional<PressureRefreshReadiness> pressure_refresh_readiness;
  if (report.pressure_refresh.has_value()) {
    pressure_refresh_readiness = evaluate_pressure_refresh_readiness(report);
  }
  const bool pressure_refresh_applied_to_candidate =
      report.pressure_refresh.has_value() && !experimental_pressure_refresh_apply;
  const bool pressure_refresh_gate_candidate =
      pressure_refresh_applied_to_candidate && pressure_refresh_readiness.has_value() &&
      pressure_refresh_readiness->ready();
  const bool gate_candidate =
      report.pressure_refresh.has_value() ? pressure_refresh_gate_candidate : true;
  const NetcdfHandle file(options.output_path, NetcdfHandle::Mode::write);
  file.check(nc_redef(file.id()), "enter define mode");
  write_double_attr(file, "DX", resolution.dx);
  write_double_attr(file, "DY", resolution.dy);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ONLY",
      bool_text(experimental_pressure_refresh_apply));
  write_text_attr(file, "TYWRF_GATE_CANDIDATE", bool_text(gate_candidate));
  write_text_attr(file, "TYWRF_INTEGRATOR_OUTPUT", bool_text(gate_candidate));
  write_text_attr(file, "TYWRF_VALIDATION_GATE_ONLY", "false");
  write_text_attr(
      file,
      "TYWRF_CANDIDATE_KIND",
      experimental_pressure_refresh_apply
          ? "selected_field_pressure_refresh_experimental_apply_v0"
          : "selected_field_integrator_v0");
  if (experimental_pressure_refresh_apply) {
    write_text_attr(file, "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY", "true");
  }
  write_text_attr(file, "TYWRF_CANDIDATE_DOMAIN", "d02");
  write_text_attr(file, "TYWRF_D02_RESOLUTION_CHECK", "d02_2km");
  write_text_attr(file, "TYWRF_CYCLE_START", options.cycle_start);
  write_text_attr(file, "TYWRF_CYCLE_END", options.cycle_end);
  write_text_attr(file, "TYWRF_D01_START_STATE_SOURCE", options.d01_start_state_path.string());
  write_text_attr(file, "TYWRF_D02_START_STATE_SOURCE", options.d02_start_state_path.string());
  write_text_attr(file, "TYWRF_TEMPLATE_SOURCE", options.template_path.string());
  write_text_attr(file, "TYWRF_STATE_VARIABLES", join_variables(options.variables));
  write_text_attr(
      file,
      "TYWRF_PARENT_INTERPOLATED_STATE_VARIABLES",
      join_variables(parent_interpolation_variable_names()));
  write_text_attr(file, "TYWRF_FROM_PARENT_START", std::to_string(options.from_parent_start.i_parent_start) + "," + std::to_string(options.from_parent_start.j_parent_start));
  write_text_attr(file, "TYWRF_TO_PARENT_START", std::to_string(options.to_parent_start.i_parent_start) + "," + std::to_string(options.to_parent_start.j_parent_start));
  write_int_attr(file, "I_PARENT_START", options.to_parent_start.i_parent_start);
  write_int_attr(file, "J_PARENT_START", options.to_parent_start.j_parent_start);
  write_double_attr(file, "CEN_LAT", report.cen_lat);
  write_double_attr(file, "CEN_LON", report.cen_lon);
  write_double_attr(file, "TYWRF_PARENT_GRID_RATIO", static_cast<double>(descriptor.parent_grid_ratio));
  write_double_attr(
      file,
      "TYWRF_SELECTED_FIELD_CHANGED_POINTS",
      static_cast<double>(report.changed_selected_points));
  write_text_attr(file, "TYWRF_STATIC_REFRESH_APPLIED", "true");
  write_text_attr(
      file,
      "TYWRF_STATIC_REFRESH_METHOD",
      "overlap_shift_xlat_xlong_extrapolate_hgt_parent_bilinear_v0");
  write_text_attr(file, "TYWRF_STATIC_REFRESH_USES_REFERENCE_END", "false");
  write_text_attr(file, "TYWRF_STATIC_REFRESH_D02_START_SOURCE", options.d02_start_state_path.string());
  write_text_attr(file, "TYWRF_STATIC_REFRESH_D01_HGT_SOURCE", options.d01_start_state_path.string());
  write_text_attr(file, "TYWRF_STATIC_REFRESH_TEMPLATE_SOURCE", options.template_path.string());
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_OVERLAP_CELLS",
      static_cast<double>(report.static_refresh.overlap_cell_count));
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_EXPOSED_CELLS",
      static_cast<double>(report.static_refresh.exposed_cell_count));
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_COORD_EXTRAPOLATED_CELLS",
      static_cast<double>(report.static_refresh.coordinate_extrapolated_cell_count));
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_HGT_PARENT_INTERPOLATED_CELLS",
      static_cast<double>(report.static_refresh.parent_hgt_interpolated_cell_count));
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_CHANGED_TEMPLATE_POINTS",
      static_cast<double>(report.changed_static_template_points));
  write_double_attr(
      file,
      "TYWRF_EXPOSED_EXCHANGE_POINTS",
      static_cast<double>(report.exchange.report.exchange_point_count));
  write_double_attr(
      file,
      "TYWRF_INTERPOLATED_POINTS",
      static_cast<double>(report.interpolation.interpolated_point_count));
  if (report.pressure_refresh.has_value()) {
    const auto& pressure = *report.pressure_refresh;
    const auto parity = pressure_refresh_report_parity(report, pressure);
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_OPT_IN", "true");
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_APPLIED", "true");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS",
        experimental_pressure_refresh_apply ? "experimental_apply_test_only"
                                            : "applied_to_candidate");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY",
        experimental_pressure_refresh_apply ? "true" : "false");
    if (pressure_refresh_readiness.has_value()) {
      write_pressure_refresh_readiness_attrs(file, *pressure_refresh_readiness);
    }
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_PROVIDER_OK", "true");
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_STAGING_OK", "true");
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED", "true");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TERRAIN_OVERRIDE_USED",
        pressure.provider_report.terrain_override_used ? "true" : "false");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE",
        pressure.provider_report.terrain_source_name);
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE",
        pressure.provider_report.terrain_provenance);
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE",
        report.pressure_refresh_metadata_source.string());
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_METADATA_TIME_INDEX",
        static_cast<double>(report.pressure_refresh_metadata_time_index));
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_HELPER_NAME", "refresh_krosa_moving_nest_pressure");
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS",
        static_cast<double>(pressure.synced_pb_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS",
        static_cast<double>(pressure.synced_mub_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS",
        static_cast<double>(pressure.synced_phb_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT",
        static_cast<double>(pressure.compute_report.target_column_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT",
        static_cast<double>(pressure.compute_report.refreshed_column_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT",
        static_cast<double>(pressure.compute_report.refreshed_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT",
        static_cast<double>(pressure.compute_report.skipped_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT",
        static_cast<double>(pressure.compute_report.invalid_point_count));
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS",
        pressure.compute_report.touched_overlap_cells ? "true" : "false");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TOUCHED_HALO_CELLS",
        pressure.compute_report.touched_halo_cells ? "true" : "false");
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS",
        static_cast<double>(pressure.compute_report.refreshed_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS",
        static_cast<double>(report.pressure_refresh_changed_p_points));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS",
        static_cast<double>(report.pressure_refresh_changed_pb_points));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS",
        static_cast<double>(report.pressure_refresh_changed_mub_points));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS",
        static_cast<double>(report.pressure_refresh_changed_phb_points));
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT",
        parity.changed_p_matches_refreshed_point_count ? "true" : "false");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO",
        parity.invalid_and_skipped_points_zero ? "true" : "false");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED",
        parity.overlap_halo_untouched ? "true" : "false");
  }
  if (experimental_pressure_refresh_apply) {
    write_text_attr(
        file,
        "TYWRF_CANDIDATE_MESSAGE",
        "Experimental selected-field pressure-refresh apply seam output for tool tests only; "
        "not a validation gate pass or normal integrator output. The apply path used "
        "moved_candidate_HGT terrain override and refreshed exposed P/PB/PHB/MUB.");
  } else if (report.pressure_refresh.has_value()) {
    write_text_attr(
        file,
        "TYWRF_CANDIDATE_MESSAGE",
        "Selected-field moving-nest candidate from start states only; U/V/T/PH/MU/QVAPOR "
        "exposed cells are parent interpolated, XLAT/XLONG/HGT are refreshed from "
        "start-state pose data, and pressure refresh applied to exposed P/PB/PHB/MUB.");
  } else {
    write_text_attr(
        file,
        "TYWRF_CANDIDATE_MESSAGE",
        "Selected-field moving-nest candidate from start states only; U/V/T/PH/MU/QVAPOR "
        "exposed cells are parent interpolated, XLAT/XLONG/HGT are refreshed from "
        "start-state pose data, and P/PB/PHB/MUB remain finite d02 start-state ownership.");
  }
  file.check(nc_enddef(file.id()), "leave define mode");
}

[[nodiscard]] std::string json_escape(const std::string_view value) {
  std::ostringstream escaped;
  for (const char c : value) {
    switch (c) {
      case '\\':
        escaped << "\\\\";
        break;
      case '"':
        escaped << "\\\"";
        break;
      case '\n':
        escaped << "\\n";
        break;
      case '\r':
        escaped << "\\r";
        break;
      case '\t':
        escaped << "\\t";
        break;
      default:
        escaped << c;
        break;
    }
  }
  return escaped.str();
}

void write_json_string(
    std::ostream& stream,
    const std::string_view name,
    const std::string_view value,
    const bool comma,
    const bool pretty) {
  stream << (pretty ? "  \"" : "\"") << name << "\": \"" << json_escape(value) << "\"";
  if (comma) {
    stream << ",";
  }
  stream << (pretty ? "\n" : "");
}

void write_json_bool(
    std::ostream& stream,
    const std::string_view name,
    const bool value,
    const bool comma,
    const bool pretty) {
  stream << (pretty ? "  \"" : "\"") << name << "\": " << (value ? "true" : "false");
  if (comma) {
    stream << ",";
  }
  stream << (pretty ? "\n" : "");
}

void write_json_number(
    std::ostream& stream,
    const std::string_view name,
    const double value,
    const bool comma,
    const bool pretty) {
  stream << (pretty ? "  \"" : "\"") << name << "\": " << value;
  if (comma) {
    stream << ",";
  }
  stream << (pretty ? "\n" : "");
}

void write_pressure_refresh_readiness_json(
    std::ostream& stream,
    const PressureRefreshReadiness& readiness,
    const bool pretty) {
  write_json_bool(stream, "pressure_refresh_readiness_ready", readiness.ready(), true, pretty);
  write_json_bool(
      stream,
      "thermodynamic_base_state_consistency_ready",
      readiness.thermodynamic_base_state_consistency_ready,
      true,
      pretty);
  write_json_bool(
      stream,
      "provider_terrain_uses_moved_candidate_hgt",
      readiness.provider_terrain_uses_moved_candidate_hgt,
      true,
      pretty);
  write_json_bool(
      stream,
      "provider_base_state_reconstruct_ok",
      readiness.provider_base_state_reconstruct_ok,
      true,
      pretty);
  write_json_bool(
      stream, "base_state_sync_contract_ok", readiness.base_state_sync_contract_ok, true, pretty);
  write_json_bool(stream, "base_state_sync_dry_run", readiness.base_state_sync_dry_run, true, pretty);
  write_json_bool(stream, "base_state_sync_applied", readiness.base_state_sync_applied, true, pretty);
  write_json_number(
      stream,
      "would_sync_pb_point_count",
      static_cast<double>(readiness.would_sync_pb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "would_sync_mub_point_count",
      static_cast<double>(readiness.would_sync_mub_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "would_sync_phb_point_count",
      static_cast<double>(readiness.would_sync_phb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "sync_overlap_write_count",
      static_cast<double>(readiness.sync_overlap_write_count),
      true,
      pretty);
  write_json_number(
      stream,
      "sync_halo_write_count",
      static_cast<double>(readiness.sync_halo_write_count),
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_refresh_dry_run_compute_called",
      readiness.pressure_refresh_compute_called,
      true,
      pretty);
  write_json_bool(
      stream, "pressure_compute_dry_run", readiness.pressure_compute_dry_run, true, pretty);
  write_json_bool(
      stream,
      "pressure_compute_dry_run_called",
      readiness.pressure_compute_dry_run_called,
      true,
      pretty);
  write_json_bool(
      stream, "pressure_compute_dry_run_ok", readiness.pressure_compute_dry_run_ok, true, pretty);
  write_json_number(
      stream,
      "would_refresh_p_point_count",
      static_cast<double>(readiness.would_refresh_p_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "dry_run_invalid_p_point_count",
      static_cast<double>(readiness.dry_run_invalid_p_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "dry_run_skipped_p_point_count",
      static_cast<double>(readiness.dry_run_skipped_p_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_compute_dry_run_report_target_column_count",
      static_cast<double>(readiness.pressure_compute_dry_run_report_target_column_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_compute_dry_run_report_refreshed_point_count",
      static_cast<double>(readiness.pressure_compute_dry_run_report_refreshed_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_compute_dry_run_report_invalid_point_count",
      static_cast<double>(readiness.pressure_compute_dry_run_report_invalid_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_compute_dry_run_report_skipped_point_count",
      static_cast<double>(readiness.pressure_compute_dry_run_report_skipped_point_count),
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_compute_dry_run_report_touched_overlap_cells",
      readiness.pressure_compute_dry_run_report_touched_overlap_cells,
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_compute_dry_run_report_touched_halo_cells",
      readiness.pressure_compute_dry_run_report_touched_halo_cells,
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_refresh_dry_run_applied",
      readiness.pressure_refresh_applied,
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_refresh_readiness_provider_terrain_source",
      readiness.provider_terrain_source_name,
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_refresh_readiness_provider_terrain_provenance",
      readiness.provider_terrain_provenance,
      true,
      pretty);
}

void print_report(
    const Options& options,
    const Resolution resolution,
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const CandidateReport& report) {
  const bool pretty = options.pretty;
  const bool experimental_pressure_refresh_apply =
      options.experimental_pressure_refresh_apply && report.pressure_refresh.has_value();
  std::optional<PressureRefreshReadiness> pressure_refresh_readiness;
  if (report.pressure_refresh.has_value()) {
    pressure_refresh_readiness = evaluate_pressure_refresh_readiness(report);
  }
  const bool pressure_refresh_applied_to_candidate =
      report.pressure_refresh.has_value() && !experimental_pressure_refresh_apply;
  const bool pressure_refresh_gate_candidate =
      pressure_refresh_applied_to_candidate && pressure_refresh_readiness.has_value() &&
      pressure_refresh_readiness->ready();
  const bool gate_candidate =
      report.pressure_refresh.has_value() ? pressure_refresh_gate_candidate : true;
  std::cout << "{" << (pretty ? "\n" : "");
  write_json_string(
      std::cout,
      "status",
      experimental_pressure_refresh_apply
          ? "selected_field_pressure_refresh_experimental_apply_generated"
          : "selected_field_candidate_generated",
      true,
      pretty);
  write_json_string(
      std::cout,
      "candidate_kind",
      experimental_pressure_refresh_apply
          ? "selected_field_pressure_refresh_experimental_apply_v0"
          : "selected_field_integrator_v0",
      true,
      pretty);
  write_json_bool(std::cout, "gate_candidate", gate_candidate, true, pretty);
  write_json_bool(std::cout, "integrator_output", gate_candidate, true, pretty);
  write_json_bool(std::cout, "validation_gate_only", false, true, pretty);
  if (experimental_pressure_refresh_apply) {
    write_json_bool(std::cout, "experimental_pressure_refresh_apply", true, true, pretty);
  }
  write_json_string(std::cout, "d01_start_state", options.d01_start_state_path.string(), true, pretty);
  write_json_string(std::cout, "d02_start_state", options.d02_start_state_path.string(), true, pretty);
  write_json_string(std::cout, "candidate", options.output_path.string(), true, pretty);
  write_json_number(std::cout, "dx_m", resolution.dx, true, pretty);
  write_json_number(std::cout, "dy_m", resolution.dy, true, pretty);
  write_json_number(
      std::cout, "parent_grid_ratio", static_cast<double>(descriptor.parent_grid_ratio), true, pretty);
  write_json_number(
      std::cout,
      "selected_field_changed_points",
      static_cast<double>(report.changed_selected_points),
      true,
      pretty);
  write_json_number(
      std::cout,
      "static_refresh_changed_template_points",
      static_cast<double>(report.changed_static_template_points),
      true,
      pretty);
  write_json_bool(std::cout, "static_refresh_uses_reference_end", false, true, pretty);
  write_json_number(
      std::cout,
      "static_refresh_exposed_cells",
      static_cast<double>(report.static_refresh.exposed_cell_count),
      true,
      pretty);
  write_json_number(
      std::cout,
      "exposed_exchange_points",
      static_cast<double>(report.exchange.report.exchange_point_count),
      true,
      pretty);
  write_json_number(
      std::cout,
      "interpolated_points",
      static_cast<double>(report.interpolation.interpolated_point_count),
      true,
      pretty);
  write_json_string(
      std::cout,
      "parent_interpolated_state_variables",
      join_variables(parent_interpolation_variable_names()),
      report.pressure_refresh.has_value(),
      pretty);
  if (report.pressure_refresh.has_value()) {
    const auto& pressure = *report.pressure_refresh;
    const auto parity = pressure_refresh_report_parity(report, pressure);
    write_json_bool(std::cout, "pressure_refresh_opt_in", true, true, pretty);
    write_json_bool(std::cout, "pressure_refresh_applied", true, true, pretty);
    write_json_bool(std::cout, "pressure_refresh_provider_ok", true, true, pretty);
    write_json_bool(std::cout, "pressure_refresh_staging_ok", true, true, pretty);
    write_json_bool(std::cout, "pressure_refresh_compute_called", true, true, pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_experimental_apply",
        experimental_pressure_refresh_apply,
        true,
        pretty);
    write_json_string(
        std::cout,
        "pressure_refresh_integration_status",
        experimental_pressure_refresh_apply ? "experimental_apply_test_only"
                                            : "applied_to_candidate",
        true,
        pretty);
    if (pressure_refresh_readiness.has_value()) {
      write_pressure_refresh_readiness_json(std::cout, *pressure_refresh_readiness, pretty);
    }
    write_json_bool(
        std::cout,
        "pressure_refresh_terrain_override_used",
        pressure.provider_report.terrain_override_used,
        true,
        pretty);
    write_json_string(
        std::cout,
        "pressure_refresh_terrain_source",
        pressure.provider_report.terrain_source_name,
        true,
        pretty);
    write_json_string(
        std::cout,
        "pressure_refresh_terrain_provenance",
        pressure.provider_report.terrain_provenance,
        true,
        pretty);
    write_json_string(
        std::cout,
        "pressure_refresh_metadata_source",
        report.pressure_refresh_metadata_source.string(),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_metadata_time_index",
        static_cast<double>(report.pressure_refresh_metadata_time_index),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_target_column_count",
        static_cast<double>(pressure.compute_report.target_column_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_refreshed_column_count",
        static_cast<double>(pressure.compute_report.refreshed_column_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_refreshed_point_count",
        static_cast<double>(pressure.compute_report.refreshed_point_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_skipped_point_count",
        static_cast<double>(pressure.compute_report.skipped_point_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_invalid_point_count",
        static_cast<double>(pressure.compute_report.invalid_point_count),
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_touched_overlap_cells",
        pressure.compute_report.touched_overlap_cells,
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_touched_halo_cells",
        pressure.compute_report.touched_halo_cells,
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_refreshed_p_points",
        static_cast<double>(pressure.compute_report.refreshed_point_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_changed_p_points",
        static_cast<double>(report.pressure_refresh_changed_p_points),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_changed_pb_points",
        static_cast<double>(report.pressure_refresh_changed_pb_points),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_changed_mub_points",
        static_cast<double>(report.pressure_refresh_changed_mub_points),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_changed_phb_points",
        static_cast<double>(report.pressure_refresh_changed_phb_points),
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_changed_p_matches_refreshed_point_count",
        parity.changed_p_matches_refreshed_point_count,
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_invalid_and_skipped_points_zero",
        parity.invalid_and_skipped_points_zero,
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_overlap_halo_untouched",
        parity.overlap_halo_untouched,
        false,
        pretty);
  }
  std::cout << "}" << (pretty ? "\n" : "\n");
}

int run(Options options) {
  require_path_exists(options.d01_start_state_path, "d01 start-state input");
  require_path_exists(options.d02_start_state_path, "d02 start-state input");
  require_path_exists(options.template_path, "template input");
  if (options.variables.empty()) {
    options.variables = default_output_variables(options.d02_start_state_path);
  }
  require_output_variables(options.variables);
  if (options.output_path.has_parent_path()) {
    std::filesystem::create_directories(options.output_path.parent_path());
  }

  const auto d01_resolution = read_resolution(options.d01_start_state_path);
  const auto d02_resolution = read_resolution(options.d02_start_state_path);
  const auto template_resolution = read_resolution(options.template_path);
  require_d02_resolution(options.d02_start_state_path, d02_resolution);
  require_d02_resolution(options.template_path, template_resolution);

  const auto d01_grid = tywrf::io::derive_grid_from_wrf_file(
      options.d01_start_state_path, tywrf::uniform_halo_3d(0));
  const auto d02_grid = tywrf::io::derive_grid_from_wrf_file(
      options.d02_start_state_path, tywrf::uniform_halo_3d(0));
  const auto descriptor = make_descriptor(d01_resolution, d02_resolution, d01_grid, d02_grid);

  const auto d01_start_hgt =
      load_hgt_field(options.d01_start_state_path, d01_grid, options.d01_time_index);
  const auto d02_start_static =
      load_static_fields(options.d02_start_state_path, d02_grid, options.d02_time_index);
  const auto template_static =
      load_static_fields(options.template_path, d02_grid, options.template_time_index);
  require_finite_storage("d01 HGT", d01_start_hgt);
  require_finite_static_fields(d02_start_static);
  require_finite_static_fields(template_static);

  tywrf::State<float> d01_start(d01_grid);
  tywrf::State<float> d02_start(d02_grid);
  tywrf::io::load_wrf_state(
      options.d01_start_state_path,
      d01_start,
      {.time_index = options.d01_time_index, .variables = parent_interpolation_variable_names()});
  tywrf::io::load_wrf_state(
      options.d02_start_state_path,
      d02_start,
      {.time_index = options.d02_time_index,
       .variables = selected_d02_read_variables(options.d02_start_state_path)});
  require_finite_strict_fields(d02_start);

  tywrf::State<float> candidate(d02_grid);
  auto report = build_candidate_state(descriptor, options, d01_start, d02_start, candidate);
  StaticFieldSet output_static(d02_grid);
  refresh_static_fields(
      descriptor,
      options,
      d02_start_static,
      template_static,
      d01_start_hgt,
      output_static,
      report);
  if (options.pressure_refresh) {
    probe_pressure_refresh_provider_readiness(options, candidate, output_static, report);
    probe_pressure_refresh_dry_run_contract(options, candidate, output_static, report);
    if (!options.experimental_pressure_refresh_apply) {
      require_pressure_refresh_ready_for_compute(report);
    }
    apply_pressure_refresh(options, candidate, output_static, report);
  }

  tywrf::io::write_wrf_state(
      options.output_path,
      candidate,
      {
          .time_index = options.output_time_index,
          .variables = options.variables,
          .template_path = options.template_path,
          .template_time_index = options.template_time_index,
          .times_value = options.times_value,
      });
  overwrite_output_static_fields(options.output_path, output_static, options.output_time_index);
  stamp_gate_metadata(options, template_resolution, descriptor, report);
  print_report(options, template_resolution, descriptor, report);
  return 0;
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    const auto options = parse_options(argc, argv);
    return run(options);
  } catch (const std::exception& error) {
    std::cerr << "tywrf_selected_field_cycle: " << error.what() << "\n\n" << usage();
    return 2;
  }
}
