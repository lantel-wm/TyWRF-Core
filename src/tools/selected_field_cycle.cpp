#include "tywrf/io/wrf_state_io.hpp"
#include "tywrf/nest/parent_child_interpolation.hpp"
#include "tywrf/nest/state_exchange.hpp"
#include "tywrf/nest/state_remap.hpp"
#include "tywrf/state.hpp"

#include <netcdf.h>

#include <algorithm>
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
  static const std::vector<std::string> names = {"U", "V", "MU", "QVAPOR"};
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
  bool pretty = false;
  std::vector<std::string> variables;
};

struct Resolution {
  double dx = 0.0;
  double dy = 0.0;
};

struct CandidateReport {
  std::uint64_t changed_selected_points = 0;
  tywrf::nest::ChildStateRemapReport remap;
  tywrf::nest::StateExchangePlan exchange;
  tywrf::nest::ParentChildInterpolationReport interpolation;
};

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
  --pretty                       Pretty-print JSON report.
  --help                         Show this help.

Oracle inputs are not accepted: --end-state, --reference-end, --d01-end-state,
and --d02-end-state are rejected. This tool remaps the d02 start state from the
old moving-nest pose to the new pose, fills newly exposed U/V/MU/QVAPOR cells
from d01 start-state parent interpolation, and preserves strict plus available
d02 start-state diagnostic fields without using reference-end truth. It is a
selected-field integrator candidate, not a WRF-exact physics result.
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

template <typename Storage>
[[nodiscard]] std::uint64_t changed_points(const Storage& before, const Storage& after) {
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
  const NetcdfHandle file(options.output_path, NetcdfHandle::Mode::write);
  file.check(nc_redef(file.id()), "enter define mode");
  write_double_attr(file, "DX", resolution.dx);
  write_double_attr(file, "DY", resolution.dy);
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ONLY", "false");
  write_text_attr(file, "TYWRF_GATE_CANDIDATE", "true");
  write_text_attr(file, "TYWRF_INTEGRATOR_OUTPUT", "true");
  write_text_attr(file, "TYWRF_VALIDATION_GATE_ONLY", "false");
  write_text_attr(file, "TYWRF_CANDIDATE_KIND", "selected_field_integrator_v0");
  write_text_attr(file, "TYWRF_CANDIDATE_DOMAIN", "d02");
  write_text_attr(file, "TYWRF_D02_RESOLUTION_CHECK", "d02_2km");
  write_text_attr(file, "TYWRF_CYCLE_START", options.cycle_start);
  write_text_attr(file, "TYWRF_CYCLE_END", options.cycle_end);
  write_text_attr(file, "TYWRF_D01_START_STATE_SOURCE", options.d01_start_state_path.string());
  write_text_attr(file, "TYWRF_D02_START_STATE_SOURCE", options.d02_start_state_path.string());
  write_text_attr(file, "TYWRF_TEMPLATE_SOURCE", options.template_path.string());
  write_text_attr(file, "TYWRF_STATE_VARIABLES", join_variables(options.variables));
  write_text_attr(file, "TYWRF_FROM_PARENT_START", std::to_string(options.from_parent_start.i_parent_start) + "," + std::to_string(options.from_parent_start.j_parent_start));
  write_text_attr(file, "TYWRF_TO_PARENT_START", std::to_string(options.to_parent_start.i_parent_start) + "," + std::to_string(options.to_parent_start.j_parent_start));
  write_double_attr(file, "TYWRF_PARENT_GRID_RATIO", static_cast<double>(descriptor.parent_grid_ratio));
  write_double_attr(
      file,
      "TYWRF_SELECTED_FIELD_CHANGED_POINTS",
      static_cast<double>(report.changed_selected_points));
  write_double_attr(
      file,
      "TYWRF_EXPOSED_EXCHANGE_POINTS",
      static_cast<double>(report.exchange.report.exchange_point_count));
  write_double_attr(
      file,
      "TYWRF_INTERPOLATED_POINTS",
      static_cast<double>(report.interpolation.interpolated_point_count));
  write_text_attr(
      file,
      "TYWRF_CANDIDATE_MESSAGE",
      "Selected-field moving-nest candidate from start states only; U/V/MU/QVAPOR exposed "
      "cells are parent interpolated, and T/PH/P are preserved from finite d02 start-state data.");
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

void print_report(
    const Options& options,
    const Resolution resolution,
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const CandidateReport& report) {
  const bool pretty = options.pretty;
  std::cout << "{" << (pretty ? "\n" : "");
  write_json_string(std::cout, "status", "selected_field_candidate_generated", true, pretty);
  write_json_string(std::cout, "candidate_kind", "selected_field_integrator_v0", true, pretty);
  write_json_bool(std::cout, "gate_candidate", true, true, pretty);
  write_json_bool(std::cout, "integrator_output", true, true, pretty);
  write_json_bool(std::cout, "validation_gate_only", false, true, pretty);
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
      "exposed_exchange_points",
      static_cast<double>(report.exchange.report.exchange_point_count),
      true,
      pretty);
  write_json_number(
      std::cout,
      "interpolated_points",
      static_cast<double>(report.interpolation.interpolated_point_count),
      false,
      pretty);
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
  const auto report = build_candidate_state(descriptor, options, d01_start, d02_start, candidate);

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
