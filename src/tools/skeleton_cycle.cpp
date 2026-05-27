#include "tywrf/field_view.hpp"
#include "tywrf/io/netcdf_schema.hpp"
#include "tywrf/io/wrf_state_io.hpp"
#include "tywrf/nest/nest_interface.hpp"
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
#include <system_error>
#include <vector>

namespace {

constexpr double kD02TargetDxMeters = 2000.0;
constexpr double kD02ResolutionToleranceMeters = 0.5;

const std::vector<std::string>& template_variable_names() {
  static const std::vector<std::string> names = {"Times", "XLAT", "XLONG", "HGT"};
  return names;
}

struct Options {
  std::filesystem::path state_path;
  std::filesystem::path template_path;
  std::filesystem::path output_path;
  std::string domain = "d02";
  std::string cycle_start;
  std::string cycle_end;
  std::string times_value;
  std::size_t state_time_index = 0;
  std::size_t template_time_index = 0;
  std::size_t output_time_index = 0;
  std::vector<std::string> variables;
  bool pretty = false;
  bool diagnostic_remap_overlap = false;
  std::optional<tywrf::nest::ParentChildPosition> from_parent_start;
  std::optional<tywrf::nest::ParentChildPosition> to_parent_start;
};

struct DiagnosticRemapResult {
  tywrf::nest::RemapPlan plan;
  tywrf::nest::ChildStateRemapReport report;
};

struct Resolution {
  double dx = 0.0;
  double dy = 0.0;
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
  tywrf_skeleton_cycle --state <cycle-start wrfout> --template <cycle-start/static wrfout> --output <candidate wrfout> [options]

Options:
  --domain d02                    Only d02 is currently supported.
  --state-time-index N            Time index read from --state; default 0.
  --template-time-index N         Time index copied from --template for static fields and default Times; default 0.
  --output-time-index N           Time index written in --output; default 0.
  --cycle-start TEXT              Optional report/metadata cycle start string.
  --cycle-end TEXT                Optional report/metadata cycle end string.
  --times TEXT                    Override output Times; default copies template Times.
  --variables A,B,C               Combined State/template variables to write.
  --diagnostic-remap-overlap      Remap d02 state overlap only; leaves exposed cells unfilled.
  --from-parent-start I,J         One-based d02 parent start before diagnostic remap.
  --to-parent-start I,J           One-based d02 parent start after diagnostic remap.
  --pretty                        Pretty-print JSON report.
  --help                          Show this help.

For d02 moving-nest persistence smoke tests, pass the cycle-start wrfout as
--template and use --times for the cycle-end timestamp. This keeps
XLAT/XLONG/HGT consistent with the cycle-start state. Only use cycle-end
template coordinates after a real integrator state has already been remapped to
the end nest pose.

This executable is a skeleton candidate writer only. It does not run dynamics
or physics and marks output metadata as skeleton/not-physical/integrator_output=false.
The diagnostic remap path is opt-in, d02 2 km only, writes NaN into cells not
covered by overlap, and is not a validation-gate candidate.
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

[[nodiscard]] std::int32_t parse_i32(
    const std::string& value,
    const std::string_view option) {
  std::size_t parsed_chars = 0;
  const auto parsed = std::stol(value, &parsed_chars);
  if (parsed_chars != value.size()) {
    throw std::invalid_argument(std::string(option) + " expects an integer");
  }
  if (parsed < std::numeric_limits<std::int32_t>::min() ||
      parsed > std::numeric_limits<std::int32_t>::max()) {
    throw std::out_of_range(std::string(option) + " is outside int32 range");
  }
  return static_cast<std::int32_t>(parsed);
}

[[nodiscard]] tywrf::nest::ParentChildPosition parse_parent_start(
    const std::string& value,
    const std::string_view option) {
  const auto comma = value.find(',');
  if (comma == std::string::npos || value.find(',', comma + 1) != std::string::npos) {
    throw std::invalid_argument(std::string(option) + " expects I,J");
  }
  return {
      parse_i32(value.substr(0, comma), option),
      parse_i32(value.substr(comma + 1), option),
      tywrf::nest::IndexBase::one_based,
  };
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

[[nodiscard]] Options parse_options(const int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string arg = argv[index];
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
    } else if (arg == "--state") {
      options.state_path = require_value(arg);
    } else if (arg == "--template") {
      options.template_path = require_value(arg);
    } else if (arg == "--output") {
      options.output_path = require_value(arg);
    } else if (arg == "--domain") {
      options.domain = require_value(arg);
    } else if (arg == "--cycle-start") {
      options.cycle_start = require_value(arg);
    } else if (arg == "--cycle-end") {
      options.cycle_end = require_value(arg);
    } else if (arg == "--times") {
      options.times_value = require_value(arg);
    } else if (arg == "--state-time-index") {
      options.state_time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--template-time-index") {
      options.template_time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--output-time-index") {
      options.output_time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--variables") {
      options.variables = split_variables(require_value(arg));
    } else if (arg == "--diagnostic-remap-overlap") {
      options.diagnostic_remap_overlap = true;
    } else if (arg == "--from-parent-start") {
      options.from_parent_start = parse_parent_start(require_value(arg), arg);
    } else if (arg == "--to-parent-start") {
      options.to_parent_start = parse_parent_start(require_value(arg), arg);
    } else if (arg == "--pretty") {
      options.pretty = true;
    } else {
      throw std::invalid_argument("unknown option: " + arg);
    }
  }

  if (options.state_path.empty()) {
    throw std::invalid_argument("--state is required");
  }
  if (options.template_path.empty()) {
    throw std::invalid_argument("--template is required");
  }
  if (options.output_path.empty()) {
    throw std::invalid_argument("--output is required");
  }
  if (options.domain != "d02") {
    throw std::invalid_argument("tywrf_skeleton_cycle currently supports d02 only");
  }
  if (options.diagnostic_remap_overlap &&
      (!options.from_parent_start.has_value() || !options.to_parent_start.has_value())) {
    throw std::invalid_argument(
        "--diagnostic-remap-overlap requires --from-parent-start and --to-parent-start");
  }
  if (options.variables.empty()) {
    options.variables = template_variable_names();
    const auto state_variables = tywrf::io::wrf_state_writable_field_names();
    options.variables.insert(options.variables.end(), state_variables.begin(), state_variables.end());
  }
  return options;
}

[[nodiscard]] bool contains(
    const std::vector<std::string>& values,
    const std::string_view value) {
  return std::find(values.begin(), values.end(), value) != values.end();
}

[[nodiscard]] std::vector<std::string> state_variables_from_selection(
    const std::vector<std::string>& variables) {
  const auto writable = tywrf::io::wrf_state_writable_field_names();
  std::vector<std::string> state_variables;
  for (const auto& variable : variables) {
    if (contains(writable, variable)) {
      state_variables.push_back(variable);
    } else if (!contains(template_variable_names(), variable)) {
      throw std::invalid_argument("unsupported skeleton output variable: " + variable);
    }
  }
  return state_variables;
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
    throw std::runtime_error("d02 input must carry DX and DY attributes: " + path.string());
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

[[nodiscard]] bool paths_refer_to_same_file(
    const std::filesystem::path& lhs,
    const std::filesystem::path& rhs) {
  if (std::filesystem::absolute(lhs).lexically_normal() ==
      std::filesystem::absolute(rhs).lexically_normal()) {
    return true;
  }
  std::error_code error;
  const bool equivalent = std::filesystem::equivalent(lhs, rhs, error);
  return !error && equivalent;
}

[[nodiscard]] std::string times_source(const Options& options) {
  if (!options.times_value.empty()) {
    return "--times:" + options.times_value;
  }
  return "template:" + options.template_path.string();
}

[[nodiscard]] std::string static_coords_match_state_source(const Options& options) {
  return paths_refer_to_same_file(options.state_path, options.template_path) ? "same_file"
                                                                            : "not_verified";
}

[[nodiscard]] std::string state_static_consistency(const Options& options) {
  if (paths_refer_to_same_file(options.state_path, options.template_path)) {
    return "static_coords_same_file_as_state_source";
  }
  return "static_coords_from_template_not_verified_against_state_source";
}

[[nodiscard]] std::string parent_start_string(
    const tywrf::nest::ParentChildPosition& position) {
  std::ostringstream value;
  value << position.i_parent_start << "," << position.j_parent_start;
  return value.str();
}

[[nodiscard]] std::string bool_string(const bool value) {
  return value ? "true" : "false";
}

void mark_skeleton_output(
    const Options& options,
    const Resolution resolution,
    const std::vector<std::string>& state_variables,
    const std::optional<DiagnosticRemapResult>& remap_result) {
  const bool diagnostic_remap = remap_result.has_value();
  const NetcdfHandle file(options.output_path, NetcdfHandle::Mode::write);
  file.check(nc_redef(file.id()), "enter define mode");
  write_double_attr(file, "DX", resolution.dx);
  write_double_attr(file, "DY", resolution.dy);
  write_text_attr(
      file,
      "TYWRF_CANDIDATE_KIND",
      diagnostic_remap ? "cpp_skeleton_remap_overlap_diagnostic" : "cpp_skeleton_candidate");
  write_text_attr(file, "TYWRF_SKELETON", "true");
  write_text_attr(file, "TYWRF_NOT_PHYSICAL", "true");
  write_text_attr(file, "TYWRF_SKELETON_NOT_PHYSICAL", "true");
  write_text_attr(file, "TYWRF_INTEGRATOR_OUTPUT", "false");
  write_text_attr(file, "TYWRF_VALIDATION_GATE_ONLY", diagnostic_remap ? "false" : "true");
  write_text_attr(file, "TYWRF_EXPECTED_TO_MEET_THRESHOLDS", "false");
  write_text_attr(file, "TYWRF_CANDIDATE_DOMAIN", options.domain);
  write_text_attr(file, "TYWRF_CANDIDATE_SOURCE", options.state_path.string());
  write_text_attr(file, "TYWRF_STATE_SOURCE", options.state_path.string());
  write_text_attr(file, "TYWRF_TEMPLATE_SOURCE", options.template_path.string());
  write_text_attr(file, "TYWRF_STATIC_SOURCE", options.template_path.string());
  write_text_attr(
      file,
      "TYWRF_TEMPLATE_SOURCE_ROLE",
      "static_fields_xlat_xlong_hgt_and_default_times_when_no_times_override");
  write_text_attr(file, "TYWRF_TIMES_SOURCE", times_source(options));
  write_text_attr(
      file,
      "TYWRF_STATIC_COORDS_MATCH_STATE_SOURCE",
      static_coords_match_state_source(options));
  write_text_attr(file, "TYWRF_STATE_STATIC_CONSISTENCY", state_static_consistency(options));
  write_text_attr(file, "TYWRF_D02_RESOLUTION_CHECK", "d02_2km");
  write_text_attr(
      file,
      "TYWRF_CANDIDATE_MESSAGE",
      diagnostic_remap
          ? "C++ diagnostic remap-overlap candidate; overlap cells copied from old d02 state, "
            "exposed cells require parent fill and may contain NaN. Not physical and not a "
            "validation-gate candidate."
          : "C++ skeleton candidate generated from cycle-start reference state; not physical "
            "and not an integrator result. For d02 moving-nest persistence, static coordinates "
            "must come from the cycle-start nest pose and Times should be overridden to the "
            "cycle end.");
  if (!options.cycle_start.empty()) {
    write_text_attr(file, "TYWRF_CYCLE_START", options.cycle_start);
  }
  if (!options.cycle_end.empty()) {
    write_text_attr(file, "TYWRF_CYCLE_END", options.cycle_end);
  }
  std::ostringstream loaded;
  for (std::size_t index = 0; index < state_variables.size(); ++index) {
    if (index != 0) {
      loaded << ",";
    }
    loaded << state_variables[index];
  }
  write_text_attr(file, "TYWRF_STATE_VARIABLES", loaded.str());
  if (diagnostic_remap) {
    const auto& report = remap_result->report;
    write_text_attr(file, "TYWRF_DIAGNOSTIC_REMAP_OVERLAP", "true");
    write_text_attr(file, "TYWRF_DIAGNOSTIC_ONLY", "true");
    write_text_attr(file, "TYWRF_GATE_CANDIDATE", "false");
    write_text_attr(file, "TYWRF_NEEDS_PARENT_FILL", bool_string(report.needs_parent_fill));
    write_text_attr(file, "TYWRF_REMAP_EXPOSED_CELLS_FILLED_BY_PARENT", "false");
    write_text_attr(file, "TYWRF_UNFILLED_EXPOSED_CELLS", "nan_sentinel_pending_parent_fill");
    write_text_attr(file, "TYWRF_REMAP_FROM_PARENT_START", parent_start_string(*options.from_parent_start));
    write_text_attr(file, "TYWRF_REMAP_TO_PARENT_START", parent_start_string(*options.to_parent_start));
    write_double_attr(
        file,
        "TYWRF_REMAP_COPIED_FIELD_COUNT",
        static_cast<double>(report.copied_field_count));
    write_double_attr(
        file,
        "TYWRF_REMAP_COPIED_POINT_COUNT",
        static_cast<double>(report.copied_point_count));
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

void write_json_uint64(
    std::ostream& stream,
    const std::string_view name,
    const std::uint64_t value,
    const bool comma,
    const bool pretty) {
  stream << (pretty ? "  \"" : "\"") << name << "\": " << value;
  if (comma) {
    stream << ",";
  }
  stream << (pretty ? "\n" : "");
}

void write_json_int32(
    std::ostream& stream,
    const std::string_view name,
    const std::int32_t value,
    const bool comma,
    const bool pretty) {
  stream << (pretty ? "  \"" : "\"") << name << "\": " << value;
  if (comma) {
    stream << ",";
  }
  stream << (pretty ? "\n" : "");
}

void write_json_array(
    std::ostream& stream,
    const std::string_view name,
    const std::vector<std::string>& values,
    const bool comma,
    const bool pretty) {
  stream << (pretty ? "  \"" : "\"") << name << "\": [";
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) {
      stream << ", ";
    }
    stream << "\"" << json_escape(values[index]) << "\"";
  }
  stream << "]";
  if (comma) {
    stream << ",";
  }
  stream << (pretty ? "\n" : "");
}

void print_report(
    const Options& options,
    const Resolution resolution,
    const std::vector<std::string>& state_variables,
    const std::optional<DiagnosticRemapResult>& remap_result) {
  const bool pretty = options.pretty;
  const bool diagnostic_remap = remap_result.has_value();
  std::cout << "{" << (pretty ? "\n" : "");
  write_json_string(
      std::cout,
      "status",
      diagnostic_remap ? "diagnostic_remap_overlap_generated" : "skeleton_candidate_generated",
      true,
      pretty);
  write_json_string(
      std::cout,
      "candidate_kind",
      diagnostic_remap ? "cpp_skeleton_remap_overlap_diagnostic" : "cpp_skeleton_candidate",
      true,
      pretty);
  write_json_bool(std::cout, "skeleton", true, true, pretty);
  write_json_bool(std::cout, "not_physical", true, true, pretty);
  write_json_bool(std::cout, "integrator_output", false, true, pretty);
  write_json_bool(std::cout, "validation_gate_only", !diagnostic_remap, true, pretty);
  write_json_string(std::cout, "domain", options.domain, true, pretty);
  write_json_string(std::cout, "source", options.state_path.string(), true, pretty);
  write_json_string(std::cout, "template", options.template_path.string(), true, pretty);
  write_json_string(std::cout, "state_source", options.state_path.string(), true, pretty);
  write_json_string(std::cout, "static_source", options.template_path.string(), true, pretty);
  write_json_string(std::cout, "times_source", times_source(options), true, pretty);
  write_json_string(
      std::cout,
      "static_coords_match_state_source",
      static_coords_match_state_source(options),
      true,
      pretty);
  write_json_string(
      std::cout,
      "state_static_consistency",
      state_static_consistency(options),
      true,
      pretty);
  write_json_string(std::cout, "candidate", options.output_path.string(), true, pretty);
  write_json_string(std::cout, "d02_resolution_check", "d02_2km", true, pretty);
  write_json_number(std::cout, "dx_m", resolution.dx, true, pretty);
  write_json_number(std::cout, "dy_m", resolution.dy, true, pretty);
  write_json_array(std::cout, "variables", options.variables, true, pretty);
  write_json_array(std::cout, "state_variables", state_variables, true, pretty);
  if (diagnostic_remap) {
    const auto& report = remap_result->report;
    write_json_bool(std::cout, "diagnostic_remap_overlap", true, true, pretty);
    write_json_bool(std::cout, "diagnostic_only", true, true, pretty);
    write_json_bool(std::cout, "gate_candidate", false, true, pretty);
    write_json_bool(std::cout, "needs_parent_fill", report.needs_parent_fill, true, pretty);
    write_json_bool(std::cout, "exposed_cells_filled_by_parent", false, true, pretty);
    write_json_string(
        std::cout,
        "unfilled_exposed_cells",
        "nan_sentinel_pending_parent_fill",
        true,
        pretty);
    write_json_string(
        std::cout,
        "from_parent_start",
        parent_start_string(*options.from_parent_start),
        true,
        pretty);
    write_json_string(
        std::cout,
        "to_parent_start",
        parent_start_string(*options.to_parent_start),
        true,
        pretty);
    write_json_int32(
        std::cout,
        "from_i_parent_start",
        options.from_parent_start->i_parent_start,
        true,
        pretty);
    write_json_int32(
        std::cout,
        "from_j_parent_start",
        options.from_parent_start->j_parent_start,
        true,
        pretty);
    write_json_int32(
        std::cout,
        "to_i_parent_start",
        options.to_parent_start->i_parent_start,
        true,
        pretty);
    write_json_int32(
        std::cout,
        "to_j_parent_start",
        options.to_parent_start->j_parent_start,
        true,
        pretty);
    write_json_uint64(
        std::cout,
        "remap_copied_field_count",
        report.copied_field_count,
        true,
        pretty);
    write_json_uint64(
        std::cout,
        "remap_copied_point_count",
        report.copied_point_count,
        true,
        pretty);
  }
  write_json_string(
      std::cout,
      "message",
      diagnostic_remap
          ? "C++ diagnostic overlap remap only; exposed cells require parent fill before any physical validation."
          : "C++ skeleton candidate only; output is not physical and not a TyWRF-Core integrator result.",
      false,
      pretty);
  std::cout << "}" << (pretty ? "\n" : "\n");
}

template <typename Storage>
void fill_storage(Storage& storage, const float value) {
  std::fill_n(storage.data(), storage.size(), value);
}

void fill_state(tywrf::State<float>& state, const float value) {
  fill_storage(state.u, value);
  fill_storage(state.v, value);
  fill_storage(state.w, value);
  fill_storage(state.ph, value);
  fill_storage(state.phb, value);
  fill_storage(state.t, value);
  fill_storage(state.p, value);
  fill_storage(state.pb, value);
  fill_storage(state.qvapor, value);
  fill_storage(state.qcloud, value);
  fill_storage(state.qrain, value);
  fill_storage(state.qice, value);
  fill_storage(state.qsnow, value);
  fill_storage(state.qgraup, value);
  fill_storage(state.qnice, value);
  fill_storage(state.qnrain, value);
  fill_storage(state.mu, value);
  fill_storage(state.mub, value);
  fill_storage(state.psfc, value);
  fill_storage(state.u10, value);
  fill_storage(state.v10, value);
  fill_storage(state.t2, value);
  fill_storage(state.q2, value);
  fill_storage(state.rainc, value);
  fill_storage(state.rainnc, value);
}

void require_krosa_d02_grid(
    const tywrf::Grid& grid,
    const tywrf::nest::ParentChildDescriptor& descriptor) {
  const auto& config = grid.config();
  if (config.mass_nx == descriptor.child.mass_nx &&
      config.mass_ny == descriptor.child.mass_ny) {
    return;
  }
  std::ostringstream message;
  message << "diagnostic remap overlap currently supports KROSA d02 active mass grid "
          << descriptor.child.mass_nx << "x" << descriptor.child.mass_ny << "; state has "
          << config.mass_nx << "x" << config.mass_ny;
  throw std::runtime_error(message.str());
}

void require_nest_result_ok(
    const tywrf::nest::NestResult result,
    const std::string_view operation) {
  if (result.ok()) {
    return;
  }
  std::ostringstream message;
  message << operation << " failed with "
          << tywrf::nest::nest_status_name(result.status) << ": " << result.message;
  throw std::runtime_error(message.str());
}

[[nodiscard]] DiagnosticRemapResult apply_diagnostic_remap_overlap(
    const Options& options,
    const tywrf::State<float>& old_state,
    tywrf::State<float>& new_state) {
  const auto descriptor = tywrf::nest::make_krosa_parent_child_descriptor();
  require_krosa_d02_grid(old_state.grid, descriptor);

  const auto from_pose = tywrf::nest::make_domain_pose(descriptor, *options.from_parent_start);
  require_nest_result_ok(from_pose.result, "from parent-start pose");
  const auto to_pose = tywrf::nest::make_domain_pose(descriptor, *options.to_parent_start);
  require_nest_result_ok(to_pose.result, "to parent-start pose");
  auto result = DiagnosticRemapResult{};
  result.plan = tywrf::nest::build_remap_plan(from_pose, to_pose);
  require_nest_result_ok(result.plan.result, "diagnostic overlap remap plan");

  fill_state(new_state, std::numeric_limits<float>::quiet_NaN());
  result.report = tywrf::nest::remap_child_state_overlap_only(result.plan, old_state, new_state);
  require_nest_result_ok(result.report.result, "diagnostic overlap remap");
  return result;
}

int run(const Options& options) {
  if (!std::filesystem::exists(options.state_path)) {
    throw std::runtime_error("state input does not exist: " + options.state_path.string());
  }
  if (!std::filesystem::exists(options.template_path)) {
    throw std::runtime_error("template input does not exist: " + options.template_path.string());
  }
  if (options.output_path.has_parent_path()) {
    std::filesystem::create_directories(options.output_path.parent_path());
  }

  const auto state_resolution = read_resolution(options.state_path);
  const auto template_resolution = read_resolution(options.template_path);
  require_d02_resolution(options.state_path, state_resolution);
  require_d02_resolution(options.template_path, template_resolution);

  const auto state_variables = state_variables_from_selection(options.variables);
  const auto grid = tywrf::io::derive_grid_from_wrf_file(options.state_path, tywrf::uniform_halo_3d(0));
  tywrf::State<float> state(grid);
  tywrf::io::load_wrf_state(
      options.state_path,
      state,
      {.time_index = options.state_time_index, .variables = state_variables});

  std::optional<DiagnosticRemapResult> remap_result;
  if (options.diagnostic_remap_overlap) {
    tywrf::State<float> remapped_state(grid);
    remap_result = apply_diagnostic_remap_overlap(options, state, remapped_state);
    tywrf::io::write_wrf_state(
        options.output_path,
        remapped_state,
        {
            .time_index = options.output_time_index,
            .variables = options.variables,
            .template_path = options.template_path,
            .template_time_index = options.template_time_index,
            .times_value = options.times_value,
        });
  } else {
    tywrf::io::write_wrf_state(
        options.output_path,
        state,
        {
            .time_index = options.output_time_index,
            .variables = options.variables,
            .template_path = options.template_path,
            .template_time_index = options.template_time_index,
            .times_value = options.times_value,
        });
  }

  mark_skeleton_output(options, template_resolution, state_variables, remap_result);
  print_report(options, template_resolution, state_variables, remap_result);
  return 0;
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    const auto options = parse_options(argc, argv);
    return run(options);
  } catch (const std::exception& error) {
    std::cerr << "tywrf_skeleton_cycle: " << error.what() << "\n\n" << usage();
    return 2;
  }
}
