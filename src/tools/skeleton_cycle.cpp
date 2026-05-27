#include "tywrf/field_view.hpp"
#include "tywrf/io/netcdf_schema.hpp"
#include "tywrf/io/wrf_state_io.hpp"
#include "tywrf/state.hpp"

#include <netcdf.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
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
  tywrf_skeleton_cycle --state <cycle-start wrfout> --template <cycle-end wrfout> --output <candidate wrfout> [options]

Options:
  --domain d02                    Only d02 is currently supported.
  --state-time-index N            Time index read from --state; default 0.
  --template-time-index N         Time index copied from --template; default 0.
  --output-time-index N           Time index written in --output; default 0.
  --cycle-start TEXT              Optional report/metadata cycle start string.
  --cycle-end TEXT                Optional report/metadata cycle end string.
  --times TEXT                    Override output Times; default copies template Times.
  --variables A,B,C               Combined State/template variables to write.
  --pretty                        Pretty-print JSON report.
  --help                          Show this help.

This executable is a skeleton candidate writer only. It does not run dynamics
or physics and marks output metadata as skeleton/not-physical/integrator_output=false.
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

void mark_skeleton_output(
    const Options& options,
    const Resolution resolution,
    const std::vector<std::string>& state_variables) {
  const NetcdfHandle file(options.output_path, NetcdfHandle::Mode::write);
  file.check(nc_redef(file.id()), "enter define mode");
  write_double_attr(file, "DX", resolution.dx);
  write_double_attr(file, "DY", resolution.dy);
  write_text_attr(file, "TYWRF_CANDIDATE_KIND", "cpp_skeleton_candidate");
  write_text_attr(file, "TYWRF_SKELETON", "true");
  write_text_attr(file, "TYWRF_NOT_PHYSICAL", "true");
  write_text_attr(file, "TYWRF_SKELETON_NOT_PHYSICAL", "true");
  write_text_attr(file, "TYWRF_INTEGRATOR_OUTPUT", "false");
  write_text_attr(file, "TYWRF_VALIDATION_GATE_ONLY", "true");
  write_text_attr(file, "TYWRF_EXPECTED_TO_MEET_THRESHOLDS", "false");
  write_text_attr(file, "TYWRF_CANDIDATE_DOMAIN", options.domain);
  write_text_attr(file, "TYWRF_CANDIDATE_SOURCE", options.state_path.string());
  write_text_attr(file, "TYWRF_TEMPLATE_SOURCE", options.template_path.string());
  write_text_attr(file, "TYWRF_D02_RESOLUTION_CHECK", "d02_2km");
  write_text_attr(
      file,
      "TYWRF_CANDIDATE_MESSAGE",
      "C++ skeleton candidate generated from cycle-start reference state; not physical and "
      "not an integrator result.");
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
    const std::vector<std::string>& state_variables) {
  const bool pretty = options.pretty;
  std::cout << "{" << (pretty ? "\n" : "");
  write_json_string(std::cout, "status", "skeleton_candidate_generated", true, pretty);
  write_json_string(std::cout, "candidate_kind", "cpp_skeleton_candidate", true, pretty);
  write_json_bool(std::cout, "skeleton", true, true, pretty);
  write_json_bool(std::cout, "not_physical", true, true, pretty);
  write_json_bool(std::cout, "integrator_output", false, true, pretty);
  write_json_bool(std::cout, "validation_gate_only", true, true, pretty);
  write_json_string(std::cout, "domain", options.domain, true, pretty);
  write_json_string(std::cout, "source", options.state_path.string(), true, pretty);
  write_json_string(std::cout, "template", options.template_path.string(), true, pretty);
  write_json_string(std::cout, "candidate", options.output_path.string(), true, pretty);
  write_json_string(std::cout, "d02_resolution_check", "d02_2km", true, pretty);
  write_json_number(std::cout, "dx_m", resolution.dx, true, pretty);
  write_json_number(std::cout, "dy_m", resolution.dy, true, pretty);
  write_json_array(std::cout, "variables", options.variables, true, pretty);
  write_json_array(std::cout, "state_variables", state_variables, true, pretty);
  write_json_string(
      std::cout,
      "message",
      "C++ skeleton candidate only; output is not physical and not a TyWRF-Core integrator result.",
      false,
      pretty);
  std::cout << "}" << (pretty ? "\n" : "\n");
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

  mark_skeleton_output(options, template_resolution, state_variables);
  print_report(options, template_resolution, state_variables);
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
