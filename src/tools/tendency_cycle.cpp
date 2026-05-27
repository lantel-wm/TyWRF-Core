#include "tywrf/dynamics/dynamics_loop.hpp"
#include "tywrf/dynamics/state_stepper.hpp"
#include "tywrf/field_view.hpp"
#include "tywrf/io/wrf_state_io.hpp"
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
  std::filesystem::path start_state_path;
  std::filesystem::path end_state_path;
  std::filesystem::path template_path;
  std::filesystem::path output_path;
  std::string domain = "d02";
  std::string cycle_start;
  std::string cycle_end;
  std::string times_value;
  std::size_t start_time_index = 0;
  std::size_t end_time_index = 0;
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
  tywrf_tendency_cycle --start-state <cycle-start wrfout> --end-state <cycle-end wrfout> --template <static/template wrfout> --output <candidate wrfout> --cycle-start TEXT --cycle-end TEXT --times TEXT [options]

Options:
  --domain d02                    Only d02 is supported; the d02 2 km guard is enforced.
  --start-time-index N            Time index read from --start-state; default 0.
  --end-time-index N              Time index read from --end-state; default 0.
  --template-time-index N         Time index copied from --template for static fields; default 0.
  --output-time-index N           Time index written in --output; default 0.
  --variables A,B,C               Combined State/template variables to write.
  --pretty                        Pretty-print JSON report.
  --help                          Show this help.

This executable is a diagnostic reference-delta oracle for the first 10 minute
KROSA d02 validation segment. It derives explicit tendencies from
(end-state - start-state) / 600 s and routes them through the skeleton
explicit-tendency hook. The result is non-physical by construction and is
stamped so strict gates reject it by default.
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
    } else if (arg == "--start-state") {
      options.start_state_path = require_value(arg);
    } else if (arg == "--end-state") {
      options.end_state_path = require_value(arg);
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
    } else if (arg == "--start-time-index") {
      options.start_time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--end-time-index") {
      options.end_time_index = parse_size(require_value(arg), arg);
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

  if (options.start_state_path.empty()) {
    throw std::invalid_argument("--start-state is required");
  }
  if (options.end_state_path.empty()) {
    throw std::invalid_argument("--end-state is required");
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
  if (options.domain != "d02") {
    throw std::invalid_argument("tywrf_tendency_cycle currently supports d02 only");
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
      throw std::invalid_argument("unsupported tendency-cycle output variable: " + variable);
    }
  }
  if (state_variables.empty()) {
    throw std::invalid_argument(
        "tendency-cycle oracle requires at least one State-backed output variable");
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

[[nodiscard]] int parse_fixed_int(
    const std::string& value,
    const std::size_t offset,
    const std::size_t count,
    const std::string_view label) {
  int parsed = 0;
  for (std::size_t index = 0; index < count; ++index) {
    const char c = value[offset + index];
    if (!std::isdigit(static_cast<unsigned char>(c))) {
      throw std::invalid_argument(std::string(label) + " must match YYYY-MM-DD_HH:MM:SS");
    }
    parsed = parsed * 10 + (c - '0');
  }
  return parsed;
}

[[nodiscard]] std::int64_t days_from_civil(
    int year,
    const unsigned month,
    const unsigned day) noexcept {
  year -= month <= 2U ? 1 : 0;
  const int era = (year >= 0 ? year : year - 399) / 400;
  const unsigned year_of_era = static_cast<unsigned>(year - era * 400);
  const unsigned month_prime = month > 2U ? month - 3U : month + 9U;
  const unsigned day_of_year = (153U * month_prime + 2U) / 5U + day - 1U;
  const unsigned day_of_era =
      year_of_era * 365U + year_of_era / 4U - year_of_era / 100U + day_of_year;
  return static_cast<std::int64_t>(era) * 146097L +
         static_cast<std::int64_t>(day_of_era) - 719468L;
}

[[nodiscard]] std::int64_t parse_wrf_timestamp_seconds(
    const std::string& value,
    const std::string_view label) {
  if (value.size() != 19 || value[4] != '-' || value[7] != '-' || value[10] != '_' ||
      value[13] != ':' || value[16] != ':') {
    throw std::invalid_argument(std::string(label) + " must match YYYY-MM-DD_HH:MM:SS");
  }
  const int year = parse_fixed_int(value, 0, 4, label);
  const int month = parse_fixed_int(value, 5, 2, label);
  const int day = parse_fixed_int(value, 8, 2, label);
  const int hour = parse_fixed_int(value, 11, 2, label);
  const int minute = parse_fixed_int(value, 14, 2, label);
  const int second = parse_fixed_int(value, 17, 2, label);
  if (month < 1 || month > 12 || day < 1 || day > 31 || hour < 0 || hour > 23 ||
      minute < 0 || minute > 59 || second < 0 || second > 60) {
    throw std::invalid_argument(std::string(label) + " has an out-of-range timestamp field");
  }

  const auto days = days_from_civil(year, static_cast<unsigned>(month), static_cast<unsigned>(day));
  return days * 86'400L + static_cast<std::int64_t>(hour) * 3600L +
         static_cast<std::int64_t>(minute) * 60L + static_cast<std::int64_t>(second);
}

[[nodiscard]] std::int64_t cycle_dt_seconds(const Options& options) {
  const auto start_seconds = parse_wrf_timestamp_seconds(options.cycle_start, "--cycle-start");
  const auto end_seconds = parse_wrf_timestamp_seconds(options.cycle_end, "--cycle-end");
  const auto dt_seconds = end_seconds - start_seconds;
  if (dt_seconds <= 0) {
    throw std::invalid_argument("--cycle-end must be after --cycle-start");
  }
  return dt_seconds;
}

template <typename Storage>
void compute_storage_tendency(
    const std::string_view name,
    const Storage& start,
    const Storage& end,
    Storage& tendency,
    const float inverse_dt) {
  if (start.size() != end.size() || start.size() != tendency.size()) {
    throw std::runtime_error("state storage size mismatch while deriving tendency for " +
                             std::string(name));
  }
  for (std::size_t index = 0; index < start.size(); ++index) {
    tendency.data()[index] = (end.data()[index] - start.data()[index]) * inverse_dt;
  }
}

void compute_state_tendencies(
    const tywrf::State<float>& start,
    const tywrf::State<float>& end,
    tywrf::State<float>& tendency,
    const std::int64_t dt_seconds) {
  if (dt_seconds <= 0 ||
      static_cast<long double>(dt_seconds) >
          static_cast<long double>(std::numeric_limits<float>::max())) {
    throw std::runtime_error("invalid tendency dt seconds");
  }
  const float inverse_dt = 1.0F / static_cast<float>(dt_seconds);
  compute_storage_tendency("U", start.u, end.u, tendency.u, inverse_dt);
  compute_storage_tendency("V", start.v, end.v, tendency.v, inverse_dt);
  compute_storage_tendency("W", start.w, end.w, tendency.w, inverse_dt);
  compute_storage_tendency("PH", start.ph, end.ph, tendency.ph, inverse_dt);
  compute_storage_tendency("PHB", start.phb, end.phb, tendency.phb, inverse_dt);
  compute_storage_tendency("T", start.t, end.t, tendency.t, inverse_dt);
  compute_storage_tendency("P", start.p, end.p, tendency.p, inverse_dt);
  compute_storage_tendency("PB", start.pb, end.pb, tendency.pb, inverse_dt);
  compute_storage_tendency("QVAPOR", start.qvapor, end.qvapor, tendency.qvapor, inverse_dt);
  compute_storage_tendency("QCLOUD", start.qcloud, end.qcloud, tendency.qcloud, inverse_dt);
  compute_storage_tendency("QRAIN", start.qrain, end.qrain, tendency.qrain, inverse_dt);
  compute_storage_tendency("QICE", start.qice, end.qice, tendency.qice, inverse_dt);
  compute_storage_tendency("QSNOW", start.qsnow, end.qsnow, tendency.qsnow, inverse_dt);
  compute_storage_tendency("QGRAUP", start.qgraup, end.qgraup, tendency.qgraup, inverse_dt);
  compute_storage_tendency("QNICE", start.qnice, end.qnice, tendency.qnice, inverse_dt);
  compute_storage_tendency("QNRAIN", start.qnrain, end.qnrain, tendency.qnrain, inverse_dt);
  compute_storage_tendency("MU", start.mu, end.mu, tendency.mu, inverse_dt);
  compute_storage_tendency("MUB", start.mub, end.mub, tendency.mub, inverse_dt);
  compute_storage_tendency("PSFC", start.psfc, end.psfc, tendency.psfc, inverse_dt);
  compute_storage_tendency("U10", start.u10, end.u10, tendency.u10, inverse_dt);
  compute_storage_tendency("V10", start.v10, end.v10, tendency.v10, inverse_dt);
  compute_storage_tendency("T2", start.t2, end.t2, tendency.t2, inverse_dt);
  compute_storage_tendency("Q2", start.q2, end.q2, tendency.q2, inverse_dt);
  compute_storage_tendency("RAINC", start.rainc, end.rainc, tendency.rainc, inverse_dt);
  compute_storage_tendency("RAINNC", start.rainnc, end.rainnc, tendency.rainnc, inverse_dt);
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

void stamp_oracle_metadata(
    const Options& options,
    const Resolution resolution,
    const std::vector<std::string>& state_variables,
    const tywrf::dynamics::SkeletonStateStepReport& report,
    const std::int64_t dt_seconds) {
  const NetcdfHandle file(options.output_path, NetcdfHandle::Mode::write);
  file.check(nc_redef(file.id()), "enter define mode");
  write_double_attr(file, "DX", resolution.dx);
  write_double_attr(file, "DY", resolution.dy);
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ONLY", "true");
  write_text_attr(file, "TYWRF_GATE_CANDIDATE", "false");
  write_text_attr(file, "TYWRF_INTEGRATOR_OUTPUT", "false");
  write_text_attr(file, "TYWRF_VALIDATION_GATE_ONLY", "true");
  write_text_attr(file, "TYWRF_CANDIDATE_KIND", "reference_delta_oracle");
  write_text_attr(file, "TYWRF_NOT_PHYSICAL", "true");
  write_text_attr(file, "TYWRF_CANDIDATE_DOMAIN", options.domain);
  write_text_attr(file, "TYWRF_D02_RESOLUTION_CHECK", "d02_2km");
  write_text_attr(file, "TYWRF_CYCLE_START", options.cycle_start);
  write_text_attr(file, "TYWRF_CYCLE_END", options.cycle_end);
  write_double_attr(file, "TYWRF_DT_SECONDS", static_cast<double>(dt_seconds));
  write_double_attr(file, "TYWRF_REFERENCE_DELTA_DT_SECONDS", static_cast<double>(dt_seconds));
  write_double_attr(
      file,
      "TYWRF_D02_EXPLICIT_TENDENCY_APPLY_COUNT",
      static_cast<double>(report.d02.explicit_tendency_apply_count));
  write_double_attr(
      file,
      "TYWRF_D02_TENDENCY_DT_SECONDS",
      static_cast<double>(report.d02.tendency_dt_seconds));
  write_text_attr(file, "TYWRF_START_STATE_SOURCE", options.start_state_path.string());
  write_text_attr(file, "TYWRF_END_STATE_SOURCE", options.end_state_path.string());
  write_text_attr(file, "TYWRF_TEMPLATE_SOURCE", options.template_path.string());
  write_text_attr(file, "TYWRF_STATE_VARIABLES", join_variables(state_variables));
  write_text_attr(
      file,
      "TYWRF_CANDIDATE_MESSAGE",
      "Diagnostic reference-delta oracle routed through the skeleton explicit-tendency hook; "
      "not physical and not a TyWRF-Core integrator output.");
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
    const std::vector<std::string>& state_variables,
    const tywrf::dynamics::SkeletonStateStepReport& step_report,
    const std::int64_t dt_seconds) {
  const bool pretty = options.pretty;
  std::cout << "{" << (pretty ? "\n" : "");
  write_json_string(std::cout, "status", "reference_delta_oracle_generated", true, pretty);
  write_json_string(std::cout, "candidate_kind", "reference_delta_oracle", true, pretty);
  write_json_bool(std::cout, "diagnostic_only", true, true, pretty);
  write_json_bool(std::cout, "gate_candidate", false, true, pretty);
  write_json_bool(std::cout, "integrator_output", false, true, pretty);
  write_json_bool(std::cout, "validation_gate_only", true, true, pretty);
  write_json_bool(std::cout, "not_physical", true, true, pretty);
  write_json_string(std::cout, "domain", options.domain, true, pretty);
  write_json_string(std::cout, "start_state", options.start_state_path.string(), true, pretty);
  write_json_string(std::cout, "end_state", options.end_state_path.string(), true, pretty);
  write_json_string(std::cout, "template", options.template_path.string(), true, pretty);
  write_json_string(std::cout, "candidate", options.output_path.string(), true, pretty);
  write_json_string(std::cout, "cycle_start", options.cycle_start, true, pretty);
  write_json_string(std::cout, "cycle_end", options.cycle_end, true, pretty);
  write_json_number(std::cout, "dt_seconds", static_cast<double>(dt_seconds), true, pretty);
  write_json_number(std::cout, "dx_m", resolution.dx, true, pretty);
  write_json_number(std::cout, "dy_m", resolution.dy, true, pretty);
  write_json_number(
      std::cout,
      "d02_explicit_tendency_apply_count",
      static_cast<double>(step_report.d02.explicit_tendency_apply_count),
      true,
      pretty);
  write_json_number(
      std::cout,
      "d02_tendency_dt_seconds",
      static_cast<double>(step_report.d02.tendency_dt_seconds),
      true,
      pretty);
  write_json_array(std::cout, "variables", options.variables, true, pretty);
  write_json_array(std::cout, "state_variables", state_variables, true, pretty);
  write_json_string(
      std::cout,
      "message",
      "Reference-delta oracle only; output is deliberately not a physical integrator candidate.",
      false,
      pretty);
  std::cout << "}" << (pretty ? "\n" : "\n");
}

void require_path_exists(const std::filesystem::path& path, const std::string_view label) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error(std::string(label) + " does not exist: " + path.string());
  }
}

int run(const Options& options) {
  require_path_exists(options.start_state_path, "start-state input");
  require_path_exists(options.end_state_path, "end-state input");
  require_path_exists(options.template_path, "template input");
  if (options.output_path.has_parent_path()) {
    std::filesystem::create_directories(options.output_path.parent_path());
  }

  const auto start_resolution = read_resolution(options.start_state_path);
  const auto end_resolution = read_resolution(options.end_state_path);
  const auto template_resolution = read_resolution(options.template_path);
  require_d02_resolution(options.start_state_path, start_resolution);
  require_d02_resolution(options.end_state_path, end_resolution);
  require_d02_resolution(options.template_path, template_resolution);

  const auto dt_seconds = cycle_dt_seconds(options);
  const auto loop_config = tywrf::dynamics::make_krosa_10min_validation_loop_config();
  if (dt_seconds != loop_config.timing.segment_seconds) {
    std::ostringstream message;
    message << "tywrf_tendency_cycle is only for the first 10 minute segment; cycle dt is "
            << dt_seconds << " s but KROSA validation config expects "
            << loop_config.timing.segment_seconds << " s";
    throw std::runtime_error(message.str());
  }

  const auto state_variables = state_variables_from_selection(options.variables);
  const auto grid =
      tywrf::io::derive_grid_from_wrf_file(options.start_state_path, tywrf::uniform_halo_3d(0));
  tywrf::State<float> start_state(grid);
  tywrf::State<float> end_state(grid);
  tywrf::State<float> tendency(grid);
  tywrf::io::load_wrf_state(
      options.start_state_path,
      start_state,
      {.time_index = options.start_time_index, .variables = state_variables});
  tywrf::io::load_wrf_state(
      options.end_state_path,
      end_state,
      {.time_index = options.end_time_index, .variables = state_variables});
  compute_state_tendencies(start_state, end_state, tendency, dt_seconds);

  tywrf::State<float> d01_dummy_state(grid);
  const tywrf::dynamics::SkeletonStateStepper stepper(loop_config);
  const auto step_report = stepper.run_with_explicit_tendencies(
      d01_dummy_state,
      start_state,
      tywrf::dynamics::ExplicitStateTendencySet{nullptr, &tendency});
  if (step_report.status != tywrf::dynamics::StateStepStatus::ok ||
      step_report.layout_or_status_failure) {
    throw std::runtime_error("SkeletonStateStepper explicit-tendency oracle path failed");
  }
  if (step_report.d02.explicit_tendency_apply_count != step_report.loop.child_steps ||
      step_report.d02.tendency_dt_seconds != dt_seconds) {
    throw std::runtime_error(
        "d02 explicit tendency did not cover the full 600 s validation segment");
  }

  tywrf::io::write_wrf_state(
      options.output_path,
      start_state,
      {
          .time_index = options.output_time_index,
          .variables = options.variables,
          .template_path = options.template_path,
          .template_time_index = options.template_time_index,
          .times_value = options.times_value,
      });

  stamp_oracle_metadata(options, template_resolution, state_variables, step_report, dt_seconds);
  print_report(options, template_resolution, state_variables, step_report, dt_seconds);
  return 0;
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    const auto options = parse_options(argc, argv);
    return run(options);
  } catch (const std::exception& error) {
    std::cerr << "tywrf_tendency_cycle: " << error.what() << "\n\n" << usage();
    return 2;
  }
}
