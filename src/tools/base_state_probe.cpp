#include "tywrf/dynamics/base_state.hpp"
#include "tywrf/grid.hpp"
#include "tywrf/state.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

struct CliOptions {
  std::filesystem::path input;
  std::string domain;
  std::int32_t mass_nx = 0;
  std::int32_t mass_ny = 0;
  std::int32_t mass_nz = 0;
  std::int32_t full_nz = 0;
  std::int32_t halo = 0;
  std::size_t time_index = 0;
  bool pretty = false;
};

struct NetcdfVariable {
  int id = -1;
  nc_type type = NC_NAT;
  std::vector<std::string> dimensions;
  std::vector<std::size_t> shape;
};

struct DiffStats {
  std::string status = "not_available";
  std::string detail;
  std::int32_t nx = 0;
  std::int32_t ny = 0;
  std::int32_t nz = 0;
  double max_abs_diff = 0.0;
  double mean_abs_diff = 0.0;
  std::uint64_t point_count = 0;
};

class NetcdfReadHandle {
 public:
  explicit NetcdfReadHandle(std::filesystem::path path) : path_(std::move(path)) {
    check(nc_open(path_.string().c_str(), NC_NOWRITE, &id_), "open");
  }

  NetcdfReadHandle(const NetcdfReadHandle&) = delete;
  NetcdfReadHandle& operator=(const NetcdfReadHandle&) = delete;

  ~NetcdfReadHandle() {
    if (id_ >= 0) {
      nc_close(id_);
    }
  }

  [[nodiscard]] int id() const noexcept {
    return id_;
  }

  [[nodiscard]] const std::filesystem::path& path() const noexcept {
    return path_;
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
  return "usage: tywrf_base_state_probe --input PATH --domain d01|d02 "
         "--mass-nx N --mass-ny N --mass-nz N --full-nz N "
         "[--halo N] [--time-index N] [--pretty]";
}

[[nodiscard]] std::int32_t parse_i32(
    const std::string_view value,
    const std::string_view name) {
  std::size_t parsed = 0;
  int parsed_value = 0;
  try {
    parsed_value = std::stoi(std::string(value), &parsed);
  } catch (const std::exception&) {
    throw std::runtime_error("invalid integer for " + std::string(name) + ": " +
                             std::string(value));
  }
  if (parsed != value.size()) {
    throw std::runtime_error("invalid integer for " + std::string(name) + ": " +
                             std::string(value));
  }
  return static_cast<std::int32_t>(parsed_value);
}

[[nodiscard]] std::size_t parse_size(
    const std::string_view value,
    const std::string_view name) {
  std::size_t parsed = 0;
  unsigned long long parsed_value = 0;
  try {
    parsed_value = std::stoull(std::string(value), &parsed);
  } catch (const std::exception&) {
    throw std::runtime_error("invalid integer for " + std::string(name) + ": " +
                             std::string(value));
  }
  if (parsed != value.size()) {
    throw std::runtime_error("invalid integer for " + std::string(name) + ": " +
                             std::string(value));
  }
  return static_cast<std::size_t>(parsed_value);
}

[[nodiscard]] CliOptions parse_args(const int argc, char** argv) {
  CliOptions options;
  for (int index = 1; index < argc; ++index) {
    const std::string_view arg(argv[index]);
    const auto require_value = [&](const std::string_view name) -> std::string_view {
      if (index + 1 >= argc) {
        throw std::runtime_error("missing value for " + std::string(name));
      }
      ++index;
      return argv[index];
    };

    if (arg == "--input") {
      options.input = std::filesystem::path(require_value(arg));
    } else if (arg == "--domain") {
      options.domain = std::string(require_value(arg));
    } else if (arg == "--mass-nx") {
      options.mass_nx = parse_i32(require_value(arg), arg);
    } else if (arg == "--mass-ny") {
      options.mass_ny = parse_i32(require_value(arg), arg);
    } else if (arg == "--mass-nz") {
      options.mass_nz = parse_i32(require_value(arg), arg);
    } else if (arg == "--full-nz") {
      options.full_nz = parse_i32(require_value(arg), arg);
    } else if (arg == "--halo") {
      options.halo = parse_i32(require_value(arg), arg);
    } else if (arg == "--time-index") {
      options.time_index = parse_size(require_value(arg), arg);
    } else if (arg == "--pretty") {
      options.pretty = true;
    } else if (arg == "--help" || arg == "-h") {
      std::cout << usage() << '\n';
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + std::string(arg));
    }
  }

  if (options.input.empty()) {
    throw std::runtime_error("missing required --input");
  }
  if (options.domain != "d01" && options.domain != "d02") {
    throw std::runtime_error("--domain must be d01 or d02");
  }
  if (options.mass_nx <= 0 || options.mass_ny <= 0 || options.mass_nz <= 0 ||
      options.full_nz <= 0) {
    throw std::runtime_error("--mass-nx, --mass-ny, --mass-nz, and --full-nz must be positive");
  }
  if (options.full_nz != options.mass_nz + 1) {
    throw std::runtime_error("--full-nz must equal --mass-nz + 1 for KROSA base state");
  }
  if (options.halo < 0) {
    throw std::runtime_error("--halo must be non-negative");
  }
  return options;
}

[[nodiscard]] bool numeric_type(const nc_type type) noexcept {
  switch (type) {
    case NC_BYTE:
    case NC_SHORT:
    case NC_INT:
    case NC_FLOAT:
    case NC_DOUBLE:
    case NC_UBYTE:
    case NC_USHORT:
    case NC_UINT:
    case NC_INT64:
    case NC_UINT64:
      return true;
    default:
      return false;
  }
}

[[nodiscard]] bool float_type(const nc_type type) noexcept {
  return type == NC_FLOAT || type == NC_DOUBLE;
}

[[nodiscard]] std::string join_strings(const std::vector<std::string>& values) {
  std::ostringstream joined;
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) {
      joined << ",";
    }
    joined << values[index];
  }
  return joined.str();
}

[[nodiscard]] std::string join_shape(const std::vector<std::size_t>& values) {
  std::ostringstream joined;
  for (const auto value : values) {
    joined << ' ' << value;
  }
  return joined.str();
}

[[nodiscard]] bool inquire_variable_if_present(
    const NetcdfReadHandle& file,
    const std::string_view name,
    NetcdfVariable& variable) {
  const int status = nc_inq_varid(file.id(), std::string(name).c_str(), &variable.id);
  if (status == NC_ENOTVAR) {
    return false;
  }
  file.check(status, "inquire variable id");

  int ndims = 0;
  file.check(nc_inq_varndims(file.id(), variable.id, &ndims), "inquire variable rank");
  std::vector<int> dimids(static_cast<std::size_t>(ndims));
  std::array<char, NC_MAX_NAME + 1> name_buffer{};
  int natts = 0;
  file.check(
      nc_inq_var(
          file.id(),
          variable.id,
          name_buffer.data(),
          &variable.type,
          &ndims,
          dimids.empty() ? nullptr : dimids.data(),
          &natts),
      "inquire variable");

  variable.dimensions.clear();
  variable.shape.clear();
  for (const int dimid : dimids) {
    std::array<char, NC_MAX_NAME + 1> dim_name{};
    std::size_t dim_size = 0;
    file.check(nc_inq_dim(file.id(), dimid, dim_name.data(), &dim_size), "inquire dimension");
    variable.dimensions.push_back(dim_name.data());
    variable.shape.push_back(dim_size);
  }
  return true;
}

[[nodiscard]] NetcdfVariable require_variable(
    const NetcdfReadHandle& file,
    const std::string_view name) {
  NetcdfVariable variable;
  if (!inquire_variable_if_present(file, name, variable)) {
    throw std::runtime_error("missing base-state variable " + std::string(name) +
                             " in " + file.path().string());
  }
  return variable;
}

[[nodiscard]] bool has_leading_time_dimension(const NetcdfVariable& variable) {
  return !variable.dimensions.empty() && variable.dimensions.front() == "Time";
}

[[nodiscard]] std::string describe_shape(const NetcdfVariable& variable) {
  std::ostringstream message;
  message << "dimensions " << join_strings(variable.dimensions) << " shape"
          << join_shape(variable.shape);
  return message.str();
}

void require_float_variable(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const NetcdfVariable& variable) {
  if (!float_type(variable.type)) {
    throw std::runtime_error("base-state variable " + std::string(name) + " in " +
                             file.path().string() + " must be float or double");
  }
}

[[nodiscard]] bool valid_vector_shape(
    const NetcdfVariable& variable,
    const std::size_t expected_count,
    const std::size_t time_index) {
  const bool vector_shape = variable.shape.size() == 1 &&
                            variable.shape.front() == expected_count;
  const bool time_vector_shape =
      variable.shape.size() == 2 && has_leading_time_dimension(variable) &&
      variable.shape.front() > time_index && variable.shape.back() == expected_count;
  return vector_shape || time_vector_shape;
}

[[nodiscard]] std::vector<float> read_vector_variable(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const std::size_t expected_count,
    const std::size_t time_index) {
  const auto variable = require_variable(file, name);
  require_float_variable(file, name, variable);
  if (!valid_vector_shape(variable, expected_count, time_index)) {
    std::ostringstream message;
    message << "base-state variable " << name << " in " << file.path()
            << " has " << describe_shape(variable) << ", expected count "
            << expected_count << " with optional leading Time";
    throw std::runtime_error(message.str());
  }

  std::vector<float> values(expected_count, 0.0F);
  std::vector<std::size_t> start(variable.shape.size(), 0);
  std::vector<std::size_t> count = variable.shape;
  if (has_leading_time_dimension(variable)) {
    start.front() = time_index;
    count.front() = 1;
  }
  file.check(
      nc_get_vara_float(file.id(), variable.id, start.data(), count.data(), values.data()),
      "read vertical coefficient");
  return values;
}

[[nodiscard]] bool has_global_attribute(
    const NetcdfReadHandle& file,
    const std::string_view name) {
  const int status = nc_inq_att(file.id(), NC_GLOBAL, std::string(name).c_str(), nullptr, nullptr);
  if (status == NC_NOERR) {
    return true;
  }
  if (status == NC_ENOTATT) {
    return false;
  }
  file.check(status, "inquire global attribute");
  return false;
}

[[nodiscard]] float checked_p_top(
    const double value,
    const std::filesystem::path& path) {
  if (!std::isfinite(value) || value < 0.0 ||
      value > static_cast<double>(std::numeric_limits<float>::max())) {
    std::ostringstream message;
    message << "P_TOP in " << path << " must be finite and non-negative, found "
            << value;
    throw std::runtime_error(message.str());
  }
  return static_cast<float>(value);
}

[[nodiscard]] float read_p_top(
    const NetcdfReadHandle& file,
    const std::size_t time_index,
    std::string& source) {
  double value = 0.0;
  if (has_global_attribute(file, "P_TOP")) {
    nc_type type = NC_NAT;
    std::size_t length = 0;
    file.check(nc_inq_att(file.id(), NC_GLOBAL, "P_TOP", &type, &length), "inquire P_TOP attr");
    if (length != 1 || !numeric_type(type)) {
      throw std::runtime_error("global attribute P_TOP must be a scalar numeric value");
    }
    file.check(nc_get_att_double(file.id(), NC_GLOBAL, "P_TOP", &value), "read P_TOP attr");
    source = "global_attribute";
    return checked_p_top(value, file.path());
  }

  const auto variable = require_variable(file, "P_TOP");
  if (!numeric_type(variable.type)) {
    throw std::runtime_error("P_TOP must be numeric");
  }
  const bool scalar = variable.shape.empty();
  const bool time_scalar = variable.shape.size() == 1 &&
                           has_leading_time_dimension(variable) &&
                           variable.shape.front() > time_index;
  if (!scalar && !time_scalar) {
    throw std::runtime_error("P_TOP must be scalar or Time scalar");
  }
  if (scalar) {
    file.check(nc_get_var_double(file.id(), variable.id, &value), "read scalar P_TOP");
    source = "scalar_variable";
  } else {
    const std::size_t start[1] = {time_index};
    const std::size_t count[1] = {1};
    file.check(nc_get_vara_double(file.id(), variable.id, start, count, &value), "read P_TOP");
    source = "time_variable";
  }
  return checked_p_top(value, file.path());
}

[[nodiscard]] bool valid_2d_shape(
    const NetcdfVariable& variable,
    const std::int32_t nx,
    const std::int32_t ny,
    const std::size_t time_index) {
  const auto nx_size = static_cast<std::size_t>(nx);
  const auto ny_size = static_cast<std::size_t>(ny);
  const bool plain = variable.shape.size() == 2 && variable.shape[0] == ny_size &&
                     variable.shape[1] == nx_size;
  const bool timed = variable.shape.size() == 3 && has_leading_time_dimension(variable) &&
                     variable.shape[0] > time_index && variable.shape[1] == ny_size &&
                     variable.shape[2] == nx_size;
  return plain || timed;
}

[[nodiscard]] bool valid_3d_shape(
    const NetcdfVariable& variable,
    const std::int32_t nx,
    const std::int32_t ny,
    const std::int32_t nz,
    const std::size_t time_index) {
  const auto nx_size = static_cast<std::size_t>(nx);
  const auto ny_size = static_cast<std::size_t>(ny);
  const auto nz_size = static_cast<std::size_t>(nz);
  const bool plain = variable.shape.size() == 3 && variable.shape[0] == nz_size &&
                     variable.shape[1] == ny_size && variable.shape[2] == nx_size;
  const bool timed = variable.shape.size() == 4 && has_leading_time_dimension(variable) &&
                     variable.shape[0] > time_index && variable.shape[1] == nz_size &&
                     variable.shape[2] == ny_size && variable.shape[3] == nx_size;
  return plain || timed;
}

void read_2d_into_view(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const NetcdfVariable& variable,
    const std::size_t time_index,
    tywrf::FieldView2D<float> view) {
  const auto active_nx = view.nx - view.halo.i_lower - view.halo.i_upper;
  const auto active_ny = view.ny - view.halo.j_lower - view.halo.j_upper;
  std::vector<float> raw(
      static_cast<std::size_t>(active_nx) * static_cast<std::size_t>(active_ny),
      0.0F);
  std::vector<std::size_t> start(variable.shape.size(), 0);
  std::vector<std::size_t> count = variable.shape;
  if (has_leading_time_dimension(variable)) {
    start.front() = time_index;
    count.front() = 1;
  }
  file.check(
      nc_get_vara_float(file.id(), variable.id, start.data(), count.data(), raw.data()),
      "read 2D field");

  for (std::int32_t j = 0; j < active_ny; ++j) {
    for (std::int32_t i = 0; i < active_nx; ++i) {
      const auto raw_index =
          static_cast<std::size_t>(j) * static_cast<std::size_t>(active_nx) +
          static_cast<std::size_t>(i);
      view(view.halo.i_lower + i, view.halo.j_lower + j) = raw[raw_index];
    }
  }
}

void read_3d_into_view(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const NetcdfVariable& variable,
    const std::size_t time_index,
    tywrf::FieldView3D<float> view) {
  const auto active_nx = view.nx - view.halo.i_lower - view.halo.i_upper;
  const auto active_ny = view.ny - view.halo.j_lower - view.halo.j_upper;
  const auto active_nz = view.nz - view.halo.k_lower - view.halo.k_upper;
  std::vector<float> raw(
      static_cast<std::size_t>(active_nx) * static_cast<std::size_t>(active_ny) *
          static_cast<std::size_t>(active_nz),
      0.0F);
  std::vector<std::size_t> start(variable.shape.size(), 0);
  std::vector<std::size_t> count = variable.shape;
  if (has_leading_time_dimension(variable)) {
    start.front() = time_index;
    count.front() = 1;
  }
  file.check(
      nc_get_vara_float(file.id(), variable.id, start.data(), count.data(), raw.data()),
      "read 3D field");

  for (std::int32_t k = 0; k < active_nz; ++k) {
    for (std::int32_t j = 0; j < active_ny; ++j) {
      for (std::int32_t i = 0; i < active_nx; ++i) {
        const auto raw_index =
            (static_cast<std::size_t>(k) * static_cast<std::size_t>(active_ny) +
             static_cast<std::size_t>(j)) *
                static_cast<std::size_t>(active_nx) +
            static_cast<std::size_t>(i);
        view(view.halo.i_lower + i, view.halo.j_lower + j, view.halo.k_lower + k) =
            raw[raw_index];
      }
    }
  }
}

void read_required_2d(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const std::size_t time_index,
    tywrf::FieldView2D<float> view) {
  const auto variable = require_variable(file, name);
  require_float_variable(file, name, variable);
  const auto active_nx = view.nx - view.halo.i_lower - view.halo.i_upper;
  const auto active_ny = view.ny - view.halo.j_lower - view.halo.j_upper;
  if (!valid_2d_shape(variable, active_nx, active_ny, time_index)) {
    std::ostringstream message;
    message << "base-state variable " << name << " in " << file.path() << " has "
            << describe_shape(variable) << ", expected surface shape ["
            << active_ny << ',' << active_nx << "] with optional leading Time";
    throw std::runtime_error(message.str());
  }
  read_2d_into_view(file, name, variable, time_index, view);
}

[[nodiscard]] tywrf::Grid make_grid(const CliOptions& options) {
  return tywrf::Grid(
      {options.mass_nx,
       options.mass_ny,
       options.mass_nz,
       options.full_nz,
       tywrf::uniform_halo_3d(options.halo)});
}

template <typename View>
[[nodiscard]] DiffStats compare_view(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const std::size_t time_index,
    View reconstructed) {
  DiffStats stats;
  NetcdfVariable variable;
  if (!inquire_variable_if_present(file, name, variable)) {
    stats.status = "not_available";
    stats.detail = "missing";
    return stats;
  }
  if (!float_type(variable.type)) {
    stats.status = "shape_mismatch";
    stats.detail = "not a float or double variable";
    return stats;
  }

  const auto active_nx =
      reconstructed.nx - reconstructed.halo.i_lower - reconstructed.halo.i_upper;
  const auto active_ny =
      reconstructed.ny - reconstructed.halo.j_lower - reconstructed.halo.j_upper;
  if constexpr (requires { reconstructed.nz; }) {
    const auto active_nz =
        reconstructed.nz - reconstructed.halo.k_lower - reconstructed.halo.k_upper;
    stats.nx = active_nx;
    stats.ny = active_ny;
    stats.nz = active_nz;
    if (!valid_3d_shape(variable, active_nx, active_ny, active_nz, time_index)) {
      stats.status = "shape_mismatch";
      stats.detail = describe_shape(variable);
      return stats;
    }
    tywrf::FieldStorage3D<float> direct(
        tywrf::make_field_layout(
            tywrf::ActiveShape3D{active_nx, active_ny, active_nz},
            reconstructed.halo));
    read_3d_into_view(file, name, variable, time_index, direct.view());
    const auto direct_view = direct.view();
    double sum = 0.0;
    double max_value = 0.0;
    std::uint64_t count = 0;
    for (std::int32_t k = 0; k < active_nz; ++k) {
      for (std::int32_t j = 0; j < active_ny; ++j) {
        for (std::int32_t i = 0; i < active_nx; ++i) {
          const auto ii = reconstructed.halo.i_lower + i;
          const auto jj = reconstructed.halo.j_lower + j;
          const auto kk = reconstructed.halo.k_lower + k;
          const auto diff = std::abs(
              static_cast<double>(reconstructed(ii, jj, kk)) -
              static_cast<double>(direct_view(ii, jj, kk)));
          max_value = std::max(max_value, diff);
          sum += diff;
          ++count;
        }
      }
    }
    stats.status = "compared";
    stats.max_abs_diff = max_value;
    stats.mean_abs_diff = count == 0 ? 0.0 : sum / static_cast<double>(count);
    stats.point_count = count;
  } else {
    stats.nx = active_nx;
    stats.ny = active_ny;
    stats.nz = 1;
    if (!valid_2d_shape(variable, active_nx, active_ny, time_index)) {
      stats.status = "shape_mismatch";
      stats.detail = describe_shape(variable);
      return stats;
    }
    tywrf::FieldStorage2D<float> direct(
        tywrf::make_field_layout(
            tywrf::ActiveShape2D{active_nx, active_ny},
            reconstructed.halo));
    read_2d_into_view(file, name, variable, time_index, direct.view());
    const auto direct_view = direct.view();
    double sum = 0.0;
    double max_value = 0.0;
    std::uint64_t count = 0;
    for (std::int32_t j = 0; j < active_ny; ++j) {
      for (std::int32_t i = 0; i < active_nx; ++i) {
        const auto ii = reconstructed.halo.i_lower + i;
        const auto jj = reconstructed.halo.j_lower + j;
        const auto diff =
            std::abs(static_cast<double>(reconstructed(ii, jj)) -
                     static_cast<double>(direct_view(ii, jj)));
        max_value = std::max(max_value, diff);
        sum += diff;
        ++count;
      }
    }
    stats.status = "compared";
    stats.max_abs_diff = max_value;
    stats.mean_abs_diff = count == 0 ? 0.0 : sum / static_cast<double>(count);
    stats.point_count = count;
  }
  return stats;
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

void write_indent(std::ostream& output, const bool pretty, const int level) {
  if (!pretty) {
    return;
  }
  output << '\n';
  for (int index = 0; index < level; ++index) {
    output << "  ";
  }
}

void write_key(std::ostream& output, const std::string_view key, const bool pretty) {
  output << '"' << key << "\":";
  if (pretty) {
    output << ' ';
  }
}

template <typename Writer>
void write_field(
    std::ostream& output,
    const std::string_view key,
    const bool pretty,
    bool& first,
    Writer writer) {
  if (!first) {
    output << ',';
  }
  first = false;
  write_indent(output, pretty, 1);
  write_key(output, key, pretty);
  writer();
}

template <typename Writer>
void write_nested_field(
    std::ostream& output,
    const std::string_view key,
    const bool pretty,
    bool& first,
    const int level,
    Writer writer) {
  if (!first) {
    output << ',';
  }
  first = false;
  write_indent(output, pretty, level);
  write_key(output, key, pretty);
  writer();
}

void write_shape2(std::ostream& output, const std::int32_t nx, const std::int32_t ny) {
  output << '[' << nx << ',' << ny << ']';
}

void write_shape3(
    std::ostream& output,
    const std::int32_t nx,
    const std::int32_t ny,
    const std::int32_t nz) {
  output << '[' << nx << ',' << ny << ',' << nz << ']';
}

void write_diff_stats(
    std::ostream& output,
    const DiffStats& stats,
    const bool pretty,
    const int level) {
  bool first = true;
  output << '{';
  write_nested_field(output, "status", pretty, first, level + 1, [&] {
    output << '"' << stats.status << '"';
  });
  write_nested_field(output, "shape", pretty, first, level + 1, [&] {
    if (stats.nz == 1) {
      write_shape2(output, stats.nx, stats.ny);
    } else {
      write_shape3(output, stats.nx, stats.ny, stats.nz);
    }
  });
  if (!stats.detail.empty()) {
    write_nested_field(output, "detail", pretty, first, level + 1, [&] {
      output << '"' << json_escape(stats.detail) << '"';
    });
  }
  write_nested_field(output, "point_count", pretty, first, level + 1, [&] {
    output << stats.point_count;
  });
  if (stats.status == "compared") {
    write_nested_field(output, "max_abs_diff", pretty, first, level + 1, [&] {
      output << stats.max_abs_diff;
    });
    write_nested_field(output, "mean_abs_diff", pretty, first, level + 1, [&] {
      output << stats.mean_abs_diff;
    });
  }
  write_indent(output, pretty, level);
  output << '}';
}

void write_error_json(
    std::ostream& output,
    const CliOptions& options,
    const std::string_view message,
    const bool pretty) {
  bool first = true;
  output << std::setprecision(std::numeric_limits<double>::max_digits10);
  output << '{';
  write_field(output, "status", pretty, first, [&] { output << "\"error\""; });
  write_field(output, "input_file", pretty, first, [&] {
    output << '"' << json_escape(options.input.string()) << '"';
  });
  write_field(output, "domain", pretty, first, [&] {
    output << '"' << json_escape(options.domain) << '"';
  });
  write_field(output, "time_index", pretty, first, [&] { output << options.time_index; });
  write_field(output, "error", pretty, first, [&] {
    output << '"' << json_escape(message) << '"';
  });
  write_field(output, "diagnostic_only", pretty, first, [&] { output << "true"; });
  write_field(output, "gate_candidate", pretty, first, [&] { output << "false"; });
  write_field(output, "integrator_output", pretty, first, [&] { output << "false"; });
  write_field(output, "writes_netcdf", pretty, first, [&] { output << "false"; });
  write_field(output, "modifies_state", pretty, first, [&] { output << "false"; });
  write_field(output, "calls_pressure_refresh_compute", pretty, first, [&] {
    output << "false";
  });
  write_field(output, "uses_later_restart_truth", pretty, first, [&] {
    output << "false";
  });
  write_indent(output, pretty, 0);
  output << '}';
  if (pretty) {
    output << '\n';
  }
}

void write_report_json(
    std::ostream& output,
    const CliOptions& options,
    const float p_top,
    const std::string_view p_top_source,
    const tywrf::dynamics::KrosaMassBaseStateReconstructionReport& report,
    const DiffStats& pb,
    const DiffStats& t_init,
    const DiffStats& alb,
    const DiffStats& mub,
    const DiffStats& phb) {
  bool first = true;
  output << std::setprecision(std::numeric_limits<double>::max_digits10);
  output << '{';
  write_field(output, "status", options.pretty, first, [&] {
    output << '"' << (report.ok() ? "ok" : "error") << '"';
  });
  write_field(output, "input_file", options.pretty, first, [&] {
    output << '"' << json_escape(options.input.string()) << '"';
  });
  write_field(output, "domain", options.pretty, first, [&] {
    output << '"' << options.domain << '"';
  });
  write_field(output, "time_index", options.pretty, first, [&] {
    output << options.time_index;
  });
  write_field(output, "p_top_pa", options.pretty, first, [&] { output << p_top; });
  write_field(output, "p_top_source", options.pretty, first, [&] {
    output << '"' << json_escape(p_top_source) << '"';
  });
  write_field(output, "expected_shapes", options.pretty, first, [&] {
    bool shapes_first = true;
    output << '{';
    write_nested_field(output, "surface", options.pretty, shapes_first, 2, [&] {
      write_shape2(output, options.mass_nx, options.mass_ny);
    });
    write_nested_field(output, "mass", options.pretty, shapes_first, 2, [&] {
      write_shape3(output, options.mass_nx, options.mass_ny, options.mass_nz);
    });
    write_nested_field(output, "full", options.pretty, shapes_first, 2, [&] {
      write_shape3(output, options.mass_nx, options.mass_ny, options.full_nz);
    });
    write_indent(output, options.pretty, 1);
    output << '}';
  });
  write_field(output, "reconstruction", options.pretty, first, [&] {
    bool reconstruction_first = true;
    output << '{';
    write_nested_field(output, "ok", options.pretty, reconstruction_first, 2, [&] {
      output << (report.ok() ? "true" : "false");
    });
    write_nested_field(output, "message", options.pretty, reconstruction_first, 2, [&] {
      output << '"' << json_escape(report.result.message) << '"';
    });
    write_nested_field(output, "active_nx", options.pretty, reconstruction_first, 2, [&] {
      output << report.active_nx;
    });
    write_nested_field(output, "active_ny", options.pretty, reconstruction_first, 2, [&] {
      output << report.active_ny;
    });
    write_nested_field(output, "active_nz", options.pretty, reconstruction_first, 2, [&] {
      output << report.active_nz;
    });
    write_nested_field(output, "reconstructed_columns", options.pretty, reconstruction_first, 2, [&] {
      output << report.reconstructed_column_count;
    });
    write_nested_field(output, "reconstructed_points", options.pretty, reconstruction_first, 2, [&] {
      output << report.reconstructed_point_count;
    });
    write_nested_field(output, "invalid_columns", options.pretty, reconstruction_first, 2, [&] {
      output << report.invalid_column_count;
    });
    write_nested_field(output, "invalid_points", options.pretty, reconstruction_first, 2, [&] {
      output << report.invalid_point_count;
    });
    write_nested_field(output, "wrote_pb", options.pretty, reconstruction_first, 2, [&] {
      output << (report.wrote_pb ? "true" : "false");
    });
    write_nested_field(output, "wrote_t_init", options.pretty, reconstruction_first, 2, [&] {
      output << (report.wrote_t_init ? "true" : "false");
    });
    write_nested_field(output, "wrote_mub", options.pretty, reconstruction_first, 2, [&] {
      output << (report.wrote_mub ? "true" : "false");
    });
    write_nested_field(output, "wrote_alb", options.pretty, reconstruction_first, 2, [&] {
      output << (report.wrote_alb ? "true" : "false");
    });
    write_nested_field(output, "wrote_phb", options.pretty, reconstruction_first, 2, [&] {
      output << (report.wrote_phb ? "true" : "false");
    });
    write_nested_field(
        output,
        "phb_full_level_reconstruction_implemented",
        options.pretty,
        reconstruction_first,
        2,
        [&] {
          output << (report.phb_full_level_reconstruction_implemented ? "true" : "false");
        });
    write_indent(output, options.pretty, 1);
    output << '}';
  });
  write_field(output, "direct_fields", options.pretty, first, [&] {
    bool direct_first = true;
    output << '{';
    write_nested_field(output, "PB", options.pretty, direct_first, 2, [&] {
      write_diff_stats(output, pb, options.pretty, 2);
    });
    write_nested_field(output, "T_INIT", options.pretty, direct_first, 2, [&] {
      write_diff_stats(output, t_init, options.pretty, 2);
    });
    write_nested_field(output, "ALB", options.pretty, direct_first, 2, [&] {
      write_diff_stats(output, alb, options.pretty, 2);
    });
    write_nested_field(output, "MUB", options.pretty, direct_first, 2, [&] {
      write_diff_stats(output, mub, options.pretty, 2);
    });
    write_nested_field(output, "PHB", options.pretty, direct_first, 2, [&] {
      write_diff_stats(output, phb, options.pretty, 2);
    });
    write_indent(output, options.pretty, 1);
    output << '}';
  });
  write_field(output, "diagnostic_only", options.pretty, first, [&] { output << "true"; });
  write_field(output, "gate_candidate", options.pretty, first, [&] { output << "false"; });
  write_field(output, "integrator_output", options.pretty, first, [&] { output << "false"; });
  write_field(output, "writes_netcdf", options.pretty, first, [&] { output << "false"; });
  write_field(output, "modifies_state", options.pretty, first, [&] { output << "false"; });
  write_field(output, "calls_pressure_refresh_compute", options.pretty, first, [&] {
    output << "false";
  });
  write_field(output, "uses_later_restart_truth", options.pretty, first, [&] {
    output << "false";
  });
  write_indent(output, options.pretty, 0);
  output << '}';
  if (options.pretty) {
    output << '\n';
  }
}

}  // namespace

int main(const int argc, char** argv) {
  CliOptions options;
  try {
    options = parse_args(argc, argv);
    const auto grid = make_grid(options);
    const NetcdfReadHandle file(options.input);

    std::string p_top_source;
    const auto p_top = read_p_top(file, options.time_index, p_top_source);
    const auto c3f = read_vector_variable(
        file,
        "C3F",
        static_cast<std::size_t>(options.full_nz),
        options.time_index);
    const auto c4f = read_vector_variable(
        file,
        "C4F",
        static_cast<std::size_t>(options.full_nz),
        options.time_index);
    const auto c3h = read_vector_variable(
        file,
        "C3H",
        static_cast<std::size_t>(options.mass_nz),
        options.time_index);
    const auto c4h = read_vector_variable(
        file,
        "C4H",
        static_cast<std::size_t>(options.mass_nz),
        options.time_index);

    tywrf::FieldStorage2D<float> hgt(grid.surface_layout());
    tywrf::FieldStorage3D<float> pb(grid.mass_layout());
    tywrf::FieldStorage3D<float> t_init(grid.mass_layout());
    tywrf::FieldStorage2D<float> mub(grid.surface_layout());
    tywrf::FieldStorage3D<float> alb(grid.mass_layout());
    tywrf::FieldStorage3D<float> phb(grid.w_layout());

    read_required_2d(file, "HGT", options.time_index, hgt.view());

    const auto hgt_view = hgt.view();
    const tywrf::dynamics::KrosaMassBaseStateReconstructionInputs inputs{
        {hgt_view.data,
         hgt_view.nx,
         hgt_view.ny,
         hgt_view.stride_i,
         hgt_view.stride_j,
         hgt_view.halo},
        {c3h.data(), static_cast<std::int32_t>(c3h.size())},
        {c4h.data(), static_cast<std::int32_t>(c4h.size())},
        {c3f.data(), static_cast<std::int32_t>(c3f.size())},
        {c4f.data(), static_cast<std::int32_t>(c4f.size())},
        p_top};
    const tywrf::dynamics::KrosaMassBaseStateReconstructionOutputs outputs{
        pb.view(), t_init.view(), mub.view(), alb.view(), phb.view()};
    const auto report =
        tywrf::dynamics::reconstruct_krosa_mass_base_state(inputs, outputs);

    const auto pb_diff = compare_view(file, "PB", options.time_index, pb.view());
    const auto t_init_diff = compare_view(file, "T_INIT", options.time_index, t_init.view());
    const auto alb_diff = compare_view(file, "ALB", options.time_index, alb.view());
    const auto mub_diff = compare_view(file, "MUB", options.time_index, mub.view());
    const auto phb_diff = compare_view(file, "PHB", options.time_index, phb.view());

    write_report_json(
        std::cout,
        options,
        p_top,
        p_top_source,
        report,
        pb_diff,
        t_init_diff,
        alb_diff,
        mub_diff,
        phb_diff);
    return report.ok() ? 0 : 1;
  } catch (const std::exception& error) {
    if (options.input.empty()) {
      std::cerr << usage() << '\n';
      std::cerr << "base-state probe error: " << error.what() << '\n';
      return 2;
    }
    write_error_json(std::cout, options, error.what(), options.pretty);
    std::cerr << "base-state probe error: " << error.what() << '\n';
    return 1;
  }
}
