#include "tywrf/io/pressure_refresh_io.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <string_view>
#include <utility>

namespace tywrf::io {
namespace {

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
    throw PressureRefreshIoError(message.str());
  }

 private:
  std::filesystem::path path_;
  int id_ = -1;
};

struct NetcdfVariable {
  int id = -1;
  nc_type type = NC_NAT;
  std::vector<std::string> dimensions;
  std::vector<std::size_t> shape;
};

[[nodiscard]] const std::vector<std::string>& required_names_ref() {
  static const std::vector<std::string> names = {
      "P_TOP", "C3F", "C4F", "C3H", "C4H", "ALB"};
  return names;
}

[[nodiscard]] std::string nc_type_name(const nc_type type) {
  switch (type) {
    case NC_BYTE:
      return "byte";
    case NC_CHAR:
      return "char";
    case NC_SHORT:
      return "short";
    case NC_INT:
      return "int";
    case NC_FLOAT:
      return "float";
    case NC_DOUBLE:
      return "double";
    case NC_UBYTE:
      return "ubyte";
    case NC_USHORT:
      return "ushort";
    case NC_UINT:
      return "uint";
    case NC_INT64:
      return "int64";
    case NC_UINT64:
      return "uint64";
    default:
      return "unknown";
  }
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

void append_missing_unique(
    std::vector<std::string>& names,
    const std::string_view name) {
  if (std::find(names.begin(), names.end(), name) == names.end()) {
    names.emplace_back(name);
  }
}

[[nodiscard]] std::size_t checked_point_count(
    const std::int32_t nx,
    const std::int32_t ny,
    const std::int32_t nz) {
  if (nx <= 0 || ny <= 0 || nz <= 0) {
    throw PressureRefreshIoError("pressure-refresh grid dimensions must be positive");
  }
  return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
         static_cast<std::size_t>(nz);
}

[[nodiscard]] std::size_t checked_point_count(
    const std::int32_t nx,
    const std::int32_t ny) {
  if (nx <= 0 || ny <= 0) {
    throw PressureRefreshIoError("pressure-refresh grid dimensions must be positive");
  }
  return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
}

void require_grid(const Grid& grid) {
  const auto& config = grid.config();
  if (config.mass_nx <= 0 || config.mass_ny <= 0 || config.mass_nz <= 0 ||
      config.full_nz <= 0) {
    throw PressureRefreshIoError(
        "cannot read pressure-refresh inputs for invalid grid dimensions");
  }
  if (config.full_nz != config.mass_nz + 1) {
    std::ostringstream message;
    message << "cannot read pressure-refresh inputs: full_nz is " << config.full_nz
            << ", expected mass_nz + 1 = " << (config.mass_nz + 1);
    throw PressureRefreshIoError(message.str());
  }
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

[[nodiscard]] bool inquire_variable_if_present(
    const NetcdfReadHandle& file,
    const std::string_view name,
    NetcdfVariable& variable) {
  const int var_status = nc_inq_varid(file.id(), std::string(name).c_str(), &variable.id);
  if (var_status == NC_ENOTVAR) {
    return false;
  }
  file.check(var_status, "inquire variable id");

  int ndims = 0;
  file.check(nc_inq_varndims(file.id(), variable.id, &ndims), "inquire variable rank");
  std::vector<int> dimids(static_cast<std::size_t>(ndims));
  std::array<char, NC_MAX_NAME + 1> variable_name{};
  int natts = 0;
  file.check(
      nc_inq_var(
          file.id(),
          variable.id,
          variable_name.data(),
          &variable.type,
          &ndims,
          dimids.empty() ? nullptr : dimids.data(),
          &natts),
      "inquire variable");

  variable.dimensions.clear();
  variable.shape.clear();
  variable.dimensions.reserve(dimids.size());
  variable.shape.reserve(dimids.size());
  for (const int dimid : dimids) {
    std::array<char, NC_MAX_NAME + 1> dimension_name{};
    std::size_t dimension_size = 0;
    file.check(
        nc_inq_dim(file.id(), dimid, dimension_name.data(), &dimension_size),
        "inquire variable dimension");
    variable.dimensions.push_back(dimension_name.data());
    variable.shape.push_back(dimension_size);
  }
  return true;
}

[[nodiscard]] NetcdfVariable require_variable(
    const NetcdfReadHandle& file,
    const std::string_view name) {
  NetcdfVariable variable;
  if (!inquire_variable_if_present(file, name, variable)) {
    throw PressureRefreshIoError(
        "missing pressure-refresh variable " + std::string(name) + " in " +
        file.path().string());
  }
  return variable;
}

void require_numeric_variable(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const NetcdfVariable& variable) {
  if (!numeric_type(variable.type)) {
    std::ostringstream message;
    message << "pressure-refresh variable " << name << " in " << file.path()
            << " has type " << nc_type_name(variable.type) << ", expected numeric";
    throw PressureRefreshIoError(message.str());
  }
}

void require_float_variable(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const NetcdfVariable& variable) {
  if (variable.type != NC_FLOAT && variable.type != NC_DOUBLE) {
    std::ostringstream message;
    message << "pressure-refresh variable " << name << " in " << file.path()
            << " has type " << nc_type_name(variable.type) << ", expected float or double";
    throw PressureRefreshIoError(message.str());
  }
}

[[nodiscard]] bool valid_scalar_or_time_shape(
    const NetcdfVariable& variable,
    const std::size_t time_index) {
  if (variable.dimensions.empty()) {
    return true;
  }
  return variable.dimensions.size() == 1 && variable.dimensions.front() == "Time" &&
         variable.shape.front() > time_index;
}

void require_p_top_shape(
    const NetcdfReadHandle& file,
    const NetcdfVariable& variable,
    const std::size_t time_index) {
  if (valid_scalar_or_time_shape(variable, time_index)) {
    return;
  }
  std::ostringstream message;
  message << "pressure-refresh variable P_TOP in " << file.path()
          << " must be scalar or Time scalar with Time length greater than " << time_index
          << ", found dimensions " << join_strings(variable.dimensions) << " shape"
          << join_shape(variable.shape);
  throw PressureRefreshIoError(message.str());
}

[[nodiscard]] PressureRefreshPTopSource inspect_p_top_source(
    const NetcdfReadHandle& file,
    const std::size_t time_index,
    std::vector<std::string>& missing_names) {
  if (has_global_attribute(file, "P_TOP")) {
    return PressureRefreshPTopSource::global_attribute;
  }

  NetcdfVariable variable;
  if (!inquire_variable_if_present(file, "P_TOP", variable)) {
    missing_names.push_back("P_TOP");
    return PressureRefreshPTopSource::missing;
  }
  require_numeric_variable(file, "P_TOP", variable);
  require_p_top_shape(file, variable, time_index);
  return variable.dimensions.empty() ? PressureRefreshPTopSource::scalar_variable
                                     : PressureRefreshPTopSource::time_variable;
}

[[nodiscard]] float checked_float_value(
    const double value,
    const std::string_view name,
    const std::filesystem::path& path) {
  if (!std::isfinite(value) || value < 0.0 ||
      value > static_cast<double>(std::numeric_limits<float>::max())) {
    std::ostringstream message;
    message << "pressure-refresh scalar " << name << " in " << path
            << " must be finite and non-negative, found " << value;
    throw PressureRefreshIoError(message.str());
  }
  return static_cast<float>(value);
}

[[nodiscard]] float read_p_top(
    const NetcdfReadHandle& file,
    const PressureRefreshPTopSource source,
    const std::size_t time_index) {
  double value = 0.0;
  if (source == PressureRefreshPTopSource::global_attribute) {
    nc_type type = NC_NAT;
    std::size_t length = 0;
    file.check(nc_inq_att(file.id(), NC_GLOBAL, "P_TOP", &type, &length), "inquire P_TOP attr");
    if (length != 1 || !numeric_type(type)) {
      throw PressureRefreshIoError("global attribute P_TOP must be a scalar numeric value");
    }
    file.check(nc_get_att_double(file.id(), NC_GLOBAL, "P_TOP", &value), "read P_TOP attr");
    return checked_float_value(value, "P_TOP", file.path());
  }

  const auto variable = require_variable(file, "P_TOP");
  require_numeric_variable(file, "P_TOP", variable);
  require_p_top_shape(file, variable, time_index);
  if (variable.dimensions.empty()) {
    file.check(nc_get_var_double(file.id(), variable.id, &value), "read scalar P_TOP");
  } else {
    const std::size_t start[1] = {time_index};
    const std::size_t count[1] = {1};
    file.check(nc_get_vara_double(file.id(), variable.id, start, count, &value), "read P_TOP");
  }
  return checked_float_value(value, "P_TOP", file.path());
}

void require_vector_shape(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const NetcdfVariable& variable,
    const std::string_view expected_dimension,
    const std::size_t expected_count,
    const std::size_t time_index) {
  require_float_variable(file, name, variable);
  const bool vector_shape =
      variable.dimensions.size() == 1 && variable.dimensions.front() == expected_dimension &&
      variable.shape.front() == expected_count;
  const bool time_vector_shape =
      variable.dimensions.size() == 2 && variable.dimensions.front() == "Time" &&
      variable.dimensions.back() == expected_dimension && variable.shape.front() > time_index &&
      variable.shape.back() == expected_count;
  if (vector_shape || time_vector_shape) {
    return;
  }

  std::ostringstream message;
  message << "pressure-refresh variable " << name << " in " << file.path()
          << " has dimensions " << join_strings(variable.dimensions) << " shape"
          << join_shape(variable.shape) << ", expected " << expected_dimension << '='
          << expected_count << " or Time>=" << (time_index + 1) << ','
          << expected_dimension << '=' << expected_count;
  throw PressureRefreshIoError(message.str());
}

[[nodiscard]] bool has_leading_time_dimension(const NetcdfVariable& variable) {
  return !variable.dimensions.empty() && variable.dimensions.front() == "Time";
}

[[nodiscard]] std::vector<float> read_vector_variable(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const std::string_view expected_dimension,
    const std::size_t expected_count,
    const std::size_t time_index,
    std::vector<std::string>& missing_names) {
  NetcdfVariable variable;
  if (!inquire_variable_if_present(file, name, variable)) {
    missing_names.push_back(std::string(name));
    return {};
  }
  require_vector_shape(file, name, variable, expected_dimension, expected_count, time_index);

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

[[nodiscard]] std::size_t inspect_vector_count(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const std::string_view expected_dimension,
    const std::size_t expected_count,
    const std::size_t time_index,
    std::vector<std::string>& missing_names) {
  NetcdfVariable variable;
  if (!inquire_variable_if_present(file, name, variable)) {
    missing_names.push_back(std::string(name));
    return 0;
  }
  require_vector_shape(file, name, variable, expected_dimension, expected_count, time_index);
  return expected_count;
}

void require_alb_shape(
    const NetcdfReadHandle& file,
    const NetcdfVariable& variable,
    const Grid& grid,
    const std::size_t time_index) {
  require_float_variable(file, "ALB", variable);
  const auto& config = grid.config();
  const auto nx = static_cast<std::size_t>(config.mass_nx);
  const auto ny = static_cast<std::size_t>(config.mass_ny);
  const auto nz = static_cast<std::size_t>(config.mass_nz);
  const bool mass_shape =
      variable.dimensions == std::vector<std::string>{"bottom_top", "south_north", "west_east"} &&
      variable.shape == std::vector<std::size_t>{nz, ny, nx};
  const bool time_mass_shape =
      variable.dimensions ==
          std::vector<std::string>{"Time", "bottom_top", "south_north", "west_east"} &&
      variable.shape.size() == 4 && variable.shape[0] > time_index &&
      variable.shape[1] == nz && variable.shape[2] == ny && variable.shape[3] == nx;
  if (mass_shape || time_mass_shape) {
    return;
  }

  std::ostringstream message;
  message << "pressure-refresh variable ALB in " << file.path() << " has dimensions "
          << join_strings(variable.dimensions) << " shape" << join_shape(variable.shape)
          << ", expected bottom_top=" << nz << ",south_north=" << ny
          << ",west_east=" << nx << " or Time>=" << (time_index + 1)
          << ",bottom_top=" << nz << ",south_north=" << ny << ",west_east=" << nx;
  throw PressureRefreshIoError(message.str());
}

void require_terrain_shape(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const NetcdfVariable& variable,
    const Grid& grid,
    const std::size_t time_index) {
  require_float_variable(file, name, variable);
  const auto& config = grid.config();
  const auto nx = static_cast<std::size_t>(config.mass_nx);
  const auto ny = static_cast<std::size_t>(config.mass_ny);
  const bool surface_shape =
      variable.dimensions == std::vector<std::string>{"south_north", "west_east"} &&
      variable.shape == std::vector<std::size_t>{ny, nx};
  const bool time_surface_shape =
      variable.dimensions ==
          std::vector<std::string>{"Time", "south_north", "west_east"} &&
      variable.shape.size() == 3 && variable.shape[0] > time_index &&
      variable.shape[1] == ny && variable.shape[2] == nx;
  if (surface_shape || time_surface_shape) {
    return;
  }

  std::ostringstream message;
  message << "pressure-refresh terrain variable " << name << " in " << file.path()
          << " has dimensions " << join_strings(variable.dimensions) << " shape"
          << join_shape(variable.shape) << ", expected south_north=" << ny
          << ",west_east=" << nx << " or Time>=" << (time_index + 1)
          << ",south_north=" << ny << ",west_east=" << nx;
  throw PressureRefreshIoError(message.str());
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(const FieldView3D<Real> field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(const FieldView3D<Real> field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nz(const FieldView3D<Real> field) noexcept {
  return field.nz - field.halo.k_lower - field.halo.k_upper;
}

[[nodiscard]] constexpr bool canonical_view(const FieldView3D<float> field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 && field.nz > 0 &&
         field.stride_i == 1 && field.stride_k == field.nx &&
         field.stride_j == field.nx * field.nz && field.halo.i_lower >= 0 &&
         field.halo.i_upper >= 0 && field.halo.j_lower >= 0 &&
         field.halo.j_upper >= 0 && field.halo.k_lower >= 0 &&
         field.halo.k_upper >= 0 && active_nx(field) > 0 && active_ny(field) > 0 &&
         active_nz(field) > 0;
}

void require_alb_target_view(const FieldView3D<float> alb, const Grid& grid) {
  if (!canonical_view(alb)) {
    throw PressureRefreshIoError(
        "pressure-refresh ALB target must be a non-null canonical 3D field view");
  }
  const auto& config = grid.config();
  if (active_nx(alb) == config.mass_nx && active_ny(alb) == config.mass_ny &&
      active_nz(alb) == config.mass_nz) {
    return;
  }

  std::ostringstream message;
  message << "pressure-refresh ALB target active shape is " << active_nx(alb) << 'x'
          << active_ny(alb) << 'x' << active_nz(alb) << ", expected " << config.mass_nx
          << 'x' << config.mass_ny << 'x' << config.mass_nz;
  throw PressureRefreshIoError(message.str());
}

[[nodiscard]] bool inspect_alb(
    const NetcdfReadHandle& file,
    const Grid& grid,
    const std::size_t time_index,
    KrosaPressureRefreshInputReport& report) {
  NetcdfVariable variable;
  if (!inquire_variable_if_present(file, "ALB", variable)) {
    report.missing_names.push_back("ALB");
    return false;
  }
  require_alb_shape(file, variable, grid, time_index);
  report.alb_nx = grid.config().mass_nx;
  report.alb_ny = grid.config().mass_ny;
  report.alb_nz = grid.config().mass_nz;
  report.alb_point_count = report.expected_alb_point_count;
  return true;
}

[[nodiscard]] bool inspect_terrain(
    const NetcdfReadHandle& file,
    const Grid& grid,
    const std::size_t time_index,
    KrosaPressureRefreshInputReport& report) {
  for (const std::string_view name : {"HGT", "HT"}) {
    NetcdfVariable variable;
    if (!inquire_variable_if_present(file, name, variable)) {
      continue;
    }
    require_terrain_shape(file, name, variable, grid, time_index);
    report.terrain_source_name = std::string(name);
    report.terrain_nx = grid.config().mass_nx;
    report.terrain_ny = grid.config().mass_ny;
    report.terrain_point_count = report.expected_terrain_point_count;
    return true;
  }

  append_missing_unique(report.missing_base_state_reconstruction_names, "HGT/HT");
  return false;
}

void read_alb(
    const NetcdfReadHandle& file,
    const Grid& grid,
    const std::size_t time_index,
    const FieldView3D<float> target) {
  const auto variable = require_variable(file, "ALB");
  require_alb_shape(file, variable, grid, time_index);

  const auto& config = grid.config();
  const auto nx = static_cast<std::size_t>(config.mass_nx);
  const auto ny = static_cast<std::size_t>(config.mass_ny);
  const auto point_count = checked_point_count(config.mass_nx, config.mass_ny, config.mass_nz);

  std::vector<std::size_t> start(variable.shape.size(), 0);
  std::vector<std::size_t> count = variable.shape;
  if (has_leading_time_dimension(variable)) {
    start.front() = time_index;
    count.front() = 1;
  }

  std::vector<float> buffer(point_count, 0.0F);
  file.check(
      nc_get_vara_float(file.id(), variable.id, start.data(), count.data(), buffer.data()),
      "read ALB");

  for (std::int32_t k = 0; k < config.mass_nz; ++k) {
    for (std::int32_t j = 0; j < config.mass_ny; ++j) {
      const auto source_row =
          (static_cast<std::size_t>(k) * ny + static_cast<std::size_t>(j)) * nx;
      for (std::int32_t i = 0; i < config.mass_nx; ++i) {
        target(
            i + target.halo.i_lower,
            j + target.halo.j_lower,
            k + target.halo.k_lower) = buffer[source_row + static_cast<std::size_t>(i)];
      }
    }
  }
}

[[nodiscard]] std::vector<float> read_terrain(
    const NetcdfReadHandle& file,
    const std::string_view name,
    const Grid& grid,
    const std::size_t time_index) {
  const auto variable = require_variable(file, name);
  require_terrain_shape(file, name, variable, grid, time_index);

  const auto& config = grid.config();
  const auto point_count = checked_point_count(config.mass_nx, config.mass_ny);

  std::vector<std::size_t> start(variable.shape.size(), 0);
  std::vector<std::size_t> count = variable.shape;
  if (has_leading_time_dimension(variable)) {
    start.front() = time_index;
    count.front() = 1;
  }

  std::vector<float> buffer(point_count, 0.0F);
  file.check(
      nc_get_vara_float(file.id(), variable.id, start.data(), count.data(), buffer.data()),
      "read terrain");
  return buffer;
}

void fill_report_shape(KrosaPressureRefreshInputReport& report, const Grid& grid) {
  const auto& config = grid.config();
  report.expected_mass_level_count = static_cast<std::size_t>(config.mass_nz);
  report.expected_full_level_count = static_cast<std::size_t>(config.full_nz);
  report.expected_alb_nx = config.mass_nx;
  report.expected_alb_ny = config.mass_ny;
  report.expected_alb_nz = config.mass_nz;
  report.expected_alb_point_count =
      checked_point_count(config.mass_nx, config.mass_ny, config.mass_nz);
  report.expected_terrain_nx = config.mass_nx;
  report.expected_terrain_ny = config.mass_ny;
  report.expected_terrain_point_count =
      checked_point_count(config.mass_nx, config.mass_ny);
}

void fill_base_state_missing_from_counts(KrosaPressureRefreshInputReport& report) {
  if (!report.p_top_present()) {
    append_missing_unique(report.missing_base_state_reconstruction_names, "P_TOP");
  }
  if (report.c3f_count != report.expected_full_level_count) {
    append_missing_unique(report.missing_base_state_reconstruction_names, "C3F");
  }
  if (report.c4f_count != report.expected_full_level_count) {
    append_missing_unique(report.missing_base_state_reconstruction_names, "C4F");
  }
  if (report.c3h_count != report.expected_mass_level_count) {
    append_missing_unique(report.missing_base_state_reconstruction_names, "C3H");
  }
  if (report.c4h_count != report.expected_mass_level_count) {
    append_missing_unique(report.missing_base_state_reconstruction_names, "C4H");
  }
}

}  // namespace

std::vector<std::string> krosa_pressure_refresh_required_names() {
  return required_names_ref();
}

KrosaPressureRefreshInputReport inspect_krosa_pressure_refresh_inputs(
    const std::filesystem::path& path,
    const Grid& grid,
    const KrosaPressureRefreshReadOptions& options) {
  require_grid(grid);

  KrosaPressureRefreshInputReport report;
  report.source_path = path;
  report.time_index = options.time_index;
  fill_report_shape(report, grid);

  const NetcdfReadHandle file(path);
  report.p_top_source = inspect_p_top_source(file, options.time_index, report.missing_names);
  report.c3f_count = inspect_vector_count(
      file,
      "C3F",
      "bottom_top_stag",
      report.expected_full_level_count,
      options.time_index,
      report.missing_names);
  report.c4f_count = inspect_vector_count(
      file,
      "C4F",
      "bottom_top_stag",
      report.expected_full_level_count,
      options.time_index,
      report.missing_names);
  report.c3h_count = inspect_vector_count(
      file,
      "C3H",
      "bottom_top",
      report.expected_mass_level_count,
      options.time_index,
      report.missing_names);
  report.c4h_count = inspect_vector_count(
      file,
      "C4H",
      "bottom_top",
      report.expected_mass_level_count,
      options.time_index,
      report.missing_names);
  (void)inspect_terrain(file, grid, options.time_index, report);
  fill_base_state_missing_from_counts(report);
  (void)inspect_alb(file, grid, options.time_index, report);
  return report;
}

KrosaPressureRefreshReadResult read_krosa_pressure_refresh_inputs(
    const std::filesystem::path& path,
    const Grid& grid,
    const FieldView3D<float> alb,
    const KrosaPressureRefreshReadOptions& options) {
  require_grid(grid);
  require_alb_target_view(alb, grid);

  KrosaPressureRefreshReadResult result;
  result.report = inspect_krosa_pressure_refresh_inputs(path, grid, options);
  result.metadata.source_path = path;
  result.metadata.time_index = options.time_index;
  result.metadata.p_top_source = result.report.p_top_source;
  result.metadata.terrain_source_name = result.report.terrain_source_name;
  result.metadata.terrain_nx = result.report.terrain_nx;
  result.metadata.terrain_ny = result.report.terrain_ny;

  const NetcdfReadHandle file(path);
  if (result.report.p_top_present()) {
    result.metadata.p_top_pa = read_p_top(file, result.report.p_top_source, options.time_index);
  }
  result.metadata.c3f = read_vector_variable(
      file,
      "C3F",
      "bottom_top_stag",
      result.report.expected_full_level_count,
      options.time_index,
      result.report.missing_names);
  result.metadata.c4f = read_vector_variable(
      file,
      "C4F",
      "bottom_top_stag",
      result.report.expected_full_level_count,
      options.time_index,
      result.report.missing_names);
  result.metadata.c3h = read_vector_variable(
      file,
      "C3H",
      "bottom_top",
      result.report.expected_mass_level_count,
      options.time_index,
      result.report.missing_names);
  result.metadata.c4h = read_vector_variable(
      file,
      "C4H",
      "bottom_top",
      result.report.expected_mass_level_count,
      options.time_index,
      result.report.missing_names);

  if (!result.report.terrain_source_name.empty()) {
    result.metadata.terrain_height_m =
        read_terrain(file, result.report.terrain_source_name, grid, options.time_index);
    result.report.terrain_loaded = true;
  }

  if (std::find(result.report.missing_names.begin(), result.report.missing_names.end(), "ALB") ==
      result.report.missing_names.end()) {
    read_alb(file, grid, options.time_index, alb);
    result.report.alb_loaded = true;
  }
  return result;
}

KrosaPressureRefreshReadResult read_krosa_pressure_refresh_inputs(
    const std::filesystem::path& path,
    const Grid& grid,
    FieldStorage3D<float>& alb,
    const KrosaPressureRefreshReadOptions& options) {
  return read_krosa_pressure_refresh_inputs(path, grid, alb.view(), options);
}

}  // namespace tywrf::io
