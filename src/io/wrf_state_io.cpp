#include "tywrf/io/wrf_state_io.hpp"

#include <netcdf.h>

#include <array>
#include <cstdint>
#include <exception>
#include <limits>
#include <sstream>
#include <string_view>
#include <vector>

namespace tywrf::io {
namespace {

class NetcdfFile {
 public:
  explicit NetcdfFile(const std::filesystem::path& path) : path_(path) {
    check(nc_open(path.string().c_str(), NC_NOWRITE, &id_), "open");
  }

  NetcdfFile(const NetcdfFile&) = delete;
  NetcdfFile& operator=(const NetcdfFile&) = delete;

  ~NetcdfFile() {
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
    throw WrfStateIoError(message.str());
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

[[nodiscard]] const std::vector<std::string>& core_field_names_ref() {
  static const std::vector<std::string> names = {
      "U",      "V",      "W",      "PH",     "PHB",   "T",    "MU",
      "MUB",    "P",      "PB",     "QVAPOR", "QCLOUD", "QRAIN",
      "QICE",   "QSNOW",  "QGRAUP", "QNICE",  "QNRAIN", "PSFC",
      "U10",    "V10",    "T2",     "Q2",     "RAINC",  "RAINNC"};
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
    case NC_STRING:
      return "string";
    default:
      return "unknown";
  }
}

[[nodiscard]] std::string join_strings(const std::vector<std::string>& values) {
  std::ostringstream joined;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      joined << ",";
    }
    joined << values[i];
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

[[nodiscard]] std::string schema_label(const DatasetSchema& schema) {
  if (schema.path.empty()) {
    return "WRF schema";
  }
  return schema.path.string();
}

[[nodiscard]] std::int32_t require_grid_dimension(
    const DatasetSchema& schema,
    const std::string_view name) {
  const auto* dimension = find_dimension(schema, name);
  if (dimension == nullptr) {
    throw WrfStateIoError(
        "WRF grid schema error for " + schema_label(schema) +
        ": missing dimension " + std::string(name));
  }
  if (dimension->size == 0 ||
      dimension->size > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    std::ostringstream message;
    message << "WRF grid schema error for " << schema_label(schema) << ": dimension " << name
            << " has unsupported size " << dimension->size;
    throw WrfStateIoError(message.str());
  }
  return static_cast<std::int32_t>(dimension->size);
}

void require_stagger_dimension(
    const DatasetSchema& schema,
    const std::string_view name,
    const std::int32_t expected) {
  const auto actual = require_grid_dimension(schema, name);
  if (actual != expected) {
    std::ostringstream message;
    message << "WRF grid schema error for " << schema_label(schema) << ": dimension " << name
            << " has size " << actual << ", expected " << expected;
    throw WrfStateIoError(message.str());
  }
}

[[nodiscard]] NetcdfVariable inquire_variable(
    const NetcdfFile& file,
    const std::string_view name) {
  NetcdfVariable variable;
  const int var_status = nc_inq_varid(file.id(), std::string(name).c_str(), &variable.id);
  if (var_status == NC_ENOTVAR) {
    throw WrfStateIoError(
        "missing WRF state variable " + std::string(name) + " in " + file.path().string());
  }
  file.check(var_status, "inquire variable id");

  int ndims = 0;
  file.check(nc_inq_varndims(file.id(), variable.id, &ndims), "inquire variable rank");

  std::vector<int> dimids(static_cast<std::size_t>(ndims));
  std::array<char, NC_MAX_NAME + 1> var_name{};
  int natts = 0;
  file.check(
      nc_inq_var(
          file.id(),
          variable.id,
          var_name.data(),
          &variable.type,
          &ndims,
          dimids.empty() ? nullptr : dimids.data(),
          &natts),
      "inquire variable");

  variable.dimensions.reserve(dimids.size());
  variable.shape.reserve(dimids.size());
  for (const int dimid : dimids) {
    std::array<char, NC_MAX_NAME + 1> dim_name{};
    std::size_t dim_size = 0;
    file.check(nc_inq_dim(file.id(), dimid, dim_name.data(), &dim_size), "inquire dimension");
    variable.dimensions.push_back(dim_name.data());
    variable.shape.push_back(dim_size);
  }
  return variable;
}

void require_variable_layout(
    const NetcdfFile& file,
    const std::string_view name,
    const NetcdfVariable& variable,
    const std::vector<std::string>& expected_dimensions,
    const std::vector<std::size_t>& expected_trailing_shape,
    const std::size_t time_index) {
  if (variable.type != NC_FLOAT) {
    std::ostringstream message;
    message << "variable " << name << " in " << file.path() << " has type "
            << nc_type_name(variable.type) << ", expected float";
    throw WrfStateIoError(message.str());
  }

  if (variable.dimensions != expected_dimensions) {
    std::ostringstream message;
    message << "variable " << name << " in " << file.path() << " has dimensions "
            << join_strings(variable.dimensions) << ", expected "
            << join_strings(expected_dimensions);
    throw WrfStateIoError(message.str());
  }

  if (variable.shape.size() != expected_trailing_shape.size() + 1) {
    std::ostringstream message;
    message << "variable " << name << " in " << file.path() << " has shape"
            << join_shape(variable.shape) << ", expected Time plus"
            << join_shape(expected_trailing_shape);
    throw WrfStateIoError(message.str());
  }

  if (variable.shape.front() <= time_index) {
    std::ostringstream message;
    message << "variable " << name << " in " << file.path() << " has Time length "
            << variable.shape.front() << ", cannot read time index " << time_index;
    throw WrfStateIoError(message.str());
  }

  bool trailing_shape_matches = true;
  for (std::size_t dim = 0; dim < expected_trailing_shape.size(); ++dim) {
    trailing_shape_matches =
        trailing_shape_matches && variable.shape[dim + 1] == expected_trailing_shape[dim];
  }
  if (!trailing_shape_matches) {
    std::ostringstream message;
    message << "variable " << name << " in " << file.path() << " has shape"
            << join_shape(variable.shape) << ", expected Time>=" << (time_index + 1)
            << join_shape(expected_trailing_shape);
    throw WrfStateIoError(message.str());
  }
}

[[nodiscard]] std::size_t checked_count_2d(const FieldLayout2D& layout) {
  if (!layout.valid() || layout.active_nx() <= 0 || layout.active_ny() <= 0) {
    throw WrfStateIoError("cannot load WRF 2D state variable into invalid field layout");
  }
  return static_cast<std::size_t>(layout.active_nx()) *
         static_cast<std::size_t>(layout.active_ny());
}

[[nodiscard]] std::size_t checked_count_3d(const FieldLayout3D& layout) {
  if (!layout.valid() || layout.active_nx() <= 0 || layout.active_ny() <= 0 ||
      layout.active_nz() <= 0) {
    throw WrfStateIoError("cannot load WRF 3D state variable into invalid field layout");
  }
  return static_cast<std::size_t>(layout.active_nx()) *
         static_cast<std::size_t>(layout.active_ny()) *
         static_cast<std::size_t>(layout.active_nz());
}

void load_2d_variable(
    const NetcdfFile& file,
    const std::string_view name,
    const std::vector<std::string>& expected_dimensions,
    FieldStorage2D<float>& field,
    const std::size_t time_index) {
  const auto layout = field.layout();
  const auto variable = inquire_variable(file, name);
  const auto nx = static_cast<std::size_t>(layout.active_nx());
  const auto ny = static_cast<std::size_t>(layout.active_ny());
  require_variable_layout(file, name, variable, expected_dimensions, {ny, nx}, time_index);

  std::vector<float> buffer(checked_count_2d(layout));
  const std::array<std::size_t, 3> start = {time_index, 0, 0};
  const std::array<std::size_t, 3> count = {1, ny, nx};
  file.check(
      nc_get_vara_float(file.id(), variable.id, start.data(), count.data(), buffer.data()),
      "read WRF 2D state variable");

  auto view = field.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    const auto source_row = static_cast<std::size_t>(j) * nx;
    for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
      view(i + layout.i_begin(), j + layout.j_begin()) =
          buffer[source_row + static_cast<std::size_t>(i)];
    }
  }
}

void load_3d_variable(
    const NetcdfFile& file,
    const std::string_view name,
    const std::vector<std::string>& expected_dimensions,
    FieldStorage3D<float>& field,
    const std::size_t time_index) {
  const auto layout = field.layout();
  const auto variable = inquire_variable(file, name);
  const auto nx = static_cast<std::size_t>(layout.active_nx());
  const auto ny = static_cast<std::size_t>(layout.active_ny());
  const auto nz = static_cast<std::size_t>(layout.active_nz());
  require_variable_layout(file, name, variable, expected_dimensions, {nz, ny, nx}, time_index);

  std::vector<float> buffer(checked_count_3d(layout));
  const std::array<std::size_t, 4> start = {time_index, 0, 0, 0};
  const std::array<std::size_t, 4> count = {1, nz, ny, nx};
  file.check(
      nc_get_vara_float(file.id(), variable.id, start.data(), count.data(), buffer.data()),
      "read WRF 3D state variable");

  auto view = field.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t k = 0; k < layout.active_nz(); ++k) {
      const auto source_plane =
          (static_cast<std::size_t>(k) * ny + static_cast<std::size_t>(j)) * nx;
      for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
        view(i + layout.i_begin(), j + layout.j_begin(), k + layout.k_begin()) =
            buffer[source_plane + static_cast<std::size_t>(i)];
      }
    }
  }
}

void load_named_variable(
    const NetcdfFile& file,
    const std::string_view name,
    State<float>& state,
    const std::size_t time_index) {
  const std::vector<std::string> mass_2d = {"Time", "south_north", "west_east"};
  const std::vector<std::string> mass_3d = {"Time", "bottom_top", "south_north", "west_east"};
  const std::vector<std::string> u_3d = {
      "Time", "bottom_top", "south_north", "west_east_stag"};
  const std::vector<std::string> v_3d = {
      "Time", "bottom_top", "south_north_stag", "west_east"};
  const std::vector<std::string> w_3d = {
      "Time", "bottom_top_stag", "south_north", "west_east"};

  if (name == "U") {
    load_3d_variable(file, name, u_3d, state.u, time_index);
  } else if (name == "V") {
    load_3d_variable(file, name, v_3d, state.v, time_index);
  } else if (name == "W") {
    load_3d_variable(file, name, w_3d, state.w, time_index);
  } else if (name == "PH") {
    load_3d_variable(file, name, w_3d, state.ph, time_index);
  } else if (name == "PHB") {
    load_3d_variable(file, name, w_3d, state.phb, time_index);
  } else if (name == "T") {
    load_3d_variable(file, name, mass_3d, state.t, time_index);
  } else if (name == "P") {
    load_3d_variable(file, name, mass_3d, state.p, time_index);
  } else if (name == "PB") {
    load_3d_variable(file, name, mass_3d, state.pb, time_index);
  } else if (name == "QVAPOR") {
    load_3d_variable(file, name, mass_3d, state.qvapor, time_index);
  } else if (name == "QCLOUD") {
    load_3d_variable(file, name, mass_3d, state.qcloud, time_index);
  } else if (name == "QRAIN") {
    load_3d_variable(file, name, mass_3d, state.qrain, time_index);
  } else if (name == "QICE") {
    load_3d_variable(file, name, mass_3d, state.qice, time_index);
  } else if (name == "QSNOW") {
    load_3d_variable(file, name, mass_3d, state.qsnow, time_index);
  } else if (name == "QGRAUP") {
    load_3d_variable(file, name, mass_3d, state.qgraup, time_index);
  } else if (name == "QNICE") {
    load_3d_variable(file, name, mass_3d, state.qnice, time_index);
  } else if (name == "QNRAIN") {
    load_3d_variable(file, name, mass_3d, state.qnrain, time_index);
  } else if (name == "MU") {
    load_2d_variable(file, name, mass_2d, state.mu, time_index);
  } else if (name == "MUB") {
    load_2d_variable(file, name, mass_2d, state.mub, time_index);
  } else if (name == "PSFC") {
    load_2d_variable(file, name, mass_2d, state.psfc, time_index);
  } else if (name == "U10") {
    load_2d_variable(file, name, mass_2d, state.u10, time_index);
  } else if (name == "V10") {
    load_2d_variable(file, name, mass_2d, state.v10, time_index);
  } else if (name == "T2") {
    load_2d_variable(file, name, mass_2d, state.t2, time_index);
  } else if (name == "Q2") {
    load_2d_variable(file, name, mass_2d, state.q2, time_index);
  } else if (name == "RAINC") {
    load_2d_variable(file, name, mass_2d, state.rainc, time_index);
  } else if (name == "RAINNC") {
    load_2d_variable(file, name, mass_2d, state.rainnc, time_index);
  } else {
    throw WrfStateIoError("unsupported WRF state variable selection: " + std::string(name));
  }
}

}  // namespace

std::vector<std::string> wrf_state_core_field_names() {
  return core_field_names_ref();
}

Grid derive_grid_from_wrf_schema(const DatasetSchema& schema, const Halo3D halo) {
  const auto mass_nx = require_grid_dimension(schema, "west_east");
  const auto mass_ny = require_grid_dimension(schema, "south_north");
  const auto mass_nz = require_grid_dimension(schema, "bottom_top");
  const auto full_nz = require_grid_dimension(schema, "bottom_top_stag");

  require_stagger_dimension(schema, "west_east_stag", mass_nx + 1);
  require_stagger_dimension(schema, "south_north_stag", mass_ny + 1);
  if (full_nz != mass_nz + 1) {
    std::ostringstream message;
    message << "WRF grid schema error for " << schema_label(schema)
            << ": dimension bottom_top_stag has size " << full_nz << ", expected "
            << (mass_nz + 1);
    throw WrfStateIoError(message.str());
  }

  return Grid({mass_nx, mass_ny, mass_nz, full_nz, halo});
}

Grid derive_grid_from_wrf_file(const std::filesystem::path& path, const Halo3D halo) {
  try {
    return derive_grid_from_wrf_schema(read_netcdf_schema(path), halo);
  } catch (const WrfStateIoError&) {
    throw;
  } catch (const std::exception& error) {
    throw WrfStateIoError(
        "failed to derive WRF grid from " + path.string() + ": " + error.what());
  }
}

void load_wrf_state(
    const std::filesystem::path& path,
    State<float>& state,
    const WrfStateReadOptions& options) {
  const NetcdfFile file(path);
  const auto& selected = options.variables.empty() ? core_field_names_ref() : options.variables;
  for (const auto& name : selected) {
    load_named_variable(file, name, state, options.time_index);
  }
}

void write_wrf_state(
    const std::filesystem::path& path,
    const State<float>& state,
    const WrfStateWriteOptions& options) {
  (void)path;
  (void)state;
  (void)options;
  throw WrfStateIoError(
      "WRF state NetCDF writer is not implemented yet; Phase 2 currently defines "
      "the read-side shape contract only");
}

}  // namespace tywrf::io
