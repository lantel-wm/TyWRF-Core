#include "tywrf/io/wrf_state_io.hpp"

#include <netcdf.h>

#include <algorithm>
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
  enum class Mode {
    read,
    replace,
  };

  explicit NetcdfFile(const std::filesystem::path& path, const Mode mode = Mode::read)
      : path_(path) {
    if (mode == Mode::read) {
      check(nc_open(path.string().c_str(), NC_NOWRITE, &id_), "open");
    } else {
      check(nc_create(path.string().c_str(), NC_CLOBBER, &id_), "create");
    }
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

[[nodiscard]] const std::vector<std::string>& writer_missing_field_names_ref() {
  static const std::vector<std::string> names = {"Times", "XLAT", "XLONG", "HGT"};
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
    throw WrfStateIoError("cannot access WRF 2D state variable with invalid field layout");
  }
  return static_cast<std::size_t>(layout.active_nx()) *
         static_cast<std::size_t>(layout.active_ny());
}

[[nodiscard]] std::size_t checked_count_3d(const FieldLayout3D& layout) {
  if (!layout.valid() || layout.active_nx() <= 0 || layout.active_ny() <= 0 ||
      layout.active_nz() <= 0) {
    throw WrfStateIoError("cannot access WRF 3D state variable with invalid field layout");
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

enum class WrfStateVariableLayout {
  mass_2d,
  mass_3d,
  u_3d,
  v_3d,
  w_3d,
};

struct OutputDimensions {
  int time = -1;
  int bottom_top = -1;
  int bottom_top_stag = -1;
  int south_north = -1;
  int south_north_stag = -1;
  int west_east = -1;
  int west_east_stag = -1;
};

struct DefinedVariable {
  std::string name;
  int id = -1;
};

[[nodiscard]] bool is_writable_state_variable(const std::string_view name) {
  const auto& names = core_field_names_ref();
  return std::find(names.begin(), names.end(), name) != names.end();
}

[[nodiscard]] std::vector<std::string> selected_write_variables(
    const WrfStateWriteOptions& options) {
  const auto selected = options.variables.empty() ? core_field_names_ref() : options.variables;
  std::vector<std::string> variables(selected.begin(), selected.end());
  for (std::size_t index = 0; index < variables.size(); ++index) {
    const auto& name = variables[index];
    if (!is_writable_state_variable(name)) {
      std::ostringstream message;
      message << "unsupported WRF state writer variable selection: " << name
              << "; writable State-backed variables are "
              << join_strings(core_field_names_ref()) << "; not yet written by this "
              << "State writer: " << join_strings(writer_missing_field_names_ref());
      throw WrfStateIoError(message.str());
    }
    const auto current = variables.begin() + static_cast<std::ptrdiff_t>(index);
    if (std::find(variables.begin(), current, name) != current) {
      throw WrfStateIoError("duplicate WRF state writer variable selection: " + name);
    }
  }
  return variables;
}

[[nodiscard]] WrfStateVariableLayout variable_layout_for_name(const std::string_view name) {
  if (name == "U") {
    return WrfStateVariableLayout::u_3d;
  }
  if (name == "V") {
    return WrfStateVariableLayout::v_3d;
  }
  if (name == "W" || name == "PH" || name == "PHB") {
    return WrfStateVariableLayout::w_3d;
  }
  if (name == "MU" || name == "MUB" || name == "PSFC" || name == "U10" || name == "V10" ||
      name == "T2" || name == "Q2" || name == "RAINC" || name == "RAINNC") {
    return WrfStateVariableLayout::mass_2d;
  }
  return WrfStateVariableLayout::mass_3d;
}

void require_output_grid(const GridConfig& config) {
  if (config.mass_nx <= 0 || config.mass_ny <= 0 || config.mass_nz <= 0 ||
      config.full_nz <= 0) {
    throw WrfStateIoError("cannot write WRF state from invalid grid dimensions");
  }
  if (config.full_nz != config.mass_nz + 1) {
    std::ostringstream message;
    message << "cannot write WRF state: full_nz is " << config.full_nz
            << ", expected mass_nz + 1 = " << (config.mass_nz + 1);
    throw WrfStateIoError(message.str());
  }
}

int define_dimension(
    const NetcdfFile& file,
    const std::string_view name,
    const std::size_t size) {
  int id = -1;
  file.check(nc_def_dim(file.id(), std::string(name).c_str(), size, &id), "define dimension");
  return id;
}

[[nodiscard]] OutputDimensions define_output_dimensions(
    const NetcdfFile& file,
    const GridConfig& config) {
  require_output_grid(config);

  OutputDimensions dimensions;
  dimensions.time = define_dimension(file, "Time", NC_UNLIMITED);
  dimensions.bottom_top =
      define_dimension(file, "bottom_top", static_cast<std::size_t>(config.mass_nz));
  dimensions.bottom_top_stag =
      define_dimension(file, "bottom_top_stag", static_cast<std::size_t>(config.full_nz));
  dimensions.south_north =
      define_dimension(file, "south_north", static_cast<std::size_t>(config.mass_ny));
  dimensions.south_north_stag =
      define_dimension(file, "south_north_stag", static_cast<std::size_t>(config.mass_ny + 1));
  dimensions.west_east =
      define_dimension(file, "west_east", static_cast<std::size_t>(config.mass_nx));
  dimensions.west_east_stag =
      define_dimension(file, "west_east_stag", static_cast<std::size_t>(config.mass_nx + 1));
  return dimensions;
}

[[nodiscard]] std::vector<int> dimension_ids_for_layout(
    const OutputDimensions& dimensions,
    const WrfStateVariableLayout layout) {
  switch (layout) {
    case WrfStateVariableLayout::mass_2d:
      return {dimensions.time, dimensions.south_north, dimensions.west_east};
    case WrfStateVariableLayout::mass_3d:
      return {
          dimensions.time,
          dimensions.bottom_top,
          dimensions.south_north,
          dimensions.west_east};
    case WrfStateVariableLayout::u_3d:
      return {
          dimensions.time,
          dimensions.bottom_top,
          dimensions.south_north,
          dimensions.west_east_stag};
    case WrfStateVariableLayout::v_3d:
      return {
          dimensions.time,
          dimensions.bottom_top,
          dimensions.south_north_stag,
          dimensions.west_east};
    case WrfStateVariableLayout::w_3d:
      return {
          dimensions.time,
          dimensions.bottom_top_stag,
          dimensions.south_north,
          dimensions.west_east};
  }
  throw WrfStateIoError("internal error: unsupported WRF output variable layout");
}

void write_text_attribute(
    const NetcdfFile& file,
    const int variable_id,
    const std::string_view name,
    const std::string_view value) {
  file.check(
      nc_put_att_text(
          file.id(),
          variable_id,
          std::string(name).c_str(),
          value.size(),
          value.data()),
      "write variable attribute");
}

[[nodiscard]] std::string stagger_for_layout(const WrfStateVariableLayout layout) {
  switch (layout) {
    case WrfStateVariableLayout::u_3d:
      return "X";
    case WrfStateVariableLayout::v_3d:
      return "Y";
    case WrfStateVariableLayout::w_3d:
      return "Z";
    case WrfStateVariableLayout::mass_2d:
    case WrfStateVariableLayout::mass_3d:
      return "";
  }
  return "";
}

[[nodiscard]] std::string memory_order_for_layout(const WrfStateVariableLayout layout) {
  return layout == WrfStateVariableLayout::mass_2d ? "XY" : "XYZ";
}

[[nodiscard]] int define_output_variable(
    const NetcdfFile& file,
    const std::string& name,
    const OutputDimensions& dimensions) {
  const auto layout = variable_layout_for_name(name);
  auto dimids = dimension_ids_for_layout(dimensions, layout);
  int variable_id = -1;
  file.check(
      nc_def_var(
          file.id(),
          name.c_str(),
          NC_FLOAT,
          static_cast<int>(dimids.size()),
          dimids.data(),
          &variable_id),
      "define WRF state variable");
  write_text_attribute(file, variable_id, "MemoryOrder", memory_order_for_layout(layout));
  write_text_attribute(file, variable_id, "stagger", stagger_for_layout(layout));
  return variable_id;
}

void require_2d_write_shape(
    const std::string_view name,
    const FieldLayout2D& layout,
    const std::int32_t expected_nx,
    const std::int32_t expected_ny) {
  (void)checked_count_2d(layout);
  if (layout.active_nx() != expected_nx || layout.active_ny() != expected_ny) {
    std::ostringstream message;
    message << "cannot write WRF state variable " << name << ": active shape is "
            << layout.active_nx() << "x" << layout.active_ny() << ", expected "
            << expected_nx << "x" << expected_ny;
    throw WrfStateIoError(message.str());
  }
}

void require_3d_write_shape(
    const std::string_view name,
    const FieldLayout3D& layout,
    const std::int32_t expected_nx,
    const std::int32_t expected_ny,
    const std::int32_t expected_nz) {
  (void)checked_count_3d(layout);
  if (layout.active_nx() != expected_nx || layout.active_ny() != expected_ny ||
      layout.active_nz() != expected_nz) {
    std::ostringstream message;
    message << "cannot write WRF state variable " << name << ": active shape is "
            << layout.active_nx() << "x" << layout.active_ny() << "x"
            << layout.active_nz() << ", expected " << expected_nx << "x"
            << expected_ny << "x" << expected_nz;
    throw WrfStateIoError(message.str());
  }
}

void write_2d_variable(
    const NetcdfFile& file,
    const std::string_view name,
    const int variable_id,
    const FieldStorage2D<float>& field,
    const std::size_t time_index) {
  const auto layout = field.layout();
  const auto nx = static_cast<std::size_t>(layout.active_nx());
  const auto ny = static_cast<std::size_t>(layout.active_ny());
  std::vector<float> buffer(checked_count_2d(layout));

  const auto view = field.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    const auto target_row = static_cast<std::size_t>(j) * nx;
    for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
      buffer[target_row + static_cast<std::size_t>(i)] =
          view(i + layout.i_begin(), j + layout.j_begin());
    }
  }

  const std::array<std::size_t, 3> start = {time_index, 0, 0};
  const std::array<std::size_t, 3> count = {1, ny, nx};
  file.check(
      nc_put_vara_float(file.id(), variable_id, start.data(), count.data(), buffer.data()),
      "write WRF 2D state variable " + std::string(name));
}

void write_3d_variable(
    const NetcdfFile& file,
    const std::string_view name,
    const int variable_id,
    const FieldStorage3D<float>& field,
    const std::size_t time_index) {
  const auto layout = field.layout();
  const auto nx = static_cast<std::size_t>(layout.active_nx());
  const auto ny = static_cast<std::size_t>(layout.active_ny());
  const auto nz = static_cast<std::size_t>(layout.active_nz());
  std::vector<float> buffer(checked_count_3d(layout));

  const auto view = field.view();
  for (std::int32_t k = 0; k < layout.active_nz(); ++k) {
    for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
      const auto target_plane =
          (static_cast<std::size_t>(k) * ny + static_cast<std::size_t>(j)) * nx;
      for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
        buffer[target_plane + static_cast<std::size_t>(i)] =
            view(i + layout.i_begin(), j + layout.j_begin(), k + layout.k_begin());
      }
    }
  }

  const std::array<std::size_t, 4> start = {time_index, 0, 0, 0};
  const std::array<std::size_t, 4> count = {1, nz, ny, nx};
  file.check(
      nc_put_vara_float(file.id(), variable_id, start.data(), count.data(), buffer.data()),
      "write WRF 3D state variable " + std::string(name));
}

void write_checked_2d_variable(
    const NetcdfFile& file,
    const std::string_view name,
    const int variable_id,
    const FieldStorage2D<float>& field,
    const std::int32_t expected_nx,
    const std::int32_t expected_ny,
    const std::size_t time_index) {
  require_2d_write_shape(name, field.layout(), expected_nx, expected_ny);
  write_2d_variable(file, name, variable_id, field, time_index);
}

void write_checked_3d_variable(
    const NetcdfFile& file,
    const std::string_view name,
    const int variable_id,
    const FieldStorage3D<float>& field,
    const std::int32_t expected_nx,
    const std::int32_t expected_ny,
    const std::int32_t expected_nz,
    const std::size_t time_index) {
  require_3d_write_shape(name, field.layout(), expected_nx, expected_ny, expected_nz);
  write_3d_variable(file, name, variable_id, field, time_index);
}

void write_named_variable(
    const NetcdfFile& file,
    const std::string_view name,
    const int variable_id,
    const State<float>& state,
    const std::size_t time_index) {
  const auto& config = state.grid.config();
  if (name == "U") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.u,
        config.mass_nx + 1,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "V") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.v,
        config.mass_nx,
        config.mass_ny + 1,
        config.mass_nz,
        time_index);
  } else if (name == "W") {
    write_checked_3d_variable(
        file, name, variable_id, state.w, config.mass_nx, config.mass_ny, config.full_nz, time_index);
  } else if (name == "PH") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.ph,
        config.mass_nx,
        config.mass_ny,
        config.full_nz,
        time_index);
  } else if (name == "PHB") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.phb,
        config.mass_nx,
        config.mass_ny,
        config.full_nz,
        time_index);
  } else if (name == "T") {
    write_checked_3d_variable(
        file, name, variable_id, state.t, config.mass_nx, config.mass_ny, config.mass_nz, time_index);
  } else if (name == "P") {
    write_checked_3d_variable(
        file, name, variable_id, state.p, config.mass_nx, config.mass_ny, config.mass_nz, time_index);
  } else if (name == "PB") {
    write_checked_3d_variable(
        file, name, variable_id, state.pb, config.mass_nx, config.mass_ny, config.mass_nz, time_index);
  } else if (name == "QVAPOR") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.qvapor,
        config.mass_nx,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "QCLOUD") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.qcloud,
        config.mass_nx,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "QRAIN") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.qrain,
        config.mass_nx,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "QICE") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.qice,
        config.mass_nx,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "QSNOW") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.qsnow,
        config.mass_nx,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "QGRAUP") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.qgraup,
        config.mass_nx,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "QNICE") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.qnice,
        config.mass_nx,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "QNRAIN") {
    write_checked_3d_variable(
        file,
        name,
        variable_id,
        state.qnrain,
        config.mass_nx,
        config.mass_ny,
        config.mass_nz,
        time_index);
  } else if (name == "MU") {
    write_checked_2d_variable(
        file, name, variable_id, state.mu, config.mass_nx, config.mass_ny, time_index);
  } else if (name == "MUB") {
    write_checked_2d_variable(
        file, name, variable_id, state.mub, config.mass_nx, config.mass_ny, time_index);
  } else if (name == "PSFC") {
    write_checked_2d_variable(
        file, name, variable_id, state.psfc, config.mass_nx, config.mass_ny, time_index);
  } else if (name == "U10") {
    write_checked_2d_variable(
        file, name, variable_id, state.u10, config.mass_nx, config.mass_ny, time_index);
  } else if (name == "V10") {
    write_checked_2d_variable(
        file, name, variable_id, state.v10, config.mass_nx, config.mass_ny, time_index);
  } else if (name == "T2") {
    write_checked_2d_variable(
        file, name, variable_id, state.t2, config.mass_nx, config.mass_ny, time_index);
  } else if (name == "Q2") {
    write_checked_2d_variable(
        file, name, variable_id, state.q2, config.mass_nx, config.mass_ny, time_index);
  } else if (name == "RAINC") {
    write_checked_2d_variable(
        file, name, variable_id, state.rainc, config.mass_nx, config.mass_ny, time_index);
  } else if (name == "RAINNC") {
    write_checked_2d_variable(
        file, name, variable_id, state.rainnc, config.mass_nx, config.mass_ny, time_index);
  } else {
    throw WrfStateIoError("unsupported WRF state writer variable selection: " + std::string(name));
  }
}

}  // namespace

std::vector<std::string> wrf_state_core_field_names() {
  return core_field_names_ref();
}

std::vector<std::string> wrf_state_writable_field_names() {
  return core_field_names_ref();
}

std::vector<std::string> wrf_state_writer_missing_field_names() {
  return writer_missing_field_names_ref();
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
  const auto variables = selected_write_variables(options);
  const NetcdfFile file(path, NetcdfFile::Mode::replace);
  const auto dimensions = define_output_dimensions(file, state.grid.config());

  std::vector<DefinedVariable> defined;
  defined.reserve(variables.size());
  for (const auto& name : variables) {
    defined.push_back({name, define_output_variable(file, name, dimensions)});
  }

  file.check(nc_enddef(file.id()), "end definitions");

  for (const auto& variable : defined) {
    write_named_variable(file, variable.name, variable.id, state, options.time_index);
  }
}

}  // namespace tywrf::io
