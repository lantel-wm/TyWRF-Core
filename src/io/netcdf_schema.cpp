#include "tywrf/io/netcdf_schema.hpp"

#include <netcdf.h>

#include <array>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>

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

[[nodiscard]] bool should_read_variable(
    const std::string_view name,
    const std::vector<std::string>& selected_variables) {
  if (selected_variables.empty()) {
    return true;
  }
  for (const auto& selected : selected_variables) {
    if (name == selected) {
      return true;
    }
  }
  return false;
}

struct DimensionExpectation {
  std::string name;
  std::size_t size = 0;
  std::optional<bool> unlimited;
};

struct VariableExpectation {
  std::string name;
  std::string type_name;
  std::vector<std::string> dimensions;
};

void require_dimensions(
    const DatasetSchema& schema,
    const std::vector<DimensionExpectation>& expected,
    std::vector<SchemaIssue>& issues) {
  for (const auto& dimension : expected) {
    const auto* actual = find_dimension(schema, dimension.name);
    if (actual == nullptr) {
      issues.push_back({SchemaSeverity::error, "missing dimension: " + std::string(dimension.name)});
      continue;
    }
    if (actual->size != dimension.size) {
      std::ostringstream message;
      message << "dimension " << dimension.name << " has size " << actual->size
              << ", expected " << dimension.size;
      issues.push_back({SchemaSeverity::error, message.str()});
    }
    if (dimension.unlimited.has_value() && actual->unlimited != *dimension.unlimited) {
      std::ostringstream message;
      message << "dimension " << dimension.name
              << (actual->unlimited ? " is unlimited" : " is fixed")
              << ", expected " << (*dimension.unlimited ? "unlimited" : "fixed");
      issues.push_back({SchemaSeverity::error, message.str()});
    }
  }
}

[[nodiscard]] std::string join_dimensions(const std::vector<std::string>& dimensions) {
  std::ostringstream joined;
  for (std::size_t i = 0; i < dimensions.size(); ++i) {
    if (i != 0) {
      joined << ",";
    }
    joined << dimensions[i];
  }
  return joined.str();
}

[[nodiscard]] std::vector<std::size_t> expected_shape_for(
    const DatasetSchema& schema,
    const std::vector<std::string>& dimensions,
    const std::string_view variable_name,
    std::vector<SchemaIssue>& issues) {
  std::vector<std::size_t> shape;
  shape.reserve(dimensions.size());
  for (const auto& dimension_name : dimensions) {
    const auto* dimension = find_dimension(schema, dimension_name);
    if (dimension == nullptr) {
      issues.push_back(
          {SchemaSeverity::error,
           "variable " + std::string(variable_name) +
               " references missing expected dimension: " + dimension_name});
      continue;
    }
    shape.push_back(dimension->size);
  }
  return shape;
}

void require_variables(
    const DatasetSchema& schema,
    const std::vector<VariableExpectation>& required,
    std::vector<SchemaIssue>& issues) {
  for (const auto& expected : required) {
    const auto* actual = find_variable(schema, expected.name);
    if (actual == nullptr) {
      issues.push_back({SchemaSeverity::error, "missing variable: " + expected.name});
      continue;
    }

    if (actual->type_name != expected.type_name) {
      std::ostringstream message;
      message << "variable " << expected.name << " has type " << actual->type_name
              << ", expected " << expected.type_name;
      issues.push_back({SchemaSeverity::error, message.str()});
    }

    if (actual->dimensions != expected.dimensions) {
      std::ostringstream message;
      message << "variable " << expected.name << " has dimensions "
              << join_dimensions(actual->dimensions) << ", expected "
              << join_dimensions(expected.dimensions);
      issues.push_back({SchemaSeverity::error, message.str()});
    }

    const auto expected_shape =
        expected_shape_for(schema, expected.dimensions, expected.name, issues);
    if (expected_shape.size() == expected.dimensions.size() &&
        actual->shape != expected_shape) {
      std::ostringstream message;
      message << "variable " << expected.name << " has shape";
      for (const auto size : actual->shape) {
        message << ' ' << size;
      }
      message << ", expected";
      for (const auto size : expected_shape) {
        message << ' ' << size;
      }
      issues.push_back({SchemaSeverity::error, message.str()});
    }
  }
}

[[nodiscard]] VariableExpectation variable(
    std::string name,
    std::string type_name,
    std::vector<std::string> dimensions) {
  return {std::move(name), std::move(type_name), std::move(dimensions)};
}

void append_boundary_field(
    std::vector<VariableExpectation>& variables,
    const std::string_view base_name,
    std::vector<std::string> x_dimensions,
    std::vector<std::string> y_dimensions) {
  const auto append = [&](const std::string_view suffix, const std::vector<std::string>& dimensions) {
    variables.push_back(variable(std::string(base_name) + std::string(suffix), "float", dimensions));
  };

  append("_BXS", x_dimensions);
  append("_BXE", x_dimensions);
  append("_BYS", y_dimensions);
  append("_BYE", y_dimensions);
  append("_BTXS", x_dimensions);
  append("_BTXE", x_dimensions);
  append("_BTYS", y_dimensions);
  append("_BTYE", y_dimensions);
}

[[nodiscard]] std::vector<VariableExpectation> wrfinput_core_variables() {
  const std::vector<std::string> mass_2d = {"Time", "south_north", "west_east"};
  const std::vector<std::string> mass_3d = {"Time", "bottom_top", "south_north", "west_east"};
  const std::vector<std::string> w_3d = {"Time", "bottom_top_stag", "south_north", "west_east"};

  return {
      variable("Times", "char", {"Time", "DateStrLen"}),
      variable("XLAT", "float", mass_2d),
      variable("XLONG", "float", mass_2d),
      variable("HGT", "float", mass_2d),
      variable("U", "float", {"Time", "bottom_top", "south_north", "west_east_stag"}),
      variable("V", "float", {"Time", "bottom_top", "south_north_stag", "west_east"}),
      variable("W", "float", w_3d),
      variable("PH", "float", w_3d),
      variable("PHB", "float", w_3d),
      variable("T", "float", mass_3d),
      variable("MU", "float", mass_2d),
      variable("MUB", "float", mass_2d),
      variable("P", "float", mass_3d),
      variable("PB", "float", mass_3d),
      variable("QVAPOR", "float", mass_3d),
      variable("QCLOUD", "float", mass_3d),
      variable("QRAIN", "float", mass_3d),
      variable("QICE", "float", mass_3d),
      variable("QSNOW", "float", mass_3d),
      variable("QGRAUP", "float", mass_3d),
      variable("QNICE", "float", mass_3d),
      variable("QNRAIN", "float", mass_3d),
      variable("PSFC", "float", mass_2d),
      variable("U10", "float", mass_2d),
      variable("V10", "float", mass_2d),
      variable("T2", "float", mass_2d),
      variable("Q2", "float", mass_2d),
  };
}

[[nodiscard]] std::vector<VariableExpectation> wrfbdy_core_variables() {
  std::vector<VariableExpectation> variables;
  variables.push_back(variable("Times", "char", {"Time", "DateStrLen"}));

  append_boundary_field(
      variables,
      "U",
      {"Time", "bdy_width", "bottom_top", "south_north"},
      {"Time", "bdy_width", "bottom_top", "west_east_stag"});
  append_boundary_field(
      variables,
      "V",
      {"Time", "bdy_width", "bottom_top", "south_north_stag"},
      {"Time", "bdy_width", "bottom_top", "west_east"});
  append_boundary_field(
      variables,
      "W",
      {"Time", "bdy_width", "bottom_top_stag", "south_north"},
      {"Time", "bdy_width", "bottom_top_stag", "west_east"});
  append_boundary_field(
      variables,
      "PH",
      {"Time", "bdy_width", "bottom_top_stag", "south_north"},
      {"Time", "bdy_width", "bottom_top_stag", "west_east"});

  const std::vector<std::string> x_mass_3d = {"Time", "bdy_width", "bottom_top", "south_north"};
  const std::vector<std::string> y_mass_3d = {"Time", "bdy_width", "bottom_top", "west_east"};
  for (const std::string_view name :
       {"T", "QVAPOR", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP", "QNICE", "QNRAIN"}) {
    append_boundary_field(variables, name, x_mass_3d, y_mass_3d);
  }

  append_boundary_field(
      variables,
      "MU",
      {"Time", "bdy_width", "south_north"},
      {"Time", "bdy_width", "west_east"});

  return variables;
}

[[nodiscard]] std::vector<VariableExpectation> wrffdda_core_variables() {
  const std::vector<std::string> nudging_3d = {"Time", "bottom_top", "south_north", "west_east"};
  const std::vector<std::string> mu_nudging = {"Time", "one_stag", "south_north", "west_east"};
  return {
      variable("Times", "char", {"Time", "DateStrLen"}),
      variable("U_NDG_OLD", "float", nudging_3d),
      variable("V_NDG_OLD", "float", nudging_3d),
      variable("T_NDG_OLD", "float", nudging_3d),
      variable("Q_NDG_OLD", "float", nudging_3d),
      variable("PH_NDG_OLD", "float", nudging_3d),
      variable("MU_NDG_OLD", "float", mu_nudging),
      variable("U_NDG_NEW", "float", nudging_3d),
      variable("V_NDG_NEW", "float", nudging_3d),
      variable("T_NDG_NEW", "float", nudging_3d),
      variable("Q_NDG_NEW", "float", nudging_3d),
      variable("PH_NDG_NEW", "float", nudging_3d),
      variable("MU_NDG_NEW", "float", mu_nudging),
  };
}

}  // namespace

DatasetSchema read_netcdf_schema(
    const std::filesystem::path& path,
    const std::vector<std::string>& selected_variables) {
  NetcdfFile file(path);

  int ndims = 0;
  int nvars = 0;
  int ngatts = 0;
  int unlimdimid = -1;
  file.check(nc_inq(file.id(), &ndims, &nvars, &ngatts, &unlimdimid), "inquire dataset");

  DatasetSchema schema;
  schema.path = path;
  schema.dimensions.reserve(static_cast<std::size_t>(ndims));

  for (int dimid = 0; dimid < ndims; ++dimid) {
    std::array<char, NC_MAX_NAME + 1> name{};
    std::size_t size = 0;
    file.check(nc_inq_dim(file.id(), dimid, name.data(), &size), "inquire dimension");
    schema.dimensions.push_back({name.data(), size, dimid == unlimdimid});
  }

  schema.variables.reserve(static_cast<std::size_t>(nvars));
  for (int varid = 0; varid < nvars; ++varid) {
    std::array<char, NC_MAX_NAME + 1> name{};
    nc_type type = NC_NAT;
    int var_ndims = 0;
    int natts = 0;
    file.check(nc_inq_varndims(file.id(), varid, &var_ndims), "inquire variable rank");

    std::vector<int> dimids(static_cast<std::size_t>(var_ndims));
    file.check(
        nc_inq_var(
            file.id(),
            varid,
            name.data(),
            &type,
            &var_ndims,
            dimids.empty() ? nullptr : dimids.data(),
            &natts),
        "inquire variable");

    if (!should_read_variable(name.data(), selected_variables)) {
      continue;
    }

    VariableInfo variable;
    variable.name = name.data();
    variable.type_name = nc_type_name(type);
    variable.dimensions.reserve(dimids.size());
    variable.shape.reserve(dimids.size());

    for (const int dimid : dimids) {
      const auto& dimension = schema.dimensions.at(static_cast<std::size_t>(dimid));
      variable.dimensions.push_back(dimension.name);
      variable.shape.push_back(dimension.size);
    }

    schema.variables.push_back(std::move(variable));
  }

  return schema;
}

const DimensionInfo* find_dimension(const DatasetSchema& schema, const std::string_view name) noexcept {
  for (const auto& dimension : schema.dimensions) {
    if (dimension.name == name) {
      return &dimension;
    }
  }
  return nullptr;
}

const VariableInfo* find_variable(const DatasetSchema& schema, const std::string_view name) noexcept {
  for (const auto& variable : schema.variables) {
    if (variable.name == name) {
      return &variable;
    }
  }
  return nullptr;
}

std::vector<SchemaIssue> validate_wrfinput_d01_schema(const DatasetSchema& schema) {
  std::vector<SchemaIssue> issues;
  require_dimensions(
      schema,
      {
          {"Time", 1, true},
          {"DateStrLen", 19, false},
          {"west_east", 265, false},
          {"south_north", 429, false},
          {"bottom_top", 59, false},
          {"bottom_top_stag", 60, false},
          {"west_east_stag", 266, false},
          {"south_north_stag", 430, false},
      },
      issues);
  require_variables(schema, wrfinput_core_variables(), issues);
  return issues;
}

std::vector<SchemaIssue> validate_wrfbdy_d01_schema(const DatasetSchema& schema) {
  std::vector<SchemaIssue> issues;
  require_dimensions(
      schema,
      {
          {"Time", 28, true},
          {"DateStrLen", 19, false},
          {"west_east", 265, false},
          {"south_north", 429, false},
          {"bottom_top", 59, false},
          {"bottom_top_stag", 60, false},
          {"west_east_stag", 266, false},
          {"south_north_stag", 430, false},
          {"bdy_width", 5, false},
      },
      issues);
  require_variables(schema, wrfbdy_core_variables(), issues);
  return issues;
}

std::vector<SchemaIssue> validate_wrffdda_d01_schema(const DatasetSchema& schema) {
  std::vector<SchemaIssue> issues;
  require_dimensions(
      schema,
      {
          {"Time", 28, true},
          {"DateStrLen", 19, false},
          {"west_east", 265, false},
          {"south_north", 429, false},
          {"bottom_top", 59, false},
          {"one_stag", 1, false},
      },
      issues);
  require_variables(schema, wrffdda_core_variables(), issues);
  return issues;
}

bool has_schema_errors(const std::vector<SchemaIssue>& issues) noexcept {
  for (const auto& issue : issues) {
    if (issue.severity == SchemaSeverity::error) {
      return true;
    }
  }
  return false;
}

}  // namespace tywrf::io
