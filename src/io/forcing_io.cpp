#include "tywrf/io/forcing_io.hpp"

#include <netcdf.h>

#include <algorithm>
#include <sstream>
#include <utility>

namespace tywrf::io {
namespace {

void check_netcdf(
    int status,
    const std::filesystem::path& path,
    std::string_view operation);

class NetcdfReadHandle {
 public:
  explicit NetcdfReadHandle(const std::filesystem::path& path) : path_(path) {
    check_netcdf(nc_open(path_.string().c_str(), NC_NOWRITE, &id_), path_, "open");
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

  void close() {
    if (id_ < 0) {
      return;
    }
    const int status = nc_close(id_);
    id_ = -1;
    check_netcdf(status, path_, "close");
  }

 private:
  std::filesystem::path path_;
  int id_ = -1;
};

void check_netcdf(
    const int status,
    const std::filesystem::path& path,
    const std::string_view operation) {
  if (status == NC_NOERR) {
    return;
  }
  std::ostringstream message;
  message << "NetCDF " << operation << " failed for " << path << ": "
          << nc_strerror(status);
  throw ForcingIoError(message.str());
}

[[nodiscard]] std::size_t product(const std::vector<std::size_t>& shape) {
  std::size_t count = 1;
  for (const auto extent : shape) {
    count *= extent;
  }
  return count;
}

[[nodiscard]] std::string trim_time_string(std::vector<char> buffer) {
  while (!buffer.empty() && (buffer.back() == '\0' || buffer.back() == ' ')) {
    buffer.pop_back();
  }
  return std::string(buffer.begin(), buffer.end());
}

void append_boundary_field_names(
    std::vector<std::string>& names,
    const std::string_view base_name) {
  for (const std::string_view suffix :
       {"_BXS", "_BXE", "_BYS", "_BYE", "_BTXS", "_BTXE", "_BTYS", "_BTYE"}) {
    names.push_back(std::string(base_name) + std::string(suffix));
  }
}

[[nodiscard]] const std::vector<std::string>& required_names_for(
    const KrosaForcingKind kind) {
  static const std::vector<std::string> wrfbdy = [] {
    std::vector<std::string> names;
    names.push_back("Times");
    for (const std::string_view base :
         {"U", "V", "W", "PH", "T", "QVAPOR", "QCLOUD", "QRAIN",
          "QICE", "QSNOW", "QGRAUP", "QNICE", "QNRAIN", "MU"}) {
      append_boundary_field_names(names, base);
    }
    return names;
  }();

  static const std::vector<std::string> wrffdda = {
      "Times",
      "U_NDG_OLD",
      "V_NDG_OLD",
      "T_NDG_OLD",
      "Q_NDG_OLD",
      "PH_NDG_OLD",
      "MU_NDG_OLD",
      "U_NDG_NEW",
      "V_NDG_NEW",
      "T_NDG_NEW",
      "Q_NDG_NEW",
      "PH_NDG_NEW",
      "MU_NDG_NEW",
  };

  return kind == KrosaForcingKind::boundary ? wrfbdy : wrffdda;
}

void require_time_index(
    const std::filesystem::path& path,
    const std::size_t time_index,
    const std::size_t time_count) {
  if (time_index < time_count) {
    return;
  }
  std::ostringstream message;
  message << "time index " << time_index << " is outside Time dimension length "
          << time_count << " for " << path;
  throw ForcingIoError(message.str());
}

[[nodiscard]] std::size_t find_dimension_length(
    const DatasetSchema& schema,
    const std::string_view name) {
  const auto* dimension = find_dimension(schema, name);
  if (dimension == nullptr) {
    return 0;
  }
  return dimension->size;
}

}  // namespace

std::vector<std::string> krosa_wrfbdy_required_variable_names() {
  return required_names_for(KrosaForcingKind::boundary);
}

std::vector<std::string> krosa_wrffdda_required_variable_names() {
  return required_names_for(KrosaForcingKind::spectral_nudging);
}

KrosaForcingReader::KrosaForcingReader(std::filesystem::path path, const KrosaForcingKind kind)
    : path_(std::move(path)), kind_(kind), schema_(read_netcdf_schema(path_)) {
  const auto* time = find_dimension(schema_, "Time");
  if (time == nullptr) {
    std::ostringstream message;
    message << "missing Time dimension in KROSA forcing file " << path_;
    throw ForcingIoError(message.str());
  }
  time_count_ = time->size;
  date_string_length_ = find_dimension_length(schema_, "DateStrLen");
}

const std::filesystem::path& KrosaForcingReader::path() const noexcept {
  return path_;
}

KrosaForcingKind KrosaForcingReader::kind() const noexcept {
  return kind_;
}

const DatasetSchema& KrosaForcingReader::schema() const noexcept {
  return schema_;
}

std::size_t KrosaForcingReader::time_count() const noexcept {
  return time_count_;
}

std::size_t KrosaForcingReader::date_string_length() const noexcept {
  return date_string_length_;
}

bool KrosaForcingReader::has_variable(const std::string_view name) const noexcept {
  return find_variable(schema_, name) != nullptr;
}

std::vector<std::string> KrosaForcingReader::missing_required_variables() const {
  std::vector<std::string> missing;
  for (const auto& name : required_names_for(kind_)) {
    if (!has_variable(name)) {
      missing.push_back(name);
    }
  }
  return missing;
}

ForcingVariableMetadata KrosaForcingReader::variable_metadata(const std::string_view name) const {
  const auto* variable = find_variable(schema_, name);
  if (variable == nullptr) {
    std::ostringstream message;
    message << "missing KROSA forcing variable " << name << " in " << path_;
    throw ForcingIoError(message.str());
  }
  if (variable->dimensions.empty() || variable->dimensions.front() != "Time") {
    std::ostringstream message;
    message << "KROSA forcing variable " << variable->name
            << " must use Time as its leading dimension in " << path_;
    throw ForcingIoError(message.str());
  }

  ForcingVariableMetadata metadata;
  metadata.name = variable->name;
  metadata.type_name = variable->type_name;
  metadata.dimensions = variable->dimensions;
  metadata.shape = variable->shape;
  metadata.slice_dimensions.assign(
      variable->dimensions.begin() + 1,
      variable->dimensions.end());
  metadata.slice_shape.assign(variable->shape.begin() + 1, variable->shape.end());
  metadata.time_count = time_count_;
  metadata.values_per_time_slice = product(metadata.slice_shape);
  return metadata;
}

std::string KrosaForcingReader::read_time_string(const std::size_t time_index) const {
  require_time_index(path_, time_index, time_count_);

  const auto metadata = variable_metadata("Times");
  if (metadata.type_name != "char" || metadata.shape.size() != 2) {
    std::ostringstream message;
    message << "Times must be a Time,DateStrLen char variable in " << path_;
    throw ForcingIoError(message.str());
  }

  const std::size_t length = metadata.shape.at(1);
  std::vector<char> buffer(length, '\0');
  NetcdfReadHandle file(path_);
  int var_id = -1;
  check_netcdf(nc_inq_varid(file.id(), "Times", &var_id), path_, "inquire Times");
  const std::size_t start[2] = {time_index, 0};
  const std::size_t count[2] = {1, length};
  check_netcdf(
      nc_get_vara_text(file.id(), var_id, start, count, buffer.data()),
      path_,
      "read Times time slice");
  file.close();
  return trim_time_string(std::move(buffer));
}

ForcingTimeSlice KrosaForcingReader::read_float_time_slice(
    const std::string_view name,
    const std::size_t time_index) const {
  require_time_index(path_, time_index, time_count_);

  auto metadata = variable_metadata(name);
  if (metadata.type_name != "float") {
    std::ostringstream message;
    message << "KROSA forcing variable " << metadata.name << " has type "
            << metadata.type_name << ", expected float in " << path_;
    throw ForcingIoError(message.str());
  }

  std::vector<std::size_t> start(metadata.shape.size(), 0);
  std::vector<std::size_t> count = metadata.shape;
  start.front() = time_index;
  count.front() = 1;

  ForcingTimeSlice slice;
  slice.time_index = time_index;
  slice.values.assign(metadata.values_per_time_slice, 0.0F);

  NetcdfReadHandle file(path_);
  int var_id = -1;
  const auto variable_name = std::string(name);
  check_netcdf(
      nc_inq_varid(file.id(), variable_name.c_str(), &var_id),
      path_,
      "inquire forcing variable");
  check_netcdf(
      nc_get_vara_float(file.id(), var_id, start.data(), count.data(), slice.values.data()),
      path_,
      "read float time slice");
  file.close();

  slice.metadata = std::move(metadata);
  return slice;
}

}  // namespace tywrf::io
