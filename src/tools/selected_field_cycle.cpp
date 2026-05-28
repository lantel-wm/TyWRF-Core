#include "tywrf/dynamics/base_state_provider.hpp"
#include "tywrf/dynamics/pressure_refresh_hook.hpp"
#include "tywrf/dynamics/wind_tendency.hpp"
#include "tywrf/io/pressure_refresh_io.hpp"
#include "tywrf/io/wrf_state_io.hpp"
#include "tywrf/nest/base_state_exchange.hpp"
#include "tywrf/nest/base_state_exchange_adapter.hpp"
#include "tywrf/nest/base_state_source_staging.hpp"
#include "tywrf/nest/parent_child_interpolation.hpp"
#include "tywrf/nest/static_fields.hpp"
#include "tywrf/nest/state_exchange.hpp"
#include "tywrf/nest/state_remap.hpp"
#include "tywrf/state.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

constexpr double kD02TargetDxMeters = 2000.0;
constexpr double kD02ResolutionToleranceMeters = 0.5;
constexpr std::size_t kMaxPressureColumnProbeColumns = 8;
constexpr std::size_t kMaxPressureColumnProbeLevels = 16;

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
  static const std::vector<std::string> names = {"U", "V", "T", "PH", "MU", "QVAPOR"};
  return names;
}

const std::vector<std::string>& pressure_column_probe_field_names() {
  static const std::vector<std::string> names = {
      "P",    "PB",  "P+PB", "MU", "MUB",    "MU+MUB",
      "PH",   "PHB", "PH+PHB", "T", "QVAPOR", "HGT"};
  return names;
}

const std::vector<std::string>& pressure_column_probe_unavailable_names() {
  static const std::vector<std::string> names = {
      "ALB", "C3F", "C4F", "C3H", "C4H", "P_TOP", "theta_m"};
  return names;
}

const std::vector<std::string>& pressure_formula_observation_field_names() {
  static const std::vector<std::string> names = {
      "status",
      "valid",
      "i",
      "j",
      "k",
      "mu_total",
      "pfu",
      "pfd",
      "phm",
      "log_ratio",
      "phi_lower",
      "phi_upper",
      "delta_phi",
      "ALB",
      "PB",
      "theta",
      "alpha_total",
      "alpha_perturbation",
      "alpha_from_wrf_branch",
      "pressure_base",
      "total_pressure",
      "perturbation_pressure_pa"};
  return names;
}

enum class WindTendencySourceKind {
  none,
  zero,
  identity,
  self_advection,
};

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
  bool pressure_refresh = false;
  bool experimental_pressure_refresh_apply = false;
  bool diagnostic_adapter_report = false;
  WindTendencySourceKind wind_tendency_source = WindTendencySourceKind::none;
  bool pretty = false;
  std::vector<std::string> variables;
  std::vector<std::pair<std::int32_t, std::int32_t>> pressure_column_probe_columns;
  std::vector<std::int32_t> pressure_column_probe_levels;
  bool has_pressure_column_probe_levels = false;
};

struct Resolution {
  double dx = 0.0;
  double dy = 0.0;
};

struct RuntimeTimelineEvent {
  std::string name;
  std::vector<std::pair<std::string, std::string>> fields;
};

struct PressureColumnObservation {
  std::string phase;
  std::int32_t i = 0;
  std::int32_t j = 0;
  std::int32_t k = 0;
  double p = 0.0;
  double pb = 0.0;
  double p_plus_pb = 0.0;
  double mu = 0.0;
  double mub = 0.0;
  double mu_plus_mub = 0.0;
  double ph = 0.0;
  double phb = 0.0;
  double ph_plus_phb = 0.0;
  double t = 0.0;
  double qvapor = 0.0;
  double hgt = 0.0;
};

struct DiagnosticAdapterSourceChildFieldDeltaReport {
  std::uint64_t compared_value_count = 0;
  std::uint64_t differing_value_count = 0;
  double max_abs_diff = 0.0;
};

struct DiagnosticAdapterSourceChildDeltaReport {
  bool diagnostic_only = true;
  bool gate_candidate = false;
  bool integrator_output = false;
  bool writes_candidate = false;
  bool writes_netcdf = false;
  bool values_identical = false;
  std::uint64_t compared_value_count = 0;
  std::uint64_t differing_value_count = 0;
  double max_abs_diff = 0.0;
  DiagnosticAdapterSourceChildFieldDeltaReport phb;
  DiagnosticAdapterSourceChildFieldDeltaReport mub;
  DiagnosticAdapterSourceChildFieldDeltaReport ht;
  DiagnosticAdapterSourceChildFieldDeltaReport pb;
  DiagnosticAdapterSourceChildFieldDeltaReport t_init;
  DiagnosticAdapterSourceChildFieldDeltaReport alb;
};

struct WindTendencyOptInReport {
  WindTendencySourceKind source_kind = WindTendencySourceKind::none;
  tywrf::dynamics::WindTendencyReport kernel;
  std::uint64_t changed_u_points = 0;
  std::uint64_t changed_v_points = 0;
};

struct CandidateReport {
  std::uint64_t changed_selected_points = 0;
  std::uint64_t changed_static_template_points = 0;
  std::vector<RuntimeTimelineEvent> timeline;
  tywrf::nest::ChildStateRemapReport remap;
  tywrf::nest::RemapPlan remap_plan;
  tywrf::nest::StateExchangePlan exchange;
  tywrf::nest::ParentChildInterpolationReport interpolation;
  tywrf::nest::MovingNestStaticRefreshReport static_refresh;
  std::optional<tywrf::dynamics::KrosaBaseStateProviderReport>
      normal_base_state_provider_source;
  std::optional<tywrf::nest::NormalCandidateBaseStateExchangeReport>
      normal_base_state_producer;
  std::optional<tywrf::dynamics::KrosaBaseStateProviderReport>
      pressure_refresh_provider_probe;
  std::optional<tywrf::dynamics::KrosaPressureRefreshHookReport>
      pressure_refresh_dry_run_contract;
  std::optional<tywrf::dynamics::KrosaPressureRefreshHookReport> pressure_refresh;
  std::optional<tywrf::nest::ExposedBaseStateExchangeAdapterReport>
      diagnostic_adapter;
  std::optional<tywrf::dynamics::KrosaBaseStateProviderReport>
      diagnostic_adapter_provider_source;
  std::optional<tywrf::nest::BaseStateSourceStagingReport>
      diagnostic_adapter_source_staging;
  bool diagnostic_adapter_source_staging_aliases_child = false;
  std::optional<DiagnosticAdapterSourceChildDeltaReport>
      diagnostic_adapter_source_child_delta;
  std::optional<WindTendencyOptInReport> wind_tendency;
  std::filesystem::path diagnostic_adapter_metadata_source;
  std::size_t diagnostic_adapter_metadata_time_index = 0;
  std::filesystem::path pressure_refresh_metadata_source;
  std::size_t pressure_refresh_metadata_time_index = 0;
  std::filesystem::path normal_base_state_metadata_source;
  std::size_t normal_base_state_metadata_time_index = 0;
  std::uint64_t normal_base_state_changed_p_points = 0;
  std::uint64_t normal_base_state_changed_pb_points = 0;
  std::uint64_t normal_base_state_changed_mub_points = 0;
  std::uint64_t normal_base_state_changed_phb_points = 0;
  std::uint64_t normal_base_state_changed_hgt_points = 0;
  std::uint64_t pressure_refresh_changed_p_points = 0;
  std::uint64_t pressure_refresh_changed_pb_points = 0;
  std::uint64_t pressure_refresh_changed_mub_points = 0;
  std::uint64_t pressure_refresh_changed_phb_points = 0;
  double cen_lat = 0.0;
  double cen_lon = 0.0;
  std::vector<PressureColumnObservation> pressure_column_observations;
  std::vector<tywrf::dynamics::PressureRefreshFormulaObservation>
      pressure_formula_observations;
};

struct PressureRefreshReportParity {
  bool changed_p_matches_refreshed_point_count = false;
  bool invalid_and_skipped_points_zero = false;
  bool overlap_halo_untouched = false;
};

struct StaticFieldSet {
  explicit StaticFieldSet(const tywrf::Grid& grid)
      : xlat(grid.surface_layout()), xlong(grid.surface_layout()), hgt(grid.surface_layout()) {}

  tywrf::FieldStorage2D<float> xlat;
  tywrf::FieldStorage2D<float> xlong;
  tywrf::FieldStorage2D<float> hgt;
};

struct DiagnosticBaseStateAdapterStaging {
  explicit DiagnosticBaseStateAdapterStaging(const tywrf::Grid& grid)
      : phb(grid.w_layout()),
        mub(grid.surface_layout()),
        pb(grid.mass_layout()),
        t_init(grid.mass_layout()),
        alb(grid.mass_layout()),
        ht(grid.surface_layout()) {}

  tywrf::FieldStorage3D<float> phb;
  tywrf::FieldStorage2D<float> mub;
  tywrf::FieldStorage3D<float> pb;
  tywrf::FieldStorage3D<float> t_init;
  tywrf::FieldStorage3D<float> alb;
  tywrf::FieldStorage2D<float> ht;
};

struct DiagnosticAdapterProviderSource {
  tywrf::dynamics::KrosaBaseStateProvider provider;
  tywrf::dynamics::KrosaBaseStateProviderReport report;
};

[[nodiscard]] const char* bool_text(const bool value) noexcept {
  return value ? "true" : "false";
}

[[nodiscard]] constexpr bool wind_tendency_enabled(
    const WindTendencySourceKind kind) noexcept {
  return kind != WindTendencySourceKind::none;
}

[[nodiscard]] std::string_view wind_tendency_source_name(
    const WindTendencySourceKind kind) noexcept {
  switch (kind) {
    case WindTendencySourceKind::none:
      return "none";
    case WindTendencySourceKind::zero:
      return "zero";
    case WindTendencySourceKind::identity:
      return "identity";
    case WindTendencySourceKind::self_advection:
      return "self_advection";
  }
  return "none";
}

[[nodiscard]] constexpr bool wind_tendency_gate_evidence(
    const WindTendencySourceKind kind) noexcept {
  return kind == WindTendencySourceKind::self_advection;
}

[[nodiscard]] constexpr bool wind_tendency_zero_or_identity_only(
    const WindTendencySourceKind kind) noexcept {
  return kind == WindTendencySourceKind::zero || kind == WindTendencySourceKind::identity;
}

[[nodiscard]] WindTendencySourceKind parse_wind_tendency_source(
    const std::string& value,
    const std::string_view option) {
  auto trimmed = value;
  trimmed.erase(
      trimmed.begin(),
      std::find_if(trimmed.begin(), trimmed.end(), [](const unsigned char c) {
        return !std::isspace(c);
      }));
  trimmed.erase(
      std::find_if(trimmed.rbegin(), trimmed.rend(), [](const unsigned char c) {
        return !std::isspace(c);
      }).base(),
      trimmed.end());
  if (trimmed == "none") {
    return WindTendencySourceKind::none;
  }
  if (trimmed == "zero") {
    return WindTendencySourceKind::zero;
  }
  if (trimmed == "identity") {
    return WindTendencySourceKind::identity;
  }
  if (trimmed == "self-advection" || trimmed == "self_advection") {
    return WindTendencySourceKind::self_advection;
  }
  throw std::invalid_argument(
      std::string(option) + " expects one of: none, zero, identity, self-advection");
}

[[nodiscard]] std::string_view wind_tendency_status_name(
    const tywrf::dynamics::WindTendencyStatus status) noexcept {
  switch (status) {
    case tywrf::dynamics::WindTendencyStatus::ok:
      return "ok";
    case tywrf::dynamics::WindTendencyStatus::null_target:
      return "null_target";
    case tywrf::dynamics::WindTendencyStatus::null_source:
      return "null_source";
    case tywrf::dynamics::WindTendencyStatus::invalid_config:
      return "invalid_config";
    case tywrf::dynamics::WindTendencyStatus::invalid_layout:
      return "invalid_layout";
    case tywrf::dynamics::WindTendencyStatus::insufficient_halo:
      return "insufficient_halo";
    case tywrf::dynamics::WindTendencyStatus::mismatched_source_layout:
      return "mismatched_source_layout";
    case tywrf::dynamics::WindTendencyStatus::mismatched_wind_layout:
      return "mismatched_wind_layout";
  }
  return "unknown";
}

template <typename Storage>
void copy_storage(
    const Storage& source,
    Storage& destination,
    const std::string_view label) {
  if (source.size() != destination.size()) {
    throw std::runtime_error(std::string(label) + " staging shape mismatch");
  }
  std::copy(source.data(), source.data() + source.size(), destination.data());
}

[[nodiscard]] tywrf::nest::ExposedBaseStateViews<const float>
diagnostic_adapter_views(const DiagnosticBaseStateAdapterStaging& staging) {
  return {
      staging.phb.view(),
      staging.mub.view(),
      staging.pb.view(),
      staging.t_init.view(),
      staging.alb.view(),
      staging.ht.view()};
}

[[nodiscard]] tywrf::nest::ExposedBaseStateViews<float>
diagnostic_adapter_views(DiagnosticBaseStateAdapterStaging& staging) {
  return {
      staging.phb.view(),
      staging.mub.view(),
      staging.pb.view(),
      staging.t_init.view(),
      staging.alb.view(),
      staging.ht.view()};
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] bool field_data_alias(
    const tywrf::FieldView2D<LhsReal>& lhs,
    const tywrf::FieldView2D<RhsReal>& rhs) noexcept {
  return lhs.data != nullptr && rhs.data != nullptr && lhs.data == rhs.data;
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] bool field_data_alias(
    const tywrf::FieldView3D<LhsReal>& lhs,
    const tywrf::FieldView3D<RhsReal>& rhs) noexcept {
  return lhs.data != nullptr && rhs.data != nullptr && lhs.data == rhs.data;
}

[[nodiscard]] bool base_state_views_alias(
    const tywrf::nest::ExposedBaseStateViews<const float>& source,
    const tywrf::nest::ExposedBaseStateViews<float>& child) noexcept {
  return field_data_alias(source.phb, child.phb) ||
         field_data_alias(source.mub, child.mub) ||
         field_data_alias(source.pb, child.pb) ||
         field_data_alias(source.t_init, child.t_init) ||
         field_data_alias(source.alb, child.alb) ||
         field_data_alias(source.ht, child.ht);
}

[[nodiscard]] std::string timeline_value(std::string value) {
  for (char& c : value) {
    if (c == '|' || c == ';' || c == '(' || c == ')' || c == ',' || c == '=') {
      c = '_';
    }
  }
  return value;
}

[[nodiscard]] std::pair<std::string, std::string> timeline_field(
    std::string name,
    std::string value) {
  return {std::move(name), timeline_value(std::move(value))};
}

[[nodiscard]] std::pair<std::string, std::string> timeline_field(
    std::string name,
    const std::int64_t value) {
  return {std::move(name), std::to_string(value)};
}

[[nodiscard]] std::pair<std::string, std::string> timeline_field(
    std::string name,
    const std::uint64_t value) {
  return {std::move(name), std::to_string(value)};
}

void append_timeline_event(
    CandidateReport& report,
    std::string name,
    std::vector<std::pair<std::string, std::string>> fields) {
  report.timeline.push_back({std::move(name), std::move(fields)});
}

[[nodiscard]] std::string join_timeline_event_names(
    const std::vector<RuntimeTimelineEvent>& events) {
  std::ostringstream joined;
  for (std::size_t index = 0; index < events.size(); ++index) {
    if (index != 0) {
      joined << ",";
    }
    joined << events[index].name;
  }
  return joined.str();
}

[[nodiscard]] std::string join_timeline_events(
    const std::vector<RuntimeTimelineEvent>& events) {
  std::ostringstream joined;
  for (std::size_t event_index = 0; event_index < events.size(); ++event_index) {
    if (event_index != 0) {
      joined << "|";
    }
    const auto& event = events[event_index];
    joined << (event_index + 1U) << ":" << event.name << "(";
    for (std::size_t field_index = 0; field_index < event.fields.size(); ++field_index) {
      if (field_index != 0) {
        joined << ",";
      }
      joined << event.fields[field_index].first << "=" << event.fields[field_index].second;
    }
    joined << ")";
  }
  return joined.str();
}

[[nodiscard]] tywrf::dynamics::KrosaBaseStateProviderTerrainOverride
make_moved_candidate_hgt_terrain_override(const StaticFieldSet& output_static) {
  return {
      .terrain_height_m = output_static.hgt.view(),
      .source_name = "moved_candidate_HGT",
      .provenance = "override:moved_candidate_HGT"};
}

[[nodiscard]] constexpr std::string_view diagnostic_adapter_provider_source_origin() noexcept {
  return "base_state_reconstruction_provider+moved_candidate_HGT";
}

[[nodiscard]] PressureRefreshReportParity pressure_refresh_report_parity(
    const CandidateReport& candidate_report,
    const tywrf::dynamics::KrosaPressureRefreshHookReport& pressure_report) {
  return {
      .changed_p_matches_refreshed_point_count =
          candidate_report.pressure_refresh_changed_p_points ==
          pressure_report.compute_report.refreshed_point_count,
      .invalid_and_skipped_points_zero =
          pressure_report.compute_report.invalid_point_count == 0 &&
          pressure_report.compute_report.skipped_point_count == 0,
      .overlap_halo_untouched =
          !pressure_report.compute_report.touched_overlap_cells &&
          !pressure_report.compute_report.touched_halo_cells};
}

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
  --wind-tendency-source KIND    Opt in to selected-field U/V horizontal wind tendency plumbing.
                                  KIND is none, zero, identity, self-advection, or self_advection;
                                  default none. zero and identity are placeholder sources only and
                                  are not validation-gate evidence.
  --pressure-refresh             Opt in to provider-backed KROSA pressure refresh readiness check.
                                  Current selected-field state aborts before output because
                                  PB/PHB/MUB/P ownership is not yet thermodynamically consistent
                                  with the moved-pose T/PH and terrain producer.
  --pressure-column-probe I,J[;I,J...]
                                  Opt in to runtime same-column pressure telemetry for zero-based
                                  d02 mass-grid columns. Bounded to 8 columns.
  --pressure-column-levels K[,K...]
                                  Zero-based mass levels for --pressure-column-probe; default is
                                  levels 0..4 capped to the active mass-level count.
  --pretty                       Pretty-print JSON report.
  --help                         Show this help.

Oracle inputs are not accepted: --end-state, --reference-end, --d01-end-state,
and --d02-end-state are rejected. This tool remaps the d02 start state from the
old moving-nest pose to the new pose, fills newly exposed U/V/MU/QVAPOR cells
from d01 start-state parent interpolation, fills newly exposed T/PH from the
same d01 start-state parent interpolation, and preserves P plus available d02
start-state diagnostic/base-state fields without using reference-end truth. It
is a selected-field integrator candidate, not a WRF-exact physics result.
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

[[nodiscard]] std::string trim_copy(std::string token) {
  token.erase(token.begin(), std::find_if(token.begin(), token.end(), [](const unsigned char c) {
                return !std::isspace(c);
              }));
  token.erase(
      std::find_if(token.rbegin(), token.rend(), [](const unsigned char c) {
        return !std::isspace(c);
      }).base(),
      token.end());
  return token;
}

[[nodiscard]] std::int32_t parse_nonnegative_int32(
    const std::string& value,
    const std::string_view option) {
  const auto trimmed = trim_copy(value);
  if (trimmed.empty()) {
    throw std::invalid_argument(std::string(option) + " expects a non-negative integer");
  }
  std::size_t parsed_chars = 0;
  const auto parsed = std::stoll(trimmed, &parsed_chars);
  if (parsed_chars != trimmed.size() || parsed < 0 ||
      parsed > static_cast<long long>(std::numeric_limits<std::int32_t>::max())) {
    throw std::invalid_argument(std::string(option) + " expects a non-negative integer");
  }
  return static_cast<std::int32_t>(parsed);
}

[[nodiscard]] std::vector<std::string> split_variables(const std::string& value) {
  std::vector<std::string> variables;
  std::size_t begin = 0;
  while (begin <= value.size()) {
    const auto end = value.find(',', begin);
    auto token = trim_copy(
        value.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
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

[[nodiscard]] std::vector<std::pair<std::int32_t, std::int32_t>> parse_pressure_probe_columns(
    const std::string& value,
    const std::string_view option) {
  std::vector<std::pair<std::int32_t, std::int32_t>> columns;
  std::size_t begin = 0;
  while (begin <= value.size()) {
    const auto end = value.find(';', begin);
    const auto token = trim_copy(
        value.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
    if (token.empty()) {
      throw std::invalid_argument(
          std::string(option) + " must contain zero-based I,J column pairs");
    }
    const auto separator = token.find(',');
    if (separator == std::string::npos || token.find(',', separator + 1) != std::string::npos) {
      throw std::invalid_argument(
          std::string(option) + " must contain zero-based I,J column pairs");
    }
    columns.push_back({
        parse_nonnegative_int32(token.substr(0, separator), option),
        parse_nonnegative_int32(token.substr(separator + 1), option)});
    if (columns.size() > kMaxPressureColumnProbeColumns) {
      throw std::invalid_argument("--pressure-column-probe is limited to 8 columns");
    }
    if (end == std::string::npos) {
      break;
    }
    begin = end + 1;
  }
  if (columns.empty()) {
    throw std::invalid_argument(std::string(option) + " requires at least one column");
  }
  return columns;
}

[[nodiscard]] std::vector<std::int32_t> parse_pressure_probe_levels(
    const std::string& value,
    const std::string_view option) {
  std::vector<std::int32_t> levels;
  std::size_t begin = 0;
  while (begin <= value.size()) {
    const auto end = value.find(',', begin);
    const auto token = trim_copy(
        value.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
    if (token.empty()) {
      throw std::invalid_argument(
          std::string(option) + " must contain zero-based mass levels");
    }
    levels.push_back(parse_nonnegative_int32(token, option));
    if (levels.size() > kMaxPressureColumnProbeLevels) {
      throw std::invalid_argument("--pressure-column-levels is limited to 16 levels");
    }
    if (end == std::string::npos) {
      break;
    }
    begin = end + 1;
  }
  if (levels.empty()) {
    throw std::invalid_argument(std::string(option) + " requires at least one level");
  }
  return levels;
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
    } else if (arg == "--wind-tendency-source") {
      options.wind_tendency_source = parse_wind_tendency_source(require_value(arg), arg);
    } else if (arg.rfind("--wind-tendency-source=", 0) == 0) {
      options.wind_tendency_source =
          parse_wind_tendency_source(arg.substr(std::string("--wind-tendency-source=").size()), "--wind-tendency-source");
    } else if (arg == "--pressure-refresh") {
      options.pressure_refresh = true;
    } else if (arg == "--diagnostic-base-state-adapter-report") {
      options.diagnostic_adapter_report = true;
    } else if (arg == "--pressure-column-probe") {
      options.pressure_column_probe_columns =
          parse_pressure_probe_columns(require_value(arg), arg);
    } else if (arg == "--pressure-column-levels") {
      options.pressure_column_probe_levels =
          parse_pressure_probe_levels(require_value(arg), arg);
      options.has_pressure_column_probe_levels = true;
    } else if (arg == "--experimental-pressure-refresh-apply") {
      options.experimental_pressure_refresh_apply = true;
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
  if (options.experimental_pressure_refresh_apply && !options.pressure_refresh) {
    throw std::invalid_argument(
        "--experimental-pressure-refresh-apply requires --pressure-refresh");
  }
  if (options.diagnostic_adapter_report && options.pressure_refresh) {
    throw std::invalid_argument(
        "--diagnostic-base-state-adapter-report cannot be combined with --pressure-refresh");
  }
  if (options.has_pressure_column_probe_levels &&
      options.pressure_column_probe_columns.empty()) {
    throw std::invalid_argument(
        "--pressure-column-levels requires --pressure-column-probe");
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

void load_static_2d_variable(
    const std::filesystem::path& path,
    const std::string_view name,
    tywrf::FieldStorage2D<float>& field,
    const std::size_t time_index) {
  const NetcdfHandle file(path, NetcdfHandle::Mode::read);
  int variable_id = -1;
  file.check(nc_inq_varid(file.id(), std::string(name).c_str(), &variable_id), "inquire static variable");
  const auto layout = field.layout();
  const std::array<std::size_t, 3> start = {time_index, 0, 0};
  const std::array<std::size_t, 3> count = {
      1,
      static_cast<std::size_t>(layout.active_ny()),
      static_cast<std::size_t>(layout.active_nx())};
  std::vector<float> buffer(count[1] * count[2]);
  file.check(
      nc_get_vara_float(file.id(), variable_id, start.data(), count.data(), buffer.data()),
      "read static variable");

  auto view = field.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
      view(layout.i_begin() + i, layout.j_begin() + j) =
          buffer[static_cast<std::size_t>(j) * count[2] + static_cast<std::size_t>(i)];
    }
  }
}

[[nodiscard]] StaticFieldSet load_static_fields(
    const std::filesystem::path& path,
    const tywrf::Grid& grid,
    const std::size_t time_index) {
  StaticFieldSet fields(grid);
  load_static_2d_variable(path, "XLAT", fields.xlat, time_index);
  load_static_2d_variable(path, "XLONG", fields.xlong, time_index);
  load_static_2d_variable(path, "HGT", fields.hgt, time_index);
  return fields;
}

[[nodiscard]] tywrf::FieldStorage2D<float> load_hgt_field(
    const std::filesystem::path& path,
    const tywrf::Grid& grid,
    const std::size_t time_index) {
  tywrf::FieldStorage2D<float> hgt(grid.surface_layout());
  load_static_2d_variable(path, "HGT", hgt, time_index);
  return hgt;
}

void write_static_2d_variable(
    const NetcdfHandle& file,
    const std::string_view name,
    const tywrf::FieldStorage2D<float>& field,
    const std::size_t time_index) {
  int variable_id = -1;
  file.check(nc_inq_varid(file.id(), std::string(name).c_str(), &variable_id), "inquire output static variable");
  const auto layout = field.layout();
  const std::array<std::size_t, 3> start = {time_index, 0, 0};
  const std::array<std::size_t, 3> count = {
      1,
      static_cast<std::size_t>(layout.active_ny()),
      static_cast<std::size_t>(layout.active_nx())};
  std::vector<float> buffer(count[1] * count[2]);
  const auto view = field.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
      buffer[static_cast<std::size_t>(j) * count[2] + static_cast<std::size_t>(i)] =
          view(layout.i_begin() + i, layout.j_begin() + j);
    }
  }
  file.check(
      nc_put_vara_float(file.id(), variable_id, start.data(), count.data(), buffer.data()),
      "write output static variable");
}

void overwrite_output_static_fields(
    const std::filesystem::path& output_path,
    const StaticFieldSet& fields,
    const std::size_t time_index) {
  const NetcdfHandle file(output_path, NetcdfHandle::Mode::write);
  write_static_2d_variable(file, "XLAT", fields.xlat, time_index);
  write_static_2d_variable(file, "XLONG", fields.xlong, time_index);
  write_static_2d_variable(file, "HGT", fields.hgt, time_index);
}

template <typename BeforeStorage, typename AfterStorage>
[[nodiscard]] std::uint64_t changed_points(
    const BeforeStorage& before,
    const AfterStorage& after) {
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
         changed_points(before.t, after.t) + changed_points(before.ph, after.ph) +
         changed_points(before.mu, after.mu) + changed_points(before.qvapor, after.qvapor);
}

template <typename Real>
[[nodiscard]] tywrf::ActiveShape3D active_shape(const tywrf::FieldView3D<Real> view) noexcept {
  return {
      view.nx - view.halo.i_lower - view.halo.i_upper,
      view.ny - view.halo.j_lower - view.halo.j_upper,
      view.nz - view.halo.k_lower - view.halo.k_upper};
}

template <typename Real>
void fill_halo_staging_from_active_clamped(
    const tywrf::FieldView3D<const Real> source,
    tywrf::FieldStorage3D<Real>& destination) {
  const auto destination_view = destination.view();
  const auto source_i_begin = source.halo.i_lower;
  const auto source_i_end = source.nx - source.halo.i_upper;
  const auto source_j_begin = source.halo.j_lower;
  const auto source_j_end = source.ny - source.halo.j_upper;
  const auto source_k_begin = source.halo.k_lower;
  const auto source_k_end = source.nz - source.halo.k_upper;
  for (std::int32_t j = 0; j < destination_view.ny; ++j) {
    const auto source_j = std::clamp(
        j - destination_view.halo.j_lower + source_j_begin,
        source_j_begin,
        source_j_end - 1);
    for (std::int32_t k = 0; k < destination_view.nz; ++k) {
      const auto source_k = std::clamp(
          k - destination_view.halo.k_lower + source_k_begin,
          source_k_begin,
          source_k_end - 1);
      for (std::int32_t i = 0; i < destination_view.nx; ++i) {
        const auto source_i = std::clamp(
            i - destination_view.halo.i_lower + source_i_begin,
            source_i_begin,
            source_i_end - 1);
        destination_view(i, j, k) = source(source_i, source_j, source_k);
      }
    }
  }
}

template <typename Real>
void fill_storage_constant(
    tywrf::FieldStorage3D<Real>& destination,
    const Real value) {
  std::fill(destination.data(), destination.data() + destination.size(), value);
}

template <typename Real>
void copy_staging_active_to_field(
    const tywrf::FieldView3D<const Real> source,
    const tywrf::FieldView3D<Real> destination) {
  const auto destination_i_begin = destination.halo.i_lower;
  const auto destination_j_begin = destination.halo.j_lower;
  const auto destination_k_begin = destination.halo.k_lower;
  for (std::int32_t j = source.halo.j_lower; j < source.ny - source.halo.j_upper; ++j) {
    const auto destination_j = j - source.halo.j_lower + destination_j_begin;
    for (std::int32_t k = source.halo.k_lower; k < source.nz - source.halo.k_upper; ++k) {
      const auto destination_k = k - source.halo.k_lower + destination_k_begin;
      for (std::int32_t i = source.halo.i_lower; i < source.nx - source.halo.i_upper; ++i) {
        const auto destination_i = i - source.halo.i_lower + destination_i_begin;
        destination(destination_i, destination_j, destination_k) = source(i, j, k);
      }
    }
  }
}

void apply_selected_field_wind_tendency(
    const Options& options,
    tywrf::State<float>& candidate,
    CandidateReport& report) {
  if (!wind_tendency_enabled(options.wind_tendency_source)) {
    return;
  }

  auto u_target =
      tywrf::FieldStorage3D<float>(active_shape(candidate.u.view()), tywrf::uniform_halo_3d(1));
  auto v_target =
      tywrf::FieldStorage3D<float>(active_shape(candidate.v.view()), tywrf::uniform_halo_3d(1));
  auto u_source =
      tywrf::FieldStorage3D<float>(active_shape(candidate.u.view()), tywrf::uniform_halo_3d(1));
  auto v_source =
      tywrf::FieldStorage3D<float>(active_shape(candidate.v.view()), tywrf::uniform_halo_3d(1));
  auto u_advect_x =
      tywrf::FieldStorage3D<float>(active_shape(candidate.u.view()), tywrf::uniform_halo_3d(1));
  auto u_advect_y =
      tywrf::FieldStorage3D<float>(active_shape(candidate.u.view()), tywrf::uniform_halo_3d(1));
  auto v_advect_x =
      tywrf::FieldStorage3D<float>(active_shape(candidate.v.view()), tywrf::uniform_halo_3d(1));
  auto v_advect_y =
      tywrf::FieldStorage3D<float>(active_shape(candidate.v.view()), tywrf::uniform_halo_3d(1));

  const tywrf::State<float>& candidate_source = candidate;
  fill_halo_staging_from_active_clamped(candidate_source.u.view(), u_target);
  fill_halo_staging_from_active_clamped(candidate_source.v.view(), v_target);
  if (options.wind_tendency_source == WindTendencySourceKind::zero) {
    fill_storage_constant(u_source, 0.0F);
    fill_storage_constant(v_source, 0.0F);
    fill_storage_constant(u_advect_x, 0.0F);
    fill_storage_constant(u_advect_y, 0.0F);
    fill_storage_constant(v_advect_x, 0.0F);
    fill_storage_constant(v_advect_y, 0.0F);
  } else {
    fill_halo_staging_from_active_clamped(candidate_source.u.view(), u_source);
    fill_halo_staging_from_active_clamped(candidate_source.v.view(), v_source);
    fill_halo_staging_from_active_clamped(candidate_source.u.view(), u_advect_x);
    fill_halo_staging_from_active_clamped(candidate_source.u.view(), u_advect_y);
    fill_halo_staging_from_active_clamped(candidate_source.v.view(), v_advect_x);
    fill_halo_staging_from_active_clamped(candidate_source.v.view(), v_advect_y);
  }

  const auto u_before = candidate.u;
  const auto v_before = candidate.v;
  const auto& u_target_const = static_cast<const tywrf::FieldStorage3D<float>&>(u_target);
  const auto& v_target_const = static_cast<const tywrf::FieldStorage3D<float>&>(v_target);
  const auto& u_source_const = static_cast<const tywrf::FieldStorage3D<float>&>(u_source);
  const auto& v_source_const = static_cast<const tywrf::FieldStorage3D<float>&>(v_source);
  const auto& u_advect_x_const =
      static_cast<const tywrf::FieldStorage3D<float>&>(u_advect_x);
  const auto& u_advect_y_const =
      static_cast<const tywrf::FieldStorage3D<float>&>(u_advect_y);
  const auto& v_advect_x_const =
      static_cast<const tywrf::FieldStorage3D<float>&>(v_advect_x);
  const auto& v_advect_y_const =
      static_cast<const tywrf::FieldStorage3D<float>&>(v_advect_y);
  WindTendencyOptInReport wind_report;
  wind_report.source_kind = options.wind_tendency_source;
  const bool gate_evidence = wind_tendency_gate_evidence(options.wind_tendency_source);
  wind_report.kernel = tywrf::dynamics::apply_horizontal_wind_tendency(
      tywrf::dynamics::WindTendencyViews<float>{
          {u_target.view(),
           u_source_const.view(),
           u_advect_x_const.view(),
           u_advect_y_const.view()},
          {v_target.view(),
           v_source_const.view(),
           v_advect_x_const.view(),
           v_advect_y_const.view()}},
      tywrf::dynamics::WindTendencyConfig<float>{
          .dt_seconds = 8.0F,
          .dx_m = static_cast<float>(kD02TargetDxMeters),
          .dy_m = static_cast<float>(kD02TargetDxMeters),
          .enable_horizontal_advection = true,
          .diagnostic_only = false,
          .gate_candidate = gate_evidence,
          .validation_gate_evidence = gate_evidence});
  if (wind_report.kernel.status != tywrf::dynamics::WindTendencyStatus::ok) {
    throw std::runtime_error(
        "selected-field wind tendency failed: " +
        std::string(wind_tendency_status_name(wind_report.kernel.status)));
  }

  copy_staging_active_to_field(u_target_const.view(), candidate.u.view());
  copy_staging_active_to_field(v_target_const.view(), candidate.v.view());
  wind_report.changed_u_points = changed_points(u_before, candidate.u);
  wind_report.changed_v_points = changed_points(v_before, candidate.v);
  report.wind_tendency = wind_report;
  append_timeline_event(
      report,
      "wind_tendency_apply",
      {
          timeline_field("opt_in", "true"),
          timeline_field("applied", "true"),
          timeline_field(
              "source_kind",
              std::string(wind_tendency_source_name(options.wind_tendency_source))),
          timeline_field("fields", "U_V"),
          timeline_field(
              "zero_or_identity_only",
              wind_tendency_zero_or_identity_only(options.wind_tendency_source) ? "true"
                                                                                : "false"),
          timeline_field("gate_evidence", gate_evidence ? "true" : "false"),
          timeline_field("validation_gate_evidence", gate_evidence ? "true" : "false"),
          timeline_field("uses_reference_end_truth", "false"),
          timeline_field("changed_u_points", wind_report.changed_u_points),
          timeline_field("changed_v_points", wind_report.changed_v_points),
      });
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

void require_finite_static_fields(const StaticFieldSet& fields) {
  require_finite_storage("XLAT", fields.xlat);
  require_finite_storage("XLONG", fields.xlong);
  require_finite_storage("HGT", fields.hgt);
}

[[nodiscard]] std::uint64_t changed_static_points(
    const StaticFieldSet& before,
    const StaticFieldSet& after) {
  return changed_points(before.xlat, after.xlat) + changed_points(before.xlong, after.xlong) +
         changed_points(before.hgt, after.hgt);
}

[[nodiscard]] std::string join_variables(const std::vector<std::string>& variables);

[[nodiscard]] bool pressure_column_probe_enabled(const Options& options) {
  return !options.pressure_column_probe_columns.empty();
}

[[nodiscard]] bool pressure_formula_observation_enabled(
    const Options& options) {
  return options.pressure_refresh && pressure_column_probe_enabled(options);
}

void finalize_pressure_column_probe_options(Options& options, const tywrf::Grid& grid) {
  if (!pressure_column_probe_enabled(options)) {
    return;
  }

  const auto mass_layout = grid.mass_layout();
  if (mass_layout.active_nx() <= 0 || mass_layout.active_ny() <= 0 ||
      mass_layout.active_nz() <= 0) {
    throw std::runtime_error("pressure column probe requires a positive d02 mass-grid shape");
  }

  if (!options.has_pressure_column_probe_levels) {
    const auto default_count = std::min<std::int32_t>(5, mass_layout.active_nz());
    options.pressure_column_probe_levels.clear();
    for (std::int32_t k = 0; k < default_count; ++k) {
      options.pressure_column_probe_levels.push_back(k);
    }
  }

  for (const auto [i, j] : options.pressure_column_probe_columns) {
    if (i >= mass_layout.active_nx() || j >= mass_layout.active_ny()) {
      std::ostringstream message;
      message << "--pressure-column-probe column " << i << "," << j
              << " is outside zero-based d02 mass-grid shape "
              << mass_layout.active_nx() << "x" << mass_layout.active_ny();
      throw std::invalid_argument(message.str());
    }
  }
  for (const auto k : options.pressure_column_probe_levels) {
    if (k >= mass_layout.active_nz()) {
      std::ostringstream message;
      message << "--pressure-column-levels level " << k
              << " is outside zero-based d02 mass-level count "
              << mass_layout.active_nz();
      throw std::invalid_argument(message.str());
    }
  }
}

[[nodiscard]] std::string join_pressure_probe_columns(
    const std::vector<std::pair<std::int32_t, std::int32_t>>& columns) {
  std::ostringstream joined;
  for (std::size_t index = 0; index < columns.size(); ++index) {
    if (index != 0) {
      joined << ";";
    }
    joined << columns[index].first << "," << columns[index].second;
  }
  return joined.str();
}

[[nodiscard]] std::string join_pressure_probe_levels(
    const std::vector<std::int32_t>& levels) {
  std::ostringstream joined;
  for (std::size_t index = 0; index < levels.size(); ++index) {
    if (index != 0) {
      joined << ",";
    }
    joined << levels[index];
  }
  return joined.str();
}

[[nodiscard]] std::vector<std::string> pressure_column_probe_phase_names(
    const std::vector<PressureColumnObservation>& observations) {
  std::vector<std::string> phases;
  for (const auto& observation : observations) {
    if (!contains(phases, observation.phase)) {
      phases.push_back(observation.phase);
    }
  }
  return phases;
}

[[nodiscard]] std::string format_probe_double(const double value) {
  std::ostringstream formatted;
  formatted << std::setprecision(9) << value;
  return formatted.str();
}

[[nodiscard]] std::string pressure_column_probe_values(
    const std::vector<PressureColumnObservation>& observations) {
  std::ostringstream joined;
  for (std::size_t index = 0; index < observations.size(); ++index) {
    if (index != 0) {
      joined << "|";
    }
    const auto& value = observations[index];
    joined << "phase=" << value.phase << ";i=" << value.i << ";j=" << value.j
           << ";k=" << value.k << ";P=" << format_probe_double(value.p)
           << ";PB=" << format_probe_double(value.pb)
           << ";P_PLUS_PB=" << format_probe_double(value.p_plus_pb)
           << ";MU=" << format_probe_double(value.mu)
           << ";MUB=" << format_probe_double(value.mub)
           << ";MU_PLUS_MUB=" << format_probe_double(value.mu_plus_mub)
           << ";PH=" << format_probe_double(value.ph)
           << ";PHB=" << format_probe_double(value.phb)
           << ";PH_PLUS_PHB=" << format_probe_double(value.ph_plus_phb)
           << ";T=" << format_probe_double(value.t)
           << ";QVAPOR=" << format_probe_double(value.qvapor)
           << ";HGT=" << format_probe_double(value.hgt);
  }
  return joined.str();
}

[[nodiscard]] const char* pressure_formula_observation_status_name(
    const tywrf::dynamics::PressureRefreshObservationStatus status) noexcept {
  switch (status) {
    case tywrf::dynamics::PressureRefreshObservationStatus::not_recorded:
      return "not_recorded";
    case tywrf::dynamics::PressureRefreshObservationStatus::recorded:
      return "recorded";
    case tywrf::dynamics::PressureRefreshObservationStatus::request_out_of_bounds:
      return "request_out_of_bounds";
    case tywrf::dynamics::PressureRefreshObservationStatus::
        request_outside_target_region:
      return "request_outside_target_region";
    case tywrf::dynamics::PressureRefreshObservationStatus::invalid_mu_total:
      return "invalid_mu_total";
    case tywrf::dynamics::PressureRefreshObservationStatus::
        invalid_pressure_levels:
      return "invalid_pressure_levels";
    case tywrf::dynamics::PressureRefreshObservationStatus::invalid_log_ratio:
      return "invalid_log_ratio";
    case tywrf::dynamics::PressureRefreshObservationStatus::invalid_delta_phi:
      return "invalid_delta_phi";
    case tywrf::dynamics::PressureRefreshObservationStatus::invalid_alb:
      return "invalid_alb";
    case tywrf::dynamics::PressureRefreshObservationStatus::invalid_pb:
      return "invalid_pb";
    case tywrf::dynamics::PressureRefreshObservationStatus::invalid_theta:
      return "invalid_theta";
    case tywrf::dynamics::PressureRefreshObservationStatus::invalid_alpha:
      return "invalid_alpha";
    case tywrf::dynamics::PressureRefreshObservationStatus::
        invalid_pressure_base:
      return "invalid_pressure_base";
    case tywrf::dynamics::PressureRefreshObservationStatus::
        invalid_total_pressure:
      return "invalid_total_pressure";
  }
  return "not_recorded";
}

[[nodiscard]] std::vector<tywrf::dynamics::PressureRefreshObservationRequest>
make_pressure_formula_observation_requests(const Options& options) {
  std::vector<tywrf::dynamics::PressureRefreshObservationRequest> requests;
  if (!pressure_formula_observation_enabled(options)) {
    return requests;
  }
  requests.reserve(
      options.pressure_column_probe_columns.size() *
      options.pressure_column_probe_levels.size());
  for (const auto [column_i, column_j] : options.pressure_column_probe_columns) {
    for (const auto level : options.pressure_column_probe_levels) {
      requests.push_back({column_i, column_j, level});
    }
  }
  return requests;
}

[[nodiscard]] std::string pressure_formula_observation_values(
    const std::vector<tywrf::dynamics::PressureRefreshFormulaObservation>& observations) {
  std::ostringstream joined;
  for (std::size_t index = 0; index < observations.size(); ++index) {
    if (index != 0) {
      joined << "|";
    }
    const auto& value = observations[index];
    joined << "status=" << pressure_formula_observation_status_name(value.status)
           << ";valid=" << static_cast<int>(value.valid) << ";i=" << value.i
           << ";j=" << value.j << ";k=" << value.k
           << ";mu_total=" << format_probe_double(value.mu_total)
           << ";pfu=" << format_probe_double(value.pfu)
           << ";pfd=" << format_probe_double(value.pfd)
           << ";phm=" << format_probe_double(value.phm)
           << ";log_ratio=" << format_probe_double(value.log_ratio)
           << ";phi_lower=" << format_probe_double(value.phi_lower)
           << ";phi_upper=" << format_probe_double(value.phi_upper)
           << ";delta_phi=" << format_probe_double(value.delta_phi)
           << ";ALB=" << format_probe_double(value.alb)
           << ";PB=" << format_probe_double(value.pb)
           << ";theta=" << format_probe_double(value.theta)
           << ";alpha_total=" << format_probe_double(value.alpha_total)
           << ";alpha_perturbation="
           << format_probe_double(value.alpha_perturbation)
           << ";alpha_from_wrf_branch="
           << format_probe_double(value.alpha_from_wrf_branch)
           << ";pressure_base=" << format_probe_double(value.pressure_base)
           << ";total_pressure=" << format_probe_double(value.total_pressure)
           << ";perturbation_pressure_pa="
           << format_probe_double(value.perturbation_pressure_pa);
  }
  return joined.str();
}

void capture_pressure_column_observations(
    const Options& options,
    const std::string_view phase,
    const tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    CandidateReport& report) {
  if (!pressure_column_probe_enabled(options)) {
    return;
  }

  const auto p_layout = candidate.p.layout();
  const auto ph_layout = candidate.ph.layout();
  const auto surface_layout = candidate.mu.layout();
  const auto static_layout = output_static.hgt.layout();
  const auto p = candidate.p.view();
  const auto pb = candidate.pb.view();
  const auto mu = candidate.mu.view();
  const auto mub = candidate.mub.view();
  const auto ph = candidate.ph.view();
  const auto phb = candidate.phb.view();
  const auto t = candidate.t.view();
  const auto qvapor = candidate.qvapor.view();
  const auto hgt = output_static.hgt.view();

  for (const auto [column_i, column_j] : options.pressure_column_probe_columns) {
    const auto i = p_layout.i_begin() + column_i;
    const auto j = p_layout.j_begin() + column_j;
    const auto surface_i = surface_layout.i_begin() + column_i;
    const auto surface_j = surface_layout.j_begin() + column_j;
    const auto static_i = static_layout.i_begin() + column_i;
    const auto static_j = static_layout.j_begin() + column_j;
    for (const auto level : options.pressure_column_probe_levels) {
      const auto k = p_layout.k_begin() + level;
      const auto ph_k = ph_layout.k_begin() + level;
      const double p_value = p(i, j, k);
      const double pb_value = pb(i, j, k);
      const double mu_value = mu(surface_i, surface_j);
      const double mub_value = mub(surface_i, surface_j);
      const double ph_value = ph(i, j, ph_k);
      const double phb_value = phb(i, j, ph_k);
      report.pressure_column_observations.push_back({
          std::string(phase),
          column_i,
          column_j,
          level,
          p_value,
          pb_value,
          p_value + pb_value,
          mu_value,
          mub_value,
          mu_value + mub_value,
          ph_value,
          phb_value,
          ph_value + phb_value,
          static_cast<double>(t(i, j, k)),
          static_cast<double>(qvapor(i, j, k)),
          static_cast<double>(hgt(static_i, static_j))});
    }
  }
}

void append_pressure_column_probe_timeline(CandidateReport& report) {
  if (report.pressure_column_observations.empty()) {
    return;
  }
  append_timeline_event(
      report,
      "pressure_column_probe",
      {
          timeline_field("enabled", "true"),
          timeline_field(
              "phase_count",
              static_cast<std::uint64_t>(
                  pressure_column_probe_phase_names(report.pressure_column_observations).size())),
          timeline_field(
              "record_count",
              static_cast<std::uint64_t>(report.pressure_column_observations.size())),
          timeline_field("evidence_only", "true"),
      });
}

struct PressureRefreshReadiness {
  bool static_refresh_applied = false;
  bool static_refresh_uses_reference_end = false;
  bool thermodynamic_base_state_consistency_ready = false;
  bool provider_terrain_uses_moved_candidate_hgt = false;
  bool provider_base_state_reconstruct_ok = false;
  bool base_state_sync_contract_ok = false;
  bool base_state_sync_dry_run = false;
  bool base_state_sync_applied = false;
  std::uint64_t would_sync_pb_point_count = 0;
  std::uint64_t would_sync_mub_point_count = 0;
  std::uint64_t would_sync_phb_point_count = 0;
  std::uint64_t sync_overlap_write_count = 0;
  std::uint64_t sync_halo_write_count = 0;
  bool pressure_refresh_compute_called = false;
  bool pressure_compute_dry_run = false;
  bool pressure_compute_dry_run_called = false;
  bool pressure_compute_dry_run_ok = false;
  std::uint64_t would_refresh_p_point_count = 0;
  std::uint64_t dry_run_invalid_p_point_count = 0;
  std::uint64_t dry_run_skipped_p_point_count = 0;
  std::uint64_t pressure_compute_dry_run_report_target_column_count = 0;
  std::uint64_t pressure_compute_dry_run_report_refreshed_point_count = 0;
  std::uint64_t pressure_compute_dry_run_report_invalid_point_count = 0;
  std::uint64_t pressure_compute_dry_run_report_skipped_point_count = 0;
  bool pressure_compute_dry_run_report_touched_overlap_cells = false;
  bool pressure_compute_dry_run_report_touched_halo_cells = false;
  bool pressure_refresh_applied = false;
  std::string provider_terrain_source_name;
  std::string provider_terrain_provenance;

  [[nodiscard]] bool ready() const noexcept {
    return static_refresh_applied && !static_refresh_uses_reference_end &&
           thermodynamic_base_state_consistency_ready &&
           provider_terrain_uses_moved_candidate_hgt &&
           provider_base_state_reconstruct_ok && base_state_sync_contract_ok &&
           base_state_sync_dry_run && !base_state_sync_applied &&
           would_sync_pb_point_count > 0 && would_sync_mub_point_count > 0 &&
           would_sync_phb_point_count > 0 && sync_overlap_write_count == 0 &&
           sync_halo_write_count == 0 && !pressure_refresh_compute_called &&
           pressure_compute_dry_run && pressure_compute_dry_run_called &&
           pressure_compute_dry_run_ok && would_refresh_p_point_count > 0 &&
           dry_run_invalid_p_point_count == 0 &&
           dry_run_skipped_p_point_count == 0 &&
           pressure_compute_dry_run_report_target_column_count > 0 &&
           pressure_compute_dry_run_report_refreshed_point_count ==
               would_refresh_p_point_count &&
           pressure_compute_dry_run_report_invalid_point_count == 0 &&
           pressure_compute_dry_run_report_skipped_point_count == 0 &&
           !pressure_compute_dry_run_report_touched_overlap_cells &&
           !pressure_compute_dry_run_report_touched_halo_cells &&
           !pressure_refresh_applied;
  }
};

enum class CandidateDispositionKind {
  selected_field,
  selected_field_pressure_refresh_experimental_apply,
  diagnostic_adapter,
  diagnostic_helper,
  diagnostic_dry_run,
  diagnostic_staging,
  diagnostic_probe,
};

struct CandidateDispositionInput {
  CandidateDispositionKind kind = CandidateDispositionKind::selected_field;
  bool pressure_refresh_applied = false;
  bool pressure_refresh_ready = false;
};

struct CandidateDisposition {
  const char* status = "selected_field_candidate_generated";
  const char* candidate_kind = "selected_field_integrator_v0";
  bool diagnostic_only = false;
  bool gate_candidate = true;
  bool integrator_output = true;
  bool experimental_pressure_refresh_apply = false;
};

[[nodiscard]] constexpr bool is_diagnostic_disposition_kind(
    const CandidateDispositionKind kind) noexcept {
  switch (kind) {
    case CandidateDispositionKind::selected_field:
    case CandidateDispositionKind::selected_field_pressure_refresh_experimental_apply:
      return false;
    case CandidateDispositionKind::diagnostic_adapter:
    case CandidateDispositionKind::diagnostic_helper:
    case CandidateDispositionKind::diagnostic_dry_run:
    case CandidateDispositionKind::diagnostic_staging:
    case CandidateDispositionKind::diagnostic_probe:
      return true;
  }
  return true;
}

[[nodiscard]] constexpr const char* disposition_status_name(
    const CandidateDispositionKind kind) noexcept {
  switch (kind) {
    case CandidateDispositionKind::selected_field:
      return "selected_field_candidate_generated";
    case CandidateDispositionKind::selected_field_pressure_refresh_experimental_apply:
      return "selected_field_pressure_refresh_experimental_apply_generated";
    case CandidateDispositionKind::diagnostic_adapter:
      return "selected_field_diagnostic_adapter_reported";
    case CandidateDispositionKind::diagnostic_helper:
      return "selected_field_diagnostic_helper_reported";
    case CandidateDispositionKind::diagnostic_dry_run:
      return "selected_field_diagnostic_dry_run_reported";
    case CandidateDispositionKind::diagnostic_staging:
      return "selected_field_diagnostic_staging_reported";
    case CandidateDispositionKind::diagnostic_probe:
      return "selected_field_diagnostic_probe_reported";
  }
  return "selected_field_diagnostic_unknown_reported";
}

[[nodiscard]] constexpr const char* disposition_candidate_kind_name(
    const CandidateDispositionKind kind) noexcept {
  switch (kind) {
    case CandidateDispositionKind::selected_field:
      return "selected_field_integrator_v0";
    case CandidateDispositionKind::selected_field_pressure_refresh_experimental_apply:
      return "selected_field_pressure_refresh_experimental_apply_v0";
    case CandidateDispositionKind::diagnostic_adapter:
      return "selected_field_diagnostic_adapter_v0";
    case CandidateDispositionKind::diagnostic_helper:
      return "selected_field_diagnostic_helper_v0";
    case CandidateDispositionKind::diagnostic_dry_run:
      return "selected_field_diagnostic_dry_run_v0";
    case CandidateDispositionKind::diagnostic_staging:
      return "selected_field_diagnostic_staging_v0";
    case CandidateDispositionKind::diagnostic_probe:
      return "selected_field_diagnostic_probe_v0";
  }
  return "selected_field_diagnostic_unknown_v0";
}

[[nodiscard]] constexpr CandidateDisposition candidate_disposition(
    const CandidateDispositionInput input) noexcept {
  const auto kind = input.kind;
  if (is_diagnostic_disposition_kind(kind)) {
    return {
        disposition_status_name(kind),
        disposition_candidate_kind_name(kind),
        true,
        false,
        false,
        false};
  }
  if (kind == CandidateDispositionKind::selected_field_pressure_refresh_experimental_apply) {
    return {
        disposition_status_name(kind),
        disposition_candidate_kind_name(kind),
        true,
        false,
        false,
        true};
  }

  const bool gate_candidate =
      input.pressure_refresh_applied ? input.pressure_refresh_ready : true;
  return {
      disposition_status_name(kind),
      disposition_candidate_kind_name(kind),
      false,
      gate_candidate,
      gate_candidate,
      false};
}

[[nodiscard]] CandidateDisposition selected_field_candidate_disposition(
    const Options& options,
    const CandidateReport& report,
    const std::optional<PressureRefreshReadiness>& pressure_refresh_readiness) {
  if (options.diagnostic_adapter_report) {
    return candidate_disposition({
        CandidateDispositionKind::diagnostic_adapter,
        report.diagnostic_adapter.has_value(),
        report.diagnostic_adapter.has_value() && report.diagnostic_adapter->ok()});
  }
  const auto kind = options.experimental_pressure_refresh_apply
                        ? CandidateDispositionKind::
                              selected_field_pressure_refresh_experimental_apply
                        : CandidateDispositionKind::selected_field;
  return candidate_disposition({
      kind,
      report.pressure_refresh.has_value(),
      pressure_refresh_readiness.has_value() && pressure_refresh_readiness->ready()});
}

void require_candidate_disposition_case(
    const std::string_view name,
    const CandidateDispositionInput input,
    const bool expected_diagnostic_only,
    const bool expected_gate_candidate,
    const bool expected_integrator_output,
    const char* expected_candidate_kind) {
  const auto disposition = candidate_disposition(input);
  if (disposition.diagnostic_only != expected_diagnostic_only ||
      disposition.gate_candidate != expected_gate_candidate ||
      disposition.integrator_output != expected_integrator_output ||
      std::string_view(disposition.candidate_kind) != expected_candidate_kind) {
    std::ostringstream message;
    message << "candidate disposition self-test failed for " << name;
    throw std::runtime_error(message.str());
  }
  if (is_diagnostic_disposition_kind(input.kind) &&
      std::string_view(disposition.candidate_kind) == "selected_field_integrator_v0") {
    std::ostringstream message;
    message << "diagnostic disposition used integrator candidate kind for " << name;
    throw std::runtime_error(message.str());
  }
}

int run_candidate_disposition_self_test() {
  require_candidate_disposition_case(
      "selected_field_no_pressure_refresh",
      {CandidateDispositionKind::selected_field, false, false},
      false,
      true,
      true,
      "selected_field_integrator_v0");
  require_candidate_disposition_case(
      "selected_field_pressure_refresh_ready",
      {CandidateDispositionKind::selected_field, true, true},
      false,
      true,
      true,
      "selected_field_integrator_v0");
  require_candidate_disposition_case(
      "selected_field_pressure_refresh_not_ready",
      {CandidateDispositionKind::selected_field, true, false},
      false,
      false,
      false,
      "selected_field_integrator_v0");
  require_candidate_disposition_case(
      "selected_field_experimental_pressure_refresh_apply",
      {CandidateDispositionKind::selected_field_pressure_refresh_experimental_apply, true, true},
      true,
      false,
      false,
      "selected_field_pressure_refresh_experimental_apply_v0");
  require_candidate_disposition_case(
      "selected_field_experimental_pressure_refresh_apply_failed",
      {CandidateDispositionKind::selected_field_pressure_refresh_experimental_apply, false, false},
      true,
      false,
      false,
      "selected_field_pressure_refresh_experimental_apply_v0");

  const std::array<std::pair<CandidateDispositionKind, const char*>, 5> diagnostic_cases = {{
      {CandidateDispositionKind::diagnostic_adapter, "selected_field_diagnostic_adapter_v0"},
      {CandidateDispositionKind::diagnostic_helper, "selected_field_diagnostic_helper_v0"},
      {CandidateDispositionKind::diagnostic_dry_run, "selected_field_diagnostic_dry_run_v0"},
      {CandidateDispositionKind::diagnostic_staging, "selected_field_diagnostic_staging_v0"},
      {CandidateDispositionKind::diagnostic_probe, "selected_field_diagnostic_probe_v0"},
  }};
  for (const auto& [kind, expected_candidate_kind] : diagnostic_cases) {
    require_candidate_disposition_case(
        expected_candidate_kind,
        {kind, true, true},
        true,
        false,
        false,
        expected_candidate_kind);
  }
  return 0;
}

[[nodiscard]] PressureRefreshReadiness evaluate_pressure_refresh_readiness(
    const CandidateReport& report) {
  PressureRefreshReadiness readiness;
  readiness.static_refresh_applied =
      report.static_refresh.ok() && report.changed_static_template_points > 0;
  readiness.static_refresh_uses_reference_end = false;
  if (report.pressure_refresh_provider_probe.has_value()) {
    const auto& provider = *report.pressure_refresh_provider_probe;
    readiness.provider_base_state_reconstruct_ok =
        provider.ok() && provider.allocated_buffers && provider.wrote_pb &&
        provider.wrote_t_init && provider.wrote_mub && provider.wrote_alb &&
        provider.wrote_phb;
    readiness.provider_terrain_uses_moved_candidate_hgt =
        provider.ok() && provider.terrain_override_used &&
        provider.terrain_source_name == "moved_candidate_HGT" &&
        provider.terrain_provenance == "override:moved_candidate_HGT";
    readiness.provider_terrain_source_name = provider.terrain_source_name;
    readiness.provider_terrain_provenance = provider.terrain_provenance;
  }
  if (report.pressure_refresh_dry_run_contract.has_value()) {
    const auto& dry_run = *report.pressure_refresh_dry_run_contract;
    readiness.base_state_sync_contract_ok = dry_run.provider_ok &&
                                            dry_run.base_state_sync_contract_ok;
    readiness.base_state_sync_dry_run = dry_run.base_state_sync_dry_run;
    readiness.base_state_sync_applied = dry_run.base_state_sync_applied;
    readiness.would_sync_pb_point_count = dry_run.would_sync_pb_point_count;
    readiness.would_sync_mub_point_count = dry_run.would_sync_mub_point_count;
    readiness.would_sync_phb_point_count = dry_run.would_sync_phb_point_count;
    readiness.sync_overlap_write_count = dry_run.sync_overlap_write_count;
    readiness.sync_halo_write_count = dry_run.sync_halo_write_count;
    readiness.pressure_refresh_compute_called = dry_run.calls_pressure_refresh_compute;
    readiness.pressure_compute_dry_run = dry_run.pressure_compute_dry_run;
    readiness.pressure_compute_dry_run_called =
        dry_run.pressure_compute_dry_run_called;
    readiness.pressure_compute_dry_run_ok = dry_run.pressure_compute_dry_run_ok;
    readiness.would_refresh_p_point_count = dry_run.would_refresh_p_point_count;
    readiness.dry_run_invalid_p_point_count =
        dry_run.dry_run_invalid_p_point_count;
    readiness.dry_run_skipped_p_point_count =
        dry_run.pressure_compute_dry_run_report.skipped_point_count;
    readiness.pressure_compute_dry_run_report_target_column_count =
        dry_run.pressure_compute_dry_run_report.target_column_count;
    readiness.pressure_compute_dry_run_report_refreshed_point_count =
        dry_run.pressure_compute_dry_run_report.refreshed_point_count;
    readiness.pressure_compute_dry_run_report_invalid_point_count =
        dry_run.pressure_compute_dry_run_report.invalid_point_count;
    readiness.pressure_compute_dry_run_report_skipped_point_count =
        dry_run.pressure_compute_dry_run_report.skipped_point_count;
    readiness.pressure_compute_dry_run_report_touched_overlap_cells =
        dry_run.pressure_compute_dry_run_report.touched_overlap_cells;
    readiness.pressure_compute_dry_run_report_touched_halo_cells =
        dry_run.pressure_compute_dry_run_report.touched_halo_cells;
    readiness.pressure_refresh_applied = dry_run.pressure_refresh_applied;
  }
  readiness.thermodynamic_base_state_consistency_ready =
      readiness.provider_terrain_uses_moved_candidate_hgt &&
      readiness.provider_base_state_reconstruct_ok &&
      readiness.base_state_sync_contract_ok &&
      readiness.base_state_sync_dry_run && !readiness.base_state_sync_applied &&
      readiness.would_sync_pb_point_count > 0 &&
      readiness.would_sync_mub_point_count > 0 &&
      readiness.would_sync_phb_point_count > 0 &&
      readiness.sync_overlap_write_count == 0 &&
      readiness.sync_halo_write_count == 0 &&
      !readiness.pressure_refresh_compute_called &&
      readiness.pressure_compute_dry_run &&
      readiness.pressure_compute_dry_run_called &&
      readiness.pressure_compute_dry_run_ok &&
      readiness.would_refresh_p_point_count > 0 &&
      readiness.dry_run_invalid_p_point_count == 0 &&
      readiness.dry_run_skipped_p_point_count == 0 &&
      readiness.pressure_compute_dry_run_report_target_column_count > 0 &&
      readiness.pressure_compute_dry_run_report_refreshed_point_count ==
          readiness.would_refresh_p_point_count &&
      readiness.pressure_compute_dry_run_report_invalid_point_count == 0 &&
      readiness.pressure_compute_dry_run_report_skipped_point_count == 0 &&
      !readiness.pressure_compute_dry_run_report_touched_overlap_cells &&
      !readiness.pressure_compute_dry_run_report_touched_halo_cells &&
      !readiness.pressure_refresh_applied;
  return readiness;
}

[[nodiscard]] std::string pressure_refresh_not_ready_message(
    const PressureRefreshReadiness readiness) {
  std::ostringstream message;
  message << "pressure_refresh_not_ready: static refreshed but thermodynamic/base-state "
             "consistency missing";
  message << "; static_refresh_applied="
          << (readiness.static_refresh_applied ? "true" : "false");
  message << "; static_refresh_uses_reference_end="
          << (readiness.static_refresh_uses_reference_end ? "true" : "false");
  message << "; thermodynamic_base_state_consistency_ready="
          << (readiness.thermodynamic_base_state_consistency_ready ? "true" : "false");
  message << "; exposed T/PH are parent interpolated from d01 start-state fields, but "
             "PB/PHB/MUB/P base-state ownership is still preserved d02 start-state data";
  message << "; provider_terrain_uses_moved_candidate_hgt="
          << (readiness.provider_terrain_uses_moved_candidate_hgt ? "true" : "false");
  message << "; provider_base_state_reconstruct_ok="
          << (readiness.provider_base_state_reconstruct_ok ? "true" : "false");
  if (!readiness.provider_terrain_source_name.empty()) {
    message << "; provider_terrain_source=" << readiness.provider_terrain_source_name;
  }
  if (!readiness.provider_terrain_provenance.empty()) {
    message << "; provider_terrain_provenance=" << readiness.provider_terrain_provenance;
  }
  message << "; base_state_sync_contract_ok="
          << (readiness.base_state_sync_contract_ok ? "true" : "false");
  message << "; base_state_sync_dry_run="
          << (readiness.base_state_sync_dry_run ? "true" : "false");
  message << "; base_state_sync_applied="
          << (readiness.base_state_sync_applied ? "true" : "false");
  message << "; would_sync_pb_point_count=" << readiness.would_sync_pb_point_count;
  message << "; would_sync_mub_point_count=" << readiness.would_sync_mub_point_count;
  message << "; would_sync_phb_point_count=" << readiness.would_sync_phb_point_count;
  message << "; sync_overlap_write_count=" << readiness.sync_overlap_write_count;
  message << "; sync_halo_write_count=" << readiness.sync_halo_write_count;
  message << "; pressure_refresh_compute_called="
          << (readiness.pressure_refresh_compute_called ? "true" : "false");
  message << "; pressure_compute_dry_run="
          << (readiness.pressure_compute_dry_run ? "true" : "false");
  message << "; pressure_compute_dry_run_called="
          << (readiness.pressure_compute_dry_run_called ? "true" : "false");
  message << "; pressure_compute_dry_run_ok="
          << (readiness.pressure_compute_dry_run_ok ? "true" : "false");
  message << "; would_refresh_p_point_count="
          << readiness.would_refresh_p_point_count;
  message << "; dry_run_invalid_p_point_count="
          << readiness.dry_run_invalid_p_point_count;
  message << "; dry_run_skipped_p_point_count="
          << readiness.dry_run_skipped_p_point_count;
  message << "; pressure_compute_dry_run_report_target_column_count="
          << readiness.pressure_compute_dry_run_report_target_column_count;
  message << "; pressure_compute_dry_run_report_refreshed_point_count="
          << readiness.pressure_compute_dry_run_report_refreshed_point_count;
  message << "; pressure_compute_dry_run_report_invalid_point_count="
          << readiness.pressure_compute_dry_run_report_invalid_point_count;
  message << "; pressure_compute_dry_run_report_skipped_point_count="
          << readiness.pressure_compute_dry_run_report_skipped_point_count;
  message << "; pressure_compute_dry_run_report_touched_overlap_cells="
          << (readiness.pressure_compute_dry_run_report_touched_overlap_cells ? "true"
                                                                              : "false");
  message << "; pressure_compute_dry_run_report_touched_halo_cells="
          << (readiness.pressure_compute_dry_run_report_touched_halo_cells ? "true" : "false");
  message << "; pressure_refresh_applied="
          << (readiness.pressure_refresh_applied ? "true" : "false");
  return message.str();
}

void require_pressure_refresh_ready_for_compute(const CandidateReport& report) {
  const auto readiness = evaluate_pressure_refresh_readiness(report);
  if (readiness.ready()) {
    return;
  }
  throw std::runtime_error(pressure_refresh_not_ready_message(readiness));
}

void require_pressure_refresh_inputs_ready(
    const tywrf::io::KrosaPressureRefreshReadResult& inputs) {
  if (inputs.ok() || inputs.base_state_reconstruction_inputs_ready()) {
    return;
  }

  std::ostringstream message;
  message << "pressure refresh metadata source is not ready for provider-backed selected-field "
             "candidate";
  if (!inputs.report.missing_names.empty()) {
    message << "; missing direct inputs=" << join_variables(inputs.report.missing_names);
  }
  if (!inputs.report.missing_base_state_reconstruction_names.empty()) {
    message << "; missing provider inputs="
            << join_variables(inputs.report.missing_base_state_reconstruction_names);
  }
  throw std::runtime_error(message.str());
}

void probe_pressure_refresh_provider_readiness(
    const Options& options,
    const tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    CandidateReport& report) {
  tywrf::FieldStorage3D<float> direct_alb(candidate.grid.mass_layout());
  const auto inputs = tywrf::io::read_krosa_pressure_refresh_inputs(
      options.template_path,
      candidate.grid,
      direct_alb,
      {.time_index = options.template_time_index});
  require_pressure_refresh_inputs_ready(inputs);

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto terrain_override = make_moved_candidate_hgt_terrain_override(output_static);
  auto provider_report =
      provider.reconstruct(candidate.grid, inputs.metadata, terrain_override);
  const bool uses_moved_candidate_hgt =
      provider_report.terrain_override_used &&
      provider_report.terrain_source_name == "moved_candidate_HGT";
  if (!provider_report.ok() || !provider_report.allocated_buffers ||
      !provider_report.wrote_pb || !provider_report.wrote_t_init ||
      !provider_report.wrote_mub || !provider_report.wrote_alb ||
      !provider_report.wrote_phb || !uses_moved_candidate_hgt) {
    std::ostringstream message;
    message << "pressure_refresh_provider_probe_failed";
    if (provider_report.result.message != nullptr &&
        provider_report.result.message[0] != '\0') {
      message << ": " << provider_report.result.message;
    }
    message << "; provider_terrain_uses_moved_candidate_hgt="
            << (uses_moved_candidate_hgt ? "true" : "false");
    message << "; provider_terrain_source=" << provider_report.terrain_source_name;
    message << "; provider_terrain_provenance=" << provider_report.terrain_provenance;
    throw std::runtime_error(message.str());
  }

  report.pressure_refresh_provider_probe = provider_report;
  report.pressure_refresh_metadata_source = options.template_path;
  report.pressure_refresh_metadata_time_index = options.template_time_index;
}

void probe_pressure_refresh_dry_run_contract(
    const Options& options,
    tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    CandidateReport& report) {
  tywrf::FieldStorage3D<float> direct_alb(candidate.grid.mass_layout());
  const auto inputs = tywrf::io::read_krosa_pressure_refresh_inputs(
      options.template_path,
      candidate.grid,
      direct_alb,
      {.time_index = options.template_time_index});
  require_pressure_refresh_inputs_ready(inputs);

  const auto terrain_override = make_moved_candidate_hgt_terrain_override(output_static);
  tywrf::dynamics::KrosaPressureRefreshHookOptions hook_options{};
  hook_options.terrain_override = &terrain_override;
  hook_options.base_state_sync_dry_run = true;
  hook_options.pressure_compute_dry_run = true;
  auto dry_run_report = tywrf::dynamics::apply_krosa_moving_nest_pressure_refresh_hook(
      report.remap_plan,
      candidate,
      inputs.metadata,
      hook_options);
  report.pressure_refresh_dry_run_contract = dry_run_report;
  report.pressure_refresh_metadata_source = options.template_path;
  report.pressure_refresh_metadata_time_index = options.template_time_index;

  const bool dry_run_contract_ok =
      dry_run_report.ok() && dry_run_report.provider_ok &&
      dry_run_report.base_state_sync_dry_run &&
      dry_run_report.base_state_sync_contract_ok &&
      !dry_run_report.base_state_sync_applied &&
      dry_run_report.would_sync_pb_point_count > 0 &&
      dry_run_report.would_sync_mub_point_count > 0 &&
      dry_run_report.would_sync_phb_point_count > 0 &&
      dry_run_report.sync_overlap_write_count == 0 &&
      dry_run_report.sync_halo_write_count == 0 &&
      !dry_run_report.calls_pressure_refresh_compute &&
      dry_run_report.pressure_compute_dry_run &&
      dry_run_report.pressure_compute_dry_run_called &&
      dry_run_report.pressure_compute_dry_run_ok &&
      dry_run_report.would_refresh_p_point_count > 0 &&
      dry_run_report.dry_run_invalid_p_point_count == 0 &&
      dry_run_report.pressure_compute_dry_run_report.target_column_count > 0 &&
      dry_run_report.pressure_compute_dry_run_report.refreshed_point_count ==
          dry_run_report.would_refresh_p_point_count &&
      dry_run_report.pressure_compute_dry_run_report.invalid_point_count == 0 &&
      dry_run_report.pressure_compute_dry_run_report.skipped_point_count == 0 &&
      !dry_run_report.pressure_compute_dry_run_report.touched_overlap_cells &&
      !dry_run_report.pressure_compute_dry_run_report.touched_halo_cells &&
      !dry_run_report.pressure_refresh_applied;
  const auto readiness = evaluate_pressure_refresh_readiness(report);
  append_timeline_event(
      report,
      "pressure_refresh_readiness",
      {
          timeline_field("opt_in", "true"),
          timeline_field("ready", readiness.ready() ? "true" : "false"),
          timeline_field("provider_ok", dry_run_report.provider_ok ? "true" : "false"),
          timeline_field(
              "base_state_source_sync_readiness_check",
              dry_run_report.base_state_sync_dry_run ? "true" : "false"),
          timeline_field(
              "source_sync_planned_pb_points", dry_run_report.would_sync_pb_point_count),
          timeline_field(
              "source_sync_planned_mub_points", dry_run_report.would_sync_mub_point_count),
          timeline_field(
              "source_sync_planned_phb_points", dry_run_report.would_sync_phb_point_count),
          timeline_field(
              "production_refresh_planned_p_points", dry_run_report.would_refresh_p_point_count),
          timeline_field(
              "readiness_invalid_p_points", dry_run_report.dry_run_invalid_p_point_count),
          timeline_field(
              "readiness_skipped_p_points",
              dry_run_report.pressure_compute_dry_run_report.skipped_point_count),
      });
  if (dry_run_contract_ok) {
    return;
  }

  std::ostringstream message;
  message << "pressure_refresh_dry_run_contract_failed";
  if (dry_run_report.result.message != nullptr &&
      dry_run_report.result.message[0] != '\0') {
    message << ": " << dry_run_report.result.message;
  }
  message << "; provider_ok=" << (dry_run_report.provider_ok ? "true" : "false")
          << "; base_state_sync_contract_ok="
          << (dry_run_report.base_state_sync_contract_ok ? "true" : "false")
          << "; base_state_sync_dry_run="
          << (dry_run_report.base_state_sync_dry_run ? "true" : "false")
          << "; base_state_sync_applied="
          << (dry_run_report.base_state_sync_applied ? "true" : "false")
          << "; would_sync_pb_point_count=" << dry_run_report.would_sync_pb_point_count
          << "; would_sync_mub_point_count=" << dry_run_report.would_sync_mub_point_count
          << "; would_sync_phb_point_count=" << dry_run_report.would_sync_phb_point_count
          << "; sync_overlap_write_count=" << dry_run_report.sync_overlap_write_count
          << "; sync_halo_write_count=" << dry_run_report.sync_halo_write_count
          << "; pressure_refresh_compute_called="
          << (dry_run_report.calls_pressure_refresh_compute ? "true" : "false")
          << "; pressure_compute_dry_run="
          << (dry_run_report.pressure_compute_dry_run ? "true" : "false")
          << "; pressure_compute_dry_run_called="
          << (dry_run_report.pressure_compute_dry_run_called ? "true" : "false")
          << "; pressure_compute_dry_run_ok="
          << (dry_run_report.pressure_compute_dry_run_ok ? "true" : "false")
          << "; would_refresh_p_point_count=" << dry_run_report.would_refresh_p_point_count
          << "; dry_run_invalid_p_point_count="
          << dry_run_report.dry_run_invalid_p_point_count
          << "; dry_run_skipped_p_point_count="
          << dry_run_report.pressure_compute_dry_run_report.skipped_point_count
          << "; pressure_compute_dry_run_report_target_column_count="
          << dry_run_report.pressure_compute_dry_run_report.target_column_count
          << "; pressure_compute_dry_run_report_refreshed_point_count="
          << dry_run_report.pressure_compute_dry_run_report.refreshed_point_count
          << "; pressure_compute_dry_run_report_invalid_point_count="
          << dry_run_report.pressure_compute_dry_run_report.invalid_point_count
          << "; pressure_compute_dry_run_report_skipped_point_count="
          << dry_run_report.pressure_compute_dry_run_report.skipped_point_count
          << "; pressure_compute_dry_run_report_touched_overlap_cells="
          << (dry_run_report.pressure_compute_dry_run_report.touched_overlap_cells ? "true"
                                                                                   : "false")
          << "; pressure_compute_dry_run_report_touched_halo_cells="
          << (dry_run_report.pressure_compute_dry_run_report.touched_halo_cells ? "true"
                                                                                : "false")
          << "; pressure_refresh_applied="
          << (dry_run_report.pressure_refresh_applied ? "true" : "false");
  throw std::runtime_error(message.str());
}

void require_pressure_refresh_hook_success(
    const tywrf::dynamics::KrosaPressureRefreshHookReport& report) {
  const bool uses_moved_candidate_hgt =
      report.provider_report.terrain_override_used &&
      report.provider_report.terrain_source_name == "moved_candidate_HGT" &&
      report.provider_report.terrain_provenance == "override:moved_candidate_HGT";
  if (report.ok() && report.provider_ok && report.staging_ok &&
      report.calls_pressure_refresh_compute && report.pressure_refresh_applied &&
      report.base_state_sync_applied && report.synced_pb_point_count > 0 &&
      report.synced_mub_point_count > 0 && report.synced_phb_point_count > 0 &&
      report.compute_report.refreshed_point_count > 0 &&
      report.compute_report.invalid_point_count == 0 &&
      report.compute_report.skipped_point_count == 0 &&
      !report.touched_overlap_cells && !report.touched_halo_cells &&
      uses_moved_candidate_hgt) {
    return;
  }

  std::ostringstream message;
  message << "selected-field pressure refresh hook failed";
  if (report.result.message != nullptr && report.result.message[0] != '\0') {
    message << ": " << report.result.message;
  }
  message << " provider_ok=" << (report.provider_ok ? "true" : "false")
          << " staging_ok=" << (report.staging_ok ? "true" : "false")
          << " compute_called=" << (report.calls_pressure_refresh_compute ? "true" : "false")
          << " applied=" << (report.pressure_refresh_applied ? "true" : "false")
          << " base_state_sync_applied="
          << (report.base_state_sync_applied ? "true" : "false")
          << " synced_pb=" << report.synced_pb_point_count
          << " synced_mub=" << report.synced_mub_point_count
          << " synced_phb=" << report.synced_phb_point_count
          << " refreshed_points=" << report.compute_report.refreshed_point_count
          << " invalid_points=" << report.compute_report.invalid_point_count
          << " skipped_points=" << report.compute_report.skipped_point_count
          << " touched_overlap=" << (report.touched_overlap_cells ? "true" : "false")
          << " touched_halo=" << (report.touched_halo_cells ? "true" : "false")
          << " terrain_override_used="
          << (report.provider_report.terrain_override_used ? "true" : "false")
          << " terrain_source=" << report.provider_report.terrain_source_name
          << " terrain_provenance=" << report.provider_report.terrain_provenance;
  throw std::runtime_error(message.str());
}

void apply_pressure_refresh(
    const Options& options,
    tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    CandidateReport& report) {
  tywrf::FieldStorage3D<float> direct_alb(candidate.grid.mass_layout());
  const auto metadata = tywrf::io::read_krosa_pressure_refresh_inputs(
      options.template_path,
      candidate.grid,
      direct_alb,
      {.time_index = options.template_time_index});
  require_pressure_refresh_inputs_ready(metadata);

  const std::vector<float> p_before(
      candidate.p.data(), candidate.p.data() + candidate.p.size());
  const std::vector<float> pb_before(
      candidate.pb.data(), candidate.pb.data() + candidate.pb.size());
  const std::vector<float> mub_before(
      candidate.mub.data(), candidate.mub.data() + candidate.mub.size());
  const std::vector<float> phb_before(
      candidate.phb.data(), candidate.phb.data() + candidate.phb.size());
  const auto terrain_override = make_moved_candidate_hgt_terrain_override(output_static);
  auto formula_observation_requests =
      make_pressure_formula_observation_requests(options);
  std::vector<tywrf::dynamics::PressureRefreshFormulaObservation>
      formula_observation_records(formula_observation_requests.size());
  tywrf::dynamics::KrosaPressureRefreshHookOptions hook_options{};
  hook_options.terrain_override = &terrain_override;
  if (!formula_observation_requests.empty()) {
    hook_options.pressure_refresh.observation = {
        formula_observation_requests.data(),
        static_cast<std::int32_t>(formula_observation_requests.size()),
        formula_observation_records.data(),
        static_cast<std::int32_t>(formula_observation_records.size())};
  }
  auto hook_report = tywrf::dynamics::apply_krosa_moving_nest_pressure_refresh_hook(
      report.remap_plan,
      candidate,
      metadata.metadata,
      hook_options);
  require_pressure_refresh_hook_success(hook_report);
  require_finite_strict_fields(candidate);

  report.pressure_refresh_changed_p_points = changed_points(p_before, candidate.p);
  report.pressure_refresh_changed_pb_points = changed_points(pb_before, candidate.pb);
  report.pressure_refresh_changed_mub_points = changed_points(mub_before, candidate.mub);
  report.pressure_refresh_changed_phb_points = changed_points(phb_before, candidate.phb);
  if (report.pressure_refresh_changed_p_points == 0) {
    throw std::runtime_error("selected-field pressure refresh did not change any P point");
  }
  report.pressure_refresh_metadata_source = options.template_path;
  report.pressure_refresh_metadata_time_index = options.template_time_index;
  report.pressure_refresh = hook_report;
  if (!formula_observation_records.empty()) {
    report.pressure_formula_observations =
        std::move(formula_observation_records);
  }
  append_timeline_event(
      report,
      "pressure_refresh_apply",
      {
          timeline_field("opt_in", "true"),
          timeline_field("applied", "true"),
          timeline_field("refreshed_points", hook_report.compute_report.refreshed_point_count),
          timeline_field("synced_pb_points", hook_report.synced_pb_point_count),
          timeline_field("synced_mub_points", hook_report.synced_mub_point_count),
          timeline_field("synced_phb_points", hook_report.synced_phb_point_count),
          timeline_field("changed_p_points", report.pressure_refresh_changed_p_points),
          timeline_field("changed_pb_points", report.pressure_refresh_changed_pb_points),
          timeline_field("changed_mub_points", report.pressure_refresh_changed_mub_points),
          timeline_field("changed_phb_points", report.pressure_refresh_changed_phb_points),
      });
}

void fill_diagnostic_adapter_staging(
    const tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    const tywrf::FieldStorage3D<float>& alb,
    DiagnosticBaseStateAdapterStaging& staging) {
  copy_storage(candidate.phb, staging.phb, "diagnostic adapter PHB");
  copy_storage(candidate.mub, staging.mub, "diagnostic adapter MUB");
  copy_storage(candidate.pb, staging.pb, "diagnostic adapter PB");
  copy_storage(candidate.t, staging.t_init, "diagnostic adapter T_INIT");
  copy_storage(alb, staging.alb, "diagnostic adapter ALB");
  copy_storage(output_static.hgt, staging.ht, "diagnostic adapter HT/HGT");
}

[[nodiscard]] DiagnosticAdapterProviderSource
build_diagnostic_adapter_provider_source(
    const tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    const tywrf::io::KrosaPressureRefreshMetadata& metadata) {
  DiagnosticAdapterProviderSource source;
  const auto terrain_override =
      make_moved_candidate_hgt_terrain_override(output_static);
  source.report = source.provider.reconstruct(
      candidate.grid,
      metadata,
      terrain_override);

  const bool terrain_ok =
      source.report.terrain_override_used &&
      source.report.terrain_source_name == "moved_candidate_HGT" &&
      source.report.terrain_provenance == "override:moved_candidate_HGT";
  if (source.report.ok() && source.report.allocated_buffers &&
      source.report.wrote_pb && source.report.wrote_t_init &&
      source.report.wrote_mub && source.report.wrote_alb &&
      source.report.wrote_phb && terrain_ok) {
    return source;
  }

  std::ostringstream message;
  message << "diagnostic_adapter_provider_source_failed";
  if (source.report.result.message != nullptr &&
      source.report.result.message[0] != '\0') {
    message << ": " << source.report.result.message;
  }
  message << "; provider_ok=" << (source.report.ok() ? "true" : "false")
          << "; allocated_buffers="
          << (source.report.allocated_buffers ? "true" : "false")
          << "; wrote_pb=" << (source.report.wrote_pb ? "true" : "false")
          << "; wrote_t_init="
          << (source.report.wrote_t_init ? "true" : "false")
          << "; wrote_mub=" << (source.report.wrote_mub ? "true" : "false")
          << "; wrote_alb=" << (source.report.wrote_alb ? "true" : "false")
          << "; wrote_phb=" << (source.report.wrote_phb ? "true" : "false")
          << "; terrain_override_used="
          << (source.report.terrain_override_used ? "true" : "false")
          << "; terrain_source=" << source.report.terrain_source_name
          << "; terrain_provenance=" << source.report.terrain_provenance;
  throw std::runtime_error(message.str());
}

[[nodiscard]] tywrf::nest::ExposedBaseStateViews<const float>
diagnostic_adapter_provider_source_views(
    const DiagnosticAdapterProviderSource& source,
    const StaticFieldSet& output_static) {
  const auto provider_views = source.provider.views();
  return {
      provider_views.phb,
      provider_views.mub,
      provider_views.pb,
      provider_views.t_init,
      provider_views.alb,
      output_static.hgt.view()};
}

[[nodiscard]] tywrf::nest::NormalCandidateBaseStateSourceViews<const float>
normal_base_state_source_views(
    const tywrf::dynamics::KrosaBaseStateProvider& provider,
    const StaticFieldSet& output_static) {
  const auto provider_views = provider.views();
  return {
      provider_views.phb,
      provider_views.mub,
      provider_views.pb,
      output_static.hgt.view()};
}

[[nodiscard]] tywrf::nest::NormalCandidateBaseStateTargetViews<float>
normal_base_state_target_views(
    tywrf::State<float>& candidate,
    StaticFieldSet& output_static) {
  return {
      candidate.phb.view(),
      candidate.mub.view(),
      candidate.pb.view(),
      output_static.hgt.view()};
}

void run_normal_base_state_producer(
    const Options& options,
    tywrf::State<float>& candidate,
    StaticFieldSet& output_static,
    CandidateReport& report) {
  tywrf::FieldStorage3D<float> metadata_alb(candidate.grid.mass_layout());
  const auto metadata = tywrf::io::read_krosa_pressure_refresh_inputs(
      options.d02_start_state_path,
      candidate.grid,
      metadata_alb,
      {.time_index = options.d02_time_index});
  require_pressure_refresh_inputs_ready(metadata);

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto terrain_override =
      make_moved_candidate_hgt_terrain_override(output_static);
  auto provider_report = provider.reconstruct(
      candidate.grid,
      metadata.metadata,
      terrain_override);

  const bool terrain_ok =
      provider_report.terrain_override_used &&
      provider_report.terrain_source_name == "moved_candidate_HGT" &&
      provider_report.terrain_provenance == "override:moved_candidate_HGT";
  if (!provider_report.ok() || !provider_report.allocated_buffers ||
      !provider_report.wrote_pb || !provider_report.wrote_mub ||
      !provider_report.wrote_phb || !terrain_ok) {
    std::ostringstream message;
    message << "normal_base_state_provider_source_failed";
    if (provider_report.result.message != nullptr &&
        provider_report.result.message[0] != '\0') {
      message << ": " << provider_report.result.message;
    }
    message << "; provider_ok=" << (provider_report.ok() ? "true" : "false")
            << "; allocated_buffers="
            << (provider_report.allocated_buffers ? "true" : "false")
            << "; wrote_pb=" << (provider_report.wrote_pb ? "true" : "false")
            << "; wrote_mub=" << (provider_report.wrote_mub ? "true" : "false")
            << "; wrote_phb=" << (provider_report.wrote_phb ? "true" : "false")
            << "; terrain_override_used="
            << (provider_report.terrain_override_used ? "true" : "false")
            << "; terrain_source=" << provider_report.terrain_source_name
            << "; terrain_provenance=" << provider_report.terrain_provenance;
    throw std::runtime_error(message.str());
  }

  const std::vector<float> p_before(
      candidate.p.data(), candidate.p.data() + candidate.p.size());
  const std::vector<float> pb_before(
      candidate.pb.data(), candidate.pb.data() + candidate.pb.size());
  const std::vector<float> mub_before(
      candidate.mub.data(), candidate.mub.data() + candidate.mub.size());
  const std::vector<float> phb_before(
      candidate.phb.data(), candidate.phb.data() + candidate.phb.size());
  const std::vector<float> hgt_before(
      output_static.hgt.data(), output_static.hgt.data() + output_static.hgt.size());

  const auto producer_report =
      tywrf::nest::apply_normal_candidate_base_state_exchange(
          report.remap_plan,
          normal_base_state_source_views(provider, output_static),
          normal_base_state_target_views(candidate, output_static));
  if (!producer_report.ok()) {
    std::ostringstream message;
    message << "normal_base_state_producer_failed";
    if (producer_report.result.message != nullptr &&
        producer_report.result.message[0] != '\0') {
      message << ": " << producer_report.result.message;
    }
    throw std::runtime_error(message.str());
  }

  report.normal_base_state_provider_source = provider_report;
  report.normal_base_state_producer = producer_report;
  report.normal_base_state_metadata_source = options.d02_start_state_path;
  report.normal_base_state_metadata_time_index = options.d02_time_index;
  report.normal_base_state_changed_p_points = changed_points(p_before, candidate.p);
  report.normal_base_state_changed_pb_points = changed_points(pb_before, candidate.pb);
  report.normal_base_state_changed_mub_points = changed_points(mub_before, candidate.mub);
  report.normal_base_state_changed_phb_points = changed_points(phb_before, candidate.phb);
  report.normal_base_state_changed_hgt_points = changed_points(hgt_before, output_static.hgt);

  if (report.normal_base_state_changed_p_points != 0) {
    throw std::runtime_error("normal base-state producer unexpectedly changed P");
  }
  require_finite_storage("normal base-state PB", candidate.pb);
  require_finite_storage("normal base-state PHB", candidate.phb);
  require_finite_storage("normal base-state MUB", candidate.mub);
  require_finite_static_fields(output_static);

  append_timeline_event(
      report,
      "normal_base_state_producer",
      {
          timeline_field("source", std::string(producer_report.source)),
          timeline_field(
              "source_origin",
              std::string(producer_report.source_origin)),
          timeline_field("ok", std::string(bool_text(producer_report.ok()))),
          timeline_field(
              "diagnostic_only",
              std::string(bool_text(producer_report.diagnostic_only))),
          timeline_field(
              "normal_candidate_producer",
              std::string(bool_text(producer_report.normal_candidate_producer))),
          timeline_field(
              "writes_candidate",
              std::string(bool_text(producer_report.writes_candidate))),
          timeline_field(
              "uses_reference_end_truth",
              std::string(bool_text(producer_report.uses_reference_end_truth))),
          timeline_field(
              "uses_direct_p_shortcut",
              std::string(bool_text(producer_report.uses_direct_p_shortcut))),
          timeline_field(
              "reads_direct_p",
              std::string(bool_text(producer_report.reads_direct_p))),
          timeline_field(
              "writes_p",
              std::string(bool_text(producer_report.writes_p))),
          timeline_field(
              "no_gate_pass_claim",
              std::string(bool_text(producer_report.no_gate_pass_claim))),
          timeline_field(
              "no_00_20_progression",
              std::string(bool_text(producer_report.no_00_20_progression))),
          timeline_field(
              "terrain_source",
              provider_report.terrain_source_name),
          timeline_field(
              "terrain_provenance",
              provider_report.terrain_provenance),
          timeline_field("metadata_source", options.d02_start_state_path.string()),
          timeline_field(
              "exposed_mass_cells",
              producer_report.exposed_mass_cell_count),
          timeline_field(
              "pb_written_points",
              producer_report.pb_written_point_count),
          timeline_field(
              "mub_written_cells",
              producer_report.mub_written_cell_count),
          timeline_field(
              "phb_written_points",
              producer_report.phb_written_point_count),
          timeline_field(
              "ht_written_cells",
              producer_report.ht_written_cell_count),
          timeline_field(
              "changed_p_points",
              report.normal_base_state_changed_p_points),
          timeline_field(
              "changed_pb_points",
              report.normal_base_state_changed_pb_points),
          timeline_field(
              "changed_mub_points",
              report.normal_base_state_changed_mub_points),
          timeline_field(
              "changed_phb_points",
              report.normal_base_state_changed_phb_points),
      });
}

struct DiagnosticAdapterExposedRegionSet {
  std::array<tywrf::nest::RemapWindow, 4> regions{};
  std::uint8_t count = 0;
};

void append_diagnostic_adapter_exposed_region(
    DiagnosticAdapterExposedRegionSet& set,
    const tywrf::nest::HorizontalStagger stagger,
    const std::int32_t i_begin,
    const std::int32_t j_begin,
    const std::int32_t extent_i,
    const std::int32_t extent_j) noexcept {
  if (extent_i <= 0 || extent_j <= 0 || set.count >= set.regions.size()) {
    return;
  }

  auto& region = set.regions[set.count];
  region.stagger = stagger;
  region.new_i_begin = i_begin;
  region.new_j_begin = j_begin;
  region.extent_i = extent_i;
  region.extent_j = extent_j;
  ++set.count;
}

[[nodiscard]] DiagnosticAdapterExposedRegionSet
diagnostic_adapter_exposed_regions_from_overlap(
    const tywrf::nest::RemapWindow& overlap,
    const std::int32_t active_nx_value,
    const std::int32_t active_ny_value) noexcept {
  DiagnosticAdapterExposedRegionSet set{};
  const auto overlap_i0 = overlap.new_i_begin;
  const auto overlap_j0 = overlap.new_j_begin;
  const auto overlap_i1 = overlap.new_i_begin + overlap.extent_i;
  const auto overlap_j1 = overlap.new_j_begin + overlap.extent_j;

  append_diagnostic_adapter_exposed_region(
      set, overlap.stagger, 0, 0, active_nx_value, overlap_j0);
  append_diagnostic_adapter_exposed_region(
      set,
      overlap.stagger,
      0,
      overlap_j1,
      active_nx_value,
      active_ny_value - overlap_j1);
  append_diagnostic_adapter_exposed_region(
      set, overlap.stagger, 0, overlap_j0, overlap_i0, overlap.extent_j);
  append_diagnostic_adapter_exposed_region(
      set,
      overlap.stagger,
      overlap_i1,
      overlap_j0,
      active_nx_value - overlap_i1,
      overlap.extent_j);
  return set;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t diagnostic_adapter_active_nx(
    const tywrf::FieldView2D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t diagnostic_adapter_active_ny(
    const tywrf::FieldView2D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t diagnostic_adapter_active_nx(
    const tywrf::FieldView3D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t diagnostic_adapter_active_ny(
    const tywrf::FieldView3D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t diagnostic_adapter_active_nz(
    const tywrf::FieldView3D<Real>& field) noexcept {
  return field.nz - field.halo.k_lower - field.halo.k_upper;
}

template <typename LhsReal, typename RhsReal>
void require_same_diagnostic_adapter_active_shape(
    const tywrf::FieldView2D<LhsReal>& lhs,
    const tywrf::FieldView2D<RhsReal>& rhs,
    const std::string_view label) {
  if (diagnostic_adapter_active_nx(lhs) != diagnostic_adapter_active_nx(rhs) ||
      diagnostic_adapter_active_ny(lhs) != diagnostic_adapter_active_ny(rhs)) {
    throw std::runtime_error(
        "diagnostic adapter source-child delta shape mismatch for " +
        std::string(label));
  }
}

template <typename LhsReal, typename RhsReal>
void require_same_diagnostic_adapter_active_shape(
    const tywrf::FieldView3D<LhsReal>& lhs,
    const tywrf::FieldView3D<RhsReal>& rhs,
    const std::string_view label) {
  if (diagnostic_adapter_active_nx(lhs) != diagnostic_adapter_active_nx(rhs) ||
      diagnostic_adapter_active_ny(lhs) != diagnostic_adapter_active_ny(rhs) ||
      diagnostic_adapter_active_nz(lhs) != diagnostic_adapter_active_nz(rhs)) {
    throw std::runtime_error(
        "diagnostic adapter source-child delta shape mismatch for " +
        std::string(label));
  }
}

void accumulate_diagnostic_adapter_delta(
    DiagnosticAdapterSourceChildFieldDeltaReport& report,
    const float source_value,
    const float child_value) noexcept {
  ++report.compared_value_count;
  double diff = std::numeric_limits<double>::max();
  if (std::isfinite(source_value) && std::isfinite(child_value)) {
    diff = std::fabs(
        static_cast<double>(source_value) - static_cast<double>(child_value));
  }
  if (diff != 0.0) {
    ++report.differing_value_count;
  }
  report.max_abs_diff = std::max(report.max_abs_diff, diff);
}

[[nodiscard]] DiagnosticAdapterSourceChildFieldDeltaReport
diagnostic_adapter_source_child_delta_2d(
    const std::string_view label,
    const DiagnosticAdapterExposedRegionSet& regions,
    const tywrf::FieldView2D<const float>& source,
    const tywrf::FieldView2D<float>& child) {
  require_same_diagnostic_adapter_active_shape(source, child, label);
  DiagnosticAdapterSourceChildFieldDeltaReport report{};
  for (std::uint8_t region_index = 0; region_index < regions.count; ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t i = region.new_i_begin;
           i < region.new_i_begin + region.extent_i;
           ++i) {
        accumulate_diagnostic_adapter_delta(
            report,
            source(source.halo.i_lower + i, source.halo.j_lower + j),
            child(child.halo.i_lower + i, child.halo.j_lower + j));
      }
    }
  }
  return report;
}

[[nodiscard]] DiagnosticAdapterSourceChildFieldDeltaReport
diagnostic_adapter_source_child_delta_3d(
    const std::string_view label,
    const DiagnosticAdapterExposedRegionSet& regions,
    const tywrf::FieldView3D<const float>& source,
    const tywrf::FieldView3D<float>& child) {
  require_same_diagnostic_adapter_active_shape(source, child, label);
  DiagnosticAdapterSourceChildFieldDeltaReport report{};
  for (std::uint8_t region_index = 0; region_index < regions.count; ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t k = 0; k < diagnostic_adapter_active_nz(child); ++k) {
        for (std::int32_t i = region.new_i_begin;
             i < region.new_i_begin + region.extent_i;
             ++i) {
          accumulate_diagnostic_adapter_delta(
              report,
              source(
                  source.halo.i_lower + i,
                  source.halo.j_lower + j,
                  source.halo.k_lower + k),
              child(
                  child.halo.i_lower + i,
                  child.halo.j_lower + j,
                  child.halo.k_lower + k));
        }
      }
    }
  }
  return report;
}

void merge_diagnostic_adapter_delta_field(
    DiagnosticAdapterSourceChildDeltaReport& aggregate,
    const DiagnosticAdapterSourceChildFieldDeltaReport& field) noexcept {
  aggregate.compared_value_count += field.compared_value_count;
  aggregate.differing_value_count += field.differing_value_count;
  aggregate.max_abs_diff = std::max(aggregate.max_abs_diff, field.max_abs_diff);
}

[[nodiscard]] DiagnosticAdapterSourceChildDeltaReport
compare_diagnostic_adapter_source_child_delta(
    const tywrf::nest::RemapPlan& remap_plan,
    const tywrf::nest::ExposedBaseStateViews<const float>& source,
    const tywrf::nest::ExposedBaseStateViews<float>& child) {
  DiagnosticAdapterSourceChildDeltaReport report{};
  const auto mass_regions = diagnostic_adapter_exposed_regions_from_overlap(
      remap_plan.mass,
      diagnostic_adapter_active_nx(child.pb),
      diagnostic_adapter_active_ny(child.pb));
  const auto surface_regions = diagnostic_adapter_exposed_regions_from_overlap(
      remap_plan.surface,
      diagnostic_adapter_active_nx(child.mub),
      diagnostic_adapter_active_ny(child.mub));
  const auto w_full_regions = diagnostic_adapter_exposed_regions_from_overlap(
      remap_plan.w_full,
      diagnostic_adapter_active_nx(child.phb),
      diagnostic_adapter_active_ny(child.phb));

  report.phb = diagnostic_adapter_source_child_delta_3d(
      "PHB", w_full_regions, source.phb, child.phb);
  report.mub = diagnostic_adapter_source_child_delta_2d(
      "MUB", surface_regions, source.mub, child.mub);
  report.ht = diagnostic_adapter_source_child_delta_2d(
      "HT", surface_regions, source.ht, child.ht);
  report.pb = diagnostic_adapter_source_child_delta_3d(
      "PB", mass_regions, source.pb, child.pb);
  report.t_init = diagnostic_adapter_source_child_delta_3d(
      "T_INIT", mass_regions, source.t_init, child.t_init);
  report.alb = diagnostic_adapter_source_child_delta_3d(
      "ALB", mass_regions, source.alb, child.alb);

  merge_diagnostic_adapter_delta_field(report, report.phb);
  merge_diagnostic_adapter_delta_field(report, report.mub);
  merge_diagnostic_adapter_delta_field(report, report.ht);
  merge_diagnostic_adapter_delta_field(report, report.pb);
  merge_diagnostic_adapter_delta_field(report, report.t_init);
  merge_diagnostic_adapter_delta_field(report, report.alb);
  report.values_identical =
      report.compared_value_count > 0 && report.differing_value_count == 0 &&
      report.max_abs_diff == 0.0;
  return report;
}

void run_diagnostic_adapter_report(
    const Options& options,
    const tywrf::State<float>& candidate,
    const StaticFieldSet& output_static,
    CandidateReport& report) {
  tywrf::FieldStorage3D<float> metadata_alb(candidate.grid.mass_layout());
  const auto metadata = tywrf::io::read_krosa_pressure_refresh_inputs(
      options.template_path,
      candidate.grid,
      metadata_alb,
      {.time_index = options.template_time_index});
  require_pressure_refresh_inputs_ready(metadata);

  DiagnosticBaseStateAdapterStaging child(candidate.grid);
  fill_diagnostic_adapter_staging(candidate, output_static, metadata_alb, child);

  const auto provider_source = build_diagnostic_adapter_provider_source(
      candidate,
      output_static,
      metadata.metadata);
  const auto source_views =
      diagnostic_adapter_provider_source_views(provider_source, output_static);
  report.diagnostic_adapter_provider_source = provider_source.report;
  append_timeline_event(
      report,
      "diagnostic_adapter_provider_source",
      {
          timeline_field(
              "origin",
              std::string(diagnostic_adapter_provider_source_origin())),
          timeline_field(
              "source_origin",
              std::string(diagnostic_adapter_provider_source_origin())),
          timeline_field(
              "provider_ok",
              std::string(bool_text(provider_source.report.ok()))),
          timeline_field(
              "diagnostic_only",
              std::string(bool_text(provider_source.report.diagnostic_only))),
          timeline_field(
              "gate_candidate",
              std::string(bool_text(provider_source.report.gate_candidate))),
          timeline_field(
              "integrator_output",
              std::string(bool_text(provider_source.report.integrator_output))),
          timeline_field("writes_candidate", "false"),
          timeline_field("reads_direct_p", "false"),
          timeline_field(
              "terrain_source",
              provider_source.report.terrain_source_name),
          timeline_field(
              "terrain_provenance",
              provider_source.report.terrain_provenance),
          timeline_field(
              "wrote_pb",
              std::string(bool_text(provider_source.report.wrote_pb))),
          timeline_field(
              "wrote_t_init",
              std::string(bool_text(provider_source.report.wrote_t_init))),
          timeline_field(
              "wrote_mub",
              std::string(bool_text(provider_source.report.wrote_mub))),
          timeline_field(
              "wrote_alb",
              std::string(bool_text(provider_source.report.wrote_alb))),
          timeline_field(
              "wrote_phb",
              std::string(bool_text(provider_source.report.wrote_phb))),
          timeline_field(
              "provider_reconstructed_phb_not_wrf_rebalance_validated",
              "true"),
      });

  tywrf::nest::BaseStateSourceStagingProvider source_provider;
  const auto source_staging_report = source_provider.stage(
      candidate.grid,
      report.remap_plan,
      source_views);
  const auto provider_source_views = source_provider.views();
  const auto child_views = diagnostic_adapter_views(child);
  report.diagnostic_adapter_source_staging = source_staging_report;
  report.diagnostic_adapter_source_staging_aliases_child =
      base_state_views_alias(provider_source_views, child_views);
  append_timeline_event(
      report,
      "diagnostic_adapter_source_staging",
      {
          timeline_field("provider", "BaseStateSourceStagingProvider"),
          timeline_field("source", std::string(source_staging_report.source)),
          timeline_field("ok", std::string(bool_text(source_staging_report.ok()))),
          timeline_field(
              "diagnostic_only",
              std::string(bool_text(source_staging_report.diagnostic_only))),
          timeline_field(
              "gate_candidate",
              std::string(bool_text(source_staging_report.gate_candidate))),
          timeline_field(
              "integrator_output",
              std::string(bool_text(source_staging_report.integrator_output))),
          timeline_field(
              "writes_candidate",
              std::string(bool_text(source_staging_report.writes_candidate))),
          timeline_field(
              "writes_netcdf",
              std::string(bool_text(source_staging_report.writes_netcdf))),
          timeline_field(
              "uses_reference_end_truth",
              std::string(bool_text(source_staging_report.uses_reference_end_truth))),
          timeline_field(
              "uses_direct_p_shortcut",
              std::string(bool_text(source_staging_report.uses_direct_p_shortcut))),
          timeline_field(
              "aliases_child",
              std::string(
                  bool_text(report.diagnostic_adapter_source_staging_aliases_child))),
          timeline_field(
              "exposed_regions",
              static_cast<std::uint64_t>(source_staging_report.exposed_region_count)),
          timeline_field(
              "exposed_mass_points",
              source_staging_report.exposed_mass_point_count),
          timeline_field(
              "masked_mass_points",
              source_staging_report.masked_mass_point_count),
          timeline_field("staged_values", source_staging_report.staged_value_count),
          timeline_field(
              "invalid_values",
              source_staging_report.invalid_exposed_value_count),
      });
  if (!source_staging_report.ok()) {
    std::ostringstream message;
    message << "diagnostic_base_state_source_staging_failed";
    if (source_staging_report.result.message != nullptr &&
        source_staging_report.result.message[0] != '\0') {
      message << ": " << source_staging_report.result.message;
    }
    throw std::runtime_error(message.str());
  }
  if (report.diagnostic_adapter_source_staging_aliases_child) {
    throw std::runtime_error(
        "diagnostic_base_state_source_staging unexpectedly aliases child staging");
  }

  const auto source_child_delta =
      compare_diagnostic_adapter_source_child_delta(
          report.remap_plan, provider_source_views, child_views);
  report.diagnostic_adapter_source_child_delta = source_child_delta;
  append_timeline_event(
      report,
      "diagnostic_adapter_source_child_delta",
      {
          timeline_field("diagnostic_only", "true"),
          timeline_field("gate_candidate", "false"),
          timeline_field("integrator_output", "false"),
          timeline_field("writes_candidate", "false"),
          timeline_field("writes_netcdf", "false"),
          timeline_field(
              "values_identical",
              std::string(bool_text(source_child_delta.values_identical))),
          timeline_field("compared_values", source_child_delta.compared_value_count),
          timeline_field("differing_values", source_child_delta.differing_value_count),
          timeline_field(
              "max_abs_diff",
              format_probe_double(source_child_delta.max_abs_diff)),
      });

  const auto adapter_report =
      tywrf::nest::apply_exposed_base_state_exchange_adapter(
          report.remap_plan,
          provider_source_views,
          child_views,
          {{metadata.metadata.c3h.data(),
            static_cast<std::int32_t>(metadata.metadata.c3h.size())},
           {metadata.metadata.c4h.data(),
            static_cast<std::int32_t>(metadata.metadata.c4h.size())},
           metadata.metadata.p_top_pa});

  report.diagnostic_adapter = adapter_report;
  report.diagnostic_adapter_metadata_source = options.template_path;
  report.diagnostic_adapter_metadata_time_index = options.template_time_index;
  append_timeline_event(
      report,
      "diagnostic_adapter_report",
      {
          timeline_field("opt_in", "true"),
          timeline_field(
              "diagnostic_only",
              std::string(bool_text(adapter_report.diagnostic_only))),
          timeline_field(
              "gate_candidate",
              std::string(bool_text(adapter_report.gate_candidate))),
          timeline_field(
              "integrator_output",
              std::string(bool_text(adapter_report.integrator_output))),
          timeline_field("ok", std::string(bool_text(adapter_report.ok()))),
          timeline_field(
              "called_d68",
              std::string(bool_text(adapter_report.called_d68_exchange))),
          timeline_field(
              "called_d69",
              std::string(bool_text(adapter_report.called_d69_recompute))),
          timeline_field(
              "writes_candidate",
              std::string(bool_text(adapter_report.writes_candidate))),
          timeline_field(
              "writes_netcdf",
              std::string(bool_text(adapter_report.writes_netcdf))),
          timeline_field(
              "exposed_regions",
              static_cast<std::uint64_t>(adapter_report.exposed_region_count)),
          timeline_field("exposed_mass_cells", adapter_report.exposed_mass_cell_count),
          timeline_field("recomputed_points", adapter_report.recomputed_point_count),
      });
  if (adapter_report.ok()) {
    return;
  }

  std::ostringstream message;
  message << "diagnostic_base_state_adapter_report_failed";
  if (adapter_report.result.message != nullptr &&
      adapter_report.result.message[0] != '\0') {
    message << ": " << adapter_report.result.message;
  }
  throw std::runtime_error(message.str());
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
  append_timeline_event(
      report,
      "cycle_start",
      {
          timeline_field("cycle_start", options.cycle_start),
          timeline_field("cycle_end", options.cycle_end),
      });
  const auto parent_delta_i =
      static_cast<std::int64_t>(options.to_parent_start.i_parent_start) -
      static_cast<std::int64_t>(options.from_parent_start.i_parent_start);
  const auto parent_delta_j =
      static_cast<std::int64_t>(options.to_parent_start.j_parent_start) -
      static_cast<std::int64_t>(options.from_parent_start.j_parent_start);
  append_timeline_event(
      report,
      "move_from_to_parent_start",
      {
          timeline_field("from_i", static_cast<std::int64_t>(options.from_parent_start.i_parent_start)),
          timeline_field("from_j", static_cast<std::int64_t>(options.from_parent_start.j_parent_start)),
          timeline_field("to_i", static_cast<std::int64_t>(options.to_parent_start.i_parent_start)),
          timeline_field("to_j", static_cast<std::int64_t>(options.to_parent_start.j_parent_start)),
          timeline_field("parent_grid_ratio", static_cast<std::int64_t>(descriptor.parent_grid_ratio)),
          timeline_field("parent_delta_i", parent_delta_i),
          timeline_field("parent_delta_j", parent_delta_j),
          timeline_field(
              "child_delta_i",
              parent_delta_i * static_cast<std::int64_t>(descriptor.parent_grid_ratio)),
          timeline_field(
              "child_delta_j",
              parent_delta_j * static_cast<std::int64_t>(descriptor.parent_grid_ratio)),
      });
  report.remap_plan = remap_plan;
  report.remap = tywrf::nest::remap_child_state_overlap_only(remap_plan, d02_start, candidate);
  if (!report.remap.ok()) {
    throw std::runtime_error("failed to remap d02 overlap: " + std::string(report.remap.result.message));
  }
  append_timeline_event(
      report,
      "overlap_remap",
      {
          timeline_field("copied_points", report.remap.copied_point_count),
          timeline_field("copied_fields", static_cast<std::uint64_t>(report.remap.copied_field_count)),
          timeline_field("child_delta_i", static_cast<std::int64_t>(report.remap_plan.delta.child_di)),
          timeline_field("child_delta_j", static_cast<std::int64_t>(report.remap_plan.delta.child_dj)),
          timeline_field("needs_parent_fill", report.remap.needs_parent_fill ? "true" : "false"),
          timeline_field(
              "needs_derived_pressure_refresh",
              report.remap.needs_derived_pressure_refresh ? "true" : "false"),
      });

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
  append_timeline_event(
      report,
      "exchange_plan_build",
      {
          timeline_field("exchange_points", report.exchange.report.exchange_point_count),
          timeline_field(
              "requires_parent_interpolation",
              report.exchange.report.requires_parent_interpolation ? "true" : "false"),
          timeline_field(
              "modifies_overlap", report.exchange.report.modifies_overlap ? "true" : "false"),
          timeline_field("modifies_halo", report.exchange.report.modifies_halo ? "true" : "false"),
      });

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
  append_timeline_event(
      report,
      "parent_interpolation",
      {
          timeline_field("interpolated_points", report.interpolation.interpolated_point_count),
          timeline_field(
              "wrote_overlap", report.interpolation.wrote_overlap ? "true" : "false"),
          timeline_field("wrote_halo", report.interpolation.wrote_halo ? "true" : "false"),
      });

  apply_selected_field_wind_tendency(options, candidate, report);

  report.changed_selected_points = changed_selected_points(d02_start, candidate);
  if (report.changed_selected_points == 0) {
    throw std::runtime_error("selected-field candidate did not change any selected-field point");
  }
  append_timeline_event(
      report,
      "selected_field_change_summary",
      {
          timeline_field("changed_points", report.changed_selected_points),
      });
  require_finite_strict_fields(candidate);
  return report;
}

void refresh_static_fields(
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const Options& options,
    const StaticFieldSet& d02_start_static,
    const StaticFieldSet& template_static,
    const tywrf::FieldStorage2D<float>& d01_start_hgt,
    StaticFieldSet& output_static,
    CandidateReport& report) {
  output_static = template_static;
  report.static_refresh = tywrf::nest::refresh_moving_nest_static_fields(
      descriptor,
      options.to_parent_start,
      report.remap_plan,
      d02_start_static.xlat.view(),
      d02_start_static.xlong.view(),
      d02_start_static.hgt.view(),
      d01_start_hgt.view(),
      output_static.xlat.view(),
      output_static.xlong.view(),
      output_static.hgt.view());
  if (!report.static_refresh.ok()) {
    throw std::runtime_error(
        "failed moving-nest static refresh: " +
        std::string(report.static_refresh.result.message));
  }
  require_finite_static_fields(output_static);
  report.changed_static_template_points = changed_static_points(template_static, output_static);
  append_timeline_event(
      report,
      "static_refresh",
      {
          timeline_field("overlap_cells", report.static_refresh.overlap_cell_count),
          timeline_field("exposed_cells", report.static_refresh.exposed_cell_count),
          timeline_field(
              "coord_extrapolated_cells",
              report.static_refresh.coordinate_extrapolated_cell_count),
          timeline_field(
              "hgt_parent_interpolated_cells",
              report.static_refresh.parent_hgt_interpolated_cell_count),
          timeline_field("changed_template_points", report.changed_static_template_points),
          timeline_field("uses_reference_end", "false"),
      });

  const auto layout = output_static.xlat.layout();
  const auto center_i = layout.i_begin() + layout.active_nx() / 2;
  const auto center_j = layout.j_begin() + layout.active_ny() / 2;
  const auto xlat_view = output_static.xlat.view();
  const auto xlong_view = output_static.xlong.view();
  report.cen_lat = xlat_view(center_i, center_j);
  report.cen_lon = xlong_view(center_i, center_j);
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

void write_int_attr(
    const NetcdfHandle& file,
    const std::string_view name,
    const std::int32_t value) {
  file.check(
      nc_put_att_int(file.id(), NC_GLOBAL, std::string(name).c_str(), NC_INT, 1, &value),
      "write global integer attribute");
}

void write_pressure_refresh_readiness_attrs(
    const NetcdfHandle& file,
    const PressureRefreshReadiness& readiness) {
  write_text_attr(file, "TYWRF_PRESSURE_REFRESH_READINESS_READY", bool_text(readiness.ready()));
  write_text_attr(
      file,
      "TYWRF_THERMODYNAMIC_BASE_STATE_CONSISTENCY_READY",
      bool_text(readiness.thermodynamic_base_state_consistency_ready));
  write_text_attr(
      file,
      "TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT",
      bool_text(readiness.provider_terrain_uses_moved_candidate_hgt));
  write_text_attr(
      file,
      "TYWRF_PROVIDER_BASE_STATE_RECONSTRUCT_OK",
      bool_text(readiness.provider_base_state_reconstruct_ok));
  write_text_attr(
      file,
      "TYWRF_BASE_STATE_SYNC_CONTRACT_OK",
      bool_text(readiness.base_state_sync_contract_ok));
  write_text_attr(
      file,
      "TYWRF_BASE_STATE_SOURCE_SYNC_READINESS_CHECK",
      bool_text(readiness.base_state_sync_dry_run));
  write_text_attr(
      file,
      "TYWRF_BASE_STATE_SOURCE_SYNC_APPLIED",
      bool_text(readiness.base_state_sync_applied));
  write_double_attr(
      file,
      "TYWRF_SOURCE_SYNC_PLANNED_PB_POINT_COUNT",
      static_cast<double>(readiness.would_sync_pb_point_count));
  write_double_attr(
      file,
      "TYWRF_SOURCE_SYNC_PLANNED_MUB_POINT_COUNT",
      static_cast<double>(readiness.would_sync_mub_point_count));
  write_double_attr(
      file,
      "TYWRF_SOURCE_SYNC_PLANNED_PHB_POINT_COUNT",
      static_cast<double>(readiness.would_sync_phb_point_count));
  write_double_attr(
      file,
      "TYWRF_SOURCE_SYNC_OVERLAP_WRITE_COUNT",
      static_cast<double>(readiness.sync_overlap_write_count));
  write_double_attr(
      file,
      "TYWRF_SOURCE_SYNC_HALO_WRITE_COUNT",
      static_cast<double>(readiness.sync_halo_write_count));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_REFRESH_READINESS_COMPUTE_CALLED",
      bool_text(readiness.pressure_refresh_compute_called));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK",
      bool_text(readiness.pressure_compute_dry_run));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK_CALLED",
      bool_text(readiness.pressure_compute_dry_run_called));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_CHECK_OK",
      bool_text(readiness.pressure_compute_dry_run_ok));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_REFRESH_PLANNED_P_POINT_COUNT",
      static_cast<double>(readiness.would_refresh_p_point_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_READINESS_INVALID_P_POINT_COUNT",
      static_cast<double>(readiness.dry_run_invalid_p_point_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_READINESS_SKIPPED_P_POINT_COUNT",
      static_cast<double>(readiness.dry_run_skipped_p_point_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_TARGET_COLUMN_COUNT",
      static_cast<double>(readiness.pressure_compute_dry_run_report_target_column_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_REFRESHED_POINT_COUNT",
      static_cast<double>(readiness.pressure_compute_dry_run_report_refreshed_point_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_INVALID_POINT_COUNT",
      static_cast<double>(readiness.pressure_compute_dry_run_report_invalid_point_count));
  write_double_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_SKIPPED_POINT_COUNT",
      static_cast<double>(readiness.pressure_compute_dry_run_report_skipped_point_count));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_TOUCHED_OVERLAP_CELLS",
      bool_text(readiness.pressure_compute_dry_run_report_touched_overlap_cells));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COMPUTE_READINESS_REPORT_TOUCHED_HALO_CELLS",
      bool_text(readiness.pressure_compute_dry_run_report_touched_halo_cells));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_REFRESH_READINESS_APPLIED",
      bool_text(readiness.pressure_refresh_applied));
  if (!readiness.provider_terrain_source_name.empty()) {
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_SOURCE",
        readiness.provider_terrain_source_name);
  }
  if (!readiness.provider_terrain_provenance.empty()) {
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_READINESS_PROVIDER_TERRAIN_PROVENANCE",
        readiness.provider_terrain_provenance);
  }
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

void write_pressure_column_probe_attrs(
    const NetcdfHandle& file,
    const Options& options,
    const CandidateReport& report) {
  if (report.pressure_column_observations.empty()) {
    return;
  }

  const auto phases = pressure_column_probe_phase_names(report.pressure_column_observations);
  write_text_attr(file, "TYWRF_PRESSURE_COLUMN_PROBE_VERSION", "runtime_v0");
  write_text_attr(file, "TYWRF_PRESSURE_COLUMN_PROBE_ENABLED", "true");
  write_text_attr(file, "TYWRF_PRESSURE_COLUMN_PROBE_EVIDENCE_ONLY", "true");
  write_text_attr(file, "TYWRF_PRESSURE_COLUMN_PROBE_INDEX_BASE", "zero_based_mass_grid");
  write_int_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_COLUMN_COUNT",
      static_cast<std::int32_t>(options.pressure_column_probe_columns.size()));
  write_int_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_LEVEL_COUNT",
      static_cast<std::int32_t>(options.pressure_column_probe_levels.size()));
  write_int_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_PHASE_COUNT",
      static_cast<std::int32_t>(phases.size()));
  write_int_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_RECORD_COUNT",
      static_cast<std::int32_t>(report.pressure_column_observations.size()));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_COLUMNS",
      join_pressure_probe_columns(options.pressure_column_probe_columns));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_LEVELS",
      join_pressure_probe_levels(options.pressure_column_probe_levels));
  write_text_attr(file, "TYWRF_PRESSURE_COLUMN_PROBE_PHASES", join_variables(phases));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_FIELDS",
      join_variables(pressure_column_probe_field_names()));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_NOT_AVAILABLE",
      join_variables(pressure_column_probe_unavailable_names()));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_COLUMN_PROBE_VALUES",
      pressure_column_probe_values(report.pressure_column_observations));
}

void write_pressure_formula_observation_attrs(
    const NetcdfHandle& file,
    const CandidateReport& report) {
  if (report.pressure_formula_observations.empty() ||
      !report.pressure_refresh.has_value()) {
    return;
  }

  const auto& compute = report.pressure_refresh->compute_report;
  write_text_attr(file, "TYWRF_PRESSURE_FORMULA_OBSERVATION_VERSION", "runtime_v0");
  write_text_attr(file, "TYWRF_PRESSURE_FORMULA_OBSERVATION_ENABLED", "true");
  write_text_attr(file, "TYWRF_PRESSURE_FORMULA_OBSERVATION_EVIDENCE_ONLY", "true");
  write_text_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_INDEX_BASE",
      "zero_based_mass_grid");
  write_int_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_REQUEST_COUNT",
      static_cast<std::int32_t>(compute.observation_request_count));
  write_int_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_RECORD_COUNT",
      static_cast<std::int32_t>(compute.observation_record_count));
  write_int_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALID_COUNT",
      static_cast<std::int32_t>(compute.observation_valid_count));
  write_int_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_INVALID_COUNT",
      static_cast<std::int32_t>(compute.observation_invalid_count));
  write_int_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_OUT_OF_BOUNDS_COUNT",
      static_cast<std::int32_t>(compute.observation_out_of_bounds_count));
  write_int_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_OUTSIDE_TARGET_REGION_COUNT",
      static_cast<std::int32_t>(
          compute.observation_outside_target_region_count));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_FIELDS",
      join_variables(pressure_formula_observation_field_names()));
  write_text_attr(
      file,
      "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALUES",
      pressure_formula_observation_values(report.pressure_formula_observations));
}

void write_normal_base_state_producer_attrs(
    const NetcdfHandle& file,
    const CandidateReport& report) {
  if (!report.normal_base_state_producer.has_value()) {
    return;
  }

  const auto& producer = *report.normal_base_state_producer;
  const auto* provider = report.normal_base_state_provider_source.has_value()
                             ? &*report.normal_base_state_provider_source
                             : nullptr;
  write_text_attr(file, "TYWRF_NORMAL_BASE_STATE_PRODUCER_VERSION", "a79_normal_v0");
  write_text_attr(file, "TYWRF_NORMAL_BASE_STATE_PRODUCER_SOURCE", producer.source);
  write_text_attr(file, "TYWRF_NORMAL_BASE_STATE_PRODUCER_DISPOSITION", producer.disposition);
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_SOURCE_ORIGIN",
      producer.source_origin);
  write_text_attr(file, "TYWRF_NORMAL_BASE_STATE_PRODUCER_OK", bool_text(producer.ok()));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_DIAGNOSTIC_ONLY",
      bool_text(producer.diagnostic_only));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_NORMAL_CANDIDATE_PRODUCER",
      bool_text(producer.normal_candidate_producer));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_GATE_CANDIDATE",
      bool_text(producer.gate_candidate));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_INTEGRATOR_OUTPUT",
      bool_text(producer.integrator_output));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_WRITES_CANDIDATE",
      bool_text(producer.writes_candidate));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_WRITES_NETCDF",
      bool_text(producer.writes_netcdf));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_USES_REFERENCE_END_TRUTH",
      bool_text(producer.uses_reference_end_truth));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_NO_REFERENCE_END_TRUTH",
      bool_text(!producer.uses_reference_end_truth));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_USES_DIRECT_P_SHORTCUT",
      bool_text(producer.uses_direct_p_shortcut));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_NO_DIRECT_P_SHORTCUT",
      bool_text(!producer.uses_direct_p_shortcut));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_READS_DIRECT_P",
      bool_text(producer.reads_direct_p));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_WRITES_P",
      bool_text(producer.writes_p));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_NO_GATE_PASS_CLAIM",
      bool_text(producer.no_gate_pass_claim));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_NO_00_20_PROGRESSION",
      bool_text(producer.no_00_20_progression));
  write_text_attr(
      file, "TYWRF_NORMAL_BASE_STATE_PRODUCER_WRITTEN_FIELDS", producer.written_fields);
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_METADATA_SOURCE",
      report.normal_base_state_metadata_source.string());
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_METADATA_TIME_INDEX",
      static_cast<double>(report.normal_base_state_metadata_time_index));
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_PROVIDER_SOURCE",
      provider == nullptr ? "" : provider->source);
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_TERRAIN_SOURCE",
      provider == nullptr ? "" : provider->terrain_source_name);
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_TERRAIN_PROVENANCE",
      provider == nullptr ? "" : provider->terrain_provenance);
  write_text_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_TERRAIN_OVERRIDE_USED",
      bool_text(provider != nullptr && provider->terrain_override_used));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_EXPOSED_MASS_CELL_COUNT",
      static_cast<double>(producer.exposed_mass_cell_count));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_EXPOSED_SURFACE_CELL_COUNT",
      static_cast<double>(producer.exposed_surface_cell_count));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_EXPOSED_W_FULL_COLUMN_COUNT",
      static_cast<double>(producer.exposed_w_full_column_count));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_PB_WRITTEN_POINT_COUNT",
      static_cast<double>(producer.pb_written_point_count));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_MUB_WRITTEN_CELL_COUNT",
      static_cast<double>(producer.mub_written_cell_count));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_PHB_WRITTEN_POINT_COUNT",
      static_cast<double>(producer.phb_written_point_count));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_HT_WRITTEN_CELL_COUNT",
      static_cast<double>(producer.ht_written_cell_count));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_DIRECT_WRITE_POINT_COUNT",
      static_cast<double>(producer.direct_write_point_count));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_P_POINTS",
      static_cast<double>(report.normal_base_state_changed_p_points));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_PB_POINTS",
      static_cast<double>(report.normal_base_state_changed_pb_points));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_MUB_POINTS",
      static_cast<double>(report.normal_base_state_changed_mub_points));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_PHB_POINTS",
      static_cast<double>(report.normal_base_state_changed_phb_points));
  write_double_attr(
      file,
      "TYWRF_NORMAL_BASE_STATE_PRODUCER_CHANGED_HGT_POINTS",
      static_cast<double>(report.normal_base_state_changed_hgt_points));
}

void write_diagnostic_adapter_provider_source_attrs(
    const NetcdfHandle& file,
    const CandidateReport& report) {
  if (!report.diagnostic_adapter_provider_source.has_value()) {
    return;
  }

  const auto& provider = *report.diagnostic_adapter_provider_source;
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_VERSION",
      "d77_provider_source_v0");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_ORIGIN",
      diagnostic_adapter_provider_source_origin());
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_SOURCE_ORIGIN",
      diagnostic_adapter_provider_source_origin());
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_PROVIDER_SOURCE",
      provider.source);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_TERRAIN_SOURCE",
      provider.terrain_source_name);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_TERRAIN_PROVENANCE",
      provider.terrain_provenance);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_TERRAIN_OVERRIDE_USED",
      bool_text(provider.terrain_override_used));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_HT_SOURCE",
      "output_static.hgt");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_HT_PROVENANCE",
      "adapter_HT_from_output_static_HGT");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_PROVIDER_OK",
      bool_text(provider.ok()));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_RESULT_MESSAGE",
      provider.result.message == nullptr ? "" : provider.result.message);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_DIAGNOSTIC_ONLY",
      bool_text(provider.diagnostic_only));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_GATE_CANDIDATE",
      bool_text(provider.gate_candidate));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_INTEGRATOR_OUTPUT",
      bool_text(provider.integrator_output));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WRITES_CANDIDATE",
      "false");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WRITES_NETCDF",
      "false");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_NO_CANDIDATE_WRITE",
      "true");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_USES_REFERENCE_END_TRUTH",
      "false");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_NO_REFERENCE_END_TRUTH",
      "true");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_USES_DIRECT_P_SHORTCUT",
      "false");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_NO_DIRECT_P_SHORTCUT",
      "true");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_READS_DIRECT_P",
      "false");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_PB",
      bool_text(provider.wrote_pb));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_T_INIT",
      bool_text(provider.wrote_t_init));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_MUB",
      bool_text(provider.wrote_mub));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_ALB",
      bool_text(provider.wrote_alb));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_WROTE_PHB",
      bool_text(provider.wrote_phb));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_PROVIDER_RECONSTRUCTED_PHB_NOT_WRF_REBALANCE_VALIDATED",
      "true");
}

void write_diagnostic_adapter_source_staging_attrs(
    const NetcdfHandle& file,
    const CandidateReport& report) {
  if (!report.diagnostic_adapter_source_staging.has_value()) {
    return;
  }

  const auto& staging = *report.diagnostic_adapter_source_staging;
  write_text_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_VERSION", "d75_provider_report_v0");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_PROVIDER_KIND",
      "BaseStateSourceStagingProvider");
  write_text_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_SOURCE", staging.source);
  write_text_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_DISPOSITION", staging.disposition);
  write_text_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_SOURCE_SHAPE", staging.source_shape);
  write_text_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_OK", bool_text(staging.ok()));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_RESULT_MESSAGE",
      staging.result.message == nullptr ? "" : staging.result.message);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_DIAGNOSTIC_ONLY",
      bool_text(staging.diagnostic_only));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_GATE_CANDIDATE",
      bool_text(staging.gate_candidate));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INTEGRATOR_OUTPUT",
      bool_text(staging.integrator_output));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_WRITES_CANDIDATE",
      bool_text(staging.writes_candidate));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_WRITES_NETCDF",
      bool_text(staging.writes_netcdf));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_CANDIDATE_BUFFERS_PRESERVED",
      bool_text(staging.candidate_buffers_preserved));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_OWNS_STAGING_BUFFERS",
      bool_text(staging.owns_staging_buffers));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_ALLOCATED_BUFFERS",
      bool_text(staging.allocated_buffers));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_USES_REFERENCE_END_TRUTH",
      bool_text(staging.uses_reference_end_truth));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_USES_DIRECT_P_SHORTCUT",
      bool_text(staging.uses_direct_p_shortcut));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_READS_DIRECT_P",
      bool_text(staging.reads_direct_p));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_ALIASES_CHILD",
      bool_text(report.diagnostic_adapter_source_staging_aliases_child));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_ACTIVE_NX",
      static_cast<double>(staging.active_nx));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_ACTIVE_NY",
      static_cast<double>(staging.active_ny));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_MASS_NZ",
      static_cast<double>(staging.mass_nz));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_FULL_NZ",
      static_cast<double>(staging.full_nz));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_REGION_COUNT",
      static_cast<double>(staging.exposed_region_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_MASS_CELL_COUNT",
      static_cast<double>(staging.exposed_mass_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_MASS_POINT_COUNT",
      static_cast<double>(staging.exposed_mass_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_SURFACE_CELL_COUNT",
      static_cast<double>(staging.exposed_surface_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_W_FULL_COLUMN_COUNT",
      static_cast<double>(staging.exposed_w_full_column_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_EXPOSED_W_FULL_POINT_COUNT",
      static_cast<double>(staging.exposed_w_full_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_MASKED_MASS_CELL_COUNT",
      static_cast<double>(staging.masked_mass_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_MASKED_MASS_POINT_COUNT",
      static_cast<double>(staging.masked_mass_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_MASKED_SURFACE_CELL_COUNT",
      static_cast<double>(staging.masked_surface_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_MASKED_W_FULL_COLUMN_COUNT",
      static_cast<double>(staging.masked_w_full_column_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_MASKED_W_FULL_POINT_COUNT",
      static_cast<double>(staging.masked_w_full_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_PHB_POINT_COUNT",
      static_cast<double>(staging.staged_phb_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_MUB_CELL_COUNT",
      static_cast<double>(staging.staged_mub_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_HT_CELL_COUNT",
      static_cast<double>(staging.staged_ht_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_PB_POINT_COUNT",
      static_cast<double>(staging.staged_pb_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_T_INIT_POINT_COUNT",
      static_cast<double>(staging.staged_t_init_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_ALB_POINT_COUNT",
      static_cast<double>(staging.staged_alb_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_STAGED_VALUE_COUNT",
      static_cast<double>(staging.staged_value_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_ACTIVE_MASKED_VALUE_COUNT",
      static_cast<double>(staging.active_masked_value_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_HALO_MASKED_VALUE_COUNT",
      static_cast<double>(staging.halo_masked_value_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INVALID_EXPOSED_VALUE_COUNT",
      static_cast<double>(staging.invalid_exposed_value_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INVALID_EXPOSED_PHB_POINT_COUNT",
      static_cast<double>(staging.invalid_exposed_phb_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INVALID_EXPOSED_MUB_CELL_COUNT",
      static_cast<double>(staging.invalid_exposed_mub_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INVALID_EXPOSED_HT_CELL_COUNT",
      static_cast<double>(staging.invalid_exposed_ht_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INVALID_EXPOSED_PB_POINT_COUNT",
      static_cast<double>(staging.invalid_exposed_pb_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INVALID_EXPOSED_T_INIT_POINT_COUNT",
      static_cast<double>(staging.invalid_exposed_t_init_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_INVALID_EXPOSED_ALB_POINT_COUNT",
      static_cast<double>(staging.invalid_exposed_alb_point_count));
}

void write_diagnostic_adapter_source_child_delta_field_attrs(
    const NetcdfHandle& file,
    const std::string_view field,
    const DiagnosticAdapterSourceChildFieldDeltaReport& delta) {
  const std::string prefix =
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_" + std::string(field);
  write_double_attr(
      file,
      prefix + "_COMPARED_VALUE_COUNT",
      static_cast<double>(delta.compared_value_count));
  write_double_attr(
      file,
      prefix + "_DIFFERING_VALUE_COUNT",
      static_cast<double>(delta.differing_value_count));
  write_double_attr(file, prefix + "_MAX_ABS_DIFF", delta.max_abs_diff);
}

void write_diagnostic_adapter_source_child_delta_attrs(
    const NetcdfHandle& file,
    const CandidateReport& report) {
  if (!report.diagnostic_adapter_source_child_delta.has_value()) {
    return;
  }

  const auto& delta = *report.diagnostic_adapter_source_child_delta;
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_VERSION",
      "a76_source_child_delta_v0");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_SOURCE",
      "BaseStateSourceStagingProvider_vs_child_staging_pre_adapter");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_SCOPE",
      "exposed_base_state_values_only");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_FIELDS",
      "PHB,MUB,HT,PB,T_INIT,ALB");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_DIAGNOSTIC_ONLY",
      bool_text(delta.diagnostic_only));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_GATE_CANDIDATE",
      bool_text(delta.gate_candidate));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_INTEGRATOR_OUTPUT",
      bool_text(delta.integrator_output));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_WRITES_CANDIDATE",
      bool_text(delta.writes_candidate));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_WRITES_NETCDF",
      bool_text(delta.writes_netcdf));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_VALUES_IDENTICAL",
      bool_text(delta.values_identical));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_COMPARED_VALUE_COUNT",
      static_cast<double>(delta.compared_value_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_DIFFERING_VALUE_COUNT",
      static_cast<double>(delta.differing_value_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_MAX_ABS_DIFF",
      delta.max_abs_diff);
  write_diagnostic_adapter_source_child_delta_field_attrs(file, "PHB", delta.phb);
  write_diagnostic_adapter_source_child_delta_field_attrs(file, "MUB", delta.mub);
  write_diagnostic_adapter_source_child_delta_field_attrs(file, "HT", delta.ht);
  write_diagnostic_adapter_source_child_delta_field_attrs(file, "PB", delta.pb);
  write_diagnostic_adapter_source_child_delta_field_attrs(file, "T_INIT", delta.t_init);
  write_diagnostic_adapter_source_child_delta_field_attrs(file, "ALB", delta.alb);
}

void write_diagnostic_adapter_attrs(
    const NetcdfHandle& file,
    const CandidateReport& report) {
  if (!report.diagnostic_adapter.has_value()) {
    return;
  }

  const auto& adapter = *report.diagnostic_adapter;
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ADAPTER_OPT_IN", "true");
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ADAPTER_VERSION", "d70_report_v0");
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE", adapter.source);
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ADAPTER_DISPOSITION", adapter.disposition);
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ADAPTER_EXCHANGE_SOURCE", adapter.exchange_source);
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ADAPTER_RECOMPUTE_SOURCE", adapter.recompute_source);
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ADAPTER_OK", bool_text(adapter.ok()));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_RESULT_MESSAGE",
      adapter.result.message == nullptr ? "" : adapter.result.message);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_DIAGNOSTIC_ONLY",
      bool_text(adapter.diagnostic_only));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_GATE_CANDIDATE",
      bool_text(adapter.gate_candidate));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_INTEGRATOR_OUTPUT",
      bool_text(adapter.integrator_output));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_SELECTED_FIELD_NUMERICS_ENABLED",
      bool_text(adapter.selected_field_numerics_enabled));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_ENABLES_SELECTED_FIELD_NUMERICS",
      bool_text(adapter.enables_selected_field_numerics));
  write_text_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_WRITES_NETCDF", bool_text(adapter.writes_netcdf));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_WRITES_CANDIDATE",
      bool_text(adapter.writes_candidate));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_CALLED_D68_EXCHANGE",
      bool_text(adapter.called_d68_exchange));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_CALLED_D69_RECOMPUTE",
      bool_text(adapter.called_d69_recompute));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_HT_SOURCE_NAME",
      adapter.ht_source_name);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_HT_DIAGNOSTIC_LABEL",
      adapter.ht_diagnostic_label);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_HT_IS_HGT_ALIAS",
      bool_text(adapter.ht_is_hgt_alias));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_CREATES_TERRAIN_OWNER",
      bool_text(adapter.creates_terrain_owner));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_TERRAIN_OWNER_CREATED",
      bool_text(adapter.terrain_owner_created));
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_METADATA_SOURCE",
      report.diagnostic_adapter_metadata_source.string());
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_METADATA_TIME_INDEX",
      static_cast<double>(report.diagnostic_adapter_metadata_time_index));
  write_double_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_ACTIVE_NX", static_cast<double>(adapter.active_nx));
  write_double_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_ACTIVE_NY", static_cast<double>(adapter.active_ny));
  write_double_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_ACTIVE_NZ", static_cast<double>(adapter.active_nz));
  write_double_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_C3H_COUNT", static_cast<double>(adapter.c3h_count));
  write_double_attr(
      file, "TYWRF_DIAGNOSTIC_ADAPTER_C4H_COUNT", static_cast<double>(adapter.c4h_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_EXPOSED_REGION_COUNT",
      static_cast<double>(adapter.exposed_region_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_EXPOSED_MASS_CELL_COUNT",
      static_cast<double>(adapter.exposed_mass_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_EXPOSED_SURFACE_CELL_COUNT",
      static_cast<double>(adapter.exposed_surface_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_EXPOSED_W_FULL_COLUMN_COUNT",
      static_cast<double>(adapter.exposed_w_full_column_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PHB_WRITTEN_POINT_COUNT",
      static_cast<double>(adapter.phb_written_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_MUB_WRITTEN_CELL_COUNT",
      static_cast<double>(adapter.mub_written_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_HT_WRITTEN_CELL_COUNT",
      static_cast<double>(adapter.ht_written_cell_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_DIRECT_WRITE_POINT_COUNT",
      static_cast<double>(adapter.direct_write_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_EXCHANGE_RECOMPUTE_MARK_COUNT",
      static_cast<double>(adapter.exchange_recompute_mark_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_RECOMPUTED_POINT_COUNT",
      static_cast<double>(adapter.recomputed_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_PB_RECOMPUTED_POINT_COUNT",
      static_cast<double>(adapter.pb_recomputed_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_T_INIT_RECOMPUTED_POINT_COUNT",
      static_cast<double>(adapter.t_init_recomputed_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_ALB_RECOMPUTED_POINT_COUNT",
      static_cast<double>(adapter.alb_recomputed_point_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_INVALID_COLUMN_COUNT",
      static_cast<double>(adapter.invalid_column_count));
  write_double_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_INVALID_POINT_COUNT",
      static_cast<double>(adapter.invalid_point_count));
  write_diagnostic_adapter_provider_source_attrs(file, report);
  write_diagnostic_adapter_source_staging_attrs(file, report);
  write_diagnostic_adapter_source_child_delta_attrs(file, report);
  write_text_attr(file, "TYWRF_DIAGNOSTIC_ADAPTER_WRITTEN_FIELDS", "PHB,MUB,HGT,PB,T_INIT,ALB");
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ADAPTER_INTEGRATION_STATUS",
      "staging_report_only_no_gate_no_integrator");
}

void write_wind_tendency_attrs(
    const NetcdfHandle& file,
    const CandidateReport& report) {
  if (!report.wind_tendency.has_value()) {
    return;
  }

  const auto& wind = *report.wind_tendency;
  const bool gate_evidence = wind_tendency_gate_evidence(wind.source_kind);
  const bool zero_or_identity_only = wind_tendency_zero_or_identity_only(wind.source_kind);
  write_text_attr(file, "TYWRF_WIND_TENDENCY_OPT_IN", "true");
  write_text_attr(file, "TYWRF_WIND_TENDENCY_APPLIED", "true");
  write_text_attr(
      file,
      "TYWRF_WIND_TENDENCY_SOURCE_KIND",
      wind_tendency_source_name(wind.source_kind));
  write_text_attr(
      file, "TYWRF_WIND_TENDENCY_GATE_EVIDENCE", gate_evidence ? "true" : "false");
  write_text_attr(
      file,
      "TYWRF_WIND_TENDENCY_VALIDATION_GATE_EVIDENCE",
      gate_evidence ? "true" : "false");
  write_text_attr(file, "TYWRF_WIND_TENDENCY_USES_REFERENCE_END_TRUTH", "false");
  write_text_attr(
      file,
      "TYWRF_WIND_TENDENCY_ZERO_OR_IDENTITY_ONLY",
      zero_or_identity_only ? "true" : "false");
  write_text_attr(file, "TYWRF_WIND_TENDENCY_WRITTEN_FIELDS", "U,V");
  write_text_attr(
      file,
      "TYWRF_WIND_TENDENCY_STATUS",
      wind_tendency_status_name(wind.kernel.status));
  write_double_attr(
      file,
      "TYWRF_WIND_TENDENCY_ACTIVE_U_POINTS",
      static_cast<double>(wind.kernel.active_u_points));
  write_double_attr(
      file,
      "TYWRF_WIND_TENDENCY_ACTIVE_V_POINTS",
      static_cast<double>(wind.kernel.active_v_points));
  write_double_attr(
      file,
      "TYWRF_WIND_TENDENCY_UPDATED_U_POINTS",
      static_cast<double>(wind.kernel.updated_u_points));
  write_double_attr(
      file,
      "TYWRF_WIND_TENDENCY_UPDATED_V_POINTS",
      static_cast<double>(wind.kernel.updated_v_points));
  write_double_attr(
      file,
      "TYWRF_WIND_TENDENCY_CHANGED_U_POINTS",
      static_cast<double>(wind.changed_u_points));
  write_double_attr(
      file,
      "TYWRF_WIND_TENDENCY_CHANGED_V_POINTS",
      static_cast<double>(wind.changed_v_points));
}

void stamp_gate_metadata(
    const Options& options,
    const Resolution resolution,
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const CandidateReport& report) {
  std::optional<PressureRefreshReadiness> pressure_refresh_readiness;
  if (report.pressure_refresh.has_value()) {
    pressure_refresh_readiness = evaluate_pressure_refresh_readiness(report);
  }
  const auto disposition =
      selected_field_candidate_disposition(options, report, pressure_refresh_readiness);
  const bool experimental_pressure_refresh_apply =
      disposition.experimental_pressure_refresh_apply;
  const NetcdfHandle file(options.output_path, NetcdfHandle::Mode::write);
  file.check(nc_redef(file.id()), "enter define mode");
  write_double_attr(file, "DX", resolution.dx);
  write_double_attr(file, "DY", resolution.dy);
  write_text_attr(
      file,
      "TYWRF_DIAGNOSTIC_ONLY",
      bool_text(disposition.diagnostic_only));
  write_text_attr(file, "TYWRF_GATE_CANDIDATE", bool_text(disposition.gate_candidate));
  write_text_attr(file, "TYWRF_INTEGRATOR_OUTPUT", bool_text(disposition.integrator_output));
  write_text_attr(file, "TYWRF_VALIDATION_GATE_ONLY", "false");
  write_text_attr(
      file,
      "TYWRF_CANDIDATE_KIND",
      disposition.candidate_kind);
  if (experimental_pressure_refresh_apply) {
    write_text_attr(file, "TYWRF_EXPERIMENTAL_PRESSURE_REFRESH_APPLY", "true");
  }
  write_text_attr(file, "TYWRF_CANDIDATE_DOMAIN", "d02");
  write_text_attr(file, "TYWRF_D02_RESOLUTION_CHECK", "d02_2km");
  write_text_attr(file, "TYWRF_CYCLE_START", options.cycle_start);
  write_text_attr(file, "TYWRF_CYCLE_END", options.cycle_end);
  write_text_attr(file, "TYWRF_D01_START_STATE_SOURCE", options.d01_start_state_path.string());
  write_text_attr(file, "TYWRF_D02_START_STATE_SOURCE", options.d02_start_state_path.string());
  write_text_attr(file, "TYWRF_TEMPLATE_SOURCE", options.template_path.string());
  write_text_attr(file, "TYWRF_STATE_VARIABLES", join_variables(options.variables));
  write_text_attr(
      file,
      "TYWRF_PARENT_INTERPOLATED_STATE_VARIABLES",
      join_variables(parent_interpolation_variable_names()));
  write_text_attr(file, "TYWRF_SELECTED_FIELD_TIMELINE_VERSION", "runtime_v0");
  write_text_attr(file, "TYWRF_SELECTED_FIELD_TIMELINE_EVIDENCE_ONLY", "true");
  write_int_attr(
      file,
      "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_COUNT",
      static_cast<std::int32_t>(report.timeline.size()));
  write_text_attr(
      file,
      "TYWRF_SELECTED_FIELD_TIMELINE_EVENT_NAMES",
      join_timeline_event_names(report.timeline));
  write_text_attr(
      file,
      "TYWRF_SELECTED_FIELD_TIMELINE_EVENTS",
      join_timeline_events(report.timeline));
  write_text_attr(file, "TYWRF_FROM_PARENT_START", std::to_string(options.from_parent_start.i_parent_start) + "," + std::to_string(options.from_parent_start.j_parent_start));
  write_text_attr(file, "TYWRF_TO_PARENT_START", std::to_string(options.to_parent_start.i_parent_start) + "," + std::to_string(options.to_parent_start.j_parent_start));
  write_int_attr(file, "I_PARENT_START", options.to_parent_start.i_parent_start);
  write_int_attr(file, "J_PARENT_START", options.to_parent_start.j_parent_start);
  write_double_attr(file, "CEN_LAT", report.cen_lat);
  write_double_attr(file, "CEN_LON", report.cen_lon);
  write_double_attr(file, "TYWRF_PARENT_GRID_RATIO", static_cast<double>(descriptor.parent_grid_ratio));
  write_double_attr(
      file,
      "TYWRF_SELECTED_FIELD_CHANGED_POINTS",
      static_cast<double>(report.changed_selected_points));
  write_text_attr(file, "TYWRF_STATIC_REFRESH_APPLIED", "true");
  write_text_attr(
      file,
      "TYWRF_STATIC_REFRESH_METHOD",
      "overlap_shift_xlat_xlong_extrapolate_hgt_parent_bilinear_v0");
  write_text_attr(file, "TYWRF_STATIC_REFRESH_USES_REFERENCE_END", "false");
  write_text_attr(file, "TYWRF_STATIC_REFRESH_D02_START_SOURCE", options.d02_start_state_path.string());
  write_text_attr(file, "TYWRF_STATIC_REFRESH_D01_HGT_SOURCE", options.d01_start_state_path.string());
  write_text_attr(file, "TYWRF_STATIC_REFRESH_TEMPLATE_SOURCE", options.template_path.string());
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_OVERLAP_CELLS",
      static_cast<double>(report.static_refresh.overlap_cell_count));
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_EXPOSED_CELLS",
      static_cast<double>(report.static_refresh.exposed_cell_count));
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_COORD_EXTRAPOLATED_CELLS",
      static_cast<double>(report.static_refresh.coordinate_extrapolated_cell_count));
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_HGT_PARENT_INTERPOLATED_CELLS",
      static_cast<double>(report.static_refresh.parent_hgt_interpolated_cell_count));
  write_double_attr(
      file,
      "TYWRF_STATIC_REFRESH_CHANGED_TEMPLATE_POINTS",
      static_cast<double>(report.changed_static_template_points));
  write_double_attr(
      file,
      "TYWRF_EXPOSED_EXCHANGE_POINTS",
      static_cast<double>(report.exchange.report.exchange_point_count));
  write_double_attr(
      file,
      "TYWRF_INTERPOLATED_POINTS",
      static_cast<double>(report.interpolation.interpolated_point_count));
  if (report.pressure_refresh.has_value()) {
    const auto& pressure = *report.pressure_refresh;
    const auto parity = pressure_refresh_report_parity(report, pressure);
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_OPT_IN", "true");
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_APPLIED", "true");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS",
        experimental_pressure_refresh_apply ? "experimental_apply_test_only"
                                            : "applied_to_candidate");
    if (experimental_pressure_refresh_apply) {
      write_text_attr(file, "TYWRF_PRESSURE_REFRESH_EXPERIMENTAL_APPLY", "true");
    }
    if (pressure_refresh_readiness.has_value()) {
      write_pressure_refresh_readiness_attrs(file, *pressure_refresh_readiness);
    }
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_PROVIDER_OK", "true");
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_SOURCE_SYNC_OK", "true");
    write_text_attr(file, "TYWRF_PRESSURE_REFRESH_COMPUTE_CALLED", "true");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TERRAIN_OVERRIDE_USED",
        pressure.provider_report.terrain_override_used ? "true" : "false");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE",
        pressure.provider_report.terrain_source_name);
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TERRAIN_PROVENANCE",
        pressure.provider_report.terrain_provenance);
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_METADATA_SOURCE",
        report.pressure_refresh_metadata_source.string());
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_METADATA_TIME_INDEX",
        static_cast<double>(report.pressure_refresh_metadata_time_index));
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_PRODUCTION_SOURCE",
        "krosa_moving_nest_pressure_refresh");
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_SYNCED_PB_POINTS",
        static_cast<double>(pressure.synced_pb_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_SYNCED_MUB_POINTS",
        static_cast<double>(pressure.synced_mub_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_SYNCED_PHB_POINTS",
        static_cast<double>(pressure.synced_phb_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT",
        static_cast<double>(pressure.compute_report.target_column_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT",
        static_cast<double>(pressure.compute_report.refreshed_column_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT",
        static_cast<double>(pressure.compute_report.refreshed_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_SKIPPED_POINT_COUNT",
        static_cast<double>(pressure.compute_report.skipped_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_INVALID_POINT_COUNT",
        static_cast<double>(pressure.compute_report.invalid_point_count));
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TOUCHED_OVERLAP_CELLS",
        pressure.compute_report.touched_overlap_cells ? "true" : "false");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_TOUCHED_HALO_CELLS",
        pressure.compute_report.touched_halo_cells ? "true" : "false");
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS",
        static_cast<double>(pressure.compute_report.refreshed_point_count));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS",
        static_cast<double>(report.pressure_refresh_changed_p_points));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_PB_POINTS",
        static_cast<double>(report.pressure_refresh_changed_pb_points));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_MUB_POINTS",
        static_cast<double>(report.pressure_refresh_changed_mub_points));
    write_double_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_PHB_POINTS",
        static_cast<double>(report.pressure_refresh_changed_phb_points));
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_MATCHES_REFRESHED_POINT_COUNT",
        parity.changed_p_matches_refreshed_point_count ? "true" : "false");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_INVALID_AND_SKIPPED_POINTS_ZERO",
        parity.invalid_and_skipped_points_zero ? "true" : "false");
    write_text_attr(
        file,
        "TYWRF_PRESSURE_REFRESH_OVERLAP_HALO_UNTOUCHED",
        parity.overlap_halo_untouched ? "true" : "false");
  }
  write_pressure_column_probe_attrs(file, options, report);
  write_pressure_formula_observation_attrs(file, report);
  write_normal_base_state_producer_attrs(file, report);
  write_diagnostic_adapter_attrs(file, report);
  write_wind_tendency_attrs(file, report);
  if (options.diagnostic_adapter_report) {
    write_text_attr(
        file,
        "TYWRF_CANDIDATE_MESSAGE",
        "Diagnostic-only selected-field base-state adapter staging report; "
        "not a validation gate pass or normal integrator output. The D70 adapter "
        "ran on staging buffers only and did not write selected-field candidate numerics.");
  } else if (experimental_pressure_refresh_apply) {
    write_text_attr(
        file,
        "TYWRF_CANDIDATE_MESSAGE",
        "Experimental selected-field pressure-refresh apply seam output for tool tests only; "
        "not a validation gate pass or normal integrator output. The apply path used "
        "moved_candidate_HGT terrain override and refreshed exposed P/PB/PHB/MUB.");
  } else if (report.pressure_refresh.has_value()) {
    write_text_attr(
        file,
        "TYWRF_CANDIDATE_MESSAGE",
        "Selected-field moving-nest candidate from start states only; U/V/T/PH/MU/QVAPOR "
        "exposed cells are parent interpolated, XLAT/XLONG/HGT are refreshed from "
        "start-state pose data, and pressure refresh applied to exposed P/PB/PHB/MUB.");
  } else {
    write_text_attr(
        file,
        "TYWRF_CANDIDATE_MESSAGE",
        "Selected-field moving-nest candidate from start states only; U/V/T/PH/MU/QVAPOR "
        "exposed cells are parent interpolated, XLAT/XLONG/HGT are refreshed from "
        "start-state pose data, exposed PB/PHB/MUB are updated by the normal "
        "non-oracle base-state producer, and P remains finite d02 start-state ownership.");
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

void write_pressure_refresh_readiness_json(
    std::ostream& stream,
    const PressureRefreshReadiness& readiness,
    const bool pretty) {
  write_json_bool(stream, "pressure_refresh_readiness_ready", readiness.ready(), true, pretty);
  write_json_bool(
      stream,
      "thermodynamic_base_state_consistency_ready",
      readiness.thermodynamic_base_state_consistency_ready,
      true,
      pretty);
  write_json_bool(
      stream,
      "provider_terrain_uses_moved_candidate_hgt",
      readiness.provider_terrain_uses_moved_candidate_hgt,
      true,
      pretty);
  write_json_bool(
      stream,
      "provider_base_state_reconstruct_ok",
      readiness.provider_base_state_reconstruct_ok,
      true,
      pretty);
  write_json_bool(
      stream, "base_state_sync_contract_ok", readiness.base_state_sync_contract_ok, true, pretty);
  write_json_bool(
      stream,
      "base_state_source_sync_readiness_check",
      readiness.base_state_sync_dry_run,
      true,
      pretty);
  write_json_bool(
      stream,
      "base_state_source_sync_applied",
      readiness.base_state_sync_applied,
      true,
      pretty);
  write_json_number(
      stream,
      "source_sync_planned_pb_point_count",
      static_cast<double>(readiness.would_sync_pb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "source_sync_planned_mub_point_count",
      static_cast<double>(readiness.would_sync_mub_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "source_sync_planned_phb_point_count",
      static_cast<double>(readiness.would_sync_phb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "source_sync_overlap_write_count",
      static_cast<double>(readiness.sync_overlap_write_count),
      true,
      pretty);
  write_json_number(
      stream,
      "source_sync_halo_write_count",
      static_cast<double>(readiness.sync_halo_write_count),
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_refresh_readiness_compute_called",
      readiness.pressure_refresh_compute_called,
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_compute_readiness_check",
      readiness.pressure_compute_dry_run,
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_compute_readiness_check_called",
      readiness.pressure_compute_dry_run_called,
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_compute_readiness_check_ok",
      readiness.pressure_compute_dry_run_ok,
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_refresh_planned_p_point_count",
      static_cast<double>(readiness.would_refresh_p_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_readiness_invalid_p_point_count",
      static_cast<double>(readiness.dry_run_invalid_p_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_readiness_skipped_p_point_count",
      static_cast<double>(readiness.dry_run_skipped_p_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_compute_readiness_report_target_column_count",
      static_cast<double>(readiness.pressure_compute_dry_run_report_target_column_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_compute_readiness_report_refreshed_point_count",
      static_cast<double>(readiness.pressure_compute_dry_run_report_refreshed_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_compute_readiness_report_invalid_point_count",
      static_cast<double>(readiness.pressure_compute_dry_run_report_invalid_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_compute_readiness_report_skipped_point_count",
      static_cast<double>(readiness.pressure_compute_dry_run_report_skipped_point_count),
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_compute_readiness_report_touched_overlap_cells",
      readiness.pressure_compute_dry_run_report_touched_overlap_cells,
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_compute_readiness_report_touched_halo_cells",
      readiness.pressure_compute_dry_run_report_touched_halo_cells,
      true,
      pretty);
  write_json_bool(
      stream,
      "pressure_refresh_readiness_applied",
      readiness.pressure_refresh_applied,
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_refresh_readiness_provider_terrain_source",
      readiness.provider_terrain_source_name,
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_refresh_readiness_provider_terrain_provenance",
      readiness.provider_terrain_provenance,
      true,
      pretty);
}

void write_pressure_column_probe_json(
    std::ostream& stream,
    const Options& options,
    const CandidateReport& report,
    const bool comma,
    const bool pretty) {
  const auto phases = pressure_column_probe_phase_names(report.pressure_column_observations);
  write_json_bool(stream, "pressure_column_probe_enabled", true, true, pretty);
  write_json_string(stream, "pressure_column_probe_version", "runtime_v0", true, pretty);
  write_json_bool(stream, "pressure_column_probe_evidence_only", true, true, pretty);
  write_json_string(
      stream, "pressure_column_probe_index_base", "zero_based_mass_grid", true, pretty);
  write_json_number(
      stream,
      "pressure_column_probe_column_count",
      static_cast<double>(options.pressure_column_probe_columns.size()),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_column_probe_level_count",
      static_cast<double>(options.pressure_column_probe_levels.size()),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_column_probe_phase_count",
      static_cast<double>(phases.size()),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_column_probe_record_count",
      static_cast<double>(report.pressure_column_observations.size()),
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_column_probe_columns",
      join_pressure_probe_columns(options.pressure_column_probe_columns),
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_column_probe_levels",
      join_pressure_probe_levels(options.pressure_column_probe_levels),
      true,
      pretty);
  write_json_string(stream, "pressure_column_probe_phases", join_variables(phases), true, pretty);
  write_json_string(
      stream,
      "pressure_column_probe_fields",
      join_variables(pressure_column_probe_field_names()),
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_column_probe_not_available",
      join_variables(pressure_column_probe_unavailable_names()),
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_column_probe_values",
      pressure_column_probe_values(report.pressure_column_observations),
      comma,
      pretty);
}

void write_normal_base_state_producer_json(
    std::ostream& stream,
    const CandidateReport& report,
    const bool comma,
    const bool pretty) {
  const auto& producer = *report.normal_base_state_producer;
  const auto* provider = report.normal_base_state_provider_source.has_value()
                             ? &*report.normal_base_state_provider_source
                             : nullptr;
  write_json_string(
      stream, "normal_base_state_producer_version", "a79_normal_v0", true, pretty);
  write_json_string(
      stream, "normal_base_state_producer_source", producer.source, true, pretty);
  write_json_string(
      stream,
      "normal_base_state_producer_disposition",
      producer.disposition,
      true,
      pretty);
  write_json_string(
      stream,
      "normal_base_state_producer_source_origin",
      producer.source_origin,
      true,
      pretty);
  write_json_bool(stream, "normal_base_state_producer_ok", producer.ok(), true, pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_diagnostic_only",
      producer.diagnostic_only,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_normal_candidate_producer",
      producer.normal_candidate_producer,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_gate_candidate",
      producer.gate_candidate,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_integrator_output",
      producer.integrator_output,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_writes_candidate",
      producer.writes_candidate,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_writes_netcdf",
      producer.writes_netcdf,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_uses_reference_end_truth",
      producer.uses_reference_end_truth,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_no_reference_end_truth",
      !producer.uses_reference_end_truth,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_uses_direct_p_shortcut",
      producer.uses_direct_p_shortcut,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_no_direct_p_shortcut",
      !producer.uses_direct_p_shortcut,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_reads_direct_p",
      producer.reads_direct_p,
      true,
      pretty);
  write_json_bool(
      stream, "normal_base_state_producer_writes_p", producer.writes_p, true, pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_no_gate_pass_claim",
      producer.no_gate_pass_claim,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_no_00_20_progression",
      producer.no_00_20_progression,
      true,
      pretty);
  write_json_string(
      stream,
      "normal_base_state_producer_written_fields",
      producer.written_fields,
      true,
      pretty);
  write_json_string(
      stream,
      "normal_base_state_producer_metadata_source",
      report.normal_base_state_metadata_source.string(),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_metadata_time_index",
      static_cast<double>(report.normal_base_state_metadata_time_index),
      true,
      pretty);
  write_json_string(
      stream,
      "normal_base_state_producer_provider_source",
      provider == nullptr ? "" : provider->source,
      true,
      pretty);
  write_json_string(
      stream,
      "normal_base_state_producer_terrain_source",
      provider == nullptr ? "" : provider->terrain_source_name,
      true,
      pretty);
  write_json_string(
      stream,
      "normal_base_state_producer_terrain_provenance",
      provider == nullptr ? "" : provider->terrain_provenance,
      true,
      pretty);
  write_json_bool(
      stream,
      "normal_base_state_producer_terrain_override_used",
      provider != nullptr && provider->terrain_override_used,
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_exposed_mass_cell_count",
      static_cast<double>(producer.exposed_mass_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_exposed_surface_cell_count",
      static_cast<double>(producer.exposed_surface_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_exposed_w_full_column_count",
      static_cast<double>(producer.exposed_w_full_column_count),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_pb_written_point_count",
      static_cast<double>(producer.pb_written_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_mub_written_cell_count",
      static_cast<double>(producer.mub_written_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_phb_written_point_count",
      static_cast<double>(producer.phb_written_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_ht_written_cell_count",
      static_cast<double>(producer.ht_written_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_direct_write_point_count",
      static_cast<double>(producer.direct_write_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_changed_p_points",
      static_cast<double>(report.normal_base_state_changed_p_points),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_changed_pb_points",
      static_cast<double>(report.normal_base_state_changed_pb_points),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_changed_mub_points",
      static_cast<double>(report.normal_base_state_changed_mub_points),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_changed_phb_points",
      static_cast<double>(report.normal_base_state_changed_phb_points),
      true,
      pretty);
  write_json_number(
      stream,
      "normal_base_state_producer_changed_hgt_points",
      static_cast<double>(report.normal_base_state_changed_hgt_points),
      comma,
      pretty);
}

void write_wind_tendency_json(
    std::ostream& stream,
    const CandidateReport& report,
    const bool comma,
    const bool pretty) {
  if (!report.wind_tendency.has_value()) {
    return;
  }

  const auto& wind = *report.wind_tendency;
  const bool gate_evidence = wind_tendency_gate_evidence(wind.source_kind);
  const bool zero_or_identity_only = wind_tendency_zero_or_identity_only(wind.source_kind);
  write_json_bool(stream, "wind_tendency_opt_in", true, true, pretty);
  write_json_bool(stream, "wind_tendency_applied", true, true, pretty);
  write_json_string(
      stream,
      "wind_tendency_source_kind",
      wind_tendency_source_name(wind.source_kind),
      true,
      pretty);
  write_json_bool(stream, "wind_tendency_gate_evidence", gate_evidence, true, pretty);
  write_json_bool(
      stream, "wind_tendency_validation_gate_evidence", gate_evidence, true, pretty);
  write_json_bool(stream, "wind_tendency_uses_reference_end_truth", false, true, pretty);
  write_json_bool(
      stream, "wind_tendency_zero_or_identity_only", zero_or_identity_only, true, pretty);
  write_json_string(stream, "wind_tendency_written_fields", "U,V", true, pretty);
  write_json_string(
      stream,
      "wind_tendency_status",
      wind_tendency_status_name(wind.kernel.status),
      true,
      pretty);
  write_json_number(
      stream,
      "wind_tendency_active_u_points",
      static_cast<double>(wind.kernel.active_u_points),
      true,
      pretty);
  write_json_number(
      stream,
      "wind_tendency_active_v_points",
      static_cast<double>(wind.kernel.active_v_points),
      true,
      pretty);
  write_json_number(
      stream,
      "wind_tendency_updated_u_points",
      static_cast<double>(wind.kernel.updated_u_points),
      true,
      pretty);
  write_json_number(
      stream,
      "wind_tendency_updated_v_points",
      static_cast<double>(wind.kernel.updated_v_points),
      true,
      pretty);
  write_json_number(
      stream,
      "wind_tendency_changed_u_points",
      static_cast<double>(wind.changed_u_points),
      true,
      pretty);
  write_json_number(
      stream,
      "wind_tendency_changed_v_points",
      static_cast<double>(wind.changed_v_points),
      comma,
      pretty);
}

void write_pressure_formula_observation_json(
    std::ostream& stream,
    const CandidateReport& report,
    const bool comma,
    const bool pretty) {
  const auto& compute = report.pressure_refresh->compute_report;
  write_json_bool(stream, "pressure_formula_observation_enabled", true, true, pretty);
  write_json_string(
      stream, "pressure_formula_observation_version", "runtime_v0", true, pretty);
  write_json_bool(
      stream, "pressure_formula_observation_evidence_only", true, true, pretty);
  write_json_string(
      stream,
      "pressure_formula_observation_index_base",
      "zero_based_mass_grid",
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_formula_observation_request_count",
      static_cast<double>(compute.observation_request_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_formula_observation_record_count",
      static_cast<double>(compute.observation_record_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_formula_observation_valid_count",
      static_cast<double>(compute.observation_valid_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_formula_observation_invalid_count",
      static_cast<double>(compute.observation_invalid_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_formula_observation_out_of_bounds_count",
      static_cast<double>(compute.observation_out_of_bounds_count),
      true,
      pretty);
  write_json_number(
      stream,
      "pressure_formula_observation_outside_target_region_count",
      static_cast<double>(compute.observation_outside_target_region_count),
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_formula_observation_fields",
      join_variables(pressure_formula_observation_field_names()),
      true,
      pretty);
  write_json_string(
      stream,
      "pressure_formula_observation_values",
      pressure_formula_observation_values(report.pressure_formula_observations),
      comma,
      pretty);
}

void write_diagnostic_adapter_provider_source_json(
    std::ostream& stream,
    const CandidateReport& report,
    const bool pretty) {
  if (!report.diagnostic_adapter_provider_source.has_value()) {
    return;
  }

  const auto& provider = *report.diagnostic_adapter_provider_source;
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_version",
      "d77_provider_source_v0",
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_origin",
      diagnostic_adapter_provider_source_origin(),
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_source_origin",
      diagnostic_adapter_provider_source_origin(),
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_provider_source",
      provider.source,
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_terrain_source",
      provider.terrain_source_name,
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_terrain_provenance",
      provider.terrain_provenance,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_terrain_override_used",
      provider.terrain_override_used,
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_ht_source",
      "output_static.hgt",
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_ht_provenance",
      "adapter_HT_from_output_static_HGT",
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_provider_ok",
      provider.ok(),
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_provider_source_result_message",
      provider.result.message == nullptr ? "" : provider.result.message,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_diagnostic_only",
      provider.diagnostic_only,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_gate_candidate",
      provider.gate_candidate,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_integrator_output",
      provider.integrator_output,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_writes_candidate",
      false,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_writes_netcdf",
      false,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_no_candidate_write",
      true,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_uses_reference_end_truth",
      false,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_no_reference_end_truth",
      true,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_uses_direct_p_shortcut",
      false,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_no_direct_p_shortcut",
      true,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_reads_direct_p",
      false,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_wrote_pb",
      provider.wrote_pb,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_wrote_t_init",
      provider.wrote_t_init,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_wrote_mub",
      provider.wrote_mub,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_wrote_alb",
      provider.wrote_alb,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_wrote_phb",
      provider.wrote_phb,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_provider_source_provider_reconstructed_phb_not_wrf_rebalance_validated",
      true,
      true,
      pretty);
}

void write_diagnostic_adapter_source_staging_json(
    std::ostream& stream,
    const CandidateReport& report,
    const bool pretty) {
  if (!report.diagnostic_adapter_source_staging.has_value()) {
    return;
  }

  const auto& staging = *report.diagnostic_adapter_source_staging;
  write_json_string(
      stream,
      "diagnostic_adapter_source_staging_version",
      "d75_provider_report_v0",
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_source_staging_provider_kind",
      "BaseStateSourceStagingProvider",
      true,
      pretty);
  write_json_string(
      stream, "diagnostic_adapter_source_staging_source", staging.source, true, pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_source_staging_disposition",
      staging.disposition,
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_source_staging_source_shape",
      staging.source_shape,
      true,
      pretty);
  write_json_bool(stream, "diagnostic_adapter_source_staging_ok", staging.ok(), true, pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_source_staging_result_message",
      staging.result.message == nullptr ? "" : staging.result.message,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_diagnostic_only",
      staging.diagnostic_only,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_gate_candidate",
      staging.gate_candidate,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_integrator_output",
      staging.integrator_output,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_writes_candidate",
      staging.writes_candidate,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_writes_netcdf",
      staging.writes_netcdf,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_candidate_buffers_preserved",
      staging.candidate_buffers_preserved,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_owns_staging_buffers",
      staging.owns_staging_buffers,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_allocated_buffers",
      staging.allocated_buffers,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_uses_reference_end_truth",
      staging.uses_reference_end_truth,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_uses_direct_p_shortcut",
      staging.uses_direct_p_shortcut,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_reads_direct_p",
      staging.reads_direct_p,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_staging_aliases_child",
      report.diagnostic_adapter_source_staging_aliases_child,
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_active_nx",
      static_cast<double>(staging.active_nx),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_active_ny",
      static_cast<double>(staging.active_ny),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_mass_nz",
      static_cast<double>(staging.mass_nz),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_full_nz",
      static_cast<double>(staging.full_nz),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_exposed_region_count",
      static_cast<double>(staging.exposed_region_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_exposed_mass_cell_count",
      static_cast<double>(staging.exposed_mass_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_exposed_mass_point_count",
      static_cast<double>(staging.exposed_mass_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_exposed_surface_cell_count",
      static_cast<double>(staging.exposed_surface_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_exposed_w_full_column_count",
      static_cast<double>(staging.exposed_w_full_column_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_exposed_w_full_point_count",
      static_cast<double>(staging.exposed_w_full_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_masked_mass_cell_count",
      static_cast<double>(staging.masked_mass_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_masked_mass_point_count",
      static_cast<double>(staging.masked_mass_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_masked_surface_cell_count",
      static_cast<double>(staging.masked_surface_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_masked_w_full_column_count",
      static_cast<double>(staging.masked_w_full_column_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_masked_w_full_point_count",
      static_cast<double>(staging.masked_w_full_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_staged_phb_point_count",
      static_cast<double>(staging.staged_phb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_staged_mub_cell_count",
      static_cast<double>(staging.staged_mub_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_staged_ht_cell_count",
      static_cast<double>(staging.staged_ht_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_staged_pb_point_count",
      static_cast<double>(staging.staged_pb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_staged_t_init_point_count",
      static_cast<double>(staging.staged_t_init_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_staged_alb_point_count",
      static_cast<double>(staging.staged_alb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_staged_value_count",
      static_cast<double>(staging.staged_value_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_active_masked_value_count",
      static_cast<double>(staging.active_masked_value_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_halo_masked_value_count",
      static_cast<double>(staging.halo_masked_value_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_invalid_exposed_value_count",
      static_cast<double>(staging.invalid_exposed_value_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_invalid_exposed_phb_point_count",
      static_cast<double>(staging.invalid_exposed_phb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_invalid_exposed_mub_cell_count",
      static_cast<double>(staging.invalid_exposed_mub_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_invalid_exposed_ht_cell_count",
      static_cast<double>(staging.invalid_exposed_ht_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_invalid_exposed_pb_point_count",
      static_cast<double>(staging.invalid_exposed_pb_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_invalid_exposed_t_init_point_count",
      static_cast<double>(staging.invalid_exposed_t_init_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_staging_invalid_exposed_alb_point_count",
      static_cast<double>(staging.invalid_exposed_alb_point_count),
      true,
      pretty);
}

void write_diagnostic_adapter_source_child_delta_field_json(
    std::ostream& stream,
    const std::string_view field,
    const DiagnosticAdapterSourceChildFieldDeltaReport& delta,
    const bool pretty) {
  const std::string prefix =
      "diagnostic_adapter_source_child_delta_" + std::string(field);
  write_json_number(
      stream,
      prefix + "_compared_value_count",
      static_cast<double>(delta.compared_value_count),
      true,
      pretty);
  write_json_number(
      stream,
      prefix + "_differing_value_count",
      static_cast<double>(delta.differing_value_count),
      true,
      pretty);
  write_json_number(
      stream, prefix + "_max_abs_diff", delta.max_abs_diff, true, pretty);
}

void write_diagnostic_adapter_source_child_delta_json(
    std::ostream& stream,
    const CandidateReport& report,
    const bool pretty) {
  if (!report.diagnostic_adapter_source_child_delta.has_value()) {
    return;
  }

  const auto& delta = *report.diagnostic_adapter_source_child_delta;
  write_json_string(
      stream,
      "diagnostic_adapter_source_child_delta_version",
      "a76_source_child_delta_v0",
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_source_child_delta_source",
      "BaseStateSourceStagingProvider_vs_child_staging_pre_adapter",
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_source_child_delta_scope",
      "exposed_base_state_values_only",
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_source_child_delta_fields",
      "PHB,MUB,HT,PB,T_INIT,ALB",
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_child_delta_diagnostic_only",
      delta.diagnostic_only,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_child_delta_gate_candidate",
      delta.gate_candidate,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_child_delta_integrator_output",
      delta.integrator_output,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_child_delta_writes_candidate",
      delta.writes_candidate,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_child_delta_writes_netcdf",
      delta.writes_netcdf,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_source_child_delta_values_identical",
      delta.values_identical,
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_child_delta_compared_value_count",
      static_cast<double>(delta.compared_value_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_child_delta_differing_value_count",
      static_cast<double>(delta.differing_value_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_source_child_delta_max_abs_diff",
      delta.max_abs_diff,
      true,
      pretty);
  write_diagnostic_adapter_source_child_delta_field_json(
      stream, "phb", delta.phb, pretty);
  write_diagnostic_adapter_source_child_delta_field_json(
      stream, "mub", delta.mub, pretty);
  write_diagnostic_adapter_source_child_delta_field_json(
      stream, "ht", delta.ht, pretty);
  write_diagnostic_adapter_source_child_delta_field_json(
      stream, "pb", delta.pb, pretty);
  write_diagnostic_adapter_source_child_delta_field_json(
      stream, "t_init", delta.t_init, pretty);
  write_diagnostic_adapter_source_child_delta_field_json(
      stream, "alb", delta.alb, pretty);
}

void write_diagnostic_adapter_json(
    std::ostream& stream,
    const CandidateReport& report,
    const bool comma,
    const bool pretty) {
  const auto& adapter = *report.diagnostic_adapter;
  write_json_bool(stream, "diagnostic_adapter_opt_in", true, true, pretty);
  write_json_string(stream, "diagnostic_adapter_version", "d70_report_v0", true, pretty);
  write_json_string(stream, "diagnostic_adapter_source", adapter.source, true, pretty);
  write_json_string(stream, "diagnostic_adapter_disposition", adapter.disposition, true, pretty);
  write_json_string(
      stream, "diagnostic_adapter_exchange_source", adapter.exchange_source, true, pretty);
  write_json_string(
      stream, "diagnostic_adapter_recompute_source", adapter.recompute_source, true, pretty);
  write_json_bool(stream, "diagnostic_adapter_ok", adapter.ok(), true, pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_result_message",
      adapter.result.message == nullptr ? "" : adapter.result.message,
      true,
      pretty);
  write_json_bool(
      stream, "diagnostic_adapter_diagnostic_only", adapter.diagnostic_only, true, pretty);
  write_json_bool(
      stream, "diagnostic_adapter_gate_candidate", adapter.gate_candidate, true, pretty);
  write_json_bool(
      stream, "diagnostic_adapter_integrator_output", adapter.integrator_output, true, pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_selected_field_numerics_enabled",
      adapter.selected_field_numerics_enabled,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_enables_selected_field_numerics",
      adapter.enables_selected_field_numerics,
      true,
      pretty);
  write_json_bool(
      stream, "diagnostic_adapter_writes_netcdf", adapter.writes_netcdf, true, pretty);
  write_json_bool(
      stream, "diagnostic_adapter_writes_candidate", adapter.writes_candidate, true, pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_called_d68_exchange",
      adapter.called_d68_exchange,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_called_d69_recompute",
      adapter.called_d69_recompute,
      true,
      pretty);
  write_json_bool(
      stream, "diagnostic_adapter_ht_is_hgt_alias", adapter.ht_is_hgt_alias, true, pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_creates_terrain_owner",
      adapter.creates_terrain_owner,
      true,
      pretty);
  write_json_bool(
      stream,
      "diagnostic_adapter_terrain_owner_created",
      adapter.terrain_owner_created,
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_metadata_source",
      report.diagnostic_adapter_metadata_source.string(),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_metadata_time_index",
      static_cast<double>(report.diagnostic_adapter_metadata_time_index),
      true,
      pretty);
  write_json_number(
      stream, "diagnostic_adapter_active_nx", static_cast<double>(adapter.active_nx), true, pretty);
  write_json_number(
      stream, "diagnostic_adapter_active_ny", static_cast<double>(adapter.active_ny), true, pretty);
  write_json_number(
      stream, "diagnostic_adapter_active_nz", static_cast<double>(adapter.active_nz), true, pretty);
  write_json_number(
      stream, "diagnostic_adapter_c3h_count", static_cast<double>(adapter.c3h_count), true, pretty);
  write_json_number(
      stream, "diagnostic_adapter_c4h_count", static_cast<double>(adapter.c4h_count), true, pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_exposed_region_count",
      static_cast<double>(adapter.exposed_region_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_exposed_mass_cell_count",
      static_cast<double>(adapter.exposed_mass_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_exposed_surface_cell_count",
      static_cast<double>(adapter.exposed_surface_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_exposed_w_full_column_count",
      static_cast<double>(adapter.exposed_w_full_column_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_phb_written_point_count",
      static_cast<double>(adapter.phb_written_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_mub_written_cell_count",
      static_cast<double>(adapter.mub_written_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_ht_written_cell_count",
      static_cast<double>(adapter.ht_written_cell_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_direct_write_point_count",
      static_cast<double>(adapter.direct_write_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_exchange_recompute_mark_count",
      static_cast<double>(adapter.exchange_recompute_mark_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_recomputed_point_count",
      static_cast<double>(adapter.recomputed_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_pb_recomputed_point_count",
      static_cast<double>(adapter.pb_recomputed_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_t_init_recomputed_point_count",
      static_cast<double>(adapter.t_init_recomputed_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_alb_recomputed_point_count",
      static_cast<double>(adapter.alb_recomputed_point_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_invalid_column_count",
      static_cast<double>(adapter.invalid_column_count),
      true,
      pretty);
  write_json_number(
      stream,
      "diagnostic_adapter_invalid_point_count",
      static_cast<double>(adapter.invalid_point_count),
      true,
      pretty);
  write_diagnostic_adapter_provider_source_json(stream, report, pretty);
  write_diagnostic_adapter_source_staging_json(stream, report, pretty);
  write_diagnostic_adapter_source_child_delta_json(stream, report, pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_written_fields",
      "PHB,MUB,HGT,PB,T_INIT,ALB",
      true,
      pretty);
  write_json_string(
      stream,
      "diagnostic_adapter_integration_status",
      "staging_report_only_no_gate_no_integrator",
      comma,
      pretty);
}

void print_report(
    const Options& options,
    const Resolution resolution,
    const tywrf::nest::ParentChildDescriptor& descriptor,
    const CandidateReport& report) {
  const bool pretty = options.pretty;
  std::optional<PressureRefreshReadiness> pressure_refresh_readiness;
  if (report.pressure_refresh.has_value()) {
    pressure_refresh_readiness = evaluate_pressure_refresh_readiness(report);
  }
  const auto disposition =
      selected_field_candidate_disposition(options, report, pressure_refresh_readiness);
  const bool experimental_pressure_refresh_apply =
      disposition.experimental_pressure_refresh_apply;
  const bool has_normal_base_state_producer =
      report.normal_base_state_producer.has_value();
  const bool has_diagnostic_adapter = report.diagnostic_adapter.has_value();
  const bool has_wind_tendency = report.wind_tendency.has_value();
  const bool has_pressure_column_probe = !report.pressure_column_observations.empty();
  const bool has_pressure_formula_observation =
      !report.pressure_formula_observations.empty();
  std::cout << "{" << (pretty ? "\n" : "");
  write_json_string(
      std::cout,
      "status",
      disposition.status,
      true,
      pretty);
  write_json_string(
      std::cout,
      "candidate_kind",
      disposition.candidate_kind,
      true,
      pretty);
  write_json_bool(std::cout, "diagnostic_only", disposition.diagnostic_only, true, pretty);
  write_json_bool(std::cout, "gate_candidate", disposition.gate_candidate, true, pretty);
  write_json_bool(std::cout, "integrator_output", disposition.integrator_output, true, pretty);
  write_json_bool(std::cout, "validation_gate_only", false, true, pretty);
  if (experimental_pressure_refresh_apply) {
    write_json_bool(std::cout, "experimental_pressure_refresh_apply", true, true, pretty);
  }
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
      "static_refresh_changed_template_points",
      static_cast<double>(report.changed_static_template_points),
      true,
      pretty);
  write_json_bool(std::cout, "static_refresh_uses_reference_end", false, true, pretty);
  write_json_number(
      std::cout,
      "static_refresh_exposed_cells",
      static_cast<double>(report.static_refresh.exposed_cell_count),
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
      true,
      pretty);
  write_json_string(
      std::cout,
      "parent_interpolated_state_variables",
      join_variables(parent_interpolation_variable_names()),
      true,
      pretty);
  write_json_number(
      std::cout,
      "selected_field_timeline_event_count",
      static_cast<double>(report.timeline.size()),
      true,
      pretty);
  write_json_string(
      std::cout,
      "selected_field_timeline_event_names",
      join_timeline_event_names(report.timeline),
      true,
      pretty);
  write_json_string(
      std::cout,
      "selected_field_timeline_events",
      join_timeline_events(report.timeline),
      has_normal_base_state_producer || has_diagnostic_adapter ||
          has_wind_tendency || report.pressure_refresh.has_value() ||
          has_pressure_column_probe,
      pretty);
  if (has_wind_tendency) {
    write_wind_tendency_json(std::cout, report, true, pretty);
  }
  if (has_normal_base_state_producer) {
    write_normal_base_state_producer_json(
        std::cout,
        report,
        has_diagnostic_adapter || report.pressure_refresh.has_value() ||
            has_pressure_column_probe || has_pressure_formula_observation,
        pretty);
  }
  if (report.pressure_refresh.has_value()) {
    const auto& pressure = *report.pressure_refresh;
    const auto parity = pressure_refresh_report_parity(report, pressure);
    write_json_bool(std::cout, "pressure_refresh_opt_in", true, true, pretty);
    write_json_bool(std::cout, "pressure_refresh_applied", true, true, pretty);
    write_json_bool(std::cout, "pressure_refresh_provider_ok", true, true, pretty);
    write_json_bool(std::cout, "pressure_refresh_source_sync_ok", true, true, pretty);
    write_json_bool(std::cout, "pressure_refresh_compute_called", true, true, pretty);
    if (experimental_pressure_refresh_apply) {
      write_json_bool(
          std::cout,
          "pressure_refresh_experimental_apply",
          true,
          true,
          pretty);
    }
    write_json_string(
        std::cout,
        "pressure_refresh_integration_status",
        experimental_pressure_refresh_apply ? "experimental_apply_test_only"
                                            : "applied_to_candidate",
        true,
        pretty);
    if (pressure_refresh_readiness.has_value()) {
      write_pressure_refresh_readiness_json(std::cout, *pressure_refresh_readiness, pretty);
    }
    write_json_bool(
        std::cout,
        "pressure_refresh_terrain_override_used",
        pressure.provider_report.terrain_override_used,
        true,
        pretty);
    write_json_string(
        std::cout,
        "pressure_refresh_terrain_source",
        pressure.provider_report.terrain_source_name,
        true,
        pretty);
    write_json_string(
        std::cout,
        "pressure_refresh_terrain_provenance",
        pressure.provider_report.terrain_provenance,
        true,
        pretty);
    write_json_string(
        std::cout,
        "pressure_refresh_metadata_source",
        report.pressure_refresh_metadata_source.string(),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_metadata_time_index",
        static_cast<double>(report.pressure_refresh_metadata_time_index),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_target_column_count",
        static_cast<double>(pressure.compute_report.target_column_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_refreshed_column_count",
        static_cast<double>(pressure.compute_report.refreshed_column_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_refreshed_point_count",
        static_cast<double>(pressure.compute_report.refreshed_point_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_skipped_point_count",
        static_cast<double>(pressure.compute_report.skipped_point_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_invalid_point_count",
        static_cast<double>(pressure.compute_report.invalid_point_count),
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_touched_overlap_cells",
        pressure.compute_report.touched_overlap_cells,
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_touched_halo_cells",
        pressure.compute_report.touched_halo_cells,
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_refreshed_p_points",
        static_cast<double>(pressure.compute_report.refreshed_point_count),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_changed_p_points",
        static_cast<double>(report.pressure_refresh_changed_p_points),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_changed_pb_points",
        static_cast<double>(report.pressure_refresh_changed_pb_points),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_changed_mub_points",
        static_cast<double>(report.pressure_refresh_changed_mub_points),
        true,
        pretty);
    write_json_number(
        std::cout,
        "pressure_refresh_changed_phb_points",
        static_cast<double>(report.pressure_refresh_changed_phb_points),
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_changed_p_matches_refreshed_point_count",
        parity.changed_p_matches_refreshed_point_count,
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_invalid_and_skipped_points_zero",
        parity.invalid_and_skipped_points_zero,
        true,
        pretty);
    write_json_bool(
        std::cout,
        "pressure_refresh_overlap_halo_untouched",
        parity.overlap_halo_untouched,
        has_pressure_column_probe || has_pressure_formula_observation,
        pretty);
  }
  if (has_diagnostic_adapter) {
    write_diagnostic_adapter_json(
        std::cout,
        report,
        has_pressure_formula_observation || has_pressure_column_probe,
        pretty);
  }
  if (has_pressure_formula_observation) {
    write_pressure_formula_observation_json(
        std::cout, report, has_pressure_column_probe, pretty);
  }
  if (has_pressure_column_probe) {
    write_pressure_column_probe_json(std::cout, options, report, false, pretty);
  }
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
  finalize_pressure_column_probe_options(options, d02_grid);

  const auto d01_start_hgt =
      load_hgt_field(options.d01_start_state_path, d01_grid, options.d01_time_index);
  const auto d02_start_static =
      load_static_fields(options.d02_start_state_path, d02_grid, options.d02_time_index);
  const auto template_static =
      load_static_fields(options.template_path, d02_grid, options.template_time_index);
  require_finite_storage("d01 HGT", d01_start_hgt);
  require_finite_static_fields(d02_start_static);
  require_finite_static_fields(template_static);

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
  auto report = build_candidate_state(descriptor, options, d01_start, d02_start, candidate);
  StaticFieldSet output_static(d02_grid);
  refresh_static_fields(
      descriptor,
      options,
      d02_start_static,
      template_static,
      d01_start_hgt,
      output_static,
      report);
  if (!options.diagnostic_adapter_report) {
    run_normal_base_state_producer(options, candidate, output_static, report);
  }
  capture_pressure_column_observations(
      options, "post_static_refresh", candidate, output_static, report);
  if (options.diagnostic_adapter_report) {
    run_diagnostic_adapter_report(options, candidate, output_static, report);
  }
  if (options.pressure_refresh) {
    probe_pressure_refresh_provider_readiness(options, candidate, output_static, report);
    probe_pressure_refresh_dry_run_contract(options, candidate, output_static, report);
    if (!options.experimental_pressure_refresh_apply) {
      require_pressure_refresh_ready_for_compute(report);
    }
    apply_pressure_refresh(options, candidate, output_static, report);
    capture_pressure_column_observations(
        options, "post_pressure_refresh", candidate, output_static, report);
  } else {
    append_timeline_event(
        report,
        "pressure_refresh_readiness",
        {
            timeline_field("opt_in", "false"),
            timeline_field("ready", "not_applicable"),
            timeline_field("status", "skipped"),
        });
    append_timeline_event(
        report,
        "pressure_refresh_apply",
        {
            timeline_field("opt_in", "false"),
            timeline_field("applied", "false"),
            timeline_field("status", "skipped"),
        });
    capture_pressure_column_observations(
        options, "pressure_refresh_skipped", candidate, output_static, report);
  }
  append_pressure_column_probe_timeline(report);
  append_timeline_event(
      report,
      "cycle_end",
      {
          timeline_field("cycle_start", options.cycle_start),
          timeline_field("cycle_end", options.cycle_end),
      });
  append_timeline_event(
      report,
      "output_write_preparation",
      {
          timeline_field("time_index", static_cast<std::uint64_t>(options.output_time_index)),
          timeline_field("variable_count", static_cast<std::uint64_t>(options.variables.size())),
          timeline_field("times", options.times_value),
          timeline_field("state_write", "pending"),
          timeline_field("metadata_write", "pending"),
      });

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
  overwrite_output_static_fields(options.output_path, output_static, options.output_time_index);
  stamp_gate_metadata(options, template_resolution, descriptor, report);
  print_report(options, template_resolution, descriptor, report);
  return 0;
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    if (argc == 2 && std::string_view(argv[1]) == "--candidate-disposition-self-test") {
      return run_candidate_disposition_self_test();
    }
    const auto options = parse_options(argc, argv);
    return run(options);
  } catch (const std::exception& error) {
    std::cerr << "tywrf_selected_field_cycle: " << error.what() << "\n\n" << usage();
    return 2;
  }
}
