#include "tywrf/dynamics/base_state_provider.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace tywrf::dynamics {
namespace {

struct ExposedRegionSet {
  std::array<nest::RemapWindow, 4> regions{};
  std::uint8_t count = 0;
  std::uint64_t horizontal_cell_count = 0;
};

[[nodiscard]] constexpr nest::NestResult ok_result() noexcept {
  return {nest::NestStatus::ok, "ok"};
}

[[nodiscard]] constexpr nest::NestResult result(
    const nest::NestStatus status,
    const char* message) noexcept {
  return {status, message};
}

[[nodiscard]] constexpr std::int32_t vector_size_or_negative(
    const std::size_t size) noexcept {
  if (size > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    return -1;
  }
  return static_cast<std::int32_t>(size);
}

[[nodiscard]] constexpr std::size_t point_count_2d(
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
}

[[nodiscard]] constexpr std::uint64_t point_count_2d_u64(
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny);
}

[[nodiscard]] constexpr bool valid_grid_config(const GridConfig& config) noexcept {
  return config.mass_nx > 0 && config.mass_ny > 0 && config.mass_nz > 0 &&
         config.full_nz == config.mass_nz + 1 && config.halo.i_lower >= 0 &&
         config.halo.i_upper >= 0 && config.halo.j_lower >= 0 &&
         config.halo.j_upper >= 0 && config.halo.k_lower >= 0 &&
         config.halo.k_upper >= 0;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(
    const FieldView2D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(
    const FieldView2D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(
    const FieldView3D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(
    const FieldView3D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nz(
    const FieldView3D<Real>& field) noexcept {
  return field.nz - field.halo.k_lower - field.halo.k_upper;
}

template <typename Real>
[[nodiscard]] constexpr bool valid_canonical_terrain_view(
    const FieldView2D<Real>& field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 &&
         field.stride_i == 1 && field.stride_j == field.nx &&
         field.halo.i_lower >= 0 && field.halo.i_upper >= 0 &&
         field.halo.j_lower >= 0 && field.halo.j_upper >= 0 &&
         active_nx(field) > 0 && active_ny(field) > 0;
}

template <typename Real>
[[nodiscard]] constexpr bool valid_canonical_view(
    const FieldView2D<Real>& field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 &&
         field.stride_i == 1 && field.stride_j == field.nx &&
         field.halo.i_lower >= 0 && field.halo.i_upper >= 0 &&
         field.halo.j_lower >= 0 && field.halo.j_upper >= 0 &&
         active_nx(field) > 0 && active_ny(field) > 0;
}

template <typename Real>
[[nodiscard]] constexpr bool valid_canonical_view(
    const FieldView3D<Real>& field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 &&
         field.nz > 0 && field.stride_i == 1 && field.stride_k == field.nx &&
         field.stride_j == field.nx * field.nz && field.halo.i_lower >= 0 &&
         field.halo.i_upper >= 0 && field.halo.j_lower >= 0 &&
         field.halo.j_upper >= 0 && field.halo.k_lower >= 0 &&
         field.halo.k_upper >= 0 && active_nx(field) > 0 &&
         active_ny(field) > 0 && active_nz(field) > 0;
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_active_shape(
    const FieldView3D<LhsReal>& lhs,
    const FieldView3D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs) &&
         active_nz(lhs) == active_nz(rhs);
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_horizontal_active_shape(
    const FieldView2D<LhsReal>& lhs,
    const FieldView3D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

[[nodiscard]] constexpr bool sane_terrain_height_m(const float value) noexcept {
  return value >= -1000.0F && value <= 20000.0F;
}

[[nodiscard]] bool valid_options(
    const KrosaMassBaseStateReconstructionOptions& options) noexcept {
  return std::isfinite(options.dry_air_gas_constant) &&
         std::isfinite(options.gravity) &&
         std::isfinite(options.sea_level_base_pressure_pa) &&
         std::isfinite(options.sea_level_base_temperature_k) &&
         std::isfinite(options.base_lapse_k) &&
         std::isfinite(options.isothermal_temperature_k) &&
         std::isfinite(options.stratosphere_base_pressure_pa) &&
         std::isfinite(options.stratosphere_lapse_k) &&
         std::isfinite(options.reference_pressure_pa) &&
         std::isfinite(options.base_potential_temperature_k) &&
         std::isfinite(options.specific_heat_cp) && std::isfinite(options.cvpm) &&
         options.dry_air_gas_constant > 0.0F && options.gravity > 0.0F &&
         options.sea_level_base_pressure_pa > 0.0F &&
         options.base_lapse_k > 0.0F && options.reference_pressure_pa > 0.0F &&
         options.specific_heat_cp > 0.0F &&
         options.stratosphere_base_pressure_pa >= 0.0F;
}

[[nodiscard]] constexpr bool valid_coefficients(
    const BaseStateVerticalCoefficientView coefficients,
    const std::int32_t expected_count) noexcept {
  return coefficients.values != nullptr && coefficients.count == expected_count &&
         expected_count > 0;
}

[[nodiscard]] bool finite_coefficients(
    const BaseStateVerticalCoefficientView coefficients) noexcept {
  for (std::int32_t k = 0; k < coefficients.count; ++k) {
    if (!std::isfinite(coefficients.values[k])) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] constexpr bool window_fits_active_shape(
    const nest::RemapWindow& window,
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return window.stagger == nest::HorizontalStagger::mass &&
         window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0 &&
         window.new_i_begin + window.extent_i <= nx &&
         window.new_j_begin + window.extent_j <= ny;
}

void append_region(
    ExposedRegionSet& set,
    const nest::HorizontalStagger stagger,
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
  set.horizontal_cell_count += point_count_2d_u64(extent_i, extent_j);
}

[[nodiscard]] ExposedRegionSet exposed_regions_from_overlap(
    const nest::RemapWindow& overlap,
    const std::int32_t active_nx_value,
    const std::int32_t active_ny_value) noexcept {
  ExposedRegionSet set{};
  const auto overlap_i0 = overlap.new_i_begin;
  const auto overlap_j0 = overlap.new_j_begin;
  const auto overlap_i1 = overlap.new_i_begin + overlap.extent_i;
  const auto overlap_j1 = overlap.new_j_begin + overlap.extent_j;

  append_region(set, overlap.stagger, 0, 0, active_nx_value, overlap_j0);
  append_region(
      set,
      overlap.stagger,
      0,
      overlap_j1,
      active_nx_value,
      active_ny_value - overlap_j1);
  append_region(set, overlap.stagger, 0, overlap_j0, overlap_i0, overlap.extent_j);
  append_region(
      set,
      overlap.stagger,
      overlap_i1,
      overlap_j0,
      active_nx_value - overlap_i1,
      overlap.extent_j);
  return set;
}

[[nodiscard]] double wrf_base_temperature_k(
    const double pb,
    const KrosaMassBaseStateReconstructionOptions& options) noexcept {
  const auto p00 = static_cast<double>(options.sea_level_base_pressure_pa);
  const auto t00 = static_cast<double>(options.sea_level_base_temperature_k);
  const auto lapse = static_cast<double>(options.base_lapse_k);
  const auto tiso = static_cast<double>(options.isothermal_temperature_k);
  const auto troposphere_temperature = t00 + lapse * std::log(pb / p00);
  return troposphere_temperature > tiso ? troposphere_temperature : tiso;
}

[[nodiscard]] double recomputed_t_init(
    const double pb,
    const KrosaMassBaseStateReconstructionOptions& options) noexcept {
  const auto p00 = static_cast<double>(options.sea_level_base_pressure_pa);
  const auto rd_over_cp =
      static_cast<double>(options.dry_air_gas_constant) /
      static_cast<double>(options.specific_heat_cp);
  return wrf_base_temperature_k(pb, options) * std::pow(p00 / pb, rd_over_cp) -
         static_cast<double>(options.base_potential_temperature_k);
}

[[nodiscard]] double recomputed_alb(
    const double pb,
    const double t_init,
    const KrosaMassBaseStateReconstructionOptions& options) noexcept {
  const auto p0 = static_cast<double>(options.reference_pressure_pa);
  const auto theta_base =
      t_init + static_cast<double>(options.base_potential_temperature_k);
  return (static_cast<double>(options.dry_air_gas_constant) / p0) *
         theta_base * std::pow(pb / p0, static_cast<double>(options.cvpm));
}

[[nodiscard]] FieldView2D<const float> terrain_view(
    const io::KrosaPressureRefreshMetadata& metadata) noexcept {
  return {
      metadata.terrain_height_m.data(),
      metadata.terrain_nx,
      metadata.terrain_ny,
      1,
      metadata.terrain_nx,
      {}};
}

[[nodiscard]] BaseStateVerticalCoefficientView coefficient_view(
    const std::vector<float>& values) noexcept {
  return {values.data(), vector_size_or_negative(values.size())};
}

[[nodiscard]] KrosaMassBaseStateReconstructionInputs make_inputs(
    const io::KrosaPressureRefreshMetadata& metadata,
    const FieldView2D<const float> terrain_height_m) noexcept {
  return {
      terrain_height_m,
      coefficient_view(metadata.c3h),
      coefficient_view(metadata.c4h),
      coefficient_view(metadata.c3f),
      coefficient_view(metadata.c4f),
      metadata.p_top_pa};
}

[[nodiscard]] KrosaMassBaseStateReconstructionOutputs make_outputs(
    FieldStorage3D<float>& pb,
    FieldStorage3D<float>& t_init,
    FieldStorage2D<float>& mub,
    FieldStorage3D<float>& alb,
    FieldStorage3D<float>& phb) noexcept {
  return {pb.view(), t_init.view(), mub.view(), alb.view(), phb.view()};
}

[[nodiscard]] KrosaBaseStateProviderReport base_report(
    const Grid& grid,
    const io::KrosaPressureRefreshMetadata& metadata,
    const FieldView2D<const float> terrain_height_m,
    const bool terrain_override_used,
    const std::string_view terrain_source_name,
    const std::string_view terrain_provenance) {
  const auto& config = grid.config();
  KrosaBaseStateProviderReport report{};
  report.active_nx = config.mass_nx;
  report.active_ny = config.mass_ny;
  report.mass_nz = config.mass_nz;
  report.full_nz = config.full_nz;
  report.expected_mass_level_count = config.mass_nz;
  report.expected_full_level_count = config.mass_nz + 1;
  report.c3f_count = vector_size_or_negative(metadata.c3f.size());
  report.c4f_count = vector_size_or_negative(metadata.c4f.size());
  report.c3h_count = vector_size_or_negative(metadata.c3h.size());
  report.c4h_count = vector_size_or_negative(metadata.c4h.size());
  report.terrain_override_used = terrain_override_used;
  if (terrain_override_used) {
    report.terrain_nx = active_nx(terrain_height_m);
    report.terrain_ny = active_ny(terrain_height_m);
    if (report.terrain_nx > 0 && report.terrain_ny > 0) {
      report.terrain_point_count =
          point_count_2d(report.terrain_nx, report.terrain_ny);
    }
  } else {
    report.terrain_nx = metadata.terrain_nx;
    report.terrain_ny = metadata.terrain_ny;
    report.terrain_point_count = metadata.terrain_height_m.size();
  }
  report.terrain_source_name = std::string(terrain_source_name);
  if (terrain_override_used) {
    report.terrain_provenance = std::string(terrain_provenance);
  } else {
    report.terrain_provenance = "metadata";
    if (!report.terrain_source_name.empty()) {
      report.terrain_provenance += ":";
      report.terrain_provenance += report.terrain_source_name;
    }
  }
  if (config.mass_nx > 0 && config.mass_ny > 0) {
    report.expected_terrain_point_count =
        point_count_2d(config.mass_nx, config.mass_ny);
  }
  report.p_top_present =
      metadata.p_top_source != io::PressureRefreshPTopSource::missing;
  return report;
}

[[nodiscard]] bool coefficient_counts_match(
    const KrosaBaseStateProviderReport& report) noexcept {
  return report.c3f_count == report.expected_full_level_count &&
         report.c4f_count == report.expected_full_level_count &&
         report.c3h_count == report.expected_mass_level_count &&
         report.c4h_count == report.expected_mass_level_count;
}

[[nodiscard]] nest::NestResult validate_terrain_override(
    const GridConfig& config,
    const FieldView2D<const float> terrain_height_m) noexcept {
  if (!valid_canonical_terrain_view(terrain_height_m)) {
    return result(
        nest::NestStatus::invalid_contract,
        "base-state provider terrain override must be a non-null canonical 2D view");
  }

  if (active_nx(terrain_height_m) != config.mass_nx ||
      active_ny(terrain_height_m) != config.mass_ny) {
    return result(
        nest::NestStatus::invalid_contract,
        "base-state provider terrain override active shape must match grid");
  }

  for (std::int32_t j = 0; j < config.mass_ny; ++j) {
    for (std::int32_t i = 0; i < config.mass_nx; ++i) {
      const auto value = terrain_height_m(
          terrain_height_m.halo.i_lower + i,
          terrain_height_m.halo.j_lower + j);
      if (!std::isfinite(value) || !sane_terrain_height_m(value)) {
        return result(
            nest::NestStatus::invalid_contract,
            "base-state provider terrain override values must be finite and sane");
      }
    }
  }

  return ok_result();
}

}  // namespace

KrosaBaseStateProviderReport KrosaBaseStateProvider::reconstruct(
    const Grid& grid,
    const io::KrosaPressureRefreshMetadata& metadata,
    const KrosaMassBaseStateReconstructionOptions options) {
  return reconstruct_with_terrain(
      grid,
      metadata,
      terrain_view(metadata),
      false,
      metadata.terrain_source_name,
      {},
      options);
}

KrosaBaseStateProviderReport KrosaBaseStateProvider::reconstruct(
    const Grid& grid,
    const io::KrosaPressureRefreshMetadata& metadata,
    const KrosaBaseStateProviderTerrainOverride& terrain_override,
    const KrosaMassBaseStateReconstructionOptions options) {
  return reconstruct_with_terrain(
      grid,
      metadata,
      terrain_override.terrain_height_m,
      true,
      terrain_override.source_name,
      terrain_override.provenance,
      options);
}

KrosaBaseStateProviderReport KrosaBaseStateProvider::reconstruct_with_terrain(
    const Grid& grid,
    const io::KrosaPressureRefreshMetadata& metadata,
    const FieldView2D<const float> terrain_height_m,
    const bool terrain_override_used,
    const std::string_view terrain_source_name,
    const std::string_view terrain_provenance,
    const KrosaMassBaseStateReconstructionOptions options) {
  clear_views();

  auto report = base_report(
      grid,
      metadata,
      terrain_height_m,
      terrain_override_used,
      terrain_source_name,
      terrain_provenance);
  const auto& config = grid.config();

  if (!valid_grid_config(config)) {
    report.result = result(
        nest::NestStatus::invalid_configuration,
        "base-state provider requires positive mass shape and full_nz = mass_nz + 1");
    return report;
  }

  if (terrain_override_used) {
    if (const auto validation =
            validate_terrain_override(config, terrain_height_m);
        !validation.ok()) {
      report.result = validation;
      return report;
    }
  } else if (
      metadata.terrain_nx != config.mass_nx ||
      metadata.terrain_ny != config.mass_ny ||
      metadata.terrain_height_m.size() != report.expected_terrain_point_count) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "base-state provider terrain metadata must match grid horizontal shape");
    return report;
  }

  if (!coefficient_counts_match(report)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "base-state provider coefficient counts must match grid vertical levels");
    return report;
  }

  if (!report.p_top_present || !std::isfinite(metadata.p_top_pa) ||
      metadata.p_top_pa < 0.0F) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "base-state provider P_TOP must be present, finite, and non-negative");
    return report;
  }

  pb_ = FieldStorage3D<float>(grid.mass_layout());
  t_init_ = FieldStorage3D<float>(grid.mass_layout());
  mub_ = FieldStorage2D<float>(grid.surface_layout());
  alb_ = FieldStorage3D<float>(grid.mass_layout());
  phb_ = FieldStorage3D<float>(grid.w_layout());
  report.allocated_buffers = true;

  report.reconstruction = reconstruct_krosa_mass_base_state(
      make_inputs(metadata, terrain_height_m),
      make_outputs(pb_, t_init_, mub_, alb_, phb_),
      options);
  report.result = report.reconstruction.result;
  report.wrote_pb = report.reconstruction.wrote_pb;
  report.wrote_t_init = report.reconstruction.wrote_t_init;
  report.wrote_mub = report.reconstruction.wrote_mub;
  report.wrote_alb = report.reconstruction.wrote_alb;
  report.wrote_phb = report.reconstruction.wrote_phb;
  if (report.ok()) {
    report.result = ok_result();
  }
  return report;
}

void KrosaBaseStateProvider::clear_views() noexcept {
  pb_ = FieldStorage3D<float>{};
  t_init_ = FieldStorage3D<float>{};
  mub_ = FieldStorage2D<float>{};
  alb_ = FieldStorage3D<float>{};
  phb_ = FieldStorage3D<float>{};
}

KrosaBaseStateProviderViews KrosaBaseStateProvider::views() const noexcept {
  return {pb_.view(), t_init_.view(), mub_.view(), alb_.view(), phb_.view()};
}

KrosaExposedBaseStateRecomputeReport recompute_exposed_base_state_from_mub(
    const nest::RemapWindow& overlap_window,
    const KrosaExposedBaseStateRecomputeInputs& inputs,
    const KrosaExposedBaseStateRecomputeOutputs& outputs,
    const KrosaMassBaseStateReconstructionOptions options) noexcept {
  KrosaExposedBaseStateRecomputeReport report{};
  report.active_nx = active_nx(outputs.pb);
  report.active_ny = active_ny(outputs.pb);
  report.active_nz = active_nz(outputs.pb);
  report.c3h_count = inputs.c3h.count;
  report.c4h_count = inputs.c4h.count;

  if (!valid_canonical_view(inputs.mub) ||
      !valid_canonical_view(outputs.pb) ||
      !valid_canonical_view(outputs.t_init) ||
      !valid_canonical_view(outputs.alb)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "exposed MUB base-state recompute requires non-null canonical views");
    return report;
  }

  if (!same_active_shape(outputs.pb, outputs.t_init) ||
      !same_active_shape(outputs.pb, outputs.alb) ||
      !same_horizontal_active_shape(inputs.mub, outputs.pb)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "exposed MUB base-state recompute views must share active mass shape");
    return report;
  }

  if (!window_fits_active_shape(
          overlap_window, active_nx(outputs.pb), active_ny(outputs.pb))) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "exposed MUB base-state recompute requires a valid mass overlap window");
    return report;
  }

  if (!valid_coefficients(inputs.c3h, active_nz(outputs.pb)) ||
      !valid_coefficients(inputs.c4h, active_nz(outputs.pb))) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "exposed MUB base-state recompute C3H/C4H counts must match mass levels");
    return report;
  }

  if (!valid_options(options) || !std::isfinite(inputs.p_top_pa) ||
      inputs.p_top_pa < 0.0F) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "exposed MUB base-state recompute constants and P_TOP must be finite and valid");
    return report;
  }

  if (!finite_coefficients(inputs.c3h) ||
      !finite_coefficients(inputs.c4h)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "exposed MUB base-state recompute vertical coefficients must be finite");
    return report;
  }

  const auto regions = exposed_regions_from_overlap(
      overlap_window, active_nx(outputs.pb), active_ny(outputs.pb));
  report.exposed_region_count = static_cast<std::uint32_t>(regions.count);
  report.exposed_mass_cell_count = regions.horizontal_cell_count;
  report.recomputed_point_count =
      regions.horizontal_cell_count *
      static_cast<std::uint64_t>(active_nz(outputs.pb));

  if (report.recomputed_point_count == 0) {
    report.result = ok_result();
    return report;
  }

  const auto p_top = static_cast<double>(inputs.p_top_pa);
  const auto p_strat =
      static_cast<double>(options.stratosphere_base_pressure_pa);

  for (std::uint8_t region_index = 0; region_index < regions.count;
       ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t i = region.new_i_begin;
           i < region.new_i_begin + region.extent_i;
           ++i) {
        const auto mub = static_cast<double>(inputs.mub(
            inputs.mub.halo.i_lower + i,
            inputs.mub.halo.j_lower + j));
        if (!std::isfinite(mub)) {
          ++report.invalid_column_count;
          report.invalid_point_count +=
              static_cast<std::uint64_t>(active_nz(outputs.pb));
          continue;
        }

        for (std::int32_t k = 0; k < active_nz(outputs.pb); ++k) {
          const auto pb =
              static_cast<double>(inputs.c3h.values[k]) * mub +
              static_cast<double>(inputs.c4h.values[k]) + p_top;
          if (!std::isfinite(pb) || pb <= 0.0) {
            ++report.invalid_point_count;
            continue;
          }
          if (p_strat > 0.0 && pb < p_strat) {
            ++report.unsupported_stratosphere_point_count;
          }
        }
      }
    }
  }

  if (report.invalid_column_count != 0 || report.invalid_point_count != 0) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "exposed MUB base-state recompute inputs produce invalid pressure");
    return report;
  }

  if (report.unsupported_stratosphere_point_count != 0) {
    report.result = result(
        nest::NestStatus::not_implemented,
        "exposed MUB base-state recompute stratosphere pressure branch is not implemented");
    return report;
  }

  for (std::uint8_t region_index = 0; region_index < regions.count;
       ++region_index) {
    const auto& region = regions.regions[region_index];
    for (std::int32_t j = region.new_j_begin;
         j < region.new_j_begin + region.extent_j;
         ++j) {
      for (std::int32_t k = 0; k < active_nz(outputs.pb); ++k) {
        for (std::int32_t i = region.new_i_begin;
             i < region.new_i_begin + region.extent_i;
             ++i) {
          const auto mub = static_cast<double>(inputs.mub(
              inputs.mub.halo.i_lower + i,
              inputs.mub.halo.j_lower + j));
          const auto pb =
              static_cast<double>(inputs.c3h.values[k]) * mub +
              static_cast<double>(inputs.c4h.values[k]) + p_top;
          const auto t_init = recomputed_t_init(pb, options);
          const auto alb = recomputed_alb(pb, t_init, options);

          outputs.pb(
              outputs.pb.halo.i_lower + i,
              outputs.pb.halo.j_lower + j,
              outputs.pb.halo.k_lower + k) = static_cast<float>(pb);
          outputs.t_init(
              outputs.t_init.halo.i_lower + i,
              outputs.t_init.halo.j_lower + j,
              outputs.t_init.halo.k_lower + k) =
              static_cast<float>(t_init);
          outputs.alb(
              outputs.alb.halo.i_lower + i,
              outputs.alb.halo.j_lower + j,
              outputs.alb.halo.k_lower + k) = static_cast<float>(alb);
        }
      }
    }
  }

  report.result = ok_result();
  report.read_interpolated_mub = true;
  report.wrote_pb = true;
  report.wrote_t_init = true;
  report.wrote_alb = true;
  report.pb_recomputed_point_count = report.recomputed_point_count;
  report.t_init_recomputed_point_count = report.recomputed_point_count;
  report.alb_recomputed_point_count = report.recomputed_point_count;
  return report;
}

}  // namespace tywrf::dynamics
