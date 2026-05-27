#include "tywrf/dynamics/base_state_provider.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace tywrf::dynamics {
namespace {

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
[[nodiscard]] constexpr bool valid_canonical_terrain_view(
    const FieldView2D<Real>& field) noexcept {
  return field.data != nullptr && field.nx > 0 && field.ny > 0 &&
         field.stride_i == 1 && field.stride_j == field.nx &&
         field.halo.i_lower >= 0 && field.halo.i_upper >= 0 &&
         field.halo.j_lower >= 0 && field.halo.j_upper >= 0 &&
         active_nx(field) > 0 && active_ny(field) > 0;
}

[[nodiscard]] constexpr bool sane_terrain_height_m(const float value) noexcept {
  return value >= -1000.0F && value <= 20000.0F;
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

}  // namespace tywrf::dynamics
