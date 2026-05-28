#include "tywrf/nest/base_state_exchange_adapter.hpp"

#include "tywrf/dynamics/base_state_provider.hpp"

#include <array>
#include <cmath>
#include <cstdint>

namespace tywrf::nest {
namespace {

struct ExposedRegionSet {
  std::array<RemapWindow, 4> regions{};
  std::uint8_t count = 0;
  std::uint64_t horizontal_cell_count = 0;
};

[[nodiscard]] constexpr NestResult ok_result() noexcept {
  return {NestStatus::ok, "ok"};
}

[[nodiscard]] constexpr NestResult result(
    const NestStatus status,
    const char* message) noexcept {
  return {status, message};
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
    const FieldView2D<LhsReal>& lhs,
    const FieldView2D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
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

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_horizontal_active_shape(
    const FieldView3D<LhsReal>& lhs,
    const FieldView2D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

template <typename LhsReal, typename RhsReal>
[[nodiscard]] constexpr bool same_horizontal_active_shape(
    const FieldView3D<LhsReal>& lhs,
    const FieldView3D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

[[nodiscard]] constexpr bool window_fits_active_shape(
    const RemapWindow& window,
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0 &&
         window.new_i_begin + window.extent_i <= nx &&
         window.new_j_begin + window.extent_j <= ny;
}

[[nodiscard]] constexpr bool mass_window_fits_active_shape(
    const RemapWindow& window,
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return window.stagger == HorizontalStagger::mass &&
         window_fits_active_shape(window, nx, ny);
}

template <typename Real>
[[nodiscard]] NestResult validate_window(
    const RemapWindow& window,
    const FieldView2D<Real>& field,
    const char* message) noexcept {
  if (!window_fits_active_shape(window, active_nx(field), active_ny(field))) {
    return result(NestStatus::invalid_contract, message);
  }
  return ok_result();
}

template <typename Real>
[[nodiscard]] NestResult validate_window(
    const RemapWindow& window,
    const FieldView3D<Real>& field,
    const char* message) noexcept {
  if (!window_fits_active_shape(window, active_nx(field), active_ny(field))) {
    return result(NestStatus::invalid_contract, message);
  }
  return ok_result();
}

[[nodiscard]] bool valid_options(
    const dynamics::KrosaMassBaseStateReconstructionOptions& options) noexcept {
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
    const dynamics::BaseStateVerticalCoefficientView coefficients,
    const std::int32_t expected_count) noexcept {
  return coefficients.values != nullptr && coefficients.count == expected_count &&
         expected_count > 0;
}

[[nodiscard]] bool finite_coefficients(
    const dynamics::BaseStateVerticalCoefficientView coefficients) noexcept {
  for (std::int32_t k = 0; k < coefficients.count; ++k) {
    if (!std::isfinite(coefficients.values[k])) {
      return false;
    }
  }
  return true;
}

void append_region(
    ExposedRegionSet& set,
    const HorizontalStagger stagger,
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
  set.horizontal_cell_count +=
      static_cast<std::uint64_t>(extent_i) * static_cast<std::uint64_t>(extent_j);
}

[[nodiscard]] ExposedRegionSet exposed_regions_from_overlap(
    const RemapWindow& overlap,
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

[[nodiscard]] bool horizontally_exposed(
    const RemapWindow& overlap,
    const std::int32_t active_i,
    const std::int32_t active_j) noexcept {
  const bool inside_overlap_i =
      active_i >= overlap.new_i_begin &&
      active_i < overlap.new_i_begin + overlap.extent_i;
  const bool inside_overlap_j =
      active_j >= overlap.new_j_begin &&
      active_j < overlap.new_j_begin + overlap.extent_j;
  return !(inside_overlap_i && inside_overlap_j);
}

[[nodiscard]] float prospective_mub_after_exchange(
    const RemapWindow& surface_overlap,
    const ExposedBaseStateViews<const float>& source,
    const ExposedBaseStateViews<float>& child,
    const std::int32_t active_i,
    const std::int32_t active_j) noexcept {
  if (horizontally_exposed(surface_overlap, active_i, active_j)) {
    return source.mub(
        source.mub.halo.i_lower + active_i,
        source.mub.halo.j_lower + active_j);
  }
  return child.mub(
      child.mub.halo.i_lower + active_i,
      child.mub.halo.j_lower + active_j);
}

[[nodiscard]] constexpr FieldView2D<const float> const_view(
    const FieldView2D<float>& field) noexcept {
  return {
      field.data,
      field.nx,
      field.ny,
      field.stride_i,
      field.stride_j,
      field.halo};
}

[[nodiscard]] NestResult validate_exchange_inputs(
    const RemapPlan& remap_plan,
    const ExposedBaseStateViews<const float>& source,
    const ExposedBaseStateViews<float>& child) noexcept {
  if (!valid_canonical_view(source.phb) || !valid_canonical_view(child.phb) ||
      !valid_canonical_view(source.mub) || !valid_canonical_view(child.mub) ||
      !valid_canonical_view(source.ht) || !valid_canonical_view(child.ht) ||
      !valid_canonical_view(child.pb) || !valid_canonical_view(child.t_init) ||
      !valid_canonical_view(child.alb)) {
    return result(
        NestStatus::invalid_contract,
        "adapter requires D68 non-null canonical exposed base-state views");
  }

  if (!same_active_shape(source.phb, child.phb) ||
      !same_active_shape(source.mub, child.mub) ||
      !same_active_shape(source.ht, child.ht)) {
    return result(
        NestStatus::invalid_contract,
        "adapter requires D68 source and child exposed views to share active shapes");
  }

  if (!same_active_shape(child.pb, child.t_init) ||
      !same_active_shape(child.pb, child.alb) ||
      !same_horizontal_active_shape(child.pb, child.mub) ||
      !same_horizontal_active_shape(child.pb, child.ht) ||
      !same_horizontal_active_shape(child.pb, child.phb)) {
    return result(
        NestStatus::invalid_contract,
        "adapter requires child mass fields to share D68/D69 active shapes");
  }

  if (const auto validation = validate_window(
          remap_plan.w_full,
          child.phb,
          "adapter D68 PHB remap window is outside child active extents");
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_window(
          remap_plan.surface,
          child.mub,
          "adapter D68 MUB remap window is outside child active extents");
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_window(
          remap_plan.surface,
          child.ht,
          "adapter D68 HT/HGT remap window is outside child active extents");
      !validation.ok()) {
    return validation;
  }
  if (const auto validation = validate_window(
          remap_plan.mass,
          child.pb,
          "adapter D68 mass remap window is outside child active extents");
      !validation.ok()) {
    return validation;
  }

  return ok_result();
}

[[nodiscard]] NestResult validate_recompute_inputs(
    const RemapPlan& remap_plan,
    const ExposedBaseStateViews<const float>& source,
    const ExposedBaseStateViews<float>& child,
    const ExposedBaseStateExchangeAdapterInputs& inputs,
    const dynamics::KrosaMassBaseStateReconstructionOptions& options,
    ExposedBaseStateExchangeAdapterReport& report) noexcept {
  report.active_nx = active_nx(child.pb);
  report.active_ny = active_ny(child.pb);
  report.active_nz = active_nz(child.pb);
  report.c3h_count = inputs.c3h.count;
  report.c4h_count = inputs.c4h.count;

  if (!mass_window_fits_active_shape(
          remap_plan.mass, active_nx(child.pb), active_ny(child.pb))) {
    return result(
        NestStatus::invalid_contract,
        "adapter D69 recompute requires a valid mass overlap window");
  }

  if (!valid_coefficients(inputs.c3h, active_nz(child.pb)) ||
      !valid_coefficients(inputs.c4h, active_nz(child.pb))) {
    return result(
        NestStatus::invalid_contract,
        "adapter D69 C3H/C4H counts must match mass levels");
  }

  if (!valid_options(options) || !std::isfinite(inputs.p_top_pa) ||
      inputs.p_top_pa < 0.0F) {
    return result(
        NestStatus::invalid_contract,
        "adapter D69 constants and P_TOP must be finite and valid");
  }

  if (!finite_coefficients(inputs.c3h) || !finite_coefficients(inputs.c4h)) {
    return result(
        NestStatus::invalid_contract,
        "adapter D69 vertical coefficients must be finite");
  }

  const auto regions = exposed_regions_from_overlap(
      remap_plan.mass, active_nx(child.pb), active_ny(child.pb));
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
        const auto mub = static_cast<double>(
            prospective_mub_after_exchange(remap_plan.surface, source, child, i, j));
        if (!std::isfinite(mub)) {
          ++report.invalid_column_count;
          report.invalid_point_count +=
              static_cast<std::uint64_t>(active_nz(child.pb));
          continue;
        }

        for (std::int32_t k = 0; k < active_nz(child.pb); ++k) {
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
    return result(
        NestStatus::invalid_contract,
        "adapter D69 preflight inputs produce invalid pressure");
  }

  if (report.unsupported_stratosphere_point_count != 0) {
    return result(
        NestStatus::not_implemented,
        "adapter D69 stratosphere pressure branch is not implemented");
  }

  return ok_result();
}

void merge_exchange_report(
    ExposedBaseStateExchangeAdapterReport& report,
    const ExposedBaseStateExchangeReport& exchange) noexcept {
  report.exchange_result = exchange.result;
  report.rebalance_requested = exchange.rebalance_requested;
  report.rebalance_applied = exchange.rebalance_applied;
  report.exchange_exposed_region_count = exchange.exposed_region_count;
  report.exposed_region_count += exchange.exposed_region_count;
  report.exposed_mass_cell_count = exchange.exposed_mass_cell_count;
  report.exposed_surface_cell_count = exchange.exposed_surface_cell_count;
  report.exposed_w_full_column_count = exchange.exposed_w_full_column_count;
  report.wrote_overlap = report.wrote_overlap || exchange.wrote_overlap;
  report.wrote_halo = report.wrote_halo || exchange.wrote_halo;
  report.wrote_phb = report.wrote_phb || exchange.wrote_phb;
  report.wrote_mub = report.wrote_mub || exchange.wrote_mub;
  report.wrote_ht = report.wrote_ht || exchange.wrote_ht;
  report.phb_written_point_count = exchange.phb_written_point_count;
  report.mub_written_cell_count = exchange.mub_written_cell_count;
  report.ht_written_cell_count = exchange.ht_written_cell_count;
  report.direct_write_point_count = exchange.direct_write_point_count;
  report.exchange_recompute_mark_count = exchange.recompute_mark_count;
}

void merge_recompute_report(
    ExposedBaseStateExchangeAdapterReport& report,
    const dynamics::KrosaExposedBaseStateRecomputeReport& recompute) noexcept {
  report.recompute_result = recompute.result;
  report.active_nx = recompute.active_nx;
  report.active_ny = recompute.active_ny;
  report.active_nz = recompute.active_nz;
  report.c3h_count = recompute.c3h_count;
  report.c4h_count = recompute.c4h_count;
  report.recompute_exposed_region_count = recompute.exposed_region_count;
  report.exposed_region_count += recompute.exposed_region_count;
  report.exposed_mass_cell_count = recompute.exposed_mass_cell_count;
  report.recomputed_point_count = recompute.recomputed_point_count;
  report.wrote_pb = report.wrote_pb || recompute.wrote_pb;
  report.wrote_t_init = report.wrote_t_init || recompute.wrote_t_init;
  report.wrote_alb = report.wrote_alb || recompute.wrote_alb;
  report.pb_recomputed_point_count = recompute.pb_recomputed_point_count;
  report.t_init_recomputed_point_count =
      recompute.t_init_recomputed_point_count;
  report.alb_recomputed_point_count = recompute.alb_recomputed_point_count;
  report.invalid_column_count += recompute.invalid_column_count;
  report.invalid_point_count += recompute.invalid_point_count;
  report.unsupported_stratosphere_point_count +=
      recompute.unsupported_stratosphere_point_count;
}

}  // namespace

ExposedBaseStateExchangeAdapterReport apply_exposed_base_state_exchange_adapter(
    const RemapPlan& remap_plan,
    const ExposedBaseStateViews<const float>& source_interpolated_base_state,
    const ExposedBaseStateViews<float>& child_base_state,
    const ExposedBaseStateExchangeAdapterInputs& inputs,
    const ExposedBaseStateExchangeAdapterOptions options) noexcept {
  ExposedBaseStateExchangeAdapterReport report{};
  report.rebalance_requested = options.exchange.rebalance;

  if (!remap_plan.ok()) {
    report.result = remap_plan.result;
    report.exchange_result = remap_plan.result;
    return report;
  }

  if (options.exchange.rebalance) {
    report.result = result(
        NestStatus::not_implemented,
        "diagnostic exposed base-state adapter does not implement rebalancing");
    report.exchange_result = report.result;
    return report;
  }

  if (const auto validation = validate_exchange_inputs(
          remap_plan, source_interpolated_base_state, child_base_state);
      !validation.ok()) {
    report.result = validation;
    report.exchange_result = validation;
    return report;
  }

  if (const auto validation = validate_recompute_inputs(
          remap_plan,
          source_interpolated_base_state,
          child_base_state,
          inputs,
          options.recompute,
          report);
      !validation.ok()) {
    report.result = validation;
    report.recompute_result = validation;
    return report;
  }

  const auto exchange_report = apply_exposed_base_state_exchange(
      remap_plan,
      source_interpolated_base_state,
      child_base_state,
      options.exchange);
  report.called_d68_exchange = true;
  merge_exchange_report(report, exchange_report);
  if (!exchange_report.ok()) {
    report.result = exchange_report.result;
    return report;
  }

  const auto recompute_report =
      dynamics::recompute_exposed_base_state_from_mub(
          remap_plan.mass,
          {const_view(child_base_state.mub),
           inputs.c3h,
           inputs.c4h,
           inputs.p_top_pa},
          {child_base_state.pb, child_base_state.t_init, child_base_state.alb},
          options.recompute);
  report.called_d69_recompute = true;
  merge_recompute_report(report, recompute_report);
  if (!recompute_report.ok()) {
    report.result = recompute_report.result;
    return report;
  }

  report.result = ok_result();
  return report;
}

}  // namespace tywrf::nest
