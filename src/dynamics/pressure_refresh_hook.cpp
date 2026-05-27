#include "tywrf/dynamics/pressure_refresh_hook.hpp"

#include <algorithm>
#include <cstdint>

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

[[nodiscard]] constexpr bool in_new_overlap_window(
    const nest::RemapWindow& window,
    const std::int32_t i,
    const std::int32_t j) noexcept {
  return i >= window.new_i_begin && i < window.new_i_begin + window.extent_i &&
         j >= window.new_j_begin && j < window.new_j_begin + window.extent_j;
}

template <typename Real>
[[nodiscard]] bool window_fits(
    const nest::RemapWindow& window,
    const FieldView2D<Real>& field,
    const nest::HorizontalStagger stagger) noexcept {
  return window.stagger == stagger && window.old_i_begin >= 0 &&
         window.old_j_begin >= 0 && window.new_i_begin >= 0 &&
         window.new_j_begin >= 0 && window.extent_i > 0 &&
         window.extent_j > 0 &&
         window.old_i_begin + window.extent_i <= active_nx(field) &&
         window.old_j_begin + window.extent_j <= active_ny(field) &&
         window.new_i_begin + window.extent_i <= active_nx(field) &&
         window.new_j_begin + window.extent_j <= active_ny(field);
}

template <typename Real>
[[nodiscard]] bool window_fits(
    const nest::RemapWindow& window,
    const FieldView3D<Real>& field,
    const nest::HorizontalStagger stagger) noexcept {
  return window.stagger == stagger && window.old_i_begin >= 0 &&
         window.old_j_begin >= 0 && window.new_i_begin >= 0 &&
         window.new_j_begin >= 0 && window.extent_i > 0 &&
         window.extent_j > 0 &&
         window.old_i_begin + window.extent_i <= active_nx(field) &&
         window.old_j_begin + window.extent_j <= active_ny(field) &&
         window.new_i_begin + window.extent_i <= active_nx(field) &&
         window.new_j_begin + window.extent_j <= active_ny(field);
}

[[nodiscard]] nest::NestResult validate_sync_contract(
    const nest::RemapPlan& plan,
    const StateView<float>& state,
    const KrosaBaseStateProviderViews& provider) noexcept {
  if (!plan.ok()) {
    return plan.result;
  }

  if (!valid_canonical_view(state.pb) || !valid_canonical_view(state.mub) ||
      !valid_canonical_view(state.phb) || !valid_canonical_view(provider.pb) ||
      !valid_canonical_view(provider.mub) || !valid_canonical_view(provider.phb) ||
      !valid_canonical_view(provider.alb)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh hook requires canonical state and provider views");
  }

  if (!same_active_shape(state.pb, provider.pb) ||
      !same_active_shape(state.mub, provider.mub) ||
      !same_active_shape(state.phb, provider.phb)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh hook provider views must match child state shape");
  }

  if (!window_fits(plan.mass, state.pb, nest::HorizontalStagger::mass) ||
      !window_fits(plan.surface, state.mub, nest::HorizontalStagger::surface) ||
      !window_fits(plan.w_full, state.phb, nest::HorizontalStagger::w_full)) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh hook remap windows are outside active field extents");
  }

  return ok_result();
}

template <typename Real>
[[nodiscard]] std::uint64_t exposed_point_count_2d(
    const nest::RemapWindow& window,
    const FieldView2D<Real>& target) noexcept {
  const auto nx = active_nx(target);
  const auto ny = active_ny(target);
  const auto active_points =
      static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny);
  const auto overlap_points =
      static_cast<std::uint64_t>(window.extent_i) *
      static_cast<std::uint64_t>(window.extent_j);
  return active_points - overlap_points;
}

template <typename Real>
[[nodiscard]] std::uint64_t exposed_point_count_3d(
    const nest::RemapWindow& window,
    const FieldView3D<Real>& target) noexcept {
  const auto nx = active_nx(target);
  const auto ny = active_ny(target);
  const auto active_columns =
      static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny);
  const auto overlap_columns =
      static_cast<std::uint64_t>(window.extent_i) *
      static_cast<std::uint64_t>(window.extent_j);
  return (active_columns - overlap_columns) *
         static_cast<std::uint64_t>(active_nz(target));
}

template <typename Real>
std::uint64_t sync_exposed_2d(
    const nest::RemapWindow& window,
    const FieldView2D<const Real>& source,
    const FieldView2D<Real>& target) noexcept {
  const auto nx = active_nx(target);
  const auto ny = active_ny(target);
  const auto source_i0 = source.halo.i_lower;
  const auto source_j0 = source.halo.j_lower;
  const auto target_i0 = target.halo.i_lower;
  const auto target_j0 = target.halo.j_lower;
  std::uint64_t synced = 0;

  for (std::int32_t j = 0; j < ny; ++j) {
    for (std::int32_t i = 0; i < nx; ++i) {
      if (in_new_overlap_window(window, i, j)) {
        continue;
      }
      target(target_i0 + i, target_j0 + j) =
          source(source_i0 + i, source_j0 + j);
      ++synced;
    }
  }

  return synced;
}

template <typename Real>
std::uint64_t sync_exposed_3d(
    const nest::RemapWindow& window,
    const FieldView3D<const Real>& source,
    const FieldView3D<Real>& target) noexcept {
  const auto nx = active_nx(target);
  const auto ny = active_ny(target);
  const auto nz = active_nz(target);
  const auto source_i0 = source.halo.i_lower;
  const auto source_j0 = source.halo.j_lower;
  const auto source_k0 = source.halo.k_lower;
  const auto target_i0 = target.halo.i_lower;
  const auto target_j0 = target.halo.j_lower;
  const auto target_k0 = target.halo.k_lower;
  std::uint64_t synced = 0;

  for (std::int32_t j = 0; j < ny; ++j) {
    for (std::int32_t k = 0; k < nz; ++k) {
      for (std::int32_t i = 0; i < nx; ++i) {
        if (in_new_overlap_window(window, i, j)) {
          continue;
        }
        target(target_i0 + i, target_j0 + j, target_k0 + k) =
            source(source_i0 + i, source_j0 + j, source_k0 + k);
        ++synced;
      }
    }
  }

  return synced;
}

template <typename Storage>
void copy_storage(const Storage& source, Storage& target) {
  std::copy(source.data(), source.data() + source.size(), target.data());
}

[[nodiscard]] constexpr nest::NestResult
pressure_compute_postcondition_result(
    const PressureRefreshReport& compute_report) noexcept {
  if (!compute_report.ok()) {
    return compute_report.result;
  }

  if (compute_report.refreshed_point_count == 0) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh hook requires at least one refreshed pressure point");
  }

  if (compute_report.invalid_point_count != 0 ||
      compute_report.skipped_point_count != 0) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh hook rejected partial or invalid pressure compute");
  }

  if (compute_report.touched_overlap_cells ||
      compute_report.touched_halo_cells) {
    return result(
        nest::NestStatus::invalid_contract,
        "pressure refresh hook rejected overlap or halo pressure writes");
  }

  return ok_result();
}

}  // namespace

KrosaPressureRefreshHookReport apply_krosa_moving_nest_pressure_refresh_hook(
    const nest::RemapPlan& plan,
    State<float>& new_child,
    const io::KrosaPressureRefreshMetadata& metadata,
    const KrosaPressureRefreshHookOptions options) {
  KrosaPressureRefreshHookReport report{};
  report.result = ok_result();
  report.diagnostic_only = true;
  report.gate_candidate = false;
  report.integrator_output = false;
  report.base_state_sync_dry_run = options.base_state_sync_dry_run;
  report.pressure_compute_dry_run = options.pressure_compute_dry_run;

  if (options.pressure_compute_dry_run && !options.base_state_sync_dry_run) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure compute dry-run requires base-state sync dry-run");
    return report;
  }

  KrosaBaseStateProvider provider;
  if (options.terrain_override == nullptr) {
    report.provider_report =
        provider.reconstruct(new_child.grid, metadata, options.base_state);
  } else {
    report.provider_report = provider.reconstruct(
        new_child.grid,
        metadata,
        *options.terrain_override,
        options.base_state);
  }
  report.provider_ok = report.provider_report.ok();
  if (!report.provider_ok) {
    report.result = report.provider_report.result;
    return report;
  }

  const auto provider_views = provider.views();
  const auto state_view = new_child.view();
  if (const auto validation =
          validate_sync_contract(plan, state_view, provider_views);
      !validation.ok()) {
    report.result = validation;
    return report;
  }
  report.base_state_sync_contract_ok = true;
  report.would_sync_pb_point_count =
      exposed_point_count_3d(plan.mass, state_view.pb);
  report.would_sync_mub_point_count =
      exposed_point_count_2d(plan.surface, state_view.mub);
  report.would_sync_phb_point_count =
      exposed_point_count_3d(plan.w_full, state_view.phb);

  if (options.base_state_sync_dry_run) {
    if (!options.pressure_compute_dry_run) {
      return report;
    }

    FieldStorage3D<float> scratch_p(new_child.p.layout());
    FieldStorage3D<float> scratch_pb(new_child.pb.layout());
    FieldStorage3D<float> scratch_phb(new_child.phb.layout());
    FieldStorage2D<float> scratch_mub(new_child.mub.layout());
    copy_storage(new_child.p, scratch_p);
    copy_storage(new_child.pb, scratch_pb);
    copy_storage(new_child.phb, scratch_phb);
    copy_storage(new_child.mub, scratch_mub);

    auto scratch_view = state_view;
    scratch_view.p = scratch_p.view();
    scratch_view.pb = scratch_pb.view();
    scratch_view.phb = scratch_phb.view();
    scratch_view.mub = scratch_mub.view();

    (void)sync_exposed_3d(plan.mass, provider_views.pb, scratch_view.pb);
    (void)sync_exposed_2d(plan.surface, provider_views.mub, scratch_view.mub);
    (void)sync_exposed_3d(plan.w_full, provider_views.phb, scratch_view.phb);

    const auto staging = make_krosa_pressure_refresh_inputs(
        scratch_view,
        provider_views.alb,
        metadata,
        PressureRefreshAlbSource::base_state_reconstruction_provider);
    report.staging_report = staging.report;
    report.staging_ok = staging.ok();
    if (!report.staging_ok) {
      report.result = report.staging_report.result;
      return report;
    }

    report.pressure_compute_dry_run_called = true;
    report.pressure_compute_dry_run_report =
        refresh_krosa_moving_nest_pressure(
            plan,
            staging.inputs,
            options.pressure_refresh);
    report.result = pressure_compute_postcondition_result(
        report.pressure_compute_dry_run_report);
    report.pressure_compute_dry_run_ok = report.result.ok();
    report.would_refresh_p_point_count =
        report.pressure_compute_dry_run_report.refreshed_point_count;
    report.dry_run_invalid_p_point_count =
        report.pressure_compute_dry_run_report.invalid_point_count;
    report.touched_overlap_cells =
        report.pressure_compute_dry_run_report.touched_overlap_cells;
    report.touched_halo_cells =
        report.pressure_compute_dry_run_report.touched_halo_cells;
    return report;
  }

  report.synced_pb_point_count =
      sync_exposed_3d(plan.mass, provider_views.pb, state_view.pb);
  report.synced_mub_point_count =
      sync_exposed_2d(plan.surface, provider_views.mub, state_view.mub);
  report.synced_phb_point_count =
      sync_exposed_3d(plan.w_full, provider_views.phb, state_view.phb);
  report.base_state_sync_applied =
      report.synced_pb_point_count > 0 || report.synced_mub_point_count > 0 ||
      report.synced_phb_point_count > 0;

  const auto staging = make_krosa_pressure_refresh_inputs(
      state_view,
      provider_views.alb,
      metadata,
      PressureRefreshAlbSource::base_state_reconstruction_provider);
  report.staging_report = staging.report;
  report.staging_ok = staging.ok();
  if (!report.staging_ok) {
    report.result = report.staging_report.result;
    return report;
  }

  report.calls_pressure_refresh_compute = true;
  report.compute_report =
      refresh_krosa_moving_nest_pressure(
          plan,
          staging.inputs,
          options.pressure_refresh);
  report.result = pressure_compute_postcondition_result(report.compute_report);
  report.pressure_refresh_applied = report.result.ok();
  report.touched_overlap_cells = report.compute_report.touched_overlap_cells;
  report.touched_halo_cells = report.compute_report.touched_halo_cells;
  return report;
}

}  // namespace tywrf::dynamics
