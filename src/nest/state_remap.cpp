#include "tywrf/nest/state_remap.hpp"

#include <cstdint>

namespace tywrf::nest {
namespace {

[[nodiscard]] constexpr NestResult ok_result() noexcept {
  return {NestStatus::ok, "ok"};
}

[[nodiscard]] constexpr NestResult result(
    const NestStatus status,
    const char* message) noexcept {
  return {status, message};
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(const FieldView2D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(const FieldView2D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(const FieldView3D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(const FieldView3D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nz(const FieldView3D<Real>& field) noexcept {
  return field.nz - field.halo.k_lower - field.halo.k_upper;
}

[[nodiscard]] constexpr bool positive_window(const RemapWindow& window) noexcept {
  return window.old_i_begin >= 0 && window.old_j_begin >= 0 &&
         window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0;
}

template <typename OldReal, typename NewReal>
[[nodiscard]] bool valid_2d_window(
    const RemapWindow& window,
    const FieldView2D<OldReal>& old_field,
    const FieldView2D<NewReal>& new_field) noexcept {
  return positive_window(window) &&
         window.old_i_begin + window.extent_i <= active_nx(old_field) &&
         window.old_j_begin + window.extent_j <= active_ny(old_field) &&
         window.new_i_begin + window.extent_i <= active_nx(new_field) &&
         window.new_j_begin + window.extent_j <= active_ny(new_field);
}

template <typename OldReal, typename NewReal>
[[nodiscard]] bool valid_3d_window(
    const RemapWindow& window,
    const FieldView3D<OldReal>& old_field,
    const FieldView3D<NewReal>& new_field) noexcept {
  return active_nz(old_field) == active_nz(new_field) && positive_window(window) &&
         window.old_i_begin + window.extent_i <= active_nx(old_field) &&
         window.old_j_begin + window.extent_j <= active_ny(old_field) &&
         window.new_i_begin + window.extent_i <= active_nx(new_field) &&
         window.new_j_begin + window.extent_j <= active_ny(new_field);
}

template <typename Real>
[[nodiscard]] bool window_covers_new(
    const RemapWindow& window,
    const FieldView2D<Real>& field) noexcept {
  return window.new_i_begin == 0 && window.new_j_begin == 0 &&
         window.extent_i == active_nx(field) && window.extent_j == active_ny(field);
}

template <typename Real>
[[nodiscard]] bool window_covers_new(
    const RemapWindow& window,
    const FieldView3D<Real>& field) noexcept {
  return window.new_i_begin == 0 && window.new_j_begin == 0 &&
         window.extent_i == active_nx(field) && window.extent_j == active_ny(field);
}

[[nodiscard]] bool all_new_active_cells_covered(
    const RemapPlan& plan,
    const StateView<float>& new_child) noexcept {
  return window_covers_new(plan.u, new_child.u) &&
         window_covers_new(plan.v, new_child.v) &&
         window_covers_new(plan.w_full, new_child.w) &&
         window_covers_new(plan.mass, new_child.t) &&
         window_covers_new(plan.surface, new_child.mu);
}

[[nodiscard]] NestResult validate_remap_inputs(
    const RemapPlan& plan,
    const StateView<const float>& old_child,
    const StateView<float>& new_child) noexcept {
  if (!plan.ok()) {
    return plan.result;
  }

  if (!valid_3d_window(plan.u, old_child.u, new_child.u) ||
      !valid_3d_window(plan.v, old_child.v, new_child.v) ||
      !valid_3d_window(plan.w_full, old_child.w, new_child.w) ||
      !valid_3d_window(plan.w_full, old_child.ph, new_child.ph) ||
      !valid_3d_window(plan.w_full, old_child.phb, new_child.phb) ||
      !valid_3d_window(plan.mass, old_child.t, new_child.t) ||
      !valid_3d_window(plan.mass, old_child.p, new_child.p) ||
      !valid_3d_window(plan.mass, old_child.pb, new_child.pb) ||
      !valid_3d_window(plan.mass, old_child.qvapor, new_child.qvapor) ||
      !valid_3d_window(plan.mass, old_child.qcloud, new_child.qcloud) ||
      !valid_3d_window(plan.mass, old_child.qrain, new_child.qrain) ||
      !valid_3d_window(plan.mass, old_child.qice, new_child.qice) ||
      !valid_3d_window(plan.mass, old_child.qsnow, new_child.qsnow) ||
      !valid_3d_window(plan.mass, old_child.qgraup, new_child.qgraup) ||
      !valid_3d_window(plan.mass, old_child.qnice, new_child.qnice) ||
      !valid_3d_window(plan.mass, old_child.qnrain, new_child.qnrain) ||
      !valid_2d_window(plan.surface, old_child.mu, new_child.mu) ||
      !valid_2d_window(plan.surface, old_child.mub, new_child.mub) ||
      !valid_2d_window(plan.surface, old_child.psfc, new_child.psfc) ||
      !valid_2d_window(plan.surface, old_child.u10, new_child.u10) ||
      !valid_2d_window(plan.surface, old_child.v10, new_child.v10) ||
      !valid_2d_window(plan.surface, old_child.t2, new_child.t2) ||
      !valid_2d_window(plan.surface, old_child.q2, new_child.q2) ||
      !valid_2d_window(plan.surface, old_child.rainc, new_child.rainc) ||
      !valid_2d_window(plan.surface, old_child.rainnc, new_child.rainnc)) {
    return result(
        NestStatus::invalid_contract,
        "child state remap window is outside active field extents");
  }

  return ok_result();
}

void copy_overlap_2d(
    const RemapWindow& window,
    const FieldView2D<const float>& old_field,
    const FieldView2D<float>& new_field,
    ChildStateRemapReport& report) noexcept {
  const auto old_i0 = old_field.halo.i_lower + window.old_i_begin;
  const auto old_j0 = old_field.halo.j_lower + window.old_j_begin;
  const auto new_i0 = new_field.halo.i_lower + window.new_i_begin;
  const auto new_j0 = new_field.halo.j_lower + window.new_j_begin;

  for (std::int32_t j = 0; j < window.extent_j; ++j) {
    const auto old_j = old_j0 + j;
    const auto new_j = new_j0 + j;
    for (std::int32_t i = 0; i < window.extent_i; ++i) {
      new_field(new_i0 + i, new_j) = old_field(old_i0 + i, old_j);
    }
  }

  ++report.copied_field_count;
  report.copied_point_count +=
      static_cast<std::uint64_t>(window.extent_i) *
      static_cast<std::uint64_t>(window.extent_j);
}

void copy_overlap_3d(
    const RemapWindow& window,
    const FieldView3D<const float>& old_field,
    const FieldView3D<float>& new_field,
    ChildStateRemapReport& report) noexcept {
  const auto old_i0 = old_field.halo.i_lower + window.old_i_begin;
  const auto old_j0 = old_field.halo.j_lower + window.old_j_begin;
  const auto old_k0 = old_field.halo.k_lower;
  const auto new_i0 = new_field.halo.i_lower + window.new_i_begin;
  const auto new_j0 = new_field.halo.j_lower + window.new_j_begin;
  const auto new_k0 = new_field.halo.k_lower;
  const auto nz = active_nz(old_field);

  for (std::int32_t j = 0; j < window.extent_j; ++j) {
    const auto old_j = old_j0 + j;
    const auto new_j = new_j0 + j;
    for (std::int32_t k = 0; k < nz; ++k) {
      const auto old_k = old_k0 + k;
      const auto new_k = new_k0 + k;
      for (std::int32_t i = 0; i < window.extent_i; ++i) {
        new_field(new_i0 + i, new_j, new_k) = old_field(old_i0 + i, old_j, old_k);
      }
    }
  }

  ++report.copied_field_count;
  report.copied_point_count +=
      static_cast<std::uint64_t>(window.extent_i) *
      static_cast<std::uint64_t>(window.extent_j) *
      static_cast<std::uint64_t>(nz);
}

}  // namespace

ChildStateRemapReport remap_child_state_overlap_only(
    const RemapPlan& plan,
    const State<float>& old_child,
    State<float>& new_child) noexcept {
  ChildStateRemapReport report{};
  report.result = ok_result();

  const auto old_view = old_child.view();
  const auto new_view = new_child.view();
  report.needs_parent_fill = true;

  if (const auto validation = validate_remap_inputs(plan, old_view, new_view);
      !validation.ok()) {
    report.result = validation;
    return report;
  }

  report.needs_parent_fill = !all_new_active_cells_covered(plan, new_view);

  copy_overlap_3d(plan.u, old_view.u, new_view.u, report);
  copy_overlap_3d(plan.v, old_view.v, new_view.v, report);
  copy_overlap_3d(plan.w_full, old_view.w, new_view.w, report);
  copy_overlap_3d(plan.w_full, old_view.ph, new_view.ph, report);
  copy_overlap_3d(plan.w_full, old_view.phb, new_view.phb, report);

  copy_overlap_3d(plan.mass, old_view.t, new_view.t, report);
  copy_overlap_3d(plan.mass, old_view.p, new_view.p, report);
  copy_overlap_3d(plan.mass, old_view.pb, new_view.pb, report);
  copy_overlap_3d(plan.mass, old_view.qvapor, new_view.qvapor, report);
  copy_overlap_3d(plan.mass, old_view.qcloud, new_view.qcloud, report);
  copy_overlap_3d(plan.mass, old_view.qrain, new_view.qrain, report);
  copy_overlap_3d(plan.mass, old_view.qice, new_view.qice, report);
  copy_overlap_3d(plan.mass, old_view.qsnow, new_view.qsnow, report);
  copy_overlap_3d(plan.mass, old_view.qgraup, new_view.qgraup, report);
  copy_overlap_3d(plan.mass, old_view.qnice, new_view.qnice, report);
  copy_overlap_3d(plan.mass, old_view.qnrain, new_view.qnrain, report);

  copy_overlap_2d(plan.surface, old_view.mu, new_view.mu, report);
  copy_overlap_2d(plan.surface, old_view.mub, new_view.mub, report);
  copy_overlap_2d(plan.surface, old_view.psfc, new_view.psfc, report);
  copy_overlap_2d(plan.surface, old_view.u10, new_view.u10, report);
  copy_overlap_2d(plan.surface, old_view.v10, new_view.v10, report);
  copy_overlap_2d(plan.surface, old_view.t2, new_view.t2, report);
  copy_overlap_2d(plan.surface, old_view.q2, new_view.q2, report);
  copy_overlap_2d(plan.surface, old_view.rainc, new_view.rainc, report);
  copy_overlap_2d(plan.surface, old_view.rainnc, new_view.rainnc, report);

  return report;
}

}  // namespace tywrf::nest
