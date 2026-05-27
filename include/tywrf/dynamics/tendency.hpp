#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/state.hpp"

#include <cstddef>
#include <cstdint>

namespace tywrf::dynamics {

enum class TendencyApplyStatus {
  ok,
  null_field,
  null_tendency,
  invalid_layout,
  mismatched_layout,
};

struct TendencyApplyReport {
  TendencyApplyStatus status = TendencyApplyStatus::ok;
  std::int64_t active_points = 0;
};

template <typename T>
[[nodiscard]] constexpr bool valid_tendency_view(const FieldView2D<T> view) noexcept {
  const auto i_begin = view.halo.i_lower;
  const auto i_end = view.nx - view.halo.i_upper;
  const auto j_begin = view.halo.j_lower;
  const auto j_end = view.ny - view.halo.j_upper;
  return view.data != nullptr && view.nx >= 0 && view.ny >= 0 && view.stride_i == 1 &&
         view.stride_j == view.nx && view.halo.i_lower >= 0 && view.halo.i_upper >= 0 &&
         view.halo.j_lower >= 0 && view.halo.j_upper >= 0 && i_begin <= i_end &&
         j_begin <= j_end;
}

template <typename T>
[[nodiscard]] constexpr bool valid_tendency_view(const FieldView3D<T> view) noexcept {
  const auto i_begin = view.halo.i_lower;
  const auto i_end = view.nx - view.halo.i_upper;
  const auto j_begin = view.halo.j_lower;
  const auto j_end = view.ny - view.halo.j_upper;
  const auto k_begin = view.halo.k_lower;
  const auto k_end = view.nz - view.halo.k_upper;
  return view.data != nullptr && view.nx >= 0 && view.ny >= 0 && view.nz >= 0 &&
         view.stride_i == 1 && view.stride_k == view.nx &&
         view.stride_j == view.nx * view.nz && view.halo.i_lower >= 0 &&
         view.halo.i_upper >= 0 && view.halo.j_lower >= 0 &&
         view.halo.j_upper >= 0 && view.halo.k_lower >= 0 &&
         view.halo.k_upper >= 0 && i_begin <= i_end && j_begin <= j_end &&
         k_begin <= k_end;
}

template <typename FieldReal, typename TendencyReal>
[[nodiscard]] constexpr bool same_tendency_layout(
    const FieldView2D<FieldReal> field,
    const FieldView2D<TendencyReal> tendency) noexcept {
  return field.nx == tendency.nx && field.ny == tendency.ny &&
         field.stride_i == tendency.stride_i && field.stride_j == tendency.stride_j &&
         field.halo.i_lower == tendency.halo.i_lower &&
         field.halo.i_upper == tendency.halo.i_upper &&
         field.halo.j_lower == tendency.halo.j_lower &&
         field.halo.j_upper == tendency.halo.j_upper;
}

template <typename FieldReal, typename TendencyReal>
[[nodiscard]] constexpr bool same_tendency_layout(
    const FieldView3D<FieldReal> field,
    const FieldView3D<TendencyReal> tendency) noexcept {
  return field.nx == tendency.nx && field.ny == tendency.ny && field.nz == tendency.nz &&
         field.stride_i == tendency.stride_i && field.stride_k == tendency.stride_k &&
         field.stride_j == tendency.stride_j &&
         field.halo.i_lower == tendency.halo.i_lower &&
         field.halo.i_upper == tendency.halo.i_upper &&
         field.halo.j_lower == tendency.halo.j_lower &&
         field.halo.j_upper == tendency.halo.j_upper &&
         field.halo.k_lower == tendency.halo.k_lower &&
         field.halo.k_upper == tendency.halo.k_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int64_t active_point_count(
    const FieldView2D<Real> view) noexcept {
  return static_cast<std::int64_t>(view.nx - view.halo.i_lower - view.halo.i_upper) *
         static_cast<std::int64_t>(view.ny - view.halo.j_lower - view.halo.j_upper);
}

template <typename Real>
[[nodiscard]] constexpr std::int64_t active_point_count(
    const FieldView3D<Real> view) noexcept {
  return static_cast<std::int64_t>(view.nx - view.halo.i_lower - view.halo.i_upper) *
         static_cast<std::int64_t>(view.ny - view.halo.j_lower - view.halo.j_upper) *
         static_cast<std::int64_t>(view.nz - view.halo.k_lower - view.halo.k_upper);
}

template <typename Real>
[[nodiscard]] constexpr TendencyApplyReport validate_tendency_apply(
    const FieldView2D<Real> field,
    const FieldView2D<const Real> tendency) noexcept {
  if (field.data == nullptr) {
    return {TendencyApplyStatus::null_field, 0};
  }
  if (tendency.data == nullptr) {
    return {TendencyApplyStatus::null_tendency, 0};
  }
  if (!valid_tendency_view(field) || !valid_tendency_view(tendency)) {
    return {TendencyApplyStatus::invalid_layout, 0};
  }
  if (!same_tendency_layout(field, tendency)) {
    return {TendencyApplyStatus::mismatched_layout, 0};
  }
  return {TendencyApplyStatus::ok, active_point_count(field)};
}

template <typename Real>
[[nodiscard]] constexpr TendencyApplyReport validate_tendency_apply(
    const FieldView3D<Real> field,
    const FieldView3D<const Real> tendency) noexcept {
  if (field.data == nullptr) {
    return {TendencyApplyStatus::null_field, 0};
  }
  if (tendency.data == nullptr) {
    return {TendencyApplyStatus::null_tendency, 0};
  }
  if (!valid_tendency_view(field) || !valid_tendency_view(tendency)) {
    return {TendencyApplyStatus::invalid_layout, 0};
  }
  if (!same_tendency_layout(field, tendency)) {
    return {TendencyApplyStatus::mismatched_layout, 0};
  }
  return {TendencyApplyStatus::ok, active_point_count(field)};
}

template <typename Real>
[[nodiscard]] constexpr TendencyApplyReport validate_zero_tendency_apply(
    const FieldView2D<Real> field) noexcept {
  if (field.data == nullptr) {
    return {TendencyApplyStatus::null_field, 0};
  }
  if (!valid_tendency_view(field)) {
    return {TendencyApplyStatus::invalid_layout, 0};
  }
  return {TendencyApplyStatus::ok, active_point_count(field)};
}

template <typename Real>
[[nodiscard]] constexpr TendencyApplyReport validate_zero_tendency_apply(
    const FieldView3D<Real> field) noexcept {
  if (field.data == nullptr) {
    return {TendencyApplyStatus::null_field, 0};
  }
  if (!valid_tendency_view(field)) {
    return {TendencyApplyStatus::invalid_layout, 0};
  }
  return {TendencyApplyStatus::ok, active_point_count(field)};
}

[[nodiscard]] constexpr TendencyApplyReport merge_reports(
    const TendencyApplyReport total,
    const TendencyApplyReport next) noexcept {
  if (total.status != TendencyApplyStatus::ok) {
    return total;
  }
  if (next.status != TendencyApplyStatus::ok) {
    return next;
  }
  return {TendencyApplyStatus::ok, total.active_points + next.active_points};
}

template <typename Real>
[[nodiscard]] TendencyApplyReport apply_tendency(
    const FieldView2D<Real> field,
    const FieldView2D<const Real> tendency,
    const Real dt) noexcept {
  const auto validation = validate_tendency_apply(field, tendency);
  if (validation.status != TendencyApplyStatus::ok) {
    return validation;
  }
  for (std::int32_t j = field.halo.j_lower; j < field.ny - field.halo.j_upper; ++j) {
    const auto row = static_cast<std::size_t>(j) * static_cast<std::size_t>(field.stride_j);
    for (std::int32_t i = field.halo.i_lower; i < field.nx - field.halo.i_upper; ++i) {
      const auto idx = row + static_cast<std::size_t>(i);
      field.data[idx] += dt * tendency.data[idx];
    }
  }

  return validation;
}

template <typename Real>
[[nodiscard]] TendencyApplyReport apply_tendency(
    const FieldView3D<Real> field,
    const FieldView3D<const Real> tendency,
    const Real dt) noexcept {
  const auto validation = validate_tendency_apply(field, tendency);
  if (validation.status != TendencyApplyStatus::ok) {
    return validation;
  }

  for (std::int32_t j = field.halo.j_lower; j < field.ny - field.halo.j_upper; ++j) {
    for (std::int32_t k = field.halo.k_lower; k < field.nz - field.halo.k_upper; ++k) {
      const auto plane =
          static_cast<std::size_t>(j) * static_cast<std::size_t>(field.stride_j) +
          static_cast<std::size_t>(k) * static_cast<std::size_t>(field.stride_k);
      for (std::int32_t i = field.halo.i_lower; i < field.nx - field.halo.i_upper; ++i) {
        const auto idx = plane + static_cast<std::size_t>(i);
        field.data[idx] += dt * tendency.data[idx];
      }
    }
  }

  return validation;
}

template <typename Real>
[[nodiscard]] TendencyApplyReport apply_zero_tendency(
    const FieldView2D<Real> field) noexcept {
  return validate_zero_tendency_apply(field);
}

template <typename Real>
[[nodiscard]] TendencyApplyReport apply_zero_tendency(
    const FieldView3D<Real> field) noexcept {
  return validate_zero_tendency_apply(field);
}

template <typename Real>
[[nodiscard]] TendencyApplyReport apply_state_tendencies(
    const StateView<Real> state,
    const StateView<const Real> tendency,
    const Real dt) noexcept {
  auto report = TendencyApplyReport{};
  report = merge_reports(report, validate_tendency_apply(state.u, tendency.u));
  report = merge_reports(report, validate_tendency_apply(state.v, tendency.v));
  report = merge_reports(report, validate_tendency_apply(state.w, tendency.w));
  report = merge_reports(report, validate_tendency_apply(state.ph, tendency.ph));
  report = merge_reports(report, validate_tendency_apply(state.phb, tendency.phb));
  report = merge_reports(report, validate_tendency_apply(state.t, tendency.t));
  report = merge_reports(report, validate_tendency_apply(state.p, tendency.p));
  report = merge_reports(report, validate_tendency_apply(state.pb, tendency.pb));
  report = merge_reports(report, validate_tendency_apply(state.qvapor, tendency.qvapor));
  report = merge_reports(report, validate_tendency_apply(state.qcloud, tendency.qcloud));
  report = merge_reports(report, validate_tendency_apply(state.qrain, tendency.qrain));
  report = merge_reports(report, validate_tendency_apply(state.qice, tendency.qice));
  report = merge_reports(report, validate_tendency_apply(state.qsnow, tendency.qsnow));
  report = merge_reports(report, validate_tendency_apply(state.qgraup, tendency.qgraup));
  report = merge_reports(report, validate_tendency_apply(state.qnice, tendency.qnice));
  report = merge_reports(report, validate_tendency_apply(state.qnrain, tendency.qnrain));
  report = merge_reports(report, validate_tendency_apply(state.mu, tendency.mu));
  report = merge_reports(report, validate_tendency_apply(state.mub, tendency.mub));
  report = merge_reports(report, validate_tendency_apply(state.psfc, tendency.psfc));
  report = merge_reports(report, validate_tendency_apply(state.u10, tendency.u10));
  report = merge_reports(report, validate_tendency_apply(state.v10, tendency.v10));
  report = merge_reports(report, validate_tendency_apply(state.t2, tendency.t2));
  report = merge_reports(report, validate_tendency_apply(state.q2, tendency.q2));
  report = merge_reports(report, validate_tendency_apply(state.rainc, tendency.rainc));
  report = merge_reports(report, validate_tendency_apply(state.rainnc, tendency.rainnc));
  if (report.status != TendencyApplyStatus::ok) {
    return report;
  }

  (void)apply_tendency(state.u, tendency.u, dt);
  (void)apply_tendency(state.v, tendency.v, dt);
  (void)apply_tendency(state.w, tendency.w, dt);
  (void)apply_tendency(state.ph, tendency.ph, dt);
  (void)apply_tendency(state.phb, tendency.phb, dt);
  (void)apply_tendency(state.t, tendency.t, dt);
  (void)apply_tendency(state.p, tendency.p, dt);
  (void)apply_tendency(state.pb, tendency.pb, dt);
  (void)apply_tendency(state.qvapor, tendency.qvapor, dt);
  (void)apply_tendency(state.qcloud, tendency.qcloud, dt);
  (void)apply_tendency(state.qrain, tendency.qrain, dt);
  (void)apply_tendency(state.qice, tendency.qice, dt);
  (void)apply_tendency(state.qsnow, tendency.qsnow, dt);
  (void)apply_tendency(state.qgraup, tendency.qgraup, dt);
  (void)apply_tendency(state.qnice, tendency.qnice, dt);
  (void)apply_tendency(state.qnrain, tendency.qnrain, dt);
  (void)apply_tendency(state.mu, tendency.mu, dt);
  (void)apply_tendency(state.mub, tendency.mub, dt);
  (void)apply_tendency(state.psfc, tendency.psfc, dt);
  (void)apply_tendency(state.u10, tendency.u10, dt);
  (void)apply_tendency(state.v10, tendency.v10, dt);
  (void)apply_tendency(state.t2, tendency.t2, dt);
  (void)apply_tendency(state.q2, tendency.q2, dt);
  (void)apply_tendency(state.rainc, tendency.rainc, dt);
  (void)apply_tendency(state.rainnc, tendency.rainnc, dt);
  return report;
}

template <typename Real>
[[nodiscard]] TendencyApplyReport apply_zero_state_tendency(
    const StateView<Real> state) noexcept {
  auto report = TendencyApplyReport{};
  report = merge_reports(report, validate_zero_tendency_apply(state.u));
  report = merge_reports(report, validate_zero_tendency_apply(state.v));
  report = merge_reports(report, validate_zero_tendency_apply(state.w));
  report = merge_reports(report, validate_zero_tendency_apply(state.ph));
  report = merge_reports(report, validate_zero_tendency_apply(state.phb));
  report = merge_reports(report, validate_zero_tendency_apply(state.t));
  report = merge_reports(report, validate_zero_tendency_apply(state.p));
  report = merge_reports(report, validate_zero_tendency_apply(state.pb));
  report = merge_reports(report, validate_zero_tendency_apply(state.qvapor));
  report = merge_reports(report, validate_zero_tendency_apply(state.qcloud));
  report = merge_reports(report, validate_zero_tendency_apply(state.qrain));
  report = merge_reports(report, validate_zero_tendency_apply(state.qice));
  report = merge_reports(report, validate_zero_tendency_apply(state.qsnow));
  report = merge_reports(report, validate_zero_tendency_apply(state.qgraup));
  report = merge_reports(report, validate_zero_tendency_apply(state.qnice));
  report = merge_reports(report, validate_zero_tendency_apply(state.qnrain));
  report = merge_reports(report, validate_zero_tendency_apply(state.mu));
  report = merge_reports(report, validate_zero_tendency_apply(state.mub));
  report = merge_reports(report, validate_zero_tendency_apply(state.psfc));
  report = merge_reports(report, validate_zero_tendency_apply(state.u10));
  report = merge_reports(report, validate_zero_tendency_apply(state.v10));
  report = merge_reports(report, validate_zero_tendency_apply(state.t2));
  report = merge_reports(report, validate_zero_tendency_apply(state.q2));
  report = merge_reports(report, validate_zero_tendency_apply(state.rainc));
  report = merge_reports(report, validate_zero_tendency_apply(state.rainnc));
  return report;
}

}  // namespace tywrf::dynamics
