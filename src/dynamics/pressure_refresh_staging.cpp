#include "tywrf/dynamics/pressure_refresh_staging.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

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
    const FieldView3D<RhsReal>& rhs) noexcept {
  return active_nx(lhs) == active_nx(rhs) && active_ny(lhs) == active_ny(rhs);
}

template <typename Real>
[[nodiscard]] constexpr FieldView3D<const Real> const_view(
    const FieldView3D<Real>& field) noexcept {
  return {
      field.data,
      field.nx,
      field.ny,
      field.nz,
      field.stride_i,
      field.stride_k,
      field.stride_j,
      field.halo};
}

template <typename Real>
[[nodiscard]] constexpr FieldView2D<const Real> const_view(
    const FieldView2D<Real>& field) noexcept {
  return {
      field.data,
      field.nx,
      field.ny,
      field.stride_i,
      field.stride_j,
      field.halo};
}

[[nodiscard]] std::int32_t vector_size_or_negative(
    const std::size_t size) noexcept {
  if (size > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    return -1;
  }
  return static_cast<std::int32_t>(size);
}

[[nodiscard]] VerticalCoefficientView coefficient_view(
    const std::vector<float>& values) noexcept {
  return {values.data(), vector_size_or_negative(values.size())};
}

[[nodiscard]] PressureRefreshStagingReport base_report(
    const StateView<float>& state,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept {
  PressureRefreshStagingReport report{};
  report.expected_mass_level_count = active_nz(state.p);
  report.expected_full_level_count = report.expected_mass_level_count + 1;
  report.c3f_count = vector_size_or_negative(metadata.c3f.size());
  report.c4f_count = vector_size_or_negative(metadata.c4f.size());
  report.c3h_count = vector_size_or_negative(metadata.c3h.size());
  report.c4h_count = vector_size_or_negative(metadata.c4h.size());
  report.alb_is_external = true;
  return report;
}

}  // namespace

PressureRefreshStagingReport validate_krosa_pressure_refresh_staging(
    const StateView<float>& state,
    const FieldView3D<const float> external_alb,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept {
  auto report = base_report(state, metadata);

  if (!valid_canonical_view(state.p) || !valid_canonical_view(state.pb) ||
      !valid_canonical_view(state.t) || !valid_canonical_view(state.ph) ||
      !valid_canonical_view(state.phb) || !valid_canonical_view(state.mu) ||
      !valid_canonical_view(state.mub)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure refresh staging requires non-null canonical state views");
    return report;
  }

  if (!same_active_shape(state.p, state.pb) ||
      !same_active_shape(state.p, state.t)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure refresh staging state mass views must share active shape");
    return report;
  }

  if (!same_horizontal_active_shape(state.p, state.ph) ||
      !same_horizontal_active_shape(state.p, state.phb) ||
      active_nz(state.ph) != active_nz(state.p) + 1 ||
      active_nz(state.phb) != active_nz(state.p) + 1) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure refresh staging PH/PHB must have one extra full level");
    return report;
  }

  if (!same_horizontal_active_shape(state.mu, state.p) ||
      !same_horizontal_active_shape(state.mub, state.p)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure refresh staging MU/MUB must match mass horizontal shape");
    return report;
  }

  if (!valid_canonical_view(external_alb)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure refresh staging requires an external canonical ALB view");
    return report;
  }

  if (!same_active_shape(state.p, external_alb)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "external ALB active shape must match state mass shape");
    return report;
  }

  if (report.c3f_count != report.expected_full_level_count ||
      report.c4f_count != report.expected_full_level_count ||
      report.c3h_count != report.expected_mass_level_count ||
      report.c4h_count != report.expected_mass_level_count) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure refresh metadata coefficient counts do not match state levels");
    return report;
  }

  if (metadata.c3f.data() == nullptr || metadata.c4f.data() == nullptr ||
      metadata.c3h.data() == nullptr || metadata.c4h.data() == nullptr) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure refresh metadata coefficient buffers are missing");
    return report;
  }

  if (!std::isfinite(metadata.p_top_pa) || metadata.p_top_pa < 0.0F ||
      metadata.p_top_source == io::PressureRefreshPTopSource::missing) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "pressure refresh metadata P_TOP must be present, finite, and non-negative");
    return report;
  }

  report.result = ok_result();
  return report;
}

PressureRefreshStagingReport validate_krosa_pressure_refresh_staging(
    State<float>& state,
    const FieldStorage3D<float>& external_alb,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept {
  return validate_krosa_pressure_refresh_staging(
      state.view(),
      external_alb.view(),
      metadata);
}

PressureRefreshStagingResult make_krosa_pressure_refresh_inputs(
    const StateView<float> state,
    const FieldView3D<const float> external_alb,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept {
  PressureRefreshStagingResult staging{};
  staging.report =
      validate_krosa_pressure_refresh_staging(state, external_alb, metadata);
  if (!staging.ok()) {
    return staging;
  }

  staging.inputs = {
      state.p,
      const_view(state.pb),
      const_view(state.t),
      external_alb,
      const_view(state.ph),
      const_view(state.phb),
      const_view(state.mu),
      const_view(state.mub),
      coefficient_view(metadata.c3f),
      coefficient_view(metadata.c4f),
      coefficient_view(metadata.c3h),
      coefficient_view(metadata.c4h),
      metadata.p_top_pa,
  };
  return staging;
}

PressureRefreshStagingResult make_krosa_pressure_refresh_inputs(
    State<float>& state,
    const FieldStorage3D<float>& external_alb,
    const io::KrosaPressureRefreshMetadata& metadata) noexcept {
  return make_krosa_pressure_refresh_inputs(
      state.view(),
      external_alb.view(),
      metadata);
}

}  // namespace tywrf::dynamics
