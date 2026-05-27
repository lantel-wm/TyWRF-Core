#include "tywrf/dynamics/base_state.hpp"

#include <cmath>

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

[[nodiscard]] bool valid_options(
    const BaseStateReconstructionOptions& options) noexcept {
  return std::isfinite(options.dry_air_gas_constant) &&
         std::isfinite(options.reference_pressure_pa) &&
         std::isfinite(options.base_potential_temperature_k) &&
         std::isfinite(options.specific_heat_cp) && std::isfinite(options.cvpm) &&
         options.dry_air_gas_constant > 0.0F &&
         options.reference_pressure_pa > 0.0F &&
         options.specific_heat_cp > 0.0F;
}

[[nodiscard]] std::uint64_t invalid_active_input_count(
    const BaseStateReconstructionInputs& inputs) noexcept {
  std::uint64_t invalid_count = 0;
  for (std::int32_t j = 0; j < active_ny(inputs.pb); ++j) {
    for (std::int32_t k = 0; k < active_nz(inputs.pb); ++k) {
      for (std::int32_t i = 0; i < active_nx(inputs.pb); ++i) {
        const auto pb = inputs.pb(
            inputs.pb.halo.i_lower + i,
            inputs.pb.halo.j_lower + j,
            inputs.pb.halo.k_lower + k);
        const auto t_init = inputs.t_init(
            inputs.t_init.halo.i_lower + i,
            inputs.t_init.halo.j_lower + j,
            inputs.t_init.halo.k_lower + k);
        if (!std::isfinite(pb) || !std::isfinite(t_init) || pb <= 0.0F) {
          ++invalid_count;
        }
      }
    }
  }
  return invalid_count;
}

}  // namespace

BaseStateReconstructionReport reconstruct_alb_from_pb_t_init(
    const BaseStateReconstructionInputs& inputs,
    const BaseStateReconstructionOptions options) noexcept {
  BaseStateReconstructionReport report{};
  report.active_nx = active_nx(inputs.pb);
  report.active_ny = active_ny(inputs.pb);
  report.active_nz = active_nz(inputs.pb);

  if (!valid_canonical_view(inputs.pb) ||
      !valid_canonical_view(inputs.t_init) ||
      !valid_canonical_view(inputs.alb)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "base-state ALB reconstruction requires non-null canonical mass views");
    return report;
  }

  if (!same_active_shape(inputs.pb, inputs.t_init) ||
      !same_active_shape(inputs.pb, inputs.alb)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "base-state ALB reconstruction views must share active shape");
    return report;
  }

  if (!valid_options(options)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "base-state ALB reconstruction options must be finite and positive");
    return report;
  }

  report.invalid_point_count = invalid_active_input_count(inputs);
  if (report.invalid_point_count != 0) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "base-state ALB reconstruction requires finite T_INIT and positive finite PB");
    return report;
  }

  const auto p0 = static_cast<double>(options.reference_pressure_pa);
  const auto rd_over_p0 =
      static_cast<double>(options.dry_air_gas_constant) / p0;
  const auto t0 = static_cast<double>(options.base_potential_temperature_k);
  const auto cvpm = static_cast<double>(options.cvpm);

  for (std::int32_t j = 0; j < active_ny(inputs.pb); ++j) {
    for (std::int32_t k = 0; k < active_nz(inputs.pb); ++k) {
      for (std::int32_t i = 0; i < active_nx(inputs.pb); ++i) {
        const auto pb_i = inputs.pb.halo.i_lower + i;
        const auto pb_j = inputs.pb.halo.j_lower + j;
        const auto pb_k = inputs.pb.halo.k_lower + k;
        const auto t_i = inputs.t_init.halo.i_lower + i;
        const auto t_j = inputs.t_init.halo.j_lower + j;
        const auto t_k = inputs.t_init.halo.k_lower + k;
        const auto alb_i = inputs.alb.halo.i_lower + i;
        const auto alb_j = inputs.alb.halo.j_lower + j;
        const auto alb_k = inputs.alb.halo.k_lower + k;

        const auto theta_base =
            static_cast<double>(inputs.t_init(t_i, t_j, t_k)) + t0;
        const auto pressure_ratio =
            static_cast<double>(inputs.pb(pb_i, pb_j, pb_k)) / p0;
        inputs.alb(alb_i, alb_j, alb_k) = static_cast<float>(
            rd_over_p0 * theta_base * std::pow(pressure_ratio, cvpm));
        ++report.reconstructed_point_count;
      }
    }
  }

  report.result = ok_result();
  report.used_wrf_inverse_base_density_alb = true;
  report.wrote_alb = true;
  return report;
}

}  // namespace tywrf::dynamics
