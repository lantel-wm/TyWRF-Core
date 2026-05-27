#include "tywrf/dynamics/base_state.hpp"

#include <cmath>
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

[[nodiscard]] constexpr bool has_optional_alb(
    const FieldView3D<float>& alb) noexcept {
  return alb.data != nullptr;
}

[[nodiscard]] constexpr bool has_optional_phb(
    const FieldView3D<float>& phb) noexcept {
  return phb.data != nullptr;
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

[[nodiscard]] double surface_pressure_pa(
    const double terrain_height_m,
    const KrosaMassBaseStateReconstructionOptions& options) noexcept {
  const auto p00 = static_cast<double>(options.sea_level_base_pressure_pa);
  const auto t00 = static_cast<double>(options.sea_level_base_temperature_k);
  const auto lapse = static_cast<double>(options.base_lapse_k);
  const auto gravity = static_cast<double>(options.gravity);
  const auto rd = static_cast<double>(options.dry_air_gas_constant);
  const auto t00_over_lapse = t00 / lapse;
  const auto discriminant =
      t00_over_lapse * t00_over_lapse -
      (2.0 * gravity * terrain_height_m) / (lapse * rd);
  if (discriminant < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return p00 * std::exp(-t00_over_lapse + std::sqrt(discriminant));
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

[[nodiscard]] BaseStateReconstructionOptions alb_options(
    const KrosaMassBaseStateReconstructionOptions& options) noexcept {
  return {
      options.dry_air_gas_constant,
      options.reference_pressure_pa,
      options.base_potential_temperature_k,
      options.specific_heat_cp,
      options.cvpm};
}

[[nodiscard]] bool valid_phb_column_pressures(
    const double mub,
    const KrosaMassBaseStateReconstructionInputs& inputs,
    const KrosaMassBaseStateReconstructionOptions& options,
    const std::int32_t mass_k) noexcept {
  const auto p_top = static_cast<double>(inputs.p_top_pa);
  const auto pfu =
      static_cast<double>(inputs.c3f.values[mass_k + 1]) * mub +
      static_cast<double>(inputs.c4f.values[mass_k + 1]) + p_top;
  const auto pfd = static_cast<double>(inputs.c3f.values[mass_k]) * mub +
                   static_cast<double>(inputs.c4f.values[mass_k]) + p_top;
  const auto phm = static_cast<double>(inputs.c3h.values[mass_k]) * mub +
                   static_cast<double>(inputs.c4h.values[mass_k]) + p_top;
  return std::isfinite(pfu) && std::isfinite(pfd) && std::isfinite(phm) &&
         pfu > 0.0 && pfd > 0.0 && phm > 0.0 && pfd > pfu &&
         std::isfinite(std::log(pfd / pfu)) &&
         std::isfinite(static_cast<double>(options.gravity));
}

[[nodiscard]] bool reconstruct_phb_column(
    const std::int32_t i,
    const std::int32_t j,
    const double terrain_height_m,
    const double mub,
    const KrosaMassBaseStateReconstructionInputs& inputs,
    const KrosaMassBaseStateReconstructionOutputs& outputs,
    const KrosaMassBaseStateReconstructionOptions& options) noexcept {
  const auto phb_i = outputs.phb.halo.i_lower + i;
  const auto phb_j = outputs.phb.halo.j_lower + j;
  outputs.phb(phb_i, phb_j, outputs.phb.halo.k_lower) =
      static_cast<float>(terrain_height_m * static_cast<double>(options.gravity));

  const auto p_top = static_cast<double>(inputs.p_top_pa);
  for (std::int32_t k = 0; k < active_nz(outputs.pb); ++k) {
    const auto pfu =
        static_cast<double>(inputs.c3f.values[k + 1]) * mub +
        static_cast<double>(inputs.c4f.values[k + 1]) + p_top;
    const auto pfd = static_cast<double>(inputs.c3f.values[k]) * mub +
                     static_cast<double>(inputs.c4f.values[k]) + p_top;
    const auto phm = static_cast<double>(inputs.c3h.values[k]) * mub +
                     static_cast<double>(inputs.c4h.values[k]) + p_top;
    const auto alb = static_cast<double>(outputs.alb(
        outputs.alb.halo.i_lower + i,
        outputs.alb.halo.j_lower + j,
        outputs.alb.halo.k_lower + k));
    if (!std::isfinite(alb)) {
      return false;
    }

    const auto previous_phb = static_cast<double>(outputs.phb(
        phb_i,
        phb_j,
        outputs.phb.halo.k_lower + k));
    const auto next_phb = previous_phb + alb * phm * std::log(pfd / pfu);
    if (!std::isfinite(next_phb)) {
      return false;
    }
    outputs.phb(phb_i, phb_j, outputs.phb.halo.k_lower + k + 1) =
        static_cast<float>(next_phb);
  }
  return true;
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

KrosaMassBaseStateReconstructionReport reconstruct_krosa_mass_base_state(
    const KrosaMassBaseStateReconstructionInputs& inputs,
    const KrosaMassBaseStateReconstructionOutputs& outputs,
    const KrosaMassBaseStateReconstructionOptions options) noexcept {
  KrosaMassBaseStateReconstructionReport report{};
  report.active_nx = active_nx(outputs.pb);
  report.active_ny = active_ny(outputs.pb);
  report.active_nz = active_nz(outputs.pb);

  if (!valid_canonical_view(inputs.terrain_height_m) ||
      !valid_canonical_view(outputs.pb) ||
      !valid_canonical_view(outputs.t_init) ||
      !valid_canonical_view(outputs.mub)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "KROSA base-state reconstruction requires canonical terrain and output views");
    return report;
  }

  if (!same_active_shape(outputs.pb, outputs.t_init) ||
      !same_horizontal_active_shape(inputs.terrain_height_m, outputs.pb) ||
      !same_active_shape(inputs.terrain_height_m, outputs.mub)) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "KROSA base-state reconstruction views must share active mass shape");
    return report;
  }

  if (has_optional_alb(outputs.alb) &&
      (!valid_canonical_view(outputs.alb) ||
       !same_active_shape(outputs.pb, outputs.alb))) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "optional KROSA ALB output must be a canonical mass view");
    return report;
  }

  if (has_optional_phb(outputs.phb)) {
    if (!valid_canonical_view(outputs.phb) ||
        !same_horizontal_active_shape(outputs.pb, outputs.phb) ||
        active_nz(outputs.phb) != active_nz(outputs.pb) + 1) {
      report.result = result(
          nest::NestStatus::invalid_contract,
          "optional KROSA PHB output must be a canonical full-level view");
      return report;
    }
    if (!has_optional_alb(outputs.alb)) {
      report.result = result(
          nest::NestStatus::invalid_contract,
          "optional KROSA PHB reconstruction requires an ALB output view");
      return report;
    }
  }

  if (!valid_coefficients(inputs.c3h, active_nz(outputs.pb)) ||
      !valid_coefficients(inputs.c4h, active_nz(outputs.pb))) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "KROSA base-state C3H/C4H counts must match mass levels");
    return report;
  }

  if (has_optional_phb(outputs.phb) &&
      (!valid_coefficients(inputs.c3f, active_nz(outputs.pb) + 1) ||
       !valid_coefficients(inputs.c4f, active_nz(outputs.pb) + 1))) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "KROSA base-state C3F/C4F counts must match full levels");
    return report;
  }

  if (!valid_options(options) || !std::isfinite(inputs.p_top_pa) ||
      inputs.p_top_pa < 0.0F) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "KROSA base-state constants and P_TOP must be finite and valid");
    return report;
  }

  if (!finite_coefficients(inputs.c3h) || !finite_coefficients(inputs.c4h) ||
      (has_optional_phb(outputs.phb) &&
       (!finite_coefficients(inputs.c3f) ||
        !finite_coefficients(inputs.c4f)))) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "KROSA base-state vertical coefficients must be finite");
    return report;
  }

  const auto p_top = static_cast<double>(inputs.p_top_pa);
  const auto p00 = static_cast<double>(options.sea_level_base_pressure_pa);
  const auto rd_over_cp =
      static_cast<double>(options.dry_air_gas_constant) /
      static_cast<double>(options.specific_heat_cp);
  const auto theta_base =
      static_cast<double>(options.base_potential_temperature_k);
  const auto p_strat =
      static_cast<double>(options.stratosphere_base_pressure_pa);

  for (std::int32_t j = 0; j < active_ny(outputs.pb); ++j) {
    for (std::int32_t i = 0; i < active_nx(outputs.pb); ++i) {
      const auto h_i = inputs.terrain_height_m.halo.i_lower + i;
      const auto h_j = inputs.terrain_height_m.halo.j_lower + j;
      const auto terrain =
          static_cast<double>(inputs.terrain_height_m(h_i, h_j));
      const auto p_surf = surface_pressure_pa(terrain, options);
      if (!std::isfinite(terrain) || !std::isfinite(p_surf) ||
          p_surf <= p_top) {
        ++report.invalid_column_count;
        report.invalid_point_count +=
            static_cast<std::uint64_t>(active_nz(outputs.pb));
        continue;
      }

      const auto mub = p_surf - p_top;
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
        if (has_optional_phb(outputs.phb) &&
            !valid_phb_column_pressures(mub, inputs, options, k)) {
          ++report.invalid_point_count;
        }
      }
    }
  }

  if (report.invalid_column_count != 0 || report.invalid_point_count != 0) {
    report.result = result(
        nest::NestStatus::invalid_contract,
        "KROSA base-state reconstruction inputs produce invalid pressure");
    return report;
  }

  if (report.unsupported_stratosphere_point_count != 0) {
    report.result = result(
        nest::NestStatus::not_implemented,
        "KROSA base-state stratosphere pressure branch is not implemented");
    return report;
  }

  for (std::int32_t j = 0; j < active_ny(outputs.pb); ++j) {
    for (std::int32_t i = 0; i < active_nx(outputs.pb); ++i) {
      const auto h_i = inputs.terrain_height_m.halo.i_lower + i;
      const auto h_j = inputs.terrain_height_m.halo.j_lower + j;
      const auto p_surf = surface_pressure_pa(
          static_cast<double>(inputs.terrain_height_m(h_i, h_j)),
          options);
      const auto mub = p_surf - p_top;
      outputs.mub(
          outputs.mub.halo.i_lower + i,
          outputs.mub.halo.j_lower + j) = static_cast<float>(mub);
      ++report.reconstructed_column_count;

      for (std::int32_t k = 0; k < active_nz(outputs.pb); ++k) {
        const auto pb_i = outputs.pb.halo.i_lower + i;
        const auto pb_j = outputs.pb.halo.j_lower + j;
        const auto pb_k = outputs.pb.halo.k_lower + k;
        const auto t_i = outputs.t_init.halo.i_lower + i;
        const auto t_j = outputs.t_init.halo.j_lower + j;
        const auto t_k = outputs.t_init.halo.k_lower + k;

        const auto pb =
            static_cast<double>(inputs.c3h.values[k]) * mub +
            static_cast<double>(inputs.c4h.values[k]) + p_top;
        const auto temperature = wrf_base_temperature_k(pb, options);
        const auto t_init =
            temperature * std::pow(p00 / pb, rd_over_cp) - theta_base;

        outputs.pb(pb_i, pb_j, pb_k) = static_cast<float>(pb);
        outputs.t_init(t_i, t_j, t_k) = static_cast<float>(t_init);
        ++report.reconstructed_point_count;
      }
    }
  }

  report.used_wrf_surface_pressure_formula = true;
  report.used_wrf_hybrid_mass_coefficients = true;
  report.used_wrf_t_init_formula = true;
  report.wrote_pb = true;
  report.wrote_t_init = true;
  report.wrote_mub = true;

  if (has_optional_alb(outputs.alb)) {
    const auto alb_report = reconstruct_alb_from_pb_t_init(
        {const_view(outputs.pb), const_view(outputs.t_init), outputs.alb},
        alb_options(options));
    if (!alb_report.ok()) {
      report.result = alb_report.result;
      report.invalid_point_count = alb_report.invalid_point_count;
      return report;
    }
    report.reused_alb_helper = true;
    report.wrote_alb = alb_report.wrote_alb;
  }

  if (has_optional_phb(outputs.phb)) {
    for (std::int32_t j = 0; j < active_ny(outputs.pb); ++j) {
      for (std::int32_t i = 0; i < active_nx(outputs.pb); ++i) {
        const auto h_i = inputs.terrain_height_m.halo.i_lower + i;
        const auto h_j = inputs.terrain_height_m.halo.j_lower + j;
        const auto terrain =
            static_cast<double>(inputs.terrain_height_m(h_i, h_j));
        const auto mub = static_cast<double>(outputs.mub(
            outputs.mub.halo.i_lower + i,
            outputs.mub.halo.j_lower + j));
        if (!reconstruct_phb_column(i, j, terrain, mub, inputs, outputs, options)) {
          report.result = result(
              nest::NestStatus::invalid_contract,
              "KROSA PHB reconstruction requires finite log-linear inputs");
          return report;
        }
      }
    }
    report.wrote_phb = true;
    report.phb_full_level_reconstruction_implemented = true;
  }

  report.result = ok_result();
  return report;
}

}  // namespace tywrf::dynamics
