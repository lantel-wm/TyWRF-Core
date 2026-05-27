#include "tywrf/dynamics/base_state.hpp"

#include "tywrf/state.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace {

constexpr float kSentinel = -9'999.0F;

tywrf::FieldLayout3D make_layout() {
  return tywrf::make_field_layout(
      tywrf::ActiveShape3D{3, 2, 2},
      tywrf::Halo3D{1, 1, 1, 1, 1, 1});
}

void fill_inputs(
    tywrf::FieldStorage3D<float>& pb,
    tywrf::FieldStorage3D<float>& t_init) {
  std::fill(pb.data(), pb.data() + pb.size(), kSentinel);
  std::fill(t_init.data(), t_init.data() + t_init.size(), kSentinel);

  const auto layout = pb.layout();
  auto pb_view = pb.view();
  auto t_view = t_init.view();
  for (std::int32_t j = layout.j_begin(); j < layout.j_end(); ++j) {
    for (std::int32_t k = layout.k_begin(); k < layout.k_end(); ++k) {
      for (std::int32_t i = layout.i_begin(); i < layout.i_end(); ++i) {
        pb_view(i, j, k) =
            90'000.0F - 5'000.0F * static_cast<float>(k - layout.k_begin()) -
            100.0F * static_cast<float>(i - layout.i_begin()) -
            50.0F * static_cast<float>(j - layout.j_begin());
        t_view(i, j, k) =
            1.0F + 0.25F * static_cast<float>(i - layout.i_begin()) +
            0.5F * static_cast<float>(j - layout.j_begin()) +
            0.125F * static_cast<float>(k - layout.k_begin());
      }
    }
  }
}

[[nodiscard]] float expected_alb(
    const float pb,
    const float t_init,
    const tywrf::dynamics::BaseStateReconstructionOptions options = {}) {
  const auto p0 = static_cast<double>(options.reference_pressure_pa);
  const auto cvpm = static_cast<double>(
      -(options.specific_heat_cp - options.dry_air_gas_constant) /
      options.specific_heat_cp);
  return static_cast<float>(
      (static_cast<double>(options.dry_air_gas_constant) / p0) *
      (static_cast<double>(t_init) +
       static_cast<double>(options.base_potential_temperature_k)) *
      std::pow(static_cast<double>(pb) / p0, cvpm));
}

template <typename Real>
[[nodiscard]] tywrf::FieldView3D<const Real> const_view(
    const tywrf::FieldView3D<Real>& field) noexcept {
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

[[nodiscard]] tywrf::dynamics::BaseStateReconstructionInputs make_inputs(
    tywrf::FieldStorage3D<float>& pb,
    tywrf::FieldStorage3D<float>& t_init,
    tywrf::FieldStorage3D<float>& alb) noexcept {
  return {
      static_cast<const tywrf::FieldStorage3D<float>&>(pb).view(),
      static_cast<const tywrf::FieldStorage3D<float>&>(t_init).view(),
      alb.view()};
}

void test_formula_reconstructs_wrf_inverse_base_density_alb() {
  const tywrf::dynamics::BaseStateReconstructionOptions default_options{};
  const auto expected_default_cvpm =
      -(default_options.specific_heat_cp -
        default_options.dry_air_gas_constant) /
      default_options.specific_heat_cp;
  assert(std::abs(default_options.cvpm - expected_default_cvpm) < 1.0e-7F);
  assert(std::abs(default_options.cvpm - -0.7142857F) < 1.0e-6F);

  tywrf::FieldStorage3D<float> pb(make_layout());
  tywrf::FieldStorage3D<float> t_init(make_layout());
  tywrf::FieldStorage3D<float> alb(make_layout());
  fill_inputs(pb, t_init);
  std::fill(alb.data(), alb.data() + alb.size(), kSentinel);

  const auto report = tywrf::dynamics::reconstruct_alb_from_pb_t_init(
      make_inputs(pb, t_init, alb));

  assert(report.ok());
  assert(report.requires_external_t_init_and_pb);
  assert(report.used_wrf_inverse_base_density_alb);
  assert(report.wrote_alb);
  assert(report.reconstructed_point_count == 12);

  const auto layout = pb.layout();
  const auto pb_view = pb.view();
  const auto t_view = t_init.view();
  const auto alb_view = alb.view();
  for (std::int32_t j = layout.j_begin(); j < layout.j_end(); ++j) {
    for (std::int32_t k = layout.k_begin(); k < layout.k_end(); ++k) {
      for (std::int32_t i = layout.i_begin(); i < layout.i_end(); ++i) {
        const auto expected = expected_alb(pb_view(i, j, k), t_view(i, j, k));
        assert(std::abs(alb_view(i, j, k) - expected) < 1.0e-7F);
      }
    }
  }
}

void test_bad_shape_returns_invalid_contract_without_writing() {
  tywrf::FieldStorage3D<float> pb(make_layout());
  tywrf::FieldStorage3D<float> t_init(
      tywrf::ActiveShape3D{4, 2, 2},
      tywrf::Halo3D{1, 1, 1, 1, 1, 1});
  tywrf::FieldStorage3D<float> alb(make_layout());
  fill_inputs(pb, t_init);
  std::fill(alb.data(), alb.data() + alb.size(), kSentinel);

  const auto report = tywrf::dynamics::reconstruct_alb_from_pb_t_init(
      make_inputs(pb, t_init, alb));

  assert(!report.ok());
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  for (std::size_t index = 0; index < alb.size(); ++index) {
    assert(alb.data()[index] == kSentinel);
  }
}

void test_nonfinite_or_nonpositive_input_returns_invalid_contract() {
  tywrf::FieldStorage3D<float> pb(make_layout());
  tywrf::FieldStorage3D<float> t_init(make_layout());
  tywrf::FieldStorage3D<float> alb(make_layout());
  fill_inputs(pb, t_init);
  std::fill(alb.data(), alb.data() + alb.size(), kSentinel);
  pb.view()(pb.layout().i_begin(), pb.layout().j_begin(), pb.layout().k_begin()) =
      std::numeric_limits<float>::quiet_NaN();

  const auto nonfinite_report = tywrf::dynamics::reconstruct_alb_from_pb_t_init(
      make_inputs(pb, t_init, alb));
  assert(!nonfinite_report.ok());
  assert(
      nonfinite_report.result.status ==
      tywrf::nest::NestStatus::invalid_contract);
  assert(nonfinite_report.invalid_point_count == 1);

  fill_inputs(pb, t_init);
  pb.view()(pb.layout().i_begin(), pb.layout().j_begin(), pb.layout().k_begin()) =
      0.0F;
  const auto pressure_report = tywrf::dynamics::reconstruct_alb_from_pb_t_init(
      make_inputs(pb, t_init, alb));
  assert(!pressure_report.ok());
  assert(pressure_report.result.status == tywrf::nest::NestStatus::invalid_contract);

  fill_inputs(pb, t_init);
  t_init.view()(t_init.layout().i_begin(), t_init.layout().j_begin(), t_init.layout().k_begin()) =
      std::numeric_limits<float>::infinity();
  const auto temperature_report =
      tywrf::dynamics::reconstruct_alb_from_pb_t_init(
          make_inputs(pb, t_init, alb));
  assert(!temperature_report.ok());
  assert(
      temperature_report.result.status ==
      tywrf::nest::NestStatus::invalid_contract);
}

void test_only_output_alb_active_cells_change() {
  tywrf::FieldStorage3D<float> pb(make_layout());
  tywrf::FieldStorage3D<float> t_init(make_layout());
  tywrf::FieldStorage3D<float> alb(make_layout());
  fill_inputs(pb, t_init);
  std::fill(alb.data(), alb.data() + alb.size(), kSentinel);
  const std::vector<float> pb_before(pb.data(), pb.data() + pb.size());
  const std::vector<float> t_before(t_init.data(), t_init.data() + t_init.size());

  const auto report = tywrf::dynamics::reconstruct_alb_from_pb_t_init(
      make_inputs(pb, t_init, alb));

  assert(report.ok());
  assert(std::equal(pb.data(), pb.data() + pb.size(), pb_before.begin()));
  assert(std::equal(t_init.data(), t_init.data() + t_init.size(), t_before.begin()));

  const auto layout = alb.layout();
  const auto alb_view = alb.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t k = 0; k < layout.nz; ++k) {
      for (std::int32_t i = 0; i < layout.nx; ++i) {
        const bool active =
            i >= layout.i_begin() && i < layout.i_end() &&
            j >= layout.j_begin() && j < layout.j_end() &&
            k >= layout.k_begin() && k < layout.k_end();
        assert(active || alb_view(i, j, k) == kSentinel);
      }
    }
  }
}

void test_requires_i_contiguous_canonical_views() {
  tywrf::FieldStorage3D<float> pb(make_layout());
  tywrf::FieldStorage3D<float> t_init(make_layout());
  tywrf::FieldStorage3D<float> alb(make_layout());
  fill_inputs(pb, t_init);
  std::fill(alb.data(), alb.data() + alb.size(), kSentinel);

  auto bad_pb = pb.view();
  bad_pb.stride_i = 2;
  const auto report = tywrf::dynamics::reconstruct_alb_from_pb_t_init(
      {const_view(bad_pb),
       static_cast<const tywrf::FieldStorage3D<float>&>(t_init).view(),
       alb.view()});

  assert(!report.ok());
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(bad_pb.index(1, 0, 0) != tywrf::canonical_index_3d(1, 0, 0, bad_pb.nx, bad_pb.nz));
  assert(pb.layout().stride_i == 1);
  assert(pb.layout().stride_k == pb.layout().nx);
  assert(pb.layout().stride_j == pb.layout().nx * pb.layout().nz);
}

void test_surface_albedo_and_restart_truth_are_not_inputs() {
  static_assert(std::is_standard_layout_v<tywrf::dynamics::BaseStateReconstructionInputs>);
  static_assert(std::is_trivially_copyable_v<tywrf::dynamics::BaseStateReconstructionInputs>);
  static_assert(std::is_standard_layout_v<tywrf::dynamics::BaseStateReconstructionOptions>);
  static_assert(std::is_trivially_copyable_v<tywrf::dynamics::BaseStateReconstructionOptions>);
  static_assert(std::is_standard_layout_v<tywrf::dynamics::BaseStateReconstructionReport>);
  static_assert(std::is_trivially_copyable_v<tywrf::dynamics::BaseStateReconstructionReport>);

  // Contract check: this helper only accepts PB/T_INIT/ALB mass views. It does
  // not consume surface ALBEDO and does not use later restart/history ALB truth.
  tywrf::FieldStorage3D<float> pb(make_layout());
  tywrf::FieldStorage3D<float> t_init(make_layout());
  tywrf::FieldStorage3D<float> alb(make_layout());
  fill_inputs(pb, t_init);
  std::fill(alb.data(), alb.data() + alb.size(), kSentinel);
  const auto report = tywrf::dynamics::reconstruct_alb_from_pb_t_init(
      make_inputs(pb, t_init, alb));
  assert(report.ok());
  assert(report.requires_external_t_init_and_pb);
}

}  // namespace

int main() {
  test_formula_reconstructs_wrf_inverse_base_density_alb();
  test_bad_shape_returns_invalid_contract_without_writing();
  test_nonfinite_or_nonpositive_input_returns_invalid_contract();
  test_only_output_alb_active_cells_change();
  test_requires_i_contiguous_canonical_views();
  test_surface_albedo_and_restart_truth_are_not_inputs();

  std::cout << "Validated narrow base-state ALB reconstruction helper\n";
  return 0;
}
