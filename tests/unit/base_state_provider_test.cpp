#include "tywrf/dynamics/base_state_provider.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

constexpr float kSentinel = -9'999.0F;

tywrf::Grid make_grid() {
  return tywrf::Grid({4, 3, 2, 3, tywrf::uniform_halo_3d(1)});
}

std::vector<float> make_terrain() {
  std::vector<float> terrain(12, 0.0F);
  for (std::int32_t j = 0; j < 3; ++j) {
    for (std::int32_t i = 0; i < 4; ++i) {
      terrain[static_cast<std::size_t>(j) * 4U + static_cast<std::size_t>(i)] =
          25.0F + 10.0F * static_cast<float>(j) + static_cast<float>(i);
    }
  }
  return terrain;
}

tywrf::io::KrosaPressureRefreshMetadata make_metadata() {
  tywrf::io::KrosaPressureRefreshMetadata metadata{};
  metadata.time_index = 0;
  metadata.p_top_pa = 5'000.0F;
  metadata.p_top_source = tywrf::io::PressureRefreshPTopSource::global_attribute;
  metadata.c3f = {1.0F, 0.5F, 0.0F};
  metadata.c4f = {0.0F, 0.0F, 0.0F};
  metadata.c3h = {0.72F, 0.28F};
  metadata.c4h = {750.0F, 125.0F};
  metadata.terrain_source_name = "HGT";
  metadata.terrain_nx = 4;
  metadata.terrain_ny = 3;
  metadata.terrain_height_m = make_terrain();
  return metadata;
}

tywrf::FieldStorage2D<float> make_override_terrain(
    const tywrf::ActiveShape2D active,
    const float offset_m) {
  tywrf::FieldStorage2D<float> terrain(
      active,
      tywrf::horizontal_halo(make_grid().config().halo));
  auto view = terrain.view();
  for (std::int32_t j = 0; j < active.ny; ++j) {
    for (std::int32_t i = 0; i < active.nx; ++i) {
      view(view.halo.i_lower + i, view.halo.j_lower + j) =
          25.0F + offset_m + 10.0F * static_cast<float>(j) +
          static_cast<float>(i);
    }
  }
  return terrain;
}

tywrf::FieldStorage2D<float> make_override_terrain(const float offset_m) {
  return make_override_terrain({4, 3}, offset_m);
}

tywrf::dynamics::KrosaBaseStateProviderTerrainOverride moved_terrain_override(
    const tywrf::FieldStorage2D<float>& terrain) {
  return {
      terrain.view(),
      "moved_candidate_HGT",
      "override:moved_candidate_HGT"};
}

template <typename Real>
std::int32_t active_nx(const tywrf::FieldView3D<Real>& view) {
  return view.nx - view.halo.i_lower - view.halo.i_upper;
}

template <typename Real>
std::int32_t active_ny(const tywrf::FieldView3D<Real>& view) {
  return view.ny - view.halo.j_lower - view.halo.j_upper;
}

template <typename Real>
std::int32_t active_nz(const tywrf::FieldView3D<Real>& view) {
  return view.nz - view.halo.k_lower - view.halo.k_upper;
}

template <typename Real>
std::int32_t active_nx(const tywrf::FieldView2D<Real>& view) {
  return view.nx - view.halo.i_lower - view.halo.i_upper;
}

template <typename Real>
std::int32_t active_ny(const tywrf::FieldView2D<Real>& view) {
  return view.ny - view.halo.j_lower - view.halo.j_upper;
}

template <typename Storage>
std::vector<float> snapshot(const Storage& storage) {
  return std::vector<float>(storage.data(), storage.data() + storage.size());
}

template <typename Storage>
void fill_storage(Storage& storage, const float value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

template <typename Storage>
void assert_storage_unchanged(
    const Storage& storage,
    const std::vector<float>& before,
    const std::string_view label) {
  if (before.size() != storage.size() ||
      !std::equal(before.begin(), before.end(), storage.data())) {
    std::cerr << label << " changed unexpectedly\n";
    assert(false);
  }
}

tywrf::nest::RemapWindow make_mass_overlap_window(
    const std::int32_t i_begin,
    const std::int32_t j_begin,
    const std::int32_t extent_i,
    const std::int32_t extent_j) {
  return {
      tywrf::nest::HorizontalStagger::mass,
      0,
      0,
      i_begin,
      j_begin,
      extent_i,
      extent_j,
      0,
      0};
}

bool horizontally_exposed(
    const tywrf::nest::RemapWindow& overlap,
    const std::int32_t active_i,
    const std::int32_t active_j) {
  const bool inside_overlap_i =
      active_i >= overlap.new_i_begin &&
      active_i < overlap.new_i_begin + overlap.extent_i;
  const bool inside_overlap_j =
      active_j >= overlap.new_j_begin &&
      active_j < overlap.new_j_begin + overlap.extent_j;
  return !(inside_overlap_i && inside_overlap_j);
}

void fill_interpolated_mub(tywrf::FieldStorage2D<float>& mub) {
  fill_storage(mub, kSentinel);
  const auto layout = mub.layout();
  auto view = mub.view();
  for (std::int32_t j = layout.j_begin(); j < layout.j_end(); ++j) {
    for (std::int32_t i = layout.i_begin(); i < layout.i_end(); ++i) {
      view(i, j) =
          85'000.0F + 100.0F * static_cast<float>(i - layout.i_begin()) +
          25.0F * static_cast<float>(j - layout.j_begin());
    }
  }
}

float expected_t_init_from_pb(
    const double pb,
    const tywrf::dynamics::KrosaMassBaseStateReconstructionOptions options =
        {}) {
  const auto p00 = static_cast<double>(options.sea_level_base_pressure_pa);
  const auto troposphere_temperature =
      static_cast<double>(options.sea_level_base_temperature_k) +
      static_cast<double>(options.base_lapse_k) * std::log(pb / p00);
  const auto base_temperature =
      troposphere_temperature >
              static_cast<double>(options.isothermal_temperature_k)
          ? troposphere_temperature
          : static_cast<double>(options.isothermal_temperature_k);
  return static_cast<float>(
      base_temperature *
          std::pow(
              p00 / pb,
              static_cast<double>(options.dry_air_gas_constant) /
                  static_cast<double>(options.specific_heat_cp)) -
      static_cast<double>(options.base_potential_temperature_k));
}

float expected_alb_from_pb_t_init(
    const double pb,
    const double t_init,
    const tywrf::dynamics::KrosaMassBaseStateReconstructionOptions options =
        {}) {
  const auto p0 = static_cast<double>(options.reference_pressure_pa);
  return static_cast<float>(
      (static_cast<double>(options.dry_air_gas_constant) / p0) *
      (t_init + static_cast<double>(options.base_potential_temperature_k)) *
      std::pow(pb / p0, static_cast<double>(options.cvpm)));
}

void assert_exposed_recompute_flags(
    const tywrf::dynamics::KrosaExposedBaseStateRecomputeReport& report) {
  assert(
      std::string_view(report.source) ==
      "exposed_mub_base_state_recompute_provider");
  assert(
      std::string_view(report.disposition) ==
      "provider_buffer_only_no_gate_no_selected_field_numerics");
  assert(report.diagnostic_only);
  assert(!report.gate_candidate);
  assert(!report.integrator_output);
  assert(!report.selected_field_numerics_enabled);
  assert(!report.calls_pressure_refresh_compute);
  assert(!report.reads_hgt);
  assert(!report.touched_phb);
  assert(!report.wrote_state_t);
  assert(!report.wrote_state_alb);
}

void assert_report_flags(const tywrf::dynamics::KrosaBaseStateProviderReport& report) {
  assert(report.source == std::string_view("base_state_reconstruction_provider"));
  assert(report.diagnostic_only);
  assert(!report.gate_candidate);
  assert(!report.integrator_output);
  assert(!report.calls_pressure_refresh_compute);
}

void assert_views_empty(const tywrf::dynamics::KrosaBaseStateProvider& provider) {
  const auto views = provider.views();
  assert(views.pb.nx == 0);
  assert(views.t_init.nx == 0);
  assert(views.mub.nx == 0);
  assert(views.alb.nx == 0);
  assert(views.phb.nx == 0);
}

void test_success_reconstructs_and_holds_views() {
  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report = provider.reconstruct(make_grid(), make_metadata());
  assert(report.ok());
  assert_report_flags(report);
  assert(report.allocated_buffers);
  assert(report.wrote_pb);
  assert(report.wrote_t_init);
  assert(report.wrote_mub);
  assert(report.wrote_alb);
  assert(report.wrote_phb);
  assert(report.reconstruction.reused_alb_helper);
  assert(report.reconstruction.phb_full_level_reconstruction_implemented);
  assert(report.active_nx == 4);
  assert(report.active_ny == 3);
  assert(report.mass_nz == 2);
  assert(report.full_nz == 3);
  assert(report.expected_terrain_point_count == 12);
  assert(report.terrain_point_count == 12);
  assert(!report.terrain_override_used);
  assert(report.terrain_source_name == "HGT");
  assert(report.terrain_provenance == "metadata:HGT");

  const auto views = provider.views();
  assert(views.pb.data != nullptr);
  assert(views.t_init.data != nullptr);
  assert(views.mub.data != nullptr);
  assert(views.alb.data != nullptr);
  assert(views.phb.data != nullptr);
  assert(active_nx(views.pb) == 4);
  assert(active_ny(views.pb) == 3);
  assert(active_nz(views.pb) == 2);
  assert(active_nx(views.t_init) == 4);
  assert(active_ny(views.t_init) == 3);
  assert(active_nz(views.t_init) == 2);
  assert(active_nx(views.alb) == 4);
  assert(active_ny(views.alb) == 3);
  assert(active_nz(views.alb) == 2);
  assert(active_nx(views.phb) == 4);
  assert(active_ny(views.phb) == 3);
  assert(active_nz(views.phb) == 3);
  assert(active_nx(views.mub) == 4);
  assert(active_ny(views.mub) == 3);

  const auto i = views.pb.halo.i_lower;
  const auto j = views.pb.halo.j_lower;
  const auto k = views.pb.halo.k_lower;
  assert(std::isfinite(views.pb(i, j, k)));
  assert(std::isfinite(views.t_init(i, j, k)));
  assert(std::isfinite(views.alb(i, j, k)));
  assert(std::isfinite(views.mub(views.mub.halo.i_lower, views.mub.halo.j_lower)));
  assert(std::isfinite(
      views.phb(views.phb.halo.i_lower, views.phb.halo.j_lower, views.phb.halo.k_lower)));
}

void test_override_terrain_reconstructs_and_changes_output() {
  tywrf::dynamics::KrosaBaseStateProvider metadata_provider;
  const auto metadata_report =
      metadata_provider.reconstruct(make_grid(), make_metadata());
  assert(metadata_report.ok());
  const auto metadata_views = metadata_provider.views();
  const auto i = metadata_views.mub.halo.i_lower;
  const auto j = metadata_views.mub.halo.j_lower;
  const auto metadata_mub = metadata_views.mub(i, j);
  const auto metadata_phb_surface = metadata_views.phb(
      metadata_views.phb.halo.i_lower,
      metadata_views.phb.halo.j_lower,
      metadata_views.phb.halo.k_lower);

  const auto moved_terrain = make_override_terrain(150.0F);
  tywrf::dynamics::KrosaBaseStateProvider override_provider;
  const auto report = override_provider.reconstruct(
      make_grid(),
      make_metadata(),
      moved_terrain_override(moved_terrain));
  assert(report.ok());
  assert_report_flags(report);
  assert(report.allocated_buffers);
  assert(report.terrain_override_used);
  assert(report.terrain_source_name == "moved_candidate_HGT");
  assert(report.terrain_provenance == "override:moved_candidate_HGT");
  assert(report.terrain_nx == 4);
  assert(report.terrain_ny == 3);
  assert(report.expected_terrain_point_count == 12);
  assert(report.terrain_point_count == 12);

  const auto override_views = override_provider.views();
  const auto override_mub = override_views.mub(
      override_views.mub.halo.i_lower,
      override_views.mub.halo.j_lower);
  const auto override_phb_surface = override_views.phb(
      override_views.phb.halo.i_lower,
      override_views.phb.halo.j_lower,
      override_views.phb.halo.k_lower);
  assert(std::isfinite(override_mub));
  assert(std::fabs(override_mub - metadata_mub) > 1.0F);
  assert(std::fabs(override_phb_surface - metadata_phb_surface) > 100.0F);
}

void test_missing_terrain_is_rejected_before_allocation() {
  auto metadata = make_metadata();
  metadata.terrain_nx = 0;
  metadata.terrain_ny = 0;
  metadata.terrain_height_m.clear();

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report = provider.reconstruct(make_grid(), metadata);
  assert(!report.ok());
  assert_report_flags(report);
  assert(!report.allocated_buffers);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(std::string_view(report.result.message).find("terrain") != std::string_view::npos);
}

void test_bad_override_shape_is_rejected_before_allocation() {
  const auto bad_terrain = make_override_terrain({5, 3}, 0.0F);

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report = provider.reconstruct(
      make_grid(),
      make_metadata(),
      moved_terrain_override(bad_terrain));
  assert(!report.ok());
  assert_report_flags(report);
  assert(report.terrain_override_used);
  assert(!report.allocated_buffers);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(std::string_view(report.result.message).find("override") !=
         std::string_view::npos);
  assert_views_empty(provider);
}

void test_bad_override_stride_is_rejected_before_allocation() {
  const auto terrain = make_override_terrain(0.0F);
  auto bad_view = terrain.view();
  bad_view.stride_j = bad_view.nx + 1;
  const tywrf::dynamics::KrosaBaseStateProviderTerrainOverride bad_override{
      bad_view,
      "moved_candidate_HGT",
      "override:moved_candidate_HGT"};

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report =
      provider.reconstruct(make_grid(), make_metadata(), bad_override);
  assert(!report.ok());
  assert_report_flags(report);
  assert(report.terrain_override_used);
  assert(!report.allocated_buffers);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(std::string_view(report.result.message).find("canonical") !=
         std::string_view::npos);
  assert_views_empty(provider);
}

void test_nonfinite_override_value_is_rejected_before_allocation() {
  auto terrain = make_override_terrain(0.0F);
  auto view = terrain.view();
  view(view.halo.i_lower + 1, view.halo.j_lower + 1) =
      std::numeric_limits<float>::quiet_NaN();

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report = provider.reconstruct(
      make_grid(),
      make_metadata(),
      moved_terrain_override(terrain));
  assert(!report.ok());
  assert_report_flags(report);
  assert(report.terrain_override_used);
  assert(!report.allocated_buffers);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(std::string_view(report.result.message).find("finite") !=
         std::string_view::npos);
  assert_views_empty(provider);
}

void test_failed_override_reconstruct_clears_previous_views() {
  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto success = provider.reconstruct(make_grid(), make_metadata());
  assert(success.ok());
  assert(provider.views().pb.nx > 0);

  auto terrain = make_override_terrain(0.0F);
  auto view = terrain.view();
  view(view.halo.i_lower, view.halo.j_lower) =
      std::numeric_limits<float>::quiet_NaN();
  const auto failure = provider.reconstruct(
      make_grid(),
      make_metadata(),
      moved_terrain_override(terrain));
  assert(!failure.ok());
  assert(!failure.allocated_buffers);
  assert(failure.terrain_override_used);
  assert_views_empty(provider);
}

void test_failed_reconstruct_clears_previous_views() {
  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto success = provider.reconstruct(make_grid(), make_metadata());
  assert(success.ok());
  assert(provider.views().alb.nx > 0);

  auto metadata = make_metadata();
  metadata.terrain_height_m.clear();
  const auto failure = provider.reconstruct(make_grid(), metadata);
  assert(!failure.ok());
  assert(!failure.allocated_buffers);

  assert_views_empty(provider);
}

void test_bad_coefficient_count_is_rejected_before_allocation() {
  auto metadata = make_metadata();
  metadata.c3h.pop_back();

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report = provider.reconstruct(make_grid(), metadata);
  assert(!report.ok());
  assert_report_flags(report);
  assert(!report.allocated_buffers);
  assert(report.c3h_count == 1);
  assert(report.expected_mass_level_count == 2);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(std::string_view(report.result.message).find("coefficient") != std::string_view::npos);
}

void test_missing_p_top_source_is_rejected() {
  auto metadata = make_metadata();
  metadata.p_top_source = tywrf::io::PressureRefreshPTopSource::missing;

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report = provider.reconstruct(make_grid(), metadata);
  assert(!report.ok());
  assert_report_flags(report);
  assert(!report.allocated_buffers);
  assert(!report.p_top_present);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
}

void test_invalid_grid_shape_is_rejected() {
  auto metadata = make_metadata();
  tywrf::Grid grid({4, 3, 2, 4, tywrf::uniform_halo_3d(1)});

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report = provider.reconstruct(grid, metadata);
  assert(!report.ok());
  assert_report_flags(report);
  assert(!report.allocated_buffers);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_configuration);
}

void test_exposed_mub_recompute_api_contract_is_cuda_ready() {
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::KrosaExposedBaseStateRecomputeInputs>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::KrosaExposedBaseStateRecomputeInputs>);
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::KrosaExposedBaseStateRecomputeOutputs>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::KrosaExposedBaseStateRecomputeOutputs>);
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::KrosaExposedBaseStateRecomputeReport>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::KrosaExposedBaseStateRecomputeReport>);
}

void test_exposed_mub_recompute_writes_only_exposed_provider_buffers() {
  const auto grid = make_grid();
  constexpr float p_top = 5'000.0F;
  const std::array<float, 2> c3h{0.72F, 0.28F};
  const std::array<float, 2> c4h{750.0F, 125.0F};

  tywrf::FieldStorage2D<float> mub(grid.surface_layout());
  tywrf::FieldStorage3D<float> pb(grid.mass_layout());
  tywrf::FieldStorage3D<float> t_init(grid.mass_layout());
  tywrf::FieldStorage3D<float> alb(grid.mass_layout());
  tywrf::FieldStorage3D<float> phb(grid.w_layout());
  tywrf::State<float> state(grid);
  fill_interpolated_mub(mub);
  fill_storage(pb, kSentinel);
  fill_storage(t_init, kSentinel);
  fill_storage(alb, kSentinel);
  fill_storage(phb, 12'345.0F);
  fill_storage(state.t, 54'321.0F);
  assert(t_init.data() != state.t.data());

  const auto mub_before = snapshot(mub);
  const auto phb_before = snapshot(phb);
  const auto state_t_before = snapshot(state.t);
  const auto overlap = make_mass_overlap_window(1, 1, 2, 1);

  const auto report =
      tywrf::dynamics::recompute_exposed_base_state_from_mub(
          overlap,
          {static_cast<const tywrf::FieldStorage2D<float>&>(mub).view(),
           {c3h.data(), static_cast<std::int32_t>(c3h.size())},
           {c4h.data(), static_cast<std::int32_t>(c4h.size())},
           p_top},
          {pb.view(), t_init.view(), alb.view()});

  assert(report.ok());
  assert_exposed_recompute_flags(report);
  assert(report.read_interpolated_mub);
  assert(report.wrote_pb);
  assert(report.wrote_t_init);
  assert(report.wrote_alb);
  assert(report.active_nx == 4);
  assert(report.active_ny == 3);
  assert(report.active_nz == 2);
  assert(report.exposed_region_count == 4);
  assert(report.exposed_mass_cell_count == 10);
  assert(report.recomputed_point_count == 20);
  assert(report.pb_recomputed_point_count == 20);
  assert(report.t_init_recomputed_point_count == 20);
  assert(report.alb_recomputed_point_count == 20);

  assert_storage_unchanged(mub, mub_before, "MUB input");
  assert_storage_unchanged(phb, phb_before, "PHB");
  assert_storage_unchanged(state.t, state_t_before, "State::t");

  const auto mass_layout = pb.layout();
  const auto surface_layout = mub.layout();
  const auto pb_view = pb.view();
  const auto t_view = t_init.view();
  const auto alb_view = alb.view();
  const auto mub_view = mub.view();
  bool saw_provider_t_init_write = false;
  for (std::int32_t j = 0; j < mass_layout.ny; ++j) {
    for (std::int32_t k = 0; k < mass_layout.nz; ++k) {
      for (std::int32_t i = 0; i < mass_layout.nx; ++i) {
        const bool active =
            i >= mass_layout.i_begin() && i < mass_layout.i_end() &&
            j >= mass_layout.j_begin() && j < mass_layout.j_end() &&
            k >= mass_layout.k_begin() && k < mass_layout.k_end();
        const auto active_i = i - mass_layout.i_begin();
        const auto active_j = j - mass_layout.j_begin();
        const auto active_k = k - mass_layout.k_begin();
        const bool exposed =
            active && horizontally_exposed(overlap, active_i, active_j);

        if (!exposed) {
          assert(pb_view(i, j, k) == kSentinel);
          assert(t_view(i, j, k) == kSentinel);
          assert(alb_view(i, j, k) == kSentinel);
          continue;
        }

        const auto mub_i = surface_layout.i_begin() + active_i;
        const auto mub_j = surface_layout.j_begin() + active_j;
        const auto expected_pb =
            static_cast<double>(c3h[static_cast<std::size_t>(active_k)]) *
                static_cast<double>(mub_view(mub_i, mub_j)) +
            static_cast<double>(c4h[static_cast<std::size_t>(active_k)]) +
            static_cast<double>(p_top);
        const auto expected_t = expected_t_init_from_pb(expected_pb);
        const auto expected_alb =
            expected_alb_from_pb_t_init(expected_pb, expected_t);
        assert(std::abs(static_cast<double>(pb_view(i, j, k)) - expected_pb) <
               1.0e-2);
        assert(std::abs(t_view(i, j, k) - expected_t) < 1.0e-5F);
        assert(std::abs(alb_view(i, j, k) - expected_alb) < 1.0e-7F);
        saw_provider_t_init_write = true;
      }
    }
  }
  assert(saw_provider_t_init_write);
}

void test_exposed_mub_recompute_no_exposed_no_writes() {
  const auto grid = make_grid();
  constexpr float p_top = 5'000.0F;
  const std::array<float, 2> c3h{0.72F, 0.28F};
  const std::array<float, 2> c4h{750.0F, 125.0F};

  tywrf::FieldStorage2D<float> mub(grid.surface_layout());
  tywrf::FieldStorage3D<float> pb(grid.mass_layout());
  tywrf::FieldStorage3D<float> t_init(grid.mass_layout());
  tywrf::FieldStorage3D<float> alb(grid.mass_layout());
  fill_interpolated_mub(mub);
  fill_storage(pb, kSentinel);
  fill_storage(t_init, kSentinel);
  fill_storage(alb, kSentinel);
  const auto pb_before = snapshot(pb);
  const auto t_init_before = snapshot(t_init);
  const auto alb_before = snapshot(alb);

  const auto report =
      tywrf::dynamics::recompute_exposed_base_state_from_mub(
          make_mass_overlap_window(0, 0, 4, 3),
          {static_cast<const tywrf::FieldStorage2D<float>&>(mub).view(),
           {c3h.data(), static_cast<std::int32_t>(c3h.size())},
           {c4h.data(), static_cast<std::int32_t>(c4h.size())},
           p_top},
          {pb.view(), t_init.view(), alb.view()});

  assert(report.ok());
  assert_exposed_recompute_flags(report);
  assert(!report.read_interpolated_mub);
  assert(!report.wrote_pb);
  assert(!report.wrote_t_init);
  assert(!report.wrote_alb);
  assert(report.exposed_region_count == 0);
  assert(report.exposed_mass_cell_count == 0);
  assert(report.recomputed_point_count == 0);
  assert_storage_unchanged(pb, pb_before, "PB");
  assert_storage_unchanged(t_init, t_init_before, "T_INIT");
  assert_storage_unchanged(alb, alb_before, "ALB");
}

void test_exposed_mub_recompute_invalid_shape_or_window_does_not_write() {
  const auto grid = make_grid();
  constexpr float p_top = 5'000.0F;
  const std::array<float, 2> c3h{0.72F, 0.28F};
  const std::array<float, 2> c4h{750.0F, 125.0F};

  {
    tywrf::FieldStorage2D<float> mub(grid.surface_layout());
    tywrf::FieldStorage3D<float> pb(grid.mass_layout());
    tywrf::FieldStorage3D<float> t_init(grid.mass_layout());
    tywrf::FieldStorage3D<float> alb(grid.mass_layout());
    fill_interpolated_mub(mub);
    fill_storage(pb, kSentinel);
    fill_storage(t_init, kSentinel);
    fill_storage(alb, kSentinel);
    const auto pb_before = snapshot(pb);
    const auto t_init_before = snapshot(t_init);
    const auto alb_before = snapshot(alb);

    const auto report =
        tywrf::dynamics::recompute_exposed_base_state_from_mub(
            make_mass_overlap_window(3, 0, 2, 3),
            {static_cast<const tywrf::FieldStorage2D<float>&>(mub).view(),
             {c3h.data(), static_cast<std::int32_t>(c3h.size())},
             {c4h.data(), static_cast<std::int32_t>(c4h.size())},
             p_top},
            {pb.view(), t_init.view(), alb.view()});

    assert(!report.ok());
    assert_exposed_recompute_flags(report);
    assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
    assert(!report.wrote_pb);
    assert(!report.wrote_t_init);
    assert(!report.wrote_alb);
    assert_storage_unchanged(pb, pb_before, "PB invalid window");
    assert_storage_unchanged(t_init, t_init_before, "T_INIT invalid window");
    assert_storage_unchanged(alb, alb_before, "ALB invalid window");
  }

  {
    tywrf::FieldStorage2D<float> mub(grid.surface_layout());
    tywrf::FieldStorage3D<float> pb(grid.mass_layout());
    tywrf::FieldStorage3D<float> t_init(
        tywrf::ActiveShape3D{5, 3, 2},
        tywrf::uniform_halo_3d(1));
    tywrf::FieldStorage3D<float> alb(grid.mass_layout());
    fill_interpolated_mub(mub);
    fill_storage(pb, kSentinel);
    fill_storage(t_init, kSentinel);
    fill_storage(alb, kSentinel);
    const auto pb_before = snapshot(pb);
    const auto t_init_before = snapshot(t_init);
    const auto alb_before = snapshot(alb);

    const auto report =
        tywrf::dynamics::recompute_exposed_base_state_from_mub(
            make_mass_overlap_window(1, 1, 2, 1),
            {static_cast<const tywrf::FieldStorage2D<float>&>(mub).view(),
             {c3h.data(), static_cast<std::int32_t>(c3h.size())},
             {c4h.data(), static_cast<std::int32_t>(c4h.size())},
             p_top},
            {pb.view(), t_init.view(), alb.view()});

    assert(!report.ok());
    assert_exposed_recompute_flags(report);
    assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
    assert(!report.wrote_pb);
    assert(!report.wrote_t_init);
    assert(!report.wrote_alb);
    assert_storage_unchanged(pb, pb_before, "PB invalid shape");
    assert_storage_unchanged(t_init, t_init_before, "T_INIT invalid shape");
    assert_storage_unchanged(alb, alb_before, "ALB invalid shape");
  }
}

void smoke_real_krosa_file_if_configured(
    const char* env_name,
    const tywrf::Grid& grid) {
  const char* path = std::getenv(env_name);
  if (path == nullptr || std::string_view(path).empty()) {
    return;
  }

  tywrf::FieldStorage3D<float> direct_alb(grid.mass_layout());
  const auto read_result = tywrf::io::read_krosa_pressure_refresh_inputs(
      std::filesystem::path(path),
      grid,
      direct_alb);
  assert(read_result.base_state_reconstruction_inputs_ready());

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto report = provider.reconstruct(grid, read_result.metadata);
  assert(report.ok());
  assert_report_flags(report);
  assert(report.wrote_pb);
  assert(report.wrote_t_init);
  assert(report.wrote_mub);
  assert(report.wrote_alb);
  assert(report.wrote_phb);
}

void test_real_krosa_smoke_if_configured() {
  smoke_real_krosa_file_if_configured(
      "TYWRF_KROSA_BASE_STATE_D01",
      tywrf::Grid({265, 429, 59, 60, tywrf::uniform_halo_3d(1)}));
  smoke_real_krosa_file_if_configured(
      "TYWRF_KROSA_BASE_STATE_D02",
      tywrf::Grid({210, 210, 59, 60, tywrf::uniform_halo_3d(1)}));
}

}  // namespace

int main() {
  try {
    test_success_reconstructs_and_holds_views();
    test_override_terrain_reconstructs_and_changes_output();
    test_missing_terrain_is_rejected_before_allocation();
    test_bad_override_shape_is_rejected_before_allocation();
    test_bad_override_stride_is_rejected_before_allocation();
    test_nonfinite_override_value_is_rejected_before_allocation();
    test_failed_override_reconstruct_clears_previous_views();
    test_failed_reconstruct_clears_previous_views();
    test_bad_coefficient_count_is_rejected_before_allocation();
    test_missing_p_top_source_is_rejected();
    test_invalid_grid_shape_is_rejected();
    test_exposed_mub_recompute_api_contract_is_cuda_ready();
    test_exposed_mub_recompute_writes_only_exposed_provider_buffers();
    test_exposed_mub_recompute_no_exposed_no_writes();
    test_exposed_mub_recompute_invalid_shape_or_window_does_not_write();
    test_real_krosa_smoke_if_configured();
  } catch (const std::exception& error) {
    std::cerr << "base-state provider test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
