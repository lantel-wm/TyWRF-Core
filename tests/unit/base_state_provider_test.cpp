#include "tywrf/dynamics/base_state_provider.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

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

void assert_report_flags(const tywrf::dynamics::KrosaBaseStateProviderReport& report) {
  assert(report.source == std::string_view("base_state_reconstruction_provider"));
  assert(report.diagnostic_only);
  assert(!report.gate_candidate);
  assert(!report.integrator_output);
  assert(!report.calls_pressure_refresh_compute);
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

  const auto views = provider.views();
  assert(views.pb.nx == 0);
  assert(views.t_init.nx == 0);
  assert(views.mub.nx == 0);
  assert(views.alb.nx == 0);
  assert(views.phb.nx == 0);
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
    test_missing_terrain_is_rejected_before_allocation();
    test_failed_reconstruct_clears_previous_views();
    test_bad_coefficient_count_is_rejected_before_allocation();
    test_missing_p_top_source_is_rejected();
    test_invalid_grid_shape_is_rejected();
    test_real_krosa_smoke_if_configured();
  } catch (const std::exception& error) {
    std::cerr << "base-state provider test failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
