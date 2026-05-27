#include "tywrf/dynamics/base_state_provider.hpp"
#include "tywrf/dynamics/pressure_refresh_hook.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

constexpr float kSentinel = -9'999.0F;
constexpr float kInitialT = 7.0F;

tywrf::Grid make_grid() {
  return tywrf::Grid({4, 3, 2, 3, tywrf::uniform_halo_3d(1)});
}

tywrf::Grid make_zero_halo_grid() {
  return tywrf::Grid({4, 3, 2, 3, tywrf::uniform_halo_3d(0)});
}

tywrf::nest::RemapWindow make_window(
    const tywrf::nest::HorizontalStagger stagger) {
  return {
      stagger,
      0,
      0,
      1,
      1,
      2,
      1,
      -1,
      -1,
  };
}

tywrf::nest::RemapPlan make_plan() {
  tywrf::nest::RemapPlan plan{};
  plan.mass = make_window(tywrf::nest::HorizontalStagger::mass);
  plan.surface = make_window(tywrf::nest::HorizontalStagger::surface);
  plan.w_full = make_window(tywrf::nest::HorizontalStagger::w_full);
  return plan;
}

std::vector<float> make_terrain() {
  std::vector<float> terrain(12, 0.0F);
  for (std::int32_t j = 0; j < 3; ++j) {
    for (std::int32_t i = 0; i < 4; ++i) {
      terrain[static_cast<std::size_t>(j) * 4U + static_cast<std::size_t>(i)] =
          20.0F + 5.0F * static_cast<float>(j) + static_cast<float>(i);
    }
  }
  return terrain;
}

tywrf::io::KrosaPressureRefreshMetadata make_metadata() {
  tywrf::io::KrosaPressureRefreshMetadata metadata{};
  metadata.p_top_pa = 5'000.0F;
  metadata.p_top_source = tywrf::io::PressureRefreshPTopSource::global_attribute;
  metadata.c3f = {1.0F, 0.5F, 0.0F};
  metadata.c4f = {0.0F, 0.0F, 0.0F};
  metadata.c3h = {0.75F, 0.25F};
  metadata.c4h = {0.0F, 0.0F};
  metadata.terrain_source_name = "HGT";
  metadata.terrain_nx = 4;
  metadata.terrain_ny = 3;
  metadata.terrain_height_m = make_terrain();
  return metadata;
}

template <typename Storage>
void fill_storage(Storage& storage, const float value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

void fill_state(tywrf::State<float>& state) {
  fill_storage(state.u, kSentinel);
  fill_storage(state.v, kSentinel);
  fill_storage(state.w, kSentinel);
  fill_storage(state.ph, 0.0F);
  fill_storage(state.phb, kSentinel);
  fill_storage(state.t, kInitialT);
  fill_storage(state.p, kSentinel);
  fill_storage(state.pb, kSentinel);
  fill_storage(state.qvapor, kSentinel);
  fill_storage(state.qcloud, kSentinel);
  fill_storage(state.qrain, kSentinel);
  fill_storage(state.qice, kSentinel);
  fill_storage(state.qsnow, kSentinel);
  fill_storage(state.qgraup, kSentinel);
  fill_storage(state.qnice, kSentinel);
  fill_storage(state.qnrain, kSentinel);
  fill_storage(state.mu, kSentinel);
  fill_storage(state.mub, kSentinel);
  fill_storage(state.psfc, kSentinel);
  fill_storage(state.u10, kSentinel);
  fill_storage(state.v10, kSentinel);
  fill_storage(state.t2, kSentinel);
  fill_storage(state.q2, kSentinel);
  fill_storage(state.rainc, kSentinel);
  fill_storage(state.rainnc, kSentinel);
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

[[nodiscard]] bool in_window(
    const tywrf::nest::RemapWindow& window,
    const std::int32_t i,
    const std::int32_t j) {
  return i >= window.new_i_begin && i < window.new_i_begin + window.extent_i &&
         j >= window.new_j_begin && j < window.new_j_begin + window.extent_j;
}

template <typename Storage>
std::vector<float> snapshot(const Storage& storage) {
  return {storage.data(), storage.data() + storage.size()};
}

template <typename Storage>
void assert_unchanged(
    const Storage& storage,
    const std::vector<float>& before,
    const std::string_view label) {
  assert(storage.size() == before.size());
  if (!std::equal(storage.data(), storage.data() + storage.size(), before.begin())) {
    std::cerr << label << " changed unexpectedly\n";
    assert(false);
  }
}

struct StateSnapshot {
  std::vector<float> p;
  std::vector<float> pb;
  std::vector<float> phb;
  std::vector<float> t;
  std::vector<float> mub;
};

StateSnapshot snapshot_state(const tywrf::State<float>& state) {
  return {
      snapshot(state.p),
      snapshot(state.pb),
      snapshot(state.phb),
      snapshot(state.t),
      snapshot(state.mub),
  };
}

void assert_state_subset_unchanged(
    const tywrf::State<float>& state,
    const StateSnapshot& before) {
  assert_unchanged(state.p, before.p, "P");
  assert_unchanged(state.pb, before.pb, "PB");
  assert_unchanged(state.phb, before.phb, "PHB");
  assert_unchanged(state.t, before.t, "T");
  assert_unchanged(state.mub, before.mub, "MUB");
}

void assert_halo_unchanged(const tywrf::FieldStorage3D<float>& field) {
  const auto layout = field.layout();
  const auto view = field.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t k = 0; k < layout.nz; ++k) {
      for (std::int32_t i = 0; i < layout.nx; ++i) {
        const bool halo =
            i < layout.i_begin() || i >= layout.i_end() ||
            j < layout.j_begin() || j >= layout.j_end() ||
            k < layout.k_begin() || k >= layout.k_end();
        if (halo) {
          assert(view(i, j, k) == kSentinel);
        }
      }
    }
  }
}

void assert_halo_unchanged(const tywrf::FieldStorage2D<float>& field) {
  const auto layout = field.layout();
  const auto view = field.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t i = 0; i < layout.nx; ++i) {
      const bool halo =
          i < layout.i_begin() || i >= layout.i_end() ||
          j < layout.j_begin() || j >= layout.j_end();
      if (halo) {
        assert(view(i, j) == kSentinel);
      }
    }
  }
}

void assert_provider_reports_are_non_gate(
    const tywrf::dynamics::KrosaPressureRefreshHookReport& report) {
  assert(report.diagnostic_only);
  assert(!report.gate_candidate);
  assert(!report.integrator_output);
  assert(report.provider_report.diagnostic_only);
  assert(!report.provider_report.gate_candidate);
  assert(!report.provider_report.integrator_output);
  assert(!report.provider_report.calls_pressure_refresh_compute);
  assert(report.staging_report.diagnostic_only);
  assert(!report.staging_report.gate_candidate);
  assert(!report.staging_report.pressure_refresh_applied);
}

void test_success_syncs_provider_base_state_and_refreshes_exposed_pressure() {
  const auto grid = make_grid();
  auto metadata = make_metadata();
  auto plan = make_plan();
  tywrf::State<float> state(grid);
  fill_state(state);

  tywrf::dynamics::KrosaBaseStateProvider provider;
  const auto provider_report = provider.reconstruct(grid, metadata);
  assert(provider_report.ok());
  const auto provider_views = provider.views();

  const auto report =
      tywrf::dynamics::apply_krosa_moving_nest_pressure_refresh_hook(
          plan,
          state,
          metadata);

  assert(report.ok());
  assert(report.provider_ok);
  assert(report.staging_ok);
  assert(report.pressure_refresh_applied);
  assert(report.calls_pressure_refresh_compute);
  assert(report.synced_pb_point_count == 20);
  assert(report.synced_mub_point_count == 10);
  assert(report.synced_phb_point_count == 30);
  assert(report.compute_report.target_column_count == 10);
  assert(report.compute_report.refreshed_point_count == 20);
  assert(!report.touched_overlap_cells);
  assert(!report.touched_halo_cells);
  assert_provider_reports_are_non_gate(report);
  assert(
      report.staging_report.alb_source ==
      tywrf::dynamics::PressureRefreshAlbSource::
          base_state_reconstruction_provider);

  const auto p_view = state.p.view();
  const auto pb_view = state.pb.view();
  const auto phb_view = state.phb.view();
  const auto t_view = state.t.view();
  const auto provider_pb = provider_views.pb;
  const auto provider_phb = provider_views.phb;
  const auto provider_t_init = provider_views.t_init;
  const auto provider_alb = provider_views.alb;

  for (std::int32_t j = 0; j < active_ny(p_view); ++j) {
    for (std::int32_t k = 0; k < active_nz(p_view); ++k) {
      for (std::int32_t i = 0; i < active_nx(p_view); ++i) {
        const auto ii = p_view.halo.i_lower + i;
        const auto jj = p_view.halo.j_lower + j;
        const auto kk = p_view.halo.k_lower + k;
        if (in_window(plan.mass, i, j)) {
          assert(p_view(ii, jj, kk) == kSentinel);
          assert(pb_view(ii, jj, kk) == kSentinel);
        } else {
          assert(std::isfinite(p_view(ii, jj, kk)));
          assert(p_view(ii, jj, kk) != kSentinel);
          assert(pb_view(ii, jj, kk) == provider_pb(ii, jj, kk));
          assert(pb_view(ii, jj, kk) != provider_alb(ii, jj, kk));
        }
        assert(t_view(ii, jj, kk) == kInitialT);
        assert(t_view(ii, jj, kk) != provider_t_init(ii, jj, kk));
      }
    }
  }

  const auto mub_view = state.mub.view();
  const auto provider_mub = provider_views.mub;
  for (std::int32_t j = 0; j < active_ny(mub_view); ++j) {
    for (std::int32_t i = 0; i < active_nx(mub_view); ++i) {
      const auto ii = mub_view.halo.i_lower + i;
      const auto jj = mub_view.halo.j_lower + j;
      if (in_window(plan.surface, i, j)) {
        assert(mub_view(ii, jj) == kSentinel);
      } else {
        assert(mub_view(ii, jj) == provider_mub(ii, jj));
      }
    }
  }

  for (std::int32_t j = 0; j < active_ny(phb_view); ++j) {
    for (std::int32_t k = 0; k < active_nz(phb_view); ++k) {
      for (std::int32_t i = 0; i < active_nx(phb_view); ++i) {
        const auto ii = phb_view.halo.i_lower + i;
        const auto jj = phb_view.halo.j_lower + j;
        const auto kk = phb_view.halo.k_lower + k;
        if (in_window(plan.w_full, i, j)) {
          assert(phb_view(ii, jj, kk) == kSentinel);
        } else {
          assert(phb_view(ii, jj, kk) == provider_phb(ii, jj, kk));
        }
      }
    }
  }

  assert_halo_unchanged(state.p);
  assert_halo_unchanged(state.pb);
  assert_halo_unchanged(state.phb);
  assert_halo_unchanged(state.mub);
  for (std::size_t index = 0; index < state.qvapor.size(); ++index) {
    assert(state.qvapor.data()[index] == kSentinel);
  }
}

void test_missing_metadata_fails_without_compute_or_partial_sync() {
  const auto grid = make_grid();
  const auto plan = make_plan();

  auto missing_terrain = make_metadata();
  missing_terrain.terrain_nx = 0;
  missing_terrain.terrain_ny = 0;
  missing_terrain.terrain_height_m.clear();

  auto missing_coefficients = make_metadata();
  missing_coefficients.c4h.clear();

  auto missing_p_top = make_metadata();
  missing_p_top.p_top_source = tywrf::io::PressureRefreshPTopSource::missing;

  for (const auto& metadata : {missing_terrain, missing_coefficients, missing_p_top}) {
    tywrf::State<float> state(grid);
    fill_state(state);
    const auto before = snapshot_state(state);

    const auto report =
        tywrf::dynamics::apply_krosa_moving_nest_pressure_refresh_hook(
            plan,
            state,
            metadata);

    assert(!report.ok());
    assert(!report.provider_ok);
    assert(!report.staging_ok);
    assert(!report.pressure_refresh_applied);
    assert(!report.calls_pressure_refresh_compute);
    assert(report.synced_pb_point_count == 0);
    assert(report.synced_mub_point_count == 0);
    assert(report.synced_phb_point_count == 0);
    assert(!report.compute_report.ok() || report.compute_report.refreshed_point_count == 0);
    assert_provider_reports_are_non_gate(report);
    assert_state_subset_unchanged(state, before);
  }
}

void test_zero_halo_grid_path_is_supported() {
  const auto grid = make_zero_halo_grid();
  auto state = tywrf::State<float>(grid);
  fill_state(state);

  const auto report =
      tywrf::dynamics::apply_krosa_moving_nest_pressure_refresh_hook(
          make_plan(),
          state,
          make_metadata());

  assert(report.ok());
  assert(report.pressure_refresh_applied);
  assert(report.synced_pb_point_count == 20);
  assert(report.synced_mub_point_count == 10);
  assert(report.synced_phb_point_count == 30);
  assert(!report.touched_overlap_cells);
  assert(!report.touched_halo_cells);
  assert_provider_reports_are_non_gate(report);
}

}  // namespace

int main() {
  test_success_syncs_provider_base_state_and_refreshes_exposed_pressure();
  test_missing_metadata_fails_without_compute_or_partial_sync();
  test_zero_halo_grid_path_is_supported();

  std::cout << "Validated KROSA pressure refresh hook helper\n";
  return 0;
}
