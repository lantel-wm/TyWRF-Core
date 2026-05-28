#include "tywrf/nest/base_state_exchange_adapter.hpp"

#include "tywrf/state.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

constexpr float kSentinel = -9'999.0F;

struct BaseStateStorage {
  explicit BaseStateStorage(const tywrf::Grid& grid)
      : phb(grid.w_layout()),
        mub(grid.surface_layout()),
        pb(grid.mass_layout()),
        t_init(grid.mass_layout()),
        alb(grid.mass_layout()),
        ht(grid.surface_layout()) {}

  tywrf::FieldStorage3D<float> phb;
  tywrf::FieldStorage2D<float> mub;
  tywrf::FieldStorage3D<float> pb;
  tywrf::FieldStorage3D<float> t_init;
  tywrf::FieldStorage3D<float> alb;
  tywrf::FieldStorage2D<float> ht;
};

tywrf::Grid make_child_grid() {
  return tywrf::Grid({
      .mass_nx = 4,
      .mass_ny = 3,
      .mass_nz = 2,
      .full_nz = 3,
      .halo = tywrf::Halo3D{1, 1, 1, 1, 1, 1},
  });
}

tywrf::nest::ParentChildDescriptor make_small_descriptor() {
  return {
      tywrf::nest::HorizontalDomainDescriptor{1, 2'000, 8, 8, 9, 9},
      tywrf::nest::HorizontalDomainDescriptor{2, 2'000, 4, 3, 5, 4},
      1,
      1,
  };
}

tywrf::nest::RemapPlan make_remap_plan(
    const tywrf::nest::ParentChildPosition from,
    const tywrf::nest::ParentChildPosition to) {
  const auto descriptor = make_small_descriptor();
  return tywrf::nest::build_remap_plan(
      tywrf::nest::make_nest_pose(descriptor, from),
      tywrf::nest::make_nest_pose(descriptor, to));
}

tywrf::nest::ExposedBaseStateViews<const float> views(
    const BaseStateStorage& storage) {
  return {
      storage.phb.view(),
      storage.mub.view(),
      storage.pb.view(),
      storage.t_init.view(),
      storage.alb.view(),
      storage.ht.view()};
}

tywrf::nest::ExposedBaseStateViews<float> views(BaseStateStorage& storage) {
  return {
      storage.phb.view(),
      storage.mub.view(),
      storage.pb.view(),
      storage.t_init.view(),
      storage.alb.view(),
      storage.ht.view()};
}

tywrf::nest::ExposedBaseStateViews<const float> source_views(
    const BaseStateStorage& storage) {
  return views(storage);
}

template <typename Storage>
void fill_storage_indexed(Storage& storage, const float base) {
  for (std::size_t index = 0; index < storage.size(); ++index) {
    storage.data()[index] = base + static_cast<float>(index);
  }
}

template <typename Storage>
void fill_storage(Storage& storage, const float value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

void fill_base_state(BaseStateStorage& storage, const float base) {
  fill_storage_indexed(storage.phb, base + 1'000.0F);
  fill_storage_indexed(storage.mub, base + 80'000.0F);
  fill_storage_indexed(storage.pb, base + 3'000.0F);
  fill_storage_indexed(storage.t_init, base + 4'000.0F);
  fill_storage_indexed(storage.alb, base + 5'000.0F);
  fill_storage_indexed(storage.ht, base + 6'000.0F);
}

template <typename Storage>
std::vector<float> snapshot(const Storage& storage) {
  return std::vector<float>(storage.data(), storage.data() + storage.size());
}

template <typename Storage>
void expect_storage_unchanged(
    const Storage& storage,
    const std::vector<float>& before,
    const std::string_view label) {
  if (before.size() != storage.size() ||
      !std::equal(before.begin(), before.end(), storage.data())) {
    std::cerr << label << " changed unexpectedly\n";
    assert(false);
  }
}

[[nodiscard]] bool horizontally_exposed(
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

void expect_3d_exposed_copy_only(
    const tywrf::FieldStorage3D<float>& source,
    const tywrf::FieldStorage3D<float>& child,
    const std::vector<float>& before,
    const tywrf::nest::RemapWindow& overlap,
    const std::string_view label) {
  const auto layout = child.layout();
  const auto source_view = source.view();
  const auto child_view = child.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t k = 0; k < layout.nz; ++k) {
      for (std::int32_t i = 0; i < layout.nx; ++i) {
        const bool active =
            i >= layout.i_begin() && i < layout.i_end() &&
            j >= layout.j_begin() && j < layout.j_end() &&
            k >= layout.k_begin() && k < layout.k_end();
        const bool exposed =
            active && horizontally_exposed(
                          overlap,
                          i - layout.i_begin(),
                          j - layout.j_begin());
        const auto expected =
            exposed ? source_view(i, j, k) : before[layout.index(i, j, k)];
        if (child_view(i, j, k) != expected) {
          std::cerr << label << " mismatch at " << i << ',' << j << ',' << k
                    << ": got " << child_view(i, j, k) << ", expected "
                    << expected << '\n';
          assert(false);
        }
      }
    }
  }
}

void expect_2d_exposed_copy_only(
    const tywrf::FieldStorage2D<float>& source,
    const tywrf::FieldStorage2D<float>& child,
    const std::vector<float>& before,
    const tywrf::nest::RemapWindow& overlap,
    const std::string_view label) {
  const auto layout = child.layout();
  const auto source_view = source.view();
  const auto child_view = child.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t i = 0; i < layout.nx; ++i) {
      const bool active =
          i >= layout.i_begin() && i < layout.i_end() &&
          j >= layout.j_begin() && j < layout.j_end();
      const bool exposed =
          active && horizontally_exposed(
                        overlap,
                        i - layout.i_begin(),
                        j - layout.j_begin());
      const auto expected =
          exposed ? source_view(i, j) : before[layout.index(i, j)];
      if (child_view(i, j) != expected) {
        std::cerr << label << " mismatch at " << i << ',' << j << ": got "
                  << child_view(i, j) << ", expected " << expected << '\n';
        assert(false);
      }
    }
  }
}

void expect_recomputed_mass_fields_only(
    const BaseStateStorage& child,
    const std::vector<float>& pb_before,
    const std::vector<float>& t_init_before,
    const std::vector<float>& alb_before,
    const tywrf::nest::RemapWindow& overlap,
    const std::array<float, 2>& c3h,
    const std::array<float, 2>& c4h,
    const float p_top) {
  const auto layout = child.pb.layout();
  const auto surface_layout = child.mub.layout();
  const auto pb_view = child.pb.view();
  const auto t_view = child.t_init.view();
  const auto alb_view = child.alb.view();
  const auto mub_view = child.mub.view();
  bool saw_recomputed_cell = false;
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t k = 0; k < layout.nz; ++k) {
      for (std::int32_t i = 0; i < layout.nx; ++i) {
        const bool active =
            i >= layout.i_begin() && i < layout.i_end() &&
            j >= layout.j_begin() && j < layout.j_end() &&
            k >= layout.k_begin() && k < layout.k_end();
        const auto active_i = i - layout.i_begin();
        const auto active_j = j - layout.j_begin();
        const auto active_k = k - layout.k_begin();
        const bool exposed =
            active && horizontally_exposed(overlap, active_i, active_j);

        if (!exposed) {
          assert(pb_view(i, j, k) == pb_before[layout.index(i, j, k)]);
          assert(t_view(i, j, k) == t_init_before[layout.index(i, j, k)]);
          assert(alb_view(i, j, k) == alb_before[layout.index(i, j, k)]);
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
        saw_recomputed_cell = true;
      }
    }
  }
  assert(saw_recomputed_cell);
}

void expect_adapter_flags(
    const tywrf::nest::ExposedBaseStateExchangeAdapterReport& report) {
  assert(report.source == "diagnostic_exposed_base_state_exchange_adapter");
  assert(
      report.disposition ==
      "diagnostic_only_staging_no_gate_no_integrator_no_selected_field_numerics");
  assert(report.exchange_source == "diagnostic_exposed_base_state_exchange");
  assert(report.recompute_source == "exposed_mub_base_state_recompute_provider");
  assert(report.diagnostic_only);
  assert(!report.gate_candidate);
  assert(!report.integrator_output);
  assert(!report.selected_field_numerics_enabled);
  assert(!report.enables_selected_field_numerics);
  assert(!report.writes_netcdf);
  assert(!report.writes_candidate);
  assert(report.ht_source_name == "HT");
  assert(report.ht_diagnostic_label == "HGT");
  assert(report.ht_is_hgt_alias);
  assert(!report.creates_terrain_owner);
  assert(!report.terrain_owner_created);
}

void expect_all_child_storage_unchanged(
    const BaseStateStorage& child,
    const std::vector<float>& phb_before,
    const std::vector<float>& mub_before,
    const std::vector<float>& pb_before,
    const std::vector<float>& t_init_before,
    const std::vector<float>& alb_before,
    const std::vector<float>& ht_before,
    const std::string_view label) {
  expect_storage_unchanged(child.phb, phb_before, label);
  expect_storage_unchanged(child.mub, mub_before, label);
  expect_storage_unchanged(child.pb, pb_before, label);
  expect_storage_unchanged(child.t_init, t_init_before, label);
  expect_storage_unchanged(child.alb, alb_before, label);
  expect_storage_unchanged(child.ht, ht_before, label);
}

void test_adapter_calls_d68_and_d69_and_aggregates_counts() {
  const auto grid = make_child_grid();
  constexpr float p_top = 5'000.0F;
  constexpr std::array<float, 2> c3h{0.72F, 0.28F};
  constexpr std::array<float, 2> c4h{750.0F, 125.0F};

  BaseStateStorage source(grid);
  BaseStateStorage child(grid);
  fill_base_state(source, 10'000.0F);
  fill_base_state(child, -10'000.0F);

  const auto phb_before = snapshot(child.phb);
  const auto mub_before = snapshot(child.mub);
  const auto pb_before = snapshot(child.pb);
  const auto t_init_before = snapshot(child.t_init);
  const auto alb_before = snapshot(child.alb);
  const auto ht_before = snapshot(child.ht);

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  const auto report = tywrf::nest::apply_exposed_base_state_exchange_adapter(
      remap,
      source_views(source),
      views(child),
      {{c3h.data(), static_cast<std::int32_t>(c3h.size())},
       {c4h.data(), static_cast<std::int32_t>(c4h.size())},
       p_top});

  assert(report.ok());
  expect_adapter_flags(report);
  assert(report.called_d68_exchange);
  assert(report.called_d69_recompute);
  assert(!report.rebalance_requested);
  assert(!report.rebalance_applied);
  assert(report.active_nx == 4);
  assert(report.active_ny == 3);
  assert(report.active_nz == 2);
  assert(report.c3h_count == 2);
  assert(report.c4h_count == 2);
  assert(report.exchange_exposed_region_count == 3);
  assert(report.recompute_exposed_region_count == 1);
  assert(report.exposed_region_count == 4);
  assert(report.exposed_mass_cell_count == 3);
  assert(report.exposed_surface_cell_count == 3);
  assert(report.exposed_w_full_column_count == 3);
  assert(!report.wrote_overlap);
  assert(!report.wrote_halo);
  assert(report.wrote_phb);
  assert(report.wrote_mub);
  assert(report.wrote_ht);
  assert(report.wrote_pb);
  assert(report.wrote_t_init);
  assert(report.wrote_alb);
  assert(report.phb_written_point_count == 9);
  assert(report.mub_written_cell_count == 3);
  assert(report.ht_written_cell_count == 3);
  assert(report.direct_write_point_count == 15);
  assert(report.exchange_recompute_mark_count == 18);
  assert(report.recomputed_point_count == 6);
  assert(report.pb_recomputed_point_count == 6);
  assert(report.t_init_recomputed_point_count == 6);
  assert(report.alb_recomputed_point_count == 6);

  expect_3d_exposed_copy_only(
      source.phb, child.phb, phb_before, remap.w_full, "PHB");
  expect_2d_exposed_copy_only(
      source.mub, child.mub, mub_before, remap.surface, "MUB");
  expect_2d_exposed_copy_only(
      source.ht, child.ht, ht_before, remap.surface, "HT/HGT");
  expect_recomputed_mass_fields_only(
      child, pb_before, t_init_before, alb_before, remap.mass, c3h, c4h, p_top);
}

void test_adapter_no_exposed_count_zero_and_no_writes() {
  const auto grid = make_child_grid();
  constexpr float p_top = 5'000.0F;
  constexpr std::array<float, 2> c3h{0.72F, 0.28F};
  constexpr std::array<float, 2> c4h{750.0F, 125.0F};

  BaseStateStorage source(grid);
  BaseStateStorage child(grid);
  fill_base_state(source, 20'000.0F);
  fill_base_state(child, -20'000.0F);

  const auto phb_before = snapshot(child.phb);
  const auto mub_before = snapshot(child.mub);
  const auto pb_before = snapshot(child.pb);
  const auto t_init_before = snapshot(child.t_init);
  const auto alb_before = snapshot(child.alb);
  const auto ht_before = snapshot(child.ht);

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {1, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  const auto report = tywrf::nest::apply_exposed_base_state_exchange_adapter(
      remap,
      source_views(source),
      views(child),
      {{c3h.data(), static_cast<std::int32_t>(c3h.size())},
       {c4h.data(), static_cast<std::int32_t>(c4h.size())},
       p_top});

  assert(report.ok());
  expect_adapter_flags(report);
  assert(report.called_d68_exchange);
  assert(report.called_d69_recompute);
  assert(report.exchange_exposed_region_count == 0);
  assert(report.recompute_exposed_region_count == 0);
  assert(report.exposed_region_count == 0);
  assert(report.exposed_mass_cell_count == 0);
  assert(report.exposed_surface_cell_count == 0);
  assert(report.exposed_w_full_column_count == 0);
  assert(!report.wrote_phb);
  assert(!report.wrote_mub);
  assert(!report.wrote_ht);
  assert(!report.wrote_pb);
  assert(!report.wrote_t_init);
  assert(!report.wrote_alb);
  assert(report.direct_write_point_count == 0);
  assert(report.exchange_recompute_mark_count == 0);
  assert(report.recomputed_point_count == 0);

  expect_all_child_storage_unchanged(
      child,
      phb_before,
      mub_before,
      pb_before,
      t_init_before,
      alb_before,
      ht_before,
      "no-exposed adapter");
}

void test_invalid_d68_or_d69_input_returns_error_without_partial_writes() {
  const auto grid = make_child_grid();
  constexpr float p_top = 5'000.0F;
  constexpr std::array<float, 2> c3h{0.72F, 0.28F};
  constexpr std::array<float, 2> c4h{750.0F, 125.0F};

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  {
    BaseStateStorage source(grid);
    BaseStateStorage child(grid);
    fill_base_state(source, 30'000.0F);
    fill_base_state(child, -30'000.0F);
    const auto phb_before = snapshot(child.phb);
    const auto mub_before = snapshot(child.mub);
    const auto pb_before = snapshot(child.pb);
    const auto t_init_before = snapshot(child.t_init);
    const auto alb_before = snapshot(child.alb);
    const auto ht_before = snapshot(child.ht);

    auto bad_source = source_views(source);
    bad_source.phb.data = nullptr;
    const auto report =
        tywrf::nest::apply_exposed_base_state_exchange_adapter(
            remap,
            bad_source,
            views(child),
            {{c3h.data(), static_cast<std::int32_t>(c3h.size())},
             {c4h.data(), static_cast<std::int32_t>(c4h.size())},
             p_top});

    assert(!report.ok());
    expect_adapter_flags(report);
    assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
    assert(!report.called_d68_exchange);
    assert(!report.called_d69_recompute);
    expect_all_child_storage_unchanged(
        child,
        phb_before,
        mub_before,
        pb_before,
        t_init_before,
        alb_before,
        ht_before,
        "invalid D68 adapter input");
  }

  {
    BaseStateStorage source(grid);
    BaseStateStorage child(grid);
    fill_base_state(source, 40'000.0F);
    fill_base_state(child, -40'000.0F);
    auto source_view = source.mub.view();
    for (std::int32_t j = source_view.halo.j_lower;
         j < source_view.ny - source_view.halo.j_upper;
         ++j) {
      for (std::int32_t i = source_view.halo.i_lower;
           i < source_view.nx - source_view.halo.i_upper;
           ++i) {
        source_view(i, j) = std::numeric_limits<float>::quiet_NaN();
      }
    }
    const auto phb_before = snapshot(child.phb);
    const auto mub_before = snapshot(child.mub);
    const auto pb_before = snapshot(child.pb);
    const auto t_init_before = snapshot(child.t_init);
    const auto alb_before = snapshot(child.alb);
    const auto ht_before = snapshot(child.ht);

    const auto report =
        tywrf::nest::apply_exposed_base_state_exchange_adapter(
            remap,
            source_views(source),
            views(child),
            {{c3h.data(), static_cast<std::int32_t>(c3h.size())},
             {c4h.data(), static_cast<std::int32_t>(c4h.size())},
             p_top});

    assert(!report.ok());
    expect_adapter_flags(report);
    assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
    assert(!report.called_d68_exchange);
    assert(!report.called_d69_recompute);
    assert(report.invalid_column_count == 3);
    assert(report.invalid_point_count == 6);
    expect_all_child_storage_unchanged(
        child,
        phb_before,
        mub_before,
        pb_before,
        t_init_before,
        alb_before,
        ht_before,
        "invalid D69 adapter input");
  }
}

}  // namespace

int main() {
  test_adapter_calls_d68_and_d69_and_aggregates_counts();
  test_adapter_no_exposed_count_zero_and_no_writes();
  test_invalid_d68_or_d69_input_returns_error_without_partial_writes();
  return 0;
}
