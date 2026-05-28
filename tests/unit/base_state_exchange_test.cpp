#include "tywrf/nest/base_state_exchange.hpp"

#include "tywrf/state.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

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
void fill_storage(Storage& storage, const float base) {
  for (std::size_t index = 0; index < storage.size(); ++index) {
    storage.data()[index] = base + static_cast<float>(index);
  }
}

void fill_base_state(BaseStateStorage& storage, const float base) {
  fill_storage(storage.phb, base + 1'000.0F);
  fill_storage(storage.mub, base + 2'000.0F);
  fill_storage(storage.pb, base + 3'000.0F);
  fill_storage(storage.t_init, base + 4'000.0F);
  fill_storage(storage.alb, base + 5'000.0F);
  fill_storage(storage.ht, base + 6'000.0F);
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

void expect_diagnostic_only_report(
    const tywrf::nest::ExposedBaseStateExchangeReport& report) {
  assert(report.source == "diagnostic_exposed_base_state_exchange");
  assert(report.disposition == "diagnostic_only_no_gate_no_selected_field_numerics");
  assert(report.diagnostic_only);
  assert(!report.gate_candidate);
  assert(!report.integrator_output);
  assert(!report.selected_field_numerics_enabled);
  assert(!report.enables_selected_field_numerics);
}

void test_rebalance_false_writes_only_exposed_phb_mub_ht_and_marks_recompute() {
  const auto grid = make_child_grid();
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

  const auto report = tywrf::nest::apply_exposed_base_state_exchange(
      remap, source_views(source), views(child));
  assert(report.ok());
  expect_diagnostic_only_report(report);
  assert(!report.rebalance_requested);
  assert(!report.rebalance_applied);
  assert(report.exposed_region_count == 3);
  assert(report.exposed_mass_cell_count == 3);
  assert(report.exposed_surface_cell_count == 3);
  assert(report.exposed_w_full_column_count == 3);
  assert(!report.wrote_overlap);
  assert(!report.wrote_halo);
  assert(report.wrote_phb);
  assert(report.wrote_mub);
  assert(report.wrote_ht);
  assert(!report.wrote_pb);
  assert(!report.wrote_t_init);
  assert(!report.wrote_alb);
  assert(report.phb_written_point_count == 9);
  assert(report.mub_written_cell_count == 3);
  assert(report.ht_written_cell_count == 3);
  assert(report.direct_write_point_count == 15);
  assert(report.pb_recompute_mark_count == 6);
  assert(report.t_init_recompute_mark_count == 6);
  assert(report.alb_recompute_mark_count == 6);
  assert(report.recompute_mark_count == 18);

  expect_3d_exposed_copy_only(
      source.phb, child.phb, phb_before, remap.w_full, "PHB");
  expect_2d_exposed_copy_only(
      source.mub, child.mub, mub_before, remap.surface, "MUB");
  expect_2d_exposed_copy_only(
      source.ht, child.ht, ht_before, remap.surface, "HT");
  expect_storage_unchanged(child.pb, pb_before, "PB");
  expect_storage_unchanged(child.t_init, t_init_before, "T_INIT");
  expect_storage_unchanged(child.alb, alb_before, "ALB");
}

void test_no_move_has_no_writes_or_recompute_marks() {
  const auto grid = make_child_grid();
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

  const auto report = tywrf::nest::apply_exposed_base_state_exchange(
      remap, source_views(source), views(child));
  assert(report.ok());
  expect_diagnostic_only_report(report);
  assert(report.exposed_region_count == 0);
  assert(report.exposed_mass_cell_count == 0);
  assert(report.exposed_surface_cell_count == 0);
  assert(report.exposed_w_full_column_count == 0);
  assert(!report.wrote_phb);
  assert(!report.wrote_mub);
  assert(!report.wrote_ht);
  assert(report.direct_write_point_count == 0);
  assert(report.pb_recompute_mark_count == 0);
  assert(report.t_init_recompute_mark_count == 0);
  assert(report.alb_recompute_mark_count == 0);
  assert(report.recompute_mark_count == 0);

  expect_storage_unchanged(child.phb, phb_before, "PHB");
  expect_storage_unchanged(child.mub, mub_before, "MUB");
  expect_storage_unchanged(child.pb, pb_before, "PB");
  expect_storage_unchanged(child.t_init, t_init_before, "T_INIT");
  expect_storage_unchanged(child.alb, alb_before, "ALB");
  expect_storage_unchanged(child.ht, ht_before, "HT");
}

void test_rebalance_option_is_blocked_without_selected_field_numerics() {
  const auto grid = make_child_grid();
  BaseStateStorage source(grid);
  BaseStateStorage child(grid);
  fill_base_state(source, 30'000.0F);
  fill_base_state(child, -30'000.0F);
  const auto phb_before = snapshot(child.phb);
  const auto mub_before = snapshot(child.mub);
  const auto ht_before = snapshot(child.ht);

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  const auto report = tywrf::nest::apply_exposed_base_state_exchange(
      remap,
      source_views(source),
      views(child),
      tywrf::nest::ExposedBaseStateExchangeOptions{.rebalance = true});
  assert(!report.ok());
  assert(report.result.status == tywrf::nest::NestStatus::not_implemented);
  expect_diagnostic_only_report(report);
  assert(report.rebalance_requested);
  assert(!report.rebalance_applied);
  assert(report.direct_write_point_count == 0);
  assert(report.recompute_mark_count == 0);
  expect_storage_unchanged(child.phb, phb_before, "PHB");
  expect_storage_unchanged(child.mub, mub_before, "MUB");
  expect_storage_unchanged(child.ht, ht_before, "HT");
}

}  // namespace

int main() {
  static_assert(
      std::is_standard_layout_v<
          tywrf::nest::ExposedBaseStateViews<float>>);
  static_assert(
      std::is_trivially_copyable_v<
          tywrf::nest::ExposedBaseStateViews<float>>);

  test_rebalance_false_writes_only_exposed_phb_mub_ht_and_marks_recompute();
  test_no_move_has_no_writes_or_recompute_marks();
  test_rebalance_option_is_blocked_without_selected_field_numerics();
  return 0;
}
