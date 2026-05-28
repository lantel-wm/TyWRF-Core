#include "tywrf/nest/base_state_source_staging.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

constexpr float kMask = -876'543.0F;
constexpr float kCandidateSentinel = -9'999.0F;

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

tywrf::nest::ExposedBaseStateViews<const float> source_views(
    const BaseStateStorage& storage) {
  return {
      storage.phb.view(),
      storage.mub.view(),
      storage.pb.view(),
      storage.t_init.view(),
      storage.alb.view(),
      storage.ht.view()};
}

template <typename Storage>
void fill_storage(Storage& storage, const float value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

template <typename Storage>
void fill_storage_indexed(Storage& storage, const float base) {
  for (std::size_t index = 0; index < storage.size(); ++index) {
    storage.data()[index] = base + static_cast<float>(index);
  }
}

void fill_base_state(BaseStateStorage& storage, const float base) {
  fill_storage_indexed(storage.phb, base + 1'000.0F);
  fill_storage_indexed(storage.mub, base + 2'000.0F);
  fill_storage_indexed(storage.pb, base + 3'000.0F);
  fill_storage_indexed(storage.t_init, base + 4'000.0F);
  fill_storage_indexed(storage.alb, base + 5'000.0F);
  fill_storage_indexed(storage.ht, base + 6'000.0F);
}

void fill_candidate(BaseStateStorage& storage) {
  fill_storage(storage.phb, kCandidateSentinel);
  fill_storage(storage.mub, kCandidateSentinel);
  fill_storage(storage.pb, kCandidateSentinel);
  fill_storage(storage.t_init, kCandidateSentinel);
  fill_storage(storage.alb, kCandidateSentinel);
  fill_storage(storage.ht, kCandidateSentinel);
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

struct CandidateSnapshots {
  std::vector<float> phb;
  std::vector<float> mub;
  std::vector<float> pb;
  std::vector<float> t_init;
  std::vector<float> alb;
  std::vector<float> ht;
};

CandidateSnapshots snapshot_candidate(const BaseStateStorage& candidate) {
  return {
      snapshot(candidate.phb),
      snapshot(candidate.mub),
      snapshot(candidate.pb),
      snapshot(candidate.t_init),
      snapshot(candidate.alb),
      snapshot(candidate.ht)};
}

void expect_candidate_unchanged(
    const BaseStateStorage& candidate,
    const CandidateSnapshots& before) {
  expect_storage_unchanged(candidate.phb, before.phb, "candidate PHB");
  expect_storage_unchanged(candidate.mub, before.mub, "candidate MUB");
  expect_storage_unchanged(candidate.pb, before.pb, "candidate PB");
  expect_storage_unchanged(candidate.t_init, before.t_init, "candidate T_INIT");
  expect_storage_unchanged(candidate.alb, before.alb, "candidate ALB");
  expect_storage_unchanged(candidate.ht, before.ht, "candidate HT");
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

void expect_staged_3d(
    const tywrf::FieldView3D<const float>& source,
    const tywrf::FieldView3D<const float>& staged,
    const tywrf::nest::RemapWindow& overlap,
    const float mask,
    const std::string_view label) {
  for (std::int32_t j = 0; j < staged.ny; ++j) {
    for (std::int32_t k = 0; k < staged.nz; ++k) {
      for (std::int32_t i = 0; i < staged.nx; ++i) {
        const bool active =
            i >= staged.halo.i_lower && i < staged.nx - staged.halo.i_upper &&
            j >= staged.halo.j_lower && j < staged.ny - staged.halo.j_upper &&
            k >= staged.halo.k_lower && k < staged.nz - staged.halo.k_upper;
        const auto active_i = i - staged.halo.i_lower;
        const auto active_j = j - staged.halo.j_lower;
        const auto active_k = k - staged.halo.k_lower;
        const bool exposed =
            active && horizontally_exposed(overlap, active_i, active_j);
        const auto expected =
            exposed ? source(
                          source.halo.i_lower + active_i,
                          source.halo.j_lower + active_j,
                          source.halo.k_lower + active_k)
                    : mask;
        if (staged(i, j, k) != expected) {
          std::cerr << label << " mismatch at " << i << ',' << j << ',' << k
                    << ": got " << staged(i, j, k) << ", expected "
                    << expected << '\n';
          assert(false);
        }
      }
    }
  }
}

void expect_staged_2d(
    const tywrf::FieldView2D<const float>& source,
    const tywrf::FieldView2D<const float>& staged,
    const tywrf::nest::RemapWindow& overlap,
    const float mask,
    const std::string_view label) {
  for (std::int32_t j = 0; j < staged.ny; ++j) {
    for (std::int32_t i = 0; i < staged.nx; ++i) {
      const bool active =
          i >= staged.halo.i_lower && i < staged.nx - staged.halo.i_upper &&
          j >= staged.halo.j_lower && j < staged.ny - staged.halo.j_upper;
      const auto active_i = i - staged.halo.i_lower;
      const auto active_j = j - staged.halo.j_lower;
      const bool exposed =
          active && horizontally_exposed(overlap, active_i, active_j);
      const auto expected =
          exposed ? source(source.halo.i_lower + active_i,
                           source.halo.j_lower + active_j)
                  : mask;
      if (staged(i, j) != expected) {
        std::cerr << label << " mismatch at " << i << ',' << j << ": got "
                  << staged(i, j) << ", expected " << expected << '\n';
        assert(false);
      }
    }
  }
}

void expect_empty_views(const tywrf::nest::BaseStateSourceStagingProvider& provider) {
  const auto views = provider.views();
  assert(views.phb.nx == 0);
  assert(views.mub.nx == 0);
  assert(views.pb.nx == 0);
  assert(views.t_init.nx == 0);
  assert(views.alb.nx == 0);
  assert(views.ht.nx == 0);
}

void expect_non_candidate_report(
    const tywrf::nest::BaseStateSourceStagingReport& report) {
  assert(report.source == "explicit_base_state_source_staging_provider");
  assert(
      report.disposition ==
      "diagnostic_only_source_staging_no_gate_no_integrator_no_candidate");
  assert(report.source_shape == "child_shaped_explicit_source_views");
  assert(report.diagnostic_only);
  assert(!report.gate_candidate);
  assert(!report.integrator_output);
  assert(!report.selected_field_numerics_enabled);
  assert(!report.enables_selected_field_numerics);
  assert(!report.writes_netcdf);
  assert(!report.writes_candidate);
  assert(report.candidate_buffers_preserved);
  assert(!report.uses_reference_end_truth);
  assert(!report.uses_direct_p_shortcut);
  assert(!report.reads_direct_p);
}

void test_stages_only_exposed_source_values_into_owned_child_buffers() {
  const auto grid = make_child_grid();
  BaseStateStorage source(grid);
  BaseStateStorage candidate(grid);
  fill_base_state(source, 10'000.0F);
  fill_candidate(candidate);
  const auto candidate_before = snapshot_candidate(candidate);

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  tywrf::nest::BaseStateSourceStagingProvider provider;
  const auto report = provider.stage(
      grid,
      remap,
      source_views(source),
      tywrf::nest::BaseStateSourceStagingOptions{.mask_value = kMask});

  assert(report.ok());
  expect_non_candidate_report(report);
  assert(report.owns_staging_buffers);
  assert(report.allocated_buffers);
  assert(report.active_nx == 4);
  assert(report.active_ny == 3);
  assert(report.mass_nz == 2);
  assert(report.full_nz == 3);
  assert(report.exposed_region_count == 3);
  assert(report.active_mass_cell_count == 12);
  assert(report.active_mass_point_count == 24);
  assert(report.active_surface_cell_count == 12);
  assert(report.active_w_full_point_count == 36);
  assert(report.exposed_mass_cell_count == 3);
  assert(report.exposed_mass_point_count == 6);
  assert(report.exposed_surface_cell_count == 3);
  assert(report.exposed_w_full_column_count == 3);
  assert(report.exposed_w_full_point_count == 9);
  assert(report.masked_mass_cell_count == 9);
  assert(report.masked_mass_point_count == 18);
  assert(report.masked_surface_cell_count == 9);
  assert(report.masked_w_full_column_count == 9);
  assert(report.masked_w_full_point_count == 27);
  assert(report.staged_phb_point_count == 9);
  assert(report.staged_mub_cell_count == 3);
  assert(report.staged_ht_cell_count == 3);
  assert(report.staged_pb_point_count == 6);
  assert(report.staged_t_init_point_count == 6);
  assert(report.staged_alb_point_count == 6);
  assert(report.staged_value_count == 33);
  assert(report.active_masked_value_count == 99);
  assert(report.halo_masked_value_count > 0);
  assert(report.invalid_exposed_value_count == 0);

  const auto staged = provider.views();
  assert(staged.phb.data != source.phb.data());
  assert(staged.mub.data != source.mub.data());
  assert(staged.pb.data != source.pb.data());
  assert(staged.t_init.data != source.t_init.data());
  assert(staged.alb.data != source.alb.data());
  assert(staged.ht.data != source.ht.data());
  assert(staged.phb.data != candidate.phb.data());
  assert(staged.mub.data != candidate.mub.data());
  assert(staged.pb.data != candidate.pb.data());
  assert(staged.t_init.data != candidate.t_init.data());
  assert(staged.alb.data != candidate.alb.data());
  assert(staged.ht.data != candidate.ht.data());
  assert(staged.pb.stride_i == 1);
  assert(staged.pb.stride_k == staged.pb.nx);
  assert(staged.pb.stride_j == staged.pb.nx * staged.pb.nz);

  const auto source_view = source_views(source);
  expect_staged_3d(source_view.phb, staged.phb, remap.w_full, kMask, "PHB");
  expect_staged_2d(source_view.mub, staged.mub, remap.surface, kMask, "MUB");
  expect_staged_2d(source_view.ht, staged.ht, remap.surface, kMask, "HT");
  expect_staged_3d(source_view.pb, staged.pb, remap.mass, kMask, "PB");
  expect_staged_3d(source_view.t_init, staged.t_init, remap.mass, kMask, "T_INIT");
  expect_staged_3d(source_view.alb, staged.alb, remap.mass, kMask, "ALB");
  expect_candidate_unchanged(candidate, candidate_before);
}

void test_no_move_masks_all_values_and_preserves_candidate() {
  const auto grid = make_child_grid();
  BaseStateStorage source(grid);
  BaseStateStorage candidate(grid);
  fill_base_state(source, 20'000.0F);
  fill_candidate(candidate);
  const auto candidate_before = snapshot_candidate(candidate);

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {1, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  tywrf::nest::BaseStateSourceStagingProvider provider;
  const auto report = provider.stage(
      grid,
      remap,
      source_views(source),
      tywrf::nest::BaseStateSourceStagingOptions{.mask_value = kMask});

  assert(report.ok());
  expect_non_candidate_report(report);
  assert(report.allocated_buffers);
  assert(report.exposed_region_count == 0);
  assert(report.exposed_mass_cell_count == 0);
  assert(report.exposed_mass_point_count == 0);
  assert(report.exposed_surface_cell_count == 0);
  assert(report.exposed_w_full_column_count == 0);
  assert(report.exposed_w_full_point_count == 0);
  assert(report.staged_value_count == 0);
  assert(report.active_masked_value_count == 132);

  const auto staged = provider.views();
  const auto source_view = source_views(source);
  expect_staged_3d(source_view.phb, staged.phb, remap.w_full, kMask, "PHB");
  expect_staged_2d(source_view.mub, staged.mub, remap.surface, kMask, "MUB");
  expect_staged_2d(source_view.ht, staged.ht, remap.surface, kMask, "HT");
  expect_staged_3d(source_view.pb, staged.pb, remap.mass, kMask, "PB");
  expect_staged_3d(source_view.t_init, staged.t_init, remap.mass, kMask, "T_INIT");
  expect_staged_3d(source_view.alb, staged.alb, remap.mass, kMask, "ALB");
  expect_candidate_unchanged(candidate, candidate_before);
}

void test_shape_mismatch_rejected_without_allocation_or_candidate_write() {
  const auto grid = make_child_grid();
  BaseStateStorage source(grid);
  BaseStateStorage candidate(grid);
  fill_base_state(source, 30'000.0F);
  fill_candidate(candidate);
  const auto candidate_before = snapshot_candidate(candidate);
  tywrf::FieldStorage3D<float> bad_pb(
      tywrf::ActiveShape3D{5, 3, 2},
      tywrf::uniform_halo_3d(1));
  fill_storage(bad_pb, 1.0F);

  auto explicit_views = source_views(source);
  explicit_views.pb =
      static_cast<const tywrf::FieldStorage3D<float>&>(bad_pb).view();
  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  tywrf::nest::BaseStateSourceStagingProvider provider;
  const auto report = provider.stage(
      grid,
      remap,
      explicit_views,
      tywrf::nest::BaseStateSourceStagingOptions{.mask_value = kMask});

  assert(!report.ok());
  expect_non_candidate_report(report);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(!report.allocated_buffers);
  assert(!report.owns_staging_buffers);
  expect_empty_views(provider);
  expect_candidate_unchanged(candidate, candidate_before);
}

void test_invalid_view_rejected_without_allocation() {
  const auto grid = make_child_grid();
  BaseStateStorage source(grid);
  fill_base_state(source, 40'000.0F);

  auto explicit_views = source_views(source);
  explicit_views.mub.stride_j = explicit_views.mub.nx + 1;
  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  tywrf::nest::BaseStateSourceStagingProvider provider;
  const auto report = provider.stage(
      grid,
      remap,
      explicit_views,
      tywrf::nest::BaseStateSourceStagingOptions{.mask_value = kMask});

  assert(!report.ok());
  expect_non_candidate_report(report);
  assert(report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(!report.allocated_buffers);
  expect_empty_views(provider);
}

void find_exposed_and_masked_mass_cells(
    const tywrf::nest::RemapWindow& overlap,
    std::int32_t& exposed_i,
    std::int32_t& exposed_j,
    std::int32_t& masked_i,
    std::int32_t& masked_j) {
  exposed_i = -1;
  exposed_j = -1;
  masked_i = -1;
  masked_j = -1;
  for (std::int32_t j = 0; j < 3; ++j) {
    for (std::int32_t i = 0; i < 4; ++i) {
      if (horizontally_exposed(overlap, i, j)) {
        exposed_i = i;
        exposed_j = j;
      } else {
        masked_i = i;
        masked_j = j;
      }
    }
  }
  assert(exposed_i >= 0);
  assert(masked_i >= 0);
}

void test_nonfinite_masked_source_is_ignored_but_exposed_source_fails() {
  const auto grid = make_child_grid();
  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  std::int32_t exposed_i = -1;
  std::int32_t exposed_j = -1;
  std::int32_t masked_i = -1;
  std::int32_t masked_j = -1;
  find_exposed_and_masked_mass_cells(
      remap.mass, exposed_i, exposed_j, masked_i, masked_j);

  BaseStateStorage masked_bad_source(grid);
  fill_base_state(masked_bad_source, 50'000.0F);
  auto masked_bad_pb = masked_bad_source.pb.view();
  masked_bad_pb(
      masked_bad_pb.halo.i_lower + masked_i,
      masked_bad_pb.halo.j_lower + masked_j,
      masked_bad_pb.halo.k_lower) =
      std::numeric_limits<float>::quiet_NaN();

  tywrf::nest::BaseStateSourceStagingProvider masked_provider;
  const auto masked_report = masked_provider.stage(
      grid,
      remap,
      source_views(masked_bad_source),
      tywrf::nest::BaseStateSourceStagingOptions{.mask_value = kMask});
  assert(masked_report.ok());
  assert(masked_report.invalid_exposed_value_count == 0);

  BaseStateStorage exposed_bad_source(grid);
  fill_base_state(exposed_bad_source, 60'000.0F);
  auto exposed_bad_pb = exposed_bad_source.pb.view();
  exposed_bad_pb(
      exposed_bad_pb.halo.i_lower + exposed_i,
      exposed_bad_pb.halo.j_lower + exposed_j,
      exposed_bad_pb.halo.k_lower) =
      std::numeric_limits<float>::quiet_NaN();

  tywrf::nest::BaseStateSourceStagingProvider exposed_provider;
  const auto exposed_report = exposed_provider.stage(
      grid,
      remap,
      source_views(exposed_bad_source),
      tywrf::nest::BaseStateSourceStagingOptions{.mask_value = kMask});
  assert(!exposed_report.ok());
  expect_non_candidate_report(exposed_report);
  assert(exposed_report.result.status == tywrf::nest::NestStatus::invalid_contract);
  assert(!exposed_report.allocated_buffers);
  assert(exposed_report.exposed_mass_cell_count == 3);
  assert(exposed_report.invalid_exposed_value_count == 1);
  assert(exposed_report.invalid_exposed_pb_point_count == 1);
  assert(exposed_report.invalid_exposed_phb_point_count == 0);
  assert(exposed_report.invalid_exposed_mub_cell_count == 0);
  assert(exposed_report.invalid_exposed_ht_cell_count == 0);
  assert(exposed_report.invalid_exposed_t_init_point_count == 0);
  assert(exposed_report.invalid_exposed_alb_point_count == 0);
  expect_empty_views(exposed_provider);
}

}  // namespace

int main() {
  static_assert(std::is_standard_layout_v<
                tywrf::nest::BaseStateSourceStagingOptions>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::nest::BaseStateSourceStagingOptions>);
  static_assert(std::is_standard_layout_v<
                tywrf::nest::BaseStateSourceStagingReport>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::nest::BaseStateSourceStagingReport>);

  test_stages_only_exposed_source_values_into_owned_child_buffers();
  test_no_move_masks_all_values_and_preserves_candidate();
  test_shape_mismatch_rejected_without_allocation_or_candidate_write();
  test_invalid_view_rejected_without_allocation();
  test_nonfinite_masked_source_is_ignored_but_exposed_source_fails();

  std::cout << "Validated base-state source staging provider\n";
  return 0;
}
