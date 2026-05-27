#include "tywrf/nest/state_exchange.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

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

template <typename Storage>
void fill_storage(Storage& storage, const float base) {
  for (std::size_t index = 0; index < storage.size(); ++index) {
    storage.data()[index] = base + static_cast<float>(index);
  }
}

void fill_state(tywrf::State<float>& state) {
  fill_storage(state.u, 1'000.0F);
  fill_storage(state.v, 2'000.0F);
  fill_storage(state.w, 3'000.0F);
  fill_storage(state.ph, 4'000.0F);
  fill_storage(state.phb, 5'000.0F);
  fill_storage(state.t, 6'000.0F);
  fill_storage(state.p, 7'000.0F);
  fill_storage(state.pb, 8'000.0F);
  fill_storage(state.qvapor, 9'000.0F);
  fill_storage(state.qcloud, 10'000.0F);
  fill_storage(state.qrain, 11'000.0F);
  fill_storage(state.qice, 12'000.0F);
  fill_storage(state.qsnow, 13'000.0F);
  fill_storage(state.qgraup, 14'000.0F);
  fill_storage(state.qnice, 15'000.0F);
  fill_storage(state.qnrain, 16'000.0F);
  fill_storage(state.mu, 17'000.0F);
  fill_storage(state.mub, 18'000.0F);
  fill_storage(state.psfc, 19'000.0F);
  fill_storage(state.u10, 20'000.0F);
  fill_storage(state.v10, 21'000.0F);
  fill_storage(state.t2, 22'000.0F);
  fill_storage(state.q2, 23'000.0F);
  fill_storage(state.rainc, 24'000.0F);
  fill_storage(state.rainnc, 25'000.0F);
}

template <typename Storage>
std::vector<float> snapshot(const Storage& storage) {
  return std::vector<float>(storage.data(), storage.data() + storage.size());
}

template <typename Storage>
void expect_unchanged(
    const Storage& storage,
    const std::vector<float>& before,
    const std::string_view label) {
  if (before.size() != storage.size() ||
      !std::equal(before.begin(), before.end(), storage.data())) {
    std::cerr << label << " was modified by state exchange planning\n";
    assert(false);
  }
}

void test_exposed_strip_accounting() {
  tywrf::State<float> child(make_child_grid());
  fill_state(child);
  const auto u_before = snapshot(child.u);
  const auto v_before = snapshot(child.v);
  const auto mu_before = snapshot(child.mu);
  const auto qvapor_before = snapshot(child.qvapor);
  const auto t_before = snapshot(child.t);
  const auto ph_before = snapshot(child.ph);

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  const auto exchange =
      tywrf::nest::build_exposed_child_state_exchange_plan(
          remap, static_cast<const tywrf::State<float>&>(child).view());
  assert(exchange.ok());
  assert(exchange.operation ==
         tywrf::nest::ExchangeOperation::parent_to_child_interpolation);
  assert(exchange.field_count == 6);
  assert(exchange.report.ok());
  assert(exchange.report.planned_field_count == 6);
  assert(exchange.report.active_field_count == 6);
  assert(exchange.report.exposed_region_count == 6);
  assert(exchange.report.exposed_horizontal_cell_count == 19);
  assert(exchange.report.exchange_point_count == 38);
  assert(exchange.report.requires_parent_interpolation);
  assert(!exchange.report.performed_interpolation);
  assert(!exchange.report.modifies_overlap);
  assert(!exchange.report.modifies_halo);

  const auto& u = exchange.fields[0];
  assert(u.field == tywrf::nest::StateExchangeField::u);
  assert(u.stagger == tywrf::nest::HorizontalStagger::u);
  assert(u.active_nx == 5);
  assert(u.active_ny == 3);
  assert(u.active_k_count == 2);
  assert(u.exposed_region_count == 1);
  assert(u.exposed_horizontal_cell_count == 3);
  assert(u.exchange_point_count == 6);
  assert(u.exposed_regions[0].child_i_begin == 4);
  assert(u.exposed_regions[0].child_j_begin == 0);
  assert(u.exposed_regions[0].extent_i == 1);
  assert(u.exposed_regions[0].extent_j == 3);
  assert(!u.exposed_regions[0].owns_overlap);
  assert(!u.exposed_regions[0].owns_halo);

  const auto& v = exchange.fields[1];
  assert(v.field == tywrf::nest::StateExchangeField::v);
  assert(v.active_nx == 4);
  assert(v.active_ny == 4);
  assert(v.active_k_count == 2);
  assert(v.exposed_region_count == 1);
  assert(v.exposed_horizontal_cell_count == 4);
  assert(v.exchange_point_count == 8);

  const auto& mu = exchange.fields[2];
  assert(mu.field == tywrf::nest::StateExchangeField::mu);
  assert(mu.stagger == tywrf::nest::HorizontalStagger::surface);
  assert(!mu.three_dimensional);
  assert(mu.active_k_count == 1);
  assert(mu.exposed_horizontal_cell_count == 3);
  assert(mu.exchange_point_count == 3);

  const auto& qvapor = exchange.fields[3];
  assert(qvapor.field == tywrf::nest::StateExchangeField::qvapor);
  assert(qvapor.stagger == tywrf::nest::HorizontalStagger::mass);
  assert(qvapor.active_nx == 4);
  assert(qvapor.active_ny == 3);
  assert(qvapor.active_k_count == 2);
  assert(qvapor.exposed_horizontal_cell_count == 3);
  assert(qvapor.exchange_point_count == 6);

  const auto& t = exchange.fields[4];
  assert(t.field == tywrf::nest::StateExchangeField::t);
  assert(t.stagger == tywrf::nest::HorizontalStagger::mass);
  assert(t.three_dimensional);
  assert(t.active_nx == 4);
  assert(t.active_ny == 3);
  assert(t.active_k_count == 2);
  assert(t.exposed_horizontal_cell_count == 3);
  assert(t.exchange_point_count == 6);

  const auto& ph = exchange.fields[5];
  assert(ph.field == tywrf::nest::StateExchangeField::ph);
  assert(ph.stagger == tywrf::nest::HorizontalStagger::w_full);
  assert(ph.three_dimensional);
  assert(ph.active_nx == 4);
  assert(ph.active_ny == 3);
  assert(ph.active_k_count == 3);
  assert(ph.exposed_horizontal_cell_count == 3);
  assert(ph.exchange_point_count == 9);

  expect_unchanged(child.u, u_before, "U");
  expect_unchanged(child.v, v_before, "V");
  expect_unchanged(child.mu, mu_before, "MU");
  expect_unchanged(child.qvapor, qvapor_before, "QVAPOR");
  expect_unchanged(child.t, t_before, "T");
  expect_unchanged(child.ph, ph_before, "PH");
}

void test_full_overlap_has_no_exchange_work() {
  tywrf::State<float> child(make_child_grid());
  fill_state(child);

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {1, 1, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  const auto exchange =
      tywrf::nest::build_exposed_child_state_exchange_plan(
          remap, static_cast<const tywrf::State<float>&>(child).view());
  assert(exchange.ok());
  assert(exchange.field_count == 6);
  assert(exchange.report.planned_field_count == 6);
  assert(exchange.report.active_field_count == 0);
  assert(exchange.report.exposed_region_count == 0);
  assert(exchange.report.exposed_horizontal_cell_count == 0);
  assert(exchange.report.exchange_point_count == 0);
  assert(!exchange.report.requires_parent_interpolation);
  assert(!exchange.report.performed_interpolation);
  assert(!exchange.report.modifies_overlap);
  assert(!exchange.report.modifies_halo);
  for (std::uint8_t field_index = 0; field_index < exchange.field_count; ++field_index) {
    assert(exchange.fields[field_index].exposed_region_count == 0);
    assert(exchange.fields[field_index].exchange_point_count == 0);
  }
}

void test_diagonal_move_decomposes_exposed_cells_without_overlap_or_halo() {
  tywrf::State<float> child(make_child_grid());
  fill_state(child);

  const auto remap = make_remap_plan(
      {1, 1, tywrf::nest::IndexBase::one_based},
      {2, 2, tywrf::nest::IndexBase::one_based});
  assert(remap.ok());

  const auto exchange =
      tywrf::nest::build_exposed_child_state_exchange_plan(
          remap, static_cast<const tywrf::State<float>&>(child).view());
  assert(exchange.ok());
  assert(exchange.report.active_field_count == 6);
  assert(exchange.report.exposed_region_count == 12);
  assert(exchange.report.exposed_horizontal_cell_count == 38);
  assert(exchange.report.exchange_point_count == 76);

  const auto& qvapor = exchange.fields[3];
  assert(qvapor.exposed_region_count == 2);
  assert(qvapor.exposed_regions[0].child_i_begin == 0);
  assert(qvapor.exposed_regions[0].child_j_begin == 2);
  assert(qvapor.exposed_regions[0].extent_i == 4);
  assert(qvapor.exposed_regions[0].extent_j == 1);
  assert(qvapor.exposed_regions[1].child_i_begin == 3);
  assert(qvapor.exposed_regions[1].child_j_begin == 0);
  assert(qvapor.exposed_regions[1].extent_i == 1);
  assert(qvapor.exposed_regions[1].extent_j == 2);
  assert(!qvapor.owns_overlap);
  assert(!qvapor.owns_halo);

  const auto& t = exchange.fields[4];
  assert(t.exposed_region_count == 2);
  assert(t.exposed_horizontal_cell_count == qvapor.exposed_horizontal_cell_count);
  assert(t.exchange_point_count == qvapor.exchange_point_count);

  const auto& ph = exchange.fields[5];
  assert(ph.stagger == tywrf::nest::HorizontalStagger::w_full);
  assert(ph.exposed_region_count == 2);
  assert(ph.exposed_horizontal_cell_count == qvapor.exposed_horizontal_cell_count);
  assert(ph.exchange_point_count == 18);
}

}  // namespace

int main() {
  test_exposed_strip_accounting();
  test_full_overlap_has_no_exchange_work();
  test_diagonal_move_decomposes_exposed_cells_without_overlap_or_halo();
  return 0;
}
