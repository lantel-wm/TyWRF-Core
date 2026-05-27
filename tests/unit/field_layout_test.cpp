#include "tywrf/field_view.hpp"
#include "tywrf/grid.hpp"
#include "tywrf/state.hpp"

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace {

tywrf::Grid make_test_grid() {
  return tywrf::Grid({
      .mass_nx = 4,
      .mass_ny = 3,
      .mass_nz = 2,
      .full_nz = 3,
      .halo = {1, 2, 1, 1, 1, 0},
  });
}

}  // namespace

int main() {
  static_assert(std::is_standard_layout_v<tywrf::FieldView3D<double>>);
  static_assert(std::is_trivially_copyable_v<tywrf::FieldView3D<double>>);
  static_assert(std::is_standard_layout_v<tywrf::FieldView2D<double>>);
  static_assert(std::is_trivially_copyable_v<tywrf::FieldView2D<double>>);
  static_assert(std::is_standard_layout_v<tywrf::StateView<double>>);
  static_assert(std::is_trivially_copyable_v<tywrf::StateView<double>>);

  const auto grid = make_test_grid();
  const auto mass = grid.mass_layout();

  assert(mass.nx == 7);
  assert(mass.ny == 5);
  assert(mass.nz == 3);
  assert(mass.active_nx() == 4);
  assert(mass.active_ny() == 3);
  assert(mass.active_nz() == 2);
  assert(mass.i_begin() == 1);
  assert(mass.i_end() == 5);
  assert(mass.j_begin() == 1);
  assert(mass.j_end() == 4);
  assert(mass.k_begin() == 1);
  assert(mass.k_end() == 3);

  assert(mass.stride_i == 1);
  assert(mass.stride_k == mass.nx);
  assert(mass.stride_j == mass.nx * mass.nz);

  const std::size_t expected_index = ((2U * static_cast<std::size_t>(mass.nz)) + 1U) *
                                         static_cast<std::size_t>(mass.nx) +
                                     3U;
  assert(mass.index(3, 2, 1) == expected_index);
  assert(tywrf::canonical_index_3d(3, 2, 1, mass.nx, mass.nz) == expected_index);

  tywrf::FieldStorage3D<double> field(mass);
  auto view = field.view();
  view(3, 2, 1) = 42.0;
  view(4, 2, 1) = 43.0;
  assert(field.data()[expected_index] == 42.0);
  assert(view.index(4, 2, 1) == view.index(3, 2, 1) + 1);
  assert(field.data()[view.index(4, 2, 1)] == 43.0);

  const auto u_layout = grid.u_layout();
  const auto v_layout = grid.v_layout();
  const auto w_layout = grid.w_layout();
  const auto surface_layout = grid.surface_layout();
  assert(u_layout.active_nx() == mass.active_nx() + 1);
  assert(v_layout.active_ny() == mass.active_ny() + 1);
  assert(w_layout.active_nz() == 3);
  assert(surface_layout.active_nx() == mass.active_nx());
  assert(surface_layout.active_ny() == mass.active_ny());
  assert(surface_layout.stride_i == 1);
  assert(surface_layout.stride_j == surface_layout.nx);

  tywrf::State<double> state(grid);
  assert(state.u.data() != state.v.data());
  assert(state.t.data() != state.qvapor.data());
  assert(state.mu.data() != state.mub.data());
  assert(state.u.size() == u_layout.allocation_size());
  assert(state.v.size() == v_layout.allocation_size());
  assert(state.w.size() == w_layout.allocation_size());
  assert(state.t.size() == mass.allocation_size());
  assert(state.mu.size() == surface_layout.allocation_size());

  auto state_view = state.view();
  state_view.t(mass.i_begin(), mass.j_begin(), mass.k_begin()) = 300.0;
  state_view.mu(surface_layout.i_begin(), surface_layout.j_begin()) = 12.0;
  assert(state.t.data()[mass.index(mass.i_begin(), mass.j_begin(), mass.k_begin())] == 300.0);
  assert(state.mu.data()[surface_layout.index(surface_layout.i_begin(), surface_layout.j_begin())] ==
         12.0);

  return 0;
}
