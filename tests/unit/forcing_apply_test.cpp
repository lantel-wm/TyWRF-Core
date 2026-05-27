#include "tywrf/io/forcing_apply.hpp"

#include "tywrf/grid.hpp"
#include "tywrf/io/forcing_frames.hpp"
#include "tywrf/state.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

int failures = 0;

void expect(const bool condition, const std::string_view message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

void expect_close(
    const float actual,
    const float expected,
    const std::string_view message) {
  expect(std::abs(actual - expected) < 1.0e-6F, message);
}

template <typename Fn>
void expect_apply_error(
    Fn&& fn,
    const tywrf::io::ForcingApplyStatus expected_status,
    const std::string_view message) {
  try {
    fn();
  } catch (const tywrf::io::ForcingApplyError& error) {
    expect(error.status() == expected_status, message);
    return;
  }
  expect(false, message);
}

tywrf::Grid make_grid() {
  return tywrf::Grid({
      .mass_nx = 4,
      .mass_ny = 3,
      .mass_nz = 2,
      .full_nz = 3,
      .halo = {1, 1, 1, 1, 0, 0},
  });
}

std::vector<float> snapshot_storage(const tywrf::FieldStorage3D<float>& storage) {
  return std::vector<float>(storage.data(), storage.data() + storage.size());
}

void zero_state(tywrf::State<float>& state) {
  auto zero_storage = [](auto& storage) {
    for (std::size_t index = 0; index < storage.size(); ++index) {
      storage.data()[index] = 0.0F;
    }
  };

  zero_storage(state.u);
  zero_storage(state.v);
  zero_storage(state.w);
  zero_storage(state.ph);
  zero_storage(state.phb);
  zero_storage(state.t);
  zero_storage(state.p);
  zero_storage(state.pb);
  zero_storage(state.qvapor);
  zero_storage(state.qcloud);
  zero_storage(state.qrain);
  zero_storage(state.qice);
  zero_storage(state.qsnow);
  zero_storage(state.qgraup);
  zero_storage(state.qnice);
  zero_storage(state.qnrain);
  zero_storage(state.mu);
  zero_storage(state.mub);
  zero_storage(state.psfc);
  zero_storage(state.u10);
  zero_storage(state.v10);
  zero_storage(state.t2);
  zero_storage(state.q2);
  zero_storage(state.rainc);
  zero_storage(state.rainnc);
}

tywrf::io::PackedForcingField make_i_side_t_boundary_pack() {
  constexpr std::size_t bdy_width = 2;
  constexpr std::size_t nz = 2;
  constexpr std::size_t ny = 3;
  std::vector<float> raw(bdy_width * nz * ny);
  for (std::size_t bdy = 0; bdy < bdy_width; ++bdy) {
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t j = 0; j < ny; ++j) {
        raw[((bdy * nz) + k) * ny + j] =
            static_cast<float>(100 * bdy + 10 * k + j);
      }
    }
  }

  const std::vector<std::size_t> shape{bdy_width, nz, ny};
  return tywrf::io::pack_boundary_x_side_raw_to_canonical(raw, shape);
}

void test_boundary_validation_reports_shape_and_does_not_write_state() {
  tywrf::State<float> state(make_grid());
  zero_state(state);
  const auto t_before = snapshot_storage(state.t);

  const auto packed = make_i_side_t_boundary_pack();
  const auto report = tywrf::io::validate_boundary_pack_for_state(
      state,
      "T",
      tywrf::io::BoundarySide::i_lower,
      packed);

  expect(report.status == tywrf::io::ForcingApplyStatus::ok,
         "boundary validation status ok");
  expect(report.field_count == 1, "boundary validation reports one field");
  expect(report.point_count == 12, "boundary validation reports packed point count");
  expect(!report.would_modify_state, "boundary validation is read-only");
  expect(!report.synthetic, "boundary validation is not synthetic apply");
  expect(snapshot_storage(state.t) == t_before, "boundary validation does not write T");
}

void test_boundary_validation_rejects_wrong_shape() {
  tywrf::State<float> state(make_grid());
  zero_state(state);

  tywrf::io::PackedForcingField wrong;
  wrong.layout = tywrf::io::make_canonical_forcing_layout(2, 4, 2);
  wrong.values.assign(wrong.layout.value_count, 1.0F);

  expect_apply_error(
      [&] {
        (void)tywrf::io::validate_boundary_pack_for_state(
            state,
            "T",
            tywrf::io::BoundarySide::i_lower,
            wrong);
      },
      tywrf::io::ForcingApplyStatus::shape_mismatch,
      "wrong boundary shape raises shape_mismatch");
}

void test_synthetic_nudging_delta_writes_only_requested_window() {
  tywrf::State<float> state(make_grid());
  zero_state(state);

  tywrf::io::PackedForcingField delta;
  delta.layout = tywrf::io::make_canonical_forcing_layout(2, 1, 1);
  delta.values = {0.5F, 1.5F};

  const auto report = tywrf::io::apply_synthetic_nudging_delta(
      state,
      "T",
      delta,
      tywrf::io::SyntheticNudgingDeltaConfig{
          .active_i_begin = 1,
          .active_j_begin = 1,
          .active_k_begin = 0,
      });

  expect(report.status == tywrf::io::ForcingApplyStatus::ok,
         "synthetic delta status ok");
  expect(report.field_count == 1, "synthetic delta reports one field");
  expect(report.point_count == 2, "synthetic delta reports written points");
  expect(report.would_modify_state, "synthetic delta reports state write");
  expect(report.synthetic, "synthetic delta metadata is explicit");

  const auto view = state.view().t;
  expect_close(view(2, 2, 0), 0.5F, "synthetic delta first target point");
  expect_close(view(3, 2, 0), 1.5F, "synthetic delta second target point");

  std::size_t changed = 0;
  for (std::size_t index = 0; index < state.t.size(); ++index) {
    if (state.t.data()[index] != 0.0F) {
      ++changed;
    }
  }
  expect(changed == 2, "synthetic delta only changes requested T points");

  for (std::size_t index = 0; index < state.u.size(); ++index) {
    expect(state.u.data()[index] == 0.0F, "synthetic T delta does not change U");
  }
}

void test_synthetic_nudging_delta_rejects_out_of_range_window() {
  tywrf::State<float> state(make_grid());
  zero_state(state);

  tywrf::io::PackedForcingField delta;
  delta.layout = tywrf::io::make_canonical_forcing_layout(3, 1, 1);
  delta.values.assign(delta.layout.value_count, 1.0F);

  expect_apply_error(
      [&] {
        (void)tywrf::io::apply_synthetic_nudging_delta(
            state,
            "T",
            delta,
            tywrf::io::SyntheticNudgingDeltaConfig{
                .active_i_begin = 2,
                .active_j_begin = 0,
                .active_k_begin = 0,
            });
      },
      tywrf::io::ForcingApplyStatus::invalid_range,
      "out-of-range synthetic window raises invalid_range");
}

}  // namespace

int main() {
  test_boundary_validation_reports_shape_and_does_not_write_state();
  test_boundary_validation_rejects_wrong_shape();
  test_synthetic_nudging_delta_writes_only_requested_window();
  test_synthetic_nudging_delta_rejects_out_of_range_window();

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated forcing apply skeleton contracts\n";
  return 0;
}
