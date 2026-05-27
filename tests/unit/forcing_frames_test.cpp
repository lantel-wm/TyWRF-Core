#include "tywrf/io/forcing_frames.hpp"
#include "tywrf/io/forcing_io.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
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
    const double actual,
    const double expected,
    const std::string_view message) {
  expect(std::abs(actual - expected) < 1.0e-12, message);
}

template <typename Fn>
void expect_frame_error(Fn&& fn, const std::string_view message) {
  try {
    fn();
  } catch (const tywrf::io::ForcingFrameError&) {
    return;
  }
  expect(false, message);
}

std::filesystem::path reference_dir() {
  if (const char* env = std::getenv("TYWRF_REFERENCE_DIR")) {
    return env;
  }
  return "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF";
}

std::vector<std::string> make_synthetic_krosa_times() {
  std::vector<std::string> times;
  times.reserve(28);
  for (int index = 0; index < 28; ++index) {
    times.push_back(tywrf::io::add_seconds_to_wrf_timestamp(
        "2025-07-26_00:00:00",
        static_cast<std::int64_t>(index) *
            tywrf::io::kKrosaForcingIntervalSeconds));
  }
  return times;
}

int run_selector_test() {
  const auto times = make_synthetic_krosa_times();
  const tywrf::io::ForcingFrameSelector selector(
      tywrf::io::ForcingFrameSelectorConfig{
          28,
          tywrf::io::kKrosaForcingIntervalSeconds,
          times});

  const auto first = selector.frame_for_cycle_index(0);
  expect(first.record_index == 0, "cycle 0 selects record 0");
  expect(first.start_seconds == 0, "cycle 0 starts at 0 s");
  expect(first.end_seconds == 21'600, "cycle 0 ends at 6 h");
  expect(first.start_time == "2025-07-26_00:00:00", "cycle 0 start timestamp");
  expect(first.end_time == "2025-07-26_06:00:00", "cycle 0 end timestamp");

  const auto ordinary = selector.frame_for_start_seconds(43'200);
  expect(ordinary.record_index == 2, "ordinary start_seconds selects record 2");
  expect(ordinary.start_time == "2025-07-26_12:00:00",
         "ordinary start_seconds start timestamp");
  expect(ordinary.end_time == "2025-07-26_18:00:00",
         "ordinary start_seconds end timestamp");

  const auto final = selector.frame_for_cycle_index(27);
  expect(final.record_index == 27, "final cycle selects record 27");
  expect(final.start_seconds == 27 * 21'600, "final cycle start seconds");
  expect(final.end_seconds == 28 * 21'600, "final cycle end seconds");
  expect(final.start_time == "2025-08-01_18:00:00", "final cycle start timestamp");
  expect(final.end_time == "2025-08-02_00:00:00",
         "final cycle end timestamp is inferred without Times[28]");

  const auto start_weights = selector.weights_for(first, first.start_seconds);
  expect_close(start_weights.old_weight, 1.0, "start old weight");
  expect_close(start_weights.new_weight, 0.0, "start new weight");
  const auto mid_weights = selector.weights_for(first, 10'800);
  expect_close(mid_weights.old_weight, 0.5, "midpoint old weight");
  expect_close(mid_weights.new_weight, 0.5, "midpoint new weight");
  const auto end_weights = selector.weights_for(first, first.end_seconds);
  expect_close(end_weights.old_weight, 0.0, "endpoint old weight");
  expect_close(end_weights.new_weight, 1.0, "endpoint new weight");
  expect(selector.frame_for_start_seconds(first.end_seconds).record_index == 1,
         "endpoint record transition is explicit, not guessed by weights_for");

  expect_frame_error(
      [&] { (void)selector.frame_for_cycle_index(28); },
      "cycle index 28 should be out of range");
  expect_frame_error(
      [&] { (void)selector.frame_for_start_seconds(28 * 21'600); },
      "168 h start_seconds should be out of range for 28 records");
  expect_frame_error(
      [&] { (void)selector.frame_for_start_seconds(1); },
      "unaligned start_seconds should be rejected");
  expect_frame_error(
      [&] { (void)selector.weights_for(first, -1); },
      "weights before frame should be rejected");
  expect_frame_error(
      [&] { (void)selector.weights_for(first, first.end_seconds + 1); },
      "weights after frame should be rejected");

  auto invalid_times = times;
  invalid_times.pop_back();
  expect_frame_error(
      [&] {
        (void)tywrf::io::ForcingFrameSelector(
            tywrf::io::ForcingFrameSelectorConfig{28, 21'600, invalid_times});
      },
      "selector should reject timestamp count mismatch");

  return failures == 0 ? 0 : 1;
}

int run_staging_test() {
  const std::array<std::size_t, 3> shape{2, 3, 5};
  const auto layout = tywrf::io::make_forcing_field_layout(shape);
  expect(layout.rank == 3, "staging layout rank");
  expect(layout.extents[0] == 2 && layout.extents[1] == 3 && layout.extents[2] == 5,
         "staging layout extents");
  expect(layout.strides[0] == 15 && layout.strides[1] == 5 && layout.strides[2] == 1,
         "staging layout raw row-major strides");
  expect(layout.value_count == 30, "staging layout value count");

  auto buffer = tywrf::io::make_forcing_field_buffer(layout);
  expect(buffer.values.size() == 30, "staging buffer allocates one slice");
  auto view = buffer.view();
  const std::array<std::size_t, 3> indices{1, 2, 4};
  view.at(indices) = 29.0F;
  expect(buffer.values.back() == 29.0F, "staging view indexes last value");

  tywrf::io::ForcingTimeSlice slice;
  slice.metadata.name = "U_NDG_NEW";
  slice.metadata.slice_shape.assign(shape.begin(), shape.end());
  slice.time_index = 3;
  slice.values.resize(30);
  for (std::size_t index = 0; index < slice.values.size(); ++index) {
    slice.values[index] = static_cast<float>(index);
  }

  const auto staged = tywrf::io::stage_forcing_time_slice(slice);
  expect(staged.variable_name == "U_NDG_NEW", "staged variable name");
  expect(staged.record_index == 3, "staged record index");
  expect(staged.field.layout.value_count == 30, "staged layout value count");
  expect(staged.field.values.back() == 29.0F, "staged values copy one slice");

  const auto mappings = tywrf::io::krosa_spectral_nudging_variable_mappings();
  expect(mappings.size() == 6, "KROSA nudging mapping count");
  expect(mappings[3].canonical_name == "QVAPOR" &&
             mappings[3].old_variable_name == "Q_NDG_OLD" &&
             mappings[3].new_variable_name == "Q_NDG_NEW",
         "KROSA QVAPOR nudging mapping");

  expect_frame_error(
      [&] {
        const std::array<std::size_t, 5> too_many{1, 1, 1, 1, 1};
        (void)tywrf::io::make_forcing_field_layout(too_many);
      },
      "staging layout should reject ranks above four");

  return failures == 0 ? 0 : 1;
}

int run_packing_test() {
  {
    const std::array<std::size_t, 3> shape{2, 3, 4};
    std::vector<float> raw(shape[0] * shape[1] * shape[2]);
    for (std::size_t k = 0; k < shape[0]; ++k) {
      for (std::size_t j = 0; j < shape[1]; ++j) {
        for (std::size_t i = 0; i < shape[2]; ++i) {
          raw[((k * shape[1]) + j) * shape[2] + i] =
              static_cast<float>(100 * k + 10 * j + i);
        }
      }
    }

    const auto packed = tywrf::io::pack_fdda_3d_raw_to_canonical(raw, shape);
    expect(packed.layout.nx == 4 && packed.layout.ny == 3 &&
               packed.layout.nz == 2,
           "FDDA packed layout extents");
    for (std::size_t k = 0; k < shape[0]; ++k) {
      for (std::size_t j = 0; j < shape[1]; ++j) {
        for (std::size_t i = 0; i < shape[2]; ++i) {
          const auto expected = static_cast<float>(100 * k + 10 * j + i);
          const auto canonical_index = ((j * shape[0]) + k) * shape[2] + i;
          expect(packed.values.at(canonical_index) == expected,
                 "FDDA raw k,j,i packs to canonical j,k,i");
          expect(packed.at(i, j, k) == expected, "FDDA packed at()");
        }
      }
    }
  }

  {
    const std::array<std::size_t, 3> shape{2, 3, 4};
    std::vector<float> raw(shape[0] * shape[1] * shape[2]);
    for (std::size_t bdy = 0; bdy < shape[0]; ++bdy) {
      for (std::size_t k = 0; k < shape[1]; ++k) {
        for (std::size_t j = 0; j < shape[2]; ++j) {
          raw[((bdy * shape[1]) + k) * shape[2] + j] =
              static_cast<float>(100 * bdy + 10 * k + j);
        }
      }
    }

    const auto packed =
        tywrf::io::pack_boundary_x_side_raw_to_canonical(raw, shape);
    expect(packed.layout.nx == 2 && packed.layout.ny == 4 &&
               packed.layout.nz == 3,
           "X-side boundary packed layout extents");
    for (std::size_t bdy = 0; bdy < shape[0]; ++bdy) {
      for (std::size_t k = 0; k < shape[1]; ++k) {
        for (std::size_t j = 0; j < shape[2]; ++j) {
          const auto expected = static_cast<float>(100 * bdy + 10 * k + j);
          const auto canonical_index = ((j * shape[1]) + k) * shape[0] + bdy;
          expect(packed.values.at(canonical_index) == expected,
                 "X-side raw bdy,k,j packs to canonical i=bdy,j,k");
          expect(packed.at(bdy, j, k) == expected, "X-side packed at()");
        }
      }
    }
  }

  {
    const std::array<std::size_t, 3> shape{3, 2, 4};
    std::vector<float> raw(shape[0] * shape[1] * shape[2]);
    for (std::size_t bdy = 0; bdy < shape[0]; ++bdy) {
      for (std::size_t k = 0; k < shape[1]; ++k) {
        for (std::size_t i = 0; i < shape[2]; ++i) {
          raw[((bdy * shape[1]) + k) * shape[2] + i] =
              static_cast<float>(100 * bdy + 10 * k + i);
        }
      }
    }

    const auto packed =
        tywrf::io::pack_boundary_y_side_raw_to_canonical(raw, shape);
    expect(packed.layout.nx == 4 && packed.layout.ny == 3 &&
               packed.layout.nz == 2,
           "Y-side boundary packed layout extents");
    for (std::size_t bdy = 0; bdy < shape[0]; ++bdy) {
      for (std::size_t k = 0; k < shape[1]; ++k) {
        for (std::size_t i = 0; i < shape[2]; ++i) {
          const auto expected = static_cast<float>(100 * bdy + 10 * k + i);
          const auto canonical_index = ((bdy * shape[1]) + k) * shape[2] + i;
          expect(packed.values.at(canonical_index) == expected,
                 "Y-side raw bdy,k,i packs to canonical i,j=bdy,k");
          expect(packed.at(i, bdy, k) == expected, "Y-side packed at()");
        }
      }
    }
  }

  expect_frame_error(
      [] {
        const std::array<std::size_t, 2> shape{2, 3};
        const std::vector<float> raw(shape[0] * shape[1], 1.0F);
        (void)tywrf::io::pack_fdda_3d_raw_to_canonical(raw, shape);
      },
      "raw packing should reject non-3D rank");
  expect_frame_error(
      [] {
        const std::array<std::size_t, 3> shape{2, 0, 3};
        const std::vector<float> raw(1, 1.0F);
        (void)tywrf::io::pack_boundary_x_side_raw_to_canonical(raw, shape);
      },
      "raw packing should reject empty extents");
  expect_frame_error(
      [] {
        const std::array<std::size_t, 3> shape{2, 3, 4};
        const std::vector<float> raw(23, 1.0F);
        (void)tywrf::io::pack_boundary_y_side_raw_to_canonical(raw, shape);
      },
      "raw packing should reject value count mismatches");

  return failures == 0 ? 0 : 1;
}

int run_reference_smoke_test() {
  const auto root = reference_dir();
  const auto wrfbdy = root / "wrfbdy_d01";
  const auto wrffdda = root / "wrffdda_d01";
  if (!std::filesystem::exists(wrfbdy) || !std::filesystem::exists(wrffdda)) {
    std::cout << "Skipping forcing frame reference smoke; KROSA inputs are not present under "
              << root << '\n';
    return 0;
  }

  const tywrf::io::KrosaForcingReader boundary(
      wrfbdy,
      tywrf::io::KrosaForcingKind::boundary);
  const auto boundary_selector = tywrf::io::make_krosa_forcing_frame_selector(boundary);
  expect(boundary.time_count() == 28, "KROSA wrfbdy Time dimension is 28");
  expect(boundary_selector.time_count() == 28, "KROSA wrfbdy selector has 28 records");
  const auto boundary_record0 = boundary_selector.frame_for_cycle_index(0);
  const auto boundary_record27 = boundary_selector.frame_for_cycle_index(27);
  expect(boundary_record0.record_index == 0 && boundary_record0.start_seconds == 0 &&
             boundary_record0.end_seconds == 21'600,
         "KROSA wrfbdy record 0 interval");
  expect(boundary_record0.start_time == boundary.read_time_string(0) &&
             boundary_record0.end_time == boundary.read_time_string(1),
         "KROSA wrfbdy record 0 timestamps");
  expect(boundary_record27.record_index == 27 &&
             boundary_record27.start_seconds == 27 * 21'600 &&
             boundary_record27.end_seconds == 28 * 21'600,
         "KROSA wrfbdy record 27 interval");
  expect(boundary_record27.start_time == boundary.read_time_string(27),
         "KROSA wrfbdy record 27 start timestamp");
  expect(boundary_record27.end_time ==
             tywrf::io::add_seconds_to_wrf_timestamp(
                 boundary.read_time_string(27),
                 tywrf::io::kKrosaForcingIntervalSeconds),
         "KROSA wrfbdy record 27 inferred end timestamp");

  const auto u_bxs_slice = boundary.read_float_time_slice("U_BXS", 0);
  expect(u_bxs_slice.metadata.slice_shape ==
             std::vector<std::size_t>{5, 59, 429},
         "KROSA wrfbdy U_BXS record 0 raw shape");
  const auto staged_u_bxs = tywrf::io::stage_forcing_time_slice(u_bxs_slice);
  expect(staged_u_bxs.variable_name == "U_BXS",
         "KROSA wrfbdy U_BXS staged variable name");
  expect(staged_u_bxs.record_index == 0,
         "KROSA wrfbdy U_BXS staged record index");
  const auto packed_u_bxs = tywrf::io::pack_boundary_x_side_raw_to_canonical(
      staged_u_bxs.field.values,
      u_bxs_slice.metadata.slice_shape);
  expect(packed_u_bxs.layout.nx == 5 && packed_u_bxs.layout.ny == 429 &&
             packed_u_bxs.layout.nz == 59,
         "KROSA wrfbdy U_BXS canonical shape");
  const std::array<std::array<std::size_t, 3>, 4> u_bxs_samples{{
      {0, 0, 0},
      {4, 428, 58},
      {2, 17, 31},
      {3, 287, 7},
  }};
  for (const auto& sample : u_bxs_samples) {
    const auto i = sample[0];
    const auto j = sample[1];
    const auto k = sample[2];
    const auto raw_index = ((i * 59) + k) * 429 + j;
    expect(packed_u_bxs.at(i, j, k) == staged_u_bxs.field.values.at(raw_index),
           "KROSA wrfbdy U_BXS raw bdy,k,j maps to canonical i,j,k");
  }

  const tywrf::io::KrosaForcingReader nudging(
      wrffdda,
      tywrf::io::KrosaForcingKind::spectral_nudging);
  const auto nudging_selector = tywrf::io::make_krosa_forcing_frame_selector(nudging);
  expect(nudging.time_count() == 28, "KROSA wrffdda Time dimension is 28");
  expect(nudging_selector.time_count() == 28, "KROSA wrffdda selector has 28 records");
  const auto nudging_record0 = nudging_selector.frame_for_cycle_index(0);
  const auto nudging_record27 = nudging_selector.frame_for_start_seconds(27 * 21'600);
  expect(nudging_record0.start_seconds == 0 && nudging_record0.end_seconds == 21'600,
         "KROSA wrffdda record 0 interval");
  expect(nudging_record27.record_index == 27 &&
             nudging_record27.start_seconds == 27 * 21'600 &&
             nudging_record27.end_seconds == 28 * 21'600,
         "KROSA wrffdda record 27 interval");
  expect(nudging_record27.end_time ==
             tywrf::io::add_seconds_to_wrf_timestamp(
                 nudging.read_time_string(27),
                 tywrf::io::kKrosaForcingIntervalSeconds),
         "KROSA wrffdda record 27 inferred end timestamp");

  const auto q_new_slice = nudging.read_float_time_slice("Q_NDG_NEW", 0);
  expect(q_new_slice.metadata.slice_shape ==
             std::vector<std::size_t>{59, 429, 265},
         "KROSA wrffdda Q_NDG_NEW record 0 raw shape");
  const auto packed_q_new = tywrf::io::pack_fdda_3d_raw_to_canonical(
      q_new_slice.values,
      q_new_slice.metadata.slice_shape);
  expect(packed_q_new.layout.nx == 265 && packed_q_new.layout.ny == 429 &&
             packed_q_new.layout.nz == 59,
         "KROSA wrffdda Q_NDG_NEW canonical shape");
  const std::array<std::array<std::size_t, 3>, 4> q_new_samples{{
      {0, 0, 0},
      {264, 428, 58},
      {19, 23, 41},
      {251, 307, 6},
  }};
  for (const auto& sample : q_new_samples) {
    const auto i = sample[0];
    const auto j = sample[1];
    const auto k = sample[2];
    const auto raw_index = ((k * 429) + j) * 265 + i;
    expect(packed_q_new.at(i, j, k) == q_new_slice.values.at(raw_index),
           "KROSA wrffdda Q_NDG_NEW raw k,j,i maps to canonical i,j,k");
  }

  std::cout << "Validated KROSA forcing frame metadata and canonical pack smoke under "
            << root << '\n';
  return failures == 0 ? 0 : 1;
}

}  // namespace

int main() {
  try {
    if (const int status = run_selector_test(); status != 0) {
      return status;
    }
    if (const int status = run_staging_test(); status != 0) {
      return status;
    }
    if (const int status = run_packing_test(); status != 0) {
      return status;
    }
    if (const int status = run_reference_smoke_test(); status != 0) {
      return status;
    }
  } catch (const std::exception& error) {
    std::cerr << "forcing frame test failed: " << error.what() << '\n';
    return 1;
  }
  return failures == 0 ? 0 : 1;
}
