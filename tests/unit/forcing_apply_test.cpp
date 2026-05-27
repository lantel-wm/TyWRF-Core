#include "tywrf/io/forcing_apply.hpp"

#include "tywrf/grid.hpp"
#include "tywrf/io/forcing_frames.hpp"
#include "tywrf/io/forcing_io.hpp"
#include "tywrf/io/wrf_state_io.hpp"
#include "tywrf/state.hpp"

#include <netcdf.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
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
    const float actual,
    const float expected,
    const std::string_view message) {
  expect(std::abs(actual - expected) < 1.0e-6F, message);
}

void check_nc(const int status, const std::string_view message) {
  if (status != NC_NOERR) {
    throw std::runtime_error(
        std::string(message) + ": " + nc_strerror(status));
  }
}

std::filesystem::path reference_dir() {
  if (const char* env = std::getenv("TYWRF_REFERENCE_DIR")) {
    return env;
  }
  return "/home/zzy/Projects/tc_sim/pgwrf_2025wp12_d0110km/PGWRF/output_gfs_analysis/2025wp12/2025072600/WRF";
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

void fill_storage(tywrf::FieldStorage3D<float>& storage, const float value) {
  for (std::size_t index = 0; index < storage.size(); ++index) {
    storage.data()[index] = value;
  }
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

std::vector<float> make_raw_boundary_values(
    const std::size_t bdy_width,
    const std::size_t nz,
    const std::size_t horizontal_extent,
    const float base) {
  std::vector<float> raw(bdy_width * nz * horizontal_extent);
  for (std::size_t bdy = 0; bdy < bdy_width; ++bdy) {
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t h = 0; h < horizontal_extent; ++h) {
        raw[((bdy * nz) + k) * horizontal_extent + h] =
            base + static_cast<float>(100 * bdy + 10 * k + h);
      }
    }
  }
  return raw;
}

int define_float_var(
    const int file_id,
    const char* name,
    const std::array<int, 4> dimensions) {
  int var_id = -1;
  check_nc(
      nc_def_var(file_id, name, NC_FLOAT, 4, dimensions.data(), &var_id),
      std::string("define ") + name);
  return var_id;
}

void create_synthetic_u_wrfbdy_for_apply(const std::filesystem::path& path) {
  int file_id = -1;
  check_nc(nc_create(path.string().c_str(), NC_CLOBBER, &file_id),
           "create synthetic wrfbdy");

  int time_dim = -1;
  int date_dim = -1;
  int bdy_width_dim = -1;
  int bottom_top_dim = -1;
  int south_north_dim = -1;
  int west_east_stag_dim = -1;
  check_nc(nc_def_dim(file_id, "Time", 1, &time_dim), "define Time");
  check_nc(nc_def_dim(file_id, "DateStrLen", 19, &date_dim),
           "define DateStrLen");
  check_nc(nc_def_dim(file_id, "bdy_width", 1, &bdy_width_dim),
           "define bdy_width");
  check_nc(nc_def_dim(file_id, "bottom_top", 2, &bottom_top_dim),
           "define bottom_top");
  check_nc(nc_def_dim(file_id, "south_north", 3, &south_north_dim),
           "define south_north");
  check_nc(nc_def_dim(file_id, "west_east_stag", 5, &west_east_stag_dim),
           "define west_east_stag");

  const std::array<int, 2> time_dims{time_dim, date_dim};
  int times_var = -1;
  check_nc(nc_def_var(file_id, "Times", NC_CHAR, 2, time_dims.data(), &times_var),
           "define Times");
  const std::array<int, 4> x_dims{
      time_dim,
      bdy_width_dim,
      bottom_top_dim,
      south_north_dim,
  };
  const std::array<int, 4> y_dims{
      time_dim,
      bdy_width_dim,
      bottom_top_dim,
      west_east_stag_dim,
  };
  const int u_bxs = define_float_var(file_id, "U_BXS", x_dims);
  const int u_bxe = define_float_var(file_id, "U_BXE", x_dims);
  const int u_bys = define_float_var(file_id, "U_BYS", y_dims);
  const int u_bye = define_float_var(file_id, "U_BYE", y_dims);
  check_nc(nc_enddef(file_id), "end synthetic wrfbdy definitions");

  const char timestamp[] = "2025-07-26_00:00:00";
  const std::size_t time_start[2] = {0, 0};
  const std::size_t time_count[2] = {1, 19};
  check_nc(
      nc_put_vara_text(file_id, times_var, time_start, time_count, timestamp),
      "write Times");

  const auto u_bxs_values = make_raw_boundary_values(1, 2, 3, 1000.0F);
  const auto u_bxe_values = make_raw_boundary_values(1, 2, 3, 2000.0F);
  const auto u_bys_values = make_raw_boundary_values(1, 2, 5, 3000.0F);
  const auto u_bye_values = make_raw_boundary_values(1, 2, 5, 4000.0F);
  check_nc(nc_put_var_float(file_id, u_bxs, u_bxs_values.data()), "write U_BXS");
  check_nc(nc_put_var_float(file_id, u_bxe, u_bxe_values.data()), "write U_BXE");
  check_nc(nc_put_var_float(file_id, u_bys, u_bys_values.data()), "write U_BYS");
  check_nc(nc_put_var_float(file_id, u_bye, u_bye_values.data()), "write U_BYE");
  check_nc(nc_close(file_id), "close synthetic wrfbdy");
}

std::size_t expected_active_boundary_point_count(
    const std::size_t nx,
    const std::size_t ny,
    const std::size_t nz,
    const std::size_t width) {
  const auto interior_nx = nx > 2 * width ? nx - 2 * width : 0;
  const auto interior_ny = ny > 2 * width ? ny - 2 * width : 0;
  return nx * ny * nz - interior_nx * interior_ny * nz;
}

template <typename Predicate>
std::size_t count_active_boundary_points_matching(
    const tywrf::FieldView3D<const float> view,
    const std::size_t width,
    Predicate predicate) {
  const auto active_nx =
      static_cast<std::size_t>(view.nx - view.halo.i_lower - view.halo.i_upper);
  const auto active_ny =
      static_cast<std::size_t>(view.ny - view.halo.j_lower - view.halo.j_upper);
  const auto active_nz =
      static_cast<std::size_t>(view.nz - view.halo.k_lower - view.halo.k_upper);
  const auto upper_i_begin = active_nx > width ? active_nx - width : 0;
  const auto upper_j_begin = active_ny > width ? active_ny - width : 0;

  std::size_t count = 0;
  for (std::size_t j = 0; j < active_ny; ++j) {
    const auto storage_j = view.halo.j_lower + static_cast<std::int32_t>(j);
    for (std::size_t k = 0; k < active_nz; ++k) {
      const auto storage_k = view.halo.k_lower + static_cast<std::int32_t>(k);
      for (std::size_t i = 0; i < active_nx; ++i) {
        const auto on_boundary =
            i < width || i >= upper_i_begin || j < width || j >= upper_j_begin;
        if (!on_boundary) {
          continue;
        }
        const auto storage_i = view.halo.i_lower + static_cast<std::int32_t>(i);
        if (predicate(view(storage_i, storage_j, storage_k))) {
          ++count;
        }
      }
    }
  }
  return count;
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
  expect(report.operation == tywrf::io::ForcingApplyOperation::validation_only,
         "boundary validation operation is explicit");
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
  expect(report.operation == tywrf::io::ForcingApplyOperation::synthetic_nudging_delta,
         "synthetic delta operation is explicit");
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

void test_direct_boundary_copy_skeleton_writes_i_upper_active_edge() {
  tywrf::State<float> state(make_grid());
  zero_state(state);

  const auto packed = make_i_side_t_boundary_pack();
  const auto report = tywrf::io::apply_boundary_copy_skeleton_to_state(
      state,
      "T",
      tywrf::io::BoundarySide::i_upper,
      packed);

  expect(report.status == tywrf::io::ForcingApplyStatus::ok,
         "i_upper direct boundary copy status ok");
  expect(report.field_count == 1, "i_upper copy reports one field");
  expect(report.point_count == 12, "i_upper copy reports copied point count");
  expect(report.operation ==
             tywrf::io::ForcingApplyOperation::direct_boundary_copy_skeleton,
         "i_upper copy operation is direct boundary skeleton");
  expect(report.would_modify_state, "i_upper copy reports state write");
  expect(!report.synthetic, "i_upper copy is real forcing skeleton, not synthetic");

  const auto view = state.view().t;
  expect_close(view(3, 1, 0), 0.0F, "i_upper first copied point");
  expect_close(view(4, 1, 0), 100.0F, "i_upper second copied point");
  expect_close(view(3, 3, 1), 12.0F, "i_upper copied y/k point");
  expect_close(view(4, 3, 1), 112.0F, "i_upper copied last point");
  expect_close(view(1, 1, 0), 0.0F, "i_upper copy leaves lower active edge unchanged");
  expect_close(view(0, 1, 0), 0.0F, "i_upper copy leaves i_lower halo unchanged");
}

void test_direct_boundary_copy_skeleton_writes_j_lower_2d_edge() {
  tywrf::State<float> state(make_grid());
  zero_state(state);

  tywrf::io::PackedForcingField packed;
  packed.layout = tywrf::io::make_canonical_forcing_layout(4, 2, 1);
  packed.values = {
      1.0F, 2.0F, 3.0F, 4.0F,
      11.0F, 12.0F, 13.0F, 14.0F,
  };

  const auto report = tywrf::io::apply_boundary_copy_skeleton_to_state(
      state,
      "MU",
      tywrf::io::BoundarySide::j_lower,
      packed);

  expect(report.status == tywrf::io::ForcingApplyStatus::ok,
         "j_lower 2D direct boundary copy status ok");
  expect(report.operation ==
             tywrf::io::ForcingApplyOperation::direct_boundary_copy_skeleton,
         "j_lower 2D copy operation is direct boundary skeleton");
  expect(report.point_count == packed.values.size(),
         "j_lower 2D copy reports copied point count");

  const auto view = state.view().mu;
  expect_close(view(1, 1), 1.0F, "j_lower first copied row");
  expect_close(view(4, 1), 4.0F, "j_lower first row last copied column");
  expect_close(view(1, 2), 11.0F, "j_lower second copied row");
  expect_close(view(4, 2), 14.0F, "j_lower second row last copied column");
  expect_close(view(1, 3), 0.0F, "j_lower copy leaves next active row unchanged");
  expect_close(view(1, 0), 0.0F, "j_lower copy leaves j_lower halo unchanged");
}

void test_krosa_boundary_copy_skeleton_applies_all_u_edges() {
  const auto path =
      std::filesystem::temp_directory_path() / "tywrf_forcing_apply_u_wrfbdy.nc";
  std::filesystem::remove(path);
  create_synthetic_u_wrfbdy_for_apply(path);

  const tywrf::io::KrosaForcingReader boundary(
      path,
      tywrf::io::KrosaForcingKind::boundary);
  tywrf::State<float> state(make_grid());
  zero_state(state);

  const auto report = tywrf::io::apply_krosa_boundary_copy_skeleton_to_state(
      state,
      boundary,
      "U",
      0);

  expect(report.status == tywrf::io::ForcingApplyStatus::ok,
         "KROSA U batch boundary copy status ok");
  expect(report.field_count == 1, "KROSA U batch reports one field");
  expect(report.side_count == 4, "KROSA U batch reports four sides");
  expect(report.point_count == 32,
         "KROSA U batch accumulates all side point counts");
  expect(report.operation ==
             tywrf::io::ForcingApplyOperation::direct_boundary_copy_skeleton,
         "KROSA U batch operation is direct boundary copy skeleton");
  expect(report.would_modify_state, "KROSA U batch reports state write");
  expect(!report.synthetic, "KROSA U batch is not synthetic nudging");
  expect(report.detail.find("direct_boundary_copy_skeleton") != std::string::npos,
         "KROSA U batch detail names direct copy skeleton operation");
  expect(report.detail.find("not WRF lateral relaxation") != std::string::npos,
         "KROSA U batch detail rejects WRF relaxation semantics");

  const auto view = state.view().u;
  expect_close(view(1, 2, 0), 1001.0F, "KROSA U_BXS copied lower i edge");
  expect_close(view(5, 2, 1), 2011.0F, "KROSA U_BXE copied upper i edge");
  expect_close(view(3, 1, 0), 3002.0F, "KROSA U_BYS copied lower j edge");
  expect_close(view(4, 3, 1), 4013.0F, "KROSA U_BYE copied upper j edge");
  expect_close(view(3, 2, 0), 0.0F, "KROSA U batch leaves interior U point");
  expect_close(view(0, 2, 0), 0.0F, "KROSA U batch leaves i halo point");

  const auto changed = count_active_boundary_points_matching(
      tywrf::FieldView3D<const float>{
          view.data,
          view.nx,
          view.ny,
          view.nz,
          view.stride_i,
          view.stride_k,
          view.stride_j,
          view.halo},
      1,
      [](const float value) { return value != 0.0F; });
  expect(changed == expected_active_boundary_point_count(5, 3, 2, 1),
         "KROSA U batch changes the expected active edge point count");

  std::filesystem::remove(path);
}

int run_krosa_reference_boundary_pack_smoke() {
  const auto root = reference_dir();
  const auto wrfbdy = root / "wrfbdy_d01";
  const auto wrfinput_d01 = root / "wrfinput_d01";
  const auto wrfout_d01 = root / "wrfout_d01_2025-07-26_00:00:00";
  if (!std::filesystem::exists(wrfbdy) ||
      (!std::filesystem::exists(wrfinput_d01) &&
       !std::filesystem::exists(wrfout_d01))) {
    std::cout << "Skipping forcing apply KROSA reference smoke; wrfbdy_d01 and d01 "
                 "state/grid source are not present under "
              << root << '\n';
    return 0;
  }

  const tywrf::io::KrosaForcingReader boundary(
      wrfbdy,
      tywrf::io::KrosaForcingKind::boundary);
  const auto u_bxs_metadata = boundary.variable_metadata("U_BXS");
  const auto u_bxe_metadata = boundary.variable_metadata("U_BXE");
  const auto u_bys_metadata = boundary.variable_metadata("U_BYS");
  const auto u_bye_metadata = boundary.variable_metadata("U_BYE");
  const auto u_bxs_slice = boundary.read_float_time_slice("U_BXS", 0);
  const auto packed =
      tywrf::io::pack_boundary_x_side_raw_to_canonical(
          u_bxs_slice.values,
          u_bxs_slice.metadata.slice_shape);

  const auto grid_source =
      std::filesystem::exists(wrfinput_d01) ? wrfinput_d01 : wrfout_d01;
  tywrf::State<float> state(tywrf::io::derive_grid_from_wrf_file(grid_source));
  fill_storage(state.u, std::numeric_limits<float>::quiet_NaN());
  const auto report = tywrf::io::validate_boundary_pack_for_state(
      state,
      "U",
      tywrf::io::BoundarySide::i_lower,
      packed);

  expect(report.status == tywrf::io::ForcingApplyStatus::ok,
         "KROSA U_BXS validates against d01 U state shape");
  expect(report.field_count == 1, "KROSA U_BXS reports one validated field");
  expect(report.point_count == u_bxs_slice.values.size(),
         "KROSA U_BXS report preserves packed point count");
  expect(!report.would_modify_state, "KROSA U_BXS validation is read-only");
  expect(!report.synthetic, "KROSA U_BXS validation is real forcing, not synthetic");

  const auto copy_report = tywrf::io::apply_krosa_boundary_copy_skeleton_to_state(
      state,
      boundary,
      "U",
      0);
  expect(copy_report.status == tywrf::io::ForcingApplyStatus::ok,
         "KROSA U four-side direct copy skeleton status ok");
  expect(copy_report.field_count == 1,
         "KROSA U four-side direct copy reports one field");
  expect(copy_report.side_count == 4,
         "KROSA U four-side direct copy reports four sides");
  expect(copy_report.operation ==
             tywrf::io::ForcingApplyOperation::direct_boundary_copy_skeleton,
         "KROSA U four-side direct copy skeleton operation");
  expect(copy_report.point_count ==
             u_bxs_metadata.values_per_time_slice +
                 u_bxe_metadata.values_per_time_slice +
                 u_bys_metadata.values_per_time_slice +
                 u_bye_metadata.values_per_time_slice,
         "KROSA U four-side direct copy accumulates point count");
  expect(copy_report.would_modify_state,
         "KROSA U four-side direct copy reports state write");
  expect(!copy_report.synthetic,
         "KROSA U four-side direct copy is not synthetic nudging");
  const auto u_view = state.view().u;
  expect_close(
      u_view(u_view.halo.i_lower, u_view.halo.j_lower, u_view.halo.k_lower),
      packed.at(0, 0, 0),
      "KROSA U four-side direct copy writes first active U boundary value");

  const auto config = state.grid.config();
  expect(config.mass_nx == 265, "KROSA d01 grid mass_nx");
  expect(config.mass_ny == 429, "KROSA d01 grid mass_ny");
  expect(config.mass_nz == 59, "KROSA d01 grid mass_nz");
  expect(packed.layout.nx == 5, "KROSA U_BXS boundary width");
  expect(packed.layout.ny == static_cast<std::size_t>(config.mass_ny),
         "KROSA U_BXS y extent matches d01 U state");
  expect(packed.layout.nz == static_cast<std::size_t>(config.mass_nz),
         "KROSA U_BXS z extent matches d01 U state");
  const auto active_u_nx =
      static_cast<std::size_t>(u_view.nx - u_view.halo.i_lower - u_view.halo.i_upper);
  const auto active_u_ny =
      static_cast<std::size_t>(u_view.ny - u_view.halo.j_lower - u_view.halo.j_upper);
  const auto active_u_nz =
      static_cast<std::size_t>(u_view.nz - u_view.halo.k_lower - u_view.halo.k_upper);
  const auto changed = count_active_boundary_points_matching(
      tywrf::FieldView3D<const float>{
          u_view.data,
          u_view.nx,
          u_view.ny,
          u_view.nz,
          u_view.stride_i,
          u_view.stride_k,
          u_view.stride_j,
          u_view.halo},
      packed.layout.nx,
      [](const float value) { return !std::isnan(value); });
  expect(changed == expected_active_boundary_point_count(
                        active_u_nx,
                        active_u_ny,
                        active_u_nz,
                        packed.layout.nx),
         "KROSA U four-side direct copy changes expected active edge points");
  return 0;
}

}  // namespace

int main() {
  test_boundary_validation_reports_shape_and_does_not_write_state();
  test_boundary_validation_rejects_wrong_shape();
  test_synthetic_nudging_delta_writes_only_requested_window();
  test_synthetic_nudging_delta_rejects_out_of_range_window();
  test_direct_boundary_copy_skeleton_writes_i_upper_active_edge();
  test_direct_boundary_copy_skeleton_writes_j_lower_2d_edge();
  test_krosa_boundary_copy_skeleton_applies_all_u_edges();
  if (const int status = run_krosa_reference_boundary_pack_smoke(); status != 0) {
    return status;
  }

  if (failures != 0) {
    return 1;
  }

  std::cout << "Validated forcing apply skeleton contracts\n";
  return 0;
}
