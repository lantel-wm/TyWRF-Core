#include "tywrf/dynamics/pressure_refresh.hpp"

#include "tywrf/grid.hpp"
#include "tywrf/state.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

int failures = 0;

constexpr float kSentinel = -9'999.0F;
constexpr float kOldOverlapP = 12'345.0F;
constexpr float kWouldBeParentP = 54'321.0F;
constexpr float kPTop = 5'000.0F;

void expect(const bool condition, const std::string_view message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

void expect_near(
    const double actual,
    const double expected,
    const double tolerance,
    const std::string_view message) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL: " << message << " actual=" << actual
              << " expected=" << expected << '\n';
    ++failures;
  }
}

template <typename Storage>
void fill_storage(Storage& storage, const float value) {
  std::fill(storage.data(), storage.data() + storage.size(), value);
}

tywrf::Grid make_grid() {
  return tywrf::Grid({
      .mass_nx = 4,
      .mass_ny = 3,
      .mass_nz = 3,
      .full_nz = 4,
      .halo = tywrf::Halo3D{1, 1, 1, 1, 1, 1},
  });
}

tywrf::nest::RemapWindow make_mass_window() {
  return {
      tywrf::nest::HorizontalStagger::mass,
      1,
      1,
      0,
      0,
      3,
      2,
      1,
      1,
  };
}

tywrf::nest::RemapPlan make_plan() {
  tywrf::nest::RemapPlan plan{};
  plan.mass = make_mass_window();
  return plan;
}

[[nodiscard]] bool in_window(
    const tywrf::nest::RemapWindow& window,
    const std::int32_t i,
    const std::int32_t j) noexcept {
  return i >= window.new_i_begin && i < window.new_i_begin + window.extent_i &&
         j >= window.new_j_begin && j < window.new_j_begin + window.extent_j;
}

struct Coefficients {
  float c3f[4] = {1.0F, 0.7F, 0.4F, 0.1F};
  float c4f[4] = {0.0F, 0.0F, 0.0F, 0.0F};
  float c3h[3] = {0.85F, 0.55F, 0.25F};
  float c4h[3] = {0.0F, 0.0F, 0.0F};
};

struct RefreshFixture {
  explicit RefreshFixture(const tywrf::Grid& grid)
      : p(grid.mass_layout()),
        pb(grid.mass_layout()),
        t(grid.mass_layout()),
        alb(grid.mass_layout()),
        ph(grid.w_layout()),
        phb(grid.w_layout()),
        mu(grid.surface_layout()),
        mub(grid.surface_layout()) {}

  tywrf::FieldStorage3D<float> p;
  tywrf::FieldStorage3D<float> pb;
  tywrf::FieldStorage3D<float> t;
  tywrf::FieldStorage3D<float> alb;
  tywrf::FieldStorage3D<float> ph;
  tywrf::FieldStorage3D<float> phb;
  tywrf::FieldStorage2D<float> mu;
  tywrf::FieldStorage2D<float> mub;
  Coefficients coefficients{};
};

[[nodiscard]] tywrf::dynamics::PressureRefreshInputs make_inputs(
    RefreshFixture& fixture) noexcept {
  return {
      fixture.p.view(),
      static_cast<const tywrf::FieldStorage3D<float>&>(fixture.pb).view(),
      static_cast<const tywrf::FieldStorage3D<float>&>(fixture.t).view(),
      static_cast<const tywrf::FieldStorage3D<float>&>(fixture.alb).view(),
      static_cast<const tywrf::FieldStorage3D<float>&>(fixture.ph).view(),
      static_cast<const tywrf::FieldStorage3D<float>&>(fixture.phb).view(),
      static_cast<const tywrf::FieldStorage2D<float>&>(fixture.mu).view(),
      static_cast<const tywrf::FieldStorage2D<float>&>(fixture.mub).view(),
      {fixture.coefficients.c3f, 4},
      {fixture.coefficients.c4f, 4},
      {fixture.coefficients.c3h, 3},
      {fixture.coefficients.c4h, 3},
      kPTop,
  };
}

void fill_mass_inputs(
    RefreshFixture& fixture,
    const tywrf::nest::RemapWindow& window) {
  fill_storage(fixture.p, kSentinel);
  fill_storage(fixture.pb, kSentinel);
  fill_storage(fixture.t, kSentinel);
  fill_storage(fixture.alb, kSentinel);

  const auto layout = fixture.p.layout();
  auto p_view = fixture.p.view();
  auto pb_view = fixture.pb.view();
  auto t_view = fixture.t.view();
  auto alb_view = fixture.alb.view();

  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t k = 0; k < layout.active_nz(); ++k) {
      for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
        const auto ii = layout.i_begin() + i;
        const auto jj = layout.j_begin() + j;
        const auto kk = layout.k_begin() + k;
        p_view(ii, jj, kk) =
            in_window(window, i, j) ? kOldOverlapP : kWouldBeParentP;
        pb_view(ii, jj, kk) =
            80'000.0F - 25'000.0F * static_cast<float>(k) -
            200.0F * static_cast<float>(i) - 100.0F * static_cast<float>(j);
        t_view(ii, jj, kk) =
            2.0F + 0.25F * static_cast<float>(i) +
            0.5F * static_cast<float>(j) + 0.125F * static_cast<float>(k);
        alb_view(ii, jj, kk) =
            0.8F + 0.05F * static_cast<float>(k) +
            0.01F * static_cast<float>(i);
      }
    }
  }
}

void fill_full_level_geopotential(RefreshFixture& fixture) {
  fill_storage(fixture.ph, kSentinel);
  fill_storage(fixture.phb, kSentinel);

  const auto layout = fixture.ph.layout();
  auto ph_view = fixture.ph.view();
  auto phb_view = fixture.phb.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t k = 0; k < layout.active_nz(); ++k) {
      for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
        const auto ii = layout.i_begin() + i;
        const auto jj = layout.j_begin() + j;
        const auto kk = layout.k_begin() + k;
        ph_view(ii, jj, kk) =
            2'500.0F * static_cast<float>(k) + 5.0F * static_cast<float>(i);
        phb_view(ii, jj, kk) =
            22'500.0F * static_cast<float>(k) +
            20.0F * static_cast<float>(i) +
            30.0F * static_cast<float>(j);
      }
    }
  }
}

void fill_surface_mass(RefreshFixture& fixture) {
  fill_storage(fixture.mu, kSentinel);
  fill_storage(fixture.mub, kSentinel);

  const auto layout = fixture.mu.layout();
  auto mu_view = fixture.mu.view();
  auto mub_view = fixture.mub.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
      const auto ii = layout.i_begin() + i;
      const auto jj = layout.j_begin() + j;
      mu_view(ii, jj) =
          1'000.0F + 10.0F * static_cast<float>(i) +
          20.0F * static_cast<float>(j);
      mub_view(ii, jj) =
          79'000.0F - 30.0F * static_cast<float>(i) -
          10.0F * static_cast<float>(j);
    }
  }
}

void fill_fixture(RefreshFixture& fixture, const tywrf::nest::RemapWindow& window) {
  fill_mass_inputs(fixture, window);
  fill_full_level_geopotential(fixture);
  fill_surface_mass(fixture);
}

void expect_halos_unchanged(const tywrf::FieldStorage3D<float>& p) {
  const auto layout = p.layout();
  const auto view = p.view();
  for (std::int32_t j = 0; j < layout.ny; ++j) {
    for (std::int32_t k = 0; k < layout.nz; ++k) {
      for (std::int32_t i = 0; i < layout.nx; ++i) {
        const bool halo =
            i < layout.i_begin() || i >= layout.i_end() ||
            j < layout.j_begin() || j >= layout.j_end() ||
            k < layout.k_begin() || k >= layout.k_end();
        if (halo) {
          expect(view(i, j, k) == kSentinel, "P halo remains unchanged");
        }
      }
    }
  }
}

void expect_storage_unchanged(
    const tywrf::FieldStorage3D<float>& field,
    const std::vector<float>& before,
    const std::string_view label) {
  expect(field.size() == before.size(), "field copy size matches");
  for (std::size_t idx = 0; idx < field.size(); ++idx) {
    expect(field.data()[idx] == before[idx], label);
  }
}

void test_exposed_cells_refresh_without_parent_pressure_copy() {
  const auto grid = make_grid();
  const auto window = make_mass_window();
  RefreshFixture fixture(grid);
  fill_fixture(fixture, window);
  const std::vector<float> pb_before(
      fixture.pb.data(), fixture.pb.data() + fixture.pb.size());
  const std::vector<float> alb_before(
      fixture.alb.data(), fixture.alb.data() + fixture.alb.size());

  const auto report = tywrf::dynamics::refresh_krosa_moving_nest_pressure(
      make_plan(), make_inputs(fixture));

  expect(report.ok(), "pressure refresh succeeds");
  expect(report.target_column_count == 6, "six exposed columns targeted");
  expect(report.refreshed_column_count == 6, "six exposed columns refreshed");
  expect(report.refreshed_point_count == 18, "all exposed mass points refreshed");
  expect(report.skipped_point_count == 0, "no pressure points skipped");
  expect(report.invalid_point_count == 0, "no invalid pressure points");
  expect(!report.full_column_mode, "default mode is exposed cells only");
  expect(!report.touched_overlap_cells, "default mode does not touch overlap");
  expect(!report.touched_halo_cells, "pressure refresh does not touch halos");
  expect(!report.modified_pb, "pressure refresh reports PB unchanged");
  expect(!report.copied_parent_pressure, "pressure refresh does not copy parent P");
  expect(report.used_krosa_hypsometric_opt2, "KROSA hypsometric branch reported");
  expect(report.used_use_theta_m, "USE_THETA_M=1 branch reported");
  expect(!report.used_qvapor, "QVAPOR is not used");
  expect(report.used_wrf_vertical_coefficients, "WRF vertical coefficients used");

  const auto layout = fixture.p.layout();
  const auto p_view = fixture.p.view();
  for (std::int32_t j = 0; j < layout.active_ny(); ++j) {
    for (std::int32_t k = 0; k < layout.active_nz(); ++k) {
      for (std::int32_t i = 0; i < layout.active_nx(); ++i) {
        const auto value = p_view(
            layout.i_begin() + i,
            layout.j_begin() + j,
            layout.k_begin() + k);
        if (in_window(window, i, j)) {
          expect(value == kOldOverlapP, "overlap P remains old value");
        } else {
          expect(std::isfinite(value), "exposed P is finite");
          expect(value != kWouldBeParentP, "exposed P is not copied parent P");
          expect(value != kOldOverlapP, "exposed P is recomputed");
        }
      }
    }
  }

  expect_halos_unchanged(fixture.p);
  expect_storage_unchanged(fixture.pb, pb_before, "PB remains unchanged");
  expect_storage_unchanged(fixture.alb, alb_before, "ALB remains unchanged");
}

void test_full_column_mode_is_explicit() {
  const auto grid = make_grid();
  const auto window = make_mass_window();
  RefreshFixture fixture(grid);
  fill_fixture(fixture, window);

  tywrf::dynamics::KrosaPressureRefreshOptions options{};
  options.region = tywrf::dynamics::PressureRefreshRegion::full_active_columns;
  const auto report = tywrf::dynamics::refresh_krosa_moving_nest_pressure(
      window, make_inputs(fixture), options);

  expect(report.ok(), "full-column pressure refresh succeeds");
  expect(report.full_column_mode, "full-column mode is reported");
  expect(report.touched_overlap_cells, "full-column mode may touch overlap");
  expect(report.target_column_count == 12, "full active domain targeted");
  expect(report.refreshed_column_count == 12, "full active domain refreshed");
  expect(report.refreshed_point_count == 36, "all active mass points refreshed");
}

void test_invalid_contract_and_unsupported_options_do_not_write() {
  const auto grid = make_grid();
  const auto window = make_mass_window();
  RefreshFixture fixture(grid);
  fill_fixture(fixture, window);

  const auto p_before = fixture.p.data()[0];
  auto null_inputs = make_inputs(fixture);
  null_inputs.p.data = nullptr;
  const auto null_report =
      tywrf::dynamics::refresh_krosa_moving_nest_pressure(window, null_inputs);
  expect(
      null_report.result.status == tywrf::nest::NestStatus::invalid_contract,
      "null P view returns invalid_contract");
  expect(fixture.p.data()[0] == p_before, "null P view does not write");

  tywrf::FieldStorage3D<float> bad_ph(grid.mass_layout());
  fill_storage(bad_ph, 0.0F);
  auto bad_shape_inputs = make_inputs(fixture);
  bad_shape_inputs.ph =
      static_cast<const tywrf::FieldStorage3D<float>&>(bad_ph).view();
  const auto bad_shape_report =
      tywrf::dynamics::refresh_krosa_moving_nest_pressure(
          window, bad_shape_inputs);
  expect(
      bad_shape_report.result.status == tywrf::nest::NestStatus::invalid_contract,
      "PH mass-level shape returns invalid_contract");

  auto missing_coeff_inputs = make_inputs(fixture);
  missing_coeff_inputs.c3f = {};
  const auto missing_coeff_report =
      tywrf::dynamics::refresh_krosa_moving_nest_pressure(
          window, missing_coeff_inputs);
  expect(
      missing_coeff_report.result.status ==
          tywrf::nest::NestStatus::invalid_contract,
      "missing coefficients return invalid_contract");

  tywrf::dynamics::KrosaPressureRefreshOptions options{};
  options.formula = tywrf::dynamics::PressureRefreshFormula::wrf_start_em_general;
  const auto unsupported_formula_report =
      tywrf::dynamics::refresh_krosa_moving_nest_pressure(
          window, make_inputs(fixture), options);
  expect(
      unsupported_formula_report.result.status ==
          tywrf::nest::NestStatus::not_implemented,
      "general WRF formula is explicit not_implemented");

  options = {};
  options.thermodynamic_mode =
      tywrf::dynamics::PressureRefreshThermodynamicMode::use_theta_m_0_qvapor;
  const auto unsupported_moisture_report =
      tywrf::dynamics::refresh_krosa_moving_nest_pressure(
          window, make_inputs(fixture), options);
  expect(
      unsupported_moisture_report.result.status ==
          tywrf::nest::NestStatus::not_implemented,
      "USE_THETA_M=0 QVAPOR path is explicit not_implemented");
}

void test_invalid_cell_data_is_counted_and_skipped() {
  const auto grid = make_grid();
  const auto window = make_mass_window();
  RefreshFixture fixture(grid);
  fill_fixture(fixture, window);

  const auto layout = fixture.mu.layout();
  auto mu_view = fixture.mu.view();
  auto mub_view = fixture.mub.view();
  const std::int32_t exposed_i = 3;
  const std::int32_t exposed_j = 2;
  mu_view(layout.i_begin() + exposed_i, layout.j_begin() + exposed_j) =
      -mub_view(layout.i_begin() + exposed_i, layout.j_begin() + exposed_j);

  const auto report = tywrf::dynamics::refresh_krosa_moving_nest_pressure(
      window, make_inputs(fixture));

  expect(report.ok(), "invalid cell data leaves contract ok");
  expect(report.target_column_count == 6, "invalid data still targets exposed columns");
  expect(report.refreshed_column_count == 5, "one invalid exposed column is skipped");
  expect(report.refreshed_point_count == 15, "valid exposed points refresh");
  expect(report.invalid_point_count == 3, "invalid pressure points counted");
  expect(report.skipped_point_count == 3, "invalid pressure points skipped");

  const auto p_layout = fixture.p.layout();
  const auto p_view = fixture.p.view();
  for (std::int32_t k = 0; k < p_layout.active_nz(); ++k) {
    expect(
        p_view(
            p_layout.i_begin() + exposed_i,
            p_layout.j_begin() + exposed_j,
            p_layout.k_begin() + k) == kWouldBeParentP,
        "invalid exposed P remains unchanged");
  }
}

void test_default_path_does_not_emit_formula_observation() {
  const auto grid = make_grid();
  const auto window = make_mass_window();
  RefreshFixture fixture(grid);
  fill_fixture(fixture, window);

  tywrf::dynamics::PressureRefreshFormulaObservation records[1]{};
  records[0].status =
      tywrf::dynamics::PressureRefreshObservationStatus::request_out_of_bounds;
  records[0].i = 99;

  const auto report = tywrf::dynamics::refresh_krosa_moving_nest_pressure(
      window, make_inputs(fixture));

  expect(report.ok(), "default pressure refresh still succeeds");
  expect(report.observation_request_count == 0, "default has no observation requests");
  expect(report.observation_record_count == 0, "default emits no observations");
  expect(report.observation_valid_count == 0, "default has no valid observations");
  expect(records[0].status ==
             tywrf::dynamics::PressureRefreshObservationStatus::
                 request_out_of_bounds,
         "unconfigured observation buffer is untouched");
  expect(records[0].i == 99, "unconfigured observation payload is untouched");
}

void test_formula_observation_records_requested_point_terms() {
  const auto grid = make_grid();
  const auto window = make_mass_window();
  RefreshFixture fixture(grid);
  fill_fixture(fixture, window);

  constexpr std::int32_t i = 3;
  constexpr std::int32_t j = 2;
  constexpr std::int32_t k = 1;
  const tywrf::dynamics::PressureRefreshObservationRequest requests[1]{
      {i, j, k},
  };
  tywrf::dynamics::PressureRefreshFormulaObservation records[1]{};

  tywrf::dynamics::KrosaPressureRefreshOptions options{};
  options.observation = {requests, 1, records, 1};
  const auto report = tywrf::dynamics::refresh_krosa_moving_nest_pressure(
      window, make_inputs(fixture), options);

  expect(report.ok(), "observed pressure refresh succeeds");
  expect(report.observation_request_count == 1, "one observation requested");
  expect(report.observation_record_count == 1, "one observation record available");
  expect(report.observation_valid_count == 1, "one valid observation recorded");
  expect(report.observation_invalid_count == 0, "no invalid observations recorded");

  const auto& observation = records[0];
  expect(
      observation.status ==
          tywrf::dynamics::PressureRefreshObservationStatus::recorded,
      "observation status is recorded");
  expect(observation.valid == 1U, "observation valid flag is set");
  expect(
      observation.i == i && observation.j == j && observation.k == k,
      "observation coordinates match request");

  const double mu_total =
      (1'000.0 + 10.0 * static_cast<double>(i) +
       20.0 * static_cast<double>(j)) +
      (79'000.0 - 30.0 * static_cast<double>(i) -
       10.0 * static_cast<double>(j));
  const double pfu =
      static_cast<double>(fixture.coefficients.c3f[k + 1]) * mu_total + kPTop;
  const double pfd =
      static_cast<double>(fixture.coefficients.c3f[k]) * mu_total + kPTop;
  const double phm =
      static_cast<double>(fixture.coefficients.c3h[k]) * mu_total + kPTop;
  const double log_ratio = std::log(pfd / pfu);
  const double phi_lower =
      (2'500.0 * static_cast<double>(k) + 5.0 * static_cast<double>(i)) +
      (22'500.0 * static_cast<double>(k) + 20.0 * static_cast<double>(i) +
       30.0 * static_cast<double>(j));
  const double phi_upper =
      (2'500.0 * static_cast<double>(k + 1) + 5.0 * static_cast<double>(i)) +
      (22'500.0 * static_cast<double>(k + 1) + 20.0 * static_cast<double>(i) +
       30.0 * static_cast<double>(j));
  const double delta_phi = phi_upper - phi_lower;
  const double alb = 0.8 + 0.05 * static_cast<double>(k) +
                     0.01 * static_cast<double>(i);
  const double pb = 80'000.0 - 25'000.0 * static_cast<double>(k) -
                    200.0 * static_cast<double>(i) -
                    100.0 * static_cast<double>(j);
  const double theta =
      300.0 + 2.0 + 0.25 * static_cast<double>(i) +
      0.5 * static_cast<double>(j) + 0.125 * static_cast<double>(k);
  const double alpha_total = delta_phi / (phm * log_ratio);
  const double alpha_perturbation = alpha_total - alb;
  const double alpha_from_wrf_branch = alpha_perturbation + alb;
  const double pressure_base = 287.0 * theta / (100'000.0 * alpha_from_wrf_branch);
  const double total_pressure =
      100'000.0 * std::pow(pressure_base, 1004.5 / (1004.5 - 287.0));
  const double perturbation_pressure = total_pressure - pb;

  expect_near(observation.mu_total, mu_total, 1.0e-9, "mu_total matches fixture");
  expect_near(observation.pfu, pfu, 1.0e-6, "pfu matches fixture");
  expect_near(observation.pfd, pfd, 1.0e-6, "pfd matches fixture");
  expect_near(observation.phm, phm, 1.0e-6, "phm matches fixture");
  expect_near(observation.log_ratio, log_ratio, 1.0e-12, "log_ratio matches");
  expect_near(observation.phi_lower, phi_lower, 1.0e-9, "phi_lower matches");
  expect_near(observation.phi_upper, phi_upper, 1.0e-9, "phi_upper matches");
  expect_near(observation.delta_phi, delta_phi, 1.0e-9, "delta_phi matches");
  expect_near(observation.alb, alb, 1.0e-7, "ALB matches");
  expect_near(observation.pb, pb, 1.0e-6, "PB matches");
  expect_near(observation.theta, theta, 1.0e-7, "theta matches");
  expect_near(observation.alpha_total, alpha_total, 1.0e-12, "alpha_total matches");
  expect_near(
      observation.alpha_perturbation,
      alpha_perturbation,
      1.0e-6,
      "alpha_perturbation matches");
  expect_near(
      observation.alpha_from_wrf_branch,
      alpha_from_wrf_branch,
      1.0e-12,
      "alpha_from_wrf_branch matches");
  expect_near(observation.pressure_base, pressure_base, 1.0e-12, "pressure_base matches");
  expect_near(observation.total_pressure, total_pressure, 1.0e-2, "total_pressure matches");
  expect_near(
      observation.perturbation_pressure_pa,
      perturbation_pressure,
      1.0e-2,
      "perturbation pressure matches");

  const auto p_layout = fixture.p.layout();
  const auto p_view = fixture.p.view();
  expect_near(
      p_view(p_layout.i_begin() + i, p_layout.j_begin() + j, p_layout.k_begin() + k),
      static_cast<float>(perturbation_pressure),
      1.0e-3,
      "P field receives observed perturbation pressure");
}

void test_formula_observation_reports_out_of_target_and_invalid_points() {
  const auto grid = make_grid();
  const auto window = make_mass_window();
  RefreshFixture fixture(grid);
  fill_fixture(fixture, window);

  const auto surface_layout = fixture.mu.layout();
  auto mu_view = fixture.mu.view();
  auto mub_view = fixture.mub.view();
  constexpr std::int32_t invalid_i = 3;
  constexpr std::int32_t invalid_j = 2;
  mu_view(surface_layout.i_begin() + invalid_i, surface_layout.j_begin() + invalid_j) =
      -mub_view(surface_layout.i_begin() + invalid_i, surface_layout.j_begin() + invalid_j);

  const tywrf::dynamics::PressureRefreshObservationRequest requests[3]{
      {1, 1, 0},
      {4, 0, 0},
      {invalid_i, invalid_j, 0},
  };
  tywrf::dynamics::PressureRefreshFormulaObservation records[3]{};

  tywrf::dynamics::KrosaPressureRefreshOptions options{};
  options.observation = {requests, 3, records, 3};
  const auto report = tywrf::dynamics::refresh_krosa_moving_nest_pressure(
      window, make_inputs(fixture), options);

  expect(report.ok(), "mixed observation pressure refresh succeeds");
  expect(report.observation_request_count == 3, "three observations requested");
  expect(report.observation_record_count == 3, "three observation records written");
  expect(report.observation_outside_target_region_count == 1, "one overlap request skipped");
  expect(report.observation_out_of_bounds_count == 1, "one out-of-bounds request reported");
  expect(report.observation_invalid_count == 1, "one invalid formula observation reported");
  expect(report.observation_valid_count == 0, "no valid observations in mixed skip test");

  expect(
      records[0].status ==
          tywrf::dynamics::PressureRefreshObservationStatus::
              request_outside_target_region,
      "overlap request reports outside target region");
  expect(
      records[1].status ==
          tywrf::dynamics::PressureRefreshObservationStatus::request_out_of_bounds,
      "out-of-bounds request reports out_of_bounds");
  expect(
      records[2].status ==
          tywrf::dynamics::PressureRefreshObservationStatus::invalid_mu_total,
      "invalid exposed point reports mu_total failure");
  expect(records[2].valid == 0U, "invalid formula observation is not valid");
  expect_near(records[2].mu_total, 0.0, 1.0e-9, "invalid observation captures mu_total");
}

}  // namespace

int main() {
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::VerticalCoefficientView>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::VerticalCoefficientView>);
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::PressureRefreshObservationRequest>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureRefreshObservationRequest>);
  static_assert(std::is_standard_layout_v<
                tywrf::dynamics::PressureRefreshFormulaObservation>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureRefreshFormulaObservation>);
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::PressureRefreshObservationView>);
  static_assert(std::is_trivially_copyable_v<
                tywrf::dynamics::PressureRefreshObservationView>);
  static_assert(
      std::is_standard_layout_v<tywrf::dynamics::KrosaPressureRefreshOptions>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::KrosaPressureRefreshOptions>);
  static_assert(std::is_standard_layout_v<tywrf::dynamics::PressureRefreshInputs>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::PressureRefreshInputs>);
  static_assert(std::is_standard_layout_v<tywrf::dynamics::PressureRefreshReport>);
  static_assert(
      std::is_trivially_copyable_v<tywrf::dynamics::PressureRefreshReport>);

  test_exposed_cells_refresh_without_parent_pressure_copy();
  test_full_column_mode_is_explicit();
  test_invalid_contract_and_unsupported_options_do_not_write();
  test_invalid_cell_data_is_counted_and_skipped();
  test_default_path_does_not_emit_formula_observation();
  test_formula_observation_records_requested_point_terms();
  test_formula_observation_reports_out_of_target_and_invalid_points();

  if (failures != 0) {
    return 1;
  }
  std::cout << "Validated KROSA HYPSOMETRIC_OPT=2 pressure refresh staging\n";
  return 0;
}
