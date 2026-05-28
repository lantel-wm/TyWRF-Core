#include "tywrf/nest/state_exchange.hpp"

#include <algorithm>
#include <cstdint>

namespace tywrf::nest {
namespace {

[[nodiscard]] constexpr NestResult ok_result() noexcept {
  return {NestStatus::ok, "ok"};
}

[[nodiscard]] constexpr NestResult result(
    const NestStatus status,
    const char* message) noexcept {
  return {status, message};
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(const FieldView2D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(const FieldView2D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nx(const FieldView3D<Real>& field) noexcept {
  return field.nx - field.halo.i_lower - field.halo.i_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_ny(const FieldView3D<Real>& field) noexcept {
  return field.ny - field.halo.j_lower - field.halo.j_upper;
}

template <typename Real>
[[nodiscard]] constexpr std::int32_t active_nz(const FieldView3D<Real>& field) noexcept {
  return field.nz - field.halo.k_lower - field.halo.k_upper;
}

[[nodiscard]] constexpr bool window_fits_active_shape(
    const RemapWindow& window,
    const std::int32_t nx,
    const std::int32_t ny) noexcept {
  return window.new_i_begin >= 0 && window.new_j_begin >= 0 &&
         window.extent_i > 0 && window.extent_j > 0 &&
         window.new_i_begin + window.extent_i <= nx &&
         window.new_j_begin + window.extent_j <= ny;
}

[[nodiscard]] constexpr bool valid_active_shape(
    const std::int32_t nx,
    const std::int32_t ny,
    const std::int32_t active_k_count) noexcept {
  return nx > 0 && ny > 0 && active_k_count > 0;
}

void append_region(
    FieldStateExchangePlan& plan,
    const std::int32_t i_begin,
    const std::int32_t j_begin,
    const std::int32_t extent_i,
    const std::int32_t extent_j) noexcept {
  if (extent_i <= 0 || extent_j <= 0 ||
      plan.exposed_region_count >= plan.exposed_regions.size()) {
    return;
  }

  auto& region = plan.exposed_regions[plan.exposed_region_count];
  region.stagger = plan.stagger;
  region.child_i_begin = i_begin;
  region.child_j_begin = j_begin;
  region.extent_i = extent_i;
  region.extent_j = extent_j;
  region.active_k_count = plan.active_k_count;
  region.three_dimensional = plan.three_dimensional;
  region.owns_overlap = false;
  region.owns_halo = false;

  ++plan.exposed_region_count;
  plan.exposed_horizontal_cell_count += region.horizontal_cell_count();
  plan.exchange_point_count += region.point_count();
}

[[nodiscard]] FieldStateExchangePlan make_field_plan(
    const StateExchangeField field,
    const HorizontalStagger stagger,
    const RemapWindow& overlap,
    const std::int32_t active_nx_value,
    const std::int32_t active_ny_value,
    const std::int32_t active_k_count,
    const bool three_dimensional) noexcept {
  FieldStateExchangePlan plan{};
  plan.field = field;
  plan.stagger = stagger;
  plan.active_nx = active_nx_value;
  plan.active_ny = active_ny_value;
  plan.active_k_count = active_k_count;
  plan.three_dimensional = three_dimensional;

  const auto overlap_i0 = overlap.new_i_begin;
  const auto overlap_j0 = overlap.new_j_begin;
  const auto overlap_i1 = overlap.new_i_begin + overlap.extent_i;
  const auto overlap_j1 = overlap.new_j_begin + overlap.extent_j;

  append_region(plan, 0, 0, active_nx_value, overlap_j0);
  append_region(plan, 0, overlap_j1, active_nx_value, active_ny_value - overlap_j1);
  append_region(plan, 0, overlap_j0, overlap_i0, overlap.extent_j);
  append_region(plan, overlap_i1, overlap_j0, active_nx_value - overlap_i1, overlap.extent_j);
  return plan;
}

void add_report_field(
    const FieldStateExchangePlan& field,
    StateExchangeReport& report) noexcept {
  ++report.planned_field_count;
  report.exposed_region_count += field.exposed_region_count;
  report.exposed_horizontal_cell_count += field.exposed_horizontal_cell_count;
  report.exchange_point_count += field.exchange_point_count;
  report.modifies_overlap = report.modifies_overlap || field.owns_overlap;
  report.modifies_halo = report.modifies_halo || field.owns_halo;
  if (field.exchange_point_count > 0) {
    ++report.active_field_count;
    report.requires_parent_interpolation = true;
  }
}

template <typename Real>
[[nodiscard]] NestResult validate_3d_window(
    const RemapWindow& window,
    const FieldView3D<Real>& field) noexcept {
  const auto nx = active_nx(field);
  const auto ny = active_ny(field);
  const auto nz = active_nz(field);
  if (!valid_active_shape(nx, ny, nz)) {
    return result(NestStatus::invalid_contract, "state exchange field shape is empty");
  }
  if (!window_fits_active_shape(window, nx, ny)) {
    return result(
        NestStatus::invalid_contract,
        "state exchange remap window is outside child active extents");
  }
  return ok_result();
}

template <typename Real>
[[nodiscard]] NestResult validate_2d_window(
    const RemapWindow& window,
    const FieldView2D<Real>& field) noexcept {
  const auto nx = active_nx(field);
  const auto ny = active_ny(field);
  if (!valid_active_shape(nx, ny, 1)) {
    return result(NestStatus::invalid_contract, "state exchange field shape is empty");
  }
  if (!window_fits_active_shape(window, nx, ny)) {
    return result(
        NestStatus::invalid_contract,
        "state exchange remap window is outside child active extents");
  }
  return ok_result();
}

}  // namespace

StateExchangePlan build_exposed_child_state_exchange_plan(
    const RemapPlan& remap_plan,
    const StateView<const float>& child_state) noexcept {
  StateExchangePlan plan{};
  plan.result = ok_result();
  plan.report.result = ok_result();
  plan.operation = ExchangeOperation::parent_to_child_interpolation;

  if (!remap_plan.ok()) {
    plan.result = remap_plan.result;
    plan.report.result = remap_plan.result;
    return plan;
  }

  if (const auto validation = validate_3d_window(remap_plan.u, child_state.u);
      !validation.ok()) {
    plan.result = validation;
    plan.report.result = validation;
    return plan;
  }
  if (const auto validation = validate_3d_window(remap_plan.v, child_state.v);
      !validation.ok()) {
    plan.result = validation;
    plan.report.result = validation;
    return plan;
  }
  if (const auto validation = validate_2d_window(remap_plan.surface, child_state.mu);
      !validation.ok()) {
    plan.result = validation;
    plan.report.result = validation;
    return plan;
  }
  if (const auto validation = validate_3d_window(remap_plan.mass, child_state.qvapor);
      !validation.ok()) {
    plan.result = validation;
    plan.report.result = validation;
    return plan;
  }
  if (const auto validation = validate_3d_window(remap_plan.mass, child_state.t);
      !validation.ok()) {
    plan.result = validation;
    plan.report.result = validation;
    return plan;
  }
  if (const auto validation = validate_3d_window(remap_plan.w_full, child_state.ph);
      !validation.ok()) {
    plan.result = validation;
    plan.report.result = validation;
    return plan;
  }

  plan.field_count = 6;
  plan.fields[0] = make_field_plan(
      StateExchangeField::u,
      HorizontalStagger::u,
      remap_plan.u,
      active_nx(child_state.u),
      active_ny(child_state.u),
      active_nz(child_state.u),
      true);
  plan.fields[1] = make_field_plan(
      StateExchangeField::v,
      HorizontalStagger::v,
      remap_plan.v,
      active_nx(child_state.v),
      active_ny(child_state.v),
      active_nz(child_state.v),
      true);
  plan.fields[2] = make_field_plan(
      StateExchangeField::mu,
      HorizontalStagger::surface,
      remap_plan.surface,
      active_nx(child_state.mu),
      active_ny(child_state.mu),
      1,
      false);
  plan.fields[3] = make_field_plan(
      StateExchangeField::qvapor,
      HorizontalStagger::mass,
      remap_plan.mass,
      active_nx(child_state.qvapor),
      active_ny(child_state.qvapor),
      active_nz(child_state.qvapor),
      true);
  plan.fields[4] = make_field_plan(
      StateExchangeField::t,
      HorizontalStagger::mass,
      remap_plan.mass,
      active_nx(child_state.t),
      active_ny(child_state.t),
      active_nz(child_state.t),
      true);
  plan.fields[5] = make_field_plan(
      StateExchangeField::ph,
      HorizontalStagger::w_full,
      remap_plan.w_full,
      active_nx(child_state.ph),
      active_ny(child_state.ph),
      active_nz(child_state.ph),
      true);

  plan.report = summarize_state_exchange_plan(plan);
  return plan;
}

StateExchangeReport summarize_state_exchange_plan(const StateExchangePlan& plan) noexcept {
  StateExchangeReport report{};
  report.result = plan.result;
  if (!plan.ok()) {
    return report;
  }

  const auto field_count =
      std::min<std::uint8_t>(plan.field_count, static_cast<std::uint8_t>(plan.fields.size()));
  for (std::uint8_t field_index = 0; field_index < field_count; ++field_index) {
    add_report_field(plan.fields[field_index], report);
  }

  return report;
}

std::string_view state_exchange_field_name(const StateExchangeField field) noexcept {
  switch (field) {
    case StateExchangeField::u:
      return "U";
    case StateExchangeField::v:
      return "V";
    case StateExchangeField::mu:
      return "MU";
    case StateExchangeField::qvapor:
      return "QVAPOR";
    case StateExchangeField::t:
      return "T";
    case StateExchangeField::ph:
      return "PH";
  }
  return "unknown";
}

std::array<StateExchangeField, 6> selected_state_exchange_fields() noexcept {
  return {
      StateExchangeField::u,
      StateExchangeField::v,
      StateExchangeField::mu,
      StateExchangeField::qvapor,
      StateExchangeField::t,
      StateExchangeField::ph,
  };
}

std::array<WrfMovingNestBaseStateExchangeCandidate, 7>
wrf_moving_nest_base_state_exchange_candidates() noexcept {
  return {
      WrfMovingNestBaseStateExchangeCandidate{"PHB", true, false, false},
      WrfMovingNestBaseStateExchangeCandidate{"MUB", true, false, false},
      WrfMovingNestBaseStateExchangeCandidate{"PB", true, false, false},
      WrfMovingNestBaseStateExchangeCandidate{"ALB", false, true, false},
      WrfMovingNestBaseStateExchangeCandidate{"T_INIT", false, true, false},
      WrfMovingNestBaseStateExchangeCandidate{"HT", false, true, false},
      WrfMovingNestBaseStateExchangeCandidate{"HGT", false, true, false},
  };
}

WrfMovingNestBaseStateExchangeContractReport
describe_wrf_moving_nest_base_state_exchange_contract() noexcept {
  WrfMovingNestBaseStateExchangeContractReport report{};
  report.active_selected_fields = selected_state_exchange_fields();
  report.active_selected_field_count =
      static_cast<std::uint8_t>(report.active_selected_fields.size());
  report.base_state_candidates =
      wrf_moving_nest_base_state_exchange_candidates();
  report.base_state_candidate_count =
      static_cast<std::uint8_t>(report.base_state_candidates.size());
  report.diagnostic_only = true;
  report.enables_selected_field_numerics = false;
  return report;
}

}  // namespace tywrf::nest
