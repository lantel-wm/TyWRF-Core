#pragma once

#include "tywrf/nest/nest_interface.hpp"
#include "tywrf/state.hpp"

#include <array>
#include <cstdint>
#include <string_view>

namespace tywrf::nest {

enum class StateExchangeField : std::uint8_t {
  u,
  v,
  mu,
  qvapor,
};

struct ExposedChildRegion {
  HorizontalStagger stagger = HorizontalStagger::mass;
  std::int32_t child_i_begin = 0;
  std::int32_t child_j_begin = 0;
  std::int32_t extent_i = 0;
  std::int32_t extent_j = 0;
  std::int32_t active_k_count = 1;
  bool three_dimensional = false;
  bool owns_overlap = false;
  bool owns_halo = false;

  [[nodiscard]] constexpr bool empty() const noexcept {
    return extent_i <= 0 || extent_j <= 0 || active_k_count <= 0;
  }

  [[nodiscard]] constexpr std::uint64_t horizontal_cell_count() const noexcept {
    return empty() ? 0U
                   : static_cast<std::uint64_t>(extent_i) *
                         static_cast<std::uint64_t>(extent_j);
  }

  [[nodiscard]] constexpr std::uint64_t point_count() const noexcept {
    return horizontal_cell_count() * static_cast<std::uint64_t>(active_k_count);
  }
};

struct FieldStateExchangePlan {
  StateExchangeField field = StateExchangeField::u;
  HorizontalStagger stagger = HorizontalStagger::mass;
  std::int32_t active_nx = 0;
  std::int32_t active_ny = 0;
  std::int32_t active_k_count = 0;
  bool three_dimensional = false;
  bool owns_overlap = false;
  bool owns_halo = false;
  std::array<ExposedChildRegion, 4> exposed_regions{};
  std::uint8_t exposed_region_count = 0;
  std::uint64_t exposed_horizontal_cell_count = 0;
  std::uint64_t exchange_point_count = 0;
};

struct StateExchangeReport {
  NestResult result{NestStatus::ok, "ok"};
  std::uint32_t planned_field_count = 0;
  std::uint32_t active_field_count = 0;
  std::uint32_t exposed_region_count = 0;
  std::uint64_t exposed_horizontal_cell_count = 0;
  std::uint64_t exchange_point_count = 0;
  bool requires_parent_interpolation = false;
  bool performed_interpolation = false;
  bool modifies_overlap = false;
  bool modifies_halo = false;

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

struct StateExchangePlan {
  NestResult result{NestStatus::ok, "ok"};
  ExchangeOperation operation = ExchangeOperation::parent_to_child_interpolation;
  std::array<FieldStateExchangePlan, 4> fields{};
  std::uint8_t field_count = 0;
  StateExchangeReport report{};

  [[nodiscard]] constexpr bool ok() const noexcept {
    return result.ok();
  }
};

// Builds a host-side description of active child cells that became exposed by
// a moving-nest remap and therefore need future parent-to-child interpolation.
// Coordinates are zero-based active child indices. This function performs no
// interpolation and does not modify overlap or halo cells.
[[nodiscard]] StateExchangePlan build_exposed_child_state_exchange_plan(
    const RemapPlan& remap_plan,
    const StateView<const float>& child_state) noexcept;

[[nodiscard]] StateExchangeReport summarize_state_exchange_plan(
    const StateExchangePlan& plan) noexcept;

[[nodiscard]] std::string_view state_exchange_field_name(
    StateExchangeField field) noexcept;

}  // namespace tywrf::nest
