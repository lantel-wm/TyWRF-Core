#pragma once

#include "tywrf/field_view.hpp"

#include <cstdint>

namespace tywrf {

struct GridConfig {
  std::int32_t mass_nx = 0;
  std::int32_t mass_ny = 0;
  std::int32_t mass_nz = 0;
  std::int32_t full_nz = 0;
  Halo3D halo{};
};

class Grid {
 public:
  constexpr Grid() = default;

  explicit constexpr Grid(const GridConfig config) noexcept : config_(config) {}

  [[nodiscard]] constexpr const GridConfig& config() const noexcept {
    return config_;
  }

  [[nodiscard]] constexpr ActiveShape3D mass_shape() const noexcept {
    return {config_.mass_nx, config_.mass_ny, config_.mass_nz};
  }

  [[nodiscard]] constexpr ActiveShape3D u_shape() const noexcept {
    return {config_.mass_nx + 1, config_.mass_ny, config_.mass_nz};
  }

  [[nodiscard]] constexpr ActiveShape3D v_shape() const noexcept {
    return {config_.mass_nx, config_.mass_ny + 1, config_.mass_nz};
  }

  [[nodiscard]] constexpr ActiveShape3D w_shape() const noexcept {
    return {config_.mass_nx, config_.mass_ny, config_.full_nz};
  }

  [[nodiscard]] constexpr ActiveShape2D surface_shape() const noexcept {
    return {config_.mass_nx, config_.mass_ny};
  }

  [[nodiscard]] constexpr FieldLayout3D mass_layout() const noexcept {
    return make_field_layout(mass_shape(), config_.halo);
  }

  [[nodiscard]] constexpr FieldLayout3D u_layout() const noexcept {
    return make_field_layout(u_shape(), config_.halo);
  }

  [[nodiscard]] constexpr FieldLayout3D v_layout() const noexcept {
    return make_field_layout(v_shape(), config_.halo);
  }

  [[nodiscard]] constexpr FieldLayout3D w_layout() const noexcept {
    return make_field_layout(w_shape(), config_.halo);
  }

  [[nodiscard]] constexpr FieldLayout2D surface_layout() const noexcept {
    return make_field_layout(surface_shape(), horizontal_halo(config_.halo));
  }

 private:
  GridConfig config_{};
};

}  // namespace tywrf
