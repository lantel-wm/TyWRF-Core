#pragma once

#include <cstddef>
#include <cstdint>

namespace tywrf {

struct Halo2D {
  std::int32_t i_lower = 0;
  std::int32_t i_upper = 0;
  std::int32_t j_lower = 0;
  std::int32_t j_upper = 0;
};

struct Halo3D {
  std::int32_t i_lower = 0;
  std::int32_t i_upper = 0;
  std::int32_t j_lower = 0;
  std::int32_t j_upper = 0;
  std::int32_t k_lower = 0;
  std::int32_t k_upper = 0;
};

[[nodiscard]] constexpr Halo2D horizontal_halo(const Halo3D halo) noexcept {
  return {halo.i_lower, halo.i_upper, halo.j_lower, halo.j_upper};
}

[[nodiscard]] constexpr Halo2D uniform_halo_2d(const std::int32_t width) noexcept {
  return {width, width, width, width};
}

[[nodiscard]] constexpr Halo3D uniform_halo_3d(const std::int32_t width) noexcept {
  return {width, width, width, width, width, width};
}

struct ActiveShape2D {
  std::int32_t nx = 0;
  std::int32_t ny = 0;
};

struct ActiveShape3D {
  std::int32_t nx = 0;
  std::int32_t ny = 0;
  std::int32_t nz = 0;
};

[[nodiscard]] constexpr std::size_t canonical_index_3d(
    const std::int32_t i,
    const std::int32_t j,
    const std::int32_t k,
    const std::int32_t nx,
    const std::int32_t nz) noexcept {
  return ((static_cast<std::size_t>(j) * static_cast<std::size_t>(nz)) +
          static_cast<std::size_t>(k)) *
             static_cast<std::size_t>(nx) +
         static_cast<std::size_t>(i);
}

struct FieldLayout2D {
  std::int32_t nx = 0;
  std::int32_t ny = 0;
  std::int32_t stride_i = 1;
  std::int32_t stride_j = 0;
  Halo2D halo{};

  [[nodiscard]] constexpr std::int32_t active_nx() const noexcept {
    return nx - halo.i_lower - halo.i_upper;
  }

  [[nodiscard]] constexpr std::int32_t active_ny() const noexcept {
    return ny - halo.j_lower - halo.j_upper;
  }

  [[nodiscard]] constexpr std::int32_t i_begin() const noexcept {
    return halo.i_lower;
  }

  [[nodiscard]] constexpr std::int32_t i_end() const noexcept {
    return nx - halo.i_upper;
  }

  [[nodiscard]] constexpr std::int32_t j_begin() const noexcept {
    return halo.j_lower;
  }

  [[nodiscard]] constexpr std::int32_t j_end() const noexcept {
    return ny - halo.j_upper;
  }

  [[nodiscard]] constexpr bool valid() const noexcept {
    return nx >= 0 && ny >= 0 && stride_i == 1 && stride_j >= 0 &&
           halo.i_lower >= 0 && halo.i_upper >= 0 && halo.j_lower >= 0 &&
           halo.j_upper >= 0 && active_nx() >= 0 && active_ny() >= 0;
  }

  [[nodiscard]] constexpr std::size_t index(
      const std::int32_t i,
      const std::int32_t j) const noexcept {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(stride_j) +
           static_cast<std::size_t>(i) * static_cast<std::size_t>(stride_i);
  }

  [[nodiscard]] constexpr std::size_t allocation_size() const noexcept {
    return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  }
};

struct FieldLayout3D {
  std::int32_t nx = 0;
  std::int32_t ny = 0;
  std::int32_t nz = 0;
  std::int32_t stride_i = 1;
  std::int32_t stride_k = 0;
  std::int32_t stride_j = 0;
  Halo3D halo{};

  [[nodiscard]] constexpr std::int32_t active_nx() const noexcept {
    return nx - halo.i_lower - halo.i_upper;
  }

  [[nodiscard]] constexpr std::int32_t active_ny() const noexcept {
    return ny - halo.j_lower - halo.j_upper;
  }

  [[nodiscard]] constexpr std::int32_t active_nz() const noexcept {
    return nz - halo.k_lower - halo.k_upper;
  }

  [[nodiscard]] constexpr std::int32_t i_begin() const noexcept {
    return halo.i_lower;
  }

  [[nodiscard]] constexpr std::int32_t i_end() const noexcept {
    return nx - halo.i_upper;
  }

  [[nodiscard]] constexpr std::int32_t j_begin() const noexcept {
    return halo.j_lower;
  }

  [[nodiscard]] constexpr std::int32_t j_end() const noexcept {
    return ny - halo.j_upper;
  }

  [[nodiscard]] constexpr std::int32_t k_begin() const noexcept {
    return halo.k_lower;
  }

  [[nodiscard]] constexpr std::int32_t k_end() const noexcept {
    return nz - halo.k_upper;
  }

  [[nodiscard]] constexpr bool valid() const noexcept {
    return nx >= 0 && ny >= 0 && nz >= 0 && stride_i == 1 && stride_k == nx &&
           stride_j == nx * nz && halo.i_lower >= 0 && halo.i_upper >= 0 &&
           halo.j_lower >= 0 && halo.j_upper >= 0 && halo.k_lower >= 0 &&
           halo.k_upper >= 0 && active_nx() >= 0 && active_ny() >= 0 &&
           active_nz() >= 0;
  }

  [[nodiscard]] constexpr std::size_t index(
      const std::int32_t i,
      const std::int32_t j,
      const std::int32_t k) const noexcept {
    return canonical_index_3d(i, j, k, nx, nz);
  }

  [[nodiscard]] constexpr std::size_t allocation_size() const noexcept {
    return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
           static_cast<std::size_t>(nz);
  }
};

[[nodiscard]] constexpr FieldLayout2D make_field_layout(
    const ActiveShape2D active,
    const Halo2D halo) noexcept {
  const std::int32_t nx = active.nx + halo.i_lower + halo.i_upper;
  const std::int32_t ny = active.ny + halo.j_lower + halo.j_upper;
  return {nx, ny, 1, nx, halo};
}

[[nodiscard]] constexpr FieldLayout3D make_field_layout(
    const ActiveShape3D active,
    const Halo3D halo) noexcept {
  const std::int32_t nx = active.nx + halo.i_lower + halo.i_upper;
  const std::int32_t ny = active.ny + halo.j_lower + halo.j_upper;
  const std::int32_t nz = active.nz + halo.k_lower + halo.k_upper;
  return {nx, ny, nz, 1, nx, nx * nz, halo};
}

template <typename T>
struct FieldView2D {
  T* data = nullptr;
  std::int32_t nx = 0;
  std::int32_t ny = 0;
  std::int32_t stride_i = 1;
  std::int32_t stride_j = 0;
  Halo2D halo{};

  [[nodiscard]] constexpr std::size_t index(
      const std::int32_t i,
      const std::int32_t j) const noexcept {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(stride_j) +
           static_cast<std::size_t>(i) * static_cast<std::size_t>(stride_i);
  }

  [[nodiscard]] constexpr T& operator()(
      const std::int32_t i,
      const std::int32_t j) const noexcept {
    return data[index(i, j)];
  }
};

template <typename T>
struct FieldView3D {
  T* data = nullptr;
  std::int32_t nx = 0;
  std::int32_t ny = 0;
  std::int32_t nz = 0;
  std::int32_t stride_i = 1;
  std::int32_t stride_k = 0;
  std::int32_t stride_j = 0;
  Halo3D halo{};

  [[nodiscard]] constexpr std::size_t index(
      const std::int32_t i,
      const std::int32_t j,
      const std::int32_t k) const noexcept {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(stride_j) +
           static_cast<std::size_t>(k) * static_cast<std::size_t>(stride_k) +
           static_cast<std::size_t>(i) * static_cast<std::size_t>(stride_i);
  }

  [[nodiscard]] constexpr T& operator()(
      const std::int32_t i,
      const std::int32_t j,
      const std::int32_t k) const noexcept {
    return data[index(i, j, k)];
  }
};

template <typename T>
using FieldView = FieldView3D<T>;

template <typename T>
[[nodiscard]] constexpr FieldView2D<T> make_field_view(
    T* data,
    const FieldLayout2D layout) noexcept {
  return {data, layout.nx, layout.ny, layout.stride_i, layout.stride_j, layout.halo};
}

template <typename T>
[[nodiscard]] constexpr FieldView3D<T> make_field_view(
    T* data,
    const FieldLayout3D layout) noexcept {
  return {
      data,
      layout.nx,
      layout.ny,
      layout.nz,
      layout.stride_i,
      layout.stride_k,
      layout.stride_j,
      layout.halo};
}

}  // namespace tywrf
