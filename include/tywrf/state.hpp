#pragma once

#include "tywrf/field_view.hpp"
#include "tywrf/grid.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace tywrf {

template <typename T>
class FieldStorage2D {
 public:
  using value_type = T;

  FieldStorage2D() = default;

  explicit FieldStorage2D(const FieldLayout2D layout)
      : layout_(layout), values_(checked_size(layout)) {}

  FieldStorage2D(const ActiveShape2D active, const Halo2D halo)
      : FieldStorage2D(make_field_layout(active, halo)) {}

  [[nodiscard]] FieldView2D<T> view() noexcept {
    return make_field_view(values_.data(), layout_);
  }

  [[nodiscard]] FieldView2D<const T> view() const noexcept {
    return make_field_view(values_.data(), layout_);
  }

  [[nodiscard]] const FieldLayout2D& layout() const noexcept {
    return layout_;
  }

  [[nodiscard]] T* data() noexcept {
    return values_.data();
  }

  [[nodiscard]] const T* data() const noexcept {
    return values_.data();
  }

  [[nodiscard]] std::size_t size() const noexcept {
    return values_.size();
  }

 private:
  [[nodiscard]] static std::size_t checked_size(const FieldLayout2D layout) {
    if (!layout.valid()) {
      throw std::invalid_argument("invalid 2D field layout");
    }
    return layout.allocation_size();
  }

  FieldLayout2D layout_{};
  std::vector<T> values_;
};

template <typename T>
class FieldStorage3D {
 public:
  using value_type = T;

  FieldStorage3D() = default;

  explicit FieldStorage3D(const FieldLayout3D layout)
      : layout_(layout), values_(checked_size(layout)) {}

  FieldStorage3D(const ActiveShape3D active, const Halo3D halo)
      : FieldStorage3D(make_field_layout(active, halo)) {}

  [[nodiscard]] FieldView3D<T> view() noexcept {
    return make_field_view(values_.data(), layout_);
  }

  [[nodiscard]] FieldView3D<const T> view() const noexcept {
    return make_field_view(values_.data(), layout_);
  }

  [[nodiscard]] const FieldLayout3D& layout() const noexcept {
    return layout_;
  }

  [[nodiscard]] T* data() noexcept {
    return values_.data();
  }

  [[nodiscard]] const T* data() const noexcept {
    return values_.data();
  }

  [[nodiscard]] std::size_t size() const noexcept {
    return values_.size();
  }

 private:
  [[nodiscard]] static std::size_t checked_size(const FieldLayout3D layout) {
    if (!layout.valid()) {
      throw std::invalid_argument("invalid 3D field layout");
    }
    return layout.allocation_size();
  }

  FieldLayout3D layout_{};
  std::vector<T> values_;
};

template <typename Real>
struct StateView {
  FieldView3D<Real> u;
  FieldView3D<Real> v;
  FieldView3D<Real> w;
  FieldView3D<Real> ph;
  FieldView3D<Real> phb;
  FieldView3D<Real> t;
  FieldView3D<Real> p;
  FieldView3D<Real> pb;
  FieldView3D<Real> qvapor;
  FieldView3D<Real> qcloud;
  FieldView3D<Real> qrain;
  FieldView3D<Real> qice;
  FieldView3D<Real> qsnow;
  FieldView3D<Real> qgraup;
  FieldView3D<Real> qnice;
  FieldView3D<Real> qnrain;

  FieldView2D<Real> mu;
  FieldView2D<Real> mub;
  FieldView2D<Real> psfc;
  FieldView2D<Real> u10;
  FieldView2D<Real> v10;
  FieldView2D<Real> t2;
  FieldView2D<Real> q2;
  FieldView2D<Real> rainc;
  FieldView2D<Real> rainnc;
};

template <typename Real = double>
struct State {
  explicit State(const Grid& grid)
      : grid(grid),
        u(grid.u_layout()),
        v(grid.v_layout()),
        w(grid.w_layout()),
        ph(grid.w_layout()),
        phb(grid.w_layout()),
        t(grid.mass_layout()),
        p(grid.mass_layout()),
        pb(grid.mass_layout()),
        qvapor(grid.mass_layout()),
        qcloud(grid.mass_layout()),
        qrain(grid.mass_layout()),
        qice(grid.mass_layout()),
        qsnow(grid.mass_layout()),
        qgraup(grid.mass_layout()),
        qnice(grid.mass_layout()),
        qnrain(grid.mass_layout()),
        mu(grid.surface_layout()),
        mub(grid.surface_layout()),
        psfc(grid.surface_layout()),
        u10(grid.surface_layout()),
        v10(grid.surface_layout()),
        t2(grid.surface_layout()),
        q2(grid.surface_layout()),
        rainc(grid.surface_layout()),
        rainnc(grid.surface_layout()) {}

  [[nodiscard]] StateView<Real> view() noexcept {
    return {
        u.view(),     v.view(),      w.view(),      ph.view(),    phb.view(),
        t.view(),     p.view(),      pb.view(),     qvapor.view(), qcloud.view(),
        qrain.view(), qice.view(),   qsnow.view(),  qgraup.view(), qnice.view(),
        qnrain.view(), mu.view(),    mub.view(),    psfc.view(),   u10.view(),
        v10.view(),   t2.view(),     q2.view(),     rainc.view(),  rainnc.view()};
  }

  [[nodiscard]] StateView<const Real> view() const noexcept {
    return {
        u.view(),     v.view(),      w.view(),      ph.view(),    phb.view(),
        t.view(),     p.view(),      pb.view(),     qvapor.view(), qcloud.view(),
        qrain.view(), qice.view(),   qsnow.view(),  qgraup.view(), qnice.view(),
        qnrain.view(), mu.view(),    mub.view(),    psfc.view(),   u10.view(),
        v10.view(),   t2.view(),     q2.view(),     rainc.view(),  rainnc.view()};
  }

  Grid grid;

  FieldStorage3D<Real> u;
  FieldStorage3D<Real> v;
  FieldStorage3D<Real> w;
  FieldStorage3D<Real> ph;
  FieldStorage3D<Real> phb;
  FieldStorage3D<Real> t;
  FieldStorage3D<Real> p;
  FieldStorage3D<Real> pb;
  FieldStorage3D<Real> qvapor;
  FieldStorage3D<Real> qcloud;
  FieldStorage3D<Real> qrain;
  FieldStorage3D<Real> qice;
  FieldStorage3D<Real> qsnow;
  FieldStorage3D<Real> qgraup;
  FieldStorage3D<Real> qnice;
  FieldStorage3D<Real> qnrain;

  FieldStorage2D<Real> mu;
  FieldStorage2D<Real> mub;
  FieldStorage2D<Real> psfc;
  FieldStorage2D<Real> u10;
  FieldStorage2D<Real> v10;
  FieldStorage2D<Real> t2;
  FieldStorage2D<Real> q2;
  FieldStorage2D<Real> rainc;
  FieldStorage2D<Real> rainnc;
};

}  // namespace tywrf
