#include "tywrf/dynamics/tendency.hpp"

namespace tywrf::dynamics {

template TendencyApplyReport apply_tendency<float>(
    FieldView2D<float>,
    FieldView2D<const float>,
    float) noexcept;
template TendencyApplyReport apply_tendency<double>(
    FieldView2D<double>,
    FieldView2D<const double>,
    double) noexcept;

template TendencyApplyReport apply_tendency<float>(
    FieldView3D<float>,
    FieldView3D<const float>,
    float) noexcept;
template TendencyApplyReport apply_tendency<double>(
    FieldView3D<double>,
    FieldView3D<const double>,
    double) noexcept;

template TendencyApplyReport apply_zero_tendency<float>(FieldView2D<float>) noexcept;
template TendencyApplyReport apply_zero_tendency<double>(FieldView2D<double>) noexcept;

template TendencyApplyReport apply_zero_tendency<float>(FieldView3D<float>) noexcept;
template TendencyApplyReport apply_zero_tendency<double>(FieldView3D<double>) noexcept;

template TendencyApplyReport apply_state_tendencies<float>(
    StateView<float>,
    StateView<const float>,
    float) noexcept;
template TendencyApplyReport apply_state_tendencies<double>(
    StateView<double>,
    StateView<const double>,
    double) noexcept;

template TendencyApplyReport apply_zero_state_tendency<float>(StateView<float>) noexcept;
template TendencyApplyReport apply_zero_state_tendency<double>(StateView<double>) noexcept;

}  // namespace tywrf::dynamics
