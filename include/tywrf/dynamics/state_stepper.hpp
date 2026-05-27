#pragma once

#include "tywrf/dynamics/dynamics_loop.hpp"
#include "tywrf/dynamics/tendency.hpp"
#include "tywrf/state.hpp"

#include <cstdint>

namespace tywrf::dynamics {

enum class StateStepStatus : std::uint8_t {
  ok,
  missing_d01_state,
  missing_d02_state,
  tendency_apply_failed,
};

struct DomainStateStepReport {
  std::int64_t zero_tendency_apply_count = 0;
  std::int64_t active_points = 0;
  TendencyApplyStatus tendency_status = TendencyApplyStatus::ok;
  bool missing_state = false;
  bool failed = false;
};

struct SkeletonStateStepReport {
  LoopSummary loop{};
  DomainStateStepReport d01{};
  DomainStateStepReport d02{};
  std::int64_t total_active_points = 0;
  bool layout_or_status_failure = false;
  StateStepStatus status = StateStepStatus::ok;
  DomainId failure_domain = DomainId::d01;
  TendencyApplyStatus tendency_status = TendencyApplyStatus::ok;
};

class SkeletonStateStepper {
 public:
  explicit SkeletonStateStepper(DynamicsLoopRunner runner);
  explicit SkeletonStateStepper(DynamicsLoopConfig config);

  [[nodiscard]] const DynamicsLoopRunner& runner() const noexcept;

  [[nodiscard]] SkeletonStateStepReport run(
      State<float>* d01_state,
      State<float>* d02_state) const;

  [[nodiscard]] SkeletonStateStepReport run(
      State<float>& d01_state,
      State<float>& d02_state) const;

 private:
  DynamicsLoopRunner runner_;
};

}  // namespace tywrf::dynamics
