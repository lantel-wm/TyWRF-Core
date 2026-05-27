#include "tywrf/dynamics/state_stepper.hpp"

#include <utility>

namespace tywrf::dynamics {
namespace {

struct StepContext {
  State<float>* d01_state = nullptr;
  State<float>* d02_state = nullptr;
  SkeletonStateStepReport report{};
};

[[nodiscard]] State<float>* state_for_domain(
    StepContext& context,
    const DomainId domain) noexcept {
  switch (domain) {
    case DomainId::d01:
      return context.d01_state;
    case DomainId::d02:
      return context.d02_state;
  }
  return nullptr;
}

[[nodiscard]] DomainStateStepReport& domain_report(
    SkeletonStateStepReport& report,
    const DomainId domain) noexcept {
  switch (domain) {
    case DomainId::d01:
      return report.d01;
    case DomainId::d02:
      return report.d02;
  }
  return report.d01;
}

[[nodiscard]] StateStepStatus missing_state_status(const DomainId domain) noexcept {
  switch (domain) {
    case DomainId::d01:
      return StateStepStatus::missing_d01_state;
    case DomainId::d02:
      return StateStepStatus::missing_d02_state;
  }
  return StateStepStatus::tendency_apply_failed;
}

void mark_failure(
    SkeletonStateStepReport& report,
    const DomainId domain,
    const StateStepStatus status,
    const TendencyApplyStatus tendency_status) noexcept {
  report.layout_or_status_failure = true;
  if (report.status == StateStepStatus::ok) {
    report.status = status;
    report.failure_domain = domain;
    report.tendency_status = tendency_status;
  }
}

void apply_zero_state_tendency_event(void* user_data, const LoopEvent& event) {
  if (event.kind != LoopEventKind::zero_dynamics_tendency) {
    return;
  }

  auto* context = static_cast<StepContext*>(user_data);
  auto& report = context->report;
  auto& domain = domain_report(report, event.domain);
  auto* state = state_for_domain(*context, event.domain);

  if (state == nullptr) {
    domain.missing_state = true;
    domain.tendency_status = TendencyApplyStatus::null_field;
    domain.failed = true;
    mark_failure(report, event.domain, missing_state_status(event.domain),
                 TendencyApplyStatus::null_field);
    return;
  }

  const auto apply_report = apply_zero_state_tendency(state->view());
  ++domain.zero_tendency_apply_count;
  if (apply_report.status != TendencyApplyStatus::ok) {
    domain.tendency_status = apply_report.status;
    domain.failed = true;
    mark_failure(report, event.domain, StateStepStatus::tendency_apply_failed,
                 apply_report.status);
    return;
  }

  domain.active_points += apply_report.active_points;
  report.total_active_points += apply_report.active_points;
}

}  // namespace

SkeletonStateStepper::SkeletonStateStepper(DynamicsLoopRunner runner)
    : runner_(std::move(runner)) {}

SkeletonStateStepper::SkeletonStateStepper(DynamicsLoopConfig config)
    : runner_(std::move(config)) {}

const DynamicsLoopRunner& SkeletonStateStepper::runner() const noexcept {
  return runner_;
}

SkeletonStateStepReport SkeletonStateStepper::run(
    State<float>* d01_state,
    State<float>* d02_state) const {
  StepContext context{d01_state, d02_state, {}};
  context.report.loop = runner_.run({&context, apply_zero_state_tendency_event});
  return context.report;
}

SkeletonStateStepReport SkeletonStateStepper::run(
    State<float>& d01_state,
    State<float>& d02_state) const {
  return run(&d01_state, &d02_state);
}

}  // namespace tywrf::dynamics
