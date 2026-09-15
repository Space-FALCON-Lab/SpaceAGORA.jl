---
id: module.mission
label: AerobrakingPolicy
kind: module
source:
  file: src/mission/operations/aerobraking_policy/policy_types.jl
  symbol: AerobrakingPolicy
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: AerobrakingStrategyKind, policy configuration and selector types, pass
    schedules, strategy_for_pass, and select_strategy.
tags:
- module
charts:
- master
origin: agent
---

# AerobrakingPolicy

## Purpose
`AerobrakingPolicy` defines the mission-level policy boundary for selecting an entry or terminal guidance strategy on each aerobraking pass. The root file declares the strategy enum, the abstract selector contract, and the keyword-configured policy record, then includes the stub selector and pass-schedule implementations. This keeps policy selection data separate from the guidance models that consume the selected strategy.

## Theory & Math
The policy is discrete: `E_EDG` and `T_EDG` are enum values selected from mission state and configuration. A pass schedule maps pass index and schedule metadata to a strategy decision. The selector therefore behaves as a finite-state policy rather than a continuous controller; no force or torque is computed in this module.

## Model & Assumptions
`AerobrakingPolicyConfig` defaults to `E_EDG`, making the E-EDG strategy the fallback when no mission override is supplied. Selectors implement `AbstractAerobrakingPolicySelector`, and the schedule helpers assume pass identifiers are valid for the configured mission. The policy layer assumes downstream guidance understands both enum values and handles invalid operational conditions separately.

## Design & Implementation
`policy_types.jl` declares `AerobrakingStrategyKind`, `AerobrakingPolicyConfig`, and the abstract selector before including `selector_stub.jl` and `pass_schedule.jl`. Those files provide `DRLPolicyAdapterStub`, `DefaultAerobrakingPolicySelector`, `AerobrakingPassSchedule`, `strategy_for_pass`, and `select_strategy`. `SimulationModel` re-exports `AerobrakingPolicy` after IO owners, so guidance and callback code can use the same public types.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | AerobrakingStrategyKind, policy configuration and selector types, pass schedules, strategy_for_pass, and select_strategy. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[misc.pass_schedule_strategy_for_pass|strategy_for_pass]] · `module_api` · call · `src/mission/operations/aerobraking_policy/pass_schedule.jl`
- `api` → [[misc.policy_types_aerobrakingpolicyconfig|AerobrakingPolicyConfig]] · `module_api` · call · `src/mission/operations/aerobraking_policy/policy_types.jl`
- `api` → [[mission.maneuver_plans_earth_firing_plan|Earth_firing_plan]] · `module_api` · call · `src/mission/operations/maneuver_plans.jl`
- `api` → [[mission.maneuver_plans_magellan_firing_plan|Magellan_firing_plan]] · `module_api` · call · `src/mission/operations/maneuver_plans.jl`
- `api` → [[mission.maneuver_plans_odyssey_firing_plan|Odyssey_firing_plan]] · `module_api` · call · `src/mission/operations/maneuver_plans.jl`
- `api` → [[mission.maneuver_plans_odyssey_firing_plan_true_beginning|Odyssey_firing_plan_true_beginning]] · `module_api` · call · `src/mission/operations/maneuver_plans.jl`
- `api` → [[mission.maneuver_plans_titan_firing_plan|titan_firing_plan]] · `module_api` · call · `src/mission/operations/maneuver_plans.jl`
- `api` → [[mission.maneuver_plans_venus_express_firing_plan|Venus_Express_firing_plan]] · `module_api` · call · `src/mission/operations/maneuver_plans.jl`
- `api` → [[mission.pass_schedule_aerobrakingpassschedule|AerobrakingPassSchedule]] · `module_api` · call · `src/mission/operations/aerobraking_policy/pass_schedule.jl`
- `api` → [[mission.policy_types_abstractaerobrakingpolicyselector|AbstractAerobrakingPolicySelector]] · `module_api` · call · `src/mission/operations/aerobraking_policy/policy_types.jl`
- `api` → [[mission.policy_types_aerobrakingpolicy|AerobrakingPolicy]] · `module_api` · call · `src/mission/operations/aerobraking_policy/policy_types.jl`
- `api` → [[mission.select_strategy|select_strategy]] · `module_api` · call · `src/mission/operations/aerobraking_policy/selector_stub.jl`
- `api` → [[mission.selector_stub_defaultaerobrakingpolicyselector|DefaultAerobrakingPolicySelector]] · `module_api` · call · `src/mission/operations/aerobraking_policy/selector_stub.jl`
- `api` → [[mission.selector_stub_drlpolicyadapterstub|DRLPolicyAdapterStub]] · `module_api` · call · `src/mission/operations/aerobraking_policy/selector_stub.jl`
- `api` → [[module.core|SimulationModel]] · `mission` · call · `src/core/simulation_model.jl:117-117`
<!-- vulcan:connections:end -->

## Limitations
The root policy file does not validate mission geometry, atmospheric density, or control authority. The default strategy is a safe configuration fallback, not evidence that E-EDG is feasible for every pass. The stub adapter cannot represent a trained policy, and schedule lookup errors remain the responsibility of callers and selector implementations.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/policy_types.jl`, `selector_stub.jl`, and `pass_schedule.jl`.
