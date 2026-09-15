---
id: misc.policy_types_aerobrakingpolicyconfig
label: AerobrakingPolicyConfig
kind: struct
source:
  file: src/mission/operations/aerobraking_policy/policy_types.jl
  symbol: AerobrakingPolicyConfig
  lines:
  - 7
  - 9
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: AerobrakingPolicy module scope defining the AerobrakingStrategyKind
    enum used as the field type.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: config
  type: AerobrakingPolicyConfig
  units: n/a
  description: Keyword-constructed policy configuration carrying the default strategy
    used when no per-pass override exists.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- mission
- aerobraking
- configuration
charts:
- misc
origin: agent
---

# AerobrakingPolicyConfig

## Purpose
`AerobrakingPolicyConfig` is the campaign-level configuration record for aerobraking guidance. It names the strategy a run uses when nothing more specific applies, giving the policy layer a typed default instead of a bare enum threaded through call signatures or read from a global. It is the counterpart to `AerobrakingPassSchedule`: the config supplies the baseline, the schedule supplies the exceptions.

## Model & Assumptions
The struct declares one field, `default_strategy::AerobrakingStrategyKind`, defaulting to `E_EDG`. The two admissible values come from the module's `@enum AerobrakingStrategyKind E_EDG T_EDG`, distinguishing the energy-based entry/descent guidance formulation from the time-based one. Because the field is enum-typed, an invalid strategy cannot be represented — construction fails at the type level rather than deep inside the guidance loop. `Base.@kwdef` means every construction site is explicit about what it overrides, and a bare `AerobrakingPolicyConfig()` is a valid, meaningful default configuration.

## Design & Implementation
The type sits in `policy_types.jl` alongside the enum and the `AbstractAerobrakingPolicySelector` abstract type, which together form the module's vocabulary before any behaviour is included. The file then includes `selector_stub.jl` and `pass_schedule.jl`, so the selector implementations and the pass schedule are defined against types that already exist. All of `AerobrakingStrategyKind`, `E_EDG`, `T_EDG`, `AbstractAerobrakingPolicySelector`, `AerobrakingPolicyConfig`, `DRLPolicyAdapterStub`, `DefaultAerobrakingPolicySelector`, `AerobrakingPassSchedule`, `strategy_for_pass` and `select_strategy` are exported, so consumers see one flat policy API.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | AerobrakingPolicy module scope defining the AerobrakingStrategyKind enum used as the field type. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `config` | AerobrakingPolicyConfig | n/a | — | Keyword-constructed policy configuration carrying the default strategy used when no per-pass override exists. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:205-205`
- [[module.mission|AerobrakingPolicy]] · `api` → `module_api` · call · `src/mission/operations/aerobraking_policy/policy_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The configuration holds only the default strategy; corridor bounds, drag limits, thermal margins and abort criteria live elsewhere, so it cannot be treated as a complete description of a campaign's policy. It carries no reference to a pass schedule, so callers must pass both objects and keep them consistent themselves. Adding a strategy to the enum is a breaking change for any consumer that exhaustively branches on the two current values.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/policy_types.jl:7-9`, in the module declared at lines 1-24.
