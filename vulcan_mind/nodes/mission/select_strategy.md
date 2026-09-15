---
id: mission.select_strategy
label: select_strategy
kind: function
source:
  file: src/mission/operations/aerobraking_policy/selector_stub.jl
  symbol: select_strategy
  lines:
  - 7
  - 14
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: AerobrakingPolicy namespace supplying selectors, configuration, and
    strategy enum values.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: strategy
  type: AerobrakingStrategyKind
  units: n/a
  description: Discrete E_EDG or T_EDG strategy selected from policy configuration
    and mission input.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- mission
- aerobraking
- policy
charts:
- mission
origin: agent
---

# select_strategy

## Purpose
`select_strategy` is the policy dispatch function used by aerobraking guidance to obtain a discrete strategy for a pass. The selector stub file defines methods for the unimplemented DRL adapter and the default selector, keeping the policy contract stable while alternative policy implementations are developed.

## Theory & Math
The output is a categorical state in `{E_EDG, T_EDG}`. Selection is a policy map `π(c,x) → s`, where `c` is `AerobrakingPolicyConfig`, `x` is pass or mission input, and `s` is the selected enum. No aerodynamic or orbital force is computed at this boundary.

## Model & Assumptions
The selector assumes the input contains the mission information required by its concrete implementation. `AerobrakingPolicyConfig` supplies a default strategy, and callers assume both enum values are understood by downstream guidance dispatch. The abstract selector contract allows a future learned policy to replace the current default without changing its callers.

## Design & Implementation
`selector_stub.jl` defines `select_strategy` for `DRLPolicyAdapterStub` as an explicit not-implemented error and for `DefaultAerobrakingPolicySelector` as the configured fallback path. `guidance/aerobraking/dispatcher.jl` calls the function before selecting a guidance maneuver, so policy errors occur before thrust commands are produced.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | AerobrakingPolicy namespace supplying selectors, configuration, and strategy enum values. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `strategy` | AerobrakingStrategyKind | n/a | — | Discrete E_EDG or T_EDG strategy selected from policy configuration and mission input. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.mission|AerobrakingPolicy]] · `api` → `module_api` · call · `src/mission/operations/aerobraking_policy/selector_stub.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The DRL adapter is a stub and cannot select a learned strategy. The default selector does not evaluate atmospheric feasibility, propellant, or control authority. An enum result therefore represents policy choice only; it is not a guarantee that the selected strategy can be executed on the current spacecraft.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/selector_stub.jl:7-14`.
