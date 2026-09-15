---
id: mission.selector_stub_defaultaerobrakingpolicyselector
label: DefaultAerobrakingPolicySelector
kind: struct
source:
  file: src/mission/operations/aerobraking_policy/selector_stub.jl
  symbol: DefaultAerobrakingPolicySelector
  lines:
  - 4
  - 4
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: DefaultAerobrakingPolicySelector
  units: n/a
  description: Constructed `DefaultAerobrakingPolicySelector`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- mission
charts:
- mission
origin: agent
---

# DefaultAerobrakingPolicySelector

## Purpose

`DefaultAerobrakingPolicySelector` is the fallback aerobraking policy selector: a fieldless struct subtyping `AbstractAerobrakingPolicySelector` whose `select_strategy` method ignores the runtime `input` and returns `config.default_strategy` from the supplied `AerobrakingPolicyConfig`.

## Design & Implementation

The type carries no state, so instances are singletons for dispatch purposes and cost nothing to construct or copy. It participates in the same `select_strategy(selector, config::AerobrakingPolicyConfig, input)` contract as any learned policy adapter, which lets the mission layer swap in a trained selector later without changing call sites. This is the selector used whenever no adaptive policy has been configured.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | DefaultAerobrakingPolicySelector | n/a | — | Constructed `DefaultAerobrakingPolicySelector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:204-204`
- [[module.mission|AerobrakingPolicy]] · `api` → `module_api` · call · `src/mission/operations/aerobraking_policy/selector_stub.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The selection is entirely open loop: corridor state, measured density, heat rate and remaining orbits carried in `input` are all discarded, so the returned strategy never adapts within a mission. It also performs no validation that `config.default_strategy` is populated or admissible for the current mission phase.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/selector_stub.jl` line 4.
