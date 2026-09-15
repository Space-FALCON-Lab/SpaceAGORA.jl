---
id: mission.policy_types_aerobrakingpolicy
label: AerobrakingPolicy
kind: module
source:
  file: src/mission/operations/aerobraking_policy/policy_types.jl
  symbol: AerobrakingPolicy
  lines:
  - 1
  - 1
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# AerobrakingPolicy

## Purpose

`AerobrakingPolicy` is the module that declares the vocabulary for choosing an aerobraking pass strategy. It defines the `AerobrakingStrategyKind` enum with the two members `E_EDG` and `T_EDG`, the abstract supertype `AbstractAerobrakingPolicySelector`, and the keyword struct `AerobrakingPolicyConfig` whose single field `default_strategy` defaults to `E_EDG`.

## Design & Implementation

The module body is deliberately thin: after the type declarations it `include`s `selector_stub.jl` and `pass_schedule.jl` from `@__DIR__`, so the concrete selector (`DefaultAerobrakingPolicySelector`, `DRLPolicyAdapterStub`) and the per-pass schedule (`AerobrakingPassSchedule`, `strategy_for_pass`) live in sibling files while their names are re-exported from here. The explicit `export` list is the module's public surface, covering both enum members, the config struct, the abstract selector type and the `select_strategy` entry point.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.mission|AerobrakingPolicy]] · `api` → `module_api` · call · `src/mission/operations/aerobraking_policy/policy_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Only two strategy kinds exist, so any richer policy space needs an enum change plus recompilation of every dispatch site. The module carries no validation of `default_strategy` beyond the enum type, performs no I/O, and holds no state, so behaviour depends entirely on the included files. Because the included files are pulled in with `@__DIR__`, the module cannot be relocated without moving its siblings.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/policy_types.jl` line 1.
