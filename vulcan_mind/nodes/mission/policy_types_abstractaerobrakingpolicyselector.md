---
id: mission.policy_types_abstractaerobrakingpolicyselector
label: AbstractAerobrakingPolicySelector
kind: struct
source:
  file: src/mission/operations/aerobraking_policy/policy_types.jl
  symbol: AbstractAerobrakingPolicySelector
  lines:
  - 5
  - 5
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
  type: AbstractAerobrakingPolicySelector
  units: n/a
  description: Abstract supertype `AbstractAerobrakingPolicySelector`; no fields.
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

# AbstractAerobrakingPolicySelector

## Purpose

`AbstractAerobrakingPolicySelector` is the abstract supertype every aerobraking policy selector subtypes. It exists so that `select_strategy` can dispatch on the selector concrete type and return an `AerobrakingStrategyKind` (`E_EDG` or `T_EDG`) without the caller knowing whether the decision came from the fixed-rule `DefaultAerobrakingPolicySelector` or the learned `DRLPolicyAdapterStub`.

## Design & Implementation

Declared with a bare `abstract type ... end` in `policy_types.jl` and exported from the `AerobrakingPolicy` module. It carries no fields and no default methods; the contract is by convention, namely that each subtype implements `select_strategy` and that the returned value is a member of the `AerobrakingStrategyKind` enum. Concrete subtypes are defined in the included `selector_stub.jl`, keeping the type declaration free of any dependency on the implementations.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractAerobrakingPolicySelector | n/a | — | Abstract supertype `AbstractAerobrakingPolicySelector`; no fields. |
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

Nothing enforces that a subtype actually implements `select_strategy`; a missing method surfaces only as a `MethodError` at call time. The abstract type also fixes no interface for selector state, warm-up or reset, so stateful learned selectors must manage their own mutability and thread-safety.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/policy_types.jl` line 5.
