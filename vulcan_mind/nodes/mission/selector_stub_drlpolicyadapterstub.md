---
id: mission.selector_stub_drlpolicyadapterstub
label: DRLPolicyAdapterStub
kind: struct
source:
  file: src/mission/operations/aerobraking_policy/selector_stub.jl
  symbol: DRLPolicyAdapterStub
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
  type: DRLPolicyAdapterStub
  units: n/a
  description: Constructed `DRLPolicyAdapterStub`.
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

# DRLPolicyAdapterStub

## Purpose

`DRLPolicyAdapterStub` is the placeholder subtype of `AbstractAerobrakingPolicySelector` reserved for a deep-reinforcement-learning aerobraking policy. Its `select_strategy` method is deliberately unimplemented and always raises `ErrorException("Not implemented: DRLPolicyAdapterStub.select_strategy")`.

## Design & Implementation

Like the default selector it is a fieldless struct, so it can be constructed and passed through configuration plumbing to prove that the selector abstraction is wired end to end before any learned policy exists. Failing loudly on the first call, rather than quietly falling back to `config.default_strategy`, makes an accidental selection of the stub impossible to miss in a mission run.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | DRLPolicyAdapterStub | n/a | — | Constructed `DRLPolicyAdapterStub`. |
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

Selecting this adapter aborts the run at the first policy decision, so it is unusable in any live simulation. The struct declares no observation-space, action-space or model-artifact fields, meaning the eventual implementation will need a different constructor signature and callers that build it today will have to change.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/selector_stub.jl` line 1.
