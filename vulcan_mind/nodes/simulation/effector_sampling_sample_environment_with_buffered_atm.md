---
id: simulation.effector_sampling_sample_environment_with_buffered_atm
label: sample_environment_with_buffered_atm
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_environment_with_buffered_atm
  lines:
  - 307
  - 307
inputs:
- id: req
  type: EffectorEnvironmentRequirements
  units: n/a
  required: true
  description: Positional argument `req`.
- id: model
  type: Any
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: EnvironmentSample
  units: n/a
  description: Return value of `sample_environment_with_buffered_atm`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# sample_environment_with_buffered_atm

## Purpose
A compatibility name for the reusable-buffer environment sampler, kept so call sites written against the earlier buffered-atmosphere API continue to work.

## Design & Implementation
A one-line `@inline` forward of all six arguments to `sample_environment_with_reusable_buffers`, returning its `EnvironmentSample`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `req` | EffectorEnvironmentRequirements | n/a | yes | Positional argument `req`. |
| in | `model` | Any | n/a | yes | Positional argument `model`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EnvironmentSample | n/a | — | Return value of `sample_environment_with_buffered_atm`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`

**Downstream**

- `callees` → [[simulation.effector_sampling_sample_environment_with_reusable_buffers|sample_environment_with_reusable_buffers]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:315-315`
<!-- vulcan:connections:end -->

## Limitations
The name promises only buffered atmosphere while the implementation now also reuses planet frame and solar buffers; a reader relying on the name underestimates what is shared.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 307.
