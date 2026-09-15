---
id: simulation.effector_sampling__sample_reusable_planet_frame
label: _sample_reusable_planet_frame
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _sample_reusable_planet_frame
  lines:
  - 258
  - 258
inputs:
- id: req
  type: EffectorEnvironmentRequirements
  units: n/a
  required: true
  description: Positional argument `req`.
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
  type: Nothing
  units: n/a
  description: Return value of `_sample_reusable_planet_frame`. Returns `p.shared_buffers.rhs_planet_frame_prefilled[]
    ?` or `nothing`.
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

# _sample_reusable_planet_frame

## Purpose
Chooses between the prefilled planet-frame buffer and a fresh computation, according to what the effector requires and whether prefill ran.

## Design & Implementation
Returns `nothing` unless the requirements ask for a planet frame or atmosphere. Otherwise returns `sample_buffered_planet_frame` when `rhs_planet_frame_prefilled[]` is set and `sample_planet_frame` when it is not. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `req` | EffectorEnvironmentRequirements | n/a | yes | Positional argument `req`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_sample_reusable_planet_frame`. Returns `p.shared_buffers.rhs_planet_frame_prefilled[] ?` or `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.effector_sampling_sample_environment_with_reusable_buffers|sample_environment_with_reusable_buffers]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:294-294`

**Downstream**

- `callees` → [[simulation.effector_sampling_sample_buffered_planet_frame|sample_buffered_planet_frame]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:261-261`
- `callees` → [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:262-262`
<!-- vulcan:connections:end -->

## Limitations
The atmosphere requirement implies a planet frame even when the effector did not ask for one, so the frame is computed and then dropped from the returned `EnvironmentSample` by the caller.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 258.
