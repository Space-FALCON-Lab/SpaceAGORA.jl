---
id: simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state
label: _spacecraft_outside_atmosphere_for_current_state
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _spacecraft_outside_atmosphere_for_current_state
  lines:
  - 1942
  - 1942
inputs:
- id: sc
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc`.
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
  type: Bool
  units: n/a
  description: Return value of `_spacecraft_outside_atmosphere_for_current_state`.
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

# _spacecraft_outside_atmosphere_for_current_state

## Purpose
Tests whether a satellite is above the entry interface under a model that vanishes there, using the staged flag when current or a fresh planet-frame sample otherwise.

## Design & Implementation
Returns false unless the model vanishes above EI; uses the buffered in-atmosphere flag if its timestamp matches, else computes altitude from prefilled or fresh planet frame. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc` | Any | n/a | yes | Positional argument `sc`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_spacecraft_outside_atmosphere_for_current_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__all_active_spacecraft_outside_atmosphere|_all_active_spacecraft_outside_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1965-1965`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1343-1343`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2006-2006`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__drag_state_buffer_current|_drag_state_buffer_current]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1951-1951`
- `callees` → [[simulation.effector_sampling_sample_buffered_planet_frame|sample_buffered_planet_frame]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1956-1956`
- `callees` → [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1957-1957`
- `callees` → [[simulation.model_selection__density_model_for_sat|_density_model_for_sat]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1948-1948`
<!-- vulcan:connections:end -->

## Limitations
Fresh sampling costs a rotation per call.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1942.
