---
id: simulation.runtime_update_density_sat_bang
label: update_density_sat!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: update_density_sat!
  lines:
  - 222
  - 222
inputs:
- id: i
  type: Int
  units: n/a
  required: true
  description: Positional argument `i`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
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
  description: Return value of `update_density_sat!`; mutates `i` in place. Returns
    `nothing`.
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

# update_density_sat!

## Purpose
Per-satellite unit of work inside the density callback: sample the atmosphere for satellite `i` at time `t` and publish it to the shared buffers.

## Design & Implementation
A closure defined inside `get_density_callback` so it captures `has_j2_effector`. It reads the run-scoped env config, resolves the density model for satellite `i`, computes kinematics from `u.sc[i]`, and calls `_density_state_from_kinematics!` with the mass, geodetic coordinates, track-cache configuration, stats flag and the J2 flag combined with whether a J2 effector is actually present. The result is written with `_write_density_buffers!` at time `t`. It is called either from a plain loop or from `threaded_foreach_persistent`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `update_density_sat!`; mutates `i` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:337-337`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:222-222`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:226-226`
- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:223-223`
- `callees` → [[simulation.model_selection__density_model_for_sat|_density_model_for_sat]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:224-224`
- `callees` → [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:227-227`
- `callees` → [[simulation.runtime__extract_mass_kg|_extract_mass_kg]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:232-232`
- `callees` → [[simulation.runtime__stage_environment_kinematics|_stage_environment_kinematics]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:226-226`
- `callees` → [[simulation.runtime__write_density_buffers_bang|_write_density_buffers!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:243-243`
<!-- vulcan:connections:end -->

## Limitations
When run threaded, every invocation shares the same `caches` vector and stats accumulator, so correctness rests on the per-satellite cache entries being disjoint and the stats updater being locked; the closure itself contains no synchronisation.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 222.
