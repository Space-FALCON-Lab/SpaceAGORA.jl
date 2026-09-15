---
id: simulation.thermal_callbacks_condition
label: condition
kind: function
source:
  file: src/simulation/callbacks/thermal_callbacks.jl
  symbol: condition
  lines:
  - 66
  - 66
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `condition`. Returns `true`.
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

# condition

## Purpose
Trigger predicate of the thermal `DiscreteCallback`. It returns `true` unconditionally, which makes the thermal update fire at the end of every accepted solver step.

## Design & Implementation
Written as the one-line closure `condition(u, t, integrator) = true` inside `get_thermal_callback`. The `DiscreteCallback` contract calls this after each successful step and invokes `affect!` whenever it returns true, so an always-true condition turns the callback into a per-step hook rather than an event detector. The signature still accepts the state `u`, time `t`, and `integrator` so that a future gating rule, for instance skipping the update above the atmospheric interface, can be added without changing the callback construction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `condition`. Returns `true`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:70-70`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:247-247`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:66-66`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because it never returns false, the thermal computation runs on every step regardless of altitude, density, or whether any link actually sees heating, which costs an atmosphere-dependent evaluation per spacecraft per step even in vacuum. The callback therefore also fires on steps where the density buffers may be stale, and it offers no way to throttle the update rate.

## Provenance
Mapped from `src/simulation/callbacks/thermal_callbacks.jl` line 66.
