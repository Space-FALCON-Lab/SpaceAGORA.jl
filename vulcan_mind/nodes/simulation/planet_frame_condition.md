---
id: simulation.planet_frame_condition
label: condition
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/planet_frame.jl
  symbol: condition
  lines:
  - 70
  - 70
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
`condition(u, t, integrator)` is the trigger predicate of the planet-frame `DiscreteCallback` built by `update_planet_frame_callback`. It exists to make the frame refresh unconditional: the planet's orientation changes continuously, so there is no event to detect and the callback must fire on every accepted integration step.

## Design & Implementation
The whole definition is `condition(u, t, integrator) = true`, a closure captured inside `update_planet_frame_callback` and handed to `DiscreteCallback(condition, affect!, initialize=init_affect!)`. Returning a literal `true` is how DifferentialEquations.jl expresses a per-step callback; because the return is a compile-time constant, the solver's branch on it is eliminated and the only per-step cost is the `affect!` body itself. The state `u`, time `t` and `integrator` arguments are accepted to satisfy the `DiscreteCallback` interface and are all unused.

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
Firing unconditionally means the 3x3 rotation is recomputed at every accepted step even when the step is very short and the planet has barely turned, which is pure overhead for a fine-stepping entry solve. Conversely, the frame is only refreshed at step boundaries, so any evaluation inside a step — including every intermediate Runge-Kutta stage — uses the orientation from the start of the step; that staleness is bounded by the step size and by the planet's rotation rate, and is not corrected anywhere. Because the predicate never returns `false`, the callback cannot be disabled at run time without rebuilding the `CallbackSet`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/planet_frame.jl` line 70.
