---
id: simulation.event_callbacks_get_quaternion_projection_callback
label: get_quaternion_projection_callback
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: get_quaternion_projection_callback
  lines:
  - 193
  - 193
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: Union{Bool, DiscreteCallback}
  units: n/a
  description: Return value of `get_quaternion_projection_callback`. Returns `true`
    or `false` or `DiscreteCallback(`.
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

# get_quaternion_projection_callback

## Purpose
Constructs a `DiscreteCallback` that projects every active satellite's attitude quaternion back onto the unit sphere after each accepted integrator step, so numerical drift in `|q|` from integrating kinematics does not accumulate into attitude error over long missions.

## Design & Implementation
Takes `num_sats::Int` and `args::SimulationConfiguration`. `correction_tol = max(32 * eps(Float64), args.integration_tolerances.abstol_quaternion)` sets the drift threshold for logging. The nested `condition` returns `true` while any satellite is active. `affect!` iterates active satellites, fetches `q = _state_quaternion(u, i)` (skipping `nothing` for layouts without attitude), computes `qnorm2 = dot(q, q)`, marks `corrected` when the norm is non-finite, below `eps`, or off by more than `correction_tol`, and unconditionally writes `u.sc[i].q .= project_unit_quaternion(q)` in place. The callback is built with `initialize = (cb, u, t, integrator) -> affect!(integrator)` so the initial state is projected before the first step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Bool, DiscreteCallback} | n/a | — | Return value of `get_quaternion_projection_callback`. Returns `true` or `false` or `DiscreteCallback(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:186-186`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Projection is applied every step even when drift is negligible, so the cost is one normalisation per satellite per step. The in-place write assumes the state layout exposes `u.sc[i].q`, which differs from the accessor `_state_quaternion` used for reading; a gravity-backbone `ArrayPartition` layout would need `u.x[1].sc[i].q` and would throw here. Projecting after the step, not during, means the integrator's error estimate was computed on the unnormalised quaternion. A zero or NaN quaternion is passed to `project_unit_quaternion` and its handling depends on that function.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 193.
