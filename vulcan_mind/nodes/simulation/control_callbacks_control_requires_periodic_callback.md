---
id: simulation.control_callbacks_control_requires_periodic_callback
label: control_requires_periodic_callback
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: control_requires_periodic_callback
  lines:
  - 25
  - 25
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
  type: Bool
  units: n/a
  description: Return value of `control_requires_periodic_callback`. Returns `true`.
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

# control_requires_periodic_callback

## Purpose
Trait predicate that decides whether a control effector needs a fixed-rate `PeriodicCallback` or can instead be driven by event-scheduled callbacks. It returns `true` for any effector and `false` specifically for `BaseThrusterModel`.

## Design & Implementation
Implemented as two `@inline` one-line methods dispatching on the effector type: `control_requires_periodic_callback(::Any) = true` is the conservative default, and `control_requires_periodic_callback(::BaseThrusterModel) = false` opts thrusters out. `get_control_callbacks` branches on the result, routing thruster models to `_thruster_schedule_callbacks`, which computes the burn schedule once at initialization and registers tstops, while every other effector gets a `PeriodicCallback` at its configured `control_rate`. Neither method inspects the effector instance, so the answer is a compile-time property of the type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `control_requires_periodic_callback`. Returns `true`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:92-92`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the decision is purely type-based, a thruster model whose burn times are recomputed during the run, for example by closed-loop guidance, is still treated as schedulable up front and will not be re-evaluated at a periodic rate. Adding a new effector type silently inherits the `::Any` method and pays for a periodic callback even when it is event-driven.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl` line 25.
