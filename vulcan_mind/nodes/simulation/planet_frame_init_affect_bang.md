---
id: simulation.planet_frame_init_affect_bang
label: init_affect!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/planet_frame.jl
  symbol: init_affect!
  lines:
  - 84
  - 84
inputs:
- id: cb
  type: Any
  units: n/a
  required: true
  description: Positional argument `cb`.
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
  type: Any
  units: n/a
  description: Return value of `init_affect!`; mutates `cb` in place.
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

# init_affect!

## Purpose
`init_affect!(cb, u, t, integrator)` is the `initialize` hook of the planet-frame `DiscreteCallback`. It establishes the absolute ephemeris epoch of the run and performs the very first frame update, so that the planet orientation is already correct before the solver evaluates the RHS at the initial condition rather than only after the first completed step.

## Design & Implementation
Two statements. First it seeds the run epoch: `p.shared_buffers.et_start[] = ephemerides_time_seconds(p.args.initial_time, p.args.environment_model.ephemerides_model)`, converting the configured start time into backend-specific ephemeris seconds and storing it in the shared `Ref` that `_planet_lpi_at` later adds the relative solver time to. Second it calls `affect!(integrator)` directly, whose in-line comment states the intent — initialise the planet frame at the start of the simulation. Because `affect!` is a closure over the same scope, the initialisation path and the per-step path are literally the same code, so the two can never drift apart. The `cb` and `u` arguments exist only to satisfy the DifferentialEquations.jl `initialize` signature.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cb` | Any | n/a | yes | Positional argument `cb`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `init_affect!`; mutates `cb` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:84-84`

**Downstream**

- `callees` → [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:86-86`
- `callees` → [[simulation.event_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:87-87`
- `callees` → [[simulation.planet_frame_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:87-87`
- `callees` → [[simulation.runtime_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:87-87`
- `callees` → [[simulation.thermal_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:87-87`
<!-- vulcan:connections:end -->

## Limitations
There is an ordering hazard: `et_start[]` is only valid once this hook has run, so any callback or preprocessing step that queries the planet frame earlier in `CallbackSet` initialisation order gets an epoch of zero and a silently wrong orientation. Conversion happens once, so a run whose `initial_time` is later mutated in `p.args` will keep the stale epoch for its entire duration. No validation confirms that the resulting `et` falls inside the loaded ephemeris coverage, so an out-of-range start time fails later, inside the backend, rather than here where the cause is obvious.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/planet_frame.jl` line 84.
