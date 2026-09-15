---
id: simulation_a.planet_frame_update_planet_frame_callback
label: update_planet_frame_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/planet_frame.jl
  symbol: update_planet_frame_callback
  lines:
  - 69
  - 91
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: solver_step
  type: DEIntegrator
  units: n/a
  required: true
  description: The integrator at an accepted step, supplying the current time and
    the parameter object holding the planet model and shared buffers.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: planet_frame_callback
  type: DiscreteCallback
  units: n/a
  description: Unconditional per-step callback that refreshes `planet.L_PI` and invalidates
    the RHS execution-plan step cache.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# update_planet_frame_callback

## Purpose
`update_planet_frame_callback` returns the `DiscreteCallback` that keeps the planet-fixed rotation matrix synchronised with simulation time. Every other callback and every atmosphere query that converts inertial state to altitude, latitude and longitude depends on this matrix being current for the step being evaluated.

## Theory & Math
The planet-fixed transform $L_{PI}$ maps inertial J2000 axes into the body-fixed frame. Position transforms directly, $r_{pp} = L_{PI} r_{ii}$, while velocity requires the transport term evaluated **after** rotation:

$$v_{pp} = L_{PI} v_{ii} - \omega_p \times r_{pp}$$

because $\omega_p$ is expressed in the planet-fixed frame, where the spin axis is $+z$. Taking the cross product in J2000 axes instead points the co-rotation velocity along the J2000 pole and introduces an error of order $|\omega_p| |r|$ — roughly 250 m/s at Mars periapsis.

## Model & Assumptions
The callback condition is unconditionally true, so the affect runs exactly once per accepted solver step. `_planet_lpi_at(p, t)` resolves the transform from the configured ephemeris model, and the result is written in place into `p.args.environment_model.planet.L_PI` so downstream consumers see the update without any handoff. The initializer additionally seeds `p.shared_buffers.et_start` from the configured initial epoch, converting it to ephemeris seconds once, and then runs the affect so the frame is valid before the first derivative evaluation.

## Design & Implementation
The callback doubles as the invalidation point for the RHS execution-plan step cache. Because it already runs exactly once per accepted step, unconditionally, on every non-backbone-mode solve, it is precisely the boundary that cache needs; clearing `p.shared_buffers.rhs_plan_step_cache` here avoids adding a second `CallbackSet` entry. When the step cache is disabled the extra work is a single boolean check. The file also defines `_planet_frame_state`, which applies the rotation and the transport term together for callers converting full inertial state.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `solver_step` | DEIntegrator | n/a | yes | The integrator at an accepted step, supplying the current time and the parameter object holding the planet model and shared buffers. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `planet_frame_callback` | DiscreteCallback | n/a | — | Unconditional per-step callback that refreshes `planet.L_PI` and invalidates the RHS execution-plan step cache. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:156-156`

**Downstream**

- `callees` → [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:86-86`
- `callees` → [[simulation.event_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:71-71`
- `callees` → [[simulation.event_callbacks_condition|condition]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:70-70`
- `callees` → [[simulation.planet_frame__planet_lpi_at|_planet_lpi_at]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:73-73`
- `callees` → [[simulation.planet_frame_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:71-71`
- `callees` → [[simulation.planet_frame_condition|condition]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:70-70`
- `callees` → [[simulation.planet_frame_init_affect_bang|init_affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:84-84`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:79-79`
- `callees` → [[simulation.runtime_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:71-71`
- `callees` → [[simulation.runtime_condition|condition]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:70-70`
- `callees` → [[simulation.thermal_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:71-71`
- `callees` → [[simulation.thermal_callbacks_condition|condition]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:70-70`
<!-- vulcan:connections:end -->

## Limitations
In `:gravity_backbone_split` solver mode `get_callbacks` omits this callback, so that path must refresh the frame itself. Writing into the shared planet model means the matrix is global run state rather than per-spacecraft, which is correct for a single central body but leaves no room for per-spacecraft frame variants. Rejected solver steps do not trigger the callback, so the matrix can briefly lag a trial step's time.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/planet_frame.jl:68-91`.
