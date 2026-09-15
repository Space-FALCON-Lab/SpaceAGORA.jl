---
id: simulation.assembly__requires_guidance_orbit_counter
label: _requires_guidance_orbit_counter
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_guidance_orbit_counter
  lines:
  - 46
  - 46
inputs:
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
  type: Bool
  units: n/a
  description: Return value of `_requires_guidance_orbit_counter`.
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

# _requires_guidance_orbit_counter

## Purpose
Detects whether any configured guidance effector schedules its manoeuvres by orbit number, in which case the simulation must maintain a running orbit count.

## Design & Implementation
Loops over `args.guidance_model.guidance_effectors` under `@inbounds` and returns `true` at the first model for which `hasproperty(guidance_model, :maneuver_orbit_number)` holds, otherwise `false`. Using `hasproperty` rather than an `isa` test against a list of guidance types is a duck-typing choice: any guidance model that exposes that field opts into orbit counting automatically, without this file needing to know its type. It feeds `_requires_orbit_end_callback`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_guidance_orbit_counter`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- [[simulation.assembly__requires_orbit_end_callback|_requires_orbit_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:55-55`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`hasproperty` reflects on the property surface at run time and is not constant-folded as reliably as an `isa` test, so this loop can cost real work on every call. The presence of the field says nothing about its value: a guidance model carrying `maneuver_orbit_number = 0` or a disabled manoeuvre still forces the orbit-end callback to be installed. A guidance model that overloads `getproperty` dynamically can report the property without meaningfully supporting it, and a model that stores the orbit number under a different field name is missed entirely.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 46.
