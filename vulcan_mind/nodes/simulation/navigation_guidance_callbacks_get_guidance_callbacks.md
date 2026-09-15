---
id: simulation.navigation_guidance_callbacks_get_guidance_callbacks
label: get_guidance_callbacks
kind: function
source:
  file: src/simulation/callbacks/navigation_guidance_callbacks.jl
  symbol: get_guidance_callbacks
  lines:
  - 18
  - 18
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
  type: Vector{DiscreteCallback}
  units: n/a
  description: Return value of `get_guidance_callbacks`.
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

# get_guidance_callbacks

## Purpose
Builds one `DiscreteCallback` per configured guidance model so each fires at its own rate during integration.

## Design & Implementation
Reads `guidance_effectors` and the matching `guidance_rates` from the simulation configuration and allocates a `Vector{DiscreteCallback}` of that length. For each model it closes over the model and rate to form a callback that loops every satellite index from 1 to `num_sats`, calling `calcGuidanceEffect!` with the integrator state, parameters and time. The per-satellite loop is `@inbounds`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{DiscreteCallback} | n/a | — | Return value of `get_guidance_callbacks`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:183-183`

**Downstream**

- `callees` → [[gnc.target_energy_bracketing_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:27-27`
- `callees` → [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:27-27`
- `callees` → [[gncy.rpo_guidance_hooks_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:27-27`
<!-- vulcan:connections:end -->

## Limitations
Model and rate are paired positionally, so the two configuration vectors must be the same length and ordering; a mismatch surfaces as a bounds error rather than a validation message.

## Provenance
Mapped from `src/simulation/callbacks/navigation_guidance_callbacks.jl` line 18.
