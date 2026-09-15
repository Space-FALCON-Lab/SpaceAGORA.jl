---
id: simulation_a.navigation_guidance_callbacks_get_navigation_callbacks
label: get_navigation_callbacks
kind: function
source:
  file: src/simulation/callbacks/navigation_guidance_callbacks.jl
  symbol: get_navigation_callbacks
  lines:
  - 1
  - 16
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: num_sats_and_config
  type: Tuple{Int, SimulationConfiguration}
  units: n/a
  required: true
  description: Spacecraft count and the run configuration whose `navigation_model`
    supplies the navigation effector vector and the matching update rates.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: navigation_callbacks
  type: Vector{DiscreteCallback}
  units: n/a
  description: One `PeriodicCallback` per navigation effector, each firing at that
    effector's configured rate and updating every spacecraft.
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
# get_navigation_callbacks

## Purpose
`get_navigation_callbacks` constructs the periodic callbacks that drive the navigation half of the guidance, navigation and control chain. Each configured navigation effector — a state estimator, a sensor model, a filter update — is given its own callback firing at its own rate, which is what allows a fast inertial propagation step and a slow absolute-fix update to coexist in one solve.

## Model & Assumptions
Navigation effectors are assumed to be rate-driven rather than event-driven, so every one becomes a `PeriodicCallback`; there is no event-triggered path here, unlike the thruster branch in `get_control_callbacks`. Each effector is assumed to update all spacecraft when fired, so the callback body loops over spacecraft indices and calls `calcNavigationEffect!` with an explicit index. The effector and rate vectors are assumed parallel; the loop indexes both with the same `eachindex` iterator.

## Design & Implementation
The result vector is preallocated with `Vector{DiscreteCallback}(undef, length(navigation_models))` and filled by index rather than pushed, giving a concretely typed return value the solver can specialise on. The per-effector closure captures the model and the spacecraft count, and the spacecraft loop is marked `@inbounds`. The sibling `get_guidance_callbacks` in the same file has an identical shape over `args.guidance_model`, calling `calcGuidanceEffect!`; keeping the two together makes the symmetry of the guidance and navigation stages obvious. `get_callbacks` appends navigation before control and control before guidance, which fixes the within-step execution order.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `num_sats_and_config` | Tuple{Int, SimulationConfiguration} | n/a | yes | Spacecraft count and the run configuration whose `navigation_model` supplies the navigation effector vector and the matching update rates. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `navigation_callbacks` | Vector{DiscreteCallback} | n/a | — | One `PeriodicCallback` per navigation effector, each firing at that effector's configured rate and updating every spacecraft. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:181-181`

**Downstream**

- `callees` → [[gncz.navigation_hooks_calcnavigationeffect_bang|calcNavigationEffect!]] · `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:10-10`
<!-- vulcan:connections:end -->

## Limitations
Unlike the control path, this loop is not threaded, so runs with many spacecraft and an expensive estimator pay the full serial cost each firing. There is no validation that the rate and effector vectors have matching lengths; a mismatch surfaces as a bounds error during construction. Both callback families are skipped entirely in gravity-backbone split mode.

## Provenance
Mapped from `src/simulation/callbacks/navigation_guidance_callbacks.jl:1-16`.
