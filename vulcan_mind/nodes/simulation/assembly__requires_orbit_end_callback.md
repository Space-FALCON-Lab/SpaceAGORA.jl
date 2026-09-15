---
id: simulation.assembly__requires_orbit_end_callback
label: _requires_orbit_end_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_orbit_end_callback
  lines:
  - 53
  - 53
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
  description: Return value of `_requires_orbit_end_callback`.
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

# _requires_orbit_end_callback

## Purpose
Determines whether the orbit-end detection callback should be installed, so the run can terminate on an orbit count or trigger orbit-scheduled guidance manoeuvres.

## Design & Implementation
Returns the disjunction `args.mission_configuration.mission_type == MissionOrbits || _requires_guidance_orbit_counter(args)`. The first term covers missions whose stopping condition is expressed directly in orbits; the second covers missions of any type whose guidance effectors schedule by orbit number and therefore need the counter maintained even when the mission itself ends on time or altitude. Notably this is one of only two predicates in `get_callbacks` that is not gated on `!backbone_mode`, so orbit-end detection survives the gravity-backbone split solver policy.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_orbit_end_callback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:168-168`

**Downstream**

- `callees` → [[simulation.assembly__requires_guidance_orbit_counter|_requires_guidance_orbit_counter]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:55-55`
<!-- vulcan:connections:end -->

## Limitations
Comparing `mission_type` by equality against the single `MissionOrbits` value means any future mission type that is also orbit-counted must be added here explicitly. The predicate does not consider whether the requested orbit count is reachable within the configured time span, so a mission may install the callback and never fire it. It also cannot distinguish a guidance model that needs the counter from one that merely exposes the field, inheriting that weakness from `_requires_guidance_orbit_counter`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 53.
