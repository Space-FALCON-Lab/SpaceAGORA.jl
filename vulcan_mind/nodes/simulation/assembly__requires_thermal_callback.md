---
id: simulation.assembly__requires_thermal_callback
label: _requires_thermal_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_thermal_callback
  lines:
  - 83
  - 83
inputs:
- id: effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `effectors`.
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
  description: Return value of `_requires_thermal_callback`.
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

# _requires_thermal_callback

## Purpose
Decides whether the per-step thermal callback must be installed to integrate heat loads on spacecraft surface links.

## Design & Implementation
Short-circuits to `false` unless `_requires_density_callback(effectors, args)` holds, since aerothermal heating needs density. It then scans `args.dynamics_model.spacecraft` under `@inbounds` and returns `true` at the first spacecraft whose `links` collection is non-empty, otherwise `false`. A link is the thermal-node abstraction that accumulates heat load, so a vehicle with no links has nothing to integrate. This predicate is also one of the four reasons `_requires_staged_density_callback` returns `true`, because the thermal callback reads `shared_buffers.densities` each step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_thermal_callback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.assembly__requires_staged_density_callback|_requires_staged_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:38-38`
- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:164-164`

**Downstream**

- `callees` → [[simulation.assembly__requires_density_callback|_requires_density_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:84-84`
<!-- vulcan:connections:end -->

## Limitations
One spacecraft with links forces the thermal callback on for the whole constellation, so in a mixed fleet the callback runs and the density staging cost is paid even for vehicles that have none. The check is on emptiness only: links configured with a zero area or a disabled heating model still count. The linear scan over all spacecraft repeats on every call, and nothing validates that the links present are thermally complete, so a partially configured link surfaces as an error inside the callback rather than at assembly time.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 83.
