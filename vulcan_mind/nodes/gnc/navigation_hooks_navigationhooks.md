---
id: gnc.navigation_hooks_navigationhooks
label: NavigationHooks
kind: module
source:
  file: src/gnc/navigation/navigation_hooks.jl
  symbol: NavigationHooks
  lines:
  - 1
  - 1
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# NavigationHooks

## Purpose
`NavigationHooks` is the extension point and namespace for navigation and sensor-estimator models driven by the simulation's typed periodic callbacks. It defines the generic `calcNavigationEffect!` entry point that every navigation effector specialises, and it aggregates the RPO relative-navigation geometry and distance-query code into one exported surface.

## Design & Implementation
The module imports `LinearAlgebra` and `StaticArrays`, then declares `calcNavigationEffect!(model, u, p, t::Float64, sat_idx::Int)` whose only body is `throw(MethodError(calcNavigationEffect!, (model, u, p, t, sat_idx)))`. That fallback exists so an unimplemented navigation model fails loudly at the callback rather than silently doing nothing. Seven `include` calls pull in the RPO reference geometry (`station_geometry.jl`, `cubesat_geometry.jl`, `rpo_reference_geometry.jl`) and the distance layer (`mesh_distance.jl`, `clearance.jl`, `surface_frames.jl`, `rpo_distance_queries.jl`), all resolved relative to `@__DIR__`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/navigation/navigation_hooks.jl`

**Downstream**

- `callees` → [[gncz.navigation_hooks_calcnavigationeffect_bang|calcNavigationEffect!]] · `callers` · call · `src/gnc/navigation/navigation_hooks.jl:12-12`
<!-- vulcan:connections:end -->

## Limitations
The catch-all method matches any argument types, so a specialisation whose signature differs only by, say, an `Int32` satellite index or a non-`Float64` time will be shadowed by the fallback and raise `MethodError` at run time instead of being caught at load time. Because everything is textually included into one module, the RPO geometry types are loaded even for simulations that never perform rendezvous, and the include order is a hard, undocumented dependency: the distance files assume the geometry types already exist.

## Provenance
Mapped from `src/gnc/navigation/navigation_hooks.jl` line 1.
