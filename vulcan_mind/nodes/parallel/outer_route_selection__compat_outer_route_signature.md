---
id: parallel.outer_route_selection__compat_outer_route_signature
label: _compat_outer_route_signature
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _compat_outer_route_signature
  lines:
  - 134
  - 134
inputs:
- id: f
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Positional argument `f`.
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
  type: String
  units: n/a
  description: Return value of `_compat_outer_route_signature`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _compat_outer_route_signature

## Purpose
Builds the legacy nine-field routing signature that older feedback histories were recorded under, so a fresh full signature with no history can still fall back to statistics gathered before the signature schema was extended.

## Design & Implementation
`@inline _compat_outer_route_signature(f::OuterRouteFeatures)::String` joins with `|` the fields `cat=`, `sat=` (`_route_sat_bucket`), `links=` (`_route_link_bucket`), `mission=` (`_route_mission_bucket`), `nbody=`, `srp=`, `harm=` (`_route_harmonics_bucket`), `ctrl=`, and `orient=`, with booleans rendered as "1"/"0". It is the third and coarsest entry returned by `_outer_route_signature_hierarchy`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_compat_outer_route_signature`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:199-199`

**Downstream**

- `callees` → [[parallel.outer_route_selection__route_harmonics_bucket|_route_harmonics_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:142-142`
- `callees` → [[parallel.outer_route_selection__route_link_bucket|_route_link_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:138-138`
- `callees` → [[parallel.outer_route_selection__route_mission_bucket|_route_mission_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:139-139`
- `callees` → [[parallel.outer_route_selection__route_sat_bucket|_route_sat_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:137-137`
- `callees` → [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:149-149`
<!-- vulcan:connections:end -->

## Limitations
Omits density family, solver, rates, GRAM flags, thermal, and effector counts, so history looked up through this key can come from workloads with very different atmosphere or effector costs. `f.category` is interpolated raw without pipe sanitisation.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 134.
