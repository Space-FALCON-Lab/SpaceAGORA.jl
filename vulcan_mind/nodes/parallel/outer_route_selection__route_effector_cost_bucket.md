---
id: parallel.outer_route_selection__route_effector_cost_bucket
label: _route_effector_cost_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_effector_cost_bucket
  lines:
  - 124
  - 124
inputs:
- id: raw
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `raw`.
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
  description: Return value of `_route_effector_cost_bucket`.
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

# _route_effector_cost_bucket

## Purpose
Sanitises the effector cost-class label ("light", "medium", "heavy") for the `eff_cost=` signature field, tolerating case and whitespace differences and refusing to let pipe characters break the signature format.

## Design & Implementation
`@inline _route_effector_cost_bucket(raw::AbstractString)::String` lowercases and strips `raw`; returns it unchanged if it is one of "light", "medium", or "heavy"; returns "unknown" for an empty string; otherwise returns the token with `|` replaced by `_`. Used by `outer_route_signature` and the mid-level signature in `_outer_route_signature_hierarchy`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | AbstractString | n/a | yes | Positional argument `raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_effector_cost_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:197-197`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:177-177`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Any non-canonical, non-empty label is accepted verbatim, so a misspelled class becomes its own bucket rather than an error. The classification itself (how a workload earns "heavy") is computed elsewhere in the feature extractor and is not validated here.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 124.
