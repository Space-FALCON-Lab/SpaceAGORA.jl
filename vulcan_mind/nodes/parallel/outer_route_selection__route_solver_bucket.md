---
id: parallel.outer_route_selection__route_solver_bucket
label: _route_solver_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_solver_bucket
  lines:
  - 80
  - 80
inputs:
- id: mode
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `mode`.
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
  description: Return value of `_route_solver_bucket`.
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

# _route_solver_bucket

## Purpose
Maps the ODE solver mode string to a compact canonical token for the `solver=` field of routing signatures, so solver aliases and split-IMEX variants share routing history.

## Design & Implementation
`@inline _route_solver_bucket(mode::AbstractString)::String` lowercases and strips, then returns "auto" for ("auto", "auto_stiff", "autostiff") or an empty string, "rodas" for ("rodas5p", "rodas"), "tsit5" for "tsit5", "split" for any token starting with "split_imex", and "mrate" for "multirate". Unrecognised modes are returned with `|` replaced by `_`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mode` | AbstractString | n/a | yes | Positional argument `mode`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_solver_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:195-195`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:167-167`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
An empty mode silently means "auto", hiding a missing configuration. All `split_imex*` variants collapse to one bucket even though their sub-modes can differ in cost. New solver names fall through unbucketed, fragmenting history until an alias is added here.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 80.
