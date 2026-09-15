---
id: gnc.hypr_utils_hypr_rrt_near_indices
label: hypr_rrt_near_indices
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_rrt_near_indices
  lines:
  - 156
  - 156
inputs:
- id: tree
  type: Any
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
- id: radius
  type: Any
  units: n/a
  required: true
  description: Positional argument `radius`.
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
  description: Return value of `hypr_rrt_near_indices`. Returns `idxs`.
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

# hypr_rrt_near_indices

## Purpose
Collects every tree node within a radius of a query state, used by rewiring to find candidate parents.

## Design & Implementation
Squares the radius once, then scans `tree.nodes` pushing each index whose squared distance is at or below it onto an `Int` vector. Returns the vector, which may be empty.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | Any | n/a | yes | Positional argument `tree`. |
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `radius` | Any | n/a | yes | Positional argument `radius`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_rrt_near_indices`. Returns `idxs`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_near_indices|rpo_rrt_near_indices]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:50-50`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:157-157`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:160-160`
<!-- vulcan:connections:end -->

## Limitations
Linear scan and a growing vector; the result is unordered, so a caller wanting the nearest first must sort it.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 156.
