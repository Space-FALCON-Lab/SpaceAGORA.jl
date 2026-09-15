---
id: gnc.hypr_utils_hyprutils
label: HYPRUtils
kind: module
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: HYPRUtils
  lines:
  - 2
  - 2
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

# HYPRUtils

## Purpose
The shared HYPR toolkit — path sampling, PSO scheduling and bookkeeping, and RRT tree operations — that both the RPO planner and the robot-arm planner build on rather than each carrying its own copy.

## Design & Implementation
Depends only on LinearAlgebra and exports three groups: curve utilities (`hypr_path_length`, the two Bezier evaluators, `hypr_sample_count_path`), swarm utilities (`hypr_iteration_weights`, `hypr_material_improvement`, `hypr_protected_particle_mask`) and tree utilities (`hypr_rrt_nearest_index`, `hypr_rrt_near_indices`, `hypr_rrt_steer`, `hypr_rrt_tree_path`, `hypr_rrt_join_paths`, `hypr_rrt_refresh_subtree_costs!`). Two private accessors tolerate the RPO and robot-arm trees naming their parent and cost vectors differently, so the tree functions are duck-typed over either.

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

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Everything operates on dense `Matrix{Float64}` columns or vectors of state vectors with no spatial index, so nearest-neighbour and near-radius queries are linear scans and the RRT routines degrade quadratically with tree size.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 2.
