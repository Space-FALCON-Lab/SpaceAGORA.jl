---
id: gnc.rrt_connect_rporrtconnecttree
label: RPORRTConnectTree
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: RPORRTConnectTree
  lines:
  - 34
  - 34
inputs:
- id: nodes
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `nodes`.
- id: parents
  type: Vector{Int}
  units: n/a
  required: true
  description: Field `parents`.
- id: costs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `costs`.
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
  type: RPORRTConnectTree
  units: n/a
  description: Constructed `RPORRTConnectTree`.
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

# RPORRTConnectTree

## Purpose
Mutable-vector tree storage for RPO RRT planners: node positions, parent indices, and accumulated path cost from the root.

## Design & Implementation
An immutable struct holding three parallel vectors: `nodes::Vector{SVector{3,Float64}}` (RTN positions in metres), `parents::Vector{Int}` (index of each node's parent, `0` for the root), and `costs::Vector{Float64}` (cumulative Euclidean path length to the root). The convenience constructor `RPORRTConnectTree(root::SVector{3,Float64})` seeds the tree with `[root]`, `[0]`, `[0.0]`. Growth happens by `push!` in `rpo_rrt_extend!` and `rpo_rrt_star_add_node!`; rewiring mutates `parents` and `costs` in place.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nodes` | Vector{SVector{3, Float64}} | n/a | yes | Field `nodes`. |
| in | `parents` | Vector{Int} | n/a | yes | Field `parents`. |
| in | `costs` | Vector{Float64} | n/a | yes | Field `costs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPORRTConnectTree | n/a | — | Constructed `RPORRTConnectTree`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:533-533`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:353-353`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The three vectors are only kept consistent by convention; nothing checks equal lengths or that parent indices are less than the child index. Nearest-neighbour search delegated to `hypr_rrt_nearest_index` is linear in `length(nodes)`, so cost grows quadratically with iterations. Not thread-safe for concurrent `push!`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 34.
