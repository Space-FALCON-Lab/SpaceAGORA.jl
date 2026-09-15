---
id: gnc.planner_comparison_rpo_comparison_station_mesh_trace
label: rpo_comparison_station_mesh_trace
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_station_mesh_trace
  lines:
  - 813
  - 813
inputs:
- id: triangles
  type: Any
  units: n/a
  required: true
  description: Positional argument `triangles`.
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
  type: PlotlyJS.mesh3d
  units: n/a
  description: Return value of `rpo_comparison_station_mesh_trace`. Returns `PlotlyJS.mesh3d(`.
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

# rpo_comparison_station_mesh_trace

## Purpose
Converts a triangle soup describing the station keep-out geometry into a Plotly `mesh3d` trace, so 3-D path plots show the actual station surface rather than a point cloud when triangle data is available.

## Design & Implementation
`rpo_comparison_station_mesh_trace(triangles)` converts input to `Matrix{Float64}` with 3 rows and `3 * ntri` columns (consecutive vertex triples), copies the vertices into a fresh array, and builds zero-based index vectors `i = src`, `j = src + 1`, `k = src + 2` for each triangle as Plotly requires. Returns `PlotlyJS.mesh3d` with grey colour `rgb(150,160,170)`, `opacity = 0.55`, `flatshading = true`, name "Gateway mesh", and hover disabled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `triangles` | Any | n/a | yes | Positional argument `triangles`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlotlyJS.mesh3d | n/a | — | Return value of `rpo_comparison_station_mesh_trace`. Returns `PlotlyJS.mesh3d(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_station_trace|rpo_comparison_station_trace]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:847-847`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Vertices are not de-duplicated, so shared vertices are emitted once per triangle, tripling the payload for closed meshes. A column count not divisible by 3 silently drops the trailing vertices. The name "Gateway mesh" is hard-coded regardless of the station modelled.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 813.
