---
id: gnc.planner_comparison_rpo_comparison_station_trace
label: rpo_comparison_station_trace
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_station_trace
  lines:
  - 845
  - 845
inputs:
- id: batch
  type: Any
  units: n/a
  required: true
  description: Positional argument `batch`.
- id: max_points
  type: Integer
  units: n/a
  required: false
  description: Keyword argument `max_points` (default `2_000`).
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
  type: PlotlyJS.scatter3d
  units: n/a
  description: Return value of `rpo_comparison_station_trace`. Returns `rpo_comparison_station_mesh_trace(batch.station_triangles)`
    or `PlotlyJS.scatter3d(`.
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

# rpo_comparison_station_trace

## Purpose
Selects the station background trace for 3-D plots: a triangle mesh when `batch.station_triangles` is present, otherwise a stride-downsampled scatter of the station point cloud limited to roughly `max_points` markers.

## Design & Implementation
`rpo_comparison_station_trace(batch; max_points::Integer = 2_000)` returns `rpo_comparison_station_mesh_trace(batch.station_triangles)` when that property exists and is not `nothing`. Otherwise reads `pts = batch.geometry.station.points_body`, computes `stride = max(1, ceil(size(pts, 2) / max(1, max_points)))`, samples `1:stride:end`, and returns a `scatter3d` of size-2 grey markers at 30 % opacity with hover disabled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `max_points` | Integer | n/a | no | Keyword argument `max_points` (default `2_000`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlotlyJS.scatter3d | n/a | — | Return value of `rpo_comparison_station_trace`. Returns `rpo_comparison_station_mesh_trace(batch.station_triangles)` or `PlotlyJS.scatter3d(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_failed_paths_plot|rpo_comparison_failed_paths_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:993-993`
- [[gnc.planner_comparison_rpo_comparison_path_family_plot|rpo_comparison_path_family_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:865-865`
- [[gnc.planner_comparison_rpo_comparison_single_path_plot|rpo_comparison_single_path_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:930-930`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpo_comparison_station_mesh_trace|rpo_comparison_station_mesh_trace]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:847-847`
<!-- vulcan:connections:end -->

## Limitations
Uniform striding preserves spatial coverage only if the point cloud is stored in a spatially uncorrelated order; a cloud sorted along one axis would be thinned unevenly. The point cloud is in the body frame (`points_body`) while paths are plotted in RTN, which is correct only because the comparison assumes the station body frame is aligned with RTN.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 845.
