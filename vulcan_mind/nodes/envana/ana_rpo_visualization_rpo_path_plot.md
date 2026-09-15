---
id: envana.ana_rpo_visualization_rpo_path_plot
label: rpo_path_plot
kind: function
source:
  file: src/analysis/visualization/rpo/rpo_visualization.jl
  symbol: rpo_path_plot
  lines:
  - 7
  - 29
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RPOVisualization namespace providing the PlotlyJS trace and layout
    constructors.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: figure
  type: PlotlyJS.Plot
  units: m
  description: Three-dimensional Plotly figure of the planned relative path with optional
    station geometry markers.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- envana
origin: agent
---
# rpo_path_plot

## Purpose
`rpo_path_plot` renders a rendezvous and proximity operations trajectory as an interactive three-dimensional Plotly figure, optionally overlaying discrete station-keeping geometry points on the same axes.

## Theory & Math
The path is supplied in the RTN frame, whose orthonormal basis is defined from the target's inertial position `r` and velocity `v` as `Rhat = r / |r|`, `Nhat = (r x v) / |r x v|`, and `That = Nhat x Rhat`. The radial axis points from the central body toward the target, the normal axis is along the orbital angular momentum, and the transverse axis completes the right-handed set. Coordinates along all three axes are relative separations in metres, so a point at the origin is co-located with the target. Because the frame rotates with the orbit, a curve that appears closed in RTN corresponds to a non-closed inertial trajectory.

## Model & Assumptions
The input `path_rtn` is converted with `Matrix{Float64}(path_rtn)` and indexed as `path[1, :]`, `path[2, :]`, and `path[3, :]`, so the argument must be three-by-N with one column per sample and the rows ordered radial, transverse, normal. The optional `station_points` argument follows the same convention and is skipped entirely when it is `nothing`. The layout sets `aspectmode="data"`, which forces equal scaling on all three axes so that a metre along the radial axis is drawn the same length as a metre along the transverse axis and relative geometry is not visually distorted.

## Design & Implementation
Traces are accumulated in a `PlotlyJS.GenericTrace[]` vector, with the path drawn as a `scatter3d` in `"lines+markers"` mode named `"planned path"` and the optional geometry drawn as markers of size 3 named `"station geometry"`. Building the vector before constructing the `Plot` keeps the optional trace a simple conditional `push!`. The title defaults to `"RPO Path"` and is accepted as an `AbstractString` keyword, so both `String` and `SubString` callers work without conversion.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | RPOVisualization namespace providing the PlotlyJS trace and layout constructors. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `figure` | PlotlyJS.Plot | m | — | Three-dimensional Plotly figure of the planned relative path with optional station geometry markers. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/visualization/rpo/rpo_visualization.jl:10-10`
<!-- vulcan:connections:end -->

## Limitations
No time information is encoded, so the figure shows geometry only and a fast approach is indistinguishable from a slow one. The function does not check that the input has three rows, so a transposed matrix produces a silently wrong plot or a bounds error. Marker size is fixed at 3 and colours are left to the Plotly defaults, giving the caller no styling control. Very long paths are pushed to the browser sample by sample with no decimation.

## Provenance
Read directly from `src/analysis/visualization/rpo/rpo_visualization.jl:7-29`, including the trace construction, the optional station-geometry branch, and the `aspectmode` layout setting.
