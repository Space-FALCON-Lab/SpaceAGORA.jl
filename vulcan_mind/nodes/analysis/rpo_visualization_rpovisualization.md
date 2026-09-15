---
id: analysis.rpo_visualization_rpovisualization
label: RPOVisualization
kind: module
source:
  file: src/analysis/visualization/rpo/rpo_visualization.jl
  symbol: RPOVisualization
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
- analysis
charts:
- analysis
origin: agent
---

# RPOVisualization

## Purpose

`RPOVisualization` is the plotting module for rendezvous and proximity operations results. It wraps PlotlyJS and exports exactly two figure builders: `rpo_path_plot`, which draws a planned relative trajectory in the RTN frame together with optional station-keeping geometry points, and `rpo_tracking_plot`, which draws the tracking-error norm against time.

## Design & Implementation

The module does `using PlotlyJS` and builds figures purely functionally: each function constructs traces and returns a `Plot` object rather than displaying or writing a file, which leaves saving and layout composition to the caller. `rpo_path_plot` converts `path_rtn` to a `Matrix{Float64}` and reads rows 1, 2 and 3 as the x, y and z channels of a `scatter3d` in `lines+markers` mode, appending a second marker-only trace named `station geometry` when `station_points` is not `nothing`, and sets `scene=attr(aspectmode="data")` so relative distances are not visually distorted.

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

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/visualization/rpo/rpo_visualization.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Both entry points assume 3xN column-major input; an Nx3 array silently plots the wrong axes rather than erroring. Coordinates are labelled only through the tracking plot's axis titles (seconds and metres), so no unit conversion or validation is performed. Loading the module requires PlotlyJS to be installed, and the returned `Plot` objects need a display backend before anything is visible.

## Provenance
Mapped from `src/analysis/visualization/rpo/rpo_visualization.jl` line 1.
