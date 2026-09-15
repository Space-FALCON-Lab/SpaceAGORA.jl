---
id: gnc.path_retiming_rpo_arc_length_params
label: rpo_arc_length_params
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_retiming.jl
  symbol: rpo_arc_length_params
  lines:
  - 2
  - 2
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
  description: Return value of `rpo_arc_length_params`. Returns `s`.
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

# rpo_arc_length_params

## Purpose
Computes the cumulative arc-length coordinate `s` (metres) of every column in a sampled RPO path so that later retiming stages (`rpo_interpolate_along_path`, `rpo_curvature_from_samples`, `rpo_retime_path`) can parameterise the path by distance travelled rather than by sample index.

## Theory & Math
Cumulative chord length: $s_1 = 0,\quad s_j = s_{j-1} + \lVert \mathbf{p}_j - \mathbf{p}_{j-1} \rVert_2$ for $j = 2..n$, where $\mathbf{p}_j \in \mathbb{R}^3$ is the $j$-th column of `pts` in metres.

## Design & Implementation
`points` is converted to a `Matrix{Float64}` (one 3-vector per column) and `s` is allocated as `zeros(size(pts, 2))`. A single `@inbounds` loop over `j in 2:n` accumulates `s[j] = s[j-1] + norm(pts[:, j] - pts[:, j-1])`, i.e. the Euclidean chord length between consecutive samples. The first entry is always `0.0`. The function is pure; it does not mutate its argument and allocates a temporary vector per iteration for the column difference.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_arc_length_params`. Returns `s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_resample_polyline_points|rpo_resample_polyline_points]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:48-48`
- [[gnc.pso_refinement_rpo_refinement_sample_params|rpo_refinement_sample_params]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:107-107`
- [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:133-133`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Chord length underestimates true arc length on curved segments; accuracy depends entirely on how densely `rpo_sample_path` sampled the path. Repeated identical samples yield zero-length segments, which downstream code must guard against (they do so via `eps(Float64)` checks). An empty input matrix returns an empty vector without error; a single-column input returns `[0.0]`. The conversion `Matrix{Float64}(points)` copies the data each call.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_retiming.jl` line 2.
