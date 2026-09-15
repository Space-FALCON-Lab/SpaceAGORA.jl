---
id: gnc.rrt_connect_rpo_rrt_shortcut_path
label: rpo_rrt_shortcut_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_shortcut_path
  lines:
  - 219
  - 219
inputs:
- id: path
  type: Any
  units: n/a
  required: true
  description: Positional argument `path`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: settings
  type: RPORRTConnectSettings
  units: n/a
  required: true
  description: Positional argument `settings`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: true
  description: Keyword argument `safe_distance_m`.
- id: rng
  type: Any
  units: n/a
  required: false
  description: Keyword argument `rng` (default `Random.default_rng()`).
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
  description: Return value of `rpo_rrt_shortcut_path`. Returns `pts`.
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

# rpo_rrt_shortcut_path

## Purpose
Randomly removes intermediate waypoints from a polyline path when the direct replacement segment is collision-free and does not lengthen the path.

## Design & Implementation
Converts `path` to `Matrix{Float64}` (3 x n). If `n <= 2` it returns immediately. For up to `settings.shortcut_iters` iterations it draws `i in 1:(n-2)` and `j in (i+2):n` from `rng`, tests `rpo_rrt_segment_is_safe(pts[:,i], pts[:,j], geometry, settings; safe_distance_m)`, and on success builds `candidate = pts[:, vcat(1:i, j:n)]`, accepting it only if `rpo_path_length(candidate) <= rpo_path_length(pts) + 1e-9`. Because removing points from a polyline never lengthens it, that guard is a safety net. Returns the reduced matrix.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | Any | n/a | yes | Positional argument `path`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `settings` | RPORRTConnectSettings | n/a | yes | Positional argument `settings`. |
| in | `safe_distance_m` | Real | n/a | yes | Keyword argument `safe_distance_m`. |
| in | `rng` | Any | n/a | no | Keyword argument `rng` (default `Random.default_rng()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_shortcut_path`. Returns `pts`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:596-596`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:398-398`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_path_length|rpo_path_length]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:236-236`
- `callees` → [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:227-227`
<!-- vulcan:connections:end -->

## Limitations
Each accepted shortcut allocates a new matrix; each iteration allocates the `keep` index vector. Random pair selection can repeatedly test the same pair, so 80 iterations do not guarantee convergence to a locally minimal path. The result depends on `rng` state, making planner output non-deterministic unless the caller seeds it.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 219.
