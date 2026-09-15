---
id: gnc.pso_helpers_rpo_pso_warmstart_bounds
label: rpo_pso_warmstart_bounds
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_helpers.jl
  symbol: rpo_pso_warmstart_bounds
  lines:
  - 14
  - 14
inputs:
- id: warmstart_path
  type: Any
  units: n/a
  required: true
  description: Positional argument `warmstart_path`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: SVector
  units: n/a
  description: Return value of `rpo_pso_warmstart_bounds`. Returns `SVector{3, Float64}(lo),
    SVector{3, Float64}(hi)`.
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

# rpo_pso_warmstart_bounds

## Purpose
Derives an axis-aligned particle swarm search box that tightly encloses a rapidly-exploring random tree warm-start path, so the swarm refines near an already feasible route instead of searching the whole relative volume.

## Design & Implementation
Converts `warmstart_path` to a `Matrix{Float64}` and validates it: `size(pts, 1) == 3` or it throws `ArgumentError("warmstart_path must have three rows.")`, and `size(pts, 2) >= 2` or it throws `ArgumentError("warmstart_path must contain at least start and goal.")`. It then takes the per-row minimum and maximum across columns, pads both by `cfg.rrt_warmstart_box_margin_m` metres, and returns the pair as `SVector{3, Float64}` lower and upper bounds in the RTN frame.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `warmstart_path` | Any | n/a | yes | Positional argument `warmstart_path`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `rpo_pso_warmstart_bounds`. Returns `SVector{3, Float64}(lo), SVector{3, Float64}(hi)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_reset_swarm_bang|reset_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:326-326`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:326-326`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The box is axis-aligned in RTN, so a long diagonal path yields a volume far larger than the corridor around it, diluting the swarm. The margin is a single isotropic scalar with no per-axis control and no lower bound, so a degenerate warm start that collapses onto one plane gives a box of zero thickness plus margin in that axis. Converting to a dense matrix copies the whole path.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_helpers.jl` line 14.
