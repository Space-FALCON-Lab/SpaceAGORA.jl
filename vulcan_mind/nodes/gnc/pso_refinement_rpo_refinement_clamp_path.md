---
id: gnc.pso_refinement_rpo_refinement_clamp_path
label: rpo_refinement_clamp_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_clamp_path
  lines:
  - 17
  - 17
inputs:
- id: path
  type: Any
  units: n/a
  required: true
  description: Positional argument `path`.
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
  type: Any
  units: n/a
  description: Return value of `rpo_refinement_clamp_path`. Returns `pts`.
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

# rpo_refinement_clamp_path

## Purpose
Keeps refined control points inside the same search box the swarm was confined to, so refinement cannot escape the region the planner was configured to trust.

## Design & Implementation
Converts the path to a dense matrix, obtains the per-axis lower and upper bounds from `rpo_pso_bounds` given the start and goal columns, and clamps every interior column axis by axis. The first and last columns are left untouched because they are the fixed endpoints.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | Any | n/a | yes | Positional argument `path`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_clamp_path`. Returns `pts`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints|rpo_fit_bezier_fixed_endpoints]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:148-148`
- [[gnc.pso_refinement_rpo_try_accept_refinement|rpo_try_accept_refinement]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:165-165`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- `callees` → [[gncy.pso_helpers_rpo_pso_bounds|rpo_pso_bounds]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:19-19`
<!-- vulcan:connections:end -->

## Limitations
Clamping moves a control point to the box face rather than rejecting it, which can create a kink and a locally worse path that the subsequent cost test then has to catch.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 17.
