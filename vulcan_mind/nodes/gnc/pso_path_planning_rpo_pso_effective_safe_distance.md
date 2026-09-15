---
id: gnc.pso_path_planning_rpo_pso_effective_safe_distance
label: rpo_pso_effective_safe_distance
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_effective_safe_distance
  lines:
  - 104
  - 104
inputs:
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: safe_distance_m
  type: Any
  units: n/a
  required: true
  description: Positional argument `safe_distance_m`.
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
  description: 'Return value of `rpo_pso_effective_safe_distance`. Returns `safe >
    0.0 || cfg.safe_distance_m <= 0.0 ? safe : cfg.safe_distance_m`.'
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

# rpo_pso_effective_safe_distance

## Purpose
Resolves which safety distance the planner actually uses when the caller may or may not have supplied one.

## Design & Implementation
Returns the config's `safe_distance_m` if the argument is `nothing`. Otherwise it converts the argument to `Float64` and returns it if it is positive or if the config value is non-positive; when the caller passed zero and the config has a positive default, the config wins. This lets a zero from a generic caller mean use the default rather than disable the margin.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Any | n/a | yes | Positional argument `safe_distance_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_effective_safe_distance`. Returns `safe > 0.0 \|\| cfg.safe_distance_m <= 0.0 ? safe : cfg.safe_distance_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:191-191`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:106-106`
<!-- vulcan:connections:end -->

## Limitations
A caller who genuinely wants a zero margin cannot get one when the config default is positive; they must rebuild the config.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 104.
