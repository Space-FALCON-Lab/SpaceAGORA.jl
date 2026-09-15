---
id: gnc.pso_parameters__rpo_pso_normalize_kwargs
label: _rpo_pso_normalize_kwargs
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: _rpo_pso_normalize_kwargs
  lines:
  - 437
  - 437
inputs:
- id: kwargs
  type: Any
  units: n/a
  required: true
  description: Positional argument `kwargs`.
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
  description: Return value of `_rpo_pso_normalize_kwargs`. Returns `(; pairs...)`.
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

# _rpo_pso_normalize_kwargs

## Purpose
Translates legacy `pso_*` keyword names into the canonical `RPOPSOConfig` field names using the `RPO_PSO_CONFIG_ALIASES` table, so callers written against older keyword spellings keep working.

## Design & Implementation
Iterates `(key, value)` over the supplied `kwargs`, pushes `get(RPO_PSO_CONFIG_ALIASES, key, key) => value` onto a `Vector{Pair{Symbol, Any}}`, and returns `(; pairs...)` as a NamedTuple. The alias dictionary maps roughly 95 names such as `:pso_waypoints => :n_waypoints`, `:pso_safe_distance => :safe_distance_m`, and `:pso_rrt_warmstart_box_margin => :rrt_warmstart_box_margin_m`; unknown keys pass through unchanged.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kwargs` | Any | n/a | yes | Positional argument `kwargs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rpo_pso_normalize_kwargs`. Returns `(; pairs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:584-584`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
If both an alias and its canonical name are supplied, the later pair wins in the NamedTuple construction without warning. Unrecognised keys are not rejected here; they surface as a `MethodError` from `RPOPSOConfig(; values...)` in `rpo_pso_config`. `Pair{Symbol, Any}` storage loses type information, but the tuple is consumed immediately.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 437.
