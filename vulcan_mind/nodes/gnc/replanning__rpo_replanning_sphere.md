---
id: gnc.replanning__rpo_replanning_sphere
label: _rpo_replanning_sphere
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: _rpo_replanning_sphere
  lines:
  - 51
  - 51
inputs:
- id: value
  type: Any
  units: n/a
  required: true
  description: Positional argument `value`.
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
  type: RPOReplanningSphere
  units: n/a
  description: Return value of `_rpo_replanning_sphere`. Returns `RPOReplanningSphere(`.
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

# _rpo_replanning_sphere

## Purpose
Normalises heterogeneous obstacle inputs (an existing `RPOReplanningSphere`, a `NamedTuple`, or any struct with recognisably named fields) into a validated `RPOReplanningSphere`, so `RPOReplanningConfig` can be built from scenario files and test fixtures without a fixed schema.

## Design & Implementation
If `value isa RPOReplanningSphere` it is returned as-is. Otherwise the centre is read via `_rpo_replanning_property` from `:center_rtn` or `:center`, and the radius from `:radius_m` or `:radius`; a missing centre or radius throws `ArgumentError` with a message naming both accepted keys. Optional fields use the same aliasing: `:appear_time_s`/`:appear_time` (default 0.0), `:disappear_time_s`/`:disappear_time` (default `Inf`), `:velocity_rtn_mps`/`:velocity` (default `zeros(3)`) and `:label` (default `"dynamic_sphere"`, wrapped in `String`). The values are passed to the validating `RPOReplanningSphere(center, radius; ...)` constructor, which enforces positive radius and ordered lifetime.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value` | Any | n/a | yes | Positional argument `value`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOReplanningSphere | n/a | — | Return value of `_rpo_replanning_sphere`. Returns `RPOReplanningSphere(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:71-71`
- `callees` → [[gnc.replanning__rpo_replanning_property|_rpo_replanning_property]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:53-53`
- `callees` → [[gnc.replanning_rporeplanningconfig|RPOReplanningConfig]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:68-68`
- `callees` → [[gnc.replanning_rporeplanningsphere|RPOReplanningSphere]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:57-57`
<!-- vulcan:connections:end -->

## Limitations
A centre supplied as a 2-vector or a 4-vector fails inside `SVector{3,Float64}(center)` with a generic `DimensionMismatch` rather than a message mentioning the obstacle. Radius supplied as a string is not coerced and raises `MethodError` from `Float64(radius_m)`. Unknown extra properties are ignored silently, so a misspelt `velocity_rtn_mps` key quietly yields a stationary sphere.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 51.
