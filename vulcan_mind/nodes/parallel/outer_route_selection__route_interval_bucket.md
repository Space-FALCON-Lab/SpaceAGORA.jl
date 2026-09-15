---
id: parallel.outer_route_selection__route_interval_bucket
label: _route_interval_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_interval_bucket
  lines:
  - 98
  - 98
inputs:
- id: interval_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `interval_s`.
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
  type: String
  units: n/a
  description: Return value of `_route_interval_bucket`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _route_interval_bucket

## Purpose
Classifies a time interval in seconds (max orbit step, control, guidance, or navigation rate) into one of five tokens for the `dt=`, `ctrl_rate=`, `guid_rate=`, and `nav_rate=` signature fields.

## Design & Implementation
`@inline _route_interval_bucket(interval_s::Float64)::String` returns "na" when `!isfinite(interval_s) || interval_s <= 0.0`, "ultra" for `<= 0.2` s, "fast" for `<= 1.0` s, "med" for `<= 10.0` s, and "slow" otherwise. It is used four times inside `outer_route_signature` but not in the mid-level or legacy signatures, which is the main reason the full signature is more specific than its fallbacks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `interval_s` | Float64 | n/a | yes | Positional argument `interval_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_interval_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:168-168`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`Inf`, `NaN`, zero, and negative values all collapse to "na", so a disabled control loop and a misconfigured negative rate are indistinguishable. The boundaries are constants; a 0.21 s and a 1.0 s loop share "fast" despite a 5x cost difference.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 98.
