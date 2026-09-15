---
id: gnc.bridge_helpers__with_time_switch
label: _with_time_switch
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _with_time_switch
  lines:
  - 205
  - 205
inputs:
- id: context
  type: NamedTuple
  units: n/a
  required: true
  description: Positional argument `context`.
- id: time_switch
  type: Any
  units: n/a
  required: true
  description: Positional argument `time_switch`.
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
  type: Tuple
  units: n/a
  description: Return value of `_with_time_switch`. Returns `(; context..., time_switch)`.
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

# _with_time_switch

## Purpose
Returns a copy of an aerobraking runtime `context::NamedTuple` extended with a `time_switch` field, used by the trajectory predictor's `asim_ctrl_targeting` to carry the control switching time into downstream simulation calls.

## Design & Implementation
Implemented as `(; context..., time_switch)`, which splats every existing field and appends `time_switch`; if the context already contains `time_switch` the new value wins because later keys override earlier ones. The function is `@inline` and wrapped in an `isdefined` guard so repeated inclusion of `bridge_helpers.jl` does not redefine it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `context` | NamedTuple | n/a | yes | Positional argument `context`. |
| in | `time_switch` | Any | n/a | yes | Positional argument `time_switch`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_with_time_switch`. Returns `(; context..., time_switch)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_predictor_out_drag_passage_affect_bang|out_drag_passage_affect!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:400-400`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:400-400`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The full context is copied on each call, allocating proportionally to the number of fields. The units of `time_switch` (seconds since `time_0`) are not checked, and a `nothing` or negative value passes through unchanged. Type instability arises if callers pass `time_switch` as `Int` in some paths and `Float64` in others.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 205.
