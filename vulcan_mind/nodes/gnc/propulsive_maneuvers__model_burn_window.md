---
id: gnc.propulsive_maneuvers__model_burn_window
label: _model_burn_window
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _model_burn_window
  lines:
  - 142
  - 142
inputs:
- id: controlModel
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `controlModel`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  description: Return value of `_model_burn_window`. Returns `NaN, NaN` or `controlModel.start_burn_time[i],
    controlModel.stop_burn_time[i]`.
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

# _model_burn_window

## Purpose
Reads the raw start and stop burn times stored on the thruster model for spacecraft `i`.

## Design & Implementation
Bounds-checks `i` against `length(controlModel.start_burn_time)` and returns the tuple `(NaN, NaN)` when out of range, otherwise `(controlModel.start_burn_time[i], controlModel.stop_burn_time[i])` in seconds of simulation time. The sentinel `-1.0` written by `calcControlEffect!` on schedule clear is returned as-is, and is rejected downstream by the `stop_time > start_time` test rather than here.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_model_burn_window`. Returns `NaN, NaN` or `controlModel.start_burn_time[i], controlModel.stop_burn_time[i]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__effective_burn_window|_effective_burn_window]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:150-150`
- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:446-446`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only `start_burn_time` is bounds-checked; a `stop_burn_time` array shorter than `start_burn_time` produces an out-of-bounds error rather than the `NaN` sentinel. No ordering or finiteness validation happens at this level, so an inverted or cleared window is returned unflagged.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 142.
