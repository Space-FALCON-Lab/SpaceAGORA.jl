---
id: simulation.rhs_calibration__make_calib_flat_plan
label: _make_calib_flat_plan
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _make_calib_flat_plan
  lines:
  - 207
  - 207
inputs:
- id: allotment
  type: Int
  units: n/a
  required: true
  description: Positional argument `allotment`.
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
  description: Return value of `_make_calib_flat_plan`. Returns `(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _make_calib_flat_plan

## Purpose
Constructs the NamedTuple execution plan describing a flat constellation-by-effector work queue with a given worker allotment, used both as a sweep candidate and when rehydrating a persisted calibration.

## Design & Implementation
`@inline` constructor taking `allotment::Int`. Returns a NamedTuple with `mode = :flat_constellation_effector_queue`, `allotment = max(1, allotment)`, `scheduler = :dynamic`, `dominant_axis = :flat_effector`, `policy_applied = true`, and `effector_decision = _CALIB_SERIAL_EFFECTOR_DECISION` (a constant tuple with `use_threads = false`, `allotment = 1`, `mode = :off`, `policy_applied = false`). The shape must match what `_rhs_execution_plan` in `setup.jl` expects to read from `SharedBuffers.rhs_plan_override[]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_make_calib_flat_plan`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__rhs_calib_lookup|_rhs_calib_lookup]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:169-169`
- [[simulation.rhs_calibration__rhs_plan_candidates|_rhs_plan_candidates]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:243-243`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The plan always pins the per-effector inner decision to serial; nested effector threading is never a calibrated candidate. An allotment below 1 is silently clamped rather than rejected, and nothing checks that `allotment` does not exceed the live thread budget (callers such as `_rhs_plan_candidates` do that clamp themselves).

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 207.
