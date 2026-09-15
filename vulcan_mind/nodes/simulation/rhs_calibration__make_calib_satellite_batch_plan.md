---
id: simulation.rhs_calibration__make_calib_satellite_batch_plan
label: _make_calib_satellite_batch_plan
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _make_calib_satellite_batch_plan
  lines:
  - 196
  - 196
inputs:
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
  description: Return value of `_make_calib_satellite_batch_plan`. Returns `(`.
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

# _make_calib_satellite_batch_plan

## Purpose
Constructs the baseline satellite-batch execution plan that every calibration sweep starts from and that is restored when a persisted entry records `mode = "satellite_batch"`.

## Design & Implementation
Zero-argument `@inline` constructor returning a NamedTuple with `mode = :satellite_batch`, `allotment = 1`, `scheduler = :static`, `dominant_axis = :satellite`, `policy_applied = true`, and `effector_decision = _CALIB_SERIAL_EFFECTOR_DECISION`. Because it allocates nothing beyond the tuple it is called unconditionally at the head of `_rhs_plan_candidates` and again on every cache hit in `_rhs_calib_lookup`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_make_calib_satellite_batch_plan`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__rhs_calib_lookup|_rhs_calib_lookup]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:168-168`
- [[simulation.rhs_calibration__rhs_plan_candidates|_rhs_plan_candidates]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:227-227`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The static scheduler and single allotment are fixed; there is no way to calibrate a satellite-batch plan with a different chunk size. The tuple field set is duplicated between this function and `_make_calib_flat_plan`, so any new field required by `_rhs_execution_plan` must be added in both places or the override will fail with a `field not found` error at RHS time.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 196.
