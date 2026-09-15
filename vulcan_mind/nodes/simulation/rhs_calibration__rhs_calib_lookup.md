---
id: simulation.rhs_calibration__rhs_calib_lookup
label: _rhs_calib_lookup
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calib_lookup
  lines:
  - 160
  - 160
inputs:
- id: sig
  type: String
  units: n/a
  required: true
  description: Positional argument `sig`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_rhs_calib_lookup`.
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

# _rhs_calib_lookup

## Purpose
Resolves a calibration signature to a ready-to-use execution plan NamedTuple, or `nothing` when no persisted or in-process result exists.

## Design & Implementation
Calls `_rhs_calib_load!()` first, then reads `_rhs_calib_cache[sig]` under `_rhs_calib_lock`. With an entry present it extracts `mode` and `allotment = max(1, Int(...))` and dispatches on the mode string: `"satellite_batch"` returns `_make_calib_satellite_batch_plan()`, `"flat_constellation_effector_queue"` returns `_make_calib_flat_plan(allotment)`, and any other string returns `nothing`. Return type is `Union{Nothing, NamedTuple}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sig` | String | n/a | yes | Positional argument `sig`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_rhs_calib_lookup`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:327-327`

**Downstream**

- `callees` → [[simulation.rhs_calibration__make_calib_flat_plan|_make_calib_flat_plan]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:169-169`
- `callees` → [[simulation.rhs_calibration__make_calib_satellite_batch_plan|_make_calib_satellite_batch_plan]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:168-168`
- `callees` → [[simulation.rhs_calibration__rhs_calib_load_bang|_rhs_calib_load!]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:161-161`
<!-- vulcan:connections:end -->

## Limitations
A stale flat allotment is honoured without checking the current `effective_inner_thread_budget()`, so a file written on the same machine with `JULIA_NUM_THREADS=16` will pin a 16-worker plan in an 8-thread session (the budget is part of the signature, so this only happens if the signature scheme changes). Unknown mode strings are dropped silently rather than warned about.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 160.
