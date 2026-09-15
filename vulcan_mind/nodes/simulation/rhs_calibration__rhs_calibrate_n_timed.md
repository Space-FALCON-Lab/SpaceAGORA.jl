---
id: simulation.rhs_calibration__rhs_calibrate_n_timed
label: _rhs_calibrate_n_timed
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calibrate_n_timed
  lines:
  - 40
  - 40
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
  type: Int
  units: n/a
  description: Return value of `_rhs_calibrate_n_timed`.
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

# _rhs_calibrate_n_timed

## Purpose
Reads how many timed `spacecraft_dynamics!` calls the calibration sweep should average per candidate plan.

## Design & Implementation
`@inline` accessor that parses `SPACEAGORA_RHS_CALIBRATE_N_TIMED` via `_engine_env_get` with default `"10"`, wrapping `parse(Int, strip(...))` in a `try`/`catch` that falls back to 10 on any parse failure, then clamps with `max(1, n)`. The value multiplies the number of candidates to give the total timed RHS evaluations performed before the solve starts.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_calibrate_n_timed`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__run_rhs_sweep_bang|_run_rhs_sweep!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:254-254`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:41-41`
<!-- vulcan:connections:end -->

## Limitations
Malformed values are silently replaced by 10 with no warning. There is no upper bound, so a very large setting makes the pre-solve sweep dominate runtime. The environment variable is re-read on every call rather than memoised, though the cost is negligible because it is called once per sweep.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 40.
