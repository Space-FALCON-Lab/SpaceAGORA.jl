---
id: simulation.rhs_calibration__rhs_calib_path
label: _rhs_calib_path
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calib_path
  lines:
  - 94
  - 94
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
  type: String
  units: n/a
  description: Return value of `_rhs_calib_path`.
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

# _rhs_calib_path

## Purpose
Computes the absolute filesystem path of the TOML file that persists RHS plan calibrations for the current machine.

## Design & Implementation
If `SPACEAGORA_RHS_CALIBRATION_PATH` is set and non-empty after `strip`, it is returned through `normpath`, joined onto `pwd()` when relative. Otherwise the default is `normpath(joinpath(pwd(), "output", "parallel_policy_state", "rhs_calibration_<label>.toml"))` where `<label>` is `_calib_machine_label()`. Pure aside from the memoised label lookup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_rhs_calib_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__rhs_calib_load_bang|_rhs_calib_load!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:109-109`
- [[simulation.rhs_calibration__rhs_calib_save_bang|_rhs_calib_save!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:131-131`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:95-95`
- `callees` → [[simulation.rhs_calibration__calib_machine_label|_calib_machine_label]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:101-101`
<!-- vulcan:connections:end -->

## Limitations
The default is anchored on `pwd()` at call time, so running the same simulation from two working directories produces two independent caches and changing directory mid-session after the first load will make `_rhs_calib_save!` write to a different file than was loaded. No check is made that the directory is writable; that failure surfaces later as a `@warn` in `_rhs_calib_save!`.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 94.
