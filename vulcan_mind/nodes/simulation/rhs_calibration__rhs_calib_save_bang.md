---
id: simulation.rhs_calibration__rhs_calib_save_bang
label: _rhs_calib_save!
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calib_save!
  lines:
  - 128
  - 128
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
  type: Nothing
  units: n/a
  description: Return value of `_rhs_calib_save!`.
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

# _rhs_calib_save!

## Purpose
Atomically writes the entire in-process calibration cache back to the machine-specific TOML file so later sessions can skip the timing sweep.

## Design & Implementation
Under `_rhs_calib_lock`, returns early if `_rhs_calib_cache` is empty. Otherwise it builds a `Vector{Dict{String,Any}}` of rows sorted by signature, each with `signature`, `mode`, `allotment::Int` and `elapsed_mean_ns::Float64`, wraps them as `{"schema_version" => 1, "calibrations" => rows}`, calls `mkpath(dirname(path))`, writes via `TOML.print` to `path * ".tmp"`, then `mv(tmp, path; force=true)` for an atomic replace. Any exception in the I/O block is caught and reported with `@warn` including the path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_rhs_calib_save!`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:344-344`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:139-139`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:135-135`
- `callees` → [[simulation.rhs_calibration__rhs_calib_path|_rhs_calib_path]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:131-131`
<!-- vulcan:connections:end -->

## Limitations
The file is rewritten wholesale, so two concurrent Julia processes on the same machine race: the last writer wins and the other's fresh calibrations are lost. `schema_version` is written but never checked by `_rhs_calib_load!`. A leftover `.tmp` file from a crash between `open` and `mv` is not cleaned up on the next run.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 128.
