---
id: simulation.rhs_calibration__rhs_calib_load_bang
label: _rhs_calib_load!
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calib_load!
  lines:
  - 105
  - 105
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
  description: Return value of `_rhs_calib_load!`.
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

# _rhs_calib_load!

## Purpose
Lazily reads the machine-specific calibration TOML file into the in-process `_rhs_calib_cache` exactly once per session.

## Design & Implementation
Runs under `lock(_rhs_calib_lock)`. The `_rhs_calib_loaded::Ref{Bool}` guard is set to `true` before any I/O so a failed load is never retried. The path comes from `_rhs_calib_path()`; a missing file, a `TOML.parsefile` exception, or a `calibrations` entry that is not an `AbstractVector` all return `nothing` silently. Each `AbstractDict` row with a non-empty `signature` is inserted into `_rhs_calib_cache[sig]` as a `Dict{String,Any}` holding `mode::String`, `allotment::Int` (default 1) and `elapsed_mean_ns::Float64` (default 0.0). Mutates `_rhs_calib_cache` and `_rhs_calib_loaded`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_rhs_calib_load!`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__rhs_calib_lookup|_rhs_calib_lookup]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:161-161`
- [[simulation.rhs_calibration__rhs_calib_store_bang|_rhs_calib_store!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:176-176`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:121-121`
- `callees` → [[simulation.rhs_calibration__rhs_calib_path|_rhs_calib_path]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:109-109`
<!-- vulcan:connections:end -->

## Limitations
Parse failures are swallowed with no log message, so a corrupt file looks identical to an absent one. Rows lacking `mode` load with an empty string and are later rejected by `_rhs_calib_lookup`. `Int(get(row, "allotment", 1))` throws `InexactError` if the TOML stores a non-integral float, and that exception is not caught. Entries stored in-process before the first load are overwritten only if `_rhs_calib_store!` did not already trigger the load, which it does deliberately.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 105.
