---
id: simulation.rhs_calibration__calib_machine_label
label: _calib_machine_label
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _calib_machine_label
  lines:
  - 47
  - 47
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
  description: Return value of `_calib_machine_label`.
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

# _calib_machine_label

## Purpose
Returns a short, stable machine fingerprint string used to key the persisted RHS calibration file and the calibration signature, so timing results from one host are never reused on another.

## Design & Implementation
The result is memoised in the global `_CALIB_MACHINE_LABEL::Ref{String}`; the body only runs when that ref is empty. If the environment variable `SPACEAGORA_PERF_MACHINE_LABEL` is non-empty its stripped value is passed through `SimulationModel.ParallelPolicy._safe_token` and used verbatim. Otherwise the label is derived from `Sys.cpu_info()[1].model` (or `"unknown"` when the vector is empty) concatenated with `Sys.CPU_THREADS`, hashed with `SHA.sha256`, and truncated to the first 8 bytes rendered as 16 hex characters via `bytes2hex`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_calib_machine_label`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__rhs_calib_path|_rhs_calib_path]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:101-101`
- [[simulation.rhs_calibration__rhs_calib_signature|_rhs_calib_signature]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:84-84`

**Downstream**

- `callees` → [[parallel.env_config__safe_token|_safe_token]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:51-51`
- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:49-49`
<!-- vulcan:connections:end -->

## Limitations
The memo is a plain `Ref` with no lock, so two threads calling it for the first time may both compute the label (the values are identical, so this is benign). Changing `SPACEAGORA_PERF_MACHINE_LABEL` after the first call has no effect until the process restarts. An 8-byte hash prefix gives 64 bits of discrimination, which is ample for machine labels but is not a full digest.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 47.
