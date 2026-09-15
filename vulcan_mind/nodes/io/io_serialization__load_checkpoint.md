---
id: io.io_serialization__load_checkpoint
label: _load_checkpoint
kind: function
source:
  file: src/io/serialization/io_serialization.jl
  symbol: _load_checkpoint
  lines:
  - 62
  - 62
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_load_checkpoint`. Returns `nothing` or `(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
charts:
- io
origin: agent
---

# _load_checkpoint

## Purpose
Reads back a previously written checkpoint so a simulation can resume from the recorded time and state, returning nothing when no checkpoint exists.

## Design & Implementation
Resolves the pair of paths through `IOConfig._checkpoint_paths` and returns `nothing` immediately if the data file is absent, which is the normal first-run case rather than an error. It deserializes the payload inside a `do` block, then requires the `:t` and `:u` keys and raises `ArgumentError` if either is missing. `solver_mode` is optional and is normalised to either `nothing` or a `String`. The returned named tuple carries the time as `Float64`, the raw state, the solver mode and both paths so the caller can report where the resume came from.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_load_checkpoint`. Returns `nothing` or `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/execution.jl:237-237`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:237-237`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/io/serialization/io_serialization.jl:75-75`
- `callees` → [[grp.src_analysis_verification|analysis/verification/]] · `members_in` · call · `src/io/serialization/io_serialization.jl:75-75`
- `callees` → [[io.io_config__checkpoint_paths|_checkpoint_paths]] · `callers` · call · `src/io/serialization/io_serialization.jl:63-63`
- `callees` → [[simulation.resume_checkpoint__checkpoint_paths|_checkpoint_paths]] · `callers` · call · `src/io/serialization/io_serialization.jl:63-63`
<!-- vulcan:connections:end -->

## Limitations
The manifest's recorded size and SHA-256 are not checked before deserializing, so a truncated or tampered data file reaches `deserialize` and fails there rather than being rejected cleanly; no schema version check is performed against the payload either.

## Provenance
Mapped from `src/io/serialization/io_serialization.jl` line 62.
