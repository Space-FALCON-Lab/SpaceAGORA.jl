---
id: simulation.resume_checkpoint__load_checkpoint
label: _load_checkpoint
kind: function
source:
  file: src/simulation/engine/resume_checkpoint.jl
  symbol: _load_checkpoint
  lines:
  - 14
  - 14
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
  type: SimulationModel.IOSerialization._load_checkpoint
  units: n/a
  description: Return value of `_load_checkpoint`. Returns `SimulationModel.IOSerialization._load_checkpoint(args)`.
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

# _load_checkpoint

## Purpose
Reads a previously written resume checkpoint back for the run described by `args`, giving the engine the saved time, state vector, and solver mode needed to restart propagation mid-flight.

## Design & Implementation
Written as `@inline _load_checkpoint(args) = SimulationModel.IOSerialization._load_checkpoint(args)`. Deserialization, schema-version checking against `CHECKPOINT_SCHEMA_VERSION`, and the decision about what a missing file means are all implemented in `SimulationModel.IOSerialization`; this engine-side alias exists so that the write path (`_write_checkpoint!`, which does pass `CHECKPOINT_SCHEMA_VERSION` and a `solver_mode` keyword) and the read path sit next to each other in one small file.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.IOSerialization._load_checkpoint | n/a | — | Return value of `_load_checkpoint`. Returns `SimulationModel.IOSerialization._load_checkpoint(args)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:237-237`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Unlike the sibling `_write_checkpoint!` in this file, no schema version is passed at the call site, so version negotiation is entirely the serializer's responsibility and a checkpoint written by an incompatible build is only detected downstream. A truncated or partially flushed checkpoint file surfaces as a deserialization error rather than a clean resume-unavailable result, and there is no locking against a concurrent writer.

## Provenance
Mapped from `src/simulation/engine/resume_checkpoint.jl` line 14.
