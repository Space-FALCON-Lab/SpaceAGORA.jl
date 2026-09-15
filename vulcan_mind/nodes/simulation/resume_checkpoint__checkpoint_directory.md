---
id: simulation.resume_checkpoint__checkpoint_directory
label: _checkpoint_directory
kind: function
source:
  file: src/simulation/engine/resume_checkpoint.jl
  symbol: _checkpoint_directory
  lines:
  - 1
  - 1
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
  type: SimulationModel.IOConfig._checkpoint_directory
  units: n/a
  description: Return value of `_checkpoint_directory`. Returns `SimulationModel.IOConfig._checkpoint_directory(args)`.
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

# _checkpoint_directory

## Purpose
Resolves the on-disk directory that holds resume checkpoints for a simulation run, given the packed `args` parameter object that carries the run's IO configuration.

## Design & Implementation
Declared `@inline _checkpoint_directory(args) = SimulationModel.IOConfig._checkpoint_directory(args)`. It is a single-expression forwarding shim: the engine layer keeps a stable local name while the real path-resolution policy lives in `SimulationModel.IOConfig`, so checkpoint location rules can change without touching the propagation engine. No filesystem access, no mutation of `args`, and no argument validation happens here.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.IOConfig._checkpoint_directory | n/a | — | Return value of `_checkpoint_directory`. Returns `SimulationModel.IOConfig._checkpoint_directory(args)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_config__checkpoint_paths|_checkpoint_paths]] · `callees` → `callers` · call · `src/io/config/io_config.jl:28-28`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/resume_checkpoint.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the body is a bare delegation with an untyped `args` parameter, any failure mode surfaces from the IOConfig implementation: a missing or malformed IO configuration field raises there, not here, and the stack trace passes through this inlined frame without adding context. The function does not create the directory or check that it exists or is writable.

## Provenance
Mapped from `src/simulation/engine/resume_checkpoint.jl` line 1.
