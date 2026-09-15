---
id: simulation.setup__typed_checkpoint_enabled
label: _typed_checkpoint_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _typed_checkpoint_enabled
  lines:
  - 33
  - 33
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
  type: Any
  units: n/a
  description: Return value of `_typed_checkpoint_enabled`. Returns `args.simulation_settings.checkpoint_enabled
    || args.simulation_settings.resume_f`.
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

# _typed_checkpoint_enabled

## Purpose
Reports whether checkpoint machinery must be set up for a run, which is the case either when periodic checkpoints are requested or when the run resumes from an existing checkpoint.

## Design & Implementation
Takes the `args` bundle and returns `args.simulation_settings.checkpoint_enabled || args.simulation_settings.resume_from_checkpoint`. Both fields are `Bool`s on `simulation_settings`; the short-circuit `||` means a resume-only run still allocates checkpoint buffers even if it never writes a new one. No environment variable is consulted, in contrast to the other `_typed_*` gates.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_typed_checkpoint_enabled`. Returns `args.simulation_settings.checkpoint_enabled \|\| args.simulation_settings.resume_f`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:229-229`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It does not verify that a checkpoint path exists or is writable, nor that `CHECKPOINT_SCHEMA_VERSION` matches an on-disk file; those checks happen in the loader. Accessing a missing field on a custom `args` type raises a `FieldError` at call time.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 33.
