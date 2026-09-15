---
id: simx.engine_resume_checkpoint_write_checkpoint__write_checkpoint_bang
label: _write_checkpoint!
kind: function
source:
  file: src/simulation/engine/resume_checkpoint.jl
  symbol: _write_checkpoint!
  lines:
  - 4
  - 12
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Configuration used to resolve the checkpoint directory and file paths
    for this run.
- id: time_t
  type: Float64
  units: s
  required: true
  description: Mission-elapsed time at the end of the completed segment, which becomes
    the resume point.
- id: state_u
  type: ComponentVector
  units: m,m/s,kg,J
  required: true
  description: Deep-copied solver state at that time, captured from the last accepted
    step of the segment.
- id: solver_mode
  type: String
  units: n/a
  required: true
  description: Stringified solver mode recorded so a resume can refuse to continue
    under a different integrator.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: checkpoint_path
  type: String
  units: n/a
  description: Path of the checkpoint file written atomically by the IO serialization
    layer.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# _write_checkpoint!

## Purpose
`_write_checkpoint!` records a resumable snapshot at the end of each integration segment. It is the engine-side wrapper that supplies `CHECKPOINT_SCHEMA_VERSION` to `SimulationModel.IOSerialization._write_checkpoint!` and forwards the solver mode as a keyword.

## Model & Assumptions
A checkpoint is only meaningful together with the integrator that produced it, so the solver mode string is part of the payload rather than an afterthought. The schema version is stamped so that a checkpoint written by an older build is rejected at load time instead of being reinterpreted against a changed state layout. The snapshot is a full state vector, not a delta, which makes resume independent of every earlier segment.

## Design & Implementation
The file is fifteen lines and contains only forwarding definitions: `_checkpoint_directory` and `_checkpoint_paths` resolve locations through `IOConfig`, while `_load_checkpoint` and `_clear_checkpoint!` complete the lifecycle through `IOSerialization`. The write helper is `@inline function` rather than a one-line alias because it inserts the schema constant between the caller's positional arguments and the keyword. `run_simulation` calls it after each checkpointed segment with `t_cursor` and a `deepcopy` of the segment's final state, and calls it again on the failure path so a crashed run still leaves a usable resume point.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `args` | SimulationConfiguration | n/a | yes | Configuration used to resolve the checkpoint directory and file paths for this run. |
| in | `time_t` | Float64 | s | yes | Mission-elapsed time at the end of the completed segment, which becomes the resume point. |
| in | `state_u` | ComponentVector | m,m/s,kg,J | yes | Deep-copied solver state at that time, captured from the last accepted step of the segment. |
| in | `solver_mode` | String | n/a | yes | Stringified solver mode recorded so a resume can refuse to continue under a different integrator. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `checkpoint_path` | String | n/a | — | Path of the checkpoint file written atomically by the IO serialization layer. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:387-387`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Checkpoint granularity is the segment length set by `simulation_settings.checkpoint_interval_s`; work done since the last segment boundary is lost on a crash. The snapshot captures the solver state but not adaptive integrator internals, so a resumed run restarts step-size control from scratch and its step sequence will not reproduce the uninterrupted run exactly. Stale checkpoints are only removed by an explicit `_clear_checkpoint!`.

## Provenance
Mapped from `src/simulation/engine/resume_checkpoint.jl:4-12`; `CHECKPOINT_SCHEMA_VERSION` is defined at `src/simulation/engine/setup.jl:19` and the segment loop calls it from `src/simulation/engine/execution.jl:387`.
