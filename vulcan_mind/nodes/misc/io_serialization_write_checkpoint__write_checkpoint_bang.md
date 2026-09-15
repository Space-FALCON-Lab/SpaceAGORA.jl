---
id: misc.io_serialization_write_checkpoint__write_checkpoint_bang
label: _write_checkpoint!
kind: function
source:
  file: src/io/serialization/io_serialization.jl
  symbol: _write_checkpoint!
  lines:
  - 34
  - 60
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: IOConfig._checkpoint_paths supplying the data and manifest destinations
    for the active results directory.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: checkpoint
  type: Serialized+TOML
  units: n/a
  description: Julia-serialized state payload plus a TOML manifest carrying schema
    version, simulation time, solver mode, byte size and SHA-256.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
- checkpoint
- serialization
charts:
- misc
origin: agent
---

# _write_checkpoint!

## Purpose
`_write_checkpoint!` captures the integrator state partway through a long aerobraking simulation so the run can be resumed after a walltime limit, a preemption, or a crash. It is the write half of a three-function contract completed by `_load_checkpoint` and `_clear_checkpoint!`, all keyed off the same `IOConfig._checkpoint_paths(args)` pair of destinations.

## Model & Assumptions
The checkpoint stores simulation time `t`, the full state vector `u`, and the solver mode string, together with a schema version supplied by the caller. The state is captured with `deepcopy(u_state)` on the assumption that the solver will keep mutating its own buffers after the call returns; without that copy the serialized payload could race with continued integration. Resumption assumes the same schema version and a compatible solver configuration — the payload records the mode but does not itself validate compatibility.

## Design & Implementation
Two artifacts are written. The data file is a Julia `Serialization.serialize` dump of a named tuple carrying `schema_version`, `created_utc` from `now(UTC)`, `t`, `solver_mode` and `u`. The manifest is a TOML file recording the same schema version and timestamp plus `time_s`, `solver_mode`, `data_path`, `data_size_bytes` from `filesize`, and `data_sha256` from `_sha256_hex`. Both go through `_atomic_write_file`, whose temporary name embeds the process id, thread id and a random `UInt` so parallel workers never collide on a scratch file, and whose `finally` block removes the temporary if the writer or the move threw. Writing data before manifest means the digest in the manifest always describes a file that already exists.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | IOConfig._checkpoint_paths supplying the data and manifest destinations for the active results directory. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `checkpoint` | Serialized+TOML | n/a | — | Julia-serialized state payload plus a TOML manifest carrying schema version, simulation time, solver mode, byte size and SHA-256. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/execution.jl:387-387`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:387-387`

**Downstream**

- `callees` → [[io.io_config__checkpoint_paths|_checkpoint_paths]] · `callers` · call · `src/io/serialization/io_serialization.jl:35-35`
- `callees` → [[io.io_serialization__atomic_write_file|_atomic_write_file]] · `callers` · call · `src/io/serialization/io_serialization.jl:43-43`
- `callees` → [[io.io_serialization__sha256_hex|_sha256_hex]] · `callers` · call · `src/io/serialization/io_serialization.jl:54-54`
- `callees` → [[simulation.persistence__sha256_hex|_sha256_hex]] · `callers` · call · `src/io/serialization/io_serialization.jl:54-54`
- `callees` → [[simulation.resume_checkpoint__checkpoint_paths|_checkpoint_paths]] · `callers` · call · `src/io/serialization/io_serialization.jl:35-35`
<!-- vulcan:connections:end -->

## Limitations
`Serialization` output is not portable across Julia versions or across changes to the types stored inside the state vector, so a checkpoint is only reliably readable by the same build that produced it; the schema version guards intent but cannot repair a layout change. `_load_checkpoint` verifies only that `:t` and `:u` are present and never recomputes the manifest digest, so corruption detection is left to the caller. The function returns `nothing`, so callers learn the written paths from `IOConfig._checkpoint_paths`.

## Provenance
Mapped from `src/io/serialization/io_serialization.jl:34-60`, with the atomic writer at lines 10-26 and the loader at lines 62-81.
