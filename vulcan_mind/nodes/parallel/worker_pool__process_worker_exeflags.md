---
id: parallel.worker_pool__process_worker_exeflags
label: _process_worker_exeflags
kind: function
source:
  file: src/parallel/process/worker_pool.jl
  symbol: _process_worker_exeflags
  lines:
  - 39
  - 39
inputs:
- id: project_path
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `project_path`.
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
  type: Cmd
  units: n/a
  description: Return value of `_process_worker_exeflags`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _process_worker_exeflags

## Purpose
Builds the command-line flags every campaign worker process is started with, fixing the three settings that make a worker reproducible and non-contending.

## Design & Implementation
Returns a `Cmd` of exactly three arguments: `--threads=1`, `--startup-file=no` and `--project=` interpolated with the pool's project path. The single thread is deliberate — process workers do not share the coordinator's thread pool, so giving each one thread means inner thread-based parallelism inside a worker's own `run_simulation` neither contends with nor is affected by how many process workers are active. Suppressing the startup file keeps a developer's personal `~/.julia/config` out of campaign results.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `project_path` | AbstractString | n/a | yes | Positional argument `project_path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Cmd | n/a | — | Return value of `_process_worker_exeflags`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.ensure_process_workers_ensure_process_workers_bang|ensure_process_workers!]] · `callees` → `callers` · call · `src/parallel/process/worker_pool.jl:198-198`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The single-thread choice is hard-coded, so a campaign wanting nested threading inside each process worker cannot express that through the pool; the project path is interpolated without quoting, so a path containing shell metacharacters relies on `Cmd` argument handling rather than escaping here.

## Provenance
Mapped from `src/parallel/process/worker_pool.jl` line 39.
