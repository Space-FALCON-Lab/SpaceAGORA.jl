---
id: parallel.parallel_profile
label: ParallelProfile
kind: function
source:
  file: src/parallel/routing/profile_definitions.jl
  symbol: ParallelProfile
  lines:
  - 2
  - 18
outputs:
- id: profile
  type: ParallelProfile
  units: n/a
  description: Immutable execution-profile record containing the selected process,
    thread, and routing settings.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
- configuration
charts:
- parallel
origin: agent
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
---

# ParallelProfile

## Purpose
`ParallelProfile` is the typed value that carries a campaign’s parallel execution policy. It gives the routing layer one record for process count, thread count, outer-route preference, and related toggles instead of passing loosely typed environment strings through simulation code.

## Theory & Math
The profile does not alter the equations of motion. It selects an execution configuration, so its output is a discrete policy state. Runtime cost is influenced by worker count and thread count, but numerical state evolution should be invariant when the route preserves deterministic ordering and floating-point reduction behavior.

## Model & Assumptions
Fields are interpreted by the companion profile parser and environment mapper. A profile is assumed to be internally consistent, with nonnegative worker settings and route names recognized by the selector. The record is configuration data, not a live worker pool and not proof that the requested resources exist.

## Design & Implementation
`profile_definitions.jl` declares the record near the top of the file and exposes it through `ParallelProfiles`. `profile_config` later turns the record into the normalized configuration consumed by route selection. `with_parallel_profile` uses the same fields to scope environment changes around a caller’s function, preserving the caller’s prior environment after execution.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `profile` | ParallelProfile | n/a | — | Immutable execution-profile record containing the selected process, thread, and routing settings. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/profile_definitions.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The record cannot detect scheduler limits, oversubscription, worker startup failure, or native-library thread safety. Invalid combinations may survive construction and fail only when a route is selected or a process pool is created. Configuration changes do not retroactively update an already-created worker pool.

## Provenance
Mapped from `src/parallel/routing/profile_definitions.jl:1-18`.
