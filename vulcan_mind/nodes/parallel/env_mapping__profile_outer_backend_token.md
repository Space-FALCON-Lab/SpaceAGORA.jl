---
id: parallel.env_mapping__profile_outer_backend_token
label: _profile_outer_backend_token
kind: function
source:
  file: src/parallel/routing/env_mapping.jl
  symbol: _profile_outer_backend_token
  lines:
  - 5
  - 5
inputs:
- id: backend
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `backend`.
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
  description: Return value of `_profile_outer_backend_token`.
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

# _profile_outer_backend_token

## Purpose

Maps the outer-parallelism backend `Symbol` held in a `ParallelProfileConfig` to the string written into `SPACEAGORA_PERF_PARALLEL_BACKEND`. The four accepted symbols are `:none`, `:threads`, `:process` and `:auto`, each rendered as its lowercase name.

## Design & Implementation

An `@inline` chain of equality tests returns the matching literal token; falling off the end throws `ArgumentError("Unsupported outer backend '$backend'.")`, so an unrecognised symbol fails at profile-expansion time rather than silently producing an empty or misspelled backend selection downstream. The explicit literals keep the environment spelling independent of how the symbols are named internally.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `backend` | Symbol | n/a | yes | Positional argument `backend`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_profile_outer_backend_token`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:77-77`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The accepted set is hard-coded, so adding a backend requires editing this function as well as the profile configuration. No aliasing or case-normalisation is performed: `:Threads` or `:THREADS` are rejected. The error is raised eagerly even when the resulting variable would have been overridden by a preserved existing value.

## Provenance
Mapped from `src/parallel/routing/env_mapping.jl` line 5.
