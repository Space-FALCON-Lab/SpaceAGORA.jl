---
id: parallel.env_mapping_with_parallel_profile
label: with_parallel_profile
kind: function
source:
  file: src/parallel/routing/env_mapping.jl
  symbol: with_parallel_profile
  lines:
  - 175
  - 175
inputs:
- id: f
  type: Function
  units: n/a
  required: true
  description: Positional argument `f`.
- id: profile_in
  type: Any
  units: n/a
  required: true
  description: Positional argument `profile_in`.
- id: preserve_existing
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `preserve_existing` (default `true`).
- id: outer_parallel_active
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `outer_parallel_active` (default `false`).
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
  description: Return value of `with_parallel_profile`. Returns `withenv(env_pairs...)
    do`.
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

# with_parallel_profile

## Purpose

Runs a caller-supplied zero-argument function with the full set of `SPACEAGORA_*` environment variables implied by a parallel profile temporarily in place, restoring the previous environment when the function returns or throws. It is the scoped entry point for profile-driven parallel routing.

## Design & Implementation

Two methods exist so both argument orders work: `with_parallel_profile(f::Function, profile_in; ...)` and `with_parallel_profile(profile_in, f::Function; ...)`, the second simply forwarding to the first. The keywords `preserve_existing` (default `true`) and `outer_parallel_active` (default `false`) are passed to `profile_env_pairs`, whose vector of `Pair{String,String}` is splatted into `withenv`, which restores prior values on exit. The result of `f()` is returned unchanged.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Function | n/a | yes | Positional argument `f`. |
| in | `profile_in` | Any | n/a | yes | Positional argument `profile_in`. |
| in | `preserve_existing` | Bool | n/a | no | Keyword argument `preserve_existing` (default `true`). |
| in | `outer_parallel_active` | Bool | n/a | no | Keyword argument `outer_parallel_active` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `with_parallel_profile`. Returns `withenv(env_pairs...) do`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:169-169`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/routing/env_mapping.jl:187-187`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/routing/env_mapping.jl:187-187`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/routing/env_mapping.jl:187-187`
- `callees` → [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callers` · feedback · `src/parallel/routing/env_mapping.jl:181-181`
<!-- vulcan:connections:end -->

## Limitations

`withenv` mutates the process-wide environment for the duration of the call, so concurrent tasks on other threads observe the profile settings too and nesting two different profiles gives the inner one precedence. Because `profile_env_pairs` is evaluated eagerly, an unsupported outer backend raises `ArgumentError` before `f` is ever called. Worker processes spawned before the call do not inherit the temporary settings.

## Provenance
Mapped from `src/parallel/routing/env_mapping.jl` line 175.
