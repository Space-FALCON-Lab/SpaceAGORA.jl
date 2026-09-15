---
id: parallel.profile_definitions_parse_parallel_profile
label: parse_parallel_profile
kind: function
source:
  file: src/parallel/routing/profile_definitions.jl
  symbol: parse_parallel_profile
  lines:
  - 79
  - 79
inputs:
- id: raw
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `raw`.
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
  type: ParallelProfile
  units: n/a
  description: Return value of `parse_parallel_profile`.
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

# parse_parallel_profile

## Purpose
Turns whatever a user, script or configuration file wrote for a parallel profile into the canonical enum value, accepting the historical spellings that are still in circulation.

## Design & Implementation
Three methods. The string method normalises through `_normalize_profile_token`, then tests membership against a tuple of accepted spellings per profile — `r0` also answers to `serial` and `true_serial`, `r1_a` to `outer_only` and `threads`, `r4` to `auto_adaptive`, and `r5` additionally absorbs the legacy `r4_full_auto` family that was renamed when `R5` was introduced. An unmatched token raises `ArgumentError` quoting the raw input and listing the seven canonical names. The `Symbol` method stringifies and delegates; the `ParallelProfile` method is the identity, which lets callers accept an already-parsed value at the same argument position.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | AbstractString | n/a | yes | Positional argument `raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ParallelProfile | n/a | — | Return value of `parse_parallel_profile`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.precompile_workload_run_spaceagora_precompile_workload|_run_spaceagora_precompile_workload]] · `callees` → `callers` · call · `src/precompile_workload.jl:42-42`
- [[parallel.profile_config|profile_config]] · `callees` → `callers` · feedback · `src/parallel/routing/profile_definitions.jl:121-121`
- [[parallel.profile_definitions__normalize_profile_token|_normalize_profile_token]] · `callees` → `callers` · feedback · `src/parallel/routing/profile_definitions.jl:74-74`

**Downstream**

- `callees` → [[parallel.profile_config|profile_config]] · `callers` · call · `src/parallel/routing/profile_definitions.jl:115-115`
- `callees` → [[parallel.profile_definitions__normalize_profile_token|_normalize_profile_token]] · `callers` · call · `src/parallel/routing/profile_definitions.jl:80-80`
<!-- vulcan:connections:end -->

## Limitations
The alias tables are literal tuples inside the branch, so every rename adds another permanent entry and nothing marks which aliases are deprecated; `R4_full_auto` remains bound as a constant alias of `R5`, so code holding that name sees the newer behaviour.

## Provenance
Mapped from `src/parallel/routing/profile_definitions.jl` line 79.
