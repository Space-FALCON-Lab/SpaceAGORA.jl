---
id: parcore.env_mapping_profile_env_pairs
label: profile_env_pairs
kind: function
source:
  file: src/parallel/routing/env_mapping.jl
  symbol: profile_env_pairs
  lines:
  - 59
  - 174
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelProfiles namespace supplying profile_config, parallel_profile_name
    and the profile enumeration.
- id: profile_in
  type: ParallelProfile or Symbol
  units: n/a
  required: true
  description: Profile record or profile selector, normalised through profile_config
    before the mapping is built.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: env_pairs
  type: Vector{Pair{String,String}}
  units: n/a
  description: Ordered list of SPACEAGORA_* variable names and values implied by the
    profile, ready to apply to ENV or to pass to a worker process.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# profile_env_pairs

## Purpose
`profile_env_pairs` translates a parallel profile into the concrete set of `SPACEAGORA_*` environment variables that realise it. It is the bridge between the profile abstraction the user selects and the environment-variable surface the policy layer actually reads.

## Model & Assumptions
Every variable is emitted through `_env_or_default`, which by default preserves a value the user has already set and otherwise writes the profile's value. That preserve-existing behaviour makes the profile a baseline rather than an override, so an explicit variable on the command line still wins. The R0 profile is the documented exception: because it is defined as truly serial, the RHS batch setting is not preserved for it, and its default token becomes off instead of auto.

## Design & Implementation
The function first normalises its argument with `profile_config`, computes the inner hint defaults for that configuration, and then builds a `Vector{Pair{String,String}}` literal covering the profile name, the outer parallel backend token, the adaptive outer-route flag, the outer-active flag, the inner adaptive policy flag and the per-callback parallel modes for density, control, thermal and multibody work. Booleans are converted with `_coerce_env_bool` so the emitted values match the parser's accepted vocabulary. The same file defines two methods of `with_parallel_profile`, which apply the pairs for the duration of a block and restore the previous environment afterwards.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelProfiles namespace supplying profile_config, parallel_profile_name and the profile enumeration. |
| in | `profile_in` | ParallelProfile or Symbol | n/a | yes | Profile record or profile selector, normalised through profile_config before the mapping is built. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `env_pairs` | Vector{Pair{String,String}} | n/a | — | Ordered list of SPACEAGORA_* variable names and values implied by the profile, ready to apply to ENV or to pass to a worker process. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.env_mapping__env_or_default|_env_or_default]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:54-54`
- [[parallel.env_mapping_with_parallel_profile|with_parallel_profile]] · `callees` → `callers` · feedback · `src/parallel/routing/env_mapping.jl:181-181`

**Downstream**

- `callees` → [[parallel.env_mapping__coerce_env_bool|_coerce_env_bool]] · `callers` · call · `src/parallel/routing/env_mapping.jl:82-82`
- `callees` → [[parallel.env_mapping__env_or_default|_env_or_default]] · `callers` · feedback · `src/parallel/routing/env_mapping.jl:70-70`
- `callees` → [[parallel.env_mapping__inner_hint_defaults|_inner_hint_defaults]] · `callers` · call · `src/parallel/routing/env_mapping.jl:65-65`
- `callees` → [[parallel.env_mapping__profile_outer_backend_token|_profile_outer_backend_token]] · `callers` · call · `src/parallel/routing/env_mapping.jl:77-77`
- `callees` → [[parallel.env_mapping_with_parallel_profile|with_parallel_profile]] · `callers` · call · `src/parallel/routing/env_mapping.jl:169-169`
- `callees` → [[parallel.profile_config|profile_config]] · `callers` · call · `src/parallel/routing/env_mapping.jl:64-64`
- `callees` → [[parallel.profile_definitions_parallel_profile_name|parallel_profile_name]] · `callers` · call · `src/parallel/routing/env_mapping.jl:72-72`
<!-- vulcan:connections:end -->

## Limitations
The mapping is one-directional: it produces variables from a profile but cannot recover a profile from an arbitrary environment, so a partially overridden environment has no canonical profile name. Because preservation is keyed on the variable already being present, an empty-string value counts as set and suppresses the profile default. The pair list is a literal, so a variable added elsewhere in the policy layer is silently absent from profile mapping until it is added here too.

## Provenance
Mapped from `src/parallel/routing/env_mapping.jl:59-174`.
