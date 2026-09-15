---
id: parallel.profile_config
label: profile_config
kind: function
source:
  file: src/parallel/routing/profile_definitions.jl
  symbol: profile_config
  lines:
  - 120
  - 219
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: profile_in
  type: ParallelProfile
  units: n/a
  required: true
  description: Profile record or profile selector that is normalized into route configuration.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: config
  type: NamedTuple
  units: n/a
  description: Normalized parallel configuration consumed by route selection and worker
    setup.
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
---

# profile_config

## Purpose
`profile_config` converts a profile name, profile record, or environment-derived setting into the normalized configuration used by the parallel campaign path. It is the boundary where defaults become explicit values, which keeps `select_outer_route!` independent from TOML and environment parsing details.

## Theory & Math
This function performs discrete configuration resolution. For each option it applies a precedence rule: explicit caller data takes priority, followed by profile defaults and then environment-backed defaults. The output is a finite configuration vector; no continuous simulation quantity is transformed here.

## Model & Assumptions
The profile vocabulary and environment keys must match the definitions in `profile_definitions.jl` and `env_mapping.jl`. A missing option is expected to receive a documented default. Values that represent counts or thresholds must be parseable and within the range accepted by downstream worker creation.

## Design & Implementation
The implementation reads profile fields, applies named defaults, and returns the configuration consumed by `select_outer_route!`. The same normalized values are used by `with_parallel_profile`, so scoped environment changes and route decisions describe one policy. The profile parser also supports TOML-backed configuration, allowing a campaign manifest to reproduce a prior route choice.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `profile_in` | ParallelProfile | n/a | yes | Profile record or profile selector that is normalized into route configuration. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `config` | NamedTuple | n/a | — | Normalized parallel configuration consumed by route selection and worker setup. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.profile_definitions_parse_parallel_profile|parse_parallel_profile]] · `callees` → `callers` · call · `src/parallel/routing/profile_definitions.jl:115-115`
- [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:64-64`

**Downstream**

- `callees` → [[parallel.profile_definitions_parallelprofileconfig|ParallelProfileConfig]] · `callers` · call · `src/parallel/routing/profile_definitions.jl:123-123`
- `callees` → [[parallel.profile_definitions_parse_parallel_profile|parse_parallel_profile]] · `callers` · feedback · `src/parallel/routing/profile_definitions.jl:121-121`
<!-- vulcan:connections:end -->

## Limitations
Normalization does not verify actual process availability or benchmark the selected route. Environment variables are external mutable state and can change behavior between calls if callers do not use `with_parallel_profile`. A syntactically valid configuration can still request an inefficient or unsupported combination of processes and threads.

## Provenance
Mapped from `src/parallel/routing/profile_definitions.jl:114-219`.
