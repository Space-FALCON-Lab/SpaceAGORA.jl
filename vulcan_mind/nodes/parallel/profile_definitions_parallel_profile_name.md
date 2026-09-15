---
id: parallel.profile_definitions_parallel_profile_name
label: parallel_profile_name
kind: function
source:
  file: src/parallel/routing/profile_definitions.jl
  symbol: parallel_profile_name
  lines:
  - 50
  - 50
inputs:
- id: profile
  type: ParallelProfile
  units: n/a
  required: true
  description: Positional argument `profile`.
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
  description: Return value of `parallel_profile_name`.
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

# parallel_profile_name

## Purpose
Gives the canonical string label for a `ParallelProfile` enum value, so telemetry, logs and result filenames all name a profile identically.

## Design & Implementation
An `@inline` chain of equality tests against `R0`, `R1_a`, `R1_b`, `R2`, `R3` and `R4`, returning the matching literal, with `R5` as the unconditional fallthrough rather than a seventh branch. Returning canonical names from one place is what makes the round trip through `parse_parallel_profile` stable, since that parser accepts many aliases but this emits only one form per profile.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `profile` | ParallelProfile | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `parallel_profile_name`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.profile_definitions_parallelprofileconfig|ParallelProfileConfig]] · `callees` → `callers` · call · `src/parallel/routing/profile_definitions.jl:46-46`
- [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:72-72`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `R5` is the fallthrough rather than a tested branch, a newly added enum variant would be mislabelled as `R5` instead of raising, so extending the enum requires editing this chain in step.

## Provenance
Mapped from `src/parallel/routing/profile_definitions.jl` line 50.
