---
id: parallel.profile_definitions_parallelprofileconfig
label: ParallelProfileConfig
kind: struct
source:
  file: src/parallel/routing/profile_definitions.jl
  symbol: ParallelProfileConfig
  lines:
  - 26
  - 26
inputs:
- id: profile
  type: ParallelProfile
  units: n/a
  required: true
  description: Field `profile`.
- id: label
  type: String
  units: n/a
  required: true
  description: Field `label`.
- id: outer_backend
  type: Symbol
  units: n/a
  required: true
  description: Field `outer_backend`.
- id: inner_adaptive
  type: Bool
  units: n/a
  required: true
  description: Field `inner_adaptive`.
- id: outer_route_adaptive
  type: Bool
  units: n/a
  required: true
  description: Field `outer_route_adaptive`.
- id: density_mode
  type: String
  units: n/a
  required: true
  description: Field `density_mode`.
- id: control_mode
  type: String
  units: n/a
  required: true
  description: Field `control_mode`.
- id: thermal_mode
  type: String
  units: n/a
  required: true
  description: Field `thermal_mode`.
- id: multibody_mode
  type: String
  units: n/a
  required: true
  description: Field `multibody_mode`.
- id: effector_mode
  type: String
  units: n/a
  required: true
  description: Field `effector_mode`.
- id: inner_scheduler
  type: String
  units: n/a
  required: false
  description: Field `inner_scheduler` (default `"static"`).
- id: adaptive_window
  type: Int
  units: n/a
  required: false
  description: Field `adaptive_window` (default `8`).
- id: adaptive_control_tail_guard
  type: Bool
  units: n/a
  required: false
  description: Field `adaptive_control_tail_guard` (default `false`).
- id: adaptive_measured_reward
  type: Bool
  units: n/a
  required: false
  description: Field `adaptive_measured_reward` (default `false`).
- id: persistent_hints
  type: Bool
  units: n/a
  required: false
  description: Field `persistent_hints` (default `false`).
- id: persistent_state_persist
  type: Bool
  units: n/a
  required: false
  description: Field `persistent_state_persist` (default `false`).
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
  type: ParallelProfileConfig
  units: n/a
  description: Constructed `ParallelProfileConfig` (keyword constructor via @kwdef).
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

# ParallelProfileConfig

## Purpose
The fully resolved parallel execution settings a profile expands into, covering both the outer routing backend and the inner callback and right-hand-side threading policy.

## Design & Implementation
A `Base.@kwdef struct` of sixteen fields. Six are mandatory and describe what the profile is: the `ParallelProfile` enum value, its string `label`, the `outer_backend` symbol and three booleans covering inner and outer adaptivity. Five string mode fields — `density_mode`, `control_mode`, `thermal_mode`, `multibody_mode`, `effector_mode` — carry `off`, `auto` or `on` per workload family. The remaining six default to a conservative static configuration: `inner_scheduler` of `static`, `adaptive_window` of 8, and the tail guard, measured reward, persistent hints and persistent state flags all false, so only the richest profile opts into learned behaviour.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `profile` | ParallelProfile | n/a | yes | Field `profile`. |
| in | `label` | String | n/a | yes | Field `label`. |
| in | `outer_backend` | Symbol | n/a | yes | Field `outer_backend`. |
| in | `inner_adaptive` | Bool | n/a | yes | Field `inner_adaptive`. |
| in | `outer_route_adaptive` | Bool | n/a | yes | Field `outer_route_adaptive`. |
| in | `density_mode` | String | n/a | yes | Field `density_mode`. |
| in | `control_mode` | String | n/a | yes | Field `control_mode`. |
| in | `thermal_mode` | String | n/a | yes | Field `thermal_mode`. |
| in | `multibody_mode` | String | n/a | yes | Field `multibody_mode`. |
| in | `effector_mode` | String | n/a | yes | Field `effector_mode`. |
| in | `inner_scheduler` | String | n/a | no | Field `inner_scheduler` (default `"static"`). |
| in | `adaptive_window` | Int | n/a | no | Field `adaptive_window` (default `8`). |
| in | `adaptive_control_tail_guard` | Bool | n/a | no | Field `adaptive_control_tail_guard` (default `false`). |
| in | `adaptive_measured_reward` | Bool | n/a | no | Field `adaptive_measured_reward` (default `false`). |
| in | `persistent_hints` | Bool | n/a | no | Field `persistent_hints` (default `false`). |
| in | `persistent_state_persist` | Bool | n/a | no | Field `persistent_state_persist` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ParallelProfileConfig | n/a | — | Constructed `ParallelProfileConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.profile_config|profile_config]] · `callees` → `callers` · call · `src/parallel/routing/profile_definitions.jl:123-123`

**Downstream**

- `callees` → [[parallel.profile_definitions_parallel_profile_name|parallel_profile_name]] · `callers` · call · `src/parallel/routing/profile_definitions.jl:46-46`
<!-- vulcan:connections:end -->

## Limitations
The mode fields are free strings with no validated vocabulary, so a typo such as `aut` is carried silently into scheduling decisions rather than rejected at construction.

## Provenance
Mapped from `src/parallel/routing/profile_definitions.jl` line 26.
