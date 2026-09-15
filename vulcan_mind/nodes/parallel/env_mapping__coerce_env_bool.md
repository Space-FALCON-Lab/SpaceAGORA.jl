---
id: parallel.env_mapping__coerce_env_bool
label: _coerce_env_bool
kind: function
source:
  file: src/parallel/routing/env_mapping.jl
  symbol: _coerce_env_bool
  lines:
  - 1
  - 1
inputs:
- id: v
  type: Bool
  units: n/a
  required: true
  description: Positional argument `v`.
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
  description: Return value of `_coerce_env_bool`.
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

# _coerce_env_bool

## Purpose

Renders a Julia `Bool` as the textual form the SpaceAGORA environment variables use: `"1"` for `true` and `"0"` for `false`. It is the single place that fixes that spelling, so every boolean flag emitted by `profile_env_pairs` is written identically.

## Design & Implementation

A one-line `@inline` function with a declared `::String` return type and a `Bool` argument, so dispatch rejects anything that is not already a boolean rather than coercing it. It is applied to `cfg.outer_route_adaptive`, `cfg.inner_adaptive`, `cfg.adaptive_control_tail_guard`, `cfg.adaptive_measured_reward`, `cfg.persistent_hints`, `cfg.persistent_state_persist` and the `outer_parallel_active` keyword before those values are paired with their variable names.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v` | Bool | n/a | yes | Positional argument `v`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_coerce_env_bool`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:82-82`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

It only produces the `"1"`/`"0"` encoding; readers that accept `"true"`/`"false"` are not consulted here, so the decoding side must agree by convention. Passing an `Integer` or a `String` is a `MethodError`, which is deliberate but means callers holding a nullable flag must resolve it first.

## Provenance
Mapped from `src/parallel/routing/env_mapping.jl` line 1.
