---
id: parallel.env_mapping__inner_hint_defaults
label: _inner_hint_defaults
kind: function
source:
  file: src/parallel/routing/env_mapping.jl
  symbol: _inner_hint_defaults
  lines:
  - 32
  - 32
inputs:
- id: cfg
  type: ParallelProfileConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: NamedTuple{(:exploration,
  units: n/a
  description: Return value of `_inner_hint_defaults`.
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

# _inner_hint_defaults

## Purpose

Supplies the default adaptive-policy hint parameters for a given `ParallelProfileConfig`, returning a `NamedTuple{(:exploration, :min_samples), Tuple{Float64, Int}}` that becomes `SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION` and `SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES`.

## Design & Implementation

Any profile other than `R5` gets the neutral pair `(exploration=1.5, min_samples=2)` without touching the machine. For `R5` the function calls `_machine_parallel_class()` and scales: `:large` hosts explore more aggressively and demand more evidence at `(1.8, 3)`, `:medium` keeps `(1.5, 2)`, and `:small` hosts back off to `(1.3, 2)`. The exploration value is later rounded to three digits by `profile_env_pairs` before stringification.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | ParallelProfileConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NamedTuple{(:exploration, | n/a | — | Return value of `_inner_hint_defaults`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:65-65`

**Downstream**

- `callees` → [[parallel.env_mapping__machine_parallel_class|_machine_parallel_class]] · `callers` · call · `src/parallel/routing/env_mapping.jl:36-36`
<!-- vulcan:connections:end -->

## Limitations

Hardware sensitivity is applied only to the `R5` profile; every other profile ignores machine size entirely. The exploration and sample constants are literals with no configuration surface of their own beyond the environment override consumed by `_machine_parallel_class`. The returned tuple is a starting point only - the adaptive policy may move away from it during a run.

## Provenance
Mapped from `src/parallel/routing/env_mapping.jl` line 32.
