---
id: parallel.env_mapping__env_or_default
label: _env_or_default
kind: function
source:
  file: src/parallel/routing/env_mapping.jl
  symbol: _env_or_default
  lines:
  - 45
  - 45
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: fallback
  type: String
  units: n/a
  required: true
  description: Positional argument `fallback`.
- id: preserve_existing
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `preserve_existing`.
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
  description: Return value of `_env_or_default`.
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

# _env_or_default

## Purpose

Decides, per environment variable, whether an operator-supplied value already present in `ENV` should win over the value implied by the parallel profile. Given a variable `name`, a `fallback` string and the mandatory keyword `preserve_existing`, it returns the string that should be set.

## Design & Implementation

When `preserve_existing` is `false` the `fallback` is returned immediately without reading `ENV`. Otherwise the current value is fetched with `get(ENV, name, "")` and `strip`ped; a non-empty result is returned, and an empty or whitespace-only one falls back. This makes an explicitly blank environment variable indistinguishable from an unset one, which is the intended behaviour for shell exports that were cleared. Every entry produced by `profile_env_pairs` passes through this gate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `fallback` | String | n/a | yes | Positional argument `fallback`. |
| in | `preserve_existing` | Bool | n/a | yes | Keyword argument `preserve_existing`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_env_or_default`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callees` → `callers` · feedback · `src/parallel/routing/env_mapping.jl:70-70`

**Downstream**

- `callees` → [[parcore.env_mapping_profile_env_pairs|profile_env_pairs]] · `callers` · call · `src/parallel/routing/env_mapping.jl:54-54`
<!-- vulcan:connections:end -->

## Limitations

Only emptiness is checked - no validation that the preserved value is a legal token for that variable, so a typo in an operator override survives into the run. The stripped value is returned rather than the raw one, so leading and trailing whitespace is dropped. `SPACEAGORA_RHS_BATCH_PARALLEL` deliberately bypasses preservation under profile `R0` by having its caller pass `preserve_existing=false`.

## Provenance
Mapped from `src/parallel/routing/env_mapping.jl` line 45.
