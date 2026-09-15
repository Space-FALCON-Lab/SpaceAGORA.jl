---
id: parallel.persistent_hints_hint_layer_stats_snapshot
label: hint_layer_stats_snapshot
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: hint_layer_stats_snapshot
  lines:
  - 384
  - 384
inputs:
- id: profile
  type: Union{Nothing, AbstractString}
  units: n/a
  required: false
  description: Keyword argument `profile` (default `nothing`).
- id: machine
  type: Union{Nothing, AbstractString}
  units: n/a
  required: false
  description: Keyword argument `machine` (default `nothing`).
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
  type: Nothing
  units: n/a
  description: Return value of `hint_layer_stats_snapshot`. Returns `nothing` or `rows`.
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

# hint_layer_stats_snapshot

## Purpose
Aggregates the whole hint history into per-profile, per-machine, per-source rows so a report can show how much the bandit has learned and how confident it is.

## Design & Implementation
Loads state, normalises optional `profile` and `machine` filters through `_safe_token`, and under the lock walks every signature. Each is parsed for its profile, machine and source tokens, filtered, and its non-empty allotment entries folded into a `_HintLayerStatsAccumulator` keyed by that triple: totals of samples, successes, failures and elapsed sums, plus the best lower-confidence score, its mean and width, and the best observed mean, from which per-signature regret and confidence are accumulated. Rows are emitted in sorted key order with derived mean, standard deviation, mean confidence and mean regret, the exploration constant, minimum samples and the state path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `profile` | Union{Nothing, AbstractString} | n/a | no | Keyword argument `profile` (default `nothing`). |
| in | `machine` | Union{Nothing, AbstractString} | n/a | no | Keyword argument `machine` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `hint_layer_stats_snapshot`. Returns `nothing` or `rows`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/persistent_hints.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:419-419`
- `callees` → [[parallel.env_config__safe_token|_safe_token]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:389-389`
- `callees` → [[parallel.env_config_persistent_hints_exploration|persistent_hints_exploration]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:391-391`
- `callees` → [[parallel.env_config_persistent_hints_min_samples|persistent_hints_min_samples]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:392-392`
- `callees` → [[parallel.persistent_hints__ensure_persistent_hint_state_loaded_bang|_ensure_persistent_hint_state_loaded!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:388-388`
- `callees` → [[parallel.persistent_hints__hint_mean_and_width|_hint_mean_and_width]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:441-441`
- `callees` → [[parallel.persistent_hints__hint_signature_value|_hint_signature_value]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:404-404`
- `callees` → [[parallel.types__hintlayerstatsaccumulator|_HintLayerStatsAccumulator]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:426-426`
<!-- vulcan:connections:end -->

## Limitations
Variance is computed as mean of squares minus square of mean, which cancels catastrophically when timings are large and tightly clustered; a source token that is not a valid symbol is still converted with `Symbol`, so garbage in the file becomes a garbage layer name rather than an error.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 384.
