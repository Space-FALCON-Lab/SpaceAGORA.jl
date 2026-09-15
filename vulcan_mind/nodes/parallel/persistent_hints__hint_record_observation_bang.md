---
id: parallel.persistent_hints__hint_record_observation_bang
label: _hint_record_observation!
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_record_observation!
  lines:
  - 341
  - 341
inputs:
- id: signature
  type: String
  units: n/a
  required: true
  description: Positional argument `signature`.
- id: allotment
  type: Int64
  units: n/a
  required: true
  description: Positional argument `allotment`.
- id: elapsed_ns
  type: Int64
  units: n/a
  required: true
  description: Positional argument `elapsed_ns`.
- id: success
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `success`.
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
  description: Return value of `_hint_record_observation!`; mutates `signature` in
    place.
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

# _hint_record_observation!

## Purpose
Adds one timed outcome to the hint history for a signature and allotment, marking the state dirty so it will be saved.

## Design & Implementation
Ensures the state is loaded, returns early if hints are disabled or `allotment` is non-positive, and clamps the elapsed time at zero. Under `_persistent_hint_lock` it obtains or creates the bucket and stats record with `get!`, increments `samples` and either `successes` or `failures`, accumulates the elapsed time and its square, and sets `state.dirty`. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `signature` | String | n/a | yes | Positional argument `signature`. |
| in | `allotment` | Int64 | n/a | yes | Positional argument `allotment`. |
| in | `elapsed_ns` | Int64 | n/a | yes | Positional argument `elapsed_ns`. |
| in | `success` | Bool | n/a | yes | Keyword argument `success`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_hint_record_observation!`; mutates `signature` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:127-127`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:352-352`
- `callees` → [[parallel.persistent_hints__ensure_persistent_hint_state_loaded_bang|_ensure_persistent_hint_state_loaded!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:347-347`
- `callees` → [[parallel.types_adaptivechoicestats|AdaptiveChoiceStats]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:359-359`
<!-- vulcan:connections:end -->

## Limitations
Every observation takes the global lock, so heavily threaded workloads that observe at high frequency serialise on it; there is no cap on history size, so a long-lived process accumulates unboundedly many signatures.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 341.
