---
id: parallel.persistent_hints__hint_samples
label: _hint_samples
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_samples
  lines:
  - 239
  - 239
inputs:
- id: bucket
  type: Any
  units: n/a
  required: true
  description: Positional argument `bucket`.
- id: candidate
  type: Int64
  units: n/a
  required: true
  description: Positional argument `candidate`.
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
  type: Int64
  units: n/a
  description: Return value of `_hint_samples`.
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

# _hint_samples

## Purpose
Reads the sample count for one candidate allotment from a signature bucket, treating absence as zero.

## Design & Implementation
Looks the candidate up with `get` and returns `max(0, stats.samples)` when the value is an `AdaptiveChoiceStats`, otherwise zero. The `isa` guard tolerates a bucket that somehow holds a foreign value rather than throwing during a decision. `@inline` with an `Int64` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `bucket` | Any | n/a | yes | Positional argument `bucket`. |
| in | `candidate` | Int64 | n/a | yes | Positional argument `candidate`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int64 | n/a | — | Return value of `_hint_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.persistent_hints__hint_choose_allotment|_hint_choose_allotment]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:281-281`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It cannot distinguish an allotment never tried from one tried but recorded with zero samples, which is fine for the exploration rule but loses information for diagnostics.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 239.
