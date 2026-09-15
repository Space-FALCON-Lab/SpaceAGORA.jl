---
id: simulation.state_access__state_has_quaternion
label: _state_has_quaternion
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _state_has_quaternion
  lines:
  - 90
  - 90
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: Bool
  units: n/a
  description: Return value of `_state_has_quaternion`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _state_has_quaternion

## Purpose
Boolean capability check telling saving and guidance code whether attitude information exists for spacecraft `sat_idx` in the current integrator state.

## Design & Implementation
Defined as `!isnothing(_state_quaternion(u, sat_idx))`, so it reuses the full extraction path rather than duplicating the layout tests. That means the predicate is false both when the gravity-backbone partitioned state is active and when the per-spacecraft record simply has no `:q` field, covering translation-only spacecraft in an otherwise attitude-capable run.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_state_has_quaternion`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`

**Downstream**

- `callees` → [[simulation.state_access__state_quaternion|_state_quaternion]] · `callers` · call · `src/simulation/engine/state_access.jl:91-91`
<!-- vulcan:connections:end -->

## Limitations
Because it calls the extractor, checking the predicate builds an `SVector{4, Float64}` that is then discarded, so the common pattern of testing and then fetching does the conversion work twice. There is no caching and no type-level dispatch, so the cost recurs on every call inside a saving loop. It also inherits the extractor's behaviour of throwing when `u.sc[sat_idx]` itself is malformed.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 90.
