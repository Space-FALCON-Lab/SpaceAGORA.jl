---
id: parallel.outer_route_state_reset_outer_route_state_bang
label: reset_outer_route_state!
kind: function
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: reset_outer_route_state!
  lines:
  - 84
  - 84
inputs:
- id: state
  type: OuterRouteState
  units: n/a
  required: true
  description: Positional argument `state`.
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
  description: Return value of `reset_outer_route_state!`; mutates `state` in place.
    Returns `nothing`.
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

# reset_outer_route_state!

## Purpose
Clears all accumulated adaptive outer-route history so subsequent routing decisions start from an empty prior. Used by tests and by callers that want to discard stale timings, for example after a hardware or configuration change.

## Design & Implementation
Takes an `OuterRouteState`, acquires `state.lock` (a `ReentrantLock`) with the `lock(f, l)` form, and calls `empty!(state.history)` on the `Dict{String, Dict{Symbol, OuterRouteStats}}`. Returns `nothing`. The lock is released even if `empty!` throws, because the do-block form uses try/finally internally.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | OuterRouteState | n/a | yes | Positional argument `state`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `reset_outer_route_state!`; mutates `state` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.outer_route_state_outerroutestate|OuterRouteState]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:80-80`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the in-memory dict is cleared; any TOML file previously written by `save_outer_route_state` is untouched, so a later `load_outer_route_state!` with `replace=true` reinstates the old history. Because `empty!` keeps the dict's allocated capacity, memory is not returned. There is no way to reset a single signature or route.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl` line 84.
