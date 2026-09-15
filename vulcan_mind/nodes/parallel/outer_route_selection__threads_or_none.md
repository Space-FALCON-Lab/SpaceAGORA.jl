---
id: parallel.outer_route_selection__threads_or_none
label: _threads_or_none
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _threads_or_none
  lines:
  - 1
  - 1
inputs:
- id: threads_available
  type: Bool
  units: n/a
  required: true
  description: Positional argument `threads_available`.
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
  type: Symbol
  units: n/a
  description: Return value of `_threads_or_none`.
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

# _threads_or_none

## Purpose
Tiny inlined predicate-to-route mapper that turns a `threads_available::Bool` flag into the route symbol `:threads` or `:none`, used wherever the router would prefer thread parallelism but must fall back to serial execution on single-threaded Julia sessions.

## Design & Implementation
Declared `@inline function _threads_or_none(threads_available::Bool)::Symbol` and returns `threads_available ? :threads : :none`. It is called by `_priority_outer_route_montecarlo`, `default_outer_route` (four branches), and indirectly through candidate construction, so the threads-or-serial fallback rule lives in one place.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `threads_available` | Bool | n/a | yes | Positional argument `threads_available`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_threads_or_none`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.default_outer_route|default_outer_route]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:348-348`
- [[parallel.outer_route_selection__priority_outer_route_montecarlo|_priority_outer_route_montecarlo]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:303-303`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It never returns `:process`; callers must decide separately whether process isolation applies. `threads_available` is a caller-supplied flag rather than a query of `Threads.nthreads()`, so a caller passing `true` on a one-thread session yields a `:threads` route that will run serially.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 1.
