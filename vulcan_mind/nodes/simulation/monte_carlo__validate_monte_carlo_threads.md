---
id: simulation.monte_carlo__validate_monte_carlo_threads
label: _validate_monte_carlo_threads
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: _validate_monte_carlo_threads
  lines:
  - 65
  - 65
inputs:
- id: threads
  type: Int
  units: n/a
  required: true
  description: Positional argument `threads`.
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
  type: Any
  units: n/a
  description: Return value of `_validate_monte_carlo_threads`. Returns `threads`.
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

# _validate_monte_carlo_threads

## Purpose
Rejects a Monte Carlo worker count that exceeds the threads Julia was actually started with, turning a silent oversubscription into an actionable error.

## Design & Implementation
Compares `threads` with `Base.Threads.nthreads()` and throws `ArgumentError` when the request is larger. The message quotes both numbers and tells the user the exact restart flag, `--threads=N`, because the runner cannot create Julia threads at runtime and a user seeing only a slow campaign would not otherwise know why. Returns `threads` unchanged so it composes inline.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `threads` | Int | n/a | yes | Positional argument `threads`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_validate_monte_carlo_threads`. Returns `threads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.run_monte_carlo|run_monte_carlo]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:259-259`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It checks only the upper bound; a request equal to `nthreads()` leaves nothing for inner parallelism and is accepted without comment. The check is against the default thread pool only and ignores interactive threads.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl` line 65.
