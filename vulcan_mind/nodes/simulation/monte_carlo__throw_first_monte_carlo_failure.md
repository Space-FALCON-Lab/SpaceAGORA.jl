---
id: simulation.monte_carlo__throw_first_monte_carlo_failure
label: _throw_first_monte_carlo_failure
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: _throw_first_monte_carlo_failure
  lines:
  - 102
  - 102
inputs:
- id: samples
  type: Vector{MonteCarloSampleResult}
  units: n/a
  required: true
  description: Positional argument `samples`.
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
  description: Return value of `_throw_first_monte_carlo_failure`.
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

# _throw_first_monte_carlo_failure

## Purpose
Implements the `fail_fast` contract by surfacing the earliest failed sample, by index order, as a thrown error.

## Design & Implementation
Uses `findfirst` on the `success` flag over the samples vector and returns `nothing` when every sample succeeded. Otherwise it throws an `ErrorException` interpolating the sample's `index`, `seed` and the stringified original `error`. Because callers pass samples that are already in index order, the first unsuccessful entry is the lowest-index failure, not the first one that happened to finish.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Vector{MonteCarloSampleResult} | n/a | yes | Positional argument `samples`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_throw_first_monte_carlo_failure`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/monte_carlo.jl`
- [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:189-189`
- [[simulation.monte_carlo__run_monte_carlo_process|_run_monte_carlo_process]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:203-203`
- [[simulation.monte_carlo__run_monte_carlo_serial|_run_monte_carlo_serial]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:115-115`
- [[simulation.monte_carlo__run_monte_carlo_threaded|_run_monte_carlo_threaded]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:151-151`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The original exception is flattened into a message string, so its type and the captured backtrace are lost to whoever catches the rethrow; only the message survives. It reports one failure even when several occurred.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl` line 102.
