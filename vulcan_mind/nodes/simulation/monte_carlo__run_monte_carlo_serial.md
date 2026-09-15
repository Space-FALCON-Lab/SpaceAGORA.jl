---
id: simulation.monte_carlo__run_monte_carlo_serial
label: _run_monte_carlo_serial
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: _run_monte_carlo_serial
  lines:
  - 109
  - 109
inputs:
- id: f
  type: Any
  units: n/a
  required: true
  description: Positional argument `f`.
- id: seeds
  type: Vector
  units: n/a
  required: true
  description: Positional argument `seeds`.
- id: spec
  type: MonteCarloSpec
  units: n/a
  required: true
  description: Positional argument `spec`.
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
  description: Return value of `_run_monte_carlo_serial`. Returns `samples`.
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

# _run_monte_carlo_serial

## Purpose
The single-worker execution path for a Monte Carlo campaign, used when one worker is requested or when there is only one seed.

## Design & Implementation
Iterates `enumerate(seeds)` in order, calling `_run_monte_carlo_sample` for each and pushing the result. When `spec.fail_fast` is set and a sample fails it immediately calls `_throw_first_monte_carlo_failure` on the samples gathered so far, which throws. Because execution is sequential, the samples vector is trivially in index order and every seed before the failure has a recorded result.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Any | n/a | yes | Positional argument `f`. |
| in | `seeds` | Vector | n/a | yes | Positional argument `seeds`. |
| in | `spec` | MonteCarloSpec | n/a | yes | Positional argument `spec`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_run_monte_carlo_serial`. Returns `samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.run_monte_carlo|run_monte_carlo]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:264-264`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:113-113`
- `callees` → [[simulation.monte_carlo__run_monte_carlo_sample|_run_monte_carlo_sample]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:112-112`
- `callees` → [[simulation.monte_carlo__throw_first_monte_carlo_failure|_throw_first_monte_carlo_failure]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:115-115`
<!-- vulcan:connections:end -->

## Limitations
There is no progress reporting or checkpointing, so a long serial campaign that is interrupted loses every result accumulated in the local vector.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl` line 109.
