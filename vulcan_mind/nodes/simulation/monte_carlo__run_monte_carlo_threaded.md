---
id: simulation.monte_carlo__run_monte_carlo_threaded
label: _run_monte_carlo_threaded
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: _run_monte_carlo_threaded
  lines:
  - 121
  - 121
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
- id: worker_count
  type: Int
  units: n/a
  required: true
  description: Positional argument `worker_count`.
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
  description: Return value of `_run_monte_carlo_threaded`. Returns `completed`.
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

# _run_monte_carlo_threaded

## Purpose
Runs Monte Carlo samples across a fixed number of thread workers pulling from a shared job queue, preserving seed order in the results.

## Design & Implementation
Fills a `Channel` with every `(index, seed)` pair up front and closes it, so workers terminate naturally when it drains. `worker_count` tasks are spawned with `Threads.@spawn` inside `@sync`; each writes its result into a preallocated `Vector{Union{Nothing, MonteCarloSampleResult}}` at the sample's own index, which is what keeps results ordered without a sort. An `Atomic{Bool}` `stop_requested` implements `fail_fast`: a failing worker sets it with `atomic_xchg!` and every worker checks it before taking the next job. After the sync the `nothing` slots are filtered out and, under `fail_fast`, the first failure is rethrown.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Any | n/a | yes | Positional argument `f`. |
| in | `seeds` | Vector | n/a | yes | Positional argument `seeds`. |
| in | `spec` | MonteCarloSpec | n/a | yes | Positional argument `spec`. |
| in | `worker_count` | Int | n/a | yes | Positional argument `worker_count`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_run_monte_carlo_threaded`. Returns `completed`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.run_monte_carlo|run_monte_carlo]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:266-266`

**Downstream**

- `callees` → [[simulation.monte_carlo__run_monte_carlo_sample|_run_monte_carlo_sample]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:136-136`
- `callees` → [[simulation.monte_carlo__throw_first_monte_carlo_failure|_throw_first_monte_carlo_failure]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:151-151`
<!-- vulcan:connections:end -->

## Limitations
`fail_fast` only stops workers at their next scheduling point, so samples already in flight run to completion and are discarded from the thrown error's perspective; the job channel is bounded to `length(seeds)` and fully populated before any worker starts, so memory scales with the seed count even for lazy iterables.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl` line 121.
