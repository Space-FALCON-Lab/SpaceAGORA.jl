---
id: simulation.run_monte_carlo
label: run_monte_carlo
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: run_monte_carlo
  lines:
  - 255
  - 300
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: sample_fn
  type: Function
  units: n/a
  required: true
  description: Callable that turns one sampled scenario into a simulation result.
- id: route_features
  type: NamedTuple
  units: n/a
  required: true
  description: Campaign measurements passed to adaptive route selection.
- id: spec
  type: MonteCarloSpec
  units: n/a
  required: true
  description: Campaign specification containing sampling and aggregation policy.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: MonteCarloResult
  units: n/a
  description: Collected samples and campaign metadata returned after execution.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
- monte-carlo
charts:
- simulation
origin: agent
---

# run_monte_carlo

## Purpose
`run_monte_carlo` executes the repeated-scenario workflow described by `MonteCarloSpec`. It selects an outer route from campaign features, prepares worker resources when needed, invokes the sample function for each draw, and returns a structured result for downstream analysis.

## Theory & Math
The campaign evaluates `Yᵢ = g(Xᵢ)` for `i = 1…N`, where `Xᵢ` comes from the configured sampler and `g` is `sample_fn`. The result stores the sample collection and metadata needed for estimates such as `Ȳ = (1/N)ΣYᵢ`. Parallel execution changes scheduling, not the intended sample definition.

## Model & Assumptions
The sample function must be callable in the selected execution environment and return values compatible with the result collector. Randomness must be seeded explicitly for repeatable studies. Route features should reflect the actual campaign size, and native simulation calls must still obey the shared runtime lock.

## Design & Implementation
The implementation reads the spec, calls adaptive route selection with `route_features`, and uses the selected route to prepare process workers. It then executes or dispatches sample calls, collects outputs, and constructs `MonteCarloResult`. The `run_simulation` public API is a natural sample function for engine-backed campaigns.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `sample_fn` | Function | n/a | yes | Callable that turns one sampled scenario into a simulation result. |
| in | `route_features` | NamedTuple | n/a | yes | Campaign measurements passed to adaptive route selection. |
| in | `spec` | MonteCarloSpec | n/a | yes | Campaign specification containing sampling and aggregation policy. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MonteCarloResult | n/a | — | Collected samples and campaign metadata returned after execution. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:174-174`
- [[simulation.monte_carlo__run_monte_carlo_process|_run_monte_carlo_process]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:209-209`
- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:196-196`

**Downstream**

- `callees` → [[simulation.adaptive_routing__run_monte_carlo_adaptive|_run_monte_carlo_adaptive]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:285-285`
- `callees` → [[simulation.monte_carlo__run_monte_carlo_serial|_run_monte_carlo_serial]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:264-264`
- `callees` → [[simulation.monte_carlo__run_monte_carlo_threaded|_run_monte_carlo_threaded]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:266-266`
- `callees` → [[simulation.monte_carlo__validate_monte_carlo_threads|_validate_monte_carlo_threads]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:259-259`
- `callees` → [[simulation.monte_carlo_result|MonteCarloResult]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:257-257`
- `callees` → [[simulation.monte_carlo_spec|MonteCarloSpec]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:299-299`
<!-- vulcan:connections:end -->

## Limitations
The function does not establish statistical convergence and cannot correct a biased sampler. Worker failures, serialization errors, solver failures, or inconsistent sample shapes can abort collection or produce a partial result. Adaptive route feedback may make runtime behavior differ between repeated executions even when the physical sample distribution is unchanged.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl:208-300`.
