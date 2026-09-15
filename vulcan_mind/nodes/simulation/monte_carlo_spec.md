---
id: simulation.monte_carlo_spec
label: MonteCarloSpec
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: MonteCarloSpec
  lines:
  - 12
  - 21
outputs:
- id: spec
  type: MonteCarloSpec
  units: n/a
  description: Typed campaign specification containing sample count, sampler, execution
    policy, and aggregation settings.
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
inputs:
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
---

# MonteCarloSpec

## Purpose
`MonteCarloSpec` describes a repeated simulation campaign without executing it. It packages the number of samples, how each sample is generated, the function that evaluates a sample, randomization settings, and aggregation or routing options consumed by `run_monte_carlo`.

## Theory & Math
For samples `Y₁,…,Yₙ` generated from parameter draws `Xᵢ`, the campaign estimates statistics such as `mean(Y)` and spread from the sample collection. The specification controls `n` and the sampling map but does not claim a confidence interval or convergence guarantee by itself.

## Model & Assumptions
The sampler must produce inputs accepted by the sample function, and the requested sample count must be nonnegative. Random seeds and worker routing should be explicit when reproducibility matters. Aggregation assumes returned samples have a compatible structure or a caller-provided reducer.

## Design & Implementation
The record is declared near the beginning of `monte_carlo.jl`. `run_monte_carlo` reads the spec, combines it with route features and configuration, optionally uses `run_simulation` as the sample function, and packages outputs in `MonteCarloResult`. The spec keeps campaign description separate from process-pool state.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `spec` | MonteCarloSpec | n/a | — | Typed campaign specification containing sample count, sampler, execution policy, and aggregation settings. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.adaptive_routing__run_campaign_adaptive|_run_campaign_adaptive]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:244-244`
- [[simulation.run_monte_carlo|run_monte_carlo]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:299-299`
- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:182-182`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A valid spec can still request an infeasible simulation or a sampler that returns heterogeneous results. The record does not estimate sampling error, detect burn-in, or decide whether the chosen number of samples is scientifically sufficient. Serialization constraints apply when the spec crosses a process boundary.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl:1-21`.
