---
id: flow.campaigns
label: Campaigns
kind: group
inputs:
- id: campaign_spec
  type: MonteCarloSpec / closure
  units: n/a
  description: Seeds and the per-seed configuration builder.
- id: thread_policy
  type: PolicyDecision
  units: n/a
  description: Outer route choice.
  required: false
outputs:
- id: sample_config
  type: SimulationConfiguration
  units: n/a
  description: Each sample's run handed to the solve loop.
- id: policy_queries
  type: decision requests
  units: n/a
  description: Route selection and feedback.
- id: campaign_results
  type: MonteCarloResult
  units: n/a
  description: Ordered per-sample results with failures recorded.
tags:
- master-flow
charts:
- master
origin: agent
opens: simulation-simulation-campaigns
---

# Campaigns

## Purpose
Runs many simulations as one unit — Monte Carlo over seeds, constellation ensembles, adaptive routing studies — while recording every sample's outcome or failure and choosing serial, thread or process execution.

## Design & Implementation
`run_monte_carlo` collects seeds, validates the thread request, and dispatches through a shared job channel to thread workers or a `CachingPool` of Distributed workers, or lets the outer-route bandit choose in `:auto` mode; each sample calls the user closure, which normally builds a configuration and calls `run_simulation`. Results preserve seed order and `fail_fast` stops at the first failure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `campaign_spec` | MonteCarloSpec / closure | n/a | — | Seeds and the per-seed configuration builder. |
| in | `thread_policy` | PolicyDecision | n/a | no | Outer route choice. |
| out | `sample_config` | SimulationConfiguration | n/a | — | Each sample's run handed to the solve loop. |
| out | `policy_queries` | decision requests | n/a | — | Route selection and feedback. |
| out | `campaign_results` | MonteCarloResult | n/a | — | Ordered per-sample results with failures recorded. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.configure|Configure a run]] · `campaign_spec` → `campaign_spec` · dataflow · `src/simulation/campaigns/monte_carlo.jl`
- [[flow.parallel|Parallel routing & thread policy]] · `thread_policy` → `thread_policy` · dataflow · `src/simulation/campaigns/adaptive_routing.jl`

**Downstream**

- `policy_queries` → [[flow.parallel|Parallel routing & thread policy]] · `policy_queries` · feedback · `src/simulation/campaigns/adaptive_routing.jl`
- `sample_config` → [[flow.solve_loop|Solve loop]] · `sample_config` · dataflow · `src/simulation/campaigns/monte_carlo.jl`
<!-- vulcan:connections:end -->

## Limitations
Returning full solutions from the per-seed closure is expensive on the process backend because every result crosses a serialisation boundary; nested campaigns yield to the enclosing thread split and record no routing feedback.
