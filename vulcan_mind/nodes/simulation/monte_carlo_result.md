---
id: simulation.monte_carlo_result
label: MonteCarloResult
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: MonteCarloResult
  lines:
  - 51
  - 63
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: samples
  type: Vector
  units: n/a
  required: true
  description: Collected outputs from the configured sample function.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
tags:
- simulation
- monte-carlo
- result
charts:
- simulation
origin: agent
outputs:
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
---

# MonteCarloResult

## Purpose
`MonteCarloResult` is the return container for a repeated simulation campaign. It preserves collected sample outputs together with route, seed, timing, and configuration metadata needed to interpret the run and to compute downstream statistics without reopening worker state.

## Theory & Math
For scalar samples the result supports ordinary estimates such as the mean, variance, quantiles, and empirical failure rate. For structured simulation results, those statistics require a caller-defined projection. The container records observations; it does not assert that the empirical distribution is converged to a physical or mission-level truth.

## Model & Assumptions
Samples are assumed to be ordered consistently with the sampler or to carry their own identifiers when execution is parallel. Metadata must preserve the seed and route policy if a study is intended to be reproducible. Empty campaigns are valid only if downstream statistics define their behavior for zero observations.

## Design & Implementation
`MonteCarloResult` is declared near the specification types in `monte_carlo.jl`. `run_monte_carlo` constructs it after sample execution and attaches the metadata produced by route selection and campaign timing. Analysis modules can consume the sample vector without depending on the process pool or worker identifiers.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `samples` | Vector | n/a | yes | Collected outputs from the configured sample function. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.adaptive_routing__run_campaign_adaptive|_run_campaign_adaptive]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:241-241`
- [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:190-190`
- [[simulation.run_monte_carlo|run_monte_carlo]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:257-257`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:62-62`
<!-- vulcan:connections:end -->

## Limitations
The result does not automatically detect outliers, nonstationarity, or correlated samples. Partial campaign termination must be interpreted from status and metadata rather than from the sample vector alone. Large structured samples can retain substantial memory because the container keeps the collected outputs for later analysis.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl:42-63`.
