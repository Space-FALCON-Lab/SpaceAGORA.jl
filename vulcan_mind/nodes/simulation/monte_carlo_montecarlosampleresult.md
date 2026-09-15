---
id: simulation.monte_carlo_montecarlosampleresult
label: MonteCarloSampleResult
kind: struct
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: MonteCarloSampleResult
  lines:
  - 32
  - 32
inputs:
- id: index
  type: Int
  units: n/a
  required: true
  description: Field `index`.
- id: seed
  type: Any
  units: n/a
  required: true
  description: Field `seed`.
- id: success
  type: Bool
  units: n/a
  required: true
  description: Field `success`.
- id: elapsed_s
  type: Float64
  units: n/a
  required: true
  description: Field `elapsed_s`.
- id: value
  type: Any
  units: n/a
  required: false
  description: Field `value` (default `nothing`).
- id: error
  type: Any
  units: n/a
  required: false
  description: Field `error` (default `nothing`).
- id: backtrace
  type: Any
  units: n/a
  required: false
  description: Field `backtrace` (default `nothing`).
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
  type: MonteCarloSampleResult
  units: n/a
  description: Constructed `MonteCarloSampleResult` (keyword constructor via @kwdef).
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

# MonteCarloSampleResult

## Purpose
The per-sample record a Monte Carlo campaign returns, capturing either the user function's value or the exception that ended the sample.

## Design & Implementation
A `Base.@kwdef struct` with `index` and `seed` identifying the sample, a `success` flag, `elapsed_s` measured with `time_ns` around the user call, and three `Any`-typed slots `value`, `error` and `backtrace` that default to `nothing`. Exactly one of `value` or the `error`/`backtrace` pair is populated depending on `success`. Keeping the payload slots untyped means the same result type serves every campaign regardless of what `f(seed)` returns, at the cost of a dynamic field access.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `index` | Int | n/a | yes | Field `index`. |
| in | `seed` | Any | n/a | yes | Field `seed`. |
| in | `success` | Bool | n/a | yes | Field `success`. |
| in | `elapsed_s` | Float64 | n/a | yes | Field `elapsed_s`. |
| in | `value` | Any | n/a | no | Field `value` (default `nothing`). |
| in | `error` | Any | n/a | no | Field `error` (default `nothing`). |
| in | `backtrace` | Any | n/a | no | Field `backtrace` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MonteCarloSampleResult | n/a | — | Constructed `MonteCarloSampleResult` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/monte_carlo.jl`
- [[simulation.monte_carlo__run_monte_carlo_sample|_run_monte_carlo_sample]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:81-81`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Nothing enforces the either-or between `value` and `error`; a hand-constructed result can carry both. `seed` and `value` are `Any`, so a large returned solution object is held alive in the results vector for the whole campaign, and the type offers no way to drop it.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl` line 32.
