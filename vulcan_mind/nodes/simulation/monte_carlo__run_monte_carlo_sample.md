---
id: simulation.monte_carlo__run_monte_carlo_sample
label: _run_monte_carlo_sample
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: _run_monte_carlo_sample
  lines:
  - 76
  - 76
inputs:
- id: f
  type: Any
  units: n/a
  required: true
  description: Positional argument `f`.
- id: index
  type: Int
  units: n/a
  required: true
  description: Positional argument `index`.
- id: seed
  type: Any
  units: n/a
  required: true
  description: Positional argument `seed`.
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
  description: Return value of `_run_monte_carlo_sample`. Returns `MonteCarloSampleResult(`.
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

# _run_monte_carlo_sample

## Purpose
Runs one Monte Carlo sample under full exception capture so a single failing seed is recorded rather than aborting the campaign.

## Design & Implementation
Takes a wall-clock stamp with `time_ns`, calls `f(seed)`, and on success builds a `MonteCarloSampleResult` with `success=true`, the elapsed seconds and the returned value. On any thrown exception it captures `stacktrace(catch_backtrace())` — which must happen inside the `catch` block before anything else can throw — records elapsed time to the failure and returns a result with `success=false`, the error and the trace. Every path returns a result; nothing propagates.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Any | n/a | yes | Positional argument `f`. |
| in | `index` | Int | n/a | yes | Positional argument `index`. |
| in | `seed` | Any | n/a | yes | Positional argument `seed`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MonteCarloSampleResult | n/a | — | Return value of `_run_monte_carlo_sample`. Returns `MonteCarloSampleResult(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/monte_carlo.jl`
- [[simulation.monte_carlo__run_monte_carlo_process|_run_monte_carlo_process]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:176-176`
- [[simulation.monte_carlo__run_monte_carlo_serial|_run_monte_carlo_serial]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:112-112`
- [[simulation.monte_carlo__run_monte_carlo_threaded|_run_monte_carlo_threaded]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:136-136`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:79-79`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:79-79`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:79-79`
- `callees` → [[simulation.monte_carlo_montecarlosampleresult|MonteCarloSampleResult]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:81-81`
<!-- vulcan:connections:end -->

## Limitations
Catching every exception includes `InterruptException`, so a Ctrl-C during a sample is recorded as a failed sample instead of stopping the campaign; the timing includes exception unwinding, so failed samples read slightly longer than the useful work they did.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl` line 76.
