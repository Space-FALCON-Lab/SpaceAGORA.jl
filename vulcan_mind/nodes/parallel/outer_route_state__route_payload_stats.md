---
id: parallel.outer_route_state__route_payload_stats
label: _route_payload_stats
kind: function
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: _route_payload_stats
  lines:
  - 101
  - 101
inputs:
- id: payload
  type: Any
  units: n/a
  required: true
  description: Positional argument `payload`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_route_payload_stats`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _route_payload_stats

## Purpose
Parses one `stats` table from a persisted outer-route TOML file back into an `OuterRouteStats`, tolerating missing keys, wrong types and legacy files that predate second-moment data. Returns `nothing` when the payload is not a dict or has no samples.

## Theory & Math
The second moment is floored at the value implied by zero variance: $$S_2 \ge \frac{S_1^2}{n}$$ where $S_1$ = `elapsed_sum_s`, $S_2$ = `elapsed_sq_sum_s`, $n$ = `samples`, since $\operatorname{Var} = S_2/n - (S_1/n)^2 \ge 0$.

## Design & Implementation
Marked `@inline`. Each field is read with `get(payload, key, default)` inside a `try` block so a non-numeric value falls back to `0`, `0.0` or `NaN` instead of throwing. Consistency is then enforced: `successes = min(samples, successes)`, `failures = min(samples - successes, failures)`. If `elapsed_sq_sum_s` is not finite (legacy schema), it is reconstructed as `elapsed_sum_s^2 / samples`, which corresponds to zero variance; otherwise it is raised to at least that same lower bound so the implied variance is never negative. The result is built with the keyword constructor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `payload` | Any | n/a | yes | Positional argument `payload`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_route_payload_stats`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_state.jl`
- [[parallel.outer_route_state_load_outer_route_state_bang|load_outer_route_state!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:219-219`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:119-119`
- `callees` → [[parallel.outer_route_state_outerroutestats|OuterRouteStats]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:137-137`
<!-- vulcan:connections:end -->

## Limitations
Silent fallbacks mean a corrupted file loads as partially zeroed statistics with no warning. `Int(get(...))` of a `Float64` such as `3.5` throws `InexactError`, which is caught and replaced by `0`, discarding the row's counts while keeping its timings. A payload with `samples > 0` but `elapsed_sum_s = 0` is accepted and will bias adaptive selection toward that route with an apparent mean of zero seconds.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl` line 101.
