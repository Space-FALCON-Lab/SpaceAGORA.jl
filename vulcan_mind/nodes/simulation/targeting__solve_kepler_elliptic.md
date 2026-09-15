---
id: simulation.targeting__solve_kepler_elliptic
label: _solve_kepler_elliptic
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _solve_kepler_elliptic
  lines:
  - 48
  - 48
inputs:
- id: M
  type: Float64
  units: n/a
  required: true
  description: Positional argument `M`.
- id: e
  type: Float64
  units: n/a
  required: true
  description: Positional argument `e`.
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
  type: Float64
  units: n/a
  description: Return value of `_solve_kepler_elliptic`.
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

# _solve_kepler_elliptic

## Purpose
Solves Kepler's equation `M = E - e sin E` for the eccentric anomaly `E` given mean anomaly `M` and eccentricity `e`, the core of the two-body propagation used to aim GRAM track-cache prefetches.

## Theory & Math
Newton iteration on $f(E) = E - e\sin E - M$: $E_{k+1} = E_k - \dfrac{E_k - e\sin E_k - M}{1 - e\cos E_k}$, starting from $E_0 = M$ if $e < 0.8$ else $E_0 = \pi$.

## Design & Implementation
Reduces `M` to `M2π = mod(M, 2pi)` and chooses the initial guess `E = M2π` for `e < 0.8` and `E = π` for higher eccentricity, a standard robust start. Up to 20 Newton-Raphson iterations compute `f = E - e sin E - M2π` and `fp = 1 - e cos E`; the loop breaks early if `|fp| < 1e-14` (degenerate derivative) or if the update `|dE| <= 1e-12`. The final `E` is returned unwrapped (in `[0, 2π)` plus any Newton overshoot). The function is `@inline` and `Float64`-only.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `M` | Float64 | n/a | yes | Positional argument `M`. |
| in | `e` | Float64 | n/a | yes | Positional argument `e`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_solve_kepler_elliptic`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:91-91`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No convergence check after 20 iterations; for `e` very close to 1 Newton can converge slowly or oscillate and the returned `E` may carry error well above 1e-12. The `|fp| < 1e-14` break returns the current estimate without warning. Eccentricities outside `[0, 1)` are not rejected here; callers gate on them. The returned `E` is not re-wrapped, so a small negative value is possible for `M2π` near 0.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 48.
