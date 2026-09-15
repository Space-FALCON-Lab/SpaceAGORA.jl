---
id: simulation.vacuum_predicted_gram__natural_cubic_spline_build_bang
label: _natural_cubic_spline_build!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _natural_cubic_spline_build!
  lines:
  - 103
  - 103
inputs:
- id: Ms
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `Ms`.
- id: ys
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `ys`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
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
  type: Nothing
  units: n/a
  description: Return value of `_natural_cubic_spline_build!`; mutates `Ms` in place.
    Returns `nothing`.
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

# _natural_cubic_spline_build!

## Purpose
Computes the second-derivative coefficients of a natural cubic spline through uniformly spaced samples, storing them into a caller-owned vector. The vacuum GRAM cache runs it twice per rebuild, once on `log(ρ)` and once on temperature, so later queries can evaluate a smooth interpolant with `_eval_natural_cubic_spline`.

## Theory & Math
For knots $t_i = t_0 + (i-1)h$ and values $y_i$, the natural cubic spline second derivatives $M_i$ satisfy $M_1 = M_n = 0$ and, for $i = 2,\dots,n-1$,
$$M_{i-1} + 4M_i + M_{i+1} = \frac{6}{h^2}\left(y_{i+1} - 2y_i + y_{i-1}\right).$$
The Thomas algorithm solves this tridiagonal system in $O(n)$ with forward coefficients $c'_1 = 1/4$, $c'_i = 1/(4 - c'_{i-1})$, $d'_i = (r_i - d'_{i-1})/(4 - c'_{i-1})$ and back substitution $M_i = d'_i - c'_i M_{i+1}$.

## Design & Implementation
Signature `_natural_cubic_spline_build!(Ms::Vector{Float64}, ys::Vector{Float64}, h::Float64)`. It `resize!`s `Ms` to `n = length(ys)`, fills it with zeros (which already encodes the natural boundary conditions `M[1] = M[n] = 0`), and returns early when `n <= 2`. For the `m = n - 2` interior unknowns it solves the tridiagonal system with constant diagonals `[1 4 1]` and right-hand side `6/h^2 * (y[i+2] - 2 y[i+1] + y[i])` using the Thomas algorithm: a forward sweep computing modified coefficients `cp[i] = 1/(4 - cp[i-1])` and `dp[i]`, then back substitution writing `Ms[2..n-1]`. Two temporary vectors of length `m` are allocated per call. Returns `nothing`; `Ms` is the only mutation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `Ms` | Vector{Float64} | n/a | yes | Positional argument `Ms`. |
| in | `ys` | Vector{Float64} | n/a | yes | Positional argument `ys`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_natural_cubic_spline_build!`; mutates `Ms` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.vacuum_predicted_gram_build_vacuum_gram_cache__build_vacuum_gram_cache_bang|_build_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:228-228`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Uniform spacing is assumed; the function has no way to accept non-uniform knots. The two `Vector{Float64}(undef, m)` allocations happen on every rebuild even though `Ms` itself is reused, contradicting the allocation-free intent of the cache. A non-positive `h` produces `Inf`/`NaN` coefficients without an error. Natural boundary conditions force zero curvature at the horizon ends, which can visibly distort `log(ρ)` near atmospheric entry where curvature is largest.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 103.
