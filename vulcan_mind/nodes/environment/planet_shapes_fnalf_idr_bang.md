---
id: environment.planet_shapes_fnalf_idr_bang
label: fnALF_IDR!
kind: function
source:
  file: src/environment/ephemerides/planet_shapes.jl
  symbol: fnALF_IDR!
  lines:
  - 11
  - 11
inputs:
- id: A
  type: AbstractArray{Float64}
  units: n/a
  required: true
  description: Positional argument `A`.
- id: x
  type: Float64
  units: n/a
  required: true
  description: Positional argument `x`.
- id: N
  type: Integer
  units: n/a
  required: true
  description: Positional argument `N`.
- id: M
  type: Integer
  units: n/a
  required: true
  description: Positional argument `M`.
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
  description: Return value of `fnALF_IDR!`; mutates `A` in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# fnALF_IDR!

## Purpose

Fills a preallocated table with fully normalized associated Legendre functions (fnALFs) from degree 0 to N and order 0 to M evaluated at a single argument x, supplying the angular basis that the topography and gravity harmonic sums are evaluated against.

## Design & Implementation

Mutates the caller's array `A` in place; nothing is returned. It seeds `A[1,1] = 1.0`, computes `ξ = √(1 - x^2)`, then walks the sectorial diagonal up to `min(N, M)` with the normalised ratio `√(((1 + δ(1,n))(2n + 1)) / 2n)` applied to the previous diagonal entry. A second doubly nested loop sweeps each order column `j` upward in degree using the increasing-degree recursion, taking the two-term form for the first off-diagonal entry (`i == j + 1`, where the degree-minus-two term does not exist) and the three-term form thereafter. Both loops are `@inbounds`, which is why the array must already be at least `(N+1) × (M+1)`.

## Theory & Math

With $x = \sin\varphi$ and $\xi = \sqrt{1 - x^2}$, the sectorial seed is

$$\bar{P}_{nn} = \sqrt{\frac{(1 + \delta_{1n})(2n+1)}{2n}}\, \xi\, \bar{P}_{n-1,n-1}$$

and the increasing-degree recursion for $n > m$ is

$$\bar{P}_{nm} = g_{nm}\, x\, \bar{P}_{n-1,m} - h_{nm}\, \bar{P}_{n-2,m}$$

$$g_{nm} = \sqrt{\frac{(2n+1)(2n-1)}{(n+m)(n-m)}}, \qquad h_{nm} = \sqrt{\frac{(2n+1)(n-m-1)(n+m-1)}{(2n-3)(n+m)(n-m)}}$$

The normalisation keeps the coefficients $O(1)$ at high degree, where the unnormalised functions overflow.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `A` | AbstractArray{Float64} | n/a | yes | Positional argument `A`. |
| in | `x` | Float64 | n/a | yes | Positional argument `x`. |
| in | `N` | Integer | n/a | yes | Positional argument `N`. |
| in | `M` | Integer | n/a | yes | Positional argument `M`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `fnALF_IDR!`; mutates `A` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.env_planet_shapes_calculate_topography_harmonics_bang|calculate_topography_harmonics!]] · `callees` → `callers` · call · `src/environment/ephemerides/planet_shapes.jl:70-70`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Entries with order greater than degree are never written, so a reused buffer keeps stale values in that triangle and callers must not read it; the routine also never zeroes `A`, so any previously computed table is only partially overwritten. `@inbounds` removes the size check entirely: an undersized `A` corrupts memory rather than raising. Passing `|x| > 1` makes `ξ` a `NaN` that spreads through the whole sectorial diagonal, and the recursion assumes `N` and `M` are non-negative. Because the buffer is mutated, one `A` may not be shared between threads evaluating different latitudes.

## Provenance
Mapped from `src/environment/ephemerides/planet_shapes.jl` line 11.
