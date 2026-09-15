---
id: gnc.pso_refinement_rpo_refinement_bernstein
label: rpo_refinement_bernstein
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_bernstein
  lines:
  - 98
  - 98
inputs:
- id: n
  type: Int
  units: n/a
  required: true
  description: Positional argument `n`.
- id: j
  type: Int
  units: n/a
  required: true
  description: Positional argument `j`.
- id: u
  type: Float64
  units: n/a
  required: true
  description: Positional argument `u`.
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
  description: Return value of `rpo_refinement_bernstein`. Returns `binomial(n, j)
    * (1.0 - u)^(n - j) * u^j`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_refinement_bernstein

## Purpose
Evaluates one Bernstein basis polynomial, the building block of the least-squares Bezier fit.

## Theory & Math
$$
B_{j,n}(u) = \binom{n}{j} (1-u)^{\,n-j} u^{\,j},\qquad u \in [0, 1]
$$

## Design & Implementation
Computes `binomial(n, j) * (1-u)^(n-j) * u^j` directly. The binomial coefficient comes from `Base.binomial` on integers and the powers are floating-point, so the value is exact in structure and accurate for the small degrees the planner uses.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n` | Int | n/a | yes | Positional argument `n`. |
| in | `j` | Int | n/a | yes | Positional argument `j`. |
| in | `u` | Float64 | n/a | yes | Positional argument `u`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_bernstein`. Returns `binomial(n, j) * (1.0 - u)^(n - j) * u^j`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints|rpo_fit_bezier_fixed_endpoints]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:132-132`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`Base.binomial` overflows `Int` for degrees above about sixty, and the direct power form loses precision near `u = 0` or `u = 1` for high `j` or `n - j`; neither matters for control polygons of a dozen points but the function is not general-purpose.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 98.
