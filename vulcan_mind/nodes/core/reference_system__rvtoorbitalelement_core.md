---
id: core.reference_system__rvtoorbitalelement_core
label: _rvtoorbitalelement_core
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _rvtoorbitalelement_core
  lines:
  - 159
  - 159
inputs:
- id: r
  type: SVector
  units: n/a
  required: true
  description: Positional argument `r`.
- id: v
  type: SVector
  units: n/a
  required: true
  description: Positional argument `v`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  description: Return value of `_rvtoorbitalelement_core`. Returns `a, e, i, Ω, ω,
    ν`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# _rvtoorbitalelement_core

## Purpose
Converts an inertial Cartesian state to classical orbital elements robustly, handling the circular and equatorial singular cases by convention.

## Theory & Math
$$
\epsilon = \frac{v^2}{2} - \frac{\mu}{r},\quad a = -\frac{\mu}{2\epsilon},\quad \vec{e} = \frac{\vec{v} \times \vec{h}}{\mu} - \frac{\vec{r}}{r},\quad \cos i = \frac{h_z}{|\vec{h}|}
$$

## Design & Implementation
Computes specific energy for `a`, the eccentricity vector, angular momentum `h` and node vector `n = k × h`. Inclination comes from `h_z / |h|` through `_safe_acos`. With tolerances of 1e-12 on `e` and `|n|`, it sets `Ω` to zero when equatorial, `ω` to zero when circular or to the angle of the eccentricity vector from x when equatorial-elliptical, and `ν` to true longitude or argument of latitude in the circular cases, with quadrant fixes from `e_z`, `r_z` and `r · v`. Returns the six elements as a tuple.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r` | SVector | n/a | yes | Positional argument `r`. |
| in | `v` | SVector | n/a | yes | Positional argument `v`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rvtoorbitalelement_core`. Returns `a, e, i, Ω, ω, ν`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:226-226`

**Downstream**

- `callees` → [[core.reference_system__safe_acos|_safe_acos]] · `callers` · call · `src/core/interfaces/reference_system.jl:179-179`
- `callees` → [[core.reference_system__wrap_2pi|_wrap_2pi]] · `callers` · call · `src/core/interfaces/reference_system.jl:189-189`
<!-- vulcan:connections:end -->

## Limitations
For parabolic orbits the energy is zero and `a` is infinite; the function does not guard against that. The 1e-12 tolerances are absolute, so a nearly circular orbit at `e = 1e-11` is treated as elliptical and its `ω` and `ν` become numerically unstable.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 159.
