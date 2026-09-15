---
id: simulation.vacuum_predicted_gram__vacuum_j2_accel
label: _vacuum_j2_accel
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _vacuum_j2_accel
  lines:
  - 63
  - 63
inputs:
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_vacuum_j2_accel`.
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

# _vacuum_j2_accel

## Purpose
Evaluates the inertial-frame gravitational acceleration of a point mass plus the J2 oblateness term, used to propagate the drag-free reference trajectory that the vacuum GRAM cache samples along. It is inlined here rather than imported so the density-callback module does not depend on a private symbol of the gravity effectors.

## Theory & Math
With $r = \|\mathbf{r}\|$, $\mathbf{r} = (x, y, z)$ in the planet-centred inertial frame whose $z$ axis is the spin axis, $\mu$ the gravitational parameter, $J_2$ the second zonal harmonic and $R_e$ the equatorial radius:
$$\mathbf{a} = -\frac{\mu}{r^3}\mathbf{r} + \frac{3 J_2 \mu R_e^2}{2 r^4}\begin{pmatrix} \frac{x}{r}\left(5\frac{z^2}{r^2} - 1\right) \\ \frac{y}{r}\left(5\frac{z^2}{r^2} - 1\right) \\ \frac{z}{r}\left(5\frac{z^2}{r^2} - 3\right)\end{pmatrix}.$$

## Design & Implementation
Signature `_vacuum_j2_accel(pos::SVector{3,Float64}, planet)::SVector{3,Float64}`, `@inline`. It reads `μ = planet.μ`, `J2 = planet.J2` and the equatorial radius `Re = planet.Rp_e` (converted to `Float64`), computes `r = norm(pos)`, the spherical term `a_sph = (-μ/r^2) * normalize(pos)`, and a scale `1.5 J2 μ Re^2 / r^4` multiplying the vector `(x/r (5 z^2/r^2 - 1), y/r (5 z^2/r^2 - 1), z/r (5 z^2/r^2 - 3))`. The sum is returned as a new `SVector`; no allocation occurs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_vacuum_j2_accel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- [[simulation.vacuum_predicted_gram__vacuum_rk4_step|_vacuum_rk4_step]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:88-88`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:65-65`
<!-- vulcan:connections:end -->

## Limitations
The J2 term assumes the inertial `z` axis coincides with the planet's rotation axis, which holds for the planet-centred frames used by the engine but not for an arbitrary inertial frame. Higher zonal and tesseral harmonics, third-body and SRP effects present in the full dynamics are omitted, so over long horizons the vacuum prediction drifts from the true trajectory and forces cache rebuilds. `r = 0` divides by zero. The `planet` argument is untyped and must expose `μ`, `J2` and `Rp_e` fields.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 63.
