---
id: gnc.heat_load_control__edg_heat_load_lambdas
label: _edg_heat_load_lambdas
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_heat_load_lambdas
  lines:
  - 260
  - 260
inputs:
- id: track
  type: Any
  units: n/a
  required: true
  description: Positional argument `track`.
- id: alpha_profile
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `alpha_profile`.
- id: coeffs
  type: Any
  units: n/a
  required: true
  description: Positional argument `coeffs`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: area
  type: Float64
  units: n/a
  required: true
  description: Positional argument `area`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: scale_height
  type: Float64
  units: n/a
  required: true
  description: Positional argument `scale_height`.
- id: k
  type: Float64
  units: n/a
  required: true
  description: Positional argument `k`.
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
  description: Return value of `_edg_heat_load_lambdas`. Returns `lambda_switch, lambda_v`.
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

# _edg_heat_load_lambdas

## Purpose
Integrates the adjoint (costate) equations of the minimum-heat-load optimal control problem backwards along a track and returns the switching threshold and velocity costate that determine the bang-bang angle-of-attack profile.

## Theory & Math
$$\dot\lambda_v = -\frac{3k\rho v^2\alpha}{\pi} + \lambda_v\frac{\rho A C_D v}{m} - \lambda_\gamma\left(\frac{\rho A C_{L}}{2m} + \frac{g}{v^2} + \frac{1}{r}\right) - \lambda_h\gamma,\qquad \dot\lambda_\gamma = \lambda_v g - \lambda_h v$$ $$\dot\lambda_h = \frac{k\rho v^3\alpha}{\pi H} - \lambda_v\left(\frac{\rho A C_D v^2}{2mH} + \frac{2g\gamma}{r}\right) + \lambda_\gamma\left(\frac{\rho A C_L v}{2mH} - \frac{2g}{rv} + \frac{v}{r^2}\right),\qquad \lambda_{\mathrm{sw}} = \frac{2 k m v}{\pi A\, \partial C_D/\partial\alpha}$$ with $H$ the scale height (m), $k$ the heating weight, $C_L$ = `cl_low`, and $r = R_e + h$.

## Design & Implementation
Signature `(track, alpha_profile, coeffs, mass, area, planet, scale_height, k)`. Terminal conditions are `lambda_v[end] = track.speed[end]` and `lambda_h[end] = μ / (Rp_e + h_end)^2`, with `lambda_gamma[end] = 0`. It then steps backward with explicit Euler (`lambda[jm] = lambda[j] - lambda_dot * dt`) over the 3-state planar entry dynamics in `(v, gamma, h)` with exponential atmosphere of scale height `scale_height`, gravity `g = g_ref Rp_e^2 / r^2`, drag `cd = cd_low + alpha cd_slope`, and heating weight `k`. Finally `lambda_switch[j] = 2 k m v_j / (A cd_slope pi)` is computed for every node. Returns `(lambda_switch, lambda_v)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `alpha_profile` | Vector{Float64} | n/a | yes | Positional argument `alpha_profile`. |
| in | `coeffs` | Any | n/a | yes | Positional argument `coeffs`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `area` | Float64 | n/a | yes | Positional argument `area`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `scale_height` | Float64 | n/a | yes | Positional argument `scale_height`. |
| in | `k` | Float64 | n/a | yes | Positional argument `k`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_heat_load_lambdas`. Returns `lambda_switch, lambda_v`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_heat_load_alpha_profile|_edg_heat_load_alpha_profile]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:308-308`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Explicit Euler with the caller's grid spacing (up to several seconds) can be unstable for stiff costates near periapsis. The dynamics are planar with small-angle approximations (`sin γ ≈ γ`) and assume a linear CD-versus-alpha law. The heating term uses the factor `alpha / pi` from a specific heat-rate law and is not general.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 260.
