---
id: simulation.vacuum_predicted_gram__vacuum_rk4_step
label: _vacuum_rk4_step
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _vacuum_rk4_step
  lines:
  - 81
  - 81
inputs:
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: dt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `_vacuum_rk4_step`.
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

# _vacuum_rk4_step

## Purpose
Advances the drag-free reference state by one fixed step of size `dt` using classical fourth-order Runge-Kutta under the two-body-plus-J2 acceleration, producing the knot positions at which the vacuum GRAM cache samples the density model.

## Theory & Math
For the state $(\mathbf{r}, \mathbf{v})$ with $\dot{\mathbf{r}} = \mathbf{v}$ and $\dot{\mathbf{v}} = \mathbf{a}(\mathbf{r})$, the stages are $k_1 = (\mathbf{v}, \mathbf{a}(\mathbf{r}))$, $k_2 = (\mathbf{v} + \tfrac{h}{2}k_{1v}, \mathbf{a}(\mathbf{r} + \tfrac{h}{2}k_{1r}))$, $k_3 = (\mathbf{v} + \tfrac{h}{2}k_{2v}, \mathbf{a}(\mathbf{r} + \tfrac{h}{2}k_{2r}))$, $k_4 = (\mathbf{v} + h k_{3v}, \mathbf{a}(\mathbf{r} + h k_{3r}))$, and the update is $(\mathbf{r}, \mathbf{v})_{n+1} = (\mathbf{r}, \mathbf{v})_n + \tfrac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$ with $h$ = `dt` in seconds.

## Design & Implementation
Signature `_vacuum_rk4_step(pos, vel::SVector{3,Float64}, planet, dt::Float64)::Tuple{SVector{3,Float64}, SVector{3,Float64}}`, `@inline`. It evaluates the four RK4 stages on the coupled first-order system `dr/dt = v`, `dv/dt = a(r)`, calling `_vacuum_j2_accel` four times (at `pos`, `pos + dt/2 k1r`, `pos + dt/2 k2r`, `pos + dt k3r`), and combines them with weights `(1, 2, 2, 1) dt/6`. Both updated vectors are returned as a tuple. `_build_vacuum_gram_cache!` calls it `n_pts - 1` times with `dt = horizon_s / (n_pts - 1)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `dt` | Float64 | n/a | yes | Positional argument `dt`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `_vacuum_rk4_step`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.vacuum_predicted_gram_build_vacuum_gram_cache__build_vacuum_gram_cache_bang|_build_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:224-224`

**Downstream**

- `callees` → [[simulation.vacuum_predicted_gram__vacuum_j2_accel|_vacuum_j2_accel]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:88-88`
<!-- vulcan:connections:end -->

## Limitations
The step size is whatever the knot spacing happens to be (default 600 s / 19 ≈ 31.6 s), with no error control; for low perigee passes or large `SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S` with few knots, RK4 local error grows as `dt^5` and the predicted positions may miss the deviation threshold for the wrong reason. There is no energy or symplectic correction, so long horizons accumulate secular drift. The function ignores atmospheric drag by design, which is exactly what makes the prediction diverge inside the atmosphere.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 81.
