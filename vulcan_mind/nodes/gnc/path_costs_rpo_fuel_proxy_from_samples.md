---
id: gnc.path_costs_rpo_fuel_proxy_from_samples
label: rpo_fuel_proxy_from_samples
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_costs.jl
  symbol: rpo_fuel_proxy_from_samples
  lines:
  - 16
  - 16
inputs:
- id: samples
  type: Any
  units: n/a
  required: true
  description: Positional argument `samples`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `rpo_fuel_proxy_from_samples`. Returns `fuel`.
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

# rpo_fuel_proxy_from_samples

## Purpose
`rpo_fuel_proxy_from_samples(samples, cfg::RPOPSOConfig)` estimates the propellant mass a candidate RPO path would consume, without running a dynamics simulation. It lets the PSO objective prefer smooth, low-acceleration trajectories over jerky ones that happen to be short, using a second-difference acceleration estimate integrated through the rocket equation's small-burn limit.

## Theory & Math
For samples $\mathbf{p}_j \in \mathbb{R}^3$ uniformly spaced in time by $\Delta t = t_f/(N-1)$, the discrete acceleration is the centred second difference

$$\mathbf{a}_j = \frac{\mathbf{p}_{j+2} - 2\mathbf{p}_{j+1} + \mathbf{p}_j}{\Delta t^{2}},$$

and the proxy integrates the small-burn mass-flow law

$$m_{\mathrm{fuel}} \;=\; \sum_{j=1}^{N-2} \frac{m \, \lVert \mathbf{a}_j \rVert}{I_{sp} g_0}\, \Delta t \;\approx\; \frac{m}{I_{sp} g_0}\int_0^{t_f} \lVert \mathbf{a}(t)\rVert \, dt,$$

with $m$ = `cfg.mass_kg` (kg), $I_{sp}$ = `cfg.isp_s` (s), $g_0$ = `cfg.g0_mps2` (m/s^2). It is the $\Delta v$ integral $\int \lVert \mathbf{a} \rVert dt$ scaled by constant mass flow, so it ignores mass depletion during the transfer.

## Design & Implementation
The samples are copied into `pts = Matrix{Float64}(samples)`; paths with fewer than three columns return `0.0` immediately, since a second difference needs three points. The sample spacing in time is `dt = max(cfg.tf_s / max(size(pts,2) - 1, 1), 1.0e-6)`, i.e. the transfer duration divided uniformly over the samples. An `@inbounds` loop over `j in 1:(N-2)` forms the centred second difference of each Cartesian component, `(p[j+2] - 2p[j+1] + p[j]) / dt^2`, takes the acceleration magnitude, and accumulates `mass_kg * |a| / max(isp_s * g0_mps2, 1e-9) * dt`. The result is a mass in kilograms.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_fuel_proxy_from_samples`. Returns `fuel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:143-143`
- [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:193-193`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The sample points are spaced by arc length, not by time — `rpo_sample_path` uses a distance step `sample_ds_m` — yet `dt` assumes uniform time spacing, so the acceleration estimate is biased wherever the path's speed is non-uniform. The second difference is $O(\Delta t^2)$-accurate but amplifies sampling noise by $1/\Delta t^2$, which for a dense path makes the term dominated by discretisation rather than by physics. Constant vehicle mass is assumed throughout, so the proxy underestimates fuel for large $\Delta v$. The unconditional `Matrix{Float64}(samples)` copy is a per-particle, per-iteration allocation inside the innermost PSO loop.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_costs.jl` line 16.
