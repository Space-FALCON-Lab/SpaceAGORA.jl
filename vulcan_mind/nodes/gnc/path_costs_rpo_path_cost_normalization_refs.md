---
id: gnc.path_costs_rpo_path_cost_normalization_refs
label: rpo_path_cost_normalization_refs
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_costs.jl
  symbol: rpo_path_cost_normalization_refs
  lines:
  - 2
  - 2
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
  type: Tuple
  units: n/a
  description: Return value of `rpo_path_cost_normalization_refs`. Returns `(straight_len=straight_len,
    len_ref=len_ref, fuel_ref=max(fuel_ref, 1.0e-12))`.
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

# rpo_path_cost_normalization_refs

## Purpose
`rpo_path_cost_normalization_refs(points, cfg::RPOPSOConfig)` computes the length and fuel scale factors that make the RPO objective dimensionless. Without them, the length term (metres) and the fuel proxy (kilograms) would be summed with hand-tuned weights that only work for one problem size; normalising both by a problem-derived reference lets `cfg.w_len`, `cfg.w_obs` and `cfg.w_fuel` stay meaningful across widely different standoff distances.

## Theory & Math
With chord length $L_0 = \lVert \mathbf{p}_N - \mathbf{p}_1 \rVert$ and transfer time $t_f$ = `cfg.tf_s`, the references are

$$L_{\mathrm{ref}} = \max\!\left(\begin{cases} d_{\mathrm{ref}} & d_{\mathrm{ref}} > 0\\ L_0 & \text{otherwise}\end{cases},\; \Delta s,\; 10^{-6}\right), \qquad v_{\mathrm{ref}} = \frac{L_{\mathrm{ref}}}{\max(t_f, 10^{-6})},$$

$$m_{\mathrm{ref}} = \max\!\left(\frac{m\, v_{\mathrm{ref}}}{\max(I_{sp} g_0,\, 10^{-9})},\; 10^{-12}\right),$$

where $d_{\mathrm{ref}}$ = `cfg.cost_ref_distance_m` (m), $\Delta s$ = `cfg.sample_ds_m` (m), $m$ = `cfg.mass_kg` (kg), $I_{sp}$ = `cfg.isp_s` (s) and $g_0$ = `cfg.g0_mps2` (m/s^2). The fuel reference is the propellant mass a single $v_{\mathrm{ref}}$ impulse would consume under the Tsiolkovsky small-burn limit $\Delta m \approx m\,\Delta v/(I_{sp} g_0)$.

## Design & Implementation
The control points are materialised as `pts = Matrix{Float64}(points)` and the straight-line chord is taken between the first and last columns, `straight_len = sqrt(dx^2 + dy^2 + dz^2)` in metres. `len_ref` prefers an explicit `cfg.cost_ref_distance_m` when it is positive and otherwise uses the chord, then is floored by `max(len_ref, cfg.sample_ds_m, 1.0e-6)` so it can never be smaller than one sampling step. A reference speed `v_ref = len_ref / max(cfg.tf_s, 1e-6)` follows, and the fuel reference is the rocket-equation mass flow `cfg.mass_kg * v_ref / max(cfg.isp_s * cfg.g0_mps2, 1e-9)`. The return is the NamedTuple `(straight_len, len_ref, fuel_ref)` with `fuel_ref` floored at `1.0e-12`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `rpo_path_cost_normalization_refs`. Returns `(straight_len=straight_len, len_ref=len_ref, fuel_ref=max(fuel_ref, 1.0e-12))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:125-125`
- [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:191-191`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The chord uses only the first and last control points, so a path whose endpoints nearly coincide — a loop or an inspection circumnavigation — gets a reference length near zero and is rescued only by the `sample_ds_m` floor, which inflates the normalised length term dramatically. The `max(..., 1e-9)` and `max(..., 1e-12)` floors prevent division by zero but silently mask a misconfigured zero `isp_s` or `mass_kg` instead of reporting it. `Matrix{Float64}(points)` copies the whole control-point array on every call, and the function assumes at least three rows and at least one column, throwing `BoundsError` otherwise. `NaN` inputs propagate through `max` unpredictably.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_costs.jl` line 2.
