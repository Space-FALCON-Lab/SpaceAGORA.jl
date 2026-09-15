---
id: gnc.path_costs_rpo_obstacle_sigmoid_penalty
label: rpo_obstacle_sigmoid_penalty
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_costs.jl
  symbol: rpo_obstacle_sigmoid_penalty
  lines:
  - 75
  - 75
inputs:
- id: clearance
  type: Real
  units: n/a
  required: true
  description: Positional argument `clearance`.
- id: threshold
  type: Real
  units: n/a
  required: true
  description: Positional argument `threshold`.
- id: k
  type: Real
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
  type: Float64
  units: n/a
  description: Return value of `rpo_obstacle_sigmoid_penalty`. Returns `y / (1.0 +
    y)` or `1.0 / (1.0 + exp(x))`.
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

# rpo_obstacle_sigmoid_penalty

## Purpose
`rpo_obstacle_sigmoid_penalty(clearance, threshold, k)` converts a signed clearance distance into a smooth penalty in `[0, 1]` that is essentially one when the vehicle is inside the keep-out envelope and essentially zero when it is safely outside. Using a sigmoid rather than a hard indicator gives the PSO objective a usable gradient near the constraint boundary, so particles are pushed away from obstacles instead of merely being flagged as infeasible.

## Theory & Math
The penalty is the decreasing logistic

$$\beta(d) \;=\; \frac{1}{1 + \exp\!\big(k\,(d - d_{\mathrm{thr}})\big)},$$

where $d$ is the clearance in metres, $d_{\mathrm{thr}} = d_{\mathrm{safe}} - \epsilon_{\mathrm{tol}}$ the penalty threshold, and $k$ the sharpness in inverse metres. It satisfies $\beta(d_{\mathrm{thr}}) = 1/2$, $\beta \to 1$ as $d \to -\infty$, $\beta \to 0$ as $d \to +\infty$, and has derivative $\beta'(d) = -k\,\beta(1-\beta)$, maximal in magnitude ($-k/4$) exactly at the threshold. As $k \to \infty$ it converges pointwise to the indicator $\mathbb{1}[d < d_{\mathrm{thr}}]$.

## Design & Implementation
The function is `@inline` and takes `Real` arguments converted to `Float64`. It forms `x = k * (clearance - threshold)` and then evaluates the logistic $1/(1+e^{x})$ in whichever of two algebraically identical branches is numerically safe: for `x >= 0` it computes `y = exp(-x)` and returns `y / (1 + y)`; for `x < 0` it returns `1.0 / (1.0 + exp(x))`. Each branch only ever exponentiates a non-positive argument, so `exp` can underflow to zero but never overflow to `Inf`, and the result is monotonically decreasing in clearance. Callers pass `threshold = safe_distance_m - obstacle_sigmoid_tol_m` and `k = cfg.obstacle_sigmoid_k`, whose default in the RPO path is `1.0e6` per metre.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `clearance` | Real | n/a | yes | Positional argument `clearance`. |
| in | `threshold` | Real | n/a | yes | Positional argument `threshold`. |
| in | `k` | Real | n/a | yes | Positional argument `k`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `rpo_obstacle_sigmoid_penalty`. Returns `y / (1.0 + y)` or `1.0 / (1.0 + exp(x))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_costs_rpo_clearance_stats_from_samples|rpo_clearance_stats_from_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:50-50`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:76-76`
<!-- vulcan:connections:end -->

## Limitations
With the default `k = 1.0e6` the transition region is about four microns wide, so for all practical clearances the function has already saturated to exactly `0.0` or `1.0` in double precision and provides no gradient at all — the smoothness that motivates the sigmoid is only available if `obstacle_sigmoid_k` is reduced to the scale of the actual clearance tolerance. `k * (clearance - threshold)` can overflow to `Inf` when the clearance is large and `k` is large, in which case `exp(-Inf)` is `0.0` and the result is a correct `0.0`, but a `NaN` clearance propagates to a `NaN` penalty that then contaminates the accumulated obstacle score. Negative `k` silently inverts the penalty, rewarding collisions.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_costs.jl` line 75.
