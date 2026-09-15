---
id: gnc.trajectory_optimizers_rpo_chomp_numeric_gradient
label: rpo_chomp_numeric_gradient
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_chomp_numeric_gradient
  lines:
  - 208
  - 208
inputs:
- id: theta
  type: Any
  units: n/a
  required: true
  description: Positional argument `theta`.
- id: start
  type: Any
  units: n/a
  required: true
  description: Positional argument `start`.
- id: goal
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal`.
- id: objective_fn
  type: Any
  units: n/a
  required: true
  description: Positional argument `objective_fn`.
- id: gradient_eps
  type: Any
  units: n/a
  required: true
  description: Positional argument `gradient_eps`.
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
  description: Return value of `rpo_chomp_numeric_gradient`. Returns `grad`.
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

# rpo_chomp_numeric_gradient

## Purpose
Estimates the gradient of the trajectory objective with respect to every internal waypoint coordinate by centred finite differences, because the soft objective's clearance queries are not analytically differentiable.

## Theory & Math
$$\frac{\partial J}{\partial \theta_{d,i}} \approx \frac{J(\theta + h\,e_{d,i}) - J(\theta - h\,e_{d,i})}{2h},\qquad h = \max\left(\varepsilon_g,\ 10^{-6}\max(|\theta_{d,i}|, 1)\right)$$ where $e_{d,i}$ is the unit perturbation of axis $d$ at waypoint $i$ and $\varepsilon_g$ is `gradient_eps`.

## Design & Implementation
Signature `rpo_chomp_numeric_gradient(theta, start, goal, objective_fn, gradient_eps)`. For each waypoint `i` and axis `d in 1:3` it chooses `h = max(gradient_eps, 1e-6 * max(abs(theta[d,i]), 1))`, perturbs `theta[d,i]` in place by `+h`, evaluates `objective_fn(rpo_trajectory_points_from_internal(theta, start, goal))`, then by `-2h`, evaluates again, restores by `+h`, and stores `(c_plus - c_minus) / (2h)`. Returns a `3 x n` gradient matrix; `theta` is left numerically restored.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `theta` | Any | n/a | yes | Positional argument `theta`. |
| in | `start` | Any | n/a | yes | Positional argument `start`. |
| in | `goal` | Any | n/a | yes | Positional argument `goal`. |
| in | `objective_fn` | Any | n/a | yes | Positional argument `objective_fn`. |
| in | `gradient_eps` | Any | n/a | yes | Positional argument `gradient_eps`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_chomp_numeric_gradient`. Returns `grad`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:273-273`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:212-212`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_points_from_internal|rpo_trajectory_points_from_internal]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:214-214`
<!-- vulcan:connections:end -->

## Limitations
Costs `6 n_waypoints` objective evaluations per call, each of which re-samples the path and queries station clearance. In-place perturbation of `theta` means the function is not safe to call concurrently on a shared matrix, and floating-point round-trip leaves `theta` restored only to machine precision. The piecewise obstacle potential has gradient discontinuities that a fixed `h` can straddle.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 208.
