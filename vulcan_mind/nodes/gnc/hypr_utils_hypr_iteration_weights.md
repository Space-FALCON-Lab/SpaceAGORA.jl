---
id: gnc.hypr_utils_hypr_iteration_weights
label: hypr_iteration_weights
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_iteration_weights
  lines:
  - 80
  - 80
inputs:
- id: schedule_enable
  type: Bool
  units: n/a
  required: true
  description: Positional argument `schedule_enable`.
- id: n_iters
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_iters`.
- id: iter
  type: Int
  units: n/a
  required: true
  description: Positional argument `iter`.
- id: w_inertia
  type: Real
  units: n/a
  required: true
  description: Positional argument `w_inertia`.
- id: c1
  type: Real
  units: n/a
  required: true
  description: Positional argument `c1`.
- id: c2
  type: Real
  units: n/a
  required: true
  description: Positional argument `c2`.
- id: transition_fraction
  type: Real
  units: n/a
  required: true
  description: Positional argument `transition_fraction`.
- id: w_min
  type: Real
  units: n/a
  required: true
  description: Positional argument `w_min`.
- id: w_end_fraction
  type: Real
  units: n/a
  required: true
  description: Positional argument `w_end_fraction`.
- id: c1_end_fraction
  type: Real
  units: n/a
  required: true
  description: Positional argument `c1_end_fraction`.
- id: c2_end_fraction
  type: Real
  units: n/a
  required: true
  description: Positional argument `c2_end_fraction`.
- id: c_min
  type: Real
  units: n/a
  required: true
  description: Positional argument `c_min`.
- id: c_max
  type: Real
  units: n/a
  required: true
  description: Positional argument `c_max`.
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
  description: Return value of `hypr_iteration_weights`. Returns `(w_inertia=Float64(w_inertia),
    c1=Float64(c1), c2=Float64(c2))` or `(`.
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

# hypr_iteration_weights

## Purpose
Produces the PSO inertia and acceleration coefficients for one iteration, smoothly tapering from the configured start values toward end values over a configurable fraction of the run.

## Theory & Math
$$
p = \operatorname{clamp}\left(\frac{k - 1}{(N - 1)\,\tau}, 0, 1\right),\qquad s = p^2 (3 - 2p),\qquad w_k = w_0 + s\,(w_{\text{end}} - w_0)
$$

with $k$ the iteration, $N$ the total and $\tau$ the transition fraction; $c_1$ and $c_2$ follow the same form.

## Design & Implementation
Returns the raw values unchanged when scheduling is off or there is only one iteration. Otherwise it maps `(iter - 1)` onto progress in the unit interval over the first `transition_fraction` of the iterations, applies the smoothstep `p²(3 - 2p)`, and interpolates each weight toward its end value: inertia toward the larger of `w_min` and `w_end_fraction * w_inertia`, and each acceleration toward its end fraction clamped into `[c_min, c_max]`. Returns a named tuple `(w_inertia, c1, c2)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `schedule_enable` | Bool | n/a | yes | Positional argument `schedule_enable`. |
| in | `n_iters` | Int | n/a | yes | Positional argument `n_iters`. |
| in | `iter` | Int | n/a | yes | Positional argument `iter`. |
| in | `w_inertia` | Real | n/a | yes | Positional argument `w_inertia`. |
| in | `c1` | Real | n/a | yes | Positional argument `c1`. |
| in | `c2` | Real | n/a | yes | Positional argument `c2`. |
| in | `transition_fraction` | Real | n/a | yes | Positional argument `transition_fraction`. |
| in | `w_min` | Real | n/a | yes | Positional argument `w_min`. |
| in | `w_end_fraction` | Real | n/a | yes | Positional argument `w_end_fraction`. |
| in | `c1_end_fraction` | Real | n/a | yes | Positional argument `c1_end_fraction`. |
| in | `c2_end_fraction` | Real | n/a | yes | Positional argument `c2_end_fraction`. |
| in | `c_min` | Real | n/a | yes | Positional argument `c_min`. |
| in | `c_max` | Real | n/a | yes | Positional argument `c_max`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `hypr_iteration_weights`. Returns `(w_inertia=Float64(w_inertia), c1=Float64(c1), c2=Float64(c2))` or `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_rpo_pso_iteration_weights|rpo_pso_iteration_weights]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:3-3`
- [[gnc.swarm_and_retiming__robot_arm_hypr_iteration_weights|_robot_arm_hypr_iteration_weights]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:3-3`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:96-96`
<!-- vulcan:connections:end -->

## Limitations
After the transition point the weights are constant, so the schedule cannot express a late-run change; the clamp of `transition_fraction` to at least 1e-6 avoids division by zero but makes a near-zero fraction jump to the end values on iteration two.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 80.
