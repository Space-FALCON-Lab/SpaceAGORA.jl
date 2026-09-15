---
id: gnc.trajectory_optimizers_rpo_chomp_obstacle_potential
label: rpo_chomp_obstacle_potential
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_chomp_obstacle_potential
  lines:
  - 157
  - 157
inputs:
- id: clearance
  type: Any
  units: n/a
  required: true
  description: Positional argument `clearance`.
- id: safe_distance
  type: Any
  units: n/a
  required: true
  description: Positional argument `safe_distance`.
- id: margin
  type: Any
  units: n/a
  required: true
  description: Positional argument `margin`.
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
  description: Return value of `rpo_chomp_obstacle_potential`. Returns `-d + 0.5 *
    eps` or `0.5 * (d - eps)^2 / eps` or `0.0`.
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

# rpo_chomp_obstacle_potential

## Purpose
Evaluates the CHOMP-style piecewise soft obstacle potential for one sample given its `clearance` to the station, the required `safe_distance`, and a smoothing `margin`, producing a cost that is linear inside the keep-out zone, quadratic within the margin band, and zero beyond it.

## Theory & Math
$$c(d) = \begin{cases} -d + \tfrac{1}{2}\varepsilon & d < 0 \\ \tfrac{1}{2}\,\dfrac{(d-\varepsilon)^2}{\varepsilon} & 0 \le d \le \varepsilon \\ 0 & d > \varepsilon \end{cases}$$ with $d = \text{clearance} - \text{safe\_distance}$ in metres and $\varepsilon = \max(\text{margin}, 10^{-6})$ the band width.

## Design & Implementation
Computes `d = clearance - safe_distance` and `eps = max(Float64(margin), 1.0e-6)`. If `d < 0` returns `-d + 0.5 * eps` (linear penetration cost with continuity at `d = 0`); if `0 <= d <= eps` returns `0.5 * (d - eps)^2 / eps`; otherwise `0.0`. The function is C1-continuous at both breakpoints, which keeps the finite-difference gradient well behaved.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `clearance` | Any | n/a | yes | Positional argument `clearance`. |
| in | `safe_distance` | Any | n/a | yes | Positional argument `safe_distance`. |
| in | `margin` | Any | n/a | yes | Positional argument `margin`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `rpo_chomp_obstacle_potential`. Returns `-d + 0.5 * eps` or `0.5 * (d - eps)^2 / eps` or `0.0`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_soft_obstacle_cost_from_samples|rpo_soft_obstacle_cost_from_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:176-176`
- [[gnc.trajectory_optimizers_rpo_stomp_waypoint_state_cost|rpo_stomp_waypoint_state_cost]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:349-349`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:159-159`
<!-- vulcan:connections:end -->

## Limitations
Units are metres; the potential has no upper bound, so a sample deep inside the station contributes a large linear cost that can dominate the objective. The `1.0e-6` floor on `eps` is hard-coded. `clearance` is assumed finite; an `Inf` clearance returns 0 but a `NaN` propagates.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 157.
