---
id: gnc.trajectory_optimizers_rpo_second_difference_metric
label: rpo_second_difference_metric
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_second_difference_metric
  lines:
  - 122
  - 122
inputs:
- id: n_waypoints
  type: Integer
  units: n/a
  required: true
  description: Positional argument `n_waypoints`.
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
  description: Return value of `rpo_second_difference_metric`. Returns `R, R_inv,
    M`.
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

# rpo_second_difference_metric

## Purpose
Constructs the discrete second-difference (acceleration) metric used by CHOMP as a gradient preconditioner and by STOMP as the noise covariance and update smoother, following the standard STOMP/CHOMP formulation for `n_waypoints` internal points.

## Theory & Math
$$A_{ii} = -2,\quad A_{i,i\pm1} = 1,\qquad R = A^{\top}A + 10^{-8} I,\qquad M_{:,j} = \frac{R^{-1}_{:,j}}{n\,\max_i |R^{-1}_{ij}|}$$ where $n$ is the number of internal waypoints, $R$ is the smoothness metric whose quadratic form $\theta^{\top} R \theta$ sums squared second differences, $R^{-1}$ is the CHOMP preconditioner, and $M$ is the STOMP column-normalised projector.

## Design & Implementation
For `n_waypoints <= 0` returns three empty matrices. Otherwise builds tridiagonal `A` with `-2` on the diagonal and `1` on the off-diagonals (a finite-difference second-derivative operator with implicit zero boundary), forms `R = A' * A + 1.0e-8 * I`, inverts it via `inv(Symmetric(R))`, and derives `M` by scaling every column `j` of `R_inv` so its maximum absolute entry equals `1 / n_waypoints` (columns whose max is below `1.0e-12` are left untouched). Returns `(R, R_inv, M)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_waypoints` | Integer | n/a | yes | Positional argument `n_waypoints`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_second_difference_metric`. Returns `R, R_inv, M`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:376-376`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:244-244`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Dense `inv` is O(n^3) and allocates three `n x n` matrices; it is computed once per plan so cost is acceptable for the tens of waypoints typical here but would not scale to thousands. The `1.0e-8` Tikhonov term and `1.0e-12` column threshold are hard-coded. Boundary handling assumes fixed endpoints outside `theta`, matching `rpo_trajectory_points_from_internal`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 122.
