---
id: gnc.swarm_and_retiming__robot_arm_path_smoothness
label: _robot_arm_path_smoothness
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_path_smoothness
  lines:
  - 122
  - 122
inputs:
- id: samples
  type: Any
  units: n/a
  required: true
  description: Positional argument `samples`.
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
  description: Return value of `_robot_arm_path_smoothness`. Returns `total / max(size(samples,
    2) - 2, 1)`.
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

# _robot_arm_path_smoothness

## Purpose
Measures joint-space path roughness as the mean squared second finite difference of the sampled path, penalising sharp direction changes in the HYPR cost.

## Theory & Math
$$S = \frac{1}{n-2}\sum_{k=2}^{n-1} \left\| \mathbf{q}_{k+1} - 2\mathbf{q}_k + \mathbf{q}_{k-1} \right\|_2^2$$ where $\mathbf{q}_k$ is column $k$ of `samples` and $n$ is the number of columns.

## Design & Implementation
Returns 0.0 when `size(samples, 2) < 3`. Otherwise it loops `k` from 2 to `n - 1` under `@inbounds`, forms the second difference `samples[:, k+1] - 2 samples[:, k] + samples[:, k-1]`, accumulates `sum(abs2, ...)`, and divides by `max(n - 2, 1)`. The result is unnormalised by sample spacing, so it scales with the square of the inverse sample count for a fixed geometric path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_path_smoothness`. Returns `total / max(size(samples, 2) - 2, 1)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_robot_arm_hypr_path_cost_components|robot_arm_hypr_path_cost_components]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:14-14`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Each iteration allocates three temporary vectors from column slicing. Because the metric is not divided by the sample spacing squared, changing `n_samples` changes the cost weighting relative to path length.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 122.
