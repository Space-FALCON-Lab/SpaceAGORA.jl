---
id: gncz.hypr_utils_hypr_sample_count_path
label: hypr_sample_count_path
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_sample_count_path
  lines:
  - 49
  - 77
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: HYPRUtils module namespace exporting the shared path sampling and swarm
    scheduling helpers.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: samples
  type: Matrix{Float64}
  units: mixed
  description: Dense matrix of sampled path states with one column per requested sample
    and endpoints pinned to the control points.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# hypr_sample_count_path

## Purpose
`hypr_sample_count_path` densifies a small set of control points into a fixed-length sample matrix so that path cost, clearance, and smoothness can be evaluated on a uniform grid. Every hybrid particle swarm and refinement path planner in this codebase scores candidates through this sampling step, which makes the sample count a direct cost-accuracy control.

## Theory & Math
In Bezier mode the sample at normalised parameter $s$ is the Bernstein combination $P(s) = \sum_{i=0}^{d} \binom{d}{i}(1-s)^{d-i}s^i\,p_i$ over the $d+1$ control points, so the curve stays inside the convex hull of the control polygon. In polyline mode the parameter is mapped to $u = (n-1)s$, the segment index is the floor of $u$, and the sample is the linear blend of the two bracketing control points with weight $\alpha = u - \lfloor u \rfloor$.

## Model & Assumptions
The control points are supplied as a matrix with one column per point and any number of rows, so the same routine serves Cartesian paths and joint-space paths. At least two samples are required and the curve type must be one of the two supported symbols; both violations raise an argument error rather than degrading silently.

## Design & Implementation
The result is preallocated and filled column by column. After sampling, the first and last columns are overwritten with the first and last control points, which removes floating-point drift at the endpoints and guarantees that a planner sees exactly the start and goal states it asked for.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | HYPRUtils module namespace exporting the shared path sampling and swarm scheduling helpers. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `samples` | Matrix{Float64} | mixed | — | Dense matrix of sampled path states with one column per requested sample and endpoints pinned to the control points. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.swarm_and_retiming_robot_arm_sample_hypr_path|robot_arm_sample_hypr_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:113-113`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Uniform parameter spacing is not uniform arc-length spacing, so a Bezier with clustered control points is sampled unevenly in space. Binomial coefficients are recomputed inside the loop and the Bernstein form loses conditioning at high control point counts.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl:2-77`.
