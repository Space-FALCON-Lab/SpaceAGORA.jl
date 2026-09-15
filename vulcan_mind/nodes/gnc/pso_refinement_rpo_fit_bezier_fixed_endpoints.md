---
id: gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints
label: rpo_fit_bezier_fixed_endpoints
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_fit_bezier_fixed_endpoints
  lines:
  - 116
  - 116
inputs:
- id: samples
  type: Any
  units: n/a
  required: true
  description: Positional argument `samples`.
- id: n_control
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_control`.
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
  description: Return value of `rpo_fit_bezier_fixed_endpoints`. Returns `rpo_refinement_clamp_path(controls,
    cfg)`.
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

# rpo_fit_bezier_fixed_endpoints

## Purpose
Finds the Bezier control polygon of a given size that best reproduces a set of sampled points while pinning the first and last control points to the sample endpoints.

## Theory & Math
With samples $p_i$ at parameters $u_i$, degree $n$, fixed endpoints $c_0$ and $c_n$, and unknown internal controls $c_1 \ldots c_{n-1}$, the fit solves

$$
\min_{c_1..c_{n-1}} \sum_i \left\| p_i - B_{0,n}(u_i)\, c_0 - B_{n,n}(u_i)\, c_n - \sum_{j=1}^{n-1} B_{j,n}(u_i)\, c_j \right\|^2
$$

via $(A^\top A + \lambda I)\, C = A^\top R$ with ridge $\lambda$.

## Design & Implementation
Forces at least two control points and returns just the endpoints when there are no internal controls. Otherwise it builds the least-squares system: row `i` of `A` holds the Bernstein values for the internal indices at parameter `u[i]`, and the right-hand side is the sample minus the known endpoint contributions. The normal equations are formed, a ridge of `1e-10` times the largest diagonal entry (at least machine epsilon) is added for conditioning, and the system is solved with `\`. The fitted internal controls are clamped through `rpo_refinement_clamp_path`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `n_control` | Int | n/a | yes | Positional argument `n_control`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_fit_bezier_fixed_endpoints`. Returns `rpo_refinement_clamp_path(controls, cfg)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_seed_control_points|seed_control_points]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:305-305`
- [[gnc.pso_refinement_rpo_refine_lower_degree|rpo_refine_lower_degree]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:246-246`
- [[gnc.pso_refinement_rpo_refine_shortcut_refit|rpo_refine_shortcut_refit]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:185-185`
- [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:453-453`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:305-305`

**Downstream**

- `callees` → [[gnc.pso_refinement_rpo_refinement_bernstein|rpo_refinement_bernstein]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:132-132`
- `callees` → [[gnc.pso_refinement_rpo_refinement_clamp_path|rpo_refinement_clamp_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:148-148`
- `callees` → [[gnc.pso_refinement_rpo_refinement_sample_params|rpo_refinement_sample_params]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:128-128`
<!-- vulcan:connections:end -->

## Limitations
Forming `A' * A` squares the condition number, which the ridge only partially offsets; with many samples and few controls this is fine, but fitting a high-degree polygon to few samples can produce wildly oscillating handles that the clamp then flattens onto the box.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 116.
