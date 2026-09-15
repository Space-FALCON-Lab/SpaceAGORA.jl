---
id: simulation.solver_policy__solve_with_gravity_backbone_solver
label: _solve_with_gravity_backbone_solver
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _solve_with_gravity_backbone_solver
  lines:
  - 548
  - 548
inputs:
- id: prob
  type: Any
  units: n/a
  required: true
  description: Positional argument `prob`.
- id: cfg
  type: SolverConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_solve_with_gravity_backbone_solver`. Returns `sol,
    solver_label`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _solve_with_gravity_backbone_solver

## Purpose
Integrates a partitioned translational problem with a fixed-step `KahanLi8` symplectic core for static gravity, wrapped by explicit half velocity-kicks for perturbations, assembling a stitched `ODESolution`.

## Theory & Math
Per macro step of length $h$: $v \leftarrow v + \tfrac{h}{2} a_{pert}(r, v, t)$, then $(r, v) \leftarrow \text{KahanLi8}_h(r, v)$ under the static gravity core, then $v \leftarrow v + \tfrac{h}{2} a_{pert}(r, v, t+h)$; $a_{pert}$ is the summed explicit perturbing acceleration.

## Design & Implementation
For an empty span it returns a one-point solution via `DiffEqBase.build_solution` with `ReturnCode.Success`. Otherwise `dt_s = _gravity_backbone_fixed_dt_s(cfg, args)`, and each macro step does: `_gravity_backbone_half_kick!(u_cursor, prob.p, t_cursor, half_dt)`; `remake(prob; u0=u_cursor, tspan=(t_cursor, t_next))` solved with `_solve_with_fixed_step_solver(core_prob, cfg, KahanLi8(), segment_dt)` (one symplectic step); the second half kick at `reached_t`; then the state is appended to `solution_ts`/`solution_us` (or the last entry overwritten if time did not advance). The loop breaks on unsuccessful retcode, on a retcode whose string is not `"Success"`, on non-finite segment, or when `_gravity_backbone_time_reached` fails. Returns `(build_solution(prob, alg, ts, us; retcode), label)` where the label notes `+Kicks` when `_gravity_backbone_has_kicks`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `prob` | Any | n/a | yes | Positional argument `prob`. |
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_solve_with_gravity_backbone_solver`. Returns `sol, solver_label`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:657-657`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/solver_policy.jl:549-549`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/solver_policy.jl:592-592`
- `callees` → [[simulation.dynamics_rhs__gravity_backbone_half_kick_bang|_gravity_backbone_half_kick!]] · `callers` · call · `src/simulation/engine/solver_policy.jl:581-581`
- `callees` → [[simulation.solver_policy__gravity_backbone_fixed_dt_s|_gravity_backbone_fixed_dt_s]] · `callers` · call · `src/simulation/engine/solver_policy.jl:565-565`
- `callees` → [[simulation.solver_policy__gravity_backbone_has_kicks|_gravity_backbone_has_kicks]] · `callers` · call · `src/simulation/engine/solver_policy.jl:552-552`
- `callees` → [[simulation.solver_policy__gravity_backbone_time_reached|_gravity_backbone_time_reached]] · `callers` · call · `src/simulation/engine/solver_policy.jl:618-618`
- `callees` → [[simulation.solver_policy__solve_with_fixed_step_solver|_solve_with_fixed_step_solver]] · `callers` · call · `src/simulation/engine/solver_policy.jl:585-585`
<!-- vulcan:connections:end -->

## Limitations
Each macro step calls `remake` and a full `solve` for a single `KahanLi8` step, which is far more overhead than a hand-rolled loop. Two `deepcopy` calls per step plus the growing `solution_us` vector store every macro state regardless of `needs_full_solution`. The half kicks use the same `dt_s` as the core, so perturbations are only second-order accurate. The `string(retcode) != "Success"` check duplicates `successful_retcode` and allocates a string per step.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 548.
