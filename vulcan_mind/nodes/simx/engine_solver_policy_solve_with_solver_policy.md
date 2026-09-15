---
id: simx.engine_solver_policy_solve_with_solver_policy
label: _solve_with_solver_policy
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _solve_with_solver_policy
  lines:
  - 631
  - 748
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: problem
  type: SciMLBase.AbstractODEProblem
  units: n/a
  required: true
  description: Typed ODE or partitioned SecondOrderODEProblem built for the current
    segment, carrying the Jacobian prototype when one was supplied.
- id: solver_config
  type: SolverConfig
  units: n/a
  required: true
  description: Solver configuration whose solver_mode selects the branch and whose
    fixed-step and switching parameters tune it.
- id: tolerances
  type: NTuple{2,Float64}
  units: n/a
  required: true
  description: Relative and absolute tolerances resolved by _build_solver_tolerances
    for the adaptive branches.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: solution
  type: SciMLBase.AbstractODESolution
  units: n/a
  description: Solution object for the segment as returned by the selected integrator.
- id: solver_meta
  type: NamedTuple
  units: n/a
  description: Trace record with solver, initial_solver, fallback_used and trigger_retcode,
    appended to the run's solver_trace.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# _solve_with_solver_policy

## Purpose
`_solve_with_solver_policy` is the dispatch point between the configured solver mode and an actual integrator call. It validates that the requested mode is legal for the problem it was handed, runs the matching solve, and returns the solution together with a metadata tuple that makes the choice auditable after the fact.

## Theory & Math
The branches implement genuinely different discretisations. `:symplectic` uses `KahanLi8`, an eighth-order composition method for separable Hamiltonians $H(\vec q,\vec p) = T(\vec p) + V(\vec q)$, whose fixed-step flow map conserves a shadow Hamiltonian and so avoids the secular energy drift of a non-symplectic scheme. `:gravity_backbone_split` splits the vector field as $f = f_{grav} + f_{pert}$ and advances the dominant gravitational part with a symplectic backbone. `:multirate` applies Strang splitting, whose second-order accuracy follows from the symmetric composition $e^{\frac{h}{2}A}e^{hB}e^{\frac{h}{2}A} = e^{h(A+B)} + \mathcal{O}(h^{3})$ per step, with independent slow and fast integrators. `:auto_stiff` wraps `AutoTsit5(Rodas5P)`, switching between the explicit fifth-order Tsitouras pair and the stiffly accurate Rosenbrock method when the estimated stiffness ratio crosses the internal threshold.

## Model & Assumptions
Mode eligibility is checked before any solve. `:symplectic` demands a single-spacecraft, translational-only inverse-squared gravity configuration with no control, guidance or navigation effectors, and additionally a partitioned `SecondOrderODEProblem`; because the typed `run_simulation` path still builds a first-order problem, the error message explicitly directs callers to `:tsit5` or `:auto_stiff`. `:gravity_backbone_split` runs `_gravity_backbone_reject_reason` and the same partitioned-problem check. Every failure is an `ArgumentError` naming the constraint.

## Design & Implementation
The linear solver choice is conditional and deliberate: `_rodas5p_alg()` returns `Rodas5P(autodiff=AutoFiniteDiff(), linsolve=KLUFactorization())` only when `prob.f.jac_prototype !== nothing`, because KLU requires a sparse W matrix, and falls back to dense LU otherwise so single-satellite runs do not get wrong Newton corrections from KLU applied to a dense matrix. In the `:auto_stiff` branch, `_auto_stiff_smooth_gravity_eligible` short-circuits to plain `Tsit5` for smooth gravity-only problems; otherwise the autoswitching algorithm is used and per-step storage is kept regardless of `needs_full_solution`, because `_auto_stiff_switched` has to inspect `sol.alg_choice` across saved steps to report whether a switch occurred. The final fallthrough is `Tsit5`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `problem` | SciMLBase.AbstractODEProblem | n/a | yes | Typed ODE or partitioned SecondOrderODEProblem built for the current segment, carrying the Jacobian prototype when one was supplied. |
| in | `solver_config` | SolverConfig | n/a | yes | Solver configuration whose solver_mode selects the branch and whose fixed-step and switching parameters tune it. |
| in | `tolerances` | NTuple{2,Float64} | n/a | yes | Relative and absolute tolerances resolved by _build_solver_tolerances for the adaptive branches. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `solution` | SciMLBase.AbstractODESolution | n/a | — | Solution object for the segment as returned by the selected integrator. |
| out | `solver_meta` | NamedTuple | n/a | — | Trace record with solver, initial_solver, fallback_used and trigger_retcode, appended to the run's solver_trace. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__rodas5p_alg|_rodas5p_alg]] · `callees` → `callers` · feedback · `src/simulation/engine/solver_policy.jl:750-750`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:342-342`

**Downstream**

- `callees` → [[simulation.solver_policy__auto_stiff_switched|_auto_stiff_switched]] · `callers` · call · `src/simulation/engine/solver_policy.jl:700-700`
- `callees` → [[simulation.solver_policy__gravity_backbone_reject_reason|_gravity_backbone_reject_reason]] · `callers` · call · `src/simulation/engine/solver_policy.jl:652-652`
- `callees` → [[simulation.solver_policy__is_partitioned_second_order_problem|_is_partitioned_second_order_problem]] · `callers` · call · `src/simulation/engine/solver_policy.jl:639-639`
- `callees` → [[simulation.solver_policy__rodas5p_alg|_rodas5p_alg]] · `callers` · call · `src/simulation/engine/solver_policy.jl:669-669`
- `callees` → [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callers` · call · `src/simulation/engine/solver_policy.jl:674-674`
- `callees` → [[simulation.solver_policy__solve_with_fixed_step_solver|_solve_with_fixed_step_solver]] · `callers` · call · `src/simulation/engine/solver_policy.jl:642-642`
- `callees` → [[simulation.solver_policy__solve_with_gravity_backbone_solver|_solve_with_gravity_backbone_solver]] · `callers` · call · `src/simulation/engine/solver_policy.jl:657-657`
- `callees` → [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callers` · call · `src/simulation/engine/solver_policy.jl:721-721`
- `callees` → [[simulation.solver_policy__split_imex_solver_spec|_split_imex_solver_spec]] · `callers` · call · `src/simulation/engine/solver_policy.jl:710-710`
- `callees` → [[simulation.solver_policy__symplectic_conservative_eligible|_symplectic_conservative_eligible]] · `callers` · call · `src/simulation/engine/solver_policy.jl:636-636`
- `callees` → [[simulation.solver_policy__symplectic_fixed_dt_s|_symplectic_fixed_dt_s]] · `callers` · call · `src/simulation/engine/solver_policy.jl:642-642`
<!-- vulcan:connections:end -->

## Limitations
The two symplectic modes are unreachable from the typed `run_simulation` path today and fail fast rather than degrading, which is correct but means the configuration is only rejected at solve time. Forcing per-step storage in the autoswitching branch costs memory on long spans even when the caller asked for a summary solution. The metadata tuple reports `fallback_used` from an internal autoswitch, so a run that switched integrators mid-span cannot be attributed to a specific epoch from the trace alone.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl:631-748`, with the tolerance builder at line 15, the integrator cache type at line 257, the multirate driver at line 406 and the two-argument convenience method at line 750 of the same file.
