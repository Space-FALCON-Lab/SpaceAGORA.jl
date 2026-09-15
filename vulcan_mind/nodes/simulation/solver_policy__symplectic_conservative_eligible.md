---
id: simulation.solver_policy__symplectic_conservative_eligible
label: _symplectic_conservative_eligible
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _symplectic_conservative_eligible
  lines:
  - 100
  - 100
inputs:
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
  type: Bool
  units: n/a
  description: Return value of `_symplectic_conservative_eligible`.
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

# _symplectic_conservative_eligible

## Purpose
Checks that the run is a pure conservative single-body two-body problem so that the symplectic `KahanLi8` mode is physically appropriate.

## Design & Implementation
Returns `false` if `orientation_sim` is on, if there is not exactly one spacecraft, if any control, guidance, or navigation effectors exist, or if the dynamics model has anything other than exactly one effector. Returns `true` only when that sole effector `isa SimulationModel.InverseSquaredGravityModel`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_symplectic_conservative_eligible`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:636-636`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
J2, harmonics, and N-body gravity are excluded even though they are conservative, so the symplectic mode is limited to point-mass gravity. The check does not inspect whether the problem is actually partitioned; that is done separately. Multiple spacecraft in a conservative field are rejected regardless.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 100.
