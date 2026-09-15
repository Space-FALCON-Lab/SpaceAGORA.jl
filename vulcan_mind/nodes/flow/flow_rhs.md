---
id: flow.rhs
label: Dynamics right-hand side
kind: group
inputs:
- id: rhs_calls
  type: ODEFunction
  units: n/a
  description: Derivative requests from the integrator.
- id: staged_environment
  type: SharedBuffers
  units: n/a
  description: Atmosphere and planet frame staged by the callbacks.
- id: thread_policy
  type: PolicyDecision
  units: n/a
  description: Execution plan for this evaluation.
  required: false
- id: control_wrench
  type: (force, torque, mass rate)
  units: n/a
  description: Direct actuation returned by control effectors.
  required: false
outputs:
- id: force_requests
  type: StateSample + EnvironmentSample
  units: n/a
  description: State and environment handed to each force model.
- id: environment_queries
  type: sample_* calls
  units: n/a
  description: Planet-frame, atmosphere and ephemeris samples requested for effectors.
tags:
- master-flow
charts:
- master
origin: agent
opens: simulation-simulation-engine-dynamics-rhs-jl
---

# Dynamics right-hand side

## Purpose
Evaluates the state derivative for every satellite: sums the dynamic effectors' forces and torques, adds control effector wrenches and mass flow, applies robot-arm coupling, and writes translational, mass, heat-load and attitude derivatives.

## Design & Implementation
`dynamics_rhs.jl` routes each evaluation through an execution plan chosen by `_rhs_execution_plan` — serial, satellite batches, per-satellite effector threading, or the flat constellation effector queue with vectorised kernels for the batchable gravity, N-body and SRP models. Effectors are evaluated through the typed `wrench` interface when they implement it and the legacy `calcForceTorque` otherwise; shared samples (planet frame, atmosphere, Sun, third bodies) are prefilled once per call.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rhs_calls` | ODEFunction | n/a | — | Derivative requests from the integrator. |
| in | `staged_environment` | SharedBuffers | n/a | — | Atmosphere and planet frame staged by the callbacks. |
| in | `thread_policy` | PolicyDecision | n/a | no | Execution plan for this evaluation. |
| in | `control_wrench` | (force, torque, mass rate) | n/a | no | Direct actuation returned by control effectors. |
| out | `force_requests` | StateSample + EnvironmentSample | n/a | — | State and environment handed to each force model. |
| out | `environment_queries` | sample_* calls | n/a | — | Planet-frame, atmosphere and ephemeris samples requested for effectors. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.callbacks|Integration callbacks]] · `staged_environment` → `staged_environment` · dataflow · `src/simulation/engine/effector_sampling.jl`
- [[flow.gnc|Guidance, navigation & control]] · `control_wrench` → `control_wrench` · dataflow · `src/gnc/control/momentum_manager.jl`
- [[flow.parallel|Parallel routing & thread policy]] · `thread_policy` → `thread_policy` · dataflow · `src/simulation/engine/setup.jl`
- [[flow.solve_loop|Solve loop]] · `rhs_calls` → `rhs_calls` · dataflow · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `force_requests` → [[flow.forces|Force & torque models]] · `force_requests` · dataflow · `src/simulation/engine/dynamics_rhs.jl`
<!-- vulcan:connections:end -->

## Limitations
The routing chain has eleven return sites and a live cost model, so which path ran is only visible through the plan's `dominant_axis`; the per-satellite tail is duplicated across the monolithic, split and flat variants.
