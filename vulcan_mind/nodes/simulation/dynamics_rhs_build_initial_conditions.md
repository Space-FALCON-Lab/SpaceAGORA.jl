---
id: simulation.dynamics_rhs_build_initial_conditions
label: build_initial_conditions
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: build_initial_conditions
  lines:
  - 2289
  - 2289
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
  type: ComponentVector
  units: n/a
  description: Return value of `build_initial_conditions`.
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

# build_initial_conditions

## Purpose
Constructs the `ComponentVector` initial state for all satellites, with per-satellite position, velocity, mass, heat loads and optionally attitude.

## Design & Implementation
Builds a named-tuple shape per satellite depending on `orientation_sim`, converts each initial condition to J2000 position and velocity, fills the shapes, and assembles the component vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ComponentVector | n/a | — | Return value of `build_initial_conditions`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:177-177`

**Downstream**

- `callees` → [[core.project_unit_quaternion|project_unit_quaternion]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2338-2338`
- `callees` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2330-2330`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_coupled_cloth_robot_arm_state_shape|coupled_cloth_robot_arm_state_shape]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2317-2317`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_initialize_coupled_cloth_robot_arm_state_bang|initialize_coupled_cloth_robot_arm_state!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2346-2346`
- `callees` → [[simulation.dynamics_rhs__robot_arm_coupling|_robot_arm_coupling]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2315-2315`
<!-- vulcan:connections:end -->

## Limitations
Heat loads are sized to the link count, so link changes require rebuilding.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 2289.
