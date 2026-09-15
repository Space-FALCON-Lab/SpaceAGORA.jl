---
id: grp.src_simulation_engine
label: simulation/engine/
kind: group
inputs:
- id: members_in
  type: call
  units: n/a
  required: false
  description: Calls into any member of this block from outside it.
outputs:
- id: members_out
  type: call
  units: n/a
  description: Calls from any member of this block to nodes outside it.
expands: module.simulation
tags:
- cluster
charts:
- simulation
origin: agent
---

# simulation/engine/

## Purpose
The simulation engine: takes a `SimulationConfiguration`, sets up buffers and caches, builds the solver problem and callbacks, routes each right-hand-side evaluation, integrates, checkpoints and persists results.

## Design & Implementation
`setup.jl` initialises everything and plans RHS execution; `dynamics_rhs.jl` is the family of right-hand sides; `effector_sampling.jl` samples environment for effectors; `solver_policy.jl` picks and configures solvers; `execution.jl` runs; `persistence.jl` and `resume_checkpoint.jl` write results and checkpoints; `rhs_calibration.jl` auto-tunes the RHS plan; `adapters/` and `config/` read engine settings; `public_api.jl` is `run_simulation`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `members_in` | call | n/a | no | Calls into any member of this block from outside it. |
| out | `members_out` | call | n/a | — | Calls from any member of this block to nodes outside it. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.io|IOConfig — internals]] · `members_out` → `members_in` · call · `src/io/serialization/io_serialization.jl:63-63`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `members_in` · call · `src/simulation/callbacks/density_callbacks/config.jl:303-303`
- [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_out` → `members_in` · call · `src/simulation/campaigns/monte_carlo.jl:223-223`
- [[module.simulation|RuntimeServices]] · `api` → `members_in` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `members_out` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:178-178`
- `members_out` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/reporting.jl:31-31`
- `members_out` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:47-47`
- `members_out` → [[core.simulation_configuration_solverconfig|SolverConfig]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:128-128`
- `members_out` → [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callers` · call · `src/simulation/engine/setup.jl:1493-1493`
- `members_out` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- `members_out` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/execution.jl:78-78`
- `members_out` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- `members_out` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- `members_out` → [[grp.cli|SpaceAGORACLI — internals]] · `members_in` · call · `src/simulation/engine/rhs_calibration.jl:261-261`
- `members_out` → [[grp.io|IOConfig — internals]] · `members_in` · call · `src/simulation/engine/execution.jl:92-92`
- `members_out` → [[grp.src_analysis_verification|analysis/verification/]] · `members_in` · call · `src/simulation/engine/dynamics_rhs.jl:178-178`
- `members_out` → [[grp.src_core_interfaces|core/interfaces/]] · `members_in` · call · `src/simulation/engine/dynamics_rhs.jl:2330-2330`
- `members_out` → [[grp.src_core_numerics|core/numerics/]] · `members_in` · call · `src/simulation/engine/dynamics_rhs.jl:2338-2338`
- `members_out` → [[grp.src_core_state|core/state/]] · `members_in` · call · `src/simulation/engine/adapters/from_env.jl:128-128`
- `members_out` → [[grp.src_core_types|core/types/]] · `members_in` · call · `src/simulation/engine/setup.jl:75-75`
- `members_out` → [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_in` · call · `src/simulation/engine/setup.jl:1481-1481`
- `members_out` → [[grp.src_dynamics_multibody_cloth|dynamics/multibody_cloth/]] · `members_in` · call · `src/simulation/engine/dynamics_rhs.jl:2317-2317`
- `members_out` → [[grp.src_environment_ephemerides|environment/ephemerides/]] · `members_in` · call · `src/simulation/engine/setup.jl:1499-1499`
- `members_out` → [[grp.src_environment_gravity|environment/gravity/]] · `members_in` · call · `src/simulation/engine/setup.jl:75-75`
- `members_out` → [[grp.src_gnc_control|gnc/control/]] · `members_in` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- `members_out` → [[grp.src_gnc_guidance|gnc/guidance/]] · `members_in` · call · `src/simulation/engine/execution.jl:78-78`
- `members_out` → [[grp.src_parallel_policy|parallel/policy/]] · `members_in` · call · `src/simulation/engine/adapters/from_env.jl:125-125`
- `members_out` → [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_in` · call · `src/simulation/engine/execution.jl:206-206`
- `members_out` → [[grp.src_vehicle_thermal|vehicle/thermal/]] · `members_in` · call · `src/simulation/engine/setup.jl:60-60`
- `members_out` → [[io.io_serialization__load_checkpoint|_load_checkpoint]] · `callers` · call · `src/simulation/engine/execution.jl:237-237`
- `members_out` → [[misc.io_serialization_write_checkpoint__write_checkpoint_bang|_write_checkpoint!]] · `callers` · call · `src/simulation/engine/execution.jl:387-387`
- `members_out` → [[parcore.effector_sampling_environmentsample|EnvironmentSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:298-298`
<!-- vulcan:connections:end -->

## Limitations
`setup.jl` alone is 1,900 lines; the engine's environment-variable surface is wide and mostly documented in comments.

## Provenance
Macro block generated by `vulcan compile` (D13) from `src/simulation/engine`; groups 337 nodes.
Its members are listed in chart `simulation-simulation-engine`.
