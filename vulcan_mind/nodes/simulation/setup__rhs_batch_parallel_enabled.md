---
id: simulation.setup__rhs_batch_parallel_enabled
label: _rhs_batch_parallel_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_batch_parallel_enabled
  lines:
  - 380
  - 380
inputs:
- id: env
  type: SimulationModel.RhsPlanEnvConfig
  units: n/a
  required: true
  description: Positional argument `env`.
- id: num_spacecraft
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_spacecraft`.
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
  description: Return value of `_rhs_batch_parallel_enabled`.
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

# _rhs_batch_parallel_enabled

## Purpose
Final yes/no decision on threading the outer satellite batch of an RHS evaluation, combining the serial-profile override, the operator mode, the satellite count, and the physical core count.

## Design & Implementation
Two methods. The core one takes `env::SimulationModel.RhsPlanEnvConfig` and `num_spacecraft::Int`: returns `false` if `env.profile_forces_serial`; then `false` for `mode == :off`, `true` for `:on`; and for `:auto` returns `num_spacecraft >= env.batch_thread_threshold && Polyester.num_cores() > 1`. The convenience method `(p, num_spacecraft)` extracts the env via `_rhs_env_config(p)` and forwards. Both are `@inline` and pure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `env` | SimulationModel.RhsPlanEnvConfig | n/a | yes | Positional argument `env`. |
| in | `num_spacecraft` | Int | n/a | yes | Positional argument `num_spacecraft`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_batch_parallel_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__gravity_backbone_half_kick_bang|_gravity_backbone_half_kick!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1529-1529`
- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2098-2098`
- [[simulation.dynamics_rhs_spacecraft_dynamics_fast_control_bang|spacecraft_dynamics_fast_control!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2206-2206`
- [[simulation.dynamics_rhs_spacecraft_dynamics_gravity_backbone_bang|spacecraft_dynamics_gravity_backbone!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1560-1560`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1999-1999`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1843-1843`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1729-1729`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`Polyester.num_cores()` reports physical cores, not the Julia thread count, so on a machine started with `-t 1` the `:auto` path can return `true` and the batch loop then runs on a single thread with scheduling overhead. `:on` ignores `num_spacecraft`, so a single-satellite run can be 'parallelised'.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 380.
