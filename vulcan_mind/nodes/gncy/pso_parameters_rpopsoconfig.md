---
id: gncy.pso_parameters_rpopsoconfig
label: RPOPSOConfig
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOConfig
  lines:
  - 188
  - 307
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: settings_groups
  type: Tuple
  units: n/a
  required: true
  description: Grouped settings structs for swarm, objective, adaptive, sampling,
    culling, schedule, stagnation, early stopping, probe, reexplore, warmstart, refinement,
    and retiming behaviour.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: flat_config
  type: RPOPSOConfig
  units: n/a
  description: Flattened immutable configuration record read directly by the planner,
    sampler, cost, refinement, and retiming hot paths.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# RPOPSOConfig

## Purpose
`RPOPSOConfig` is the flattened configuration record that the entire RPO HYPR planning stack reads. It is deliberately a single flat immutable struct rather than a tree of nested settings, because the planner hot path touches these fields inside per-particle and per-sample loops where a nested field access chain would cost real time.

## Model & Assumptions
The file defines thirteen grouped keyword structs that express the configuration the way an operator thinks about it, covering swarm parameters, objective weights, adaptive policy bands, adaptive sampling, particle culling, weight scheduling, stagnation handling, early stopping, geometry probing, re-exploration, RRT-Connect warm start, refinement, and retiming. `RPOPSOConfigurator` groups them, and a constructor on line 446 flattens a configurator plus keyword overrides into this record. That gives a structured authoring surface and a flat consumption surface from one source of truth.

## Design & Implementation
Defaults encode the tuned operating point of the planner. The swarm runs 200 particles for 55 iterations over 5 interior waypoints with inertia 0.7 and both acceleration coefficients at 1.4. Objective weights are 1.0 for length, 1.0e6 for obstacles, and 1.0 for fuel, with the obstacle sigmoid sharpness at 1.0e6. Sampling defaults to 0.05 metre spacing on Bezier curves. Vehicle properties for the fuel proxy are 12 kg mass, 120 second transfer time, and 60 second specific impulse at standard gravity. Retiming defaults to 0.02 m/s^2 maximum acceleration, 0.25 second reaction time, and a 0.5 speed scale. Adaptive bands allow 3 to 8 waypoints, 140 to 320 particles, and an effort fraction from 0.75 to 1.5, with downscaling disabled by default. The functions `rpo_pso_config`, `_rpo_pso_sync_sample_ds_with_safe_distance`, and `validate_rpo_pso_config` handle derived updates and consistency checks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `settings_groups` | Tuple | n/a | yes | Grouped settings structs for swarm, objective, adaptive, sampling, culling, schedule, stagnation, early stopping, probe, reexplore, warmstart, refinement, and retiming behaviour. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `flat_config` | RPOPSOConfig | n/a | — | Flattened immutable configuration record read directly by the planner, sampler, cost, refinement, and retiming hot paths. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_740_mpc_final_pso_config|rpo_740_mpc_final_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:94-94`
- [[gnc.pso_parameters_push_bang|push!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:446-446`
- [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:583-583`
- [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:18-18`
- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:110-110`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Flattening duplicates every grouped field, so a new setting must be added in both the group struct and this record or it silently never reaches the planner. The struct is immutable, so every adjustment allocates a new record; the adaptive policy and the planner both do this several times per plan. Because sample spacing is synchronised against safe distance by a separate helper, setting the two independently can produce a record whose spacing does not resolve the keep-out shell it is meant to enforce.

## Provenance
Mapped from pso_parameters.jl lines 188-307; include site observed at guidance_hooks.jl line 65.
