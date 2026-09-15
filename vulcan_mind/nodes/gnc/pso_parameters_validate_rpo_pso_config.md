---
id: gnc.pso_parameters_validate_rpo_pso_config
label: validate_rpo_pso_config
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: validate_rpo_pso_config
  lines:
  - 593
  - 593
inputs:
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
  description: Return value of `validate_rpo_pso_config`. Returns `cfg`.
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

# validate_rpo_pso_config

## Purpose
Range-checks every numerically constrained field of an `RPOPSOConfig` and throws a descriptive `ArgumentError` on the first violation, guaranteeing the planner hot path never sees non-positive step sizes, inverted min/max pairs, or invalid curve types.

## Design & Implementation
Signature `validate_rpo_pso_config(cfg::RPOPSOConfig)`; returns `cfg` unchanged on success. Checks are sequential `cond || throw(ArgumentError(msg))` lines covering roughly 90 fields: counts such as `n_waypoints >= 0`, `n_particles > 0`, `n_iters >= 0`; positive spacings `sample_ds_m`, `refinement_sample_ds_m`, `probe_sample_ds_m`, `adaptive_sampling_max_ds_m`; `curve_type in (:bezier, :polyline)`; physical constants `mass_kg`, `tf_s`, `isp_s`, `g0_mps2 > 0`; retiming limits including `retime_max_speed_mps >= retime_min_speed_mps`; ordered pairs for every adaptive `*_min <= *_max`; fractions such as `cull_fraction_max`, `stagnation_learning_elite_fraction`, `rrt_warmstart_goal_sample_rate` in [0, 1] and `adaptive_sampling_obstacle_guard_fraction` in (0, 1]; and `schedule_c_min <= schedule_c_max`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `validate_rpo_pso_config`. Returns `cfg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:586-586`
- [[gncy.pso_adaptive_policy_rpo_adaptive_pso_config|rpo_adaptive_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:28-28`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the first failing check is reported. The cost weights `w_len`, `w_obs`, `w_fuel`, PSO coefficients `w_inertia`, `c1`, `c2`, `spread_scale`, `search_margin_m`, and the schedule end fractions are not validated at all, so negative or zero values pass. `adaptive_n_waypoints_min <= adaptive_n_waypoints_max` (and the particle/iteration analogues) are not checked, only individual non-negativity.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 593.
