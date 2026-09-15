---
id: gnc.pso_parameters_rpopsoobjectivesettings
label: RPOPSOObjectiveSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOObjectiveSettings
  lines:
  - 13
  - 13
inputs:
- id: w_len
  type: Float64
  units: n/a
  required: false
  description: Field `w_len` (default `1.0`).
- id: w_obs
  type: Float64
  units: n/a
  required: false
  description: Field `w_obs` (default `1.0e6`).
- id: w_fuel
  type: Float64
  units: n/a
  required: false
  description: Field `w_fuel` (default `1.0`).
- id: obstacle_sigmoid_k
  type: Float64
  units: n/a
  required: false
  description: Field `obstacle_sigmoid_k` (default `1.0e6`).
- id: obstacle_sigmoid_tol_m
  type: Float64
  units: n/a
  required: false
  description: Field `obstacle_sigmoid_tol_m` (default `0.0`).
- id: w_inertia
  type: Float64
  units: n/a
  required: false
  description: Field `w_inertia` (default `0.7`).
- id: c1
  type: Float64
  units: n/a
  required: false
  description: Field `c1` (default `1.4`).
- id: c2
  type: Float64
  units: n/a
  required: false
  description: Field `c2` (default `1.4`).
- id: cost_ref_distance_m
  type: Float64
  units: n/a
  required: false
  description: Field `cost_ref_distance_m` (default `20.0`).
- id: mass_kg
  type: Float64
  units: n/a
  required: false
  description: Field `mass_kg` (default `12.0`).
- id: tf_s
  type: Float64
  units: n/a
  required: false
  description: Field `tf_s` (default `120.0`).
- id: isp_s
  type: Float64
  units: n/a
  required: false
  description: Field `isp_s` (default `60.0`).
- id: g0_mps2
  type: Float64
  units: n/a
  required: false
  description: Field `g0_mps2` (default `9.80665`).
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
  type: RPOPSOObjectiveSettings
  units: n/a
  description: Constructed `RPOPSOObjectiveSettings` (keyword constructor via @kwdef).
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

# RPOPSOObjectiveSettings

## Purpose
Grouped struct bundling the HYPR path-cost weights, the obstacle sigmoid shaping, the base PSO velocity coefficients, and the physical constants used by the fuel proxy (mass, transfer time, specific impulse, standard gravity).

## Design & Implementation
Fields and defaults: `w_len = 1.0`, `w_obs = 1.0e6`, `w_fuel = 1.0` (cost weights); `obstacle_sigmoid_k = 1.0e6` and `obstacle_sigmoid_tol_m = 0.0` (steepness and tolerance of the soft collision indicator); `w_inertia = 0.7`, `c1 = 1.4`, `c2 = 1.4` (PSO inertia, cognitive, social coefficients); `cost_ref_distance_m = 20.0` (length normalisation reference); `mass_kg = 12.0`, `tf_s = 120.0`, `isp_s = 60.0`, `g0_mps2 = 9.80665` for the rocket-equation fuel proxy. All are `Float64` and flattened one-to-one into `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `w_len` | Float64 | n/a | no | Field `w_len` (default `1.0`). |
| in | `w_obs` | Float64 | n/a | no | Field `w_obs` (default `1.0e6`). |
| in | `w_fuel` | Float64 | n/a | no | Field `w_fuel` (default `1.0`). |
| in | `obstacle_sigmoid_k` | Float64 | n/a | no | Field `obstacle_sigmoid_k` (default `1.0e6`). |
| in | `obstacle_sigmoid_tol_m` | Float64 | n/a | no | Field `obstacle_sigmoid_tol_m` (default `0.0`). |
| in | `w_inertia` | Float64 | n/a | no | Field `w_inertia` (default `0.7`). |
| in | `c1` | Float64 | n/a | no | Field `c1` (default `1.4`). |
| in | `c2` | Float64 | n/a | no | Field `c2` (default `1.4`). |
| in | `cost_ref_distance_m` | Float64 | n/a | no | Field `cost_ref_distance_m` (default `20.0`). |
| in | `mass_kg` | Float64 | n/a | no | Field `mass_kg` (default `12.0`). |
| in | `tf_s` | Float64 | n/a | no | Field `tf_s` (default `120.0`). |
| in | `isp_s` | Float64 | n/a | no | Field `isp_s` (default `60.0`). |
| in | `g0_mps2` | Float64 | n/a | no | Field `g0_mps2` (default `9.80665`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOObjectiveSettings | n/a | — | Constructed `RPOPSOObjectiveSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:171-171`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct does not validate; `validate_rpo_pso_config` later requires `obstacle_sigmoid_k > 0`, `mass_kg > 0`, `tf_s > 0`, `isp_s > 0`, `g0_mps2 > 0`, and `cost_ref_distance_m >= 0`, throwing `ArgumentError` otherwise. Negative cost weights are not rejected anywhere. The `w_obs` default of 1e6 makes any collision dominate the cost, which is intentional but means the length and fuel terms only differentiate feasible paths.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 13.
