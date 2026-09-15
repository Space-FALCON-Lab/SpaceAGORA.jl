---
id: gnc.pso_parameters_rpopsocullsettings
label: RPOPSOCullSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOCullSettings
  lines:
  - 69
  - 69
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `true`).
- id: fraction_max
  type: Float64
  units: n/a
  required: false
  description: Field `fraction_max` (default `0.35`).
- id: start_iter
  type: Int
  units: n/a
  required: false
  description: Field `start_iter` (default `8`).
- id: noise_scale
  type: Float64
  units: n/a
  required: false
  description: Field `noise_scale` (default `0.25`).
- id: arc_velocity_scale
  type: Float64
  units: n/a
  required: false
  description: Field `arc_velocity_scale` (default `0.12`).
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
  type: RPOPSOCullSettings
  units: n/a
  description: Constructed `RPOPSOCullSettings` (keyword constructor via @kwdef).
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

# RPOPSOCullSettings

## Purpose
Grouped struct for the particle culling and reseeding mechanism that periodically replaces the weakest fraction of the swarm with perturbed copies of better particles to maintain diversity.

## Design & Implementation
Fields: `enabled::Bool = true`; `fraction_max::Float64 = 0.35`, the largest share of particles that may be culled in one pass; `start_iter::Int = 8`, first iteration at which culling is permitted; `noise_scale::Float64 = 0.25`, relative magnitude of positional noise added to reseeded particles; `arc_velocity_scale::Float64 = 0.12`, scale of the velocity assigned along the path arc to reseeded particles. Mapped to `cull_enable`, `cull_fraction_max`, `cull_start_iter`, `cull_noise_scale`, and `cull_arc_velocity_scale` in `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `true`). |
| in | `fraction_max` | Float64 | n/a | no | Field `fraction_max` (default `0.35`). |
| in | `start_iter` | Int | n/a | no | Field `start_iter` (default `8`). |
| in | `noise_scale` | Float64 | n/a | no | Field `noise_scale` (default `0.25`). |
| in | `arc_velocity_scale` | Float64 | n/a | no | Field `arc_velocity_scale` (default `0.12`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOCullSettings | n/a | — | Constructed `RPOPSOCullSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:174-174`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Range enforcement is deferred to `validate_rpo_pso_config` (`0 <= fraction_max <= 1`, `start_iter >= 0`, non-negative scales). `noise_scale` and `arc_velocity_scale` are dimensionless multipliers whose reference quantity is defined in the swarm loop, not documented here.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 69.
