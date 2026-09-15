---
id: gnc.pso_parameters_rpopsostagnationsettings
label: RPOPSOStagnationSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOStagnationSettings
  lines:
  - 90
  - 90
inputs:
- id: stagnation_learning_enable
  type: Bool
  units: n/a
  required: false
  description: Field `stagnation_learning_enable` (default `true`).
- id: stagnation_learning_threshold
  type: Int
  units: n/a
  required: false
  description: Field `stagnation_learning_threshold` (default `8`).
- id: stagnation_learning_elite_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `stagnation_learning_elite_fraction` (default `0.10`).
- id: stagnation_learning_max_blocks
  type: Int
  units: n/a
  required: false
  description: Field `stagnation_learning_max_blocks` (default `2`).
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
  type: RPOPSOStagnationSettings
  units: n/a
  description: Constructed `RPOPSOStagnationSettings` (keyword constructor via @kwdef).
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

# RPOPSOStagnationSettings

## Purpose
Grouped struct for HYPR's stagnation-learning response: when the global best has not improved for a number of iterations, the swarm learns from an elite subset by adjusting a bounded number of waypoint blocks.

## Design & Implementation
Fields: `stagnation_learning_enable::Bool = true`; `stagnation_learning_threshold::Int = 8`, consecutive non-improving iterations before learning triggers; `stagnation_learning_elite_fraction::Float64 = 0.10`, share of top particles used as the elite pool; `stagnation_learning_max_blocks::Int = 2`, maximum waypoint blocks modified per event. Unlike the other grouped structs, its field names already carry the `stagnation_learning_` prefix and are copied verbatim into `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `stagnation_learning_enable` | Bool | n/a | no | Field `stagnation_learning_enable` (default `true`). |
| in | `stagnation_learning_threshold` | Int | n/a | no | Field `stagnation_learning_threshold` (default `8`). |
| in | `stagnation_learning_elite_fraction` | Float64 | n/a | no | Field `stagnation_learning_elite_fraction` (default `0.10`). |
| in | `stagnation_learning_max_blocks` | Int | n/a | no | Field `stagnation_learning_max_blocks` (default `2`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOStagnationSettings | n/a | — | Constructed `RPOPSOStagnationSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:176-176`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Validation is external: threshold and max_blocks must be non-negative and elite_fraction within [0, 1] or `validate_rpo_pso_config` throws `ArgumentError`. A threshold of 0 would trigger learning every iteration, which is permitted but not guarded.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 90.
