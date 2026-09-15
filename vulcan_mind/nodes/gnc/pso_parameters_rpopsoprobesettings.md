---
id: gnc.pso_parameters_rpopsoprobesettings
label: RPOPSOProbeSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOProbeSettings
  lines:
  - 108
  - 108
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `true`).
- id: max_depth
  type: Int
  units: n/a
  required: false
  description: Field `max_depth` (default `2`).
- id: candidates
  type: Int
  units: n/a
  required: false
  description: Field `candidates` (default `24`).
- id: offset_scale
  type: Float64
  units: n/a
  required: false
  description: Field `offset_scale` (default `1.0`).
- id: sample_ds_m
  type: Float64
  units: n/a
  required: false
  description: Field `sample_ds_m` (default `0.25`).
- id: seed
  type: Int
  units: n/a
  required: false
  description: Field `seed` (default `1234`).
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
  type: RPOPSOProbeSettings
  units: n/a
  description: Constructed `RPOPSOProbeSettings` (keyword constructor via @kwdef).
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

# RPOPSOProbeSettings

## Purpose
Grouped struct configuring the quick geometry probes that adaptive sizing runs before the main swarm to estimate scene complexity by casting a small set of candidate detours around the station.

## Design & Implementation
Fields: `enabled::Bool = true`; `max_depth::Int = 2` recursion depth of detour probing; `candidates::Int = 24` candidate offsets tried per level; `offset_scale::Float64 = 1.0` multiplier on the probe offset distance; `sample_ds_m::Float64 = 0.25` metres between collision samples during probing (coarser than the main `sample_ds_m`); `seed::Int = 1234` RNG seed so probe results are deterministic. Copied to `probe_*` fields of `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `true`). |
| in | `max_depth` | Int | n/a | no | Field `max_depth` (default `2`). |
| in | `candidates` | Int | n/a | no | Field `candidates` (default `24`). |
| in | `offset_scale` | Float64 | n/a | no | Field `offset_scale` (default `1.0`). |
| in | `sample_ds_m` | Float64 | n/a | no | Field `sample_ds_m` (default `0.25`). |
| in | `seed` | Int | n/a | no | Field `seed` (default `1234`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOProbeSettings | n/a | — | Constructed `RPOPSOProbeSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:178-178`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`validate_rpo_pso_config` requires `max_depth`, `candidates`, and `offset_scale` non-negative and `sample_ds_m` positive. The fixed seed makes probing reproducible but also means every planner instance explores the same probe pattern regardless of the caller's RNG.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 108.
