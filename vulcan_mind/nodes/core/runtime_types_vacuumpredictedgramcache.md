---
id: core.runtime_types_vacuumpredictedgramcache
label: VacuumPredictedGRAMCache
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: VacuumPredictedGRAMCache
  lines:
  - 530
  - 530
inputs:
- id: valid
  type: Bool
  units: n/a
  required: true
  description: Field `valid`.
- id: t0
  type: Float64
  units: n/a
  required: true
  description: Field `t0`.
- id: t1
  type: Float64
  units: n/a
  required: true
  description: Field `t1`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Field `h`.
- id: log_rhos
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `log_rhos`.
- id: Ms_rho
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `Ms_rho`.
- id: Ts
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `Ts`.
- id: Ms_T
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `Ms_T`.
- id: winds
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `winds`.
- id: vac_alts
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `vac_alts`.
- id: vac_positions
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `vac_positions`.
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
  type: VacuumPredictedGRAMCache
  units: n/a
  description: Constructed `VacuumPredictedGRAMCache`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# VacuumPredictedGRAMCache

## Purpose
Per-satellite look-ahead cache of GRAM density along a drag-free predicted trajectory, letting the RHS interpolate smooth log-density splines instead of sampling the noisy native model at every stage.

## Design & Implementation
A mutable struct with a `valid` flag, knot span `t0` to `t1` with uniform spacing `h`, the log-density and temperature knot values with their natural cubic spline second derivatives, wind vectors for linear interpolation, and the vacuum-predicted altitude and inertial position at each knot so a query can measure how far the real trajectory has deviated from the prediction. Rebuilt when deviation exceeds the configured tolerance.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `valid` | Bool | n/a | yes | Field `valid`. |
| in | `t0` | Float64 | n/a | yes | Field `t0`. |
| in | `t1` | Float64 | n/a | yes | Field `t1`. |
| in | `h` | Float64 | n/a | yes | Field `h`. |
| in | `log_rhos` | Vector{Float64} | n/a | yes | Field `log_rhos`. |
| in | `Ms_rho` | Vector{Float64} | n/a | yes | Field `Ms_rho`. |
| in | `Ts` | Vector{Float64} | n/a | yes | Field `Ts`. |
| in | `Ms_T` | Vector{Float64} | n/a | yes | Field `Ms_T`. |
| in | `winds` | Vector{SVector{3, Float64}} | n/a | yes | Field `winds`. |
| in | `vac_alts` | Vector{Float64} | n/a | yes | Field `vac_alts`. |
| in | `vac_positions` | Vector{SVector{3, Float64}} | n/a | yes | Field `vac_positions`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VacuumPredictedGRAMCache | n/a | — | Constructed `VacuumPredictedGRAMCache`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`
- [[simulation.vacuum_predicted_gram__vacuum_gram_cache_for_sat_bang|_vacuum_gram_cache_for_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:53-53`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The prediction is two-body plus J2 with no drag, so during a deep pass the real trajectory diverges and the cache invalidates often, paying a full rebuild each time; uniform knot spacing cannot concentrate resolution near periapsis.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 530.
