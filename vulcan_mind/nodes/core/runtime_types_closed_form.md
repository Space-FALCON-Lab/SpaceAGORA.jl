---
id: core.runtime_types_closed_form
label: Closed_form
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Closed_form
  lines:
  - 461
  - 461
inputs:
- id: t_cf
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `t_cf` (default `[]`).
- id: h_cf
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `h_cf` (default `[]`).
- id: gamma_cf
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `γ_cf` (default `[]`).
- id: v_cf
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `v_cf` (default `[]`).
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
  type: Closed_form
  units: n/a
  description: Constructed `Closed_form` (keyword constructor via @kwdef).
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

# Closed_form

## Purpose
Legacy result container for the closed-form (analytic) aerobraking pass prediction, storing the predicted altitude, flight-path angle and speed histories alongside their times.

## Design & Implementation
`@kwdef mutable struct Closed_form` with four `Vector{Float64}` fields defaulting to empty: `t_cf` (s), `h_cf` (m), `γ_cf` (rad), and `v_cf` (m/s). It is embedded in `Solution.closed_form` and appended to by legacy heat-load and drag-passage predictors that mirror the original Python aerobraking tool.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_cf` | Vector{Float64} | n/a | no | Field `t_cf` (default `[]`). |
| in | `h_cf` | Vector{Float64} | n/a | no | Field `h_cf` (default `[]`). |
| in | `gamma_cf` | Vector{Float64} | n/a | no | Field `γ_cf` (default `[]`). |
| in | `v_cf` | Vector{Float64} | n/a | no | Field `v_cf` (default `[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Closed_form | n/a | — | Constructed `Closed_form` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_solution|Solution]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:508-508`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Vectors are grown by `push!` with no capacity hints and no invariant that the four fields stay the same length. The naming uses a non-idiomatic underscore-capital style retained for parity with the Python source. Modern predictors return NamedTuples instead and do not populate this struct.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 461.
