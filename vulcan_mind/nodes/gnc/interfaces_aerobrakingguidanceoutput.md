---
id: gnc.interfaces_aerobrakingguidanceoutput
label: AerobrakingGuidanceOutput
kind: struct
source:
  file: src/gnc/guidance/aerobraking/interfaces.jl
  symbol: AerobrakingGuidanceOutput
  lines:
  - 21
  - 21
inputs:
- id: time_switch_1
  type: Float64
  units: n/a
  required: false
  description: Field `time_switch_1` (default `0.0`).
- id: time_switch_2
  type: Float64
  units: n/a
  required: false
  description: Field `time_switch_2` (default `0.0`).
- id: security_mode
  type: Bool
  units: n/a
  required: false
  description: Field `security_mode` (default `false`).
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
  type: AerobrakingGuidanceOutput
  units: n/a
  description: Constructed `AerobrakingGuidanceOutput` (keyword constructor via @kwdef).
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

# AerobrakingGuidanceOutput

## Purpose
Return record of an aerobraking guidance evaluation. It carries the two angle-of-attack switch times that bound the controlled arc of the drag passage, plus a flag reporting whether the law fell back to its conservative security mode.

## Design & Implementation
A `Base.@kwdef` immutable struct with three concretely typed fields: `time_switch_1::Float64` and `time_switch_2::Float64`, both seconds measured on the same clock as `AerobrakingGuidanceInput.t`, defaulting to `0.0`, and `security_mode::Bool = false`. Concrete field types make the struct isbits, so it is returned without heap allocation and can be stored cheaply in per-satellite guidance buffers between drag passes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `time_switch_1` | Float64 | n/a | no | Field `time_switch_1` (default `0.0`). |
| in | `time_switch_2` | Float64 | n/a | no | Field `time_switch_2` (default `0.0`). |
| in | `security_mode` | Bool | n/a | no | Field `security_mode` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerobrakingGuidanceOutput | n/a | — | Constructed `AerobrakingGuidanceOutput` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.e_edg_strategy_compute_e_edg_guidance_window_bang|compute_e_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:59-59`
- [[gncy.t_edg_strategy_compute_t_edg_guidance_window_bang|compute_t_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl:8-8`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
There is no validity flag distinguishing a genuine solution of `(0.0, 0.0)` — which the solvers do return when the predicted heat load never reaches the limit — from a default-constructed instance that was never filled in. Nothing enforces `time_switch_1 <= time_switch_2`, and the struct has no room for the diagnostics (root-solver residual, iteration count) that would let a caller judge solution quality.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/interfaces.jl` line 21.
