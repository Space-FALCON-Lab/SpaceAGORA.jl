---
id: gnc.rpo_plan_buffer_rpoplanbuffer
label: RPOPlanBuffer
kind: struct
source:
  file: src/gnc/guidance/rpo/rpo_plan_buffer.jl
  symbol: RPOPlanBuffer
  lines:
  - 13
  - 13
inputs:
- id: valid
  type: Bool
  units: n/a
  required: false
  description: Field `valid` (default `false`).
- id: plan
  type: RPOPlan
  units: n/a
  required: false
  description: Field `plan` (default `RPOPlan()`).
- id: updated_at_s
  type: Float64
  units: n/a
  required: false
  description: Field `updated_at_s` (default `NaN`).
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
  type: RPOPlanBuffer
  units: n/a
  description: Constructed `RPOPlanBuffer` (keyword constructor via @kwdef).
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

# RPOPlanBuffer

## Purpose
Per-vehicle slot that carries the active RPO plan across guidance updates, so a plan computed on one call remains available on later calls without re-solving. It also records when the plan was installed, letting guidance decide whether the reference is stale.

## Design & Implementation
A `Base.@kwdef mutable struct` with three fields: `valid::Bool = false`, `plan::RPOPlan = RPOPlan()` and `updated_at_s::Float64 = NaN`. The `NaN` default is deliberate — any age computation against an untouched buffer yields `NaN` rather than a plausible-looking zero. `update_rpo_plan_buffer!(buffer, plan, t)` is the sole mutator: it assigns `buffer.plan = plan`, mirrors the plan's own flag into `buffer.valid = plan.valid`, stores `Float64(t)` into `updated_at_s`, and returns the buffer.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `valid` | Bool | n/a | no | Field `valid` (default `false`). |
| in | `plan` | RPOPlan | n/a | no | Field `plan` (default `RPOPlan()`). |
| in | `updated_at_s` | Float64 | n/a | no | Field `updated_at_s` (default `NaN`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPlanBuffer | n/a | — | Constructed `RPOPlanBuffer` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_control_types_rpompccontrolmodel|RPOMPCControlModel]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_control_types.jl:15-15`
- [[gncy.rpo_guidance_models_rpoguidancemodel|RPOGuidanceModel]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_models.jl:7-7`

**Downstream**

- `callees` → [[gnc.rpo_plan_buffer_rpoplan|RPOPlan]] · `callers` · call · `src/gnc/guidance/rpo/rpo_plan_buffer.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
Despite the docstring promising the previous plan is preserved for diagnostics, the struct has no field for it and `update_rpo_plan_buffer!` overwrites `plan` outright, so the prior plan is dropped. The buffer stores the plan by reference, so later mutation of that `RPOPlan` is visible through the buffer. There is no lock, so concurrent satellite propagation touching a shared buffer races on all three fields.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_plan_buffer.jl` line 13.
