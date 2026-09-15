---
id: core.reference_system_config_clock
label: clock
kind: struct
source:
  file: src/core/state/reference_system_config.jl
  symbol: clock
  lines:
  - 38
  - 38
inputs:
- id: year
  type: Int64
  units: n/a
  required: true
  description: Field `year`.
- id: month
  type: Int64
  units: n/a
  required: true
  description: Field `month`.
- id: day
  type: Int64
  units: n/a
  required: true
  description: Field `day`.
- id: hour
  type: Int64
  units: n/a
  required: true
  description: Field `hour`.
- id: minute
  type: Int64
  units: n/a
  required: true
  description: Field `minute`.
- id: second
  type: Float64
  units: n/a
  required: true
  description: Field `second`.
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
  type: clock
  units: n/a
  description: Constructed `clock`.
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

# clock

## Purpose

`clock` is the mutable calendar-date-and-time container used to carry an epoch alongside a reference-system state. It stores `year`, `month`, `day`, `hour` and `minute` as `Int64`, and `second` as a `Float64` so sub-second precision survives.

## Design & Implementation

It is declared as a bare `mutable struct` with no inner constructor and no validation, so the positional constructor `clock(year, month, day, hour, minute, second)` accepts any integers at all. The split between integer calendar fields and a floating-point seconds field is deliberate: it keeps the date exact while letting fractional seconds be represented for time-of-epoch conversions such as Julian date computation, which live outside this module.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `year` | Int64 | n/a | yes | Field `year`. |
| in | `month` | Int64 | n/a | yes | Field `month`. |
| in | `day` | Int64 | n/a | yes | Field `day`. |
| in | `hour` | Int64 | n/a | yes | Field `hour`. |
| in | `minute` | Int64 | n/a | yes | Field `minute`. |
| in | `second` | Float64 | n/a | yes | Field `second`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | clock | n/a | — | Constructed `clock`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:69-69`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:74-74`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:91-91`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:56-56`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:69-69`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:74-74`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:91-91`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:56-56`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/state/reference_system_config.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Nothing enforces field ranges, so a `month` of 13 or a negative `hour` is constructible and only fails downstream. The type records no time scale, so UTC, TAI, TT and UT1 epochs are indistinguishable, which matters for leap-second-sensitive conversions. Representing seconds as `Float64` limits precision to roughly the microsecond level over a long mission span, and the lower-case type name reads like a function at call sites.

## Provenance
Mapped from `src/core/state/reference_system_config.jl` line 38.
