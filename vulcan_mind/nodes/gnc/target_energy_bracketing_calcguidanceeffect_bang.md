---
id: gnc.target_energy_bracketing_calcguidanceeffect_bang
label: calcGuidanceEffect!
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: calcGuidanceEffect!
  lines:
  - 191
  - 191
inputs:
- id: model
  type: AerobrakingEnergyDepletionGuidanceModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  type: Nothing
  units: n/a
  description: Return value of `calcGuidanceEffect!`; mutates `model` in place. Returns
    `nothing`.
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

# calcGuidanceEffect!

## Purpose
`calcGuidanceEffect!` is the guidance hook the engine invokes for the `AerobrakingEnergyDepletionGuidanceModel` at each guidance update. It decides, per spacecraft, whether to run target-energy bracketing, run open-loop maximum energy depletion, or fall back to a safe low-drag attitude, and records that decision in the mutable state.

## Design & Implementation
Signature `(model::AerobrakingEnergyDepletionGuidanceModel, u, p::ODEParams, t::Float64, i::Int64)`. It returns `nothing` immediately unless `_edg_state_index_ok(state, i)`. If `:targeting in config.guidance_modes` it delegates to `_edg_run_target_energy_bracketing!(model, u, p, Float64(t), i)`, which performs the low-drag/max-depletion bracket propagation once per drag passage and sets `selected_mode`, `targeting_active`, `safe_low_drag`, target and bracket energies. Otherwise, if `:max_energy_depletion` is enabled it writes `selected_mode[i] = :max_energy_depletion` and clears `targeting_active`, `safe_low_drag`, `energy_bracketing_evaluated`; if neither mode is enabled it selects `:safe_low_drag` with `safe_low_drag[i] = true`. Only `model.state` is mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerobrakingEnergyDepletionGuidanceModel | n/a | yes | Positional argument `model`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `calcGuidanceEffect!`; mutates `model` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- [[simulation.control_callbacks__run_guidance_for_thruster_schedule_bang|_run_guidance_for_thruster_schedule!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:30-30`
- [[simulation.navigation_guidance_callbacks_get_guidance_callbacks|get_guidance_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:27-27`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:203-203`
- `callees` → [[gnc.target_energy_bracketing__edg_state_index_ok|_edg_state_index_ok]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:199-199`
- `callees` → [[gncz.target_energy_bracketing__edg_run_target_energy_bracketing_bang|_edg_run_target_energy_bracketing!]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:203-203`
<!-- vulcan:connections:end -->

## Limitations
The branch order gives `:targeting` precedence whenever present, so `guidance_modes=(:max_energy_depletion, :targeting)` behaves the same as `(:targeting,)` until bracketing declares the target unreachable. In targeting mode the bracketing is evaluated at most once per passage (`energy_bracketing_evaluated` gate) and nothing here clears that flag, so passage-to-passage reset depends on the control module. The `i::Int64` annotation makes the method fail to match on 32-bit `Int` platforms. No time-based throttling exists; the function runs on every guidance call.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 191.
