---
id: gnc.target_energy_bracketing__edg_set_targeting_fallback_bang
label: _edg_set_targeting_fallback!
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_set_targeting_fallback!
  lines:
  - 250
  - 250
inputs:
- id: config
  type: AerobrakingEnergyDepletionConfig
  units: n/a
  required: true
  description: Positional argument `config`.
- id: state
  type: AerobrakingEnergyDepletionState
  units: n/a
  required: true
  description: Positional argument `state`.
- id: i
  type: Int
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
  description: Return value of `_edg_set_targeting_fallback!`; mutates `config` in
    place. Returns `nothing`.
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

# _edg_set_targeting_fallback!

## Purpose
`_edg_set_targeting_fallback!` puts spacecraft `i` into the non-targeting fallback mode when bracketing cannot run, which happens when `_edg_run_target_energy_bracketing!` finds the vehicle outside a drag passage. It chooses maximum energy depletion if that mode is enabled and safe low-drag otherwise.

## Design & Implementation
Signature `(config::AerobrakingEnergyDepletionConfig, state::AerobrakingEnergyDepletionState, i::Int)`. It sets `state.targeting_active[i] = false`, `state.safe_low_drag[i] = !(:max_energy_depletion in config.guidance_modes)`, and `state.selected_mode[i]` to `:max_energy_depletion` when that symbol is in `config.guidance_modes` and to `:safe_low_drag` otherwise. Only these three vectors are mutated; `energy_bracketing_evaluated`, energies and counters are left untouched so bracketing can still run once the drag passage begins. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `state` | AerobrakingEnergyDepletionState | n/a | yes | Positional argument `state`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_edg_set_targeting_fallback!`; mutates `config` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.target_energy_bracketing__edg_run_target_energy_bracketing_bang|_edg_run_target_energy_bracketing!]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:277-277`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The function does not bounds-check `i`; it relies on `calcGuidanceEffect!` having already passed `_edg_state_index_ok`. Because `energy_bracketing_evaluated[i]` is not reset here, a spacecraft that was previously bracketed keeps its stale `target_energy_jkg` while in fallback. The membership test `:max_energy_depletion in config.guidance_modes` is a linear scan over the tuple on every call, negligible for the two-element vocabulary.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 250.
