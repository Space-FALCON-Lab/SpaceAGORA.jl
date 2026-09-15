---
id: gnc.targeting_control__edg_disable_uncertified_targeting_bang
label: _edg_disable_uncertified_targeting!
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_disable_uncertified_targeting!
  lines:
  - 933
  - 933
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
- id: prefer_max_energy_depletion
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `prefer_max_energy_depletion`.
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
  type: Any
  units: n/a
  description: Return value of `_edg_disable_uncertified_targeting!`; mutates `config`
    in place. Returns `Inf`.
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

# _edg_disable_uncertified_targeting!

## Purpose
Falls back from targeting when the target cannot be certified as reachable, choosing between max-energy-depletion and safe low drag.

## Design & Implementation
Clears `targeting_active[i]`. With `prefer_max_energy_depletion` and that mode configured it selects max-energy-depletion and clears the safe flag; otherwise it selects safe-low-drag and sets the flag. Returns `Inf` as the switch time so the caller stores a value that compares as never-switch.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `state` | AerobrakingEnergyDepletionState | n/a | yes | Positional argument `state`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `prefer_max_energy_depletion` | Bool | n/a | yes | Keyword argument `prefer_max_energy_depletion`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_disable_uncertified_targeting!`; mutates `config` in place. Returns `Inf`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_evaluate_candidate|evaluate_candidate]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:994-994`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The fallback is permanent for the pass; targeting is not re-attempted even if conditions later change.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 933.
