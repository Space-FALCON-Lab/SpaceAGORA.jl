---
id: gnc.propulsive_maneuvers__available_propellant_kg
label: _available_propellant_kg
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _available_propellant_kg
  lines:
  - 183
  - 183
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
- id: current_mass_kg
  type: Float64
  units: n/a
  required: true
  description: Positional argument `current_mass_kg`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_available_propellant_kg`.
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

# _available_propellant_kg

## Purpose
Estimates the propellant mass in kilograms still usable by spacecraft `i`, so a planned burn can be rejected when it would consume more than the tank holds.

## Design & Implementation
Returns `nothing` unless both `hasproperty(p, :args)` and `hasproperty(p.args, :dynamics_model)` hold. It indexes `p.args.dynamics_model.spacecraft`, returning `nothing` for `i` outside `1:length(spacecraft)`. It then returns `nothing` when `sc.prop_mass` is non-finite or non-positive, treating a propellant-less spacecraft as unconstrained rather than as having zero margin. Otherwise it returns `max(0.0, current_mass_kg - sc.dry_mass)`, the current mass less dry mass, floored at zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `current_mass_kg` | Float64 | n/a | yes | Positional argument `current_mass_kg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_available_propellant_kg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__validated_burn_plan|_validated_burn_plan]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:232-232`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A `nothing` return means the caller skips the propellant check altogether, so a spacecraft configured with `prop_mass = 0` can be scheduled for arbitrarily large burns. The estimate assumes current mass minus dry mass is entirely usable propellant, ignoring residuals, ullage and any non-propellant consumables, and it never returns a negative margin because of the `max(0.0, ...)` floor.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 183.
