---
id: gnc.propulsive_maneuvers__validated_burn_plan
label: _validated_burn_plan
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _validated_burn_plan
  lines:
  - 198
  - 198
inputs:
- id: controlModel
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `controlModel`.
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
- id: mass_kg
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass_kg`.
- id: maneuver
  type: Any
  units: n/a
  required: true
  description: Positional argument `maneuver`.
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
  description: Return value of `_validated_burn_plan`.
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

# _validated_burn_plan

## Purpose
Converts a commanded maneuver into a fully sized `PropulsiveBurnPlan` — impulse, propellant and feasibility — or returns `nothing` when the maneuver cannot be flown.

## Theory & Math
Exhaust velocity is $v_e = I_{sp}\,g_0$ with $g_0 = 9.80665\ \mathrm{m/s^2}$ and $I_{sp}$ in seconds. The rocket equation gives the mass fraction remaining after the burn, $f = e^{-\Delta v / v_e}$, so the propellant required from an initial wet mass $m$ is $m_p = m\,(1 - f)$. The commanded impulse follows as $J = m_p v_e$ in newton-seconds, and at constant thrust $F$ the burn lasts $t_b = J / F$ seconds. Here $\Delta v$ is the commanded velocity increment in m/s and $m$ is `mass_kg`.

## Design & Implementation
Takes thrust and specific impulse from `_effective_thrust_isp` and delta-v and direction from the `maneuver` named tuple. A single compound guard requires `mass_kg`, `delta_v_mps`, `thrust_n` and `isp_s` all finite and strictly positive and `direction_rad` finite, returning `nothing` otherwise. It then computes the mass fraction and rejects it unless it lies strictly in the open interval `(0, 1)`, which screens both an infinite delta-v and an underflowed exponential. Propellant, impulse and duration are each re-checked for finiteness and positivity. Finally `_available_propellant_kg` is consulted and the plan is rejected when `propellant_required_kg` exceeds the available mass by more than `1e-9` kg. On success it returns a `PropulsiveBurnPlan` with `valid=true`, no burn window yet, and `source_orbit` from the maneuver.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `mass_kg` | Float64 | n/a | yes | Positional argument `mass_kg`. |
| in | `maneuver` | Any | n/a | yes | Positional argument `maneuver`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_validated_burn_plan`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:519-519`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:206-206`
- `callees` → [[gnc.command_types_propulsiveburnplan|PropulsiveBurnPlan]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:237-237`
- `callees` → [[gnc.propulsive_maneuvers__available_propellant_kg|_available_propellant_kg]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:232-232`
- `callees` → [[gnc.propulsive_maneuvers__effective_thrust_isp|_effective_thrust_isp]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:205-205`
<!-- vulcan:connections:end -->

## Limitations
Every rejection path returns the same bare `nothing`, so a caller cannot distinguish an insufficient-propellant veto from a malformed command. The returned plan carries no start or stop time, leaving `calcControlEffect!` to fill the window, so a plan can be stored valid but unscheduled. The 1e-9 kg tolerance on the propellant check is absolute rather than relative and is meaningless for very large or very small vehicles.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 198.
