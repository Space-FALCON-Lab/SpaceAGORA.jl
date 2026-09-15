---
id: gnc.propulsive_maneuvers__constant_thrust_burn_duration_s
label: _constant_thrust_burn_duration_s
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _constant_thrust_burn_duration_s
  lines:
  - 249
  - 249
inputs:
- id: mass_kg
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass_kg`.
- id: delta_v_mps
  type: Float64
  units: n/a
  required: true
  description: Positional argument `delta_v_mps`.
- id: thrust_n
  type: Float64
  units: n/a
  required: true
  description: Positional argument `thrust_n`.
- id: isp_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `isp_s`.
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
  type: Float64
  units: n/a
  description: Return value of `_constant_thrust_burn_duration_s`.
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

# _constant_thrust_burn_duration_s

## Purpose
Computes how long a constant-thrust burn must last to deliver a commanded delta-v from a given wet mass, in seconds.

## Theory & Math
With $v_e = I_{sp} g_0$, $g_0 = 9.80665\ \mathrm{m/s^2}$, the propellant mass fraction consumed is $\mu = 1 - e^{-\Delta v / v_e}$, and constant thrust $F$ expels propellant at $\dot m = F / v_e$, so the burn duration is $t_b = \dfrac{m\,v_e}{F}\,\mu = \dfrac{m\,v_e}{F}\left(1 - e^{-\Delta v/v_e}\right)$, with $m$ the initial mass in kg and $\Delta v$ in m/s.

## Design & Implementation
Applies four sequential guards before any arithmetic: non-finite or non-positive `mass_kg` returns `NaN`; non-finite or negative `delta_v_mps` returns `NaN`; an exactly zero delta-v short-circuits to `0.0`; then non-finite or non-positive `thrust_n` and `isp_s` each return `NaN`. The standard gravity constant is re-declared locally as `g0 = 9.80665` rather than reusing the file-level `_STANDARD_GRAVITY_MPS2`. `calcControlEffect!` centres the burn window symmetrically about predicted apoapsis using the returned duration.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mass_kg` | Float64 | n/a | yes | Positional argument `mass_kg`. |
| in | `delta_v_mps` | Float64 | n/a | yes | Positional argument `delta_v_mps`. |
| in | `thrust_n` | Float64 | n/a | yes | Positional argument `thrust_n`. |
| in | `isp_s` | Float64 | n/a | yes | Positional argument `isp_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_constant_thrust_burn_duration_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:523-523`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Duplicating the gravity constant locally means the two definitions can drift apart under edit. The model is a single continuous constant-thrust arc: it captures no ignition or shutdown transient, no minimum impulse bit, and no throttling. Because the vehicle keeps thrusting while it moves along its orbit, centring a long burn on apoapsis still incurs finite-burn gravity losses that this duration does not account for.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 249.
