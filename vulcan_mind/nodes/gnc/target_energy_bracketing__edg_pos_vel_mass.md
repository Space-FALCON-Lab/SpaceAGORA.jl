---
id: gnc.target_energy_bracketing__edg_pos_vel_mass
label: _edg_pos_vel_mass
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_pos_vel_mass
  lines:
  - 184
  - 184
inputs:
- id: sc
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc`.
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
  description: Return value of `_edg_pos_vel_mass`. Returns `pos, vel, mass`.
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

# _edg_pos_vel_mass

## Purpose
`_edg_pos_vel_mass` pulls inertial position, velocity and mass out of a per-spacecraft state object regardless of whether it exposes named fields (`pos`, `vel`, `mass`) or is a raw numeric vector with the conventional layout `[x, y, z, vx, vy, vz, m, ...]`. It feeds the bracket propagation with `SVector`s.

## Design & Implementation
Declared `@inline` with signature `(sc)`. For each quantity it uses `hasproperty` to choose between the named field and positional indexing: `pos = SVector{3,Float64}(sc.pos)` or `SVector{3,Float64}(sc[1], sc[2], sc[3])`; `vel` likewise from `sc.vel` or indices 4:6; `mass = Float64(sc.mass)` or, if `length(sc) >= 7`, `Float64(sc[7])`, else `NaN`. It returns the tuple `(pos, vel, mass)` in SI units (m, m/s, kg). Conversion to `SVector{3,Float64}` copies the three components, so the result does not alias the state.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc` | Any | n/a | yes | Positional argument `sc`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_pos_vel_mass`. Returns `pos, vel, mass`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.target_energy_bracketing__edg_run_target_energy_bracketing_bang|_edg_run_target_energy_bracketing!]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:282-282`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:187-187`
<!-- vulcan:connections:end -->

## Limitations
A vector shorter than six elements raises a `BoundsError` rather than a descriptive error. A `NaN` mass is returned silently for six-element states, and downstream ballistic-coefficient calculations then produce `NaN` energies. Units are assumed, not verified; a normalised state would be misinterpreted. The `hasproperty` checks are runtime and not type-stable when the state type varies.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 184.
