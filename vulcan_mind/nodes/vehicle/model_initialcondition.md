---
id: vehicle.model_initialcondition
label: InitialCondition
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: InitialCondition
  lines:
  - 17
  - 17
inputs:
- id: a
  type: Float64
  units: n/a
  required: true
  description: Field `a`.
- id: e
  type: Float64
  units: n/a
  required: true
  description: Field `e`.
- id: i
  type: Float64
  units: n/a
  required: true
  description: Field `i`.
- id: omega
  type: Float64
  units: n/a
  required: true
  description: Field `ω`.
- id: Omega
  type: Float64
  units: n/a
  required: true
  description: Field `Ω`.
- id: nu
  type: Float64
  units: n/a
  required: true
  description: Field `ν`.
- id: q
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `q`.
- id: ang_vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `ang_vel`.
- id: a_2
  type: Float64,
  units: n/a
  required: true
  description: Field `a`.
- id: e_2
  type: Float64,
  units: n/a
  required: true
  description: Field `e`.
- id: i_2
  type: Float64,
  units: n/a
  required: true
  description: Field `i`.
- id: omega_2
  type: Float64,
  units: n/a
  required: true
  description: Field `ω`.
- id: Omega_2
  type: Float64,
  units: n/a
  required: true
  description: Field `Ω`.
- id: nu_2
  type: Float64,
  units: n/a
  required: true
  description: Field `ν`.
- id: q_2
  type: SVector{4, Float64},
  units: n/a
  required: true
  description: Field `q`.
- id: ang_vel_2
  type: SVector{3, Float64},
  units: n/a
  required: true
  description: Field `ang_vel`.
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
  type: InitialCondition
  units: n/a
  description: Constructed `InitialCondition`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# InitialCondition

## Purpose
Keplerian initial condition for a spacecraft: semimajor axis, eccentricity, inclination, argument of periapsis, RAAN and true anomaly, plus initial attitude quaternion and body angular velocity, with constructors accepting degrees, apsis radii, or apsis altitudes above an oblate planet.

## Design & Implementation
Fields: `a` (m), `e`, `i`, `ω`, `Ω`, `ν` (rad), `q::SVector{4,Float64}` in `(x, y, z, w)` order, `ang_vel::SVector{3,Float64}` (rad/s). The inner constructor takes radians and a `Val(:radians)` tag. The positional outer constructor accepts degrees for the four angles and converts with `deg2rad`. The keyword constructor accepts either `a`/`e` or an `ra`/`rp` radius pair (throwing `ArgumentError` if only one apsis is given), computing `a = (ra+rp)/2`, `e = (ra-rp)/(ra+rp)`, with `ν` defaulting to 180 deg (apoapsis) for apsis input and 0 deg otherwise. The `InitialCondition(planet; ra, hp, ...)` method treats `ra` and `hp` as altitudes above the oblate ellipsoid: it resolves `L_PI`, rotates the apsis directions into the planet frame, bisects for the radii, requires `ra_radius > rp_radius`, and delegates to the keyword constructor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Float64 | n/a | yes | Field `a`. |
| in | `e` | Float64 | n/a | yes | Field `e`. |
| in | `i` | Float64 | n/a | yes | Field `i`. |
| in | `omega` | Float64 | n/a | yes | Field `ω`. |
| in | `Omega` | Float64 | n/a | yes | Field `Ω`. |
| in | `nu` | Float64 | n/a | yes | Field `ν`. |
| in | `q` | SVector{4, Float64} | n/a | yes | Field `q`. |
| in | `ang_vel` | SVector{3, Float64} | n/a | yes | Field `ang_vel`. |
| in | `a_2` | Float64, | n/a | yes | Field `a`. |
| in | `e_2` | Float64, | n/a | yes | Field `e`. |
| in | `i_2` | Float64, | n/a | yes | Field `i`. |
| in | `omega_2` | Float64, | n/a | yes | Field `ω`. |
| in | `Omega_2` | Float64, | n/a | yes | Field `Ω`. |
| in | `nu_2` | Float64, | n/a | yes | Field `ν`. |
| in | `q_2` | SVector{4, Float64}, | n/a | yes | Field `q`. |
| in | `ang_vel_2` | SVector{3, Float64}, | n/a | yes | Field `ang_vel`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | InitialCondition | n/a | — | Constructed `InitialCondition`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__initial_condition_from_time_aligned_telemetry|_initial_condition_from_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:124-124`
- [[analysis.scenario_builders__scenario_initial_condition|_scenario_initial_condition]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:538-538`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`
- [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callees` → `callers` · call · `src/precompile_workload.jl:10-10`
- [[vehicle.model__initial_condition_radius_for_oblate_altitude|_initial_condition_radius_for_oblate_altitude]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:172-172`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/vehicle/spacecraft/model.jl:53-53`
<!-- vulcan:connections:end -->

## Limitations
No range checks on `e` (negative or >= 1 accepted) or on `a > 0` in the plain constructors. The oblate constructor depends on `_initial_condition_oblate_altitude`, whose formula deviates from the standard geodetic altitude. Mixing conventions (angles in degrees for outer constructors, radians inside) is easy to get wrong when calling the `Val(:radians)` inner constructor directly. The oblate constructor requires `planet.Rp_e` and `planet.Rp_p` and may trigger SPICE loading.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 17.
