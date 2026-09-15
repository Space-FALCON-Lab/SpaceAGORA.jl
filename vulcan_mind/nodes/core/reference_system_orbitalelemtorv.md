---
id: core.reference_system_orbitalelemtorv
label: orbitalelemtorv
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: orbitalelemtorv
  lines:
  - 122
  - 122
inputs:
- id: oe
  type: SVector{7, Float64}
  units: n/a
  required: true
  description: Positional argument `oe`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  description: Return value of `orbitalelemtorv`. Returns `collect(R), collect(V)`.
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

# orbitalelemtorv

## Purpose
Converts classical orbital elements to an inertial position and velocity, the standard first step in building an initial condition.

## Theory & Math
$$
p = a(1 - e^2),\quad h = \sqrt{\mu p},\quad r_{pf} = \frac{h^2}{\mu}\frac{1}{1 + e\cos\nu}\begin{bmatrix}\cos\nu\\ \sin\nu\\ 0\end{bmatrix},\quad v_{pf} = \frac{\mu}{h}\begin{bmatrix}-\sin\nu\\ e + \cos\nu\\ 0\end{bmatrix}
$$

then $r = Q^\top r_{pf}$, $v = Q^\top v_{pf}$ with $Q = R_3(\omega) R_1(i) R_3(\Omega)$.

## Design & Implementation
Two methods. The `SVector{7}` form reads `a, e, i, Ω, ω, ν`, computes the semi-latus rectum and specific angular momentum, forms the perifocal position and velocity, builds the perifocal-to-inertial rotation `Q` from the three angles, and returns `Q' * r` and `Q' * v` as plain vectors. The second method unpacks an `InitialCondition` struct into the seven-vector with a zero seventh slot.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `oe` | SVector{7, Float64} | n/a | yes | Positional argument `oe`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `orbitalelemtorv`. Returns `collect(R), collect(V)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__initial_condition_in_j2000|_initial_condition_in_j2000]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:514-514`
- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:117-117`
- [[gnc.targeting_solver__target_planning_impl|_target_planning_impl]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:17-17`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:25-25`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:25-25`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:26-26`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:25-25`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:117-117`
- [[simulation.dynamics_rhs_build_initial_conditions|build_initial_conditions]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2330-2330`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:103-103`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Returns `collect`ed `Vector{Float64}` rather than static vectors, so callers on a hot path pay allocations; hyperbolic orbits with `e > 1` produce a negative `p` and a `sqrt` of a negative number.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 122.
