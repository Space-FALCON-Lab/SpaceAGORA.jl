---
id: parcore.reference_system_rvtoorbitalelement
label: rvtoorbitalelement
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: rvtoorbitalelement
  lines:
  - 225
  - 233
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: SimulationModel namespace supplying the planet model type and the StaticArrays
    state representation used by the conversion.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: orbital_elements
  type: SVector{6,Float64} or SVector{7,Float64}
  units: m, -, rad
  description: Classical orbital element vector (a, e, i, RAAN, argument of periapsis,
    true anomaly), optionally carrying spacecraft mass as a seventh entry.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- parcore
origin: agent
---

# rvtoorbitalelement

## Purpose
`rvtoorbitalelement` converts an inertial Cartesian state (position and velocity) into classical Keplerian orbital elements for a given planet. Two methods are defined in this file: one appends the spacecraft mass to produce a seven-element vector consumed by the aerobraking state pipeline, and one returns the bare six-element set. Both delegate the arithmetic to `_rvtoorbitalelement_core`, which owns the singularity handling, so the exported methods only shape the return type.

## Theory & Math
The conversion uses the standard two-body relations. The specific angular momentum and eccentricity vectors are

$$\vec{h} = \vec{r} \times \vec{v}, \qquad \vec{e} = \frac{1}{\mu}\left[\left(v^2 - \frac{\mu}{r}\right)\vec{r} - (\vec{r}\cdot\vec{v})\,\vec{v}\right],$$

with semi-major axis from the vis-viva energy relation

$$a = \left(\frac{2}{r} - \frac{v^2}{\mu}\right)^{-1},$$

inclination from $\cos i = h_z/\lVert\vec{h}\rVert$, and the node vector $\vec{n} = \hat{z} \times \vec{h}$ fixing the right ascension of the ascending node. True anomaly follows from $\cos\nu = (\vec{e}\cdot\vec{r})/(e\,r)$ with the quadrant resolved by the sign of $\vec{r}\cdot\vec{v}$.

## Model & Assumptions
The transformation assumes a point-mass central body with gravitational parameter taken from the planet model, and that the input state is expressed in the same inertial frame the planet parameter describes. It is a purely instantaneous, unperturbed mapping: oblateness, drag and third-body terms are not inverted out of the state. Circular and equatorial orbits are geometrically degenerate because the eccentricity and node vectors collapse, so the core routine substitutes reference directions in those cases.

## Design & Implementation
The file keeps element conversion, latitude and longitude conversion, LVLH quaternion construction and relative RTN transformations together as the reference-frame interface layer. `rvtoorbitalelement` is written against `SVector` inputs so it allocates no heap memory inside the integrator right-hand side. The mass-carrying method exists because downstream aerobraking configuration records store the propellant-dependent mass alongside the geometry in a single static vector, which avoids a second tuple unpack in the callback path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | SimulationModel namespace supplying the planet model type and the StaticArrays state representation used by the conversion. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `orbital_elements` | SVector{6,Float64} or SVector{7,Float64} | m, -, rad | — | Classical orbital element vector (a, e, i, RAAN, argument of periapsis, true anomaly), optionally carrying spacecraft mass as a seventh entry. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system_latlongtooe|latlongtoOE]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:347-347`
- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:161-161`
- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:92-92`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:104-104`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:121-121`
- [[gnc.heat_load_control__edg_closed_form_heat_load_trajectory|_edg_closed_form_heat_load_trajectory]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:148-148`
- [[gnc.heat_load_control__edg_drag_passage_duration|_edg_drag_passage_duration]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:89-89`
- [[gnc.targeting_control__edg_orbit_metrics_from_rv|_edg_orbit_metrics_from_rv]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:299-299`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:90-90`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:92-92`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:104-104`
- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:502-502`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:121-121`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:90-90`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:92-92`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:161-161`
- [[simulation.save_fields__save_periapsis_altitude|_save_periapsis_altitude]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:79-79`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:80-80`
- [[simulation.targeting__gram_orbit_period_target|_gram_orbit_period_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:151-151`
- [[simulation.targeting__gram_periapsis_target|_gram_periapsis_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:119-119`

**Downstream**

- `callees` → [[core.reference_system__rvtoorbitalelement_core|_rvtoorbitalelement_core]] · `callers` · call · `src/core/interfaces/reference_system.jl:226-226`
<!-- vulcan:connections:end -->

## Limitations
Accuracy degrades near-parabolic, where the semi-major axis diverges, and the routine does not branch to universal-variable elements. Equatorial and circular orbits return convention-dependent angles rather than raising an error, so consumers must not compare RAAN or argument of periapsis across such cases. The gravitational parameter is read from the planet argument, so mixing a state expressed about one body with another body's model silently produces wrong elements.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl:225-233`.
