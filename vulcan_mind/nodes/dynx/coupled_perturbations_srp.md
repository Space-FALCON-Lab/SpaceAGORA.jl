---
id: dynx.coupled_perturbations_srp
label: srp
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: srp
  lines:
  - 1257
  - 1281
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace providing the SolarRadiationPressureModel
    dispatch context.
- id: orbit_state
  type: SVector{3,Float64}
  units: m
  required: true
  description: Planet-centred inertial position of the spacecraft.
- id: epoch_et
  type: Float64
  units: s
  required: true
  description: SPICE ephemeris time used to resolve the primary-body-to-Sun vector.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: accel_ii
  type: SVector{3,Float64}
  units: m/s^2
  description: Cannonball solar radiation pressure acceleration in inertial axes,
    including the eclipse shadow factor.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynx
origin: agent
---

# srp

## Purpose
`srp` resolves the Sun's position relative to the primary body through SPICE and then evaluates the cannonball solar radiation pressure acceleration for the spacecraft. It is the ephemeris-aware wrapper around `srp_cannonball_accel`, keeping the SPICE query in one place so the pure geometric kernel stays testable without kernel files loaded.

## Theory & Math
The cannonball acceleration is

$$\vec{a}_{SRP} = -\nu \, P_{\odot} \left(1 + \varrho\right) \frac{A}{m} \left(\frac{\mathrm{AU}}{d}\right)^{2} \hat{u}_{sc\to\odot}$$

where $P_{\odot} \approx 4.56\times 10^{-6}$ N/m^2 is the solar radiation pressure at one astronomical unit, $\varrho$ is the dimensionless reflection coefficient (0 for a perfect absorber, 1 for a perfect specular reflector), $A$ is the reference area in m^2, $m$ is the spacecraft mass in kg, $d$ is the spacecraft-Sun distance in m, $\mathrm{AU} = 1.495978707\times 10^{11}$ m is the astronomical unit hard-coded as the `AU_m` keyword default, and $\nu \in [0,1]$ is the eclipse shadow factor. The inverse-square factor follows from the solar irradiance falling as

$$P(d) = P_{\odot}\left(\frac{\mathrm{AU}}{d}\right)^{2}$$

and $\nu$ is computed from the illuminated fraction of the solar disc using the planetary radius $R_p$ in m, giving $\nu = 1$ in full sunlight, $\nu = 0$ in umbra, and an area-ratio value in penumbra.

## Model & Assumptions
The model treats the spacecraft as a sphere of uniform optical properties, so attitude does not change the SRP force and no SRP torque is produced. It assumes the Sun is a source at the SPICE-reported barycentric position and ignores light-time correction at the spacecraft. Accuracy is around ten percent for compact bodies and degrades sharply for large flexible arrays, where a panel model is required.

## Design & Implementation
The wrapper converts the incoming position to an `SVector{3,Float64}`, maps the planet name to its SPICE query name through `_spice_query_name`, and requests the primary-to-Sun vector with `spice_position_j2000_m("sun", et, primary_body_name)`. It then forwards the geometry, the planet equatorial radius `planet.Rp_e`, the unscaled pressure, the reflection coefficient, the area, the mass and `AU_m` to `srp_cannonball_accel`, which applies the eclipse factor from `eclipse_area_calc`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace providing the SolarRadiationPressureModel dispatch context. |
| in | `orbit_state` | SVector{3,Float64} | m | yes | Planet-centred inertial position of the spacecraft. |
| in | `epoch_et` | Float64 | s | yes | SPICE ephemeris time used to resolve the primary-body-to-Sun vector. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `accel_ii` | SVector{3,Float64} | m/s^2 | — | Cannonball solar radiation pressure acceleration in inertial axes, including the eclipse shadow factor. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:224-224`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:244-244`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:256-256`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:306-306`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:224-224`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:244-244`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:256-256`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:306-306`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:224-224`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:256-256`

**Downstream**

- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1268-1268`
- `callees` → [[dynamics.perturbations_srp_cannonball_accel|srp_cannonball_accel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1271-1271`
- `callees` → [[environment.simple_ephemerides_spice_position_j2000_m|spice_position_j2000_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1269-1269`
<!-- vulcan:connections:end -->

## Limitations
Every call performs a SPICE lookup, which must be serialised through the module's SPICE lock and dominates the cost when many spacecraft are propagated. Albedo and planetary infrared pressure are handled by separate functions and are not included here. Self-shadowing between spacecraft, specular versus diffuse split, and thermal re-radiation (Yarkovsky-like) effects are not modelled.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl:1257-1281`, with the dispatching effector method at line 1283 and the eclipse geometry at line 2277 of the same file.
