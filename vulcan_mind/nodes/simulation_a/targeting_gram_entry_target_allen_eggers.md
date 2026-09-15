---
id: simulation_a.targeting_gram_entry_target_allen_eggers
label: _gram_entry_target_allen_eggers
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_entry_target_allen_eggers
  lines:
  - 224
  - 313
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: entry_state
  type: Tuple{SVector{3,Float64}, SVector{3,Float64}, Float64, Float64, Float64}
  units: m, m/s, s, kg, m^2
  required: true
  description: Inertial position and velocity, the propagation interval, and the spacecraft
    mass and aerodynamic reference area used to form the ballistic coefficient.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: entry_target
  type: Union{Nothing, Tuple{Float64, Float64, Float64}}
  units: m, rad, rad
  description: Predicted altitude, latitude and longitude at the end of the interval,
    or `nothing` when the inputs or the planet atmosphere reference are unusable.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# _gram_entry_target_allen_eggers

## Purpose
`_gram_entry_target_allen_eggers` predicts where a decelerating vehicle will be at the end of a track-cache segment during atmospheric flight. Keplerian targeting is useless here because drag dominates, so the refresh path substitutes this drag-aware predictor whenever the entry target mode is enabled and the vehicle is inside the atmospheric band.

## Theory & Math
The predictor integrates three-degree-of-freedom point-mass entry over a rotating spherical planet. With ballistic coefficient $\beta = m/(C_D A)$ and exponential reference density $\rho(h) = \rho_{\mathrm{ref}} \exp\!\left((h_{\mathrm{ref}} - h)/H\right)$, the drag deceleration is $a_D = \tfrac{1}{2}\rho v^2/\beta$ and the state advances by explicit Euler steps of size $\Delta\tau$ on

$$\dot v = -a_D - g\sin\gamma, \qquad \dot\gamma = \left(\frac{v}{r} - \frac{g}{v}\right)\cos\gamma, \qquad \dot\chi = \frac{v\cos\gamma\sin\chi\tan\varphi}{r}$$

$$\dot h = v\sin\gamma, \qquad \dot\varphi = \frac{v\cos\gamma\cos\chi}{r}, \qquad \dot\lambda = \frac{v\cos\gamma\sin\chi}{r\cos\varphi} - \omega_p$$

with $g = \mu/r^2$ and $r = R_{p,e} + h$. Flight-path angle $\gamma$ and heading $\chi$ are initialised from the north, east and up velocity components in the local NED frame, $\gamma = \operatorname{atan2}(v_U, v_h)$ and $\chi = \operatorname{atan2}(v_E, v_N)$. This is the Allen-Eggers formulation without the closed-form ballistic solution: keeping the differential form admits the Coriolis-free rotating-planet longitude rate and a non-constant $\gamma$.

## Model & Assumptions
The model treats the vehicle as a point mass with a constant drag coefficient and no lift, flying over a spherical rotating planet with an exponential atmosphere characterised by the planet's reference density, reference altitude and scale height. Because the result only seeds a cache endpoint rather than the simulated trajectory, this fidelity is sufficient: an error in the predicted endpoint costs extra samples or an earlier cache miss, not a wrong physical answer. Mass defaults to the summed dry and propellant mass when a live value is unavailable, and to 100 kg if even that lookup fails.

## Design & Implementation
Input validation rejects non-finite or non-positive intervals, masses, areas and reference densities up front by returning `nothing`, so the caller falls back to solver-endpoint propagation. State is rotated into the planet frame with `r_intor_p!` and converted by `rtolatlong` and `latlongtoNED`. The step count is `ceil(dt / _gram_entry_target_dt())` clamped between two and the configured maximum. Inside the loop, flight-path angle is clamped to plus or minus 89 degrees and latitude to plus or minus 89.9 degrees to keep the `tan` and `sec` terms finite near the poles, speed is floored at one millimetre per second, and heading and longitude are wrapped into the principal interval with conditional additions rather than trigonometric renormalisation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `entry_state` | Tuple{SVector{3,Float64}, SVector{3,Float64}, Float64, Float64, Float64} | m, m/s, s, kg, m^2 | yes | Inertial position and velocity, the propagation interval, and the spacecraft mass and aerodynamic reference area used to form the ballistic coefficient. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `entry_target` | Union{Nothing, Tuple{Float64, Float64, Float64}} | m, rad, rad | — | Predicted altitude, latitude and longitude at the end of the interval, or `nothing` when the inputs or the planet atmosphere reference are unusable. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:108-108`

**Downstream**

- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:242-242`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:240-240`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:241-241`
- `callees` → [[simulation.config__gram_entry_target_cd|_gram_entry_target_cd]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:231-231`
- `callees` → [[simulation.config__gram_entry_target_dt|_gram_entry_target_dt]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:260-260`
- `callees` → [[simulation.config__gram_entry_target_max_steps|_gram_entry_target_max_steps]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:261-261`
- `callees` → [[simulation.targeting__gram_entry_reference_density|_gram_entry_reference_density]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:277-277`
<!-- vulcan:connections:end -->

## Limitations
Explicit Euler with a default half-second step accumulates error over long intervals, and the fixed drag coefficient ignores Mach and angle-of-attack dependence, so predictions degrade for lifting or high-angle entries. The exponential atmosphere is a single-scale-height approximation that will not match the GRAM profile the cache subsequently samples. Polar clamping means trajectories passing very near a pole return a bounded but biased latitude.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl:223-313`.
