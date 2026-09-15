---
id: dynx.ftm_perturbation_effectors_perturbationeffectors
label: PerturbationEffectors
kind: struct
source:
  file: src/dynamics/coupled/force_torque_models/perturbation_effectors.jl
  symbol: PerturbationEffectors
  lines:
  - 1
  - 33
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors parent namespace supplying calcForceTorque, wrench
    and environment_requirements for the perturbation models.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: perturbation_models
  type: Module
  units: n/a
  description: Namespace exporting NBodyGravityModel, GravitationalHarmonicsModel,
    SolarRadiationPressureModel, MagneticTorqueRodModel, EddyCurrentDampingModel and
    LVLHCascadeAttitudeControlModel.
- id: debug_flag
  type: Ref{Bool}
  units: n/a
  description: Load-time SPACEAGORA_DEBUG_COMPARE_J2 switch that enables the J2-versus-harmonics
    comparison path.
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

# PerturbationEffectors

## Purpose
`PerturbationEffectors` scopes the high-fidelity orbital and attitude perturbation models: third-body gravity from SPICE ephemerides, spherical-harmonic gravity, solar radiation pressure with albedo and infrared terms, magnetic torque rods, eddy-current damping and an LVLH cascade attitude controller. It wires the SPICE lock, the scratch workspaces and the ephemeris caches before including the 2342-line `perturbations.jl` implementation.

## Theory & Math
Third-body gravity uses the differential form that cancels the common frame acceleration,

$$\vec{a}_{3B} = \mu_3 \left( \frac{\vec{r}_{3/sc}}{\lVert \vec{r}_{3/sc} \rVert^{3}} - \frac{\vec{r}_{3/p}}{\lVert \vec{r}_{3/p} \rVert^{3}} \right)$$

with $\mu_3$ the perturbing body's gravitational parameter in m^3/s^2 and the position vectors in m. Spherical-harmonic gravity is the gradient of

$$U = \frac{\mu}{r}\left[ 1 + \sum_{n=2}^{N}\sum_{m=0}^{n} \left(\frac{R_e}{r}\right)^{n} P_{nm}(\sin\phi)\left(C_{nm}\cos m\lambda + S_{nm}\sin m\lambda\right) \right]$$

with $P_{nm}$ the associated Legendre functions, $\phi$ the geocentric latitude and $\lambda$ the longitude in rad, and $C_{nm}$, $S_{nm}$ the normalised coefficients. Solar radiation pressure follows the cannonball law $\vec{a}_{SRP} = -\nu P_{\odot}(1+\varrho) (A/m) (\mathrm{AU}/d)^2 \hat{u}_{\odot}$, and the magnetic torques follow $\vec{\tau} = \vec{m} \times \vec{B}$ in N*m.

## Model & Assumptions
The module assumes SPICE kernels are loaded and serialised through `SPICE_LOCK`, because the SPICE Fortran core is not thread-safe. Harmonic evaluation assumes the field coefficients are fully normalised and the evaluation point lies outside the Brillouin sphere, where the series converges. Magnetic models assume a dipole or IGRF field valid within the magnetosphere and degrade beyond roughly six planetary radii.

## Design & Implementation
The `using` block imports scratch workspaces (`NBodyScratchWorkspace`, `HarmonicsScratchWorkspace`) so the batched harmonics kernel can run without allocating inside the ODE right-hand side. It deliberately reuses `AerodynamicEffectors.rtolatlong` and `_multibody_thread_decision` instead of re-including `reference_system.jl`, avoiding duplicate module-scoped bindings. The `_DEBUG_COMPARE_J2` environment flag is read once at load time to keep the harmonics hot path free of syscalls.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors parent namespace supplying calcForceTorque, wrench and environment_requirements for the perturbation models. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `perturbation_models` | Module | n/a | — | Namespace exporting NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel, MagneticTorqueRodModel, EddyCurrentDampingModel and LVLHCascadeAttitudeControlModel. |
| out | `debug_flag` | Ref{Bool} | n/a | — | Load-time SPACEAGORA_DEBUG_COMPARE_J2 switch that enables the J2-versus-harmonics comparison path. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models/perturbation_effectors.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Every perturbation that needs an ephemeris takes the global SPICE lock, so heavy third-body configurations serialise across threads. The harmonics series is truncated at the configured degree, and the truncation error grows near the surface where the neglected terms scale as $(R_e/r)^{n}$. Eclipse modelling for SRP is a conical approximation and does not resolve penumbral limb darkening.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models/perturbation_effectors.jl:1-33` and the included `src/dynamics/coupled/perturbations.jl`.
