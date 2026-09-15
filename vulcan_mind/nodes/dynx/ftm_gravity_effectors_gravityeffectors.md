---
id: dynx.ftm_gravity_effectors_gravityeffectors
label: GravityEffectors
kind: struct
source:
  file: src/dynamics/coupled/force_torque_models/gravity_effectors.jl
  symbol: GravityEffectors
  lines:
  - 1
  - 13
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors parent namespace supplying the calcForceTorque, wrench
    and environment_requirements generics extended here.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: gravity_models
  type: Module
  units: n/a
  description: Namespace exporting ConstantGravityModel, InverseSquaredGravityModel,
    InverseSquaredJ2GravityModel, aerobraking_gravity_force_ii and j2_secular_rates.
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

# GravityEffectors

## Purpose
`GravityEffectors` scopes the three baseline gravity effectors used by the coupled dynamics right-hand side. It imports the parent effector generics plus the gravity-backbone sampling hooks, then includes the environment gravity model file so the constant, inverse-square and J2 models are defined in one namespace and share the effector method tables.

## Theory & Math
The exported models implement, in increasing fidelity,

$$\vec{a} = \vec{g}_0, \qquad \vec{a} = -\frac{\mu}{r^{3}}\vec{r}, \qquad \vec{a} = -\frac{\mu}{r^{3}}\vec{r} + \vec{a}_{J_2}$$

where $\mu$ is the gravitational parameter in m^3/s^2, $\vec{r}$ the planet-centred inertial position in m and $r = \lVert \vec{r} \rVert$. The oblateness term has components

$$a_{J_2,x} = -\frac{3 J_2 \mu R_e^{2}}{2 r^{5}} x \left(1 - \frac{5 z^{2}}{r^{2}}\right), \quad a_{J_2,z} = -\frac{3 J_2 \mu R_e^{2}}{2 r^{5}} z \left(3 - \frac{5 z^{2}}{r^{2}}\right)$$

with $J_2$ the dimensionless second zonal harmonic (1.0826e-3 at Earth, 1.9555e-3 at Mars) and $R_e$ the equatorial radius in m. The exported `j2_secular_rates` returns the induced mean rates $\dot{\Omega} = -\tfrac{3}{2} n J_2 (R_e/p)^2 \cos i$ and $\dot{\omega} = \tfrac{3}{4} n J_2 (R_e/p)^2 (5\cos^2 i - 1)$ in rad/s, with $n$ the mean motion, $p$ the semi-latus rectum in m and $i$ the inclination in rad.

## Model & Assumptions
Point-mass and J2 models assume the spacecraft is outside the planet's Brillouin sphere and treat the vehicle as a point mass, so gravity-gradient torque is not produced here. The constant model is only valid over spans short enough that $\mu/r^2$ barely changes, typically bench tests and near-surface scenarios. J2 truncation leaves the tesseral and higher zonal terms to the harmonics effector in `perturbations.jl`.

## Design & Implementation
The submodule imports `gravity_backbone_structure` and `gravity_backbone_acceleration_ii` from `EffectorSampling` so a solver can split gravity into a stiff analytic backbone integrated separately from the perturbing wrenches. Including `environment/gravity/gravity_models.jl` from here keeps a single definition of the models while letting the effector-side methods live under the parent generics.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors parent namespace supplying the calcForceTorque, wrench and environment_requirements generics extended here. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `gravity_models` | Module | n/a | — | Namespace exporting ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel, aerobraking_gravity_force_ii and j2_secular_rates. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models/gravity_effectors.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No self-gravity, no gravity-gradient torque, and no relativistic correction are included; those belong to the harmonics and perturbation effectors. The models do not guard against $r \to 0$, so a state that passes through the planet centre produces a non-finite acceleration. Secular J2 rates assume small eccentricity and are inaccurate for near-critical inclinations close to 63.4 degrees where $5\cos^2 i - 1$ vanishes.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models/gravity_effectors.jl:1-13` and the included `src/environment/gravity/gravity_models.jl`.
