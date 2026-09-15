---
id: dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench
label: _aero_pure_wrench
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _aero_pure_wrench
  lines:
  - 356
  - 498
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace supplying the AerodynamicEffectors sampling
    contract.
- id: state_sample
  type: StateSample
  units: m,m/s,rad/s
  required: true
  description: Inertial position, inertial velocity and body attitude of the spacecraft
    at the evaluation epoch.
- id: environment_sample
  type: EnvironmentSample
  units: kg/m^3,m/s,K
  required: true
  description: Free-stream density, planet-relative wind velocity and translational
    temperature from the atmosphere model.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: force_ii
  type: SVector{3,Float64}
  units: N
  description: Total aerodynamic force resolved in inertial axes.
- id: torque_body
  type: SVector{3,Float64}
  units: N*m
  description: Total aerodynamic torque about the spacecraft mass centre in body axes.
- id: component_accels
  type: NTuple{3,SVector{3,Float64}}
  units: N
  description: Drag, lift and cross-wind force components retained for telemetry and
    cache storage.
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

# _aero_pure_wrench

## Purpose
`_aero_pure_wrench` is the allocation-free kernel that turns a state sample and an atmosphere sample into the aerodynamic force and torque acting on a multi-body spacecraft. Every public aerodynamic effector method in this file funnels into it, so drag, lift and cross-wind bookkeeping exists in exactly one place. It walks the body tree, evaluates a panel or coefficient model per link, rotates each link contribution into inertial axes and accumulates the resulting wrench about the spacecraft mass centre.

## Theory & Math
The free-stream dynamic pressure is

$$q_\infty = \tfrac{1}{2}\,\rho\,\lVert \vec{v}_{rel} \rVert^2$$

with $\rho$ the local mass density in kg/m^3 and $\vec{v}_{rel} = \vec{v}_{ii} - \vec{v}_{atm}$ the planet-relative velocity in m/s (the co-rotating atmosphere contributes $\vec{\omega}_p \times \vec{r}$). Each link of reference area $S$ produces

$$\vec{F} = -q_\infty S \left( C_D \hat{v}_{rel} + C_L \hat{l} + C_Y \hat{y} \right), \qquad \vec{\tau}_b = \sum_k \vec{r}_{k/cm} \times \vec{F}_k$$

where $C_D$, $C_L$ and $C_Y$ are dimensionless drag, lift and side-force coefficients, $\hat{v}_{rel}$ is the unit relative-velocity vector, and $\vec{r}_{k/cm}$ is the link centre-of-pressure offset from the spacecraft mass centre in m. In free-molecular flow the coefficient model is written through the speed ratio

$$s = \frac{\lVert \vec{v}_{rel} \rVert}{\sqrt{2 R T}}$$

which is why the file caches $\sqrt{\pi}$ and $1/\sqrt{\pi}$: Schaaf-Chambre free-molecular coefficients contain $\operatorname{erf}(s\sin\alpha)$ and $e^{-s^2\sin^2\alpha}/(s\sqrt{\pi})$ terms, with $\alpha$ the local panel incidence in rad and $R$ the specific gas constant in J/(kg*K).

## Model & Assumptions
The kernel assumes rigid links with fixed reference areas, a quasi-steady flow field (no unsteady wake memory), and additive per-link wrenches with no shadowing between panels unless the coefficient model itself encodes it. Density is sampled once for the whole spacecraft unless per-link atmosphere queries are enabled, which matters only when link offsets are a non-trivial fraction of the density scale height (roughly 7 km at Earth, 11 km at Mars). Validity spans free-molecular altitudes above roughly 150 km down through the transitional regime where the supplied coefficient tables remain calibrated.

## Design & Implementation
The function is dispatched on a `Symbol` coefficient mode (`:constant`, `:fM`, `:no_ballistic_flight`) so the three public `calcForceTorque` methods share one body-tree traversal. Link wrenches are accumulated with `StaticArrays` values to avoid heap traffic inside the ODE right-hand side, and `collect_and_reset_link_wrenches!` drains the per-link buffer between calls. An optional `link_atmosphere_fn` closure lets callers that hold the live density model resolve density per link; callers without a model pass `nothing` and inherit the single-sample behaviour.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace supplying the AerodynamicEffectors sampling contract. |
| in | `state_sample` | StateSample | m,m/s,rad/s | yes | Inertial position, inertial velocity and body attitude of the spacecraft at the evaluation epoch. |
| in | `environment_sample` | EnvironmentSample | kg/m^3,m/s,K | yes | Free-stream density, planet-relative wind velocity and translational temperature from the atmosphere model. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `force_ii` | SVector{3,Float64} | N | — | Total aerodynamic force resolved in inertial axes. |
| out | `torque_body` | SVector{3,Float64} | N*m | — | Total aerodynamic torque about the spacecraft mass centre in body axes. |
| out | `component_accels` | NTuple{3,SVector{3,Float64}} | N | — | Drag, lift and cross-wind force components retained for telemetry and cache storage. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:518-518`
- [[dynamics.aerodynamic_wrench_models_wrench_caching_bang|wrench_caching!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:552-552`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:426-426`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:391-391`
- `callees` → [[dynamics.aerodynamic_wrench_models__aero_link_angles|_aero_link_angles]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:419-419`
- `callees` → [[dynamics.aerodynamic_wrench_models__aero_link_area|_aero_link_area]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:420-420`
- `callees` → [[dynamics.aerodynamic_wrench_models__constant_drag_coefficient|_constant_drag_coefficient]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:461-461`
- `callees` → [[dynamics.aerodynamic_wrench_models__fold_constant_incidence|_fold_constant_incidence]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:461-461`
- `callees` → [[dynamics.aerodynamic_wrench_models__validate_fm_incidence|_validate_fm_incidence]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:363-363`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:458-458`
<!-- vulcan:connections:end -->

## Limitations
Panel shadowing, wake interference and rarefied-transitional bridging are not modelled: coefficients are taken at face value from the supplied model. The kernel does not validate that the density sample is finite or positive, so a failed atmosphere query silently yields a zero wrench. Torque accuracy degrades when the centre-of-pressure offsets are stale relative to a deploying or articulating structure, because link geometry is read once per call rather than re-derived from the live multibody state.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl:356-498`, with the calling `calcForceTorque` methods at lines 588, 699 and 894 of the same file.
