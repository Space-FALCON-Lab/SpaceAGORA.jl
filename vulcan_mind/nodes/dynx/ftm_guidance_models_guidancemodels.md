---
id: dynx.ftm_guidance_models_guidancemodels
label: GuidanceModels
kind: struct
source:
  file: src/dynamics/coupled/force_torque_models/guidance_models.jl
  symbol: GuidanceModels
  lines:
  - 1
  - 5
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors parent namespace under which the guidance effector
    aliases are re-exported.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: guidance_types
  type: Module
  units: n/a
  description: Namespace re-exporting AerobrakingCampaignPropulsiveManeuverGuidanceModel
    and ApoapsisTargetPeriapsisRaiseGuidanceModel as effector-visible types.
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

# GuidanceModels

## Purpose
`GuidanceModels` is a re-export shim inside the effector namespace. It pulls the two campaign-level guidance types out of the top-level guidance package and makes them visible to code that constructs effector lists, so a scenario can attach a guidance model with the same syntax used for gravity or aerodynamic effectors without reaching across package boundaries.

## Theory & Math
The two re-exported types implement impulsive-approximation manoeuvre guidance. A periapsis raise executed at apoapsis costs, to first order,

$$\Delta v = \sqrt{\frac{2\mu}{r_a} - \frac{\mu}{a_{new}}} - \sqrt{\frac{2\mu}{r_a} - \frac{\mu}{a_{old}}}$$

with $\mu$ in m^3/s^2, $r_a$ the apoapsis radius in m and $a$ the semi-major axis in m, where $a_{new} = (r_a + r_{p,target})/2$. The finite-burn realisation converts this to a thrust duration through the rocket equation

$$\Delta v = I_{sp} g_0 \ln\!\left(\frac{m_0}{m_f}\right), \qquad \dot{m} = -\frac{T}{I_{sp} g_0}$$

with $I_{sp}$ the specific impulse in s, $g_0 = 9.80665$ m/s^2, $T$ the thrust in N, and $m$ the wet mass in kg. Aerobraking campaign guidance uses the same relations to schedule corridor-control burns that hold the periapsis density within a heat-rate limit.

## Model & Assumptions
The shim assumes the underlying guidance implementations already handle their own state and epoch bookkeeping; it adds no behaviour and no dispatch. Impulsive sizing is valid when the burn arc is short compared with the orbital period, roughly under a few percent, otherwise gravity losses accumulate and the achieved periapsis undershoots the target.

## Design & Implementation
The whole module body is a `using ...GuidanceModels: ...` line plus a matching `export`. Naming the submodule identically to the parent package is deliberate: effector-side code refers to `GuidanceModels.ApoapsisTargetPeriapsisRaiseGuidanceModel` and resolves it inside the effector namespace, while the definitions stay single-sourced upstream.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors parent namespace under which the guidance effector aliases are re-exported. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `guidance_types` | Module | n/a | — | Namespace re-exporting AerobrakingCampaignPropulsiveManeuverGuidanceModel and ApoapsisTargetPeriapsisRaiseGuidanceModel as effector-visible types. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models/guidance_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because it re-exports rather than defines, any change to the upstream constructor signatures propagates here with no local compatibility layer. The shim does not extend `calcForceTorque`, so guidance types must be wired to a thruster effector to produce an actual wrench. Only two of the guidance models are re-exported; others remain reachable only through the parent package.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models/guidance_models.jl:1-5`.
