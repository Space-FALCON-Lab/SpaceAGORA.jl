---
id: dynx.ftm_aerodynamic_effectors_aerodynamiceffectors
label: AerodynamicEffectors
kind: struct
source:
  file: src/dynamics/coupled/force_torque_models/aerodynamic_effectors.jl
  symbol: AerodynamicEffectors
  lines:
  - 1
  - 18
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors parent namespace whose calcForceTorque, wrench and
    environment_requirements generics this submodule extends.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: aero_models
  type: Module
  units: n/a
  description: Namespace exporting AerodynamicCoefficientConstant, AerodynamicCoefficientfM
    and AerodynamicCoefficientNoBallisticFlight together with their wrench methods.
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

# AerodynamicEffectors

## Purpose
`AerodynamicEffectors` is the submodule that scopes the aerodynamic force and torque implementations. It imports the parent effector generics, pulls in the geodetic reference-system helpers, and includes the large `aerodynamic_wrench_models.jl` body exactly once so the coefficient structs and their methods live in a single namespace with a single set of bindings.

## Theory & Math
The submodule carries no equations of its own; it defines the scope in which the aerodynamic wrench relations are evaluated. Those relations reduce to

$$\vec{F}_{aero} = -\tfrac{1}{2}\rho \lVert \vec{v}_{rel} \rVert^2 S \, C_D \, \hat{v}_{rel} + \text{lift and side terms}, \qquad \vec{\tau}_{aero} = \vec{r}_{cp/cm} \times \vec{F}_{aero}$$

with $\rho$ in kg/m^3, $\vec{v}_{rel}$ the planet-relative velocity in m/s, $S$ the reference area in m^2, $C_D$ the dimensionless drag coefficient, and $\vec{r}_{cp/cm}$ the centre-of-pressure offset in m. The three exported coefficient structs differ in how $C_D$, $C_L$ and $C_Y$ are produced: a fixed constant, a free-molecular function of the speed ratio $s = \lVert \vec{v}_{rel} \rVert / \sqrt{2RT}$, and a no-ballistic-flight variant that suppresses the lift contribution.

## Model & Assumptions
The submodule assumes the parent `DynamicEffectors` module has already declared the `calcForceTorque`, `wrench`, `wrench_caching!`, `environment_requirements` and `solver_partition` generics; it imports rather than redefines them, which keeps method tables shared. It also assumes `reference_system.jl` is safe to include here, and other submodules deliberately reuse that binding (`rtolatlong`) instead of re-including the file, because a second include would create a duplicate module-scoped copy.

## Design & Implementation
Ordering matters in the module body: the `using`/`import` block resolves the parent generics before either include runs, so the included method definitions extend the parent functions rather than creating fresh local ones. Only the three coefficient structs are exported; the wrench kernels and cache helpers stay internal. The include path is built with `joinpath(@__DIR__, ...)` so the module loads identically regardless of the process working directory.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors parent namespace whose calcForceTorque, wrench and environment_requirements generics this submodule extends. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `aero_models` | Module | n/a | — | Namespace exporting AerodynamicCoefficientConstant, AerodynamicCoefficientfM and AerodynamicCoefficientNoBallisticFlight together with their wrench methods. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models/aerodynamic_effectors.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the include is textual, every method defined in `aerodynamic_wrench_models.jl` becomes a method of the parent generics, and a symbol collision there is only visible at load time. The submodule provides no runtime guard against an atmosphere model that returns a zero or negative density. Its exports do not include the sampling helpers, so callers that need `wrench` must reach through the parent `DynamicEffectors` namespace.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models/aerodynamic_effectors.jl:1-18` and the file it includes, `src/dynamics/coupled/aerodynamic_wrench_models.jl`.
