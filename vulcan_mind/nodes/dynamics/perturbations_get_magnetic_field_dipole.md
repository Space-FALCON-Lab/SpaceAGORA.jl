---
id: dynamics.perturbations_get_magnetic_field_dipole
label: get_magnetic_field_dipole
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: get_magnetic_field_dipole
  lines:
  - 1787
  - 1787
inputs:
- id: r_ecef
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `r_ecef`.
- id: L_PI
  type: MMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Positional argument `L_PI`.
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
  description: Return value of `get_magnetic_field_dipole`. Returns `L_PI' * B_ecef`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# get_magnetic_field_dipole

## Purpose
Evaluates Earth's tilted-dipole magnetic field at a planet-fixed position and rotates it to inertial.

## Theory & Math
$$
\vec{B} = B_0 \left(\frac{R_E}{r}\right)^3 \left( \hat{m} - 3(\hat{m} \cdot \hat{r})\hat{r} \right)
$$

## Design & Implementation
Normalises the position, takes the cosine of the magnetic colatitude against the fixed dipole axis `M_HAT_ECEF`, evaluates `B0 (R/r)³ (m̂ - 3 cosθ r̂)`, and rotates through `L_PI'`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_ecef` | AbstractVector | n/a | yes | Positional argument `r_ecef`. |
| in | `L_PI` | MMatrix{3, 3, Float64} | n/a | yes | Positional argument `L_PI`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `get_magnetic_field_dipole`. Returns `L_PI' * B_ecef`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__magnetic_field_inertial|_magnetic_field_inertial]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2020-2020`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2058-2058`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1913-1913`
- [[dynamics.perturbations_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1770-1770`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynx.coupled_perturbations_calculate_magnetic_torque|calculate_magnetic_torque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1808-1808`
<!-- vulcan:connections:end -->

## Limitations
Uses 2020 Earth dipole constants regardless of planet; it is wrong for Mars or Venus.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1787.
