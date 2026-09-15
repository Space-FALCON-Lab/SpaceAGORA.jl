---
id: dynamics.aerodynamic_wrench_models__aero_link_angles
label: _aero_link_angles
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _aero_link_angles
  lines:
  - 285
  - 285
inputs:
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: body
  type: Any
  units: n/a
  required: true
  description: Positional argument `body`.
- id: root_index
  type: Int
  units: n/a
  required: true
  description: Positional argument `root_index`.
- id: orientation_sim
  type: Bool
  units: n/a
  required: true
  description: Positional argument `orientation_sim`.
- id: vel_pi
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel_pi`.
- id: theta_body
  type: Float64
  units: n/a
  required: true
  description: Positional argument `θ_body`.
- id: fixed_attitude_incidence
  type: Symbol
  units: n/a
  required: false
  description: Positional argument `fixed_attitude_incidence` (default `:max_drag`).
- id: q_root_ib
  type: Union{Nothing, SVector{4, Float64}}
  units: n/a
  required: false
  description: Positional argument `q_root_ib` (default `nothing`).
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_aero_link_angles`.
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

# _aero_link_angles

## Purpose
Computes a link's angle of attack, sideslip, and (when attitude is propagated) its body-to-inertial rotation matrix for the free-molecular coefficient evaluation.

## Design & Implementation
Signature `(spacecraft, body, root_index, orientation_sim, vel_pi, θ_body, fixed_attitude_incidence=:max_drag, q_root_ib=nothing)` returning `(α, β, R)`. With `orientation_sim`, it requires the propagated `q_root_ib` (throws `ArgumentError` if `nothing`), builds `R_root = rot(q_root_ib)'`, composes `R = R_root * rot(body.q)'` for children, rotates `body_frame_velocity = R' * vel_pi`, and returns `α = atan(v1, v3)`, `β = atan(v2, hypot(v1, v3))`, and `R`. Otherwise `α` is `_attitude_link_alpha` (`:attitude`), `π/2` (`:tumbling_average`), or the `:max_drag` root/child rule, with `β = 0` and `R = nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `body` | Any | n/a | yes | Positional argument `body`. |
| in | `root_index` | Int | n/a | yes | Positional argument `root_index`. |
| in | `orientation_sim` | Bool | n/a | yes | Positional argument `orientation_sim`. |
| in | `vel_pi` | SVector{3, Float64} | n/a | yes | Positional argument `vel_pi`. |
| in | `theta_body` | Float64 | n/a | yes | Positional argument `θ_body`. |
| in | `fixed_attitude_incidence` | Symbol | n/a | no | Positional argument `fixed_attitude_incidence` (default `:max_drag`). |
| in | `q_root_ib` | Union{Nothing, SVector{4, Float64}} | n/a | no | Positional argument `q_root_ib` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_aero_link_angles`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:419-419`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:303-303`
- `callees` → [[dynamics.aerodynamic_wrench_models__attitude_link_alpha|_attitude_link_alpha]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:312-312`
- `callees` → [[dynamics.aerodynamic_wrench_models__quaternion_link_alpha|_quaternion_link_alpha]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:316-316`
<!-- vulcan:connections:end -->

## Limitations
`root_index` and `θ_body` are accepted but unused. The comment documents a past bug where the static `Link.q` was used instead of the propagated attitude, freezing aero geometry; the fix relies on callers passing `x.q_ib`. Constructs `SVector{4}(body.q...)` per link per call, which splats a possibly non-static field.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 285.
