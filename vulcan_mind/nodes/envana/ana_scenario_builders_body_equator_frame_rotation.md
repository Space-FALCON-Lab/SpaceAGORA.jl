---
id: envana.ana_scenario_builders_body_equator_frame_rotation
label: _body_equator_frame_rotation
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _body_equator_frame_rotation
  lines:
  - 484
  - 492
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace supplying StaticArrays and the linear
    algebra norm and cross product.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: rotation
  type: SMatrix{3,3,Float64,9}
  units: n/a
  description: Orthonormal rotation matrix from the body equatorial frame to J2000.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- envana
origin: agent
---
# _body_equator_frame_rotation

## Purpose
`_body_equator_frame_rotation` constructs the rotation matrix of a planet's equatorial frame from its J2000 pole direction, so that scenario initial conditions expressed in body-equatorial elements can be converted into inertial J2000 states.

## Theory & Math
Given the pole vector `p` in J2000, the frame is built as `zhat = p / |p|`, `n = (-p_y, p_x, 0)`, `xhat = n / |n|`, and `yhat = zhat x xhat`. The vector `n` is the cross product of the J2000 z-axis with the pole, so it lies in the J2000 equatorial plane and points along the ascending node of the body equator on that plane; this is the standard node-line construction. The returned matrix `[xhat yhat zhat]` is orthonormal by construction with determinant `+1`, so its transpose is its inverse and no explicit matrix inversion is ever needed. All three columns are dimensionless unit vectors.

## Model & Assumptions
The construction degenerates when the pole is parallel to the J2000 z-axis, because then `n` vanishes and the node line is undefined. The code detects this with `node_mag <= 1e-12` and returns the identity matrix, which is the natural choice since a pole already aligned with z means the body equatorial frame and the J2000 frame share an equator and only a prime-meridian convention distinguishes them. The input pole need not be normalised, since the first operation divides by its norm.

## Design & Implementation
The function is nine lines and works entirely with `SVector{3,Float64}`, so it allocates nothing on the heap and is safe to call inside scenario construction loops. Unicode hat identifiers `ẑ`, `x̂`, and `ŷ` mirror the mathematical notation directly. The columns are assembled with `hcat`, and the declared return type `SMatrix{3,3,Float64,9}` pins the static size at the type level so downstream matrix products are fully specialised.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace supplying StaticArrays and the linear algebra norm and cross product. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `rotation` | SMatrix{3,3,Float64,9} | n/a | — | Orthonormal rotation matrix from the body equatorial frame to J2000. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__initial_condition_in_j2000|_initial_condition_in_j2000]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:513-513`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `1e-12` degeneracy threshold is a fixed absolute tolerance on a normalised quantity and is not configurable, so a pole within roughly a picoradian of the z-axis silently falls back to identity. The frame captures only the pole orientation; the prime meridian is not represented, so this matrix alone cannot place a body-fixed longitude. No check confirms that the supplied pole is non-zero, and a zero vector yields a matrix of NaNs.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/scenario_builders.jl:484-492`, including the node-line construction and the degeneracy guard.
