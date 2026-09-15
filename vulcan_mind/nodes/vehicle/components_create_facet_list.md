---
id: vehicle.components_create_facet_list
label: create_facet_list
kind: function
source:
  file: src/vehicle/spacecraft/components.jl
  symbol: create_facet_list
  lines:
  - 74
  - 74
inputs:
- id: area_list
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `area_list`.
- id: attitude_list
  type: Vector{SVector{4, Float64}}
  units: n/a
  required: true
  description: Positional argument `attitude_list`.
- id: normal_vector_list
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Positional argument `normal_vector_list`.
- id: cp_loc_list
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Positional argument `cp_loc_list`.
- id: diffuse_coeffs_list
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `diffuse_coeffs_list`.
- id: specular_coeffs_list
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `specular_coeffs_list`.
- id: facet_names_list
  type: Vector{String}
  units: n/a
  required: true
  description: Positional argument `facet_names_list`.
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
  description: Return value of `create_facet_list`. Returns `facet_vector`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# create_facet_list

## Purpose
Builds a vector of `Facet` objects from parallel per-property arrays, which is how facet geometry arrives from configuration files and mesh exports.

## Design & Implementation
Takes seven vectors — areas, attitudes, normals, centre-of-pressure locations, diffuse coefficients, specular coefficients and names — and asserts through `@assert` that all of them match the length of `area_list`, raising with a message naming the mismatch. It then preallocates a `Vector{Facet}` of that length and fills it positionally, so index `i` of every input describes the same patch.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `area_list` | Vector{Float64} | n/a | yes | Positional argument `area_list`. |
| in | `attitude_list` | Vector{SVector{4, Float64}} | n/a | yes | Positional argument `attitude_list`. |
| in | `normal_vector_list` | Vector{SVector{3, Float64}} | n/a | yes | Positional argument `normal_vector_list`. |
| in | `cp_loc_list` | Vector{SVector{3, Float64}} | n/a | yes | Positional argument `cp_loc_list`. |
| in | `diffuse_coeffs_list` | Vector{Float64} | n/a | yes | Positional argument `diffuse_coeffs_list`. |
| in | `specular_coeffs_list` | Vector{Float64} | n/a | yes | Positional argument `specular_coeffs_list`. |
| in | `facet_names_list` | Vector{String} | n/a | yes | Positional argument `facet_names_list`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `create_facet_list`. Returns `facet_vector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/components.jl`

**Downstream**

- `callees` → [[vehicle.components_facet|Facet]] · `callers` · call · `src/vehicle/spacecraft/components.jl:82-82`
<!-- vulcan:connections:end -->

## Limitations
The length assertion covers six of the seven inputs: `facet_names_list` is not included in the checked collection, so a short name vector is only caught by the bounds error at construction. `@assert` may be elided under `--check-bounds=no`, in which case a mismatched input reaches the loop unchecked.

## Provenance
Mapped from `src/vehicle/spacecraft/components.jl` line 74.
