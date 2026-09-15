---
id: analysis.types_spacecraftconfig
label: SpacecraftConfig
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: SpacecraftConfig
  lines:
  - 61
  - 61
inputs:
- id: bus_dims
  type: NTuple{3, Float64}
  units: n/a
  required: true
  description: Field `bus_dims`.
- id: panel_dims
  type: NTuple{3, Float64}
  units: n/a
  required: true
  description: Field `panel_dims`.
- id: bus_mass_kg
  type: Float64
  units: n/a
  required: true
  description: Field `bus_mass_kg`.
- id: panel_mass_each_kg
  type: Float64
  units: n/a
  required: true
  description: Field `panel_mass_each_kg`.
- id: panel_offset_y_m
  type: Float64
  units: n/a
  required: true
  description: Field `panel_offset_y_m`.
- id: prop_mass_kg
  type: Float64
  units: n/a
  required: true
  description: Field `prop_mass_kg`.
- id: id
  type: Int64
  units: n/a
  required: true
  description: Field `id`.
- id: bus_ram_face
  type: Symbol
  units: n/a
  required: false
  description: Field `bus_ram_face` (default `:legacy`).
- id: bus_attitude_q
  type: Union{Nothing, NTuple{4, Float64}}
  units: n/a
  required: false
  description: Field `bus_attitude_q` (default `nothing`).
- id: panel_attitude_q_left
  type: Union{Nothing, NTuple{4, Float64}}
  units: n/a
  required: false
  description: Field `panel_attitude_q_left` (default `nothing`).
- id: panel_attitude_q_right
  type: Union{Nothing, NTuple{4, Float64}}
  units: n/a
  required: false
  description: Field `panel_attitude_q_right` (default `nothing`).
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
  type: SpacecraftConfig
  units: n/a
  description: Constructed `SpacecraftConfig` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# SpacecraftConfig

## Purpose
Geometric and mass description of the bus-plus-two-panel spacecraft used by telemetry verification scenarios, feeding `make_three_body_spacecraft` and the free-molecular aerodynamic coefficient model.

## Design & Implementation
`Base.@kwdef struct` with required `bus_dims::NTuple{3,Float64}` and `panel_dims::NTuple{3,Float64}` (metres), `bus_mass_kg`, `panel_mass_each_kg`, `panel_offset_y_m`, `prop_mass_kg`, and an integer `id`. `bus_ram_face::Symbol=:legacy` selects the reference area: `:legacy` uses `dims[1]*dims[3]`, `:frontal` uses the flow-normal `dims[2]*dims[3]` face matching the Hart coefficient normalisation. Optional attitude quaternions `bus_attitude_q`, `panel_attitude_q_left`, `panel_attitude_q_right::Union{Nothing,NTuple{4,Float64}}` are scalar-last `(x,y,z,w)`, normalised at parse time; the bus quaternion is in the flow-aligned frame and panel quaternions are relative to the bus frame. `nothing` leaves a link at identity.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `bus_dims` | NTuple{3, Float64} | n/a | yes | Field `bus_dims`. |
| in | `panel_dims` | NTuple{3, Float64} | n/a | yes | Field `panel_dims`. |
| in | `bus_mass_kg` | Float64 | n/a | yes | Field `bus_mass_kg`. |
| in | `panel_mass_each_kg` | Float64 | n/a | yes | Field `panel_mass_each_kg`. |
| in | `panel_offset_y_m` | Float64 | n/a | yes | Field `panel_offset_y_m`. |
| in | `prop_mass_kg` | Float64 | n/a | yes | Field `prop_mass_kg`. |
| in | `id` | Int64 | n/a | yes | Field `id`. |
| in | `bus_ram_face` | Symbol | n/a | no | Field `bus_ram_face` (default `:legacy`). |
| in | `bus_attitude_q` | Union{Nothing, NTuple{4, Float64}} | n/a | no | Field `bus_attitude_q` (default `nothing`). |
| in | `panel_attitude_q_left` | Union{Nothing, NTuple{4, Float64}} | n/a | no | Field `panel_attitude_q_left` (default `nothing`). |
| in | `panel_attitude_q_right` | Union{Nothing, NTuple{4, Float64}} | n/a | no | Field `panel_attitude_q_right` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpacecraftConfig | n/a | — | Constructed `SpacecraftConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__parse_spacecraft_config|_parse_spacecraft_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:473-473`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Quaternion normalisation happens in the parser, not the constructor, so programmatically built configs can carry unnormalised quaternions. The struct does not enforce the manifest-layer rule that any non-`nothing` quaternion requires `aero_fixed_attitude_incidence = :attitude`; violating it silently changes default-mode physics. Only a symmetric two-panel layout is representable. `id` is `Int64` with no uniqueness guarantee across scenarios.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 61.
