---
id: dynamics.cloth_multibody_build_rectangular_compliant_grid
label: build_rectangular_compliant_grid
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: build_rectangular_compliant_grid
  lines:
  - 342
  - 342
inputs:
- id: rows
  type: Integer
  units: n/a
  required: true
  description: Positional argument `rows`.
- id: cols
  type: Integer
  units: n/a
  required: true
  description: Positional argument `cols`.
- id: spacing_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `spacing_m` (default `(1.0, 1.0)`).
- id: tile_size_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `tile_size_m` (default `spacing_m`).
- id: thickness_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `thickness_m` (default `1.0e-3`).
- id: mass_kg
  type: Any
  units: n/a
  required: false
  description: Keyword argument `mass_kg` (default `1.0`).
- id: k_translation_n_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `k_translation_n_m` (default `5.0e3`).
- id: c_translation_n_s_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `c_translation_n_s_m` (default `30.0`).
- id: k_rotation_n_m_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `k_rotation_n_m_rad` (default `15.0`).
- id: c_rotation_n_m_s_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `c_rotation_n_m_s_rad` (default `0.5`).
- id: anchor_index
  type: Union{Nothing, Integer}
  units: n/a
  required: false
  description: Keyword argument `anchor_index` (default `nothing`).
- id: anchor_k_translation_n_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `anchor_k_translation_n_m` (default `k_translation_n_m`).
- id: anchor_c_translation_n_s_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `anchor_c_translation_n_s_m` (default `c_translation_n_s_m`).
- id: anchor_k_rotation_n_m_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `anchor_k_rotation_n_m_rad` (default `k_rotation_n_m_rad`).
- id: anchor_c_rotation_n_m_s_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `anchor_c_rotation_n_m_s_rad` (default `c_rotation_n_m_s_rad`).
- id: base_position
  type: Any
  units: n/a
  required: false
  description: Keyword argument `base_position` (default `(0.0, 0.0, 0.0)`).
- id: base_quaternion
  type: Any
  units: n/a
  required: false
  description: Keyword argument `base_quaternion` (default `_Q_IDENTITY`).
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
  type: CompliantTopologyBuild
  units: n/a
  description: Return value of `build_rectangular_compliant_grid`.
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

# build_rectangular_compliant_grid

## Purpose
Constructs a rows-by-columns grid of thin panels joined by spring-damper edges, optionally anchored to the base, as the standard cloth test fixture.

## Design & Implementation
Validates positive counts, spacing and tile size, then places tiles centred on the origin in the xy plane with `thin_panel_inertia` from the tile size and thickness and mass from `_grid_value`. Each tile gets a `right` edge to its column neighbour and a `down` edge to its row neighbour, with attachment points at the half-spacing on the facing sides and the shared stiffness and damping defaults of 5 kN/m, 30 N·s/m, 15 N·m/rad and 0.5 N·m·s/rad. If `anchor_index` is given, a `base_anchor` edge from body zero to that tile is added with its own stiffness set. Delegates to `build_compliant_topology`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rows` | Integer | n/a | yes | Positional argument `rows`. |
| in | `cols` | Integer | n/a | yes | Positional argument `cols`. |
| in | `spacing_m` | Any | n/a | no | Keyword argument `spacing_m` (default `(1.0, 1.0)`). |
| in | `tile_size_m` | Any | n/a | no | Keyword argument `tile_size_m` (default `spacing_m`). |
| in | `thickness_m` | Real | n/a | no | Keyword argument `thickness_m` (default `1.0e-3`). |
| in | `mass_kg` | Any | n/a | no | Keyword argument `mass_kg` (default `1.0`). |
| in | `k_translation_n_m` | Any | n/a | no | Keyword argument `k_translation_n_m` (default `5.0e3`). |
| in | `c_translation_n_s_m` | Any | n/a | no | Keyword argument `c_translation_n_s_m` (default `30.0`). |
| in | `k_rotation_n_m_rad` | Any | n/a | no | Keyword argument `k_rotation_n_m_rad` (default `15.0`). |
| in | `c_rotation_n_m_s_rad` | Any | n/a | no | Keyword argument `c_rotation_n_m_s_rad` (default `0.5`). |
| in | `anchor_index` | Union{Nothing, Integer} | n/a | no | Keyword argument `anchor_index` (default `nothing`). |
| in | `anchor_k_translation_n_m` | Any | n/a | no | Keyword argument `anchor_k_translation_n_m` (default `k_translation_n_m`). |
| in | `anchor_c_translation_n_s_m` | Any | n/a | no | Keyword argument `anchor_c_translation_n_s_m` (default `c_translation_n_s_m`). |
| in | `anchor_k_rotation_n_m_rad` | Any | n/a | no | Keyword argument `anchor_k_rotation_n_m_rad` (default `k_rotation_n_m_rad`). |
| in | `anchor_c_rotation_n_m_s_rad` | Any | n/a | no | Keyword argument `anchor_c_rotation_n_m_s_rad` (default `c_rotation_n_m_s_rad`). |
| in | `base_position` | Any | n/a | no | Keyword argument `base_position` (default `(0.0, 0.0, 0.0)`). |
| in | `base_quaternion` | Any | n/a | no | Keyword argument `base_quaternion` (default `_Q_IDENTITY`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantTopologyBuild | n/a | — | Return value of `build_rectangular_compliant_grid`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:374-374`
- `callees` → [[dynamics.cloth_multibody__grid_value|_grid_value]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:374-374`
- `callees` → [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:432-432`
- `callees` → [[dynamics.cloth_multibody_complianttopologyedge|CompliantTopologyEdge]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:388-388`
- `callees` → [[dynamics.cloth_multibody_complianttopologynode|CompliantTopologyNode]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:375-375`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:375-375`
<!-- vulcan:connections:end -->

## Limitations
Only the four-neighbour lattice is built, with no diagonal shear edges, so the grid has a soft in-plane shear mode governed entirely by the rotational springs.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 342.
