---
id: simx.engine_state_access_state_position_ii
label: _state_position_ii
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _state_position_ii
  lines:
  - 32
  - 44
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: state_u
  type: AbstractVector
  units: m,m/s
  required: true
  description: Solver state in either the standard ComponentVector layout or one of
    the gravity-backbone partitioned layouts.
- id: sat_index
  type: Int
  units: index
  required: true
  description: One-based spacecraft index selecting which satellite block to read.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: position_ii
  type: SVector{3,Float64}
  units: m
  description: Planet-centred inertial position of the selected spacecraft, copied
    out of the state block.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# _state_position_ii

## Purpose
`_state_position_ii` reads one spacecraft's inertial position out of the solver state without the caller needing to know which state layout is in play. It is the accessor every effector and callback uses when it wants a position rather than a raw view.

## Model & Assumptions
Two layouts coexist. The standard path stores a `ComponentVector` whose `u.sc[i]` block begins with the three position components. The gravity-backbone solver modes instead hand the integrator a partitioned second-order state, where position and velocity live in separate halves and the spacecraft block may or may not still be wrapped in a named `sc` component. All three shapes place position in the first three slots of a spacecraft block, which is what makes a single accessor possible.

## Design & Implementation
The function branches on `_is_gravity_backbone_state(u)`. In the backbone case it extracts the position half with `_gravity_backbone_spacecraft_state(u)` and then tests `hasproperty(spacecraft_state, :sc)` to decide whether the satellites are reached through a named component or by direct indexing, covering both partitioned shapes. Every path constructs an explicit `SVector{3,Float64}` from elements one through three, so callers receive an immutable stack-allocated value rather than a view into the integrator's working memory, which prevents accidental writes and keeps the result usable inside `@batch` workers. The whole file follows this pattern with `@inline` accessors for velocity, mass, heat loads and quaternion, each paired with a `_state_has_*` predicate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `state_u` | AbstractVector | m,m/s | yes | Solver state in either the standard ComponentVector layout or one of the gravity-backbone partitioned layouts. |
| in | `sat_index` | Int | index | yes | One-based spacecraft index selecting which satellite block to read. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `position_ii` | SVector{3,Float64} | m | — | Planet-centred inertial position of the selected spacecraft, copied out of the state block. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.setup__initialize_in_atmosphere_flags_bang|_initialize_in_atmosphere_flags!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1350-1350`

**Downstream**

- `callees` → [[simulation.state_access__gravity_backbone_spacecraft_state|_gravity_backbone_spacecraft_state]] · `callers` · call · `src/simulation/engine/state_access.jl:34-34`
- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:33-33`
<!-- vulcan:connections:end -->

## Limitations
Layout detection is structural rather than typed, so a future state layout that is neither shape would be silently read as the standard one and return the wrong three numbers. `hasproperty` is resolved at compile time for concrete types but degrades if the state type is not inferable at the call site. No bounds check is performed on `sat_idx` beyond what the underlying container does.

## Provenance
Mapped from `src/simulation/engine/state_access.jl:32-44`, with the layout predicate at line 5 and the sibling velocity accessor at line 46 of the same file.
