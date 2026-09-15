---
id: simx.engine_setup_build_nbody_ephemeris_cache
label: _build_nbody_ephemeris_cache
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _build_nbody_ephemeris_cache
  lines:
  - 1547
  - 1573
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
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: SPICE name of the central body all third-body positions are expressed
    relative to.
- id: body_query_names
  type: Vector{String}
  units: n/a
  required: true
  description: SPICE query names of the perturbing bodies, in the column order the
    cache will preserve.
- id: mission_window
  type: NTuple{2,Float64}
  units: s
  required: true
  description: Ephemeris start time et_start and mission duration mission_end_s that
    bound the sampled interval.
- id: dt_s
  type: Float64
  units: s
  required: true
  description: Uniform sample spacing of the tabulated ephemeris.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: nbody_cache
  type: SimulationModel.NBodyEphemerisCache
  units: m
  description: Tabulated third-body positions in J2000 metres with their sample times,
    body query names and a name-to-column index map.
- id: spice_call_counter
  type: Base.Threads.Atomic{Int}
  units: count
  description: Optional counter incremented once per SPICE position query for prewarm
    progress reporting.
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
# _build_nbody_ephemeris_cache

## Purpose
`_build_nbody_ephemeris_cache` tabulates third-body positions across a mission window so the right-hand side can interpolate instead of calling SPICE. It is the construction half of the n-body ephemeris cache whose reuse and prewarm layers sit around it in the same file.

## Theory & Math
Sample times are laid out on a clamped uniform grid,

$$t_i = t_0 + \min\!\big((i-1)\,\Delta t,\ T_{end}\big), \qquad i = 1,\dots,N$$

so the final sample lands exactly on the mission end even when $T_{end}$ is not an integer multiple of $\Delta t$, and no query is ever made past the window. The sample count $N$ comes from `_nbody_ephemeris_cache_sample_count`. The storage is a dense $N \times M$ matrix of `SVector{3,Float64}` for $M$ bodies, so the cost model is $N M$ SPICE calls up front against one interpolation per body per RHS evaluation afterwards. Accuracy is bounded by the interpolation error over $\Delta t$, which for the slow planetary motions involved falls roughly as $\mathcal{O}(\Delta t^{2})$ for the linear scheme.

## Model & Assumptions
Positions are queried in J2000 and stored in metres relative to `primary_body_name`, matching the frame the gravity effectors work in, so no transform is needed at use time. The grid is uniform, which makes lookup an index computation rather than a search.

## Design & Implementation
The entire sample loop runs inside `lock(RuntimeServices.SPICE_LOCK) do ... end`, taking the shared lock once rather than once per query: this is both a correctness requirement, because CSPICE state is global, and the main reason prewarming is faster than lazy lookup. Inside the lock it calls `SimulationModel.EphemeridesModels._spice_position_j2000_m_unlocked`, the deliberately unlocked variant, since the lock is already held. `_increment_atomic_counter!` bumps the optional progress counter with `Base.Threads.atomic_add!` and returns immediately when no counter was supplied. The finished arrays are handed to `_nbody_ephemeris_cache_from_samples`, which builds the name-to-index map and takes a `copy` of the body-name vector so the cache does not alias the caller's array.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `primary_body_name` | String | n/a | yes | SPICE name of the central body all third-body positions are expressed relative to. |
| in | `body_query_names` | Vector{String} | n/a | yes | SPICE query names of the perturbing bodies, in the column order the cache will preserve. |
| in | `mission_window` | NTuple{2,Float64} | s | yes | Ephemeris start time et_start and mission duration mission_end_s that bound the sampled interval. |
| in | `dt_s` | Float64 | s | yes | Uniform sample spacing of the tabulated ephemeris. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `nbody_cache` | SimulationModel.NBodyEphemerisCache | m | — | Tabulated third-body positions in J2000 metres with their sample times, body query names and a name-to-column index map. |
| out | `spice_call_counter` | Base.Threads.Atomic{Int} | count | — | Optional counter incremented once per SPICE position query for prewarm progress reporting. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1834-1834`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1731-1731`

**Downstream**

- `callees` → [[environment.simple_ephemerides__spice_position_j2000_m_unlocked|_spice_position_j2000_m_unlocked]] · `callers` · call · `src/simulation/engine/setup.jl:1566-1566`
- `callees` → [[simulation.setup__increment_atomic_counter_bang|_increment_atomic_counter!]] · `callers` · call · `src/simulation/engine/setup.jl:1567-1567`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_from_samples|_nbody_ephemeris_cache_from_samples]] · `callers` · call · `src/simulation/engine/setup.jl:1572-1572`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_sample_count|_nbody_ephemeris_cache_sample_count]] · `callers` · call · `src/simulation/engine/setup.jl:1555-1555`
<!-- vulcan:connections:end -->

## Limitations
Memory grows as the product of sample count and body count, so a long mission with a small `dt_s` and several perturbing bodies can dominate the run's footprint. Because the whole build holds the SPICE lock, no other thread can perform an ephemeris or frame query while a cache is being constructed. The clamped grid repeats the final epoch when the window is not an exact multiple of the step, which any consumer differentiating the table must tolerate.

## Provenance
Mapped from `src/simulation/engine/setup.jl:1547-1573`, with the cache constructor at line 1525, the prewarm driver at line 1681 and the module-level reuse caches at lines 25-28 of the same file.
