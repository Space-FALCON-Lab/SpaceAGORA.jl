---
id: simx.runtime_services_runtimeservices
label: RuntimeServices
kind: struct
source:
  file: src/simulation/runtime_services.jl
  symbol: RuntimeServices
  lines:
  - 1
  - 17
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
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: shared_native_lock
  type: ReentrantLock
  units: n/a
  description: The single reentrant lock exported as both SPICE_LOCK and GRAM_LOCK,
    guarding every native CSPICE-backed call in the process.
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
# RuntimeServices

## Purpose
`RuntimeServices` exists to own one object: the reentrant lock that serialises every native call which ultimately reaches CSPICE. It defines `SPICE_LOCK` and then binds `GRAM_LOCK` to the very same lock, so the two names are aliases rather than two locks of the same shape.

## Model & Assumptions
The module's comment records the concurrency defect that forced this design. `libGRAM.dylib` statically links its own copy of CSPICE, but the copies are not isolated at the operating system's symbol-resolution level: `nm -gU libGRAM.dylib` shows it globally exporting CSPICE internals such as `chkin_`, `chkout_`, `trcpkg_` and `subslr_c` under exactly the names SPICE.jl's libcspice also exports. Running a native GRAM atmosphere call on one thread while a SpaceAGORA ephemerides or frame-transform call runs on another therefore corrupts CSPICE's internal call-trace stack, surfacing as `SPICE(NAMESDONOTMATCH)` and CHKOUT errors.

## Design & Implementation
When the two locks were distinct, each guarded only half of what is actually a single C-level critical section, so the failure was intermittent and thread-count dependent. Aliasing them makes that critical section indivisible by construction. The lock is a `ReentrantLock` because call chains re-enter it: `_build_nbody_ephemeris_cache` takes the lock around its whole sample loop and then calls the deliberately unlocked `_spice_position_j2000_m_unlocked` variant inside, and other paths acquire it at outer and inner levels. Keeping the lock in a tiny module with no other content means every consumer can `import ..RuntimeServices` without dragging in the engine.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `shared_native_lock` | ReentrantLock | n/a | — | The single reentrant lock exported as both SPICE_LOCK and GRAM_LOCK, guarding every native CSPICE-backed call in the process. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/runtime_services.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A single process-wide lock makes every SPICE and GRAM query a serial section, which is the dominant scaling limit for constellations with per-satellite atmosphere or third-body queries; the ephemeris caches exist precisely to amortise it. The alias is a source-level convention with no enforcement, so a future edit that gives `GRAM_LOCK` its own `ReentrantLock()` would reintroduce the corruption silently. Code that calls an unlocked variant without holding the lock is equally unguarded.

## Provenance
Mapped from `src/simulation/runtime_services.jl:1-17`; the lock is consumed by the ephemeris cache builder at `src/simulation/engine/setup.jl:1560` and imported by the engine at `src/simulation/engine/simulation_engine.jl:5`.
