---
id: core.simulation_configuration_filepaths
label: FilePaths
kind: struct
source:
  file: src/core/state/simulation_configuration.jl
  symbol: FilePaths
  lines:
  - 116
  - 116
inputs:
- id: results
  type: String
  units: n/a
  required: false
  description: Field `results` (default `"Results"`).
- id: GRAM
  type: String
  units: n/a
  required: false
  description: Field `GRAM` (default `"data/GRAMSuite.jl/GRAM Suite 2.0"`).
- id: SPICE
  type: String
  units: n/a
  required: false
  description: Field `SPICE` (default `"data/GRAMSuite.jl/GRAM Suite 2.0/SPICE"`).
- id: topography_harmonics
  type: String
  units: n/a
  required: false
  description: Field `topography_harmonics` (default `"data/Topography_harmonics_data"`).
- id: gravity_harmonics
  type: String
  units: n/a
  required: false
  description: Field `gravity_harmonics` (default `"data/Gravity_harmonics_data"`).
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
  type: FilePaths
  units: n/a
  description: Constructed `FilePaths` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# FilePaths

## Purpose
`FilePaths` holds the on-disk locations of external data assets (GRAM atmosphere suite, SPICE kernels, topography and gravity harmonics) and the default results directory. It is the first field of `SimulationConfiguration` and is consulted by the environment and I/O layers.

## Design & Implementation
A `@kwdef struct` with five `String` fields: `results = "Results"`, `GRAM = "data/GRAMSuite.jl/GRAM Suite 2.0"`, `SPICE = "data/GRAMSuite.jl/GRAM Suite 2.0/SPICE"`, `topography_harmonics = "data/Topography_harmonics_data"` and `gravity_harmonics = "data/Gravity_harmonics_data"`. The defaults are relative paths resolved against the working directory at run time. Source comments flag the two harmonics directories as candidates for moving onto the planet model.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `results` | String | n/a | no | Field `results` (default `"Results"`). |
| in | `GRAM` | String | n/a | no | Field `GRAM` (default `"data/GRAMSuite.jl/GRAM Suite 2.0"`). |
| in | `SPICE` | String | n/a | no | Field `SPICE` (default `"data/GRAMSuite.jl/GRAM Suite 2.0/SPICE"`). |
| in | `topography_harmonics` | String | n/a | no | Field `topography_harmonics` (default `"data/Topography_harmonics_data"`). |
| in | `gravity_harmonics` | String | n/a | no | Field `gravity_harmonics` (default `"data/Gravity_harmonics_data"`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | FilePaths | n/a | — | Constructed `FilePaths` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:236-236`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Relative defaults make behaviour depend on `pwd()`, which the telemetry runner works around by `cd`-ing into a temp directory. Existence of the directories is not checked at construction; missing assets fail later when GRAM or SPICE is initialised. The `SPICE` default is nested under the GRAM path, coupling two independent datasets. `results` overlaps in meaning with `SimulationSettings.results_directory` and the engine must decide which wins.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 116.
