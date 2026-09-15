---
id: parcore.compat_model_codes_legacydensitymodelcode
label: LegacyDensityModelCode
kind: struct
source:
  file: src/core/types/compat_model_codes.jl
  symbol: LegacyDensityModelCode
  lines:
  - 17
  - 23
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Base enum machinery; the module is standalone and carries no simulation
    dependencies.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: code
  type: LegacyDensityModelCode
  units: n/a
  description: Int8-backed enumeration selecting the historical atmosphere model identifier
    used by legacy configuration files and comparison scripts.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- parcore
origin: agent
---

# LegacyDensityModelCode

## Purpose
`LegacyDensityModelCode` is the enumeration of historical atmosphere-model identifiers kept for backward compatibility with the numeric model codes used by earlier versions of the tool and by archived input decks. It lets old configuration files and regression comparisons be read without reintroducing the integer literals into current code.

## Model & Assumptions
The enumeration is backed by `Int8` and pins the exact historical numbering: constant density is 0, exponential is 1, no-density is 2, GRAM is 3 and NRLMSISE is 4. Preserving the numeric values matters because archived files store the integer, not the name, so any renumbering would silently reinterpret existing data.

## Design & Implementation
The file declares the `LegacyModelCodes` module with five parallel `@enum` blocks, one per model family: gravity, density, aerodynamics, thermal and thrust control. Each block is `::Int8` and each member name is prefixed with `Legacy` so that flattening the module into the `SimulationModel` namespace cannot shadow a current model type. Two export lines list the enum types, and the remaining export lines list every member so that a caller can reference `LegacyDensityExponential` directly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Base enum machinery; the module is standalone and carries no simulation dependencies. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `code` | LegacyDensityModelCode | n/a | — | Int8-backed enumeration selecting the historical atmosphere model identifier used by legacy configuration files and comparison scripts. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/compat_model_codes.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The codes describe intent only; they carry no parameters, so a legacy exponential density code says nothing about the scale height that accompanied it in the original deck. The thermal family starts at one rather than zero, so code that maps enum values to array indices must handle the two numberings separately. Nothing in this module converts a legacy code into a current model instance; that mapping lives with the reader that consumes the archived configuration.

## Provenance
Mapped from `src/core/types/compat_model_codes.jl:17-23`.
