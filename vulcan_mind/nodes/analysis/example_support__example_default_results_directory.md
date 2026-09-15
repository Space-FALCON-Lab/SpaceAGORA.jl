---
id: analysis.example_support__example_default_results_directory
label: _example_default_results_directory
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: _example_default_results_directory
  lines:
  - 5
  - 5
inputs:
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
  type: String
  units: n/a
  description: Return value of `_example_default_results_directory`.
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

# _example_default_results_directory

## Purpose
Resolves where example runs write their outputs, honouring a CLI-provided override and otherwise defaulting into the repository.

## Design & Implementation
Reads `SPACEAGORA_CLI_OUTPUT_DIR`, strips it, and if non-empty returns its `abspath` so relative overrides are anchored to the current working directory at call time. When unset it returns `joinpath(REPO_ROOT, "output")`. Declared `@inline` with a `::String` return, and used as the default value of `make_example_config`'s `results_directory` keyword so every example shares the same resolution rule.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_example_default_results_directory`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:132-132`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The directory is resolved but not created, so a caller pointing the override at a non-existent path gets a failure at first write rather than here; `abspath` is evaluated at call time, so two calls from different working directories with the same relative override yield different locations.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/example_support.jl` line 5.
