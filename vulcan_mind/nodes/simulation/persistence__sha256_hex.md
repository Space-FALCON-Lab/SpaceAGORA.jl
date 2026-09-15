---
id: simulation.persistence__sha256_hex
label: _sha256_hex
kind: function
source:
  file: src/simulation/engine/persistence.jl
  symbol: _sha256_hex
  lines:
  - 8
  - 8
inputs:
- id: path
  type: String
  units: n/a
  required: true
  description: Positional argument `path`.
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
  type: SimulationModel.IOSerialization._sha256_hex
  units: n/a
  description: Return value of `_sha256_hex`. Returns `SimulationModel.IOSerialization._sha256_hex(path)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _sha256_hex

## Purpose
Computes the lowercase hexadecimal SHA-256 digest of a file on disk, used to stamp the data and CSV members of a results bundle into its manifest for integrity checking.

## Design & Implementation
An `@inline` forwarder to `SimulationModel.IOSerialization._sha256_hex(path)`. The implementation opens `path` in read mode, calls `read(io)` to slurp the entire contents into a byte vector, hashes it with `SHA.sha256`, and converts the 32-byte digest with `bytes2hex`. The `do` block guarantees the handle is closed even if hashing throws.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | String | n/a | yes | Positional argument `path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.IOSerialization._sha256_hex | n/a | — | Return value of `_sha256_hex`. Returns `SimulationModel.IOSerialization._sha256_hex(path)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.io_outputs_write_results_bundle__write_results_bundle_bang|_write_results_bundle!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:137-137`
- [[misc.io_serialization_write_checkpoint__write_checkpoint_bang|_write_checkpoint!]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:54-54`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the whole file is materialised in memory before hashing, a large Arrow results file costs its own size in transient allocation; a streaming, chunked digest would avoid that. A missing or unreadable path throws `SystemError` rather than returning a sentinel, so manifest writing must be ordered after the data file is fully flushed. The digest is taken at one instant and carries no locking, so a file still being written yields a hash of a partial state.

## Provenance
Mapped from `src/simulation/engine/persistence.jl` line 8.
