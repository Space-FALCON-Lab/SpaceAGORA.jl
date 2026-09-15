---
id: io.io_serialization__sha256_hex
label: _sha256_hex
kind: function
source:
  file: src/io/serialization/io_serialization.jl
  symbol: _sha256_hex
  lines:
  - 28
  - 28
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
  type: String
  units: n/a
  description: Return value of `_sha256_hex`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
charts:
- io
origin: agent
---

# _sha256_hex

## Purpose
Produces the lowercase hexadecimal SHA-256 digest of a file, used to record checkpoint integrity in the manifest.

## Design & Implementation
Opens the path for reading inside a `do` block so the handle closes even on error, reads the entire contents into a byte vector, hashes it with `SHA.sha256` and converts the digest through `bytes2hex`. Reading whole rather than streaming keeps the implementation to one expression, which is acceptable because it is called once per checkpoint rather than in a loop.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | String | n/a | yes | Positional argument `path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_sha256_hex`. |
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
The whole file is held in memory at once, so a very large checkpoint payload doubles peak memory during the hash; there is no streaming or chunked variant.

## Provenance
Mapped from `src/io/serialization/io_serialization.jl` line 28.
