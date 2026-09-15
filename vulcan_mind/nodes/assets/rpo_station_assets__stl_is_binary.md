---
id: assets.rpo_station_assets__stl_is_binary
label: _stl_is_binary
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: _stl_is_binary
  lines:
  - 44
  - 44
inputs:
- id: path
  type: AbstractString
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
  type: Any
  units: n/a
  description: Return value of `_stl_is_binary`. Returns `stat(path).size == 84 +
    50 * Int(ntri)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- assets
charts:
- assets
origin: agent
---

# _stl_is_binary

## Purpose
Distinguishes binary from ASCII STL by the one structural property that is reliable: whether the declared triangle count matches the file size.

## Design & Implementation
Returns false for any file under 84 bytes, the minimum binary header plus count. It reads exactly 84 bytes, reinterprets bytes 81 to 84 as a little-endian `UInt32` triangle count, and returns whether the file size equals `84 + 50 * ntri`, since each binary triangle occupies fifty bytes. This avoids the classic failure of sniffing for the ASCII keyword `solid`, which many binary exporters also write into the header.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | AbstractString | n/a | yes | Positional argument `path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_stl_is_binary`. Returns `stat(path).size == 84 + 50 * Int(ntri)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[assets.rpo_station_assets__load_stl_triangles|_load_stl_triangles]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:55-55`
- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A binary STL with trailing bytes, or one whose count field is wrong, is misclassified as ASCII and then parsed as text, yielding zero vertices rather than an error; the header is read into a fresh 84-byte buffer on every call.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 44.
