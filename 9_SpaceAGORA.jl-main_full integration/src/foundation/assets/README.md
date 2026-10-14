# Asset responsibilities

`asset_catalog.jl` reads the repository asset manifest, resolves its paths, and
checks local file or directory readiness: missing paths, empty directories,
Git LFS pointers, zero-byte payloads and payload byte presence. The service
returns field records and performs no downloads. The CLI retains its public report types, command
parsing, and presentation, while delegating these reusable operations.

`odyssey_surrogate_assets.jl` owns the Odyssey artifact lookup and verification flow.
`rpo_station_assets.jl` owns station asset lookup and its existing CAD sampling
behavior. Its binary STL decoding uses the geometry reader in
`SimulationModel.Structure`; ASCII parsing and sampling retain their separate
legacy behavior.

Atmosphere preset schema and scientific validation remain in the atmosphere
model. A generic asset readiness check does not establish that a model or
native runtime is scientifically qualified.
