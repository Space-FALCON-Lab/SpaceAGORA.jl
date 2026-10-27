# Result and checkpoint files

| File | Module | What it does |
|---|---|---|
| `paths/io_paths.jl` | `IOConfig` | Derives result and checkpoint file names from the configuration. It does not parse simulation settings or decide when to write a checkpoint. |
| `checkpoints/checkpoint_storage.jl` | `IOSerialization` | Reads and writes checkpoint payloads and their manifests, with atomic file replacement and checksums. |
| `outputs/io_outputs.jl` | `IOOutputs` | Builds the result tables from the saved samples and writes them in their output formats. |

The simulation engine, not this folder, decides when to save a checkpoint and
whether to resume from one.
