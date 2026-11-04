# I/O responsibilities

`paths/io_paths.jl` derives result and checkpoint filenames from configuration.
It retains the `IOConfig` module identity. It does not parse simulation settings
or decide when a checkpoint should be written.

`checkpoints/checkpoint_storage.jl` reads and writes checkpoint payloads and their
manifests, including atomic file replacement and checksums. Its module remains
`IOSerialization`, so moving the source does not rename Julia bindings or change
the stored format. Checkpoint scheduling and restart decisions belong to the
simulation engine.

`outputs/` contains result writers and their output-format handling.
