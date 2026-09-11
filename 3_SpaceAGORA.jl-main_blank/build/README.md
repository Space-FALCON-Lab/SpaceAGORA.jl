# Build folder

This folder is used to create a custom Julia system image for the project.

## What it does

- Builds a precompiled Julia sysimage with the project packages baked in.
- Speeds up subsequent Julia startup times by avoiding repeated package compilation.
- Is defined by [build_sysimage.jl](build_sysimage.jl).

## Main script

- [build_sysimage.jl](build_sysimage.jl)

This script creates a file such as `SpaceAGORA.so` in the repository root.

## Is it required?

No. It is optional.

- Without it, the project can still run normally.
- It will just take longer to compile packages on first run.
- The build folder is mainly a performance optimization, not a required part of the application itself.
