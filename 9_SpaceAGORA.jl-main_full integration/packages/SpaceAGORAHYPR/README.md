# SpaceAGORAHYPR

Compatibility entry point for the separately maintained HYPR package. The shim
contains no search implementation. Use the SpaceAGORA checkout's supported setup
helper; it installs the exact HYPR Git revision from `HYPRSource.toml`:

```sh
julia --startup-file=no scripts/setup_hypr.jl /path/to/my-rpo-project
julia --project=/path/to/my-rpo-project -e 'using SpaceAGORAHYPR'
```

This checkout pins HYPR 0.1.1 at `39e4cb9648fa94b15712d00017f12c3aad451e69`,
paired with SpaceAGORA 0.2.0, HYPRServices 1.0.0 and this 0.2.0 shim.
`HYPRSource.toml` is the authoritative source pin; the package compatibility range
does not certify every HYPR revision in that range.
Package sources do not propagate from a dependency's Project.toml,
so developing the shim alone is not the supported installation recipe. Both
loading orders are supported. An incomplete or failed extension receives a
compatibility error before aliases are used; discard that process after failure.

The source pin records the exact HYPR revision used by this checkout. The
`SPACEAGORA_HYPR_PATH` override is only for explicit coordinated development.
