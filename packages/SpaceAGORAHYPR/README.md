# SpaceAGORAHYPR

Compatibility entry point for the separately maintained HYPR package. The shim
contains no search implementation. Use the SpaceAGORA checkout's supported setup
helper; it installs the exact HYPR Git revision from `HYPRSource.toml`:

```sh
julia --startup-file=no scripts/setup_hypr.jl /path/to/my-rpo-project
julia --project=/path/to/my-rpo-project -e 'using SpaceAGORAHYPR'
```

The proposed pair is HYPR 0.1.0, SpaceAGORA 0.2.0 with HYPRServices 1.0.0 and this
0.2.0 shim. Package sources do not propagate from a dependency's Project.toml,
so developing the shim alone is not the supported installation recipe. Both
loading orders are supported. An incomplete or failed extension receives a
compatibility error before aliases are used; discard that process after failure.

The candidates are local until publication is approved and verified. The
`SPACEAGORA_HYPR_PATH` override is only for explicit coordinated development.
