# Public RPO planner pilot

From the repository root, prepare this separate environment once:

```sh
julia --project=examples/rpo_planner_env examples/rpo_planner_env/setup.jl
```

Setup creates a local `Project.toml` from the tracked `Project.template.toml`
and seeds `Manifest.toml` from the repository dependency manifest. It adds this
checkout and the `SpaceAGORAHYPR` compatibility shim by path, plus HYPR at the
immutable revision in `packages/SpaceAGORAHYPR/HYPRSource.toml`, preserving the
existing external dependency versions. Both generated installation files are
ignored; the tracked template stays unchanged. Repeating setup reuses the local
environment. After setup the run needs neither native assets nor network access:

```sh
JULIA_PKG_OFFLINE=true julia --project=examples/rpo_planner_env examples/rpo_planner_env/smoke.jl
```

The HYPR RPO/Cloth examples use this same environment, for example:

```sh
julia --project=examples/rpo_planner_env examples/Earth_RPO_CubeSat_MPC.jl
```

Those research examples require the assets described in the examples catalog and
load the companion before activating the root environment for their other imports.
Baseline examples still use the root project without loading the companion.

`smoke.jl` explicitly loads `SpaceAGORAHYPR`, then uses the public SpaceAGORA planner surface and its declared dependencies.
It runs the same two-second Earth corridor case with `DirectRPOPlanner` and
`HYPRRPOPlanner`, then defines a small user-owned planner and force. It does not
switch projects, include package source, or access internal modules. Test imports
are verification tools, not part of the planner API.

The short run checks integration and finite, bounded controller outputs. It does
not finish the reference, certify collision avoidance or accept physical tracking
for arbitrary missions. See the [planner guide](../../docs/src/user/rpo_planner_pilot.md).
