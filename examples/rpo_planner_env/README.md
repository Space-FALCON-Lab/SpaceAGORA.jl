# Public RPO planner pilot

From the repository root, prepare this separate environment once:

```sh
julia --project=examples/rpo_planner_env examples/rpo_planner_env/setup.jl
```

It seeds a fresh environment from the repository dependency manifest, adds this
checkout and its optional `SpaceAGORAHYPR` companion by path and preserves those versions. The generated local manifest is
ignored. After setup the run needs neither native assets nor network access:

```sh
JULIA_PKG_OFFLINE=true julia --project=examples/rpo_planner_env examples/rpo_planner_env/smoke.jl
```

`smoke.jl` explicitly loads `SpaceAGORAHYPR`, then uses the public SpaceAGORA planner surface and its declared dependencies.
It runs the same two-second Earth corridor case with `DirectRPOPlanner` and
`HYPRRPOPlanner`, then defines a small user-owned planner and force. It does not
switch projects, include package source, or access internal modules. Test imports
are verification tools, not part of the planner API.

The short run checks integration and finite, bounded controller outputs. It does
not finish the reference, certify collision avoidance or accept physical tracking
for arbitrary missions. See the [planner guide](../../docs/src/user/rpo_planner_pilot.md).
