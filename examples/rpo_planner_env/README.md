# RPO planner pilot

A separate environment for running rendezvous and proximity operations (RPO)
with interchangeable planners, including HYPR. From the repository root,
prepare it once:

```sh
julia --project=examples/rpo_planner_env examples/rpo_planner_env/setup.jl
```

Setup copies the tracked `Project.template.toml` to a local `Project.toml`,
starts `Manifest.toml` from the repository manifest, adds this checkout by
path, and adds HYPR at the revision pinned in `packages/HYPRSource.toml`.
Packages the repository manifest records keep their versions; only the
libraries bundled with Julia, such as OpenSSL_jll, follow your Julia install.
Git ignores both generated files, and running setup again reuses them. After
setup the pilot needs no native assets and no network:

```sh
JULIA_PKG_OFFLINE=true julia --project=examples/rpo_planner_env examples/rpo_planner_env/smoke.jl
```

`smoke.jl` loads `SpaceAGORA` and `HYPR` and uses only the public planner API.
It runs the same two-second Earth case with `DirectRPOPlanner` and with
`HYPRRPOPlanner`, then defines a small planner and a force of its own. It does
not switch projects, include package source or reach into internal modules;
the `Test` imports are only there to check the results. The run checks that
the pieces fit together and that the controller outputs are finite and
bounded. It does not finish the reference, prove collision avoidance or show
that the tracking is good enough for any other mission. The
[planner guide](../../docs/src/user/rpo_planner_pilot.md) explains the setup.

The HYPR examples also start in this environment: `Earth_RPO_CubeSat_MPC.jl`
and its batch, planner-comparison and replanning variants, and
`Robot_Arm_Planner_Cloth_Demo.jl`, all under `examples/`. They load HYPR here,
then continue in the examples environment for their other packages, so prepare
that one too:

```sh
julia --startup-file=no examples/setup.jl
julia --project=examples/rpo_planner_env examples/Earth_RPO_CubeSat_MPC.jl
```

Without the examples environment they stop at once and print its setup command.
They also need the assets the
[examples catalog](../../docs/src/user/examples_catalog.md) lists for them; the
smoke run above needs none. Examples that do not use HYPR run in the examples
environment alone.
