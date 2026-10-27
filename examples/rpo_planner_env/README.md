# RPO planner pilot

A separate environment for running rendezvous and proximity operations (RPO)
with interchangeable planners, including HYPR. From the repository root,
prepare it once:

```sh
julia --startup-file=no examples/rpo_planner_env/setup.jl
```

Setup uses HYPR's standalone installer and the same source-pair and destination
safeguards as EDG and MuJoCo. The managed environment lives under the first Julia
depot's `environments/` directory, with a name derived from this checkout. It is
outside both source checkouts. No Project or Manifest is written into this example
folder, and the tracked template remains unchanged. Set `SPACEAGORA_RPO_ENV` to
choose another separate destination. The setup script and launcher use the same
setting; repeat setup reuses only an environment owned by this source pair.

The launcher below selects that environment. After setup the smoke run needs no
native assets and no network:

```sh
JULIA_PKG_OFFLINE=true julia --startup-file=no examples/rpo_planner_env/run.jl examples/rpo_planner_env/smoke.jl
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
julia --startup-file=no examples/rpo_planner_env/run.jl examples/Earth_RPO_CubeSat_MPC.jl
```

Without the examples environment they stop at once and print its setup command.
They also need the assets the
[examples catalog](../../docs/src/user/examples_catalog.md) lists for them; the
smoke run above needs none. Examples that do not use HYPR run in the examples
environment alone.
