# Robotics and motion-planning packages

Two optional companions are pinned here as submodules. Their code and releases
belong to their own repositories. The folder name is only a grouping: HYPR
plans spacecraft rendezvous and proximity operations as well as robot-arm
motion.

| Package | Pinned revision | Role |
| --- | --- | --- |
| [HYPR.jl](https://github.com/Space-FALCON-Lab/HYPR.jl/tree/2816df11c86206cc348ec2c4f82732f958f18ad2) | `2816df11c86206cc348ec2c4f82732f958f18ad2` | Search and trajectory planning. Its SpaceAGORA extension uses the `HYPRServices` contract 3.0.0. |
| [SpaceAGORAMuJoCo.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/tree/7b5a00318f6883d7cdc7c7e7c0c18830f96ce57f) | `7b5a00318f6883d7cdc7c7e7c0c18830f96ce57f` | Orbital proximity scenes in MuJoCo, and spacecraft it propagates itself through the `ExternalPropagationServices` contract 2.0.0. |

## HYPR

From the SpaceAGORA repository root:

```sh
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
julia --startup-file=no --project=/path/to/hypr-env
```

Then load `SpaceAGORA` and `HYPR`. Loading both activates the adapter that HYPR
implements; there is no separate loader package. The
[HYPR user guide](../../docs/src/user/hypr_planning.md) explains the planner,
and the [RPO guide](../../docs/src/user/rpo_planner_pilot.md) shows how to run
it in a simulation.

## MuJoCo

```sh
julia --startup-file=no scripts/optional_packages.jl setup mujoco /path/to/mujoco-env
julia --startup-file=no --project=/path/to/mujoco-env
```

Then load `SpaceAGORA` and `SpaceAGORAMuJoCo`. This route hands the companion's
installer the SpaceAGORA checkout you run it from, in place of the core
revision the companion records for its standalone installation. Native scenes
need Linux x86_64 and use the MuJoCo 3.11.0 artifact, which is downloaded the
first time a native operation needs it.
[Orbital scenes with MuJoCo](../../docs/src/user/mujoco_scenes.md) covers
frames, ownership and the limits of the coupling.

Each `setup` initializes only its own submodule and checks the source before it
installs anything. Neither one loads a companion into a Julia session that is
already running. The [package overview](../README.md) describes `status`,
`checkout` and `check`, and why a Git clone is needed.
