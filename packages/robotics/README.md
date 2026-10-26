# Robotics and motion-planning packages

These repositories provide optional planning and mechanical-scene capabilities.
They are checked out here at reviewed revisions, while their code and releases
remain owned by their separate repositories. The directory is an organizational
home, not a restriction on where the algorithms can be used: HYPR serves both
spacecraft rendezvous/proximity operations and robot-arm planning.

| Package | Revision recorded by this layout | Role |
| --- | --- | --- |
| [HYPR.jl](https://github.com/Space-FALCON-Lab/HYPR.jl/tree/2816df11c86206cc348ec2c4f82732f958f18ad2) | `2816df11c86206cc348ec2c4f82732f958f18ad2` | Search and trajectory planning, with an optional SpaceAGORA extension using `HYPRServices` contract 3.0.0. |
| [SpaceAGORAMuJoCo.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/tree/7b5a00318f6883d7cdc7c7e7c0c18830f96ce57f) | `7b5a00318f6883d7cdc7c7e7c0c18830f96ce57f` | Standalone orbital proximity scenes and spacecraft propagation through `ExternalPropagationServices` contract 2.0.0. |

From the SpaceAGORA repository root, prepare the integration you need:

```sh
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
julia --startup-file=no --project=/path/to/hypr-env
```

Then load `SpaceAGORA` and `HYPR` directly. Loading both activates the adapter
implemented by HYPR; no separate loader package is needed. The [HYPR user guide](../../docs/src/user/hypr_planning.md)
explains the planner and the [RPO guide](../../docs/src/user/rpo_planner_pilot.md)
shows its supported simulation entry points.

For MuJoCo, select `mujoco` and a separate environment:

```sh
julia --startup-file=no scripts/optional_packages.jl setup mujoco /path/to/mujoco-env
julia --startup-file=no --project=/path/to/mujoco-env
```

Then load `SpaceAGORA` and `SpaceAGORAMuJoCo`. This route gives the companion
installer the current core checkout instead of its standalone default core pin.
Native scenes currently require Linux x86_64 and use the lazy MuJoCo 3.11.0
artifact. See [Orbital scenes with MuJoCo](../../docs/src/user/mujoco_scenes.md)
for frames, ownership and coupling limits.

Both setup commands acquire only their selected submodule and verify the source
before installing dependencies. Neither imports a provider into an existing
simulation process. The [package overview](../README.md) describes `status`,
`checkout` and `check`, Git-clone requirements and source-archive behavior.
