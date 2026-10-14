# Optional packages

This directory makes SpaceAGORA's optional integrations easy to find without
making them required dependencies. The domain folders organize source checkouts;
each external repository still owns its implementation, releases and tests.
Loading SpaceAGORA alone does not activate these packages.

| Location | Role |
| --- | --- |
| [robotics/](robotics/README.md) | Pinned HYPR and MuJoCo repositories for motion planning and mechanical scenes. HYPR supports both spacecraft RPO and robot-arm planning. |
| [gnc/](gnc/README.md) | Pinned SpaceAGORAEDG repository for experimental energy-depletion guidance; core retains simulation lifecycle and actuation. |
| [VisualizerSource.toml](VisualizerSource.toml) | Source pin for the separate SpaceAGORAVisualizer installer; the renderer is not stored in this tree. |
| [SpaceAGORAHYPR/](SpaceAGORAHYPR/README.md) | Existing compatibility package that loads SpaceAGORA and HYPR together. Its location and role are unchanged. |

GRAMSuite keeps its submodule at `data/GRAMSuite.jl`. The viewer implementation
lives in SpaceAGORAVisualizer; core retains the scene contract and compatibility
entrypoints.

## Choose an integration

Run these commands from a SpaceAGORA Git clone, with Julia 1.12 and Git available:

```sh
julia --startup-file=no scripts/optional_packages.jl status
julia --startup-file=no scripts/optional_packages.jl checkout hypr
julia --startup-file=no scripts/optional_packages.jl check hypr
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
```

Replace `hypr` with `mujoco` or `edg` to select another companion. `status` reports the
optional checkouts. `checkout` acquires only the selected repository at the
recorded Git revision, without recursively initializing other submodules.
`check` verifies its revision, package identity, version and clean state; it is
a source check, not a simulation or native-library test. `setup` performs the
selected checkout and verification, then delegates to that package's existing
installer with this SpaceAGORA checkout. The destination must be a separate
Julia environment, outside the source checkouts. Package installation may need
network access.

The helper does not reset a modified or mismatched checkout. Inspect and preserve
local work before changing it. Keep the source checkouts in place while using an
environment that refers to their paths. See the [installation guide](../docs/src/user/installation_environment.md)
for loading each integration and for standalone installation routes.

## Pins and source archives

The parent Git commit records each external source revision as a submodule
gitlink. Review updates to these pins together with the compatibility evidence;
the helper does not follow a moving branch to choose a newer package.

GitHub's source ZIP does not contain submodule contents or the Git metadata used
to verify these pins. Use a Git clone for the shared commands. Alternatively,
follow a companion's standalone installation guide. A folder's presence alone
does not mean that the package has been installed or loaded in Julia.
