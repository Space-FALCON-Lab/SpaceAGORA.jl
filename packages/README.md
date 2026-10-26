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
| [HYPRSource.toml](HYPRSource.toml) | Immutable HYPR source pin for the standalone installer; it must agree with the external HYPR gitlink. |

GRAMSuite keeps its submodule at `data/GRAMSuite.jl`. Browser rendering lives in
[SpaceAGORAVisualizer.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl),
which retains the renderer's attribution and third-party notices. SpaceAGORA
owns saved simulation data, scene metadata and portable viewer JSON; the
companion turns that data into HTML. See the
[mission workflow](../docs/src/user/mission_workflows.md) for a new run and the
[visualization guide](../docs/src/user/visualization.md) for saved results.

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

The source records serve different installation routes. `OptionalPackages.toml`
describes the HYPR, MuJoCo and EDG checkouts; their Git revisions are recorded by
gitlinks. `HYPRSource.toml` also supplies HYPR's revision to its
standalone installer, while `VisualizerSource.toml` pins the separate
visualizer installer. Keep these routes compatible when updating a companion.
Combining the records requires migrating their readers and validating both clone
and standalone installation; removing a record alone would break a supported
route.

The parent Git commit records each external source revision as a submodule
gitlink. Review updates to these pins together with the compatibility evidence;
the helper does not follow a moving branch to choose a newer package.

GitHub's source ZIP does not contain submodule contents or the Git metadata used
to verify these pins. Use a Git clone for the shared commands. Alternatively,
follow a companion's standalone installation guide. A folder's presence alone
does not mean that the package has been installed or loaded in Julia.
