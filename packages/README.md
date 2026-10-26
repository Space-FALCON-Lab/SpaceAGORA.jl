# Optional packages

The optional companions live in their own repositories, with their own
releases and tests. This folder pins the revisions of them that this SpaceAGORA
commit works with, as Git submodules grouped by domain, and records where the
separate viewer comes from. `using SpaceAGORA` loads none of them.

| Location | What it holds |
| --- | --- |
| [robotics/](robotics/README.md) | HYPR (motion planning for both spacecraft rendezvous and robot arms) and SpaceAGORAMuJoCo (contact scenes). |
| [gnc/](gnc/README.md) | SpaceAGORAEDG, the energy-depletion guidance kernels. |
| `OptionalPackages.toml` | Name, path, repository, UUID and version of the three companions above. Their revisions are the submodule gitlinks. |
| `HYPRSource.toml` | HYPR's repository and revision for the optional-package helper and its compatibility shortcut `scripts/setup_hypr.jl`. `scripts/optional_packages.jl` refuses to run if it differs from the HYPR gitlink. |
| `VisualizerSource.toml` | The SpaceAGORAVisualizer revision that `scripts/setup_visualizer.jl` installs. The viewer is not a submodule here. |

The GRAMSuite wrapper is a submodule too, at `data/GRAMSuite.jl`. SpaceAGORA
writes the saved results, scene metadata and viewer JSON;
[SpaceAGORAVisualizer.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl)
turns them into an HTML page and carries the renderer's attribution and
third-party notices. See [Mission Workflows](../docs/src/user/mission_workflows.md)
to set up a run and [Interactive Visualization](../docs/src/user/visualization.md)
to view saved results.

## Install a companion

From the root of a SpaceAGORA Git clone, with Julia 1.12 and Git:

```sh
julia --startup-file=no scripts/optional_packages.jl status
julia --startup-file=no scripts/optional_packages.jl checkout hypr
julia --startup-file=no scripts/optional_packages.jl check hypr
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
```

Use `mujoco` or `edg` in place of `hypr` for the other companions.

- `status` lists each companion's path and pinned revision and says whether it
  is checked out. It makes no network requests and writes nothing.
- `checkout` initializes only the selected submodule, at its pinned revision,
  without touching the others.
- `check` confirms that the checkout is at the pinned revision, comes from the
  declared repository, has the declared name, UUID and version, and has no local
  changes. It does not run the package or its tests.
- `setup` does `checkout` and `check`, then runs the selected companion's
  `scripts/setup_env.jl` with this SpaceAGORA checkout. HYPR, MuJoCo and EDG
  use the same generated installation policy. The environment directory must
  be outside both source checkouts. Installing can download packages.

The helper never resets a checkout with local changes or at another revision;
it stops and asks you to sort that out. The environment it creates points at
the checkout by path, so leave the checkout in place while you use it. The
[installation guide](../docs/src/user/installation_environment.md) shows how to
load each companion and the standalone installation routes.

## Moving a pin

A companion moves to a new revision when its gitlink changes in a SpaceAGORA
commit; for HYPR, `HYPRSource.toml` must change with it. Review such a change
together with the companion's compatibility evidence. The helper never follows
a branch to pick up a newer revision.

GitHub's source ZIP leaves out submodule contents and the Git metadata needed
to check these pins, so the commands above need a Git clone. Without one,
follow a companion's own installation guide. A folder that exists here is not
necessarily installed or loaded in your Julia environment.
