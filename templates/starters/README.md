# Contributor starters

Three starting points for a repository of your own built on SpaceAGORA. Core CI
tests all three: its `contributor-kit` job generates each starter from the
commit under test and runs it through the reusable workflow
`.github/workflows/contributor-ci.yml`.

| Kind | For | Example it ships |
|---|---|---|
| `project` | a research study on the public API | a force model (`wrench`) and a short case |
| `package` | an algorithm package for an existing family | an RPO reference planner |
| `companion` | a companion that serves one of core's service facades | an external propagator (`ExternalPropagationServices`) |

Use Julia 1.12 and run the generator from the root of a SpaceAGORA checkout
with no uncommitted changes to tracked files, at the commit you want to build
on. The target directory must be empty or not exist yet, and must be outside
the checkout:

```bash
julia --startup-file=no scripts/new_contributor_repository.jl package ../MyPlanner --name MyPlanner --holder "Your Lab"
```

The generator refuses a name that is not a valid package name, `SpaceAGORA`
itself, a name ending in `jl`, and any name containing the placeholder
`StarterName`. It fills in these placeholders:

| Placeholder | Becomes |
|---|---|
| `StarterName` (in text and file names) | the `--name` |
| `00000000-0000-0000-0000-000000000000` | a fresh UUID |
| `0000000000000000000000000000000000000000` | the checkout's `HEAD` commit (pin file, README checkout command and workflow call) |
| `SpaceAGORA = "0.0"` | the core's release series in `[compat]` |
| `version = "0.0.0"` (pin file) | the core's version at that commit |
| `CORE_CONTRACT_VERSION = v"0.0.0"` | the facade's current contract version (companion) |
| `Copyright (c) 0000 Your name or organization` (`LICENSE`) | this year and `--holder`; without `--holder` the holder is left for the author to fill in |

Only the package and companion starters ship a `LICENSE`. A project gets a copy
of the core commit's `Manifest.toml` instead, so every result can be traced to
exact package versions.

The starters never name a hook that SpaceAGORA has removed;
`test/gates/ci_removed_hook_names_gate.jl` checks every file under
`templates/starters/`, README files and workflows included. They use no GNC
tick hook. The project and package starters use only exported names and
`SpaceAGORA.Conformance`; the companion also uses its facade,
`SpaceAGORA.ExternalPropagationServices`. For projects and packages the
contributor workflow flags `SpaceAGORA.SimulationModel`,
`SpaceAGORA.SimulationEngine` and names that start with an underscore.

Passing the conformance checks shows that a repository meets SpaceAGORA's
interfaces. Validating the study or method is the author's work.
