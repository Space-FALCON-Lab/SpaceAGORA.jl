# Contributor starters

Three repository starters for work built on SpaceAGORA, each tested by core CI
(the `contributor-kit` job generates every starter from the commit under test and
runs it through `.github/workflows/contributor-ci.yml`):

| Kind | For | Example it ships |
|---|---|---|
| `project` | a research study on the public API | a force model (`wrench`) and a short case |
| `package` | an algorithm package for an existing family | an RPO reference planner |
| `companion` | a companion that serves a core service facade | an external propagator (`ExternalPropagationServices`) |

Use Julia 1.12 and run the generator from the root of a clean SpaceAGORA
checkout at the commit you want to build on. The output directory must not
already contain files:

```bash
julia --startup-file=no scripts/new_contributor_repository.jl package ../MyPlanner --name MyPlanner --holder "Your Lab"
```

The generator fills in the placeholders below and refuses a checkout with tracked
changes, a non-empty or in-checkout directory, and a name that is not a package name.

| Placeholder | Becomes |
|---|---|
| `StarterName` (in text and file names) | the `--name` |
| `00000000-0000-0000-0000-000000000000` | a fresh UUID |
| `0000000000000000000000000000000000000000` | the checkout's `HEAD` commit (pin file, README checkout command and workflow call) |
| `SpaceAGORA = "0.0"` | the core's release series in `[compat]` |
| `version = "0.0.0"` (pin file) | the core's version at that commit |
| `CORE_CONTRACT_VERSION = v"0.0.0"` | the facade's current contract version |
| `Copyright (c) 0000 Your name or organization` (`LICENSE`) | this year and `--holder`; without `--holder`, the holder stays for the author to fill |

The starters use only current hook names (the deprecated-hook gate refuses any
deprecated name under `templates/starters/`) and only root names or
`SpaceAGORA.Conformance` (the contributor workflow scans projects and packages for
internal names). They use no GNC tick hook. A conformance pass certifies interface
behaviour; validation of a study or a method stays with its author.
