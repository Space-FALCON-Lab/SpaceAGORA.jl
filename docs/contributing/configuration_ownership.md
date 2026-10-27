# Configuration source ownership

The one-run configuration types live in `src/orchestration/simulation/config/`.
Applications use the documented types exported by `SpaceAGORA`; internally they
are defined in `SpaceAGORA.SimulationModel` submodules, and no wrapper type
sits between the two. `configuration.jl` assembles the `SimConfig` module from
four of the files below.

| Source file under `src/orchestration/simulation/config/` | Responsibility |
| --- | --- |
| `run_settings.jl` | Mission type, epoch (`InitialTime`), file paths, `SimulationSettings` (output, checkpoint and resume settings) and `MissionConfiguration` (duration or orbit count, sampling) |
| `solver_settings.jl` | `SolverConfig` and `IntegrationTolerances` |
| `environment_settings.jl` | `EnvironmentModel`, which selects and composes the environment models without evaluating their physics |
| `simulation_configuration.jl` | The final `SimulationConfiguration` container and the `_with_fields` and `_with_configuration` copy helpers |
| `constellation_configuration.jl` | `DynamicsModel` (the member spacecraft and the selected dynamic effectors) and `scoped`, for an effector that acts on some spacecraft only |
| `example_configuration.jl` | The `ExampleConfiguration` module (`make_example_config`, `make_three_body_spacecraft`) |
| `presets/no_gram_presets.jl` | The `NoGramPresets` module: planets, density models and environments for runs without GRAM |

`constellation_configuration.jl` is included inside the `SpacecraftModels`
module, so that `DynamicsModel` keeps its module identity; the four files that
`configuration.jl` includes define `SimConfig`. The composition loader
(`src/orchestration/simulation/composition/simulation_model.jl`) includes
`example_configuration.jl` and the presets. A one-spacecraft run and a
constellation use the same collection. The `DynamicsModel` constructor keeps
the spacecraft vector it is given: it does not generate orbital layouts (that
is `build_constellation`, in `src/models/vehicle/constellations/`), schedule
activities or add collection checks of its own.

Execution policy is a separate responsibility, in
`src/orchestration/simulation/engine/config/`: `SimulationEngineConfig`
composes the parallel, solver, runtime-policy and artifact settings, and uses
the same `SolverConfig` type rather than a second solver type. Output and
checkpoint paths are derived in `src/interfaces/io/paths/`, and the solver's
environment variables are parsed in
`src/orchestration/simulation/engine/adapters/from_env.jl`. A configuration's
own `solver_config` takes precedence over those environment variables; they
apply only when it is `nothing`.

`_with_configuration` makes a shallow copy with named fields replaced and
keeps the references of every other field. Any other configuration struct is
copied the same way with `_with_fields(x; changes...)`: name only the fields
that change, and every other field passes through the type's own constructor,
so no copy lists another instance's fields one by one. Isolating mutable run
state is the engine's job (`run_simulation`'s `isolate_state`), not the copy
helpers'.
