# Configuration source ownership

The configuration implementation is assembled by
`src/orchestration/simulation/config/configuration.jl`. Applications use the
documented types exported by `SpaceAGORA`; their existing internal
`SpaceAGORA.SimulationModel` bindings retain their identities. This source-file
organization introduces no new configuration wrapper.

| Source file under `src/orchestration/simulation/config/` | Responsibility |
| --- | --- |
| `run_settings.jl` | Epoch, mission duration/orbits, sampling, paths, output and checkpoint settings |
| `solver_settings.jl` | `SolverConfig` and `IntegrationTolerances` |
| `environment_settings.jl` | Select and compose environmental models |
| `constellation_configuration.jl` | Existing `DynamicsModel`: member spacecraft and selected dynamic effectors |
| `simulation_configuration.jl` | Final container and `_with_configuration` helper |

`constellation_configuration.jl` is included inside the existing `SpacecraftModels`
module to preserve the identity of `DynamicsModel`; the other definitions remain
inside `SimConfig`. A one-spacecraft run and a constellation use the same collection.
The constructor retains the supplied spacecraft vector. It does not generate orbital
layouts, schedule activities or introduce additional collection validation.

Execution policy remains in `src/orchestration/simulation/engine/config/`: `SimulationEngineConfig`
composes parallel, solver, runtime-policy and artifact settings. It uses the same
`SolverConfig` definition, not a second solver type. Output and checkpoint path
derivation stays in `src/interfaces/io/paths/`; solver environment settings are parsed in
`src/orchestration/simulation/engine/adapters/from_env.jl`. These are distinct responsibilities
from assembling a scenario.

`_with_configuration` makes a shallow update and preserves unspecified references.
Any other configuration struct is copied the same way, with `_with_fields(x; changes...)`,
which names only the fields that change and passes every other field through the type's
constructor, so no copy lists another instance's fields one by one.
Runtime state isolation remains the engine's responsibility. Moving the definitions
changes neither those semantics nor typed-solver precedence over environment settings.
