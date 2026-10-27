# Simulation

| Folder or file | What it does |
|---|---|
| `config/` | The physical scenario and its declared run settings, including the onboarding presets in `config/presets/`. |
| `engine/` | Prepares and integrates a run. `engine/config/` assembles the execution settings and reuses the scenario's solver configuration type. |
| `callbacks/` | Scheduled and event-driven updates during a run. |
| `campaigns/` | Repeated runs: Monte Carlo, constellation ensembles and their results. |
| `composition/` | Builds the `SimulationModel` module from files across `src/`. |
| `state/` | Concrete run buffers and caches (`ConfigTypes`). |
| `callback_services.jl` | `CallbackServices`: the internal state and sampling operations callbacks may use. `engine/callback_services.jl` implements them with the engine's own functions, so callbacks never look up the engine module at run time. Solver policy, frame conventions and mutation rules stay with the engine. |
| `external_propagation_services.jl` | The public `ExternalPropagationServices` facade for companions that integrate some spacecraft themselves. |
| `runtime_services.jl` | `RuntimeServices`: bookkeeping for finding the GRAMSuite package on the load path. |

`SimulationModel`, `ConfigTypes` and the other modules defined here keep stable
names, so callers and Julia serialization see the same types.

Mission policies are in `../mission/`; parallel routing, cost estimates and
worker pools in `../parallel/`. Result files and checkpoint storage belong to
`src/interfaces/io/`.
