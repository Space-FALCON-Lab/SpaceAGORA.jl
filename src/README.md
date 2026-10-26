# Source directory

SpaceAGORA organizes its implementation into four responsibility groups. The
child directories keep their domain ownership; the parent folders make the
command surfaces, execution machinery, models and shared utilities easier to
find.

```text
src/
├── SpaceAGORA.jl
├── precompile_workload.jl
├── interfaces/
│   ├── cli/
│   ├── io/
│   └── analysis/
├── orchestration/
│   ├── simulation/
│   ├── mission/
│   └── parallel/
├── models/
│   ├── vehicle/
│   ├── environment/
│   ├── dynamics/
│   └── gnc/
└── foundation/
    ├── core/
    └── assets/
```

**Interfaces** contains the CLI, output paths, result formats and checkpoint
storage, and analysis and visualization support. It includes the operations used to enter a run or inspect
its results; it is not a registry of every abstract interface in the package.
Model-specific hooks stay with their model owners.

**Orchestration** prepares and executes simulations, dispatches callbacks and
campaigns, selects mission policy and routes parallel work. **Models** holds
vehicle construction, environmental models, dynamical equations and the GNC
contracts and algorithms used by those runs. **Foundation** holds shared model
contracts, frame and numerical utilities, and reusable asset lookup. Concrete
run buffers and caches belong to `orchestration/simulation/state/`.

These are folders, not new Julia namespace modules or strict dependency layers.
`SpaceAGORA.jl` remains the package entry point, with its package-level
`precompile_workload.jl` helper alongside it. The composition loader
`orchestration/simulation/composition/simulation_model.jl` builds
`SimulationModel` by including
definitions from several groups. The grouping preserves module identities,
public names and simulation behavior; it does not change which module owns a
type or function. Prefer the supported package API over including source files
from a script.

Within simulation, `config/` owns scenario settings and presets, `engine/`
owns execution, and `callbacks/` consumes the internal operations declared in
`callback_services.jl` and implemented by the engine. Result naming lives in
`interfaces/io/paths/`; checkpoint file storage lives in
`interfaces/io/checkpoints/`. The engine still decides when to save and resume.

Planning sources are grouped under `models/gnc/planning/`: RPO, robot-arm
references, and shared geometry and tree operations. The core-owned adapters of
the optional companions are in `models/gnc/companions/`.
Navigation retains its estimator hook; reference and collision geometry belongs
to planning. Completed planner-comparison reports are implemented in
`interfaces/analysis/planning/`, while historical function bindings remain
available from GNC. HYPR search and the typed EDG numerical kernels live in their
optional companion packages. Core retains their supported services and simulation
integration, including EDG model types, lifecycle and actuation.

Environment sources supply atmospheric, ephemeris and gravitational-field
quantities. Dynamics sources adapt those quantities to spacecraft forces and
torques. Their existing module identities and mathematical kernels are shared
across the physical source boundary.

Julia package extensions remain in the repository-root [`ext/`](../ext/), where
Julia discovers them. Optional companion repositories such as HYPR and
SpaceAGORAMuJoCo keep their own source trees outside this directory. Loading
SpaceAGORA and HYPR directly activates HYPR's extension. Runnable scenarios remain in [`examples/`](../examples/),
operational tools in [`scripts/`](../scripts/) and browser code in
[SpaceAGORAVisualizer.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl).

See [Architecture and Responsibilities](../docs/src/user/architecture.md) for
the runtime flow, detailed source owners and guidance on where new work belongs.
